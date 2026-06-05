#!/usr/bin/env python3
"""
ChEMBL latest SQLite -> unified ChEMBL drug-disease indication CSV

What this script does
---------------------
1. Reads the official ChEMBL latest FTP directory.
2. Finds the latest ChEMBL SQLite archive, e.g. chembl_36_sqlite.tar.gz.
3. Downloads:
   - chembl_XX_sqlite.tar.gz
   - checksums.txt
   - LICENSE
   - README
   - REQUIRED.ATTRIBUTION
   - chembl_XX_release_notes.txt if available
   - schema_documentation.txt if available
4. Verifies the downloaded archive using checksums.txt when possible.
5. Extracts the SQLite database safely.
6. Inspects key ChEMBL tables.
7. Extracts records from the drug_indication table.
8. Joins molecule_dictionary and disease_dictionary where available.
9. Exports a provenance-aware CSV using the unified drug-disease schema.

Source interpretation
---------------------
ChEMBL drug_indication records represent curated drug indication / clinical-phase evidence.
They should be interpreted separately from AACT/ClinicalTrials.gov trial co-occurrence records.

Input source
------------
ChEMBL FTP latest directory:
    https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/

Main table used:
    drug_indication

Optional joined tables:
    molecule_dictionary
    disease_dictionary

Unified output columns
----------------------
drug_name
drug_identifier
drug_identifier_type
disease_or_condition_name
disease_or_condition_identifier
disease_or_condition_identifier_type
relationship_type
evidence_type
source
internal_source

ChEMBL-specific metadata columns
--------------------------------
drug_id
disease_id
disease_name
mesh_id
mesh_heading
max_phase

Requirements
------------
pip install pandas requests
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import sqlite3
import tarfile
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import urljoin

import pandas as pd
import requests


# ============================================================
# CONFIG
# ============================================================

BASE_URL = "https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/"


UNIFIED_OUTPUT_COLUMNS = [
    "drug_name",
    "drug_identifier",
    "drug_identifier_type",
    "disease_or_condition_name",
    "disease_or_condition_identifier",
    "disease_or_condition_identifier_type",
    "relationship_type",
    "evidence_type",
    "source",
    "internal_source",
    "drug_id",
    "disease_id",
    "disease_name",
    "mesh_id",
    "mesh_heading",
    "max_phase",
]


# ============================================================
# GENERAL HELPERS
# ============================================================

def now_utc() -> str:
    return datetime.now(timezone.utc).isoformat()


def print_section(title: str) -> None:
    print("\n" + "=" * 100)
    print(title)
    print("=" * 100)


def human_size(num_bytes: int) -> str:
    size = float(num_bytes)
    for unit in ["B", "KB", "MB", "GB", "TB"]:
        if size < 1024:
            return f"{size:.2f} {unit}"
        size /= 1024
    return f"{size:.2f} PB"


def sha256_file(path: Path, chunk_size: int = 1024 * 1024) -> str:
    h = hashlib.sha256()

    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(chunk_size), b""):
            h.update(chunk)

    return h.hexdigest()


def download_file(
    url: str,
    out_path: Path,
    force: bool = False,
    timeout: int = 1200,
) -> Path:
    out_path.parent.mkdir(parents=True, exist_ok=True)

    if out_path.exists() and out_path.stat().st_size > 0 and not force:
        print(f"[CACHE] {out_path}")
        print(f"        Size:   {human_size(out_path.stat().st_size)}")
        print(f"        SHA256: {sha256_file(out_path)}")
        return out_path

    print(f"[DOWNLOAD] {url}")

    tmp = out_path.with_suffix(out_path.suffix + ".part")

    with requests.get(url, stream=True, timeout=timeout) as r:
        r.raise_for_status()

        with open(tmp, "wb") as f:
            for chunk in r.iter_content(chunk_size=1024 * 1024):
                if chunk:
                    f.write(chunk)

    tmp.replace(out_path)

    print(f"[DOWNLOADED] {out_path}")
    print(f"             Size:   {human_size(out_path.stat().st_size)}")
    print(f"             SHA256: {sha256_file(out_path)}")

    return out_path


def fetch_text(url: str, timeout: int = 120) -> str:
    r = requests.get(url, timeout=timeout)
    r.raise_for_status()
    return r.text


def clean_text_series(s: pd.Series) -> pd.Series:
    return s.fillna("").astype(str).str.strip()


# ============================================================
# ChEMBL LATEST DETECTION
# ============================================================

def find_latest_sqlite_filename(base_url: str = BASE_URL) -> tuple[str, str]:
    """
    Finds the latest chembl_XX_sqlite.tar.gz from the ChEMBL latest directory.

    Returns:
        filename, version
        e.g. ("chembl_36_sqlite.tar.gz", "36")
    """
    print_section("1. Detecting latest ChEMBL SQLite archive")

    html = fetch_text(base_url)
    matches = re.findall(r'href="(chembl_(\d+)_sqlite\.tar\.gz)"', html)

    if not matches:
        raise RuntimeError(
            f"Could not find chembl_XX_sqlite.tar.gz in {base_url}. "
            "Check the FTP directory format."
        )

    matches = sorted(matches, key=lambda x: int(x[1]), reverse=True)
    filename, version = matches[0]

    print(f"[LATEST SQLITE] {filename}")
    print(f"[VERSION] ChEMBL {version}")
    print(f"[URL] {urljoin(base_url, filename)}")

    return filename, version


def expected_files_for_version(version: str) -> list[str]:
    return [
        "checksums.txt",
        "LICENSE",
        "README",
        "REQUIRED.ATTRIBUTION",
        f"chembl_{version}_release_notes.txt",
        "schema_documentation.txt",
    ]


def parse_checksums(checksums_path: Path) -> dict[str, str]:
    """
    Parse ChEMBL checksums file.

    Expected patterns include:
        <sha256>  <filename>
    """
    checksums: dict[str, str] = {}

    if not checksums_path.exists():
        return checksums

    with open(checksums_path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.strip()

            if not line or line.startswith("#"):
                continue

            parts = re.split(r"\s+", line)

            if len(parts) < 2:
                continue

            sha = None
            fname = None

            for p in parts:
                if re.fullmatch(r"[a-fA-F0-9]{64}", p):
                    sha = p.lower()
                elif (
                    "chembl_" in p
                    or p in {
                        "LICENSE",
                        "README",
                        "REQUIRED.ATTRIBUTION",
                        "schema_documentation.txt",
                    }
                ):
                    fname = Path(p).name

            if sha and fname:
                checksums[fname] = sha

    return checksums


def verify_checksum_if_available(path: Path, checksums: dict[str, str]) -> dict:
    actual = sha256_file(path)
    expected = checksums.get(path.name)

    result = {
        "file": path.name,
        "actual_sha256": actual,
        "expected_sha256": expected,
        "verified": None,
    }

    if expected is None:
        print(f"[CHECKSUM] No expected SHA256 found for {path.name}; actual={actual}")
        result["verified"] = "not_available"
        return result

    if actual.lower() != expected.lower():
        result["verified"] = False

        raise RuntimeError(
            f"Checksum mismatch for {path.name}\n"
            f"Expected: {expected}\n"
            f"Actual:   {actual}"
        )

    print(f"[CHECKSUM PASS] {path.name}")
    result["verified"] = True

    return result


# ============================================================
# TAR EXTRACTION
# ============================================================

def inspect_tar_and_find_sqlite(tar_path: Path) -> tarfile.TarInfo:
    print_section("3. Inspecting ChEMBL SQLite tar.gz")

    with tarfile.open(tar_path, "r:gz") as tar:
        members = tar.getmembers()

    print(f"[TAR MEMBERS] {len(members):,}")
    print("[FIRST 30 MEMBERS]")

    for m in members[:30]:
        print(f"  - {m.name} | {human_size(m.size)}")

    candidates = []

    for m in members:
        if not m.isfile():
            continue

        lname = m.name.lower()

        if lname.endswith((".db", ".sqlite", ".sqlite3")):
            candidates.append(m)

    if not candidates:
        raise RuntimeError("No .db/.sqlite/.sqlite3 file found in ChEMBL archive.")

    candidates.sort(key=lambda x: x.size, reverse=True)
    chosen = candidates[0]

    print(f"[SELECTED SQLITE MEMBER] {chosen.name}")
    print(f"[SQLITE SIZE IN TAR] {human_size(chosen.size)}")

    return chosen


def extract_sqlite_member(
    tar_path: Path,
    member: tarfile.TarInfo,
    out_dir: Path,
    force: bool = False,
) -> Path:
    print_section("4. Extracting SQLite database")

    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / Path(member.name).name

    if out_path.exists() and out_path.stat().st_size > 0 and not force:
        print(f"[CACHE] Using extracted DB: {out_path}")
        print(f"        Size:   {human_size(out_path.stat().st_size)}")
        print(f"        SHA256: {sha256_file(out_path)}")
        return out_path

    with tarfile.open(tar_path, "r:gz") as tar:
        src = tar.extractfile(member)

        if src is None:
            raise RuntimeError(f"Could not open tar member: {member.name}")

        tmp = out_path.with_suffix(out_path.suffix + ".part")

        with open(tmp, "wb") as f:
            while True:
                chunk = src.read(1024 * 1024)

                if not chunk:
                    break

                f.write(chunk)

        tmp.replace(out_path)

    print(f"[EXTRACTED] {out_path}")
    print(f"            Size:   {human_size(out_path.stat().st_size)}")
    print(f"            SHA256: {sha256_file(out_path)}")

    return out_path


# ============================================================
# SQLITE HELPERS
# ============================================================

def table_exists(conn: sqlite3.Connection, table: str) -> bool:
    q = "SELECT 1 FROM sqlite_master WHERE type='table' AND name=? LIMIT 1;"
    return conn.execute(q, (table,)).fetchone() is not None


def get_columns(conn: sqlite3.Connection, table: str) -> list[str]:
    df = pd.read_sql_query(f"PRAGMA table_info({table});", conn)
    return df["name"].astype(str).tolist()


def choose_first_present(cols: set[str], *candidates: str) -> str | None:
    for c in candidates:
        if c in cols:
            return c
    return None


def inspect_database(conn: sqlite3.Connection, head_n: int = 5) -> dict:
    print_section("5. Inspecting ChEMBL database tables")

    tables = pd.read_sql_query(
        "SELECT name FROM sqlite_master WHERE type='table' ORDER BY name;",
        conn,
    )["name"].astype(str).tolist()

    print(f"[TABLE COUNT] {len(tables):,}")
    print("[FIRST 60 TABLES]")

    for t in tables[:60]:
        print(f"  - {t}")

    required_tables = ["drug_indication"]
    useful_tables = [
        "drug_indication",
        "molecule_dictionary",
        "disease_dictionary",
    ]

    for t in required_tables:
        if t not in tables:
            raise RuntimeError(f"Required table missing: {t}")

    metadata = {
        "table_count": len(tables),
        "tables_preview": tables[:100],
        "inspected_tables": {},
    }

    for t in useful_tables:
        if t not in tables:
            print(f"[WARN] Table not found: {t}")
            continue

        cols = get_columns(conn, t)

        row_count = int(
            pd.read_sql_query(
                f"SELECT COUNT(*) AS n FROM {t};",
                conn,
            )["n"].iloc[0]
        )

        head = pd.read_sql_query(f"SELECT * FROM {t} LIMIT {head_n};", conn)

        print("\n" + "-" * 100)
        print(f"[TABLE] {t}")
        print(f"[ROWS] {row_count:,}")
        print(f"[COLUMNS] {cols}")
        print("[HEAD]")
        print(head.to_string(index=False))

        metadata["inspected_tables"][t] = {
            "row_count": row_count,
            "columns": cols,
            "head_preview": head.to_dict(orient="records"),
        }

    return metadata


# ============================================================
# SQL BUILDER
# ============================================================

def build_chembl_drug_indication_sql(conn: sqlite3.Connection) -> str:
    """
    Build a schema-aware SQL query for ChEMBL drug_indication.

    Handles version differences where possible.
    """
    if not table_exists(conn, "drug_indication"):
        raise RuntimeError("Missing required table: drug_indication")

    di_cols = set(get_columns(conn, "drug_indication"))

    md_exists = table_exists(conn, "molecule_dictionary")
    md_cols = set(get_columns(conn, "molecule_dictionary")) if md_exists else set()

    dd_exists = table_exists(conn, "disease_dictionary")
    dd_cols = set(get_columns(conn, "disease_dictionary")) if dd_exists else set()

    # -----------------------------
    # Drug identifier and join logic
    # -----------------------------
    drug_join = None
    drug_id_expr = None

    if "molecule_chembl_id" in di_cols:
        drug_id_expr = "di.molecule_chembl_id"

        if md_exists and "molecule_chembl_id" in md_cols:
            drug_join = "di.molecule_chembl_id = md.molecule_chembl_id"

    if (
        drug_id_expr is None
        and "molecule_id" in di_cols
        and md_exists
        and "molecule_id" in md_cols
    ):
        drug_join = "di.molecule_id = md.molecule_id"

        if "molecule_chembl_id" in md_cols:
            drug_id_expr = "md.molecule_chembl_id"
        else:
            drug_id_expr = "CAST(md.molecule_id AS TEXT)"

    if (
        drug_id_expr is None
        and "molregno" in di_cols
        and md_exists
        and "molregno" in md_cols
    ):
        drug_join = "di.molregno = md.molregno"

        if "chembl_id" in md_cols:
            drug_id_expr = "md.chembl_id"
        elif "molecule_chembl_id" in md_cols:
            drug_id_expr = "md.molecule_chembl_id"
        else:
            drug_id_expr = "CAST(md.molregno AS TEXT)"

    if drug_id_expr is None:
        if "molecule_chembl_id" in di_cols:
            drug_id_expr = "di.molecule_chembl_id"
        elif "molecule_id" in di_cols:
            drug_id_expr = "CAST(di.molecule_id AS TEXT)"
        elif "molregno" in di_cols:
            drug_id_expr = "CAST(di.molregno AS TEXT)"
        else:
            drug_id_expr = "''"

    # -----------------------------
    # Drug name
    # -----------------------------
    if md_exists and "pref_name" in md_cols:
        drug_name_expr = "COALESCE(md.pref_name, '')"
    elif "molecule_pref_name" in di_cols:
        drug_name_expr = "COALESCE(di.molecule_pref_name, '')"
    else:
        drug_name_expr = "''"

    # -----------------------------
    # Disease identifier and name
    # -----------------------------
    disease_id_expr = "di.efo_id" if "efo_id" in di_cols else "''"
    di_efo_term = "di.efo_term" if "efo_term" in di_cols else "''"

    disease_join = None

    if dd_exists and "efo_id" in di_cols and "efo_id" in dd_cols:
        disease_join = "di.efo_id = dd.efo_id"

        if "pref_name" in dd_cols:
            disease_name_expr = f"COALESCE(dd.pref_name, {di_efo_term}, '')"
        else:
            disease_name_expr = f"COALESCE({di_efo_term}, '')"
    else:
        disease_name_expr = f"COALESCE({di_efo_term}, '')"

    # -----------------------------
    # MeSH and phase fields
    # -----------------------------
    mesh_id_expr = "di.mesh_id" if "mesh_id" in di_cols else "''"
    mesh_heading_expr = "di.mesh_heading" if "mesh_heading" in di_cols else "''"

    phase_col = choose_first_present(
        di_cols,
        "max_phase_for_ind",
        "max_phase",
        "max_phase_for_indication",
    )

    phase_expr = f"di.{phase_col}" if phase_col else "0"

    sql = f"""
        SELECT
            {drug_id_expr} AS drug_id,
            {drug_name_expr} AS drug_name,
            {disease_id_expr} AS disease_id,
            {disease_name_expr} AS disease_name,
            {mesh_id_expr} AS mesh_id,
            {mesh_heading_expr} AS mesh_heading,
            {phase_expr} AS max_phase
        FROM drug_indication di
    """

    if md_exists and drug_join:
        sql += f"\nLEFT JOIN molecule_dictionary md ON {drug_join}\n"

    if dd_exists and disease_join:
        sql += f"LEFT JOIN disease_dictionary dd ON {disease_join}\n"

    return sql


# ============================================================
# EXPORT
# ============================================================

def export_chembl_drug_disease(
    conn: sqlite3.Connection,
    out_csv: Path,
    internal_source: str,
    chunksize: int = 200_000,
) -> dict:
    print_section("6. Exporting ChEMBL drug-disease indication CSV")

    sql = build_chembl_drug_indication_sql(conn)

    print("[SQL USED]")
    print(sql)

    preview = pd.read_sql_query(sql + " LIMIT 10", conn)

    print("\n[OUTPUT PREVIEW BEFORE CLEANING]")
    print(preview.to_string(index=False))

    out_csv.parent.mkdir(parents=True, exist_ok=True)

    if out_csv.exists():
        out_csv.unlink()

    total_before_filter = 0
    total_after_filter = 0
    wrote_header = False

    missing_drug_name_before_filtering = 0
    missing_disease_id_before_filtering = 0
    missing_mesh_id_before_filtering = 0

    for chunk in pd.read_sql_query(sql, conn, chunksize=chunksize):
        total_before_filter += len(chunk)

        for c in [
            "drug_id",
            "drug_name",
            "disease_id",
            "disease_name",
            "mesh_id",
            "mesh_heading",
        ]:
            chunk[c] = clean_text_series(chunk[c])

        chunk["max_phase"] = (
            pd.to_numeric(chunk["max_phase"], errors="coerce")
            .fillna(0)
            .astype(int)
        )

        missing_drug_name_before_filtering += int((chunk["drug_name"] == "").sum())
        missing_disease_id_before_filtering += int((chunk["disease_id"] == "").sum())
        missing_mesh_id_before_filtering += int((chunk["mesh_id"] == "").sum())

        # Keep records with a drug identifier and a disease/condition name.
        # We do not require disease_id because some versions/records may have usable disease names.
        chunk = chunk[
            (chunk["drug_id"] != "")
            & (chunk["disease_name"] != "")
        ].copy()

        if len(chunk) == 0:
            continue

        # ------------------------------------------------------------
        # Unified schema columns
        # ------------------------------------------------------------
        chunk["drug_identifier"] = chunk["drug_id"]
        chunk["drug_identifier_type"] = "ChEMBL"

        chunk["disease_or_condition_name"] = chunk["disease_name"]
        chunk["disease_or_condition_identifier"] = chunk["disease_id"]
        chunk["disease_or_condition_identifier_type"] = chunk["disease_id"].apply(
            lambda x: "EFO" if str(x).strip() else ""
        )

        chunk["relationship_type"] = "drug_indication"
        chunk["evidence_type"] = "curated_drug_indication"

        chunk["source"] = "ChEMBL"
        chunk["internal_source"] = internal_source

        chunk = chunk[UNIFIED_OUTPUT_COLUMNS]

        # Remove exact duplicate rows inside chunk.
        # Global exact duplicates are checked after writing.
        chunk = chunk.drop_duplicates()

        chunk.to_csv(
            out_csv,
            mode="a",
            header=(not wrote_header),
            index=False,
        )

        wrote_header = True
        total_after_filter += len(chunk)

        print(f"[WRITE] {total_after_filter:,} rows -> {out_csv}")

    if not out_csv.exists() or out_csv.stat().st_size == 0:
        raise RuntimeError("Output CSV was not created or is empty.")

    output_head = pd.read_csv(out_csv, nrows=10, dtype=str)

    full = pd.read_csv(
        out_csv,
        dtype=str,
        low_memory=False,
    ).fillna("")

    exact_duplicate_rows = int(full.duplicated().sum())

    unique_drug_disease_pairs = int(
        full[
            [
                "drug_identifier",
                "disease_or_condition_name",
            ]
        ]
        .drop_duplicates()
        .shape[0]
    )

    summary = {
        "output_file": str(out_csv),
        "output_size_bytes": out_csv.stat().st_size,
        "output_sha256": sha256_file(out_csv),
        "rows_before_filtering": int(total_before_filter),
        "rows_after_filtering": int(len(full)),
        "exact_duplicate_rows_in_output": exact_duplicate_rows,
        "unique_drug_identifiers": int(full["drug_identifier"].nunique(dropna=True)),
        "unique_drug_names": int(full["drug_name"].nunique(dropna=True)),
        "unique_disease_identifiers": int(
            full["disease_or_condition_identifier"].nunique(dropna=True)
        ),
        "unique_disease_names": int(
            full["disease_or_condition_name"].nunique(dropna=True)
        ),
        "unique_drug_disease_pairs": unique_drug_disease_pairs,
        "missing_drug_name_before_filtering": int(
            missing_drug_name_before_filtering
        ),
        "missing_disease_id_before_filtering": int(
            missing_disease_id_before_filtering
        ),
        "missing_mesh_id_before_filtering": int(
            missing_mesh_id_before_filtering
        ),
        "relationship_type_counts": (
            full["relationship_type"].value_counts().to_dict()
            if len(full)
            else {}
        ),
        "evidence_type_counts": (
            full["evidence_type"].value_counts().to_dict()
            if len(full)
            else {}
        ),
        "max_phase_counts": (
            full["max_phase"].value_counts().to_dict()
            if len(full)
            else {}
        ),
        "output_columns": list(full.columns),
        "output_head": output_head.to_dict(orient="records"),
    }

    print("\n[OUTPUT HEAD]")
    print(output_head.to_string(index=False))

    print("\n[OUTPUT SUMMARY]")
    for k, v in summary.items():
        if k != "output_head":
            print(f"{k}: {v}")

    return summary


# ============================================================
# MAIN
# ============================================================

def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Download latest ChEMBL SQLite and export unified ChEMBL "
            "drug-disease indication CSV."
        )
    )

    parser.add_argument(
        "--base_url",
        default=BASE_URL,
        help="ChEMBL latest FTP directory URL.",
    )

    parser.add_argument(
        "--workdir",
        default="./data/raw/chembl",
        help="Directory for downloaded ChEMBL files.",
    )

    parser.add_argument(
        "--out",
        default="./data/processed/chembl_drug_disease.csv",
        help="Output CSV path.",
    )

    parser.add_argument(
        "--metadata",
        default="./data/metadata/chembl_processing_metadata.json",
        help="Metadata JSON path.",
    )

    parser.add_argument(
        "--input_archive",
        default=None,
        help=(
            "Optional existing chembl_XX_sqlite.tar.gz. "
            "If provided, download is skipped for the archive."
        ),
    )

    parser.add_argument(
        "--keep_db",
        action="store_true",
        help=(
            "Keep extracted SQLite DB under workdir/extracted. "
            "Otherwise use a temporary directory."
        ),
    )

    parser.add_argument(
        "--force_download",
        action="store_true",
        help="Force re-download even if files exist.",
    )

    parser.add_argument(
        "--force_extract",
        action="store_true",
        help="Force re-extraction of SQLite database.",
    )

    parser.add_argument(
        "--chunksize",
        type=int,
        default=200_000,
        help="Rows per SQL chunk.",
    )

    parser.add_argument(
        "--head_n",
        type=int,
        default=5,
        help="Number of rows to print from table heads.",
    )

    args = parser.parse_args()

    workdir = Path(args.workdir).resolve()
    out_csv = Path(args.out).resolve()
    metadata_path = Path(args.metadata).resolve()

    workdir.mkdir(parents=True, exist_ok=True)

    metadata = {
        "script": "chembl_to_drug_disease.py",
        "source": "ChEMBL",
        "relationship_type": "drug_indication",
        "evidence_type": "curated_drug_indication",
        "base_url": args.base_url,
        "started_at_utc": now_utc(),
        "processing_interpretation": (
            "ChEMBL drug_indication records are curated drug indication "
            "records with clinical phase metadata. They are not equivalent "
            "to AACT/ClinicalTrials.gov intervention-condition co-occurrence."
        ),
        "unified_schema_columns": UNIFIED_OUTPUT_COLUMNS,
    }

    try:
        print_section("ChEMBL latest -> unified drug-disease indication CSV")

        if args.input_archive:
            archive_path = Path(args.input_archive).expanduser().resolve()

            if not archive_path.exists():
                raise FileNotFoundError(
                    f"Input archive does not exist: {archive_path}"
                )

            match = re.search(
                r"chembl_(\d+)_sqlite\.tar\.gz",
                archive_path.name,
            )

            if not match:
                raise ValueError(
                    "Input archive name should look like chembl_36_sqlite.tar.gz "
                    "so the ChEMBL version can be inferred."
                )

            version = match.group(1)
            sqlite_filename = archive_path.name

            print(f"[MANUAL ARCHIVE MODE] {archive_path}")

        else:
            sqlite_filename, version = find_latest_sqlite_filename(args.base_url)
            archive_url = urljoin(args.base_url, sqlite_filename)
            archive_path = workdir / sqlite_filename

            # Download support files first.
            print_section("2. Downloading ChEMBL support files")

            support_files = expected_files_for_version(version)
            downloaded_support = []

            for fname in support_files:
                url = urljoin(args.base_url, fname)
                out_path = workdir / fname

                try:
                    download_file(
                        url,
                        out_path,
                        force=args.force_download,
                    )
                    downloaded_support.append(str(out_path))

                except Exception as e:
                    print(f"[WARN] Could not download support file {fname}: {e}")

            metadata["downloaded_support_files"] = downloaded_support

            # Download main archive.
            archive_path = download_file(
                archive_url,
                workdir / sqlite_filename,
                force=args.force_download,
                timeout=3600,
            )

        internal_source = f"ChEMBL{version}"

        metadata["internal_source"] = internal_source
        metadata["chembl_version"] = version
        metadata["sqlite_archive"] = {
            "filename": sqlite_filename,
            "path": str(archive_path),
            "size_bytes": archive_path.stat().st_size,
            "sha256": sha256_file(archive_path),
        }

        # Checksum verification if checksums.txt exists.
        checksum_results = []
        checksums_path = workdir / "checksums.txt"

        if checksums_path.exists():
            checksums = parse_checksums(checksums_path)

            print_section("Checksum verification")

            checksum_results.append(
                verify_checksum_if_available(
                    archive_path,
                    checksums,
                )
            )
        else:
            print(
                "[WARN] checksums.txt not available; "
                "skipping checksum verification against expected values."
            )

        metadata["checksum_verification"] = checksum_results

        member = inspect_tar_and_find_sqlite(archive_path)

        metadata["selected_sqlite_member"] = {
            "name": member.name,
            "size_bytes": member.size,
        }

        if args.keep_db:
            extract_dir = workdir / "extracted"
            temp_ctx = None
        else:
            temp_ctx = tempfile.TemporaryDirectory(prefix="chembl_sqlite_")
            extract_dir = Path(temp_ctx.name)

        db_path = extract_sqlite_member(
            archive_path,
            member,
            extract_dir,
            force=args.force_extract,
        )

        metadata["extracted_sqlite_db"] = {
            "path": str(db_path),
            "size_bytes": db_path.stat().st_size,
            "sha256": sha256_file(db_path),
        }

        conn = sqlite3.connect(str(db_path))

        try:
            metadata["database_inspection"] = inspect_database(
                conn,
                head_n=args.head_n,
            )

            metadata["processing"] = export_chembl_drug_disease(
                conn=conn,
                out_csv=out_csv,
                internal_source=internal_source,
                chunksize=args.chunksize,
            )

        finally:
            conn.close()

        if temp_ctx is not None:
            temp_ctx.cleanup()

        metadata["finished_at_utc"] = now_utc()
        metadata["status"] = "success"

    except Exception as e:
        metadata["finished_at_utc"] = now_utc()
        metadata["status"] = "failed"
        metadata["error"] = repr(e)
        raise

    finally:
        metadata_path.parent.mkdir(parents=True, exist_ok=True)

        with open(metadata_path, "w", encoding="utf-8") as f:
            json.dump(metadata, f, indent=2)

        print_section("METADATA WRITTEN")
        print(metadata_path)


if __name__ == "__main__":
    main()