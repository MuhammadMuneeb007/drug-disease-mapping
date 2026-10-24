
#!/usr/bin/env python3
"""ChEMBL 37 drug-indication extractor (one-file, merge-compatible).

Run: python chembl.py

Reuses data/raw/chembl/chembl_37_sqlite.tar.gz.
Does NOT download by default.

Overwrites data/processed/chembl_drug_disease.csv only
after successful validation.

The ten common drug--disease fields remain first;
source evidence follows.

Indication-specific clinical phase is NOT approval of
every drug--disease pair.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import sqlite3
import tarfile
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import urljoin

import pandas as pd
import requests


# ============================================================
# PATHS AND SOURCE CONFIGURATION
# ============================================================

ROOT = Path(__file__).resolve().parent

RAW = ROOT / "data/raw/chembl"

OUT = ROOT / "data/processed/chembl_drug_disease.csv"

META = (
    ROOT
    / "data/metadata/chembl_processing_metadata.json"
)

BASE_URL = (
    "https://ftp.ebi.ac.uk/pub/databases/chembl/"
    "ChEMBLdb/releases/chembl_37/"
)

ARCHIVE_NAME = "chembl_37_sqlite.tar.gz"


# ============================================================
# UNIFIED MERGE-COMPATIBLE SCHEMA
# ============================================================

CORE = [
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
]

EXTRA = [
    "drug_id",
    "disease_id",
    "disease_name",
    "mesh_id",
    "mesh_heading",
    "max_phase_for_indication",
    "molecule_max_phase",
    "drugind_id",
    "record_id",
    "molregno",
    "efo_term",
    "indication_ref_types",
    "indication_ref_count",
    "indication_refs_original_json",
    "drug_name_lookup_status",
    "disease_name_lookup_status",
]

for entity in ("drug", "disease", "mesh"):
    EXTRA += [
        f"{entity}_identifier_{k}"
        for k in (
            "original",
            "namespace",
            "canonical",
            "validation",
        )
    ]

COLS = CORE + EXTRA


# ============================================================
# GENERAL UTILITIES
# ============================================================

def now():
    return datetime.now(timezone.utc).isoformat()


def section(title):
    print(
        "\n" + "=" * 92
        + "\n" + title
        + "\n" + "=" * 92,
        flush=True,
    )


def sha256(path):
    h = hashlib.sha256()

    with Path(path).open("rb") as fh:
        for part in iter(
            lambda: fh.read(4 * 1024 * 1024),
            b"",
        ):
            h.update(part)

    return h.hexdigest()


def clean(x):
    if x is None or pd.isna(x):
        return ""

    return str(x).strip()


def phase(x):
    s = clean(x)

    if not s:
        return ""

    try:
        n = float(s)

        return (
            str(int(n))
            if n.is_integer()
            else str(n)
        )

    except ValueError:
        return s


def sql_cols(conn, table):
    return [
        x[1]
        for x in conn.execute(
            f"PRAGMA table_info({table})"
        )
    ]


def tables(conn):
    return {
        x[0]
        for x in conn.execute(
            "SELECT name FROM sqlite_master "
            "WHERE type='table'"
        )
    }


# ============================================================
# IDENTIFIER NORMALISATION
# ============================================================

def normalise_id(value, declared=""):
    """
    Validate identifier syntax and optionally create
    a canonical CURIE.

    This does not establish ontology equivalence.
    Original identifier values are never overwritten.
    """

    raw = clean(value)

    if not raw:
        return "", "", "missing"

    token = raw.upper().replace(" ", "")
    declared = declared.upper()

    if (
        declared == "CHEMBL"
        or token.startswith("CHEMBL")
    ):
        m = re.fullmatch(
            r"(?:CHEMBL:)?(CHEMBL\d+)",
            token,
        )

        if m:
            return (
                "ChEMBL",
                f"CHEMBL:{m[1]}",
                "syntax_valid",
            )

        return (
            "ChEMBL",
            "",
            "invalid_syntax",
        )

    if (
        declared == "MESH"
        or token.startswith("MESH:")
    ):
        token = token.removeprefix("MESH:")

        if re.fullmatch(
            r"[CD]\d{6,9}",
            token,
        ):
            return (
                "MeSH",
                "MESH:" + token,
                "syntax_valid",
            )

        return (
            "MeSH",
            "",
            "invalid_syntax",
        )

    patterns = [
        (
            "EFO",
            "EFO",
            r"EFO[:_](\d+)",
            7,
        ),
        (
            "MONDO",
            "MONDO",
            r"MONDO[:_](\d+)",
            7,
        ),
        (
            "HP",
            "HPO",
            r"HP[:_](\d+)",
            7,
        ),
        (
            "DOID",
            "DOID",
            r"DOID[:_](\d+)",
            0,
        ),
        (
            "ORPHA",
            "Orphanet",
            r"(?:ORPHA|ORPHANET)[:_](\d+)",
            0,
        ),
        (
            "OMIM",
            "OMIM",
            r"OMIM[:_](\d+)",
            0,
        ),
        (
            "UMLS",
            "UMLS",
            r"(?:UMLS:)?(C\d{7,9})",
            0,
        ),
    ]

    for prefix, label, pat, width in patterns:

        candidate = token.startswith(
            (
                prefix + ":",
                prefix + "_",
            )
        )

        if (
            prefix == "UMLS"
            and re.fullmatch(
                r"C\d{7,9}",
                token,
            )
        ):
            candidate = True

        if (
            prefix == "ORPHA"
            and token.startswith("ORPHANET")
        ):
            candidate = True

        if not candidate:
            continue

        m = re.fullmatch(pat, token)

        if m:
            digits = (
                m[1].zfill(width)
                if width
                else m[1]
            )

            return (
                label,
                f"{prefix}:{digits}",
                "syntax_valid",
            )

        return (
            label,
            "",
            "invalid_syntax",
        )

    return (
        "unresolved",
        "",
        "unrecognised_namespace",
    )


def annotate(
    df,
    column,
    entity,
    declared="",
):
    mapped = [
        normalise_id(v, declared)
        for v in df[column].tolist()
    ]

    df[f"{entity}_identifier_original"] = (
        df[column]
    )

    df[f"{entity}_identifier_namespace"] = [
        m[0] for m in mapped
    ]

    df[f"{entity}_identifier_canonical"] = [
        m[1] for m in mapped
    ]

    df[f"{entity}_identifier_validation"] = [
        m[2] for m in mapped
    ]

    return df


# ============================================================
# DOWNLOAD AND ARCHIVE HANDLING
# ============================================================

def ensure_archive(
    path,
    url,
    allow_download,
    force_download=False,
):
    if (
        path.is_file()
        and path.stat().st_size
        and not force_download
    ):
        print(
            f"[REUSE PINNED ARCHIVE] {path}"
        )

        return

    if not allow_download:
        raise FileNotFoundError(
            f"Missing {path}; "
            "use --download_missing if needed"
        )

    path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    tmp = path.with_name(
        path.name + ".part"
    )

    try:
        with requests.get(
            url,
            stream=True,
            timeout=3600,
        ) as response:

            response.raise_for_status()

            with tmp.open("wb") as output:
                for part in response.iter_content(
                    4 * 1024 * 1024
                ):
                    if part:
                        output.write(part)

        if not tmp.stat().st_size:
            raise ValueError(
                "Downloaded archive is empty"
            )

        tmp.replace(path)

    finally:
        tmp.unlink(missing_ok=True)


def extract_db(
    archive,
    directory,
    force=False,
):
    section(
        "2. COMPRESSED ARCHIVE CONTENTS "
        "/ SQLITE EXTRACTION"
    )

    with tarfile.open(
        archive,
        "r:gz",
    ) as tar:

        members = tar.getmembers()

        print(
            f"TAR member count: {len(members):,}"
        )

        for member in members:
            if member.isfile():
                print(
                    f"  {member.name} | "
                    f"{member.size:,} bytes"
                )

        candidates = [
            m for m in members
            if m.isfile()
            and m.name.lower().endswith(
                (
                    ".db",
                    ".sqlite",
                    ".sqlite3",
                )
            )
        ]

        if not candidates:
            raise RuntimeError(
                "Archive has no SQLite database file"
            )

        member = max(
            candidates,
            key=lambda m: m.size,
        )

        directory.mkdir(
            parents=True,
            exist_ok=True,
        )

        dest = (
            directory
            / Path(member.name).name
        )

        if (
            dest.exists()
            and dest.stat().st_size == member.size
            and not force
        ):
            print(
                f"[REUSE EXTRACTED SQLITE] {dest}"
            )

            return dest, member.name

        tmp = dest.with_name(
            dest.name + ".part"
        )

        try:
            src = tar.extractfile(member)

            if src is None:
                raise RuntimeError(
                    "Cannot extract selected "
                    "SQLite member"
                )

            with src, tmp.open("wb") as target:

                while True:
                    block = src.read(
                        4 * 1024 * 1024
                    )

                    if not block:
                        break

                    target.write(block)

            if (
                tmp.stat().st_size
                != member.size
            ):
                raise RuntimeError(
                    "Extracted database "
                    "size mismatch"
                )

            tmp.replace(dest)

        finally:
            tmp.unlink(missing_ok=True)

    print(f"[EXTRACTED] {dest}")

    return dest, member.name


# ============================================================
# SQLITE SCHEMA INSPECTION
# ============================================================

def inspect_schema(conn, n=3):

    section(
        "3. SQLITE TABLES, ORIGINAL "
        "COLUMNS AND ROW COUNTS"
    )

    available = tables(conn)

    print(
        f"SQLite tables: {len(available):,}"
    )

    required = {
        "drug_indication": {
            "drugind_id",
            "molregno",
            "efo_id",
            "efo_term",
            "max_phase_for_ind",
        },
        "molecule_dictionary": {
            "molregno",
            "chembl_id",
            "pref_name",
        },
    }

    for table, cols in required.items():

        if table not in available:
            raise RuntimeError(
                "Missing required ChEMBL 37 "
                f"table: {table}"
            )

        missing = (
            cols
            - set(sql_cols(conn, table))
        )

        if missing:
            raise RuntimeError(
                f"{table} missing columns: "
                f"{sorted(missing)}"
            )

    counts = {}
    headers = {}

    selected_tables = (
        "drug_indication",
        "molecule_dictionary",
        "indication_refs",
        "disease_dictionary",
    )

    for table in selected_tables:

        if table not in available:
            print(
                f"  {table}: NOT PRESENT"
            )
            continue

        headers[table] = sql_cols(
            conn,
            table,
        )

        counts[table] = conn.execute(
            f"SELECT COUNT(*) FROM {table}"
        ).fetchone()[0]

        print(
            f"  {table}: "
            f"{counts[table]:,} rows | "
            f"{len(headers[table])} columns"
        )

        print(
            "    "
            + " | ".join(headers[table])
        )

        if n:
            rows = conn.execute(
                f"SELECT * FROM {table} LIMIT ?",
                (n,),
            ).fetchall()

            if rows:
                print(
                    "    Sample:",
                    dict(
                        zip(
                            headers[table],
                            rows[0],
                        )
                    ),
                )

    return counts, headers


# ============================================================
# SQL JOINS
# ============================================================

def sql_builder(headers):

    di = set(
        headers["drug_indication"]
    )

    md = set(
        headers["molecule_dictionary"]
    )

    has_dd = (
        "disease_dictionary" in headers
        and "efo_id"
        in headers["disease_dictionary"]
    )

    dd_cols = set(
        headers.get(
            "disease_dictionary",
            [],
        )
    )

    ref_id = (
        "di.record_id"
        if "record_id" in di
        else "NULL"
    )

    mesh_id = (
        "di.mesh_id"
        if "mesh_id" in di
        else "NULL"
    )

    mesh_label = (
        "di.mesh_heading"
        if "mesh_heading" in di
        else "NULL"
    )

    mol_phase = (
        "md.max_phase"
        if "max_phase" in md
        else "NULL"
    )

    dd_label = (
        "dd.pref_name"
        if (
            has_dd
            and "pref_name" in dd_cols
        )
        else "NULL"
    )

    dd_join = (
        "LEFT JOIN disease_dictionary dd "
        "ON dd.efo_id = di.efo_id"
        if has_dd
        else ""
    )

    sql = f"""
        SELECT
            di.drugind_id AS drugind_id,
            {ref_id} AS record_id,
            di.molregno AS molregno,

            md.chembl_id AS drug_id,
            md.pref_name AS preferred_name,

            di.efo_id AS disease_id,
            di.efo_term AS efo_term,

            {dd_label}
                AS dictionary_disease_name,

            {mesh_id} AS mesh_id,
            {mesh_label} AS mesh_heading,

            di.max_phase_for_ind
                AS max_phase_for_indication,

            {mol_phase} AS molecule_max_phase

        FROM drug_indication di

        LEFT JOIN molecule_dictionary md
            ON di.molregno = md.molregno

        {dd_join}
    """

    return sql, {
        "molecule_join":
            "di.molregno = md.molregno",

        "disease_dictionary_joined":
            bool(has_dd),

        "indication_phase":
            "drug_indication.max_phase_for_ind",

        "molecule_phase":
            mol_phase,
    }


def validate_joins(
    conn,
    sql,
    input_rows,
):

    section(
        "4. JOIN CARDINALITY AND "
        "SOURCE PROVENANCE CHECKS"
    )

    checks = {}

    for table, key in (
        (
            "drug_indication",
            "drugind_id",
        ),
        (
            "molecule_dictionary",
            "molregno",
        ),
    ):

        checks[
            table + "_duplicate_keys"
        ] = conn.execute(
            f"""
            SELECT COUNT(*)
            FROM (
                SELECT {key}
                FROM {table}
                GROUP BY {key}
                HAVING COUNT(*) > 1
            )
            """
        ).fetchone()[0]

        if checks[
            table + "_duplicate_keys"
        ]:
            raise RuntimeError(
                "Non-unique primary join key: "
                f"{table}.{key}"
            )

    checks["unmatched_molregno"] = (
        conn.execute(
            """
            SELECT COUNT(*)
            FROM drug_indication di

            LEFT JOIN molecule_dictionary md
                ON md.molregno = di.molregno

            WHERE md.molregno IS NULL
            """
        ).fetchone()[0]
    )

    checks["joined_rows"] = (
        conn.execute(
            f"SELECT COUNT(*) FROM ({sql})"
        ).fetchone()[0]
    )

    for name, value in checks.items():
        print(
            f"  {name}: {value:,}"
        )

    if checks["joined_rows"] != input_rows:
        raise RuntimeError(
            "ChEMBL SQL join multiplied "
            "or dropped source rows"
        )

    return checks


# ============================================================
# INDICATION REFERENCES
# ============================================================

def references(
    conn,
    output_path,
):

    available = tables(conn)

    if (
        "indication_refs" not in available
        or not {
            "drugind_id",
            "ref_type",
        }.issubset(
            sql_cols(
                conn,
                "indication_refs",
            )
        )
    ):
        print(
            "[REFS] No compatible "
            "indication_refs table"
        )

        return {}, 0, None

    ref_cols = sql_cols(
        conn,
        "indication_refs",
    )

    query = (
        "SELECT "
        + ", ".join(
            '"' + c + '"'
            for c in ref_cols
        )
        + " FROM indication_refs"
    )

    grouped = defaultdict(list)
    rows = 0

    output_path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    tmp = output_path.with_name(
        output_path.name + ".tmp"
    )

    try:
        with tmp.open(
            "w",
            encoding="utf-8",
            newline="",
        ) as handle:

            import csv

            writer = csv.writer(handle)
            writer.writerow(ref_cols)

            cur = conn.execute(query)

            while batch := cur.fetchmany(10000):

                for record in batch:

                    writer.writerow(record)

                    datum = dict(
                        zip(
                            ref_cols,
                            record,
                        )
                    )

                    key = clean(
                        datum.get("drugind_id")
                    )

                    grouped[key].append(
                        datum
                    )

                    rows += 1

    except Exception:
        tmp.unlink(missing_ok=True)
        raise

    print(
        f"[REFS] {rows:,} source reference "
        f"records, {len(grouped):,} "
        "indication IDs"
    )

    return grouped, rows, tmp


# ============================================================
# TRANSFORM ONE SQLITE CHUNK
# ============================================================

def transform(
    chunk,
    grouped,
    internal_source,
):

    for field in chunk.columns:
        chunk[field] = (
            chunk[field].map(clean)
        )

    for field in (
        "drugind_id",
        "record_id",
        "molregno",
        "max_phase_for_indication",
        "molecule_max_phase",
    ):
        chunk[field] = (
            chunk[field].map(phase)
        )

    # Preferred drug name, when supplied.
    original_name = chunk["preferred_name"]

    chunk["drug_name_lookup_status"] = (
        "preferred_name"
    )

    no_preferred = original_name.eq("")

    chunk.loc[
        no_preferred,
        "drug_name_lookup_status",
    ] = "identifier_placeholder"

    # An identifier fallback is not a resolved name.
    chunk["drug_name"] = (
        original_name.where(
            ~no_preferred,
            chunk["drug_id"],
        )
    )

    # Disease name selection and provenance.
    chunk["disease_name_lookup_status"] = (
        "efo_term"
    )

    disease = chunk["efo_term"].copy()

    empty = disease.eq("")

    dictionary = chunk[
        "dictionary_disease_name"
    ]

    use_dict = (
        empty
        & dictionary.ne("")
    )

    disease = disease.mask(
        use_dict,
        dictionary,
    )

    chunk.loc[
        use_dict,
        "disease_name_lookup_status",
    ] = "disease_dictionary"

    empty = disease.eq("")

    use_mesh = (
        empty
        & chunk["mesh_heading"].ne("")
    )

    disease = disease.mask(
        use_mesh,
        chunk["mesh_heading"],
    )

    chunk.loc[
        use_mesh,
        "disease_name_lookup_status",
    ] = "mesh_heading_fallback"

    empty = disease.eq("")

    use_id = (
        empty
        & chunk["disease_id"].ne("")
    )

    disease = disease.mask(
        use_id,
        chunk["disease_id"],
    )

    chunk.loc[
        use_id,
        "disease_name_lookup_status",
    ] = "identifier_placeholder"

    chunk["disease_name"] = disease

    # Unified identity fields.
    chunk["drug_identifier"] = (
        chunk["drug_id"]
    )

    chunk[
        "disease_or_condition_identifier"
    ] = chunk["disease_id"]

    chunk[
        "disease_or_condition_name"
    ] = chunk["disease_name"]

    chunk["relationship_type"] = (
        "drug_indication"
    )

    chunk["evidence_type"] = (
        "curated_drug_indication"
    )

    chunk["source"] = "ChEMBL"

    chunk["internal_source"] = (
        internal_source
    )

    # Validate namespaces separately from raw IDs.
    for entity, field, namespace in (
        (
            "drug",
            "drug_id",
            "ChEMBL",
        ),
        (
            "disease",
            "disease_id",
            "",
        ),
        (
            "mesh",
            "mesh_id",
            "MeSH",
        ),
    ):

        annotate(
            chunk,
            field,
            entity,
            namespace,
        )

    chunk["drug_identifier_type"] = (
        chunk[
            "drug_identifier_namespace"
        ]
    )

    chunk[
        "disease_or_condition_identifier_type"
    ] = chunk[
        "disease_identifier_namespace"
    ]

    # Add all indication references.
    types = []
    nums = []
    refs_json = []

    for row_id in chunk[
        "drugind_id"
    ].tolist():

        items = grouped.get(
            row_id,
            [],
        )

        types.append(
            "|".join(
                sorted(
                    {
                        clean(
                            x.get("ref_type")
                        )
                        for x in items
                        if clean(
                            x.get("ref_type")
                        )
                    }
                )
            )
        )

        nums.append(
            str(len(items))
        )

        refs_json.append(
            json.dumps(
                items,
                ensure_ascii=False,
                default=str,
            )
        )

    chunk[
        "indication_ref_types"
    ] = types

    chunk[
        "indication_ref_count"
    ] = nums

    chunk[
        "indication_refs_original_json"
    ] = refs_json

    # Preserve all eligible source indication events.
    eligible = (
        chunk["drug_name"].ne("")
        & chunk["disease_name"].ne("")
    )

    accepted = chunk.loc[
        eligible,
        COLS,
    ].copy()

    rejected = chunk.loc[
        ~eligible
    ].copy()

    return accepted, rejected


# ============================================================
# MAIN PROCESSING
# ============================================================

def main():

    p = argparse.ArgumentParser(
        description=__doc__
    )

    p.add_argument(
        "--workdir",
        default=str(RAW),
    )

    p.add_argument(
        "--out",
        default=str(OUT),
    )

    p.add_argument(
        "--metadata",
        default=str(META),
    )

    p.add_argument(
        "--input_archive",
        default=None,
    )

    p.add_argument(
        "--base_url",
        default=BASE_URL,
    )

    p.add_argument(
        "--download_missing",
        action="store_true",
    )

    p.add_argument(
        "--force_download",
        action="store_true",
    )

    p.add_argument(
        "--force_extract",
        action="store_true",
    )

    p.add_argument(
        "--keep_db",
        action="store_true",
        help=(
            "Accepted for compatibility; "
            "extracted DB is always cached"
        ),
    )

    p.add_argument(
        "--chunksize",
        type=int,
        default=50000,
    )

    p.add_argument(
        "--head_n",
        type=int,
        default=2,
    )

    args = p.parse_args()

    if args.chunksize < 1:
        p.error(
            "--chunksize must be positive"
        )

    workdir = Path(
        args.workdir
    ).expanduser().resolve()

    out = Path(
        args.out
    ).expanduser().resolve()

    metadata_path = Path(
        args.metadata
    ).expanduser().resolve()

    archive = (
        Path(
            args.input_archive
        ).expanduser().resolve()
        if args.input_archive
        else workdir / ARCHIVE_NAME
    )

    m = re.fullmatch(
        r"chembl_(\d+)_sqlite\.tar\.gz",
        archive.name,
    )

    if not m:
        raise ValueError(
            "Archive must be named "
            "chembl_XX_sqlite.tar.gz"
        )

    version = m[1]

    if version != "37":
        raise ValueError(
            "This reproducible run is "
            "locked to ChEMBL 37"
        )

    report = {
        "started_at_utc": now(),
        "source": "ChEMBL",
        "chembl_release": version,
        "status": "running",
        "output_fields": COLS,
    }

    temp_csv = out.with_name(
        out.name + ".processing.tmp"
    )

    ref_path = out.with_name(
        "chembl_indication_refs.csv"
    )

    ref_temp = None

    try:

        # ------------------------------------------------
        # STEP 1: SOURCE ARCHIVE
        # ------------------------------------------------

        section(
            "1. INPUT FILE INVENTORY — "
            "ChEMBL 37 LOCKED SNAPSHOT"
        )

        ensure_archive(
            archive,
            urljoin(
                args.base_url.rstrip("/")
                + "/",
                archive.name,
            ),
            args.download_missing
            or args.force_download,
            args.force_download,
        )

        report[
            "raw_file_inventory"
        ] = []

        for path in sorted(
            x for x in workdir.rglob("*")
            if x.is_file()
        ):

            report[
                "raw_file_inventory"
            ].append(
                {
                    "name": str(
                        path.relative_to(workdir)
                    ),
                    "bytes": (
                        path.stat().st_size
                    ),
                }
            )

            print(
                f"  {path.relative_to(workdir)}"
                f" | {path.stat().st_size:,} bytes"
            )

        report["archive"] = {
            "path": str(archive),
            "bytes": (
                archive.stat().st_size
            ),
            "sha256": sha256(archive),
        }

        # ------------------------------------------------
        # STEP 2: EXTRACT SQLITE
        # ------------------------------------------------

        db, member = extract_db(
            archive,
            workdir / "extracted",
            args.force_extract,
        )

        report["sqlite_member"] = member
        report["sqlite_path"] = str(db)

        conn = sqlite3.connect(
            f"file:{db}?mode=ro",
            uri=True,
        )

        try:

            # --------------------------------------------
            # STEP 3: VERIFY SOURCE SCHEMA
            # --------------------------------------------

            counts, headers = inspect_schema(
                conn,
                args.head_n,
            )

            report[
                "original_table_headers"
            ] = headers

            report[
                "original_table_rows"
            ] = counts

            sql, mapping = sql_builder(
                headers
            )

            report[
                "sql_join_design"
            ] = mapping

            print(
                "\n[SQL JOIN DESIGN]",
                mapping,
            )

            print(
                "FIELD MAPPING: "
                "molecule_dictionary.chembl_id "
                "-> drug_identifier"
            )

            print(
                "FIELD MAPPING: "
                "molecule_dictionary.pref_name "
                "-> drug_name"
            )

            print(
                "FIELD MAPPING: "
                "drug_indication.efo_id/efo_term "
                "-> disease ID/name"
            )

            print(
                "FIELD MAPPING: "
                "drug_indication.max_phase_for_ind "
                "-> max_phase_for_indication"
            )

            # --------------------------------------------
            # STEP 4: VERIFY JOIN
            # --------------------------------------------

            report["join_checks"] = (
                validate_joins(
                    conn,
                    sql,
                    counts["drug_indication"],
                )
            )

            # --------------------------------------------
            # STEP 5: PRESERVE REFERENCES
            # --------------------------------------------

            grouped, ref_rows, ref_temp = (
                references(
                    conn,
                    ref_path,
                )
            )

            report[
                "reference_rows"
            ] = ref_rows

            out.parent.mkdir(
                parents=True,
                exist_ok=True,
            )

            total = 0
            exported = 0
            excluded = 0

            reason_counts = Counter()
            name_counts = Counter()
            disease_counts = Counter()
            id_types = Counter()
            phase_counts = Counter()

            wrote = False

            section(
                "5. EXTRACTING SOURCE "
                "INDICATION EVENTS"
            )

            for chunk in pd.read_sql_query(
                sql,
                conn,
                chunksize=args.chunksize,
            ):

                total += len(chunk)

                transformed, invalid = (
                    transform(
                        chunk,
                        grouped,
                        f"ChEMBL{version}",
                    )
                )

                exported += len(transformed)
                excluded += len(invalid)

                for name in (
                    invalid["drug_id"]
                    .fillna("")
                    .astype(str)
                    .tolist()
                ):
                    reason_counts[
                        "missing_drug_id"
                    ] += (
                        not name.strip()
                    )

                for name in (
                    invalid["disease_name"]
                    .fillna("")
                    .astype(str)
                    .tolist()
                ):
                    reason_counts[
                        "missing_disease_name"
                    ] += (
                        not name.strip()
                    )

                name_counts.update(
                    transformed[
                        "drug_name_lookup_status"
                    ].value_counts().to_dict()
                )

                disease_counts.update(
                    transformed[
                        "disease_name_lookup_status"
                    ].value_counts().to_dict()
                )

                id_types.update(
                    transformed[
                        "disease_or_condition_identifier_type"
                    ].value_counts().to_dict()
                )

                phase_counts.update(
                    transformed[
                        "max_phase_for_indication"
                    ]
                    .replace("", "NULL")
                    .value_counts()
                    .to_dict()
                )

                transformed.to_csv(
                    temp_csv,
                    mode="a",
                    index=False,
                    header=not wrote,
                )

                wrote = True

                print(
                    f"  Read {total:,}"
                    f" -> eligible {exported:,}"
                    f" | excluded {excluded:,}",
                    flush=True,
                )

            # --------------------------------------------
            # STEP 6: VALIDATE OUTPUT
            # --------------------------------------------

            if not wrote:
                raise RuntimeError(
                    "No eligible ChEMBL "
                    "source records"
                )

            if (
                total
                != counts["drug_indication"]
                or total != exported + excluded
            ):
                raise RuntimeError(
                    "ChEMBL processing row "
                    "accounting mismatch"
                )

            header = pd.read_csv(
                temp_csv,
                nrows=0,
            ).columns.tolist()

            if header != COLS:
                raise RuntimeError(
                    "Output fields differ "
                    "from declared "
                    "merge-compatible schema"
                )

            preview = pd.read_csv(
                temp_csv,
                nrows=min(
                    3,
                    args.head_n,
                ),
                dtype=str,
            ).fillna("")

            if exported == 0:
                raise RuntimeError(
                    "Zero records written"
                )

            section(
                "6. VALIDATED OUTPUT ROWS, "
                "SCHEMA AND SAMPLE"
            )

            print(
                f"Original drug_indication: "
                f"{total:,}"
            )

            print(
                "Eligible merged-schema records: "
                f"{exported:,}"
            )

            print(
                "Excluded missing identity: "
                f"{excluded:,}"
            )

            print(
                f"Output fields ({len(COLS)}): "
                + " | ".join(COLS)
            )

            print(
                "[INDICATION PHASE]",
                dict(phase_counts),
            )

            print(
                "[DISEASE IDENTIFIER TYPES]",
                dict(id_types),
            )

            print(
                "[DRUG NAME RESOLUTION]",
                dict(name_counts),
            )

            print(
                "[DISEASE NAME RESOLUTION]",
                dict(disease_counts),
            )

            if len(preview):
                print(
                    "[ONE REAL OUTPUT ROW]"
                )

                for field, value in (
                    preview.iloc[0].items()
                ):
                    print(
                        f"  {field:38} "
                        f"{str(value)[:160]}"
                    )

            report.update(
                {
                    "rows_input": total,
                    "rows_output": exported,
                    "rows_excluded": excluded,

                    "exclusion_reasons":
                        dict(reason_counts),

                    "drug_name_resolution":
                        dict(name_counts),

                    "disease_name_resolution":
                        dict(disease_counts),

                    "disease_namespace_counts":
                        dict(id_types),

                    "indication_phase_counts":
                        dict(phase_counts),

                    "relationship_type":
                        "drug_indication",

                    "evidence_type":
                        "curated_drug_indication",

                    "row_semantics":
                        "Source indication evidence, "
                        "not universal regulatory approval",
                }
            )

        finally:
            conn.close()

        # ------------------------------------------------
        # STEP 7: OVERWRITE IN THE SAME DIRECTORY
        # ------------------------------------------------

        temp_csv.replace(out)

        if ref_temp:
            ref_temp.replace(
                ref_path
            )

        report["output_file"] = str(out)

        report[
            "output_sha256"
        ] = sha256(out)

        report[
            "references_csv"
        ] = (
            str(ref_path)
            if ref_rows
            else ""
        )

        report[
            "finished_at_utc"
        ] = now()

        report["status"] = "success"

        metadata_path.parent.mkdir(
            parents=True,
            exist_ok=True,
        )

        summary_tmp = (
            metadata_path.with_name(
                metadata_path.name + ".tmp"
            )
        )

        summary_tmp.write_text(
            json.dumps(
                report,
                indent=2,
                ensure_ascii=False,
            ),
            encoding="utf-8",
        )

        summary_tmp.replace(
            metadata_path
        )

        print(
            f"\nSUCCESS: {out}"
            f"\nSUMMARY: {metadata_path}"
        )

    except Exception as exc:

        failure = (
            metadata_path.with_name(
                "chembl_failed_run.json"
            )
        )

        failure.parent.mkdir(
            parents=True,
            exist_ok=True,
        )

        report.update(
            {
                "status": "failed",
                "failure": repr(exc),
                "finished_at_utc": now(),
            }
        )

        failure.write_text(
            json.dumps(
                report,
                indent=2,
                ensure_ascii=False,
            ),
            encoding="utf-8",
        )

        raise

    finally:
        temp_csv.unlink(
            missing_ok=True
        )

        if ref_temp is not None:
            ref_temp.unlink(
                missing_ok=True
            )


if __name__ == "__main__":
    main()
