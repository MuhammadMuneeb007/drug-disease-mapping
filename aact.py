#!/usr/bin/env python3
"""
AACT flat files -> unified drug-disease clinical trial co-occurrence CSV

Purpose
-------
Download/process the AACT flat-file export and create a publication-ready,
provenance-aware drug-disease mapping table using the same unified schema
as the other drug-disease source processors.

Source
------
AACT / ClinicalTrials.gov flat files

Input tables used
-----------------
interventions.txt
conditions.txt
studies.txt

Interpretation
--------------
AACT/ClinicalTrials.gov records represent clinical trial intervention-condition
co-occurrence. They do NOT imply:
    - drug approval
    - positive trial outcome
    - clinical efficacy
    - regulatory indication

Each output row means:
    this ClinicalTrials.gov study listed this drug intervention together
    with this disease/condition.

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

AACT-specific metadata columns
------------------------------
nct_id
intervention_type
study_phase
overall_status
brief_title
official_title
clinicaltrials_url
source_url

Default output paths
--------------------
/data/ascher02/uqmmune1/ANNOVAR/drug_disease_data/drug-disease-mapping/data/raw/aact_20260114
/data/ascher02/uqmmune1/ANNOVAR/drug_disease_data/drug-disease-mapping/data/processed/aact/aact_drug_disease.csv
/data/ascher02/uqmmune1/ANNOVAR/drug_disease_data/drug-disease-mapping/data/processed/aact/aact_drug_disease_deduplicated.csv
/data/ascher02/uqmmune1/ANNOVAR/drug_disease_data/drug-disease-mapping/data/processed/aact/aact_summary.json

Requirements
------------
pip install duckdb requests pandas
"""

from __future__ import annotations

import argparse
import hashlib
import json
import zipfile
from datetime import datetime, timezone
from pathlib import Path

import duckdb
import pandas as pd
import requests


# ============================================================
# CONFIG
# ============================================================

PROJECT_ROOT = Path(
    "/data/ascher02/uqmmune1/ANNOVAR/drug_disease_data/"
    "drug-disease-mapping"
)

RAW_DIR = PROJECT_ROOT / "data" / "raw" / "aact_20260114"
OUTPUT_DIR = PROJECT_ROOT / "data" / "processed" / "aact"

ZIP_FILE = RAW_DIR / "aact_export.zip"
EXTRACT_DIR = RAW_DIR / "extract"

OUTPUT_FILE = OUTPUT_DIR / "aact_drug_disease.csv"
DEDUP_FILE = OUTPUT_DIR / "aact_drug_disease_deduplicated.csv"
SUMMARY_FILE = OUTPUT_DIR / "aact_summary.json"

AACT_FLAT_ZIP_URL_20260114 = (
    "https://ctti-aact.nyc3.digitaloceanspaces.com/pw7s52pt0ighmd1qcb5hbl8hikcs"
)

REQUIRED_FILES = {
    "interventions.txt",
    "conditions.txt",
    "studies.txt",
}

REQUIRED_COLUMNS = {
    "interventions.txt": {"nct_id", "name", "intervention_type"},
    "conditions.txt": {"nct_id", "name"},
    "studies.txt": {
        "nct_id",
        "phase",
        "overall_status",
        "brief_title",
        "official_title",
    },
}

UNIFIED_COLUMNS = [
    # Unified schema
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

    # AACT-specific metadata
    "nct_id",
    "intervention_type",
    "study_phase",
    "overall_status",
    "brief_title",
    "official_title",
    "clinicaltrials_url",
    "source_url",
]


# ============================================================
# HELPERS
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


def clean_text(x) -> str:
    if pd.isna(x):
        return ""
    return str(x).strip()


def download_file(
    url: str,
    out_path: Path,
    force: bool = False,
    timeout: int = 1200,
) -> Path:
    out_path.parent.mkdir(parents=True, exist_ok=True)

    if out_path.exists() and out_path.stat().st_size > 0 and not force:
        print(f"[CACHE] Using existing file: {out_path}")
        print(f"        Size:   {human_size(out_path.stat().st_size)}")
        print(f"        SHA256: {sha256_file(out_path)}")

        try:
            with zipfile.ZipFile(out_path, "r") as z:
                bad = z.testzip()
                if bad is not None:
                    raise RuntimeError(f"ZIP integrity check failed at member: {bad}")

            print("[CACHE] ZIP integrity check: PASS")
            return out_path

        except zipfile.BadZipFile:
            print("[WARN] Existing file is not a valid ZIP. Re-downloading.")

        except Exception as e:
            print(f"[WARN] Existing ZIP failed validation: {e}. Re-downloading.")

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

    with zipfile.ZipFile(out_path, "r") as z:
        bad = z.testzip()
        if bad is not None:
            raise RuntimeError(f"Downloaded ZIP failed integrity check at member: {bad}")

    print("[DOWNLOADED] ZIP integrity check: PASS")

    return out_path


def inspect_zip(zip_path: Path, show_n: int = 30) -> list[str]:
    print_section("1. Inspecting AACT ZIP archive")

    with zipfile.ZipFile(zip_path, "r") as z:
        members = z.namelist()

    print(f"[ZIP] Number of files/folders inside archive: {len(members):,}")
    print(f"[ZIP] First {min(show_n, len(members))} members:")

    for m in members[:show_n]:
        print(f"  - {m}")

    basenames = {Path(m).name for m in members}
    missing = REQUIRED_FILES - basenames

    if missing:
        raise RuntimeError(f"Required files missing from AACT ZIP: {sorted(missing)}")

    print(f"[ZIP] Required files found: {sorted(REQUIRED_FILES)}")

    return members


def extract_selected(
    zip_path: Path,
    out_dir: Path,
    required: set[str],
    force: bool = False,
) -> dict[str, Path]:
    print_section("2. Extracting selected AACT files")

    out_dir.mkdir(parents=True, exist_ok=True)
    found: dict[str, Path] = {}

    with zipfile.ZipFile(zip_path, "r") as z:
        for member in z.namelist():
            base = Path(member).name

            if base not in required:
                continue

            dest = out_dir / base

            if dest.exists() and dest.stat().st_size > 0 and not force:
                print(f"[CACHE] Already extracted: {dest}")
                found[base] = dest
                continue

            print(f"[EXTRACT] {member} -> {dest}")

            temp_extract_dir = out_dir / "_tmp_extract"
            temp_extract_dir.mkdir(parents=True, exist_ok=True)

            z.extract(member, temp_extract_dir)

            extracted_matches = list(temp_extract_dir.rglob(base))

            if not extracted_matches:
                raise RuntimeError(f"Could not find extracted file for member: {member}")

            extracted_matches[0].replace(dest)

            # Clean temporary directory contents
            for p in sorted(temp_extract_dir.rglob("*"), reverse=True):
                try:
                    if p.is_file():
                        p.unlink()
                    elif p.is_dir():
                        p.rmdir()
                except Exception:
                    pass

            try:
                temp_extract_dir.rmdir()
            except Exception:
                pass

            found[base] = dest

    missing = required - set(found)

    if missing:
        raise RuntimeError(f"Missing required extracted files: {sorted(missing)}")

    for name, path in found.items():
        print(
            f"[FILE] {name}: {path} | "
            f"{human_size(path.stat().st_size)} | "
            f"SHA256={sha256_file(path)}"
        )

    return found


def inspect_pipe_file(
    path: Path,
    required_cols: set[str],
    n: int = 5,
) -> dict:
    print()
    print("-" * 100)
    print(f"Inspecting file: {path.name}")
    print("-" * 100)

    df_head = pd.read_csv(
        path,
        sep="|",
        nrows=n,
        dtype=str,
        engine="python",
        keep_default_na=False,
    )

    columns = list(df_head.columns)
    missing_cols = required_cols - set(columns)

    print(f"[COLUMNS] {len(columns)} columns")
    print(columns)
    print(f"[REQUIRED] {sorted(required_cols)}")

    if missing_cols:
        raise RuntimeError(
            f"{path.name} is missing required columns: {sorted(missing_cols)}"
        )

    print("[HEAD]")
    print(df_head.head(n).to_string(index=False))

    return {
        "file": path.name,
        "path": str(path),
        "size_bytes": path.stat().st_size,
        "sha256": sha256_file(path),
        "columns": columns,
        "required_columns": sorted(required_cols),
        "missing_required_columns": sorted(missing_cols),
        "head_preview": df_head.head(n).to_dict(orient="records"),
    }


def validate_output_schema(df: pd.DataFrame, required_columns: list[str]) -> None:
    missing = [c for c in required_columns if c not in df.columns]

    if missing:
        raise RuntimeError(f"Output is missing required columns: {missing}")


def write_csv_and_summary_stats(
    df: pd.DataFrame,
    output_file: Path,
    dedup_file: Path,
) -> dict:
    output_file.parent.mkdir(parents=True, exist_ok=True)

    df = df[UNIFIED_COLUMNS].copy()

    df.to_csv(output_file, index=False)

    dedup_cols = [
        "drug_name",
        "disease_or_condition_name",
        "relationship_type",
        "evidence_type",
        "source",
        "internal_source",
        "nct_id",
    ]

    dedup_df = df.drop_duplicates(subset=dedup_cols).copy()
    dedup_df.to_csv(dedup_file, index=False)

    stats = {
        "output_file": str(output_file),
        "output_rows": int(len(df)),
        "output_sha256": sha256_file(output_file),
        "deduplicated_file": str(dedup_file),
        "deduplicated_rows": int(len(dedup_df)),
        "deduplicated_sha256": sha256_file(dedup_file),
        "unique_trials": int(df["nct_id"].nunique()),
        "unique_drug_names": int(df["drug_name"].nunique()),
        "unique_disease_or_condition_names": int(
            df["disease_or_condition_name"].nunique()
        ),
        "relationship_type_counts": (
            df["relationship_type"]
            .value_counts(dropna=False)
            .to_dict()
        ),
        "evidence_type_counts": (
            df["evidence_type"]
            .value_counts(dropna=False)
            .to_dict()
        ),
        "overall_status_counts_top30": (
            df["overall_status"]
            .value_counts(dropna=False)
            .head(30)
            .to_dict()
        ),
        "study_phase_counts_top30": (
            df["study_phase"]
            .value_counts(dropna=False)
            .head(30)
            .to_dict()
        ),
        "output_head": df.head(10).to_dict(orient="records"),
    }

    return stats


# ============================================================
# CORE PROCESSING
# ============================================================

def run_duckdb_processing(
    interventions: Path,
    conditions: Path,
    studies: Path,
    output_file: Path,
    dedup_file: Path,
    keep_all_intervention_types: bool = False,
) -> dict:
    print_section("4. Joining AACT tables and creating unified output")

    con = duckdb.connect(database=":memory:")

    try:
        con.execute(
            f"""
            CREATE VIEW interventions AS
            SELECT *
            FROM read_csv_auto(
                '{interventions.as_posix()}',
                delim='|',
                header=True,
                quote='"',
                escape='"',
                ignore_errors=false
            );
            """
        )

        con.execute(
            f"""
            CREATE VIEW conditions AS
            SELECT *
            FROM read_csv_auto(
                '{conditions.as_posix()}',
                delim='|',
                header=True,
                quote='"',
                escape='"',
                ignore_errors=false
            );
            """
        )

        con.execute(
            f"""
            CREATE VIEW studies AS
            SELECT *
            FROM read_csv_auto(
                '{studies.as_posix()}',
                delim='|',
                header=True,
                quote='"',
                escape='"',
                ignore_errors=false
            );
            """
        )

        print("[ROW COUNTS BEFORE JOIN]")

        n_interventions = con.execute(
            "SELECT COUNT(*) FROM interventions"
        ).fetchone()[0]

        n_conditions = con.execute(
            "SELECT COUNT(*) FROM conditions"
        ).fetchone()[0]

        n_studies = con.execute(
            "SELECT COUNT(*) FROM studies"
        ).fetchone()[0]

        print(f"  interventions: {n_interventions:,}")
        print(f"  conditions:    {n_conditions:,}")
        print(f"  studies:       {n_studies:,}")

        print()
        print("[INTERVENTION TYPE COUNTS]")

        intervention_type_counts = con.execute(
            """
            SELECT
                COALESCE(intervention_type, '') AS intervention_type,
                COUNT(*) AS n
            FROM interventions
            GROUP BY COALESCE(intervention_type, '')
            ORDER BY n DESC
            LIMIT 50
            """
        ).fetchdf()

        print(intervention_type_counts.to_string(index=False))

        if keep_all_intervention_types:
            intervention_filter_sql = "TRUE"
            relationship_type_sql = """
                CASE
                    WHEN lower(i.intervention_type) = 'drug'
                        THEN 'clinical_trial_drug_condition'
                    ELSE 'clinical_trial_intervention_condition'
                END
            """
            evidence_type_sql = """
                CASE
                    WHEN lower(i.intervention_type) = 'drug'
                        THEN 'aact_clinical_trial_drug_condition_cooccurrence'
                    ELSE 'aact_clinical_trial_intervention_condition_cooccurrence'
                END
            """
        else:
            intervention_filter_sql = "lower(i.intervention_type) = 'drug'"
            relationship_type_sql = "'clinical_trial_drug_condition'"
            evidence_type_sql = "'aact_clinical_trial_drug_condition_cooccurrence'"

        query = f"""
            SELECT DISTINCT
                trim(i.name) AS drug_name,
                '' AS drug_identifier,
                '' AS drug_identifier_type,

                trim(c.name) AS disease_or_condition_name,
                '' AS disease_or_condition_identifier,
                '' AS disease_or_condition_identifier_type,

                {relationship_type_sql} AS relationship_type,
                {evidence_type_sql} AS evidence_type,

                'AACT' AS source,
                'ClinicalTrials.gov' AS internal_source,

                i.nct_id AS nct_id,
                COALESCE(i.intervention_type, '') AS intervention_type,
                COALESCE(s.phase, '') AS study_phase,
                COALESCE(s.overall_status, '') AS overall_status,
                COALESCE(s.brief_title, '') AS brief_title,
                COALESCE(s.official_title, '') AS official_title,

                concat('https://clinicaltrials.gov/study/', i.nct_id) AS clinicaltrials_url,
                '{AACT_FLAT_ZIP_URL_20260114}' AS source_url

            FROM interventions i
            JOIN conditions c
                USING (nct_id)
            LEFT JOIN studies s
                USING (nct_id)
            WHERE
                {intervention_filter_sql}
                AND i.name IS NOT NULL
                AND c.name IS NOT NULL
                AND trim(i.name) <> ''
                AND trim(c.name) <> ''
                AND i.nct_id IS NOT NULL
                AND trim(i.nct_id) <> ''
        """

        print()
        print("[OUTPUT PREVIEW BEFORE WRITING]")

        preview = con.execute(query + " LIMIT 10").fetchdf()
        print(preview.to_string(index=False))

        df = con.execute(query).fetchdf()

        for col in UNIFIED_COLUMNS:
            if col not in df.columns:
                df[col] = ""

        df = df[UNIFIED_COLUMNS].copy()

        for col in df.columns:
            df[col] = df[col].fillna("").astype(str).str.strip()

        validate_output_schema(df, UNIFIED_COLUMNS)

        print()
        print("[FINAL OUTPUT COUNTS]")
        print(f"  Rows before exact deduplication: {len(df):,}")
        print(f"  Unique trials:                 {df['nct_id'].nunique():,}")
        print(f"  Unique drug names:             {df['drug_name'].nunique():,}")
        print(
            f"  Unique disease/condition names: "
            f"{df['disease_or_condition_name'].nunique():,}"
        )

        stats = write_csv_and_summary_stats(
            df=df,
            output_file=output_file,
            dedup_file=dedup_file,
        )

        print()
        print("[OUTPUT SAVED]")
        print(f"  Full CSV:         {output_file}")
        print(f"  Full rows:        {stats['output_rows']:,}")
        print(f"  Full SHA256:      {stats['output_sha256']}")
        print(f"  Dedup CSV:        {dedup_file}")
        print(f"  Deduplicated rows:{stats['deduplicated_rows']:,}")
        print(f"  Dedup SHA256:     {stats['deduplicated_sha256']}")

        print()
        print("[OUTPUT HEAD]")
        print(df.head(10).to_string(index=False))

        return {
            "input_row_counts": {
                "interventions": int(n_interventions),
                "conditions": int(n_conditions),
                "studies": int(n_studies),
            },
            "intervention_type_counts_top50": intervention_type_counts.to_dict(
                orient="records"
            ),
            "keep_all_intervention_types": bool(keep_all_intervention_types),
            **stats,
        }

    finally:
        con.close()


# ============================================================
# MAIN
# ============================================================

def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Build a unified AACT drug-disease clinical trial co-occurrence table."
        )
    )

    parser.add_argument(
        "--url",
        default=AACT_FLAT_ZIP_URL_20260114,
        help="AACT flat-file ZIP URL.",
    )

    parser.add_argument(
        "--workdir",
        default=str(RAW_DIR),
        help="Raw working directory for AACT ZIP and extracted files.",
    )

    parser.add_argument(
        "--out",
        default=str(OUTPUT_FILE),
        help="Output full unified AACT CSV.",
    )

    parser.add_argument(
        "--dedup_out",
        default=str(DEDUP_FILE),
        help="Output deduplicated unified AACT CSV.",
    )

    parser.add_argument(
        "--summary",
        default=str(SUMMARY_FILE),
        help="Output summary JSON.",
    )

    parser.add_argument(
        "--force_download",
        action="store_true",
        help="Force re-download of AACT ZIP.",
    )

    parser.add_argument(
        "--force_extract",
        action="store_true",
        help="Force re-extraction of selected AACT files.",
    )

    parser.add_argument(
        "--keep_all_intervention_types",
        action="store_true",
        help=(
            "Keep all intervention types, not only DRUG. "
            "If used, drug_name should be interpreted as intervention_name."
        ),
    )

    parser.add_argument(
        "--head_n",
        type=int,
        default=5,
        help="Number of rows to preview from input files.",
    )

    args = parser.parse_args()

    workdir = Path(args.workdir).resolve()
    zip_path = workdir / "aact_export.zip"
    extract_dir = workdir / "extract"

    output_file = Path(args.out).resolve()
    dedup_file = Path(args.dedup_out).resolve()
    summary_file = Path(args.summary).resolve()

    metadata = {
        "script": "aact.py",
        "source": "AACT",
        "internal_source": "ClinicalTrials.gov",
        "source_url": args.url,
        "relationship_type_default": "clinical_trial_drug_condition",
        "evidence_type_default": "aact_clinical_trial_drug_condition_cooccurrence",
        "interpretation": (
            "AACT/ClinicalTrials.gov records represent clinical trial "
            "intervention-condition co-occurrence. They do not imply approval, "
            "efficacy, or positive trial outcome."
        ),
        "project_root": str(PROJECT_ROOT),
        "raw_dir": str(workdir),
        "extract_dir": str(extract_dir),
        "output_file": str(output_file),
        "deduplicated_output_file": str(dedup_file),
        "summary_file": str(summary_file),
        "unified_columns": UNIFIED_COLUMNS,
        "started_at_utc": now_utc(),
        "keep_all_intervention_types": bool(args.keep_all_intervention_types),
    }

    try:
        print_section("AACT -> unified drug-disease clinical trial mapping")

        print(f"[SOURCE URL] {args.url}")
        print(f"[WORKDIR]    {workdir}")
        print(f"[OUTPUT]     {output_file}")
        print(f"[DEDUP]      {dedup_file}")
        print(f"[SUMMARY]    {summary_file}")

        zip_path = download_file(
            url=args.url,
            out_path=zip_path,
            force=args.force_download,
        )

        metadata["downloaded_zip"] = {
            "path": str(zip_path),
            "size_bytes": zip_path.stat().st_size,
            "sha256": sha256_file(zip_path),
        }

        members = inspect_zip(zip_path)
        metadata["zip_member_count"] = len(members)
        metadata["zip_members_preview"] = members[:50]

        extracted = extract_selected(
            zip_path=zip_path,
            out_dir=extract_dir,
            required=REQUIRED_FILES,
            force=args.force_extract,
        )

        print_section("3. Inspecting selected AACT source files")

        metadata["selected_files"] = {}

        for filename in sorted(REQUIRED_FILES):
            metadata["selected_files"][filename] = inspect_pipe_file(
                path=extracted[filename],
                required_cols=REQUIRED_COLUMNS[filename],
                n=args.head_n,
            )

        processing_metadata = run_duckdb_processing(
            interventions=extracted["interventions.txt"],
            conditions=extracted["conditions.txt"],
            studies=extracted["studies.txt"],
            output_file=output_file,
            dedup_file=dedup_file,
            keep_all_intervention_types=args.keep_all_intervention_types,
        )

        metadata["processing"] = processing_metadata
        metadata["finished_at_utc"] = now_utc()
        metadata["status"] = "success"

    except Exception as e:
        metadata["finished_at_utc"] = now_utc()
        metadata["status"] = "failed"
        metadata["error"] = repr(e)
        raise

    finally:
        summary_file.parent.mkdir(parents=True, exist_ok=True)

        with open(summary_file, "w", encoding="utf-8") as f:
            json.dump(metadata, f, indent=2)

        print_section("SUMMARY WRITTEN")
        print(summary_file)


if __name__ == "__main__":
    main()