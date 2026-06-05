#!/usr/bin/env python3
"""
SIDER 4.1 -> unified drug-disease / drug-indication mapping CSV

Purpose
-------
Download and process SIDER 4.1 drug indication records into a
publication-ready, provenance-aware CSV using the same unified
multi-source drug-disease schema used by the other source processors.

Source
------
SIDER 4.1:
    http://sideeffects.embl.de/

Files used
----------
drug_names.tsv
meddra_all_indications.tsv.gz

Correct observed SIDER indication schema
----------------------------------------
meddra_all_indications.tsv.gz is parsed as a 7-column file:

    raw_col_0 = drug_id, e.g. CID100000085
    raw_col_1 = source_concept_id
    raw_col_2 = extraction_method, e.g. text_mention / NLP_indication
    raw_col_3 = source_concept_name
    raw_col_4 = meddra_concept_type, e.g. LLT / PT
    raw_col_5 = disease_id, e.g. C0015544
    raw_col_6 = disease_name, e.g. Failure to Thrive

Interpretation
--------------
SIDER indications are label-derived drug-indication / phenotype records.
They should not be interpreted as direct evidence of clinical efficacy or
regulatory approval without additional validation.

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

SIDER-specific metadata columns
-------------------------------
sider_drug_id
source_concept_id
source_concept_name
meddra_concept_type
extraction_method
relationship
raw_col_0
raw_col_1
raw_col_2
raw_col_3
raw_col_4
raw_col_5
raw_col_6

Default output
--------------
/data/ascher02/uqmmune1/ANNOVAR/drug_disease_data/drug-disease-mapping/data/processed/sider/sider_drug_disease.csv
/data/ascher02/uqmmune1/ANNOVAR/drug_disease_data/drug-disease-mapping/data/processed/sider/sider_drug_disease_deduplicated.csv
/data/ascher02/uqmmune1/ANNOVAR/drug_disease_data/drug-disease-mapping/data/processed/sider/sider_summary.json

Requirements
------------
pip install pandas requests
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Tuple

import pandas as pd
import requests


# ============================================================
# CONFIG
# ============================================================

PROJECT_ROOT = Path(
    "/data/ascher02/uqmmune1/ANNOVAR/drug_disease_data/"
    "drug-disease-mapping"
)

SOURCE_NAME = "SIDER4.1"

BASE_URL_CANDIDATES = [
    "https://sideeffects.embl.de/media/download",
    "http://sideeffects.embl.de/media/download",
]

FILES = {
    "drug_names": "drug_names.tsv",
    "indications": "meddra_all_indications.tsv.gz",
}

RAW_DIR = PROJECT_ROOT / "data" / "raw" / "sider"
OUTPUT_DIR = PROJECT_ROOT / "data" / "processed" / "sider"

OUTPUT_FILE = OUTPUT_DIR / "sider_drug_disease.csv"
DEDUP_FILE = OUTPUT_DIR / "sider_drug_disease_deduplicated.csv"
SUMMARY_FILE = OUTPUT_DIR / "sider_summary.json"

UA = "Mozilla/5.0 (SIDERDrugDiseaseExporter/5.0)"


OUTPUT_COLUMNS = [
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

    # SIDER-specific metadata
    "sider_drug_id",
    "source_concept_id",
    "source_concept_name",
    "meddra_concept_type",
    "extraction_method",
    "relationship",
    "evidence_text",

    # Raw retained columns for traceability
    "raw_col_0",
    "raw_col_1",
    "raw_col_2",
    "raw_col_3",
    "raw_col_4",
    "raw_col_5",
    "raw_col_6",
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


def clean_cell(x: Any) -> str:
    if x is None:
        return ""

    try:
        if pd.isna(x):
            return ""
    except Exception:
        pass

    return str(x).strip()


def sha256_file(path: Path, chunk_size: int = 1024 * 1024) -> str:
    h = hashlib.sha256()

    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(chunk_size), b""):
            h.update(chunk)

    return h.hexdigest()


def human_size(num_bytes: int) -> str:
    size = float(num_bytes)

    for unit in ["B", "KB", "MB", "GB", "TB"]:
        if size < 1024:
            return f"{size:.2f} {unit}"
        size /= 1024

    return f"{size:.2f} PB"


def count_missing_or_empty(df: pd.DataFrame) -> pd.DataFrame:
    rows = []

    for col in df.columns:
        s = df[col]
        missing_count = int(s.isna().sum())
        empty_count = int(s.fillna("").astype(str).str.strip().eq("").sum())

        rows.append(
            {
                "column": col,
                "missing_count": missing_count,
                "missing_percent": round(
                    missing_count / max(len(df), 1) * 100,
                    3,
                ),
                "empty_string_count": empty_count,
                "empty_string_percent": round(
                    empty_count / max(len(df), 1) * 100,
                    3,
                ),
            }
        )

    return pd.DataFrame(rows)


def download_one_url(url: str, out_path: Path, timeout: int = 300) -> Path:
    out_path.parent.mkdir(parents=True, exist_ok=True)

    print(f"  [TRY] {url}")

    with requests.get(
        url,
        headers={"User-Agent": UA},
        stream=True,
        timeout=timeout,
    ) as r:
        r.raise_for_status()

        tmp = out_path.with_suffix(out_path.suffix + ".part")

        with open(tmp, "wb") as f:
            for chunk in r.iter_content(chunk_size=1024 * 1024):
                if chunk:
                    f.write(chunk)

        tmp.replace(out_path)

    return out_path


def download_with_fallback(
    filename: str,
    out_path: Path,
    timeout: int = 300,
    force: bool = False,
) -> Tuple[Path, str]:
    if out_path.exists() and out_path.stat().st_size > 0 and not force:
        print(f"  [CACHED] {out_path}")
        print(f"           Size:   {human_size(out_path.stat().st_size)}")
        print(f"           SHA256: {sha256_file(out_path)}")
        return out_path, "cached"

    if out_path.exists() and force:
        out_path.unlink()

    last_error = None

    for base_url in BASE_URL_CANDIDATES:
        url = f"{base_url.rstrip('/')}/{filename}"

        try:
            path = download_one_url(url, out_path, timeout=timeout)

            print(f"  [DOWNLOADED] {path}")
            print(f"               Size:   {human_size(path.stat().st_size)}")
            print(f"               SHA256: {sha256_file(path)}")

            return path, url

        except Exception as e:
            last_error = e
            print(f"  [FAILED] {url}")
            print(f"           {repr(e)}")

    raise RuntimeError(
        f"Could not download {filename} from SIDER endpoints. "
        f"Last error: {repr(last_error)}"
    )


def read_tsv(path: Path, gz: bool = False) -> pd.DataFrame:
    return pd.read_csv(
        path,
        sep="\t",
        header=None,
        compression="gzip" if gz else None,
        dtype=str,
        low_memory=False,
        keep_default_na=False,
    )


def infer_sider_drug_identifier_type(drug_id: str) -> str:
    drug_id = clean_cell(drug_id)

    if not drug_id:
        return ""

    if drug_id.startswith("CID"):
        return "SIDER_STITCH_CID"

    return "SIDER_DrugID"


def infer_disease_identifier_type(disease_id: str) -> str:
    disease_id = clean_cell(disease_id)

    if not disease_id:
        return ""

    upper = disease_id.upper()

    # SIDER MedDRA indications commonly use UMLS-style CUI disease IDs.
    if upper.startswith("C") and upper[1:].isdigit():
        return "UMLS"

    if upper.startswith("UMLS:"):
        return "UMLS"

    if upper.startswith("MESH:"):
        return "MeSH"

    if upper.startswith("MEDDRA:"):
        return "MedDRA"

    return "SIDER_DiseaseID"


# ============================================================
# LOADERS
# ============================================================

def load_drug_names(path: Path) -> Tuple[pd.DataFrame, Dict[str, str]]:
    print_section("Loading SIDER drug_names.tsv")
    print(f"[FILE] {path}")

    df = read_tsv(path, gz=False)

    if df.shape[1] < 2:
        raise RuntimeError(
            f"drug_names.tsv has unexpected number of columns: {df.shape[1]}"
        )

    df = df.iloc[:, :2].copy()
    df.columns = ["drug_id", "drug_name"]

    df["drug_id"] = df["drug_id"].map(clean_cell)
    df["drug_name"] = df["drug_name"].map(clean_cell)

    df = df[
        (df["drug_id"] != "")
        & (df["drug_name"] != "")
    ].copy()

    drug_name_map = (
        df.drop_duplicates(subset=["drug_id"], keep="first")
        .set_index("drug_id")["drug_name"]
        .to_dict()
    )

    print(f"[ROWS] {len(df):,}")
    print(f"[UNIQUE DRUG IDS] {df['drug_id'].nunique():,}")
    print(f"[PRIMARY NAME MAP] {len(drug_name_map):,}")

    print("\n[HEAD]")
    print(df.head(10).to_string(index=False))

    return df, drug_name_map


def load_indications(path: Path) -> pd.DataFrame:
    print_section("Loading SIDER meddra_all_indications.tsv.gz")
    print(f"[FILE] {path}")

    raw = read_tsv(path, gz=True)

    print(f"[RAW ROWS] {len(raw):,}")
    print(f"[RAW COLUMNS] {raw.shape[1]:,}")

    if raw.shape[1] < 7:
        raise RuntimeError(
            f"Expected at least 7 columns in meddra_all_indications.tsv.gz, "
            f"but found {raw.shape[1]}."
        )

    df = raw.copy()
    df.columns = [f"raw_col_{i}" for i in range(raw.shape[1])]

    # Correct SIDER 7-column schema.
    df["drug_id"] = df["raw_col_0"].map(clean_cell)
    df["source_concept_id"] = df["raw_col_1"].map(clean_cell)
    df["extraction_method"] = df["raw_col_2"].map(clean_cell)
    df["source_concept_name"] = df["raw_col_3"].map(clean_cell)
    df["meddra_concept_type"] = df["raw_col_4"].map(clean_cell)
    df["disease_id"] = df["raw_col_5"].map(clean_cell)
    df["disease_name"] = df["raw_col_6"].map(clean_cell)

    parsed_cols = [
        "drug_id",
        "source_concept_id",
        "extraction_method",
        "source_concept_name",
        "meddra_concept_type",
        "disease_id",
        "disease_name",
    ]

    print("[PARSED COLUMNS]")
    for col in parsed_cols:
        print(f"  - {col}")

    print("\n[HEAD]")
    print(df[parsed_cols].head(10).to_string(index=False))

    return df


# ============================================================
# PROCESSING
# ============================================================

def build_output(
    indications: pd.DataFrame,
    drug_name_map: Dict[str, str],
) -> pd.DataFrame:
    print_section("Building unified SIDER drug-disease output")

    df = indications.copy()

    out = pd.DataFrame()

    sider_drug_id = df["drug_id"].map(clean_cell)
    disease_id = df["disease_id"].map(clean_cell)
    disease_name = df["disease_name"].map(clean_cell)

    out["drug_name"] = sider_drug_id.map(lambda x: drug_name_map.get(x, ""))
    out["drug_identifier"] = sider_drug_id
    out["drug_identifier_type"] = sider_drug_id.map(infer_sider_drug_identifier_type)

    out["disease_or_condition_name"] = disease_name
    out["disease_or_condition_identifier"] = disease_id
    out["disease_or_condition_identifier_type"] = disease_id.map(
        infer_disease_identifier_type
    )

    out["relationship_type"] = "indication"
    out["evidence_type"] = "sider_label_derived_indication"

    out["source"] = SOURCE_NAME
    out["internal_source"] = "SIDER_meddra_all_indications"

    # SIDER-specific metadata.
    out["sider_drug_id"] = sider_drug_id
    out["source_concept_id"] = df["source_concept_id"].map(clean_cell)
    out["source_concept_name"] = df["source_concept_name"].map(clean_cell)
    out["meddra_concept_type"] = df["meddra_concept_type"].map(clean_cell)
    out["extraction_method"] = df["extraction_method"].map(clean_cell)
    out["relationship"] = "indication"

    out["evidence_text"] = (
        "source_concept_id="
        + out["source_concept_id"].astype(str)
        + "; source_concept_name="
        + out["source_concept_name"].astype(str)
        + "; meddra_concept_type="
        + out["meddra_concept_type"].astype(str)
        + "; extraction_method="
        + out["extraction_method"].astype(str)
    )

    # Retain raw columns for traceability.
    for i in range(7):
        raw_col = f"raw_col_{i}"

        if raw_col in df.columns:
            out[raw_col] = df[raw_col].map(clean_cell)
        else:
            out[raw_col] = ""

    before = len(out)

    out = out[
        (out["drug_identifier"] != "")
        & (out["disease_or_condition_identifier"] != "")
        & (out["disease_or_condition_name"] != "")
    ].copy()

    removed_empty = before - len(out)

    before_dup = len(out)

    out = out.drop_duplicates(
        subset=[
            "drug_identifier",
            "disease_or_condition_identifier",
            "disease_or_condition_name",
            "relationship_type",
            "meddra_concept_type",
            "extraction_method",
        ]
    ).reset_index(drop=True)

    removed_duplicates = before_dup - len(out)

    out = out[OUTPUT_COLUMNS]

    print(f"[INPUT ROWS] {len(df):,}")
    print(f"[OUTPUT ROWS] {len(out):,}")
    print(f"[REMOVED EMPTY CORE RECORDS] {removed_empty:,}")
    print(f"[REMOVED DUPLICATES] {removed_duplicates:,}")
    print(f"[UNIQUE DRUG IDS] {out['drug_identifier'].nunique():,}")
    print(
        f"[UNIQUE DRUG NAMES] "
        f"{out['drug_name'].replace('', pd.NA).nunique(dropna=True):,}"
    )
    print(
        f"[UNIQUE DISEASE IDS] "
        f"{out['disease_or_condition_identifier'].nunique():,}"
    )
    print(
        f"[UNIQUE DISEASE NAMES] "
        f"{out['disease_or_condition_name'].nunique():,}"
    )
    print(
        f"[UNIQUE DRUG-DISEASE PAIRS] "
        f"{out[['drug_identifier', 'disease_or_condition_identifier']].drop_duplicates().shape[0]:,}"
    )

    # ------------------------------------------------------------
    # Sanity checks for the known SIDER parsing bug.
    # ------------------------------------------------------------
    bad_disease_ids = int(
        out["disease_or_condition_identifier"].isin(["LLT", "PT"]).sum()
    )

    bad_drug_ids = int(
        out["drug_identifier"].str.match(r"^C\d+", na=False).sum()
    )

    print("\n[SANITY CHECKS]")
    print(f"[BAD disease identifier is LLT/PT] {bad_disease_ids:,}")
    print(f"[BAD drug identifier looks like UMLS C-code] {bad_drug_ids:,}")

    if bad_disease_ids > 0:
        raise RuntimeError(
            "Parsing error: disease identifier contains LLT/PT. "
            "The MedDRA concept type column was incorrectly used as disease ID."
        )

    if bad_drug_ids > 0:
        raise RuntimeError(
            "Parsing error: drug identifier looks like a UMLS disease ID. "
            "The wrong column was used as drug ID."
        )

    print("\n[OUTPUT HEAD]")
    preview_cols = [
        "drug_name",
        "drug_identifier",
        "drug_identifier_type",
        "disease_or_condition_name",
        "disease_or_condition_identifier",
        "disease_or_condition_identifier_type",
        "relationship_type",
        "evidence_type",
        "meddra_concept_type",
        "extraction_method",
        "source",
    ]

    print(out[preview_cols].head(20).to_string(index=False))

    return out


def write_csv(df: pd.DataFrame, out_csv: Path) -> None:
    print_section("Writing SIDER CSV")

    out_csv.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_csv, index=False, quoting=csv.QUOTE_ALL)

    print(f"[OUTPUT] {out_csv}")
    print(f"[ROWS] {len(df):,}")
    print(f"[SIZE] {human_size(out_csv.stat().st_size)}")
    print(f"[SHA256] {sha256_file(out_csv)}")


def write_deduplicated_csv(df: pd.DataFrame, dedup_csv: Path) -> pd.DataFrame:
    print_section("Writing deduplicated SIDER CSV")

    dedup_cols = [
        "drug_name",
        "drug_identifier",
        "disease_or_condition_name",
        "disease_or_condition_identifier",
        "relationship_type",
        "evidence_type",
        "meddra_concept_type",
        "extraction_method",
    ]

    dedup = df.drop_duplicates(subset=dedup_cols).copy()

    dedup_csv.parent.mkdir(parents=True, exist_ok=True)
    dedup.to_csv(dedup_csv, index=False, quoting=csv.QUOTE_ALL)

    print(f"[DEDUP OUTPUT] {dedup_csv}")
    print(f"[DEDUP ROWS] {len(dedup):,}")
    print(f"[DEDUP SIZE] {human_size(dedup_csv.stat().st_size)}")
    print(f"[DEDUP SHA256] {sha256_file(dedup_csv)}")

    return dedup


def summarise_output(out: pd.DataFrame) -> Dict[str, Any]:
    print_section("Output validation and summary")

    missing = count_missing_or_empty(out)

    print("[MISSING / EMPTY SUMMARY]")
    print(missing.to_string(index=False))

    relationship_counts = (
        out["relationship_type"]
        .replace("", "Missing")
        .value_counts()
        .reset_index()
    )
    relationship_counts.columns = ["relationship_type", "n"]

    evidence_type_counts = (
        out["evidence_type"]
        .replace("", "Missing")
        .value_counts()
        .reset_index()
    )
    evidence_type_counts.columns = ["evidence_type", "n"]

    meddra_type_counts = (
        out["meddra_concept_type"]
        .replace("", "Missing")
        .value_counts()
        .reset_index()
    )
    meddra_type_counts.columns = ["meddra_concept_type", "n"]

    extraction_counts = (
        out["extraction_method"]
        .replace("", "Missing")
        .value_counts()
        .reset_index()
    )
    extraction_counts.columns = ["extraction_method", "n"]

    drug_identifier_type_counts = (
        out["drug_identifier_type"]
        .replace("", "Missing")
        .value_counts()
        .reset_index()
    )
    drug_identifier_type_counts.columns = ["drug_identifier_type", "n"]

    disease_identifier_type_counts = (
        out["disease_or_condition_identifier_type"]
        .replace("", "Missing")
        .value_counts()
        .reset_index()
    )
    disease_identifier_type_counts.columns = [
        "disease_or_condition_identifier_type",
        "n",
    ]

    print("\n[RELATIONSHIP TYPE COUNTS]")
    print(relationship_counts.to_string(index=False))

    print("\n[EVIDENCE TYPE COUNTS]")
    print(evidence_type_counts.to_string(index=False))

    print("\n[MEDDRA CONCEPT TYPE COUNTS]")
    print(meddra_type_counts.to_string(index=False))

    print("\n[EXTRACTION METHOD COUNTS]")
    print(extraction_counts.to_string(index=False))

    summary = {
        "rows": int(len(out)),
        "columns": list(out.columns),
        "exact_duplicate_rows": int(out.duplicated().sum()),
        "unique_drug_identifiers": int(
            out["drug_identifier"].replace("", pd.NA).nunique(dropna=True)
        ),
        "unique_drug_names": int(
            out["drug_name"].replace("", pd.NA).nunique(dropna=True)
        ),
        "unique_disease_identifiers": int(
            out["disease_or_condition_identifier"]
            .replace("", pd.NA)
            .nunique(dropna=True)
        ),
        "unique_disease_names": int(
            out["disease_or_condition_name"]
            .replace("", pd.NA)
            .nunique(dropna=True)
        ),
        "unique_drug_disease_pairs": int(
            out[
                [
                    "drug_identifier",
                    "disease_or_condition_identifier",
                ]
            ]
            .drop_duplicates()
            .shape[0]
        ),
        "missing_or_empty_summary": missing.to_dict(orient="records"),
        "relationship_type_counts": relationship_counts.to_dict(orient="records"),
        "evidence_type_counts": evidence_type_counts.to_dict(orient="records"),
        "meddra_concept_type_counts": meddra_type_counts.to_dict(orient="records"),
        "extraction_method_counts": extraction_counts.to_dict(orient="records"),
        "drug_identifier_type_counts": drug_identifier_type_counts.to_dict(
            orient="records"
        ),
        "disease_identifier_type_counts": disease_identifier_type_counts.to_dict(
            orient="records"
        ),
        "output_head": out.head(20).to_dict(orient="records"),
    }

    print("\n[SUMMARY]")
    for key in [
        "rows",
        "exact_duplicate_rows",
        "unique_drug_identifiers",
        "unique_drug_names",
        "unique_disease_identifiers",
        "unique_disease_names",
        "unique_drug_disease_pairs",
    ]:
        print(f"{key}: {summary[key]}")

    return summary


def write_metadata(metadata_path: Path, metadata: Dict[str, Any]) -> None:
    metadata_path.parent.mkdir(parents=True, exist_ok=True)

    with open(metadata_path, "w", encoding="utf-8") as f:
        json.dump(metadata, f, indent=2, default=str)

    print_section("METADATA WRITTEN")
    print(metadata_path)


# ============================================================
# MAIN
# ============================================================

def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Download and process SIDER 4.1 drug-indication mappings "
            "into the unified multi-source drug-disease schema."
        )
    )

    parser.add_argument(
        "--workdir",
        default=str(RAW_DIR),
        help="Directory for downloaded SIDER files.",
    )

    parser.add_argument(
        "--out",
        default=str(OUTPUT_FILE),
        help="Output unified CSV path.",
    )

    parser.add_argument(
        "--dedup_out",
        default=str(DEDUP_FILE),
        help="Output deduplicated unified CSV path.",
    )

    parser.add_argument(
        "--metadata",
        default=str(SUMMARY_FILE),
        help="Output metadata JSON path.",
    )

    parser.add_argument(
        "--force_download",
        action="store_true",
        help="Force re-download of SIDER files.",
    )

    args = parser.parse_args()

    workdir = Path(args.workdir).resolve()
    out_csv = Path(args.out).resolve()
    dedup_csv = Path(args.dedup_out).resolve()
    metadata_path = Path(args.metadata).resolve()

    metadata: Dict[str, Any] = {
        "script": "sider.py",
        "source": SOURCE_NAME,
        "started_at_utc": now_utc(),
        "status": "started",
        "project_root": str(PROJECT_ROOT),
        "raw_dir": str(workdir),
        "output_file": str(out_csv),
        "deduplicated_output_file": str(dedup_csv),
        "metadata_file": str(metadata_path),
        "download_base_url_candidates": BASE_URL_CANDIDATES,
        "files": FILES,
        "relationship_type": "indication",
        "evidence_type": "sider_label_derived_indication",
        "interpretation": (
            "SIDER 4.1 is derived from marketed drug labels/package inserts. "
            "The meddra_all_indications file is interpreted as label-derived "
            "drug-indication/phenotype evidence, not direct therapeutic efficacy "
            "or proof of regulatory approval."
        ),
        "correct_schema": {
            "raw_col_0": "drug_id",
            "raw_col_1": "source_concept_id",
            "raw_col_2": "extraction_method",
            "raw_col_3": "source_concept_name",
            "raw_col_4": "meddra_concept_type",
            "raw_col_5": "disease_id",
            "raw_col_6": "disease_name",
        },
        "unified_output_columns": OUTPUT_COLUMNS,
    }

    try:
        print_section("SIDER 4.1 -> unified drug-disease / indication CSV")
        print(f"[WORKDIR]  {workdir}")
        print(f"[OUTPUT]   {out_csv}")
        print(f"[DEDUP]    {dedup_csv}")
        print(f"[METADATA] {metadata_path}")

        print_section("Downloading SIDER files")

        drug_names_path, drug_names_url = download_with_fallback(
            FILES["drug_names"],
            workdir / FILES["drug_names"],
            force=args.force_download,
        )

        indications_path, indications_url = download_with_fallback(
            FILES["indications"],
            workdir / FILES["indications"],
            force=args.force_download,
        )

        metadata["downloads"] = {
            "drug_names": {
                "url": drug_names_url,
                "path": str(drug_names_path),
                "size_bytes": drug_names_path.stat().st_size,
                "sha256": sha256_file(drug_names_path),
            },
            "indications": {
                "url": indications_url,
                "path": str(indications_path),
                "size_bytes": indications_path.stat().st_size,
                "sha256": sha256_file(indications_path),
            },
        }

        drug_names_df, drug_name_map = load_drug_names(drug_names_path)
        indications_df = load_indications(indications_path)

        metadata["input_summary"] = {
            "drug_names_rows": int(len(drug_names_df)),
            "drug_names_unique_drug_ids": int(
                drug_names_df["drug_id"].nunique()
            ),
            "indications_rows": int(len(indications_df)),
            "indications_columns": list(indications_df.columns),
        }

        out = build_output(indications_df, drug_name_map)

        if out.empty:
            raise RuntimeError("SIDER output is empty after processing.")

        write_csv(out, out_csv)
        dedup = write_deduplicated_csv(out, dedup_csv)

        metadata["output"] = {
            "path": str(out_csv),
            "rows": int(len(out)),
            "size_bytes": out_csv.stat().st_size,
            "sha256": sha256_file(out_csv),
        }

        metadata["deduplicated_output"] = {
            "path": str(dedup_csv),
            "rows": int(len(dedup)),
            "size_bytes": dedup_csv.stat().st_size,
            "sha256": sha256_file(dedup_csv),
        }

        metadata["output_summary"] = summarise_output(out)

        metadata["finished_at_utc"] = now_utc()
        metadata["status"] = "success"

        print_section("DONE")
        print("[SUCCESS] Unified SIDER drug-disease / indication mapping generated.")
        print(f"[CSV]      {out_csv}")
        print(f"[DEDUP]    {dedup_csv}")
        print(f"[METADATA] {metadata_path}")

    except Exception as e:
        metadata["finished_at_utc"] = now_utc()
        metadata["status"] = "failed"
        metadata["error"] = repr(e)

        print_section("FAILED")
        print(repr(e))

        raise

    finally:
        write_metadata(metadata_path, metadata)


if __name__ == "__main__":
    main()