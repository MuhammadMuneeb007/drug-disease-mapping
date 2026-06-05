#!/usr/bin/env python3
"""
Open Targets Platform -> unified drug-disease mapping CSV

Purpose
-------
Download Open Targets Platform clinical_indication, drug_molecule,
and disease parquet datasets, then export a publication-ready,
provenance-aware drug-disease mapping CSV using the unified schema.

Main dataset used
-----------------
clinical_indication

Supporting datasets
-------------------
drug_molecule
disease

Default Open Targets release
----------------------------
26.03

Input base URL
--------------
https://ftp.ebi.ac.uk/pub/databases/opentargets/platform

Example release URLs
--------------------
https://ftp.ebi.ac.uk/pub/databases/opentargets/platform/26.03/output/clinical_indication/
https://ftp.ebi.ac.uk/pub/databases/opentargets/platform/26.03/output/drug_molecule/
https://ftp.ebi.ac.uk/pub/databases/opentargets/platform/26.03/output/disease/

Interpretation
--------------
Open Targets clinical_indication rows represent drug-disease clinical
indication evidence with clinical stage metadata. These records should not
all be interpreted as approved indications; use max_phase / clinical_phase_raw
to distinguish clinical stage.

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

Open Targets-specific metadata columns
--------------------------------------
drug_id
disease_id
disease_name
max_phase
clinical_phase_raw
clinical_report_ids
opentargets_indication_id
dataset
release

Requirements
------------
pip install pandas pyarrow requests
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import pandas as pd
import pyarrow.parquet as pq
import requests


# =============================================================================
# Configuration
# =============================================================================

BASE = "https://ftp.ebi.ac.uk/pub/databases/opentargets/platform"
DEFAULT_RELEASE = "26.03"
UA = "Mozilla/5.0 (OpenTargetsDrugDiseaseExporter/2.0)"


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

    # Open Targets-specific metadata
    "drug_id",
    "disease_id",
    "disease_name",
    "max_phase",
    "clinical_phase_raw",
    "clinical_report_ids",
    "opentargets_indication_id",
    "dataset",
    "release",
]


# =============================================================================
# Helper functions
# =============================================================================

def now_utc() -> str:
    return datetime.now(timezone.utc).isoformat()


def print_section(title: str) -> None:
    print("\n" + "=" * 100)
    print(title)
    print("=" * 100)


def clean_cell(x: Any) -> str:
    """
    Clean scalar values, lists, numpy arrays, and pandas/pyarrow-like arrays.

    Open Targets clinicalReportIds can appear as an array. This converts it
    into a clean semicolon-separated string.
    """
    if x is None:
        return ""

    try:
        if pd.isna(x):
            return ""
    except Exception:
        pass

    if isinstance(x, (list, tuple, set)) or hasattr(x, "tolist"):
        try:
            values = x.tolist() if hasattr(x, "tolist") else list(x)
            cleaned = [clean_cell(v) for v in values]
            cleaned = [v for v in cleaned if v]
            return ";".join(sorted(set(cleaned)))
        except Exception:
            pass

    if isinstance(x, dict):
        return json.dumps(x, ensure_ascii=False, sort_keys=True)

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


def safe_phase_to_int(x: Any) -> int:
    """
    Convert Open Targets maxClinicalStage into integer max_phase.

    Examples:
        APPROVAL -> 4
        Phase 3  -> 3
        2        -> 2
        Unknown  -> 0
    """
    s = clean_cell(x)

    if not s:
        return 0

    low = s.lower()

    if low in {"approval", "approved", "launched", "marketed"}:
        return 4

    nums = re.findall(r"\d+(?:\.\d+)?", low)

    if nums:
        try:
            return int(float(nums[-1]))
        except Exception:
            return 0

    return 0


def infer_drug_identifier_type(drug_id: str) -> str:
    """
    Open Targets drug IDs are commonly ChEMBL IDs for small molecules/drugs.
    """
    x = clean_cell(drug_id)

    if not x:
        return ""

    if x.upper().startswith("CHEMBL"):
        return "ChEMBL"

    return "OpenTargets_DrugID"


def infer_disease_identifier_type(disease_id: str) -> str:
    """
    Open Targets disease IDs are usually EFO, MONDO, HP, Orphanet, etc.
    """
    x = clean_cell(disease_id)

    if not x:
        return ""

    upper = x.upper()

    if upper.startswith("EFO_") or upper.startswith("EFO:"):
        return "EFO"

    if upper.startswith("MONDO_") or upper.startswith("MONDO:"):
        return "MONDO"

    if upper.startswith("HP_") or upper.startswith("HP:"):
        return "HPO"

    if upper.startswith("ORPHANET_") or upper.startswith("ORPHANET:"):
        return "Orphanet"

    if upper.startswith("DOID_") or upper.startswith("DOID:"):
        return "DOID"

    return "OpenTargets_DiseaseID"


def fetch_text(url: str, timeout: int = 90) -> str:
    r = requests.get(url, headers={"User-Agent": UA}, timeout=timeout)
    r.raise_for_status()
    return r.text


def list_parquet_files(dir_url: str) -> List[str]:
    html = fetch_text(dir_url)
    files = re.findall(r'href="([^"]+\.parquet)"', html)
    files = sorted(set(files))
    return [dir_url.rstrip("/") + "/" + f for f in files]


def download(url: str, out_path: Path, timeout: int = 600) -> Path:
    out_path.parent.mkdir(parents=True, exist_ok=True)

    if out_path.exists() and out_path.stat().st_size > 0:
        print(f"  [CACHED] {out_path}")
        print(f"           Size:   {human_size(out_path.stat().st_size)}")
        print(f"           SHA256: {sha256_file(out_path)}")
        return out_path

    print(f"  [DOWNLOAD] {url}")

    tmp = out_path.with_suffix(out_path.suffix + ".part")

    with requests.get(url, headers={"User-Agent": UA}, stream=True, timeout=timeout) as r:
        r.raise_for_status()

        with open(tmp, "wb") as f:
            for chunk in r.iter_content(chunk_size=1024 * 1024):
                if chunk:
                    f.write(chunk)

    tmp.replace(out_path)

    print(f"           Size:   {human_size(out_path.stat().st_size)}")
    print(f"           SHA256: {sha256_file(out_path)}")

    return out_path


def download_dataset(
    release: str,
    dataset: str,
    workdir: Path,
) -> Tuple[str, List[Path]]:
    print_section(f"Downloading Open Targets dataset: {dataset}")

    dir_url = f"{BASE}/{release}/output/{dataset}/"
    print(f"[URL] {dir_url}")

    parquet_urls = list_parquet_files(dir_url)

    if not parquet_urls:
        raise RuntimeError(f"No parquet files found at: {dir_url}")

    print(f"[PARQUET FILES] {len(parquet_urls):,}")

    local_paths = []

    for url in parquet_urls:
        fname = url.split("/")[-1]
        local_path = workdir / release / dataset / fname
        local_paths.append(download(url, local_path))

    return dir_url, local_paths


def read_parquets(paths: List[Path], dataset: str) -> pd.DataFrame:
    print_section(f"Reading dataset: {dataset}")

    dfs = []

    for path in paths:
        print(f"[READ] {path}")
        table = pq.read_table(path)
        df = table.to_pandas()

        print(f"  rows={len(df):,}, cols={len(df.columns):,}")

        dfs.append(df)

    if not dfs:
        raise RuntimeError(f"No dataframes read for dataset: {dataset}")

    out = pd.concat(dfs, ignore_index=True)

    print(f"[COMBINED] {dataset}: rows={len(out):,}, cols={len(out.columns):,}")
    print("[COLUMNS]")

    for col in out.columns:
        print(f"  - {col}")

    return out


# =============================================================================
# Mapping functions
# =============================================================================

def build_drug_name_map(drug_df: pd.DataFrame) -> Dict[str, str]:
    """
    drug_molecule expected schema includes:
        id, name
    """
    required = {"id", "name"}
    missing = required - set(drug_df.columns)

    if missing:
        raise RuntimeError(f"drug_molecule is missing required columns: {missing}")

    mapping = {}

    for _, row in drug_df[["id", "name"]].iterrows():
        drug_id = clean_cell(row["id"])
        drug_name = clean_cell(row["name"])

        if drug_id and drug_name:
            mapping[drug_id] = drug_name

    print(f"[DRUG NAME MAP] {len(mapping):,} entries")

    return mapping


def build_disease_name_map(disease_df: pd.DataFrame) -> Dict[str, str]:
    """
    disease expected schema includes:
        id, name
    """
    required = {"id", "name"}
    missing = required - set(disease_df.columns)

    if missing:
        raise RuntimeError(f"disease is missing required columns: {missing}")

    mapping = {}

    for _, row in disease_df[["id", "name"]].iterrows():
        disease_id = clean_cell(row["id"])
        disease_name = clean_cell(row["name"])

        if disease_id and disease_name:
            mapping[disease_id] = disease_name

    print(f"[DISEASE NAME MAP] {len(mapping):,} entries")

    return mapping


def load_optional_chembl_mapping(chembl_csv: Optional[Path]) -> Dict[str, str]:
    """
    Optional fallback only. Open Targets drug_molecule already provides names.
    """
    if chembl_csv is None:
        return {}

    if not chembl_csv.exists():
        print(f"[WARNING] ChEMBL mapping file not found, skipping: {chembl_csv}")
        return {}

    print_section("Loading optional ChEMBL mapping")
    print(f"[ChEMBL CSV] {chembl_csv}")

    df = pd.read_csv(chembl_csv, dtype=str, low_memory=False).fillna("")

    # Support both old and unified ChEMBL outputs.
    if "drug_id" in df.columns and "drug_name" in df.columns:
        id_col = "drug_id"
    elif "drug_identifier" in df.columns and "drug_name" in df.columns:
        id_col = "drug_identifier"
    else:
        print("[WARNING] ChEMBL CSV must contain drug_id/drug_identifier and drug_name. Skipping.")
        return {}

    mapping = {}

    for _, row in df.iterrows():
        drug_id = clean_cell(row[id_col])
        drug_name = clean_cell(row["drug_name"])

        if not drug_id or not drug_name:
            continue

        mapping[drug_id] = drug_name

        if not drug_id.startswith("CHEMBL"):
            mapping[f"CHEMBL{drug_id}"] = drug_name

    print(f"[OPTIONAL ChEMBL MAP] {len(mapping):,} entries")

    return mapping


# =============================================================================
# Main extraction
# =============================================================================

def extract_opentargets_drug_disease(
    clinical_df: pd.DataFrame,
    drug_name_map: Dict[str, str],
    disease_name_map: Dict[str, str],
    chembl_name_map: Dict[str, str],
    release: str,
) -> pd.DataFrame:
    """
    clinical_indication expected schema:
        id
        maxClinicalStage
        clinicalReportIds
        diseaseId
        drugId
    """
    print_section("Extracting Open Targets drug-disease mappings")

    required = {
        "id",
        "maxClinicalStage",
        "clinicalReportIds",
        "diseaseId",
        "drugId",
    }

    missing = required - set(clinical_df.columns)

    if missing:
        raise RuntimeError(
            f"clinical_indication is missing expected columns: {missing}\n"
            f"Available columns: {list(clinical_df.columns)}"
        )

    out = pd.DataFrame()

    out["opentargets_indication_id"] = clinical_df["id"].map(clean_cell)
    out["drug_id"] = clinical_df["drugId"].map(clean_cell)
    out["disease_id"] = clinical_df["diseaseId"].map(clean_cell)
    out["clinical_phase_raw"] = clinical_df["maxClinicalStage"].map(clean_cell)
    out["max_phase"] = clinical_df["maxClinicalStage"].map(safe_phase_to_int).astype(int)
    out["clinical_report_ids"] = clinical_df["clinicalReportIds"].map(clean_cell)

    out["drug_name"] = out["drug_id"].map(
        lambda x: drug_name_map.get(x, chembl_name_map.get(x, x))
    )

    out["disease_name"] = out["disease_id"].map(
        lambda x: disease_name_map.get(x, x)
    )

    # -------------------------------------------------------------------------
    # Unified schema columns
    # -------------------------------------------------------------------------
    out["drug_identifier"] = out["drug_id"]
    out["drug_identifier_type"] = out["drug_id"].map(infer_drug_identifier_type)

    out["disease_or_condition_name"] = out["disease_name"]
    out["disease_or_condition_identifier"] = out["disease_id"]
    out["disease_or_condition_identifier_type"] = out["disease_id"].map(
        infer_disease_identifier_type
    )

    out["relationship_type"] = "clinical_indication"
    out["evidence_type"] = "opentargets_clinical_indication"

    out["source"] = "OpenTargets"
    out["internal_source"] = f"Open Targets Platform {release} clinical_indication"

    out["dataset"] = "clinical_indication"
    out["release"] = release

    # Clean strings
    for col in out.columns:
        if col != "max_phase":
            out[col] = out[col].map(clean_cell)

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
            "max_phase",
            "clinical_phase_raw",
            "opentargets_indication_id",
        ]
    ).reset_index(drop=True)

    removed_duplicates = before_dup - len(out)

    out = out[OUTPUT_COLUMNS]

    print(f"[INPUT ROWS] {len(clinical_df):,}")
    print(f"[OUTPUT ROWS] {len(out):,}")
    print(f"[REMOVED EMPTY DRUG/DISEASE] {removed_empty:,}")
    print(f"[REMOVED DUPLICATES] {removed_duplicates:,}")

    print("\n[OUTPUT PREVIEW]")
    preview_cols = [
        "drug_name",
        "drug_identifier",
        "disease_or_condition_name",
        "disease_or_condition_identifier",
        "relationship_type",
        "evidence_type",
        "max_phase",
        "clinical_phase_raw",
        "release",
    ]

    print(out[preview_cols].head(20).to_string(index=False))

    return out


def summarise_output(df: pd.DataFrame) -> Dict[str, Any]:
    print_section("Output summary")

    summary = {
        "rows": int(len(df)),
        "columns": list(df.columns),
        "unique_drug_identifiers": int(df["drug_identifier"].nunique()),
        "unique_drug_names": int(df["drug_name"].nunique()),
        "unique_disease_identifiers": int(
            df["disease_or_condition_identifier"].nunique()
        ),
        "unique_disease_names": int(
            df["disease_or_condition_name"].nunique()
        ),
        "unique_drug_disease_pairs": int(
            df[
                [
                    "drug_identifier",
                    "disease_or_condition_identifier",
                ]
            ]
            .drop_duplicates()
            .shape[0]
        ),
        "relationship_type_counts": (
            df["relationship_type"].value_counts().to_dict()
            if len(df)
            else {}
        ),
        "evidence_type_counts": (
            df["evidence_type"].value_counts().to_dict()
            if len(df)
            else {}
        ),
        "max_phase_counts": (
            df["max_phase"]
            .value_counts()
            .sort_index()
            .reset_index()
            .rename(columns={"index": "max_phase", "max_phase": "n"})
            .to_dict(orient="records")
        ),
        "clinical_phase_raw_counts": (
            df["clinical_phase_raw"]
            .value_counts()
            .reset_index()
            .rename(columns={"index": "clinical_phase_raw", "clinical_phase_raw": "n"})
            .to_dict(orient="records")
        ),
        "head": df.head(20).to_dict(orient="records"),
    }

    for key, value in summary.items():
        if key not in {
            "head",
            "max_phase_counts",
            "clinical_phase_raw_counts",
            "columns",
        }:
            print(f"{key}: {value}")

    print("\n[MAX PHASE COUNTS]")
    print(df["max_phase"].value_counts().sort_index().to_string())

    print("\n[CLINICAL PHASE RAW COUNTS]")
    print(df["clinical_phase_raw"].value_counts().head(50).to_string())

    return summary


def write_csv(df: pd.DataFrame, out_csv: Path) -> None:
    print_section("Writing CSV")

    out_csv.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_csv, index=False, quoting=csv.QUOTE_ALL)

    print(f"[OUTPUT] {out_csv}")
    print(f"[ROWS] {len(df):,}")
    print(f"[SIZE] {human_size(out_csv.stat().st_size)}")
    print(f"[SHA256] {sha256_file(out_csv)}")


# =============================================================================
# Main
# =============================================================================

def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Generate Open Targets drug-disease mapping from clinical_indication "
            "using the unified multi-source schema."
        )
    )

    parser.add_argument(
        "--release",
        default=DEFAULT_RELEASE,
        help="Open Targets release. Default: 26.03",
    )

    parser.add_argument(
        "--out",
        default="./data/processed/opentargets_drug_disease.csv",
        help="Output CSV path.",
    )

    parser.add_argument(
        "--metadata",
        default="./data/metadata/opentargets_processing_metadata.json",
        help="Output metadata JSON path.",
    )

    parser.add_argument(
        "--workdir",
        default="./data/raw/opentargets",
        help="Download/cache directory.",
    )

    parser.add_argument(
        "--chembl",
        default="",
        help=(
            "Optional local ChEMBL mapping CSV with drug_id/drug_identifier "
            "and drug_name columns."
        ),
    )

    args = parser.parse_args()

    release = args.release
    out_csv = Path(args.out).resolve()
    metadata_path = Path(args.metadata).resolve()
    workdir = Path(args.workdir).resolve()
    chembl_path = Path(args.chembl).resolve() if args.chembl else None

    metadata: Dict[str, Any] = {
        "script": "opentargets.py",
        "source": "Open Targets Platform",
        "release": release,
        "base_url": BASE,
        "started_at_utc": now_utc(),
        "status": "started",
        "datasets": {},
        "relationship_type": "clinical_indication",
        "evidence_type": "opentargets_clinical_indication",
        "interpretation": (
            "Open Targets clinical_indication rows represent drug-disease "
            "clinical indication relationships with clinical stage metadata. "
            "maxClinicalStage is converted into max_phase; APPROVAL is treated "
            "as phase 4. These rows should not all be interpreted as approved "
            "indications unless filtered by phase/stage."
        ),
        "unified_output_columns": OUTPUT_COLUMNS,
    }

    try:
        print_section("Open Targets -> unified drug-disease mapping CSV")
        print(f"[RELEASE] {release}")
        print(f"[OUTPUT] {out_csv}")
        print(f"[METADATA] {metadata_path}")
        print(f"[WORKDIR] {workdir}")

        # Download datasets
        clinical_url, clinical_paths = download_dataset(
            release=release,
            dataset="clinical_indication",
            workdir=workdir,
        )

        drug_url, drug_paths = download_dataset(
            release=release,
            dataset="drug_molecule",
            workdir=workdir,
        )

        disease_url, disease_paths = download_dataset(
            release=release,
            dataset="disease",
            workdir=workdir,
        )

        metadata["datasets"]["clinical_indication"] = {
            "url": clinical_url,
            "files": [str(p) for p in clinical_paths],
        }

        metadata["datasets"]["drug_molecule"] = {
            "url": drug_url,
            "files": [str(p) for p in drug_paths],
        }

        metadata["datasets"]["disease"] = {
            "url": disease_url,
            "files": [str(p) for p in disease_paths],
        }

        # Read datasets
        clinical_df = read_parquets(clinical_paths, "clinical_indication")
        drug_df = read_parquets(drug_paths, "drug_molecule")
        disease_df = read_parquets(disease_paths, "disease")

        metadata["datasets"]["clinical_indication"]["rows"] = int(len(clinical_df))
        metadata["datasets"]["clinical_indication"]["columns"] = list(clinical_df.columns)

        metadata["datasets"]["drug_molecule"]["rows"] = int(len(drug_df))
        metadata["datasets"]["drug_molecule"]["columns"] = list(drug_df.columns)

        metadata["datasets"]["disease"]["rows"] = int(len(disease_df))
        metadata["datasets"]["disease"]["columns"] = list(disease_df.columns)

        # Build maps
        drug_name_map = build_drug_name_map(drug_df)
        disease_name_map = build_disease_name_map(disease_df)
        chembl_name_map = load_optional_chembl_mapping(chembl_path)

        metadata["name_mapping"] = {
            "drug_molecule_entries": len(drug_name_map),
            "disease_entries": len(disease_name_map),
            "optional_chembl_entries": len(chembl_name_map),
            "optional_chembl_path": str(chembl_path) if chembl_path else "",
        }

        # Extract output
        out_df = extract_opentargets_drug_disease(
            clinical_df=clinical_df,
            drug_name_map=drug_name_map,
            disease_name_map=disease_name_map,
            chembl_name_map=chembl_name_map,
            release=release,
        )

        if out_df.empty:
            raise RuntimeError("Final Open Targets output is empty.")

        write_csv(out_df, out_csv)

        metadata["output"] = {
            "path": str(out_csv),
            "rows": int(len(out_df)),
            "size_bytes": int(out_csv.stat().st_size),
            "sha256": sha256_file(out_csv),
        }

        metadata["summary"] = summarise_output(out_df)
        metadata["finished_at_utc"] = now_utc()
        metadata["status"] = "success"

        print_section("DONE")
        print("[SUCCESS] Open Targets drug-disease mapping generated.")
        print(f"[CSV] {out_csv}")
        print(f"[METADATA] {metadata_path}")

    except Exception as e:
        metadata["finished_at_utc"] = now_utc()
        metadata["status"] = "failed"
        metadata["error"] = repr(e)

        print_section("FAILED")
        print(repr(e))

        raise

    finally:
        metadata_path.parent.mkdir(parents=True, exist_ok=True)

        with open(metadata_path, "w", encoding="utf-8") as f:
            json.dump(metadata, f, indent=2, default=str)

        print_section("METADATA WRITTEN")
        print(metadata_path)


if __name__ == "__main__":
    main()