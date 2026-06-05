#!/usr/bin/env python3
"""
repoDB -> unified drug-disease indication/status CSV

Purpose
-------
Process repoDB drug repositioning / indication-status records into a
publication-ready, provenance-aware CSV using the unified multi-source
drug-disease schema.

Input
-----
Expected local file:
    /data/ascher02/uqmmune1/ANNOVAR/drug_disease_data/
    drug-disease-mapping/data/raw/repodb/repodb_full_download.tsv

Source
------
repoDB web app:
    https://unmtid-shinyapps.net/shiny/repodb/

Citation
--------
Brown AS, Patel CJ.
A standard database for drug repositioning.
Scientific Data. 2017;4:170029.
doi:10.1038/sdata.2017.29

Important interpretation
------------------------
repoDB records represent drug repositioning / drug-indication status records.
Rows include approved, failed, terminated, withdrawn, suspended and other
status categories. Therefore, do not collapse all repoDB rows into approved
indications.

For status = approved:
    relationship_type = indication
    evidence_type = repodb_approved_indication

For failed/terminated/withdrawn/suspended:
    relationship_type = failed_or_discontinued_indication
    evidence_type = repodb_failed_or_discontinued_indication

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

repoDB-specific metadata columns
--------------------------------
status
phase
nct_id
detailed_status
evidence_text

Requirements
------------
pip install pandas
"""

from __future__ import annotations

import hashlib
import json
import os
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd


# ============================================================
# CONFIG
# ============================================================

INPUT_FILE = (
    "/data/ascher02/uqmmune1/ANNOVAR/drug_disease_data/"
    "drug-disease-mapping/data/raw/repodb/repodb_full_download.tsv"
)

OUTPUT_DIR = (
    "/data/ascher02/uqmmune1/ANNOVAR/drug_disease_data/"
    "drug-disease-mapping/data/processed/repodb"
)

OUTPUT_FILE = os.path.join(OUTPUT_DIR, "repodb_drug_disease.csv")
DEDUP_FILE = os.path.join(OUTPUT_DIR, "repodb_drug_disease_deduplicated.csv")
SUMMARY_FILE = os.path.join(OUTPUT_DIR, "repodb_summary.json")


SOURCE_URL = "https://unmtid-shinyapps.net/shiny/repodb/"
CITATION_DOI = "10.1038/sdata.2017.29"


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

    # repoDB-specific metadata
    "status",
    "phase",
    "nct_id",
    "detailed_status",
    "evidence_text",
]


# ============================================================
# HELPERS
# ============================================================

def now_utc() -> str:
    return datetime.now(timezone.utc).isoformat()


def clean_text(x) -> str:
    if pd.isna(x):
        return ""
    return str(x).strip()


def sha256_file(path: str | Path, chunk_size: int = 1024 * 1024) -> str:
    path = Path(path)
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


def detect_separator(file_path: str) -> str:
    with open(file_path, "r", encoding="utf-8", errors="ignore") as f:
        first_line = f.readline()

    if "\t" in first_line:
        return "\t"

    return ","


def classify_relationship(status: str) -> tuple[str, str]:
    """
    Convert repoDB status into harmonised relationship/evidence labels.

    Important:
    Do not relabel approved repoDB records as approved_indication.
    Keep approved status in evidence_type/status metadata.
    """
    status_clean = clean_text(status).lower()

    if status_clean == "approved":
        return "indication", "repodb_approved_indication"

    if status_clean in {
        "failed",
        "terminated",
        "withdrawn",
        "suspended",
    }:
        return (
            "failed_or_discontinued_indication",
            "repodb_failed_or_discontinued_indication",
        )

    if status_clean:
        safe_status = status_clean.replace(" ", "_").replace("/", "_")
        return (
            "drug_indication_with_status",
            f"repodb_{safe_status}",
        )

    return "drug_indication", "repodb_indication_status_unspecified"


def normalise_disease_identifier_type(ind_id: str) -> str:
    """
    repoDB ind_id values are commonly UMLS CUIs, but keep this defensive.
    """
    ind_id = clean_text(ind_id)

    if not ind_id:
        return ""

    upper = ind_id.upper()

    if upper.startswith("C") and upper[1:].isdigit():
        return "UMLS"

    if upper.startswith("UMLS:"):
        return "UMLS"

    if upper.startswith("MESH:"):
        return "MeSH"

    if upper.startswith("OMIM:"):
        return "OMIM"

    if upper.startswith("MONDO:"):
        return "MONDO"

    return "repoDB_IndicationID"


def normalise_drug_identifier_type(drugbank_id: str) -> str:
    drugbank_id = clean_text(drugbank_id)

    if not drugbank_id:
        return ""

    if drugbank_id.upper().startswith("DB"):
        return "DrugBank"

    return "repoDB_DrugID"


# ============================================================
# MAIN
# ============================================================

def main() -> None:
    print("=" * 100)
    print("repoDB PROCESSING")
    print("=" * 100)
    print("[SOURCE] repoDB")
    print(f"[SOURCE URL] {SOURCE_URL}")
    print(f"[DOI] {CITATION_DOI}")
    print("[INTERPRETATION] Drug repositioning / indication-status evidence")
    print("=" * 100)

    os.makedirs(OUTPUT_DIR, exist_ok=True)

    if not os.path.exists(INPUT_FILE):
        raise FileNotFoundError(
            f"Input file not found: {INPUT_FILE}\n"
            "Download repoDB full file first and save it at this path."
        )

    sep = detect_separator(INPUT_FILE)

    print(f"[INFO] Reading file: {INPUT_FILE}")
    print(f"[INFO] Detected separator: {repr(sep)}")
    print(f"[INFO] Input size: {human_size(Path(INPUT_FILE).stat().st_size)}")
    print(f"[INFO] Input SHA256: {sha256_file(INPUT_FILE)}")

    df = pd.read_csv(
        INPUT_FILE,
        sep=sep,
        dtype=str,
        keep_default_na=False,
        low_memory=False,
    )

    print(f"[INFO] Rows loaded: {len(df):,}")
    print(f"[INFO] Columns: {list(df.columns)}")

    required_cols = [
        "drug_name",
        "drugbank_id",
        "ind_name",
        "ind_id",
        "NCT",
        "status",
        "phase",
    ]

    missing = [c for c in required_cols if c not in df.columns]

    if missing:
        raise ValueError(
            f"Missing required repoDB columns: {missing}\n"
            f"Available columns: {list(df.columns)}"
        )

    rows = []
    skipped_missing_drug = 0
    skipped_missing_disease = 0

    for _, row in df.iterrows():
        drug_name = clean_text(row["drug_name"])
        drugbank_id = clean_text(row["drugbank_id"])

        disease_name = clean_text(row["ind_name"])
        disease_id = clean_text(row["ind_id"])

        nct_id = clean_text(row["NCT"])
        status = clean_text(row["status"])
        phase = clean_text(row["phase"])

        detailed_status = (
            clean_text(row["DetailedStatus"])
            if "DetailedStatus" in df.columns
            else ""
        )

        if not drug_name:
            skipped_missing_drug += 1
            continue

        if not disease_name:
            skipped_missing_disease += 1
            continue

        relationship_type, evidence_type = classify_relationship(status)

        evidence_parts = []

        if status:
            evidence_parts.append(f"status={status}")

        if phase:
            evidence_parts.append(f"phase={phase}")

        if nct_id and nct_id.upper() != "NA":
            evidence_parts.append(f"NCT={nct_id}")

        if detailed_status and detailed_status.upper() != "NA":
            evidence_parts.append(f"DetailedStatus={detailed_status}")

        rows.append(
            {
                "drug_name": drug_name,
                "drug_identifier": drugbank_id,
                "drug_identifier_type": normalise_drug_identifier_type(
                    drugbank_id
                ),
                "disease_or_condition_name": disease_name,
                "disease_or_condition_identifier": disease_id,
                "disease_or_condition_identifier_type": normalise_disease_identifier_type(
                    disease_id
                ),
                "relationship_type": relationship_type,
                "evidence_type": evidence_type,
                "source": "repoDB",
                "internal_source": "repoDB full database",
                "status": status,
                "phase": phase,
                "nct_id": nct_id,
                "detailed_status": detailed_status,
                "evidence_text": "; ".join(evidence_parts),
            }
        )

    out = pd.DataFrame(rows, columns=OUTPUT_COLUMNS)

    print(f"[INFO] Standardised rows: {len(out):,}")
    print(f"[INFO] Skipped missing drug: {skipped_missing_drug:,}")
    print(f"[INFO] Skipped missing disease: {skipped_missing_disease:,}")

    if len(out) == 0:
        raise RuntimeError("No repoDB rows were produced.")

    print()
    print("[OUTPUT PREVIEW]")
    print(out.head(10).to_string(index=False))

    out.to_csv(OUTPUT_FILE, index=False)

    dedup_cols = [
        "drug_name",
        "drug_identifier",
        "disease_or_condition_name",
        "disease_or_condition_identifier",
        "relationship_type",
        "evidence_type",
        "status",
        "phase",
    ]

    dedup = out.drop_duplicates(subset=dedup_cols)
    dedup.to_csv(DEDUP_FILE, index=False)

    status_counts = out["status"].value_counts().to_dict()
    phase_counts = out["phase"].value_counts().to_dict()
    relationship_type_counts = out["relationship_type"].value_counts().to_dict()
    evidence_type_counts = out["evidence_type"].value_counts().to_dict()

    summary = {
        "timestamp_utc": now_utc(),
        "source": "repoDB",
        "source_url": SOURCE_URL,
        "citation_doi": CITATION_DOI,
        "interpretation": (
            "repoDB records represent drug repositioning / indication-status "
            "evidence. Approved, failed, terminated, withdrawn, suspended and "
            "other statuses are preserved rather than collapsed into a single "
            "approved-indication class."
        ),
        "input_file": INPUT_FILE,
        "input_size_bytes": int(Path(INPUT_FILE).stat().st_size),
        "input_sha256": sha256_file(INPUT_FILE),
        "rows_loaded": int(len(df)),
        "rows_written": int(len(out)),
        "deduplicated_rows_written": int(len(dedup)),
        "skipped_missing_drug": int(skipped_missing_drug),
        "skipped_missing_disease": int(skipped_missing_disease),
        "unique_drugs": int(out["drug_identifier"].nunique()),
        "unique_drug_names": int(out["drug_name"].nunique()),
        "unique_diseases": int(out["disease_or_condition_identifier"].nunique()),
        "unique_disease_names": int(out["disease_or_condition_name"].nunique()),
        "status_counts": status_counts,
        "phase_counts": phase_counts,
        "relationship_type_counts": relationship_type_counts,
        "evidence_type_counts": evidence_type_counts,
        "drug_identifier_type_counts": (
            out["drug_identifier_type"].value_counts().to_dict()
        ),
        "disease_identifier_type_counts": (
            out["disease_or_condition_identifier_type"].value_counts().to_dict()
        ),
        "output_file": OUTPUT_FILE,
        "deduplicated_file": DEDUP_FILE,
        "output_file_sha256": sha256_file(OUTPUT_FILE),
        "deduplicated_file_sha256": sha256_file(DEDUP_FILE),
        "output_columns": OUTPUT_COLUMNS,
        "deduplication_columns": dedup_cols,
    }

    with open(SUMMARY_FILE, "w") as f:
        json.dump(summary, f, indent=2)

    print("=" * 100)
    print("DONE")
    print("=" * 100)
    print(f"[OUTPUT]  {OUTPUT_FILE}")
    print(f"[DEDUP]   {DEDUP_FILE}")
    print(f"[SUMMARY] {SUMMARY_FILE}")
    print("=" * 100)


if __name__ == "__main__":
    main()