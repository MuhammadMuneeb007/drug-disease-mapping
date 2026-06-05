#!/usr/bin/env python3
"""
MEDI-2 (MEDication-Indication) ensemble resource processing.

Purpose
-------
Process MEDI-2 ICD-coded and UMLS-coded medication-indication files into
a publication-ready, provenance-aware CSV using the unified drug-disease
schema.

Source page:
    https://www.vumc.org/wei-lab/medi

Input files:
    MEDI-2.csv
    MEDI-2_UMLS.csv

Input URLs:
    https://www.vumc.org/wei-lab/sites/default/files/public_files/MEDI-2.csv
    https://www.vumc.org/wei-lab/sites/default/files/public_files/MEDI-2_UMLS.csv

Citation:
    Wei WQ et al. Development and evaluation of an ensemble resource linking
    medications to their indications. Journal of the American Medical Informatics
    Association. 2013. doi:10.1136/amiajnl-2012-001431

Licence:
    CC BY-NC-SA 3.0

Important licence warning
-------------------------
MEDI is non-commercial and share-alike. Including MEDI-derived rows in a
redistributed merged dataset may impose non-commercial/share-alike obligations
on the merged work. If the final dataset is intended to be CC0 or CC BY,
consider distributing the processing script and metadata only, and let users
rebuild the MEDI-derived table locally.

Important interpretation
------------------------
MEDI records represent medication-indication ensemble evidence. They should
not be interpreted as regulatory approval evidence unless externally verified.

HIGH_PRECISION_SUBSET is treated as an evidence-quality flag, not a different
biological relationship type. Therefore:
    relationship_type = indication

Evidence type is used to distinguish:
    medi_high_precision_indication
    medi_ensemble_indication

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

MEDI-specific metadata columns
------------------------------
high_precision_subset
number_of_resources
resources_mentioning
source_record_type
evidence_text

Requirements
------------
pip install pandas
"""

from __future__ import annotations

import hashlib
import json
import os
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd


# ============================================================
# CONFIG
# ============================================================

MEDI_BASE_URL = "https://www.vumc.org/wei-lab/sites/default/files/public_files"

MEDI_ICD_URL = f"{MEDI_BASE_URL}/MEDI-2.csv"
MEDI_UMLS_URL = f"{MEDI_BASE_URL}/MEDI-2_UMLS.csv"

RAW_DIR = (
    "/data/ascher02/uqmmune1/ANNOVAR/drug_disease_data/"
    "drug-disease-mapping/data/raw/medi"
)

OUTPUT_DIR = (
    "/data/ascher02/uqmmune1/ANNOVAR/drug_disease_data/"
    "drug-disease-mapping/data/processed/medi"
)

MEDI_ICD_FILE = os.path.join(RAW_DIR, "MEDI-2.csv")
MEDI_UMLS_FILE = os.path.join(RAW_DIR, "MEDI-2_UMLS.csv")

OUTPUT_FILE = os.path.join(OUTPUT_DIR, "medi_drug_disease.csv")
DEDUP_FILE = os.path.join(OUTPUT_DIR, "medi_drug_disease_deduplicated.csv")
SUMMARY_FILE = os.path.join(OUTPUT_DIR, "medi_summary.json")


RESOURCE_FLAG_COLS = [
    "IN_RXNORM",
    "IN_MAYO",
    "IN_MEDLINE",
    "IN_SIDER",
    "IN_WEBMD",
    "IN_WIKIPEDIA",
]


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

    # MEDI-specific metadata
    "high_precision_subset",
    "number_of_resources",
    "resources_mentioning",
    "source_record_type",
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


def is_truthy(x) -> bool:
    """
    MEDI commonly stores booleans as TRUE/FALSE strings.
    """
    if pd.isna(x):
        return False

    s = str(x).strip().upper()

    return s in {"TRUE", "T", "1", "YES", "Y"}


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


def download_if_missing(url: str, path: str, force: bool = False) -> None:
    path_obj = Path(path)

    if path_obj.exists() and path_obj.stat().st_size > 0 and not force:
        print(f"[CACHE] Already present: {path}")
        print(f"        Size:   {human_size(path_obj.stat().st_size)}")
        print(f"        SHA256: {sha256_file(path_obj)}")
        return

    os.makedirs(path_obj.parent, exist_ok=True)

    print(f"[DOWNLOAD] {url}")
    urllib.request.urlretrieve(url, path)

    if not path_obj.exists() or path_obj.stat().st_size == 0:
        raise RuntimeError(f"Downloaded file is missing or empty: {path}")

    print(f"[DOWNLOADED] {path}")
    print(f"             Size:   {human_size(path_obj.stat().st_size)}")
    print(f"             SHA256: {sha256_file(path_obj)}")


def get_resource_list(row: pd.Series) -> list[str]:
    """
    Return list of MEDI resources that mention this drug-indication pair.
    """
    resources = []

    for col in RESOURCE_FLAG_COLS:
        if col in row.index and is_truthy(row[col]):
            resources.append(col.replace("IN_", "").lower())

    return resources


def validate_columns(df: pd.DataFrame, required_cols: list[str], file_label: str) -> None:
    missing = [c for c in required_cols if c not in df.columns]

    if missing:
        raise ValueError(f"{file_label} missing required columns: {missing}")


# ============================================================
# PROCESSORS
# ============================================================

def process_medi_icd(file_path: str) -> list[dict]:
    """
    Process MEDI-2.csv.

    This file contains ICD-9-CM / ICD-10-CM coded indications and
    provenance/resource flags.
    """
    print("=" * 100)
    print("Processing MEDI-2 ICD-coded file")
    print("=" * 100)

    df = pd.read_csv(file_path, dtype=str, keep_default_na=False)

    print(f"[INFO] Rows loaded: {len(df):,}")
    print(f"[INFO] Columns: {list(df.columns)}")

    required_cols = [
        "RXCUI",
        "DRUG_DESCRIPTION",
        "VOCABULARY",
        "CODE",
        "INDICATION_DESCRIPTION",
        "HIGH_PRECISION_SUBSET",
        "NUMBER_OF_RESOURCES_MENTIONED",
    ]

    validate_columns(df, required_cols, "MEDI-2.csv")

    rows = []
    skipped_missing_drug = 0
    skipped_missing_disease = 0

    for _, row in df.iterrows():
        rxcui = clean_text(row["RXCUI"])
        drug_name = clean_text(row["DRUG_DESCRIPTION"])

        vocabulary = clean_text(row["VOCABULARY"])
        code = clean_text(row["CODE"])
        disease_name = clean_text(row["INDICATION_DESCRIPTION"])

        hps = is_truthy(row["HIGH_PRECISION_SUBSET"])
        n_resources = clean_text(row["NUMBER_OF_RESOURCES_MENTIONED"])
        resources = get_resource_list(row)

        if not drug_name:
            skipped_missing_drug += 1
            continue

        if not disease_name:
            skipped_missing_disease += 1
            continue

        # Important:
        # HPS is an evidence-quality flag, not a different relationship type.
        relationship_type = "indication"

        if hps:
            evidence_type = "medi_high_precision_indication"
        else:
            evidence_type = "medi_ensemble_indication"

        evidence_parts = [
            f"source_record_type=ICD-coded",
            f"HIGH_PRECISION_SUBSET={hps}",
            f"NUMBER_OF_RESOURCES_MENTIONED={n_resources}",
        ]

        if resources:
            evidence_parts.append(f"resources={'|'.join(resources)}")

        rows.append(
            {
                "drug_name": drug_name,
                "drug_identifier": rxcui,
                "drug_identifier_type": "RxCUI",
                "disease_or_condition_name": disease_name,
                "disease_or_condition_identifier": code,
                "disease_or_condition_identifier_type": vocabulary,
                "relationship_type": relationship_type,
                "evidence_type": evidence_type,
                "source": "MEDI",
                "internal_source": "MEDI-2 ICD-coded",
                "high_precision_subset": str(hps),
                "number_of_resources": n_resources,
                "resources_mentioning": "|".join(resources),
                "source_record_type": "ICD-coded",
                "evidence_text": "; ".join(evidence_parts),
            }
        )

    print(f"[INFO] Standardised ICD-coded rows: {len(rows):,}")
    print(f"[INFO] Skipped missing drug:        {skipped_missing_drug:,}")
    print(f"[INFO] Skipped missing disease:     {skipped_missing_disease:,}")

    return rows


def process_medi_umls(file_path: str) -> list[dict]:
    """
    Process MEDI-2_UMLS.csv.

    This file contains UMLS CUI-coded indications and is useful for
    cross-source disease concept integration.
    """
    print("=" * 100)
    print("Processing MEDI-2 UMLS-coded file")
    print("=" * 100)

    df = pd.read_csv(file_path, dtype=str, keep_default_na=False)

    print(f"[INFO] Rows loaded: {len(df):,}")
    print(f"[INFO] Columns: {list(df.columns)}")

    required_cols = [
        "RXCUI",
        "DRUG_DESC",
        "CUI",
        "STR",
    ]

    validate_columns(df, required_cols, "MEDI-2_UMLS.csv")

    rows = []
    skipped_missing_drug = 0
    skipped_missing_disease = 0

    for _, row in df.iterrows():
        rxcui = clean_text(row["RXCUI"])
        drug_name = clean_text(row["DRUG_DESC"])

        umls_cui = clean_text(row["CUI"])
        disease_name = clean_text(row["STR"])

        if not drug_name:
            skipped_missing_drug += 1
            continue

        if not disease_name:
            skipped_missing_disease += 1
            continue

        rows.append(
            {
                "drug_name": drug_name,
                "drug_identifier": rxcui,
                "drug_identifier_type": "RxCUI",
                "disease_or_condition_name": disease_name,
                "disease_or_condition_identifier": umls_cui,
                "disease_or_condition_identifier_type": "UMLS",
                "relationship_type": "indication",
                "evidence_type": "medi_ensemble_indication",
                "source": "MEDI",
                "internal_source": "MEDI-2 UMLS-coded",
                "high_precision_subset": "",
                "number_of_resources": "",
                "resources_mentioning": "",
                "source_record_type": "UMLS-coded",
                "evidence_text": "source_record_type=UMLS-coded; source=MEDI-2_UMLS",
            }
        )

    print(f"[INFO] Standardised UMLS-coded rows: {len(rows):,}")
    print(f"[INFO] Skipped missing drug:         {skipped_missing_drug:,}")
    print(f"[INFO] Skipped missing disease:      {skipped_missing_disease:,}")

    return rows


# ============================================================
# MAIN
# ============================================================

def main() -> None:
    print("=" * 100)
    print("MEDI-2 PROCESSING")
    print("=" * 100)
    print("[SOURCE] MEDI-2 / MEDication-Indication resource")
    print(f"[ICD URL]  {MEDI_ICD_URL}")
    print(f"[UMLS URL] {MEDI_UMLS_URL}")
    print("[DOI] 10.1136/amiajnl-2012-001431")
    print("[LICENSE] CC BY-NC-SA 3.0")
    print("[WARNING] Non-commercial and share-alike licence.")
    print("=" * 100)

    os.makedirs(RAW_DIR, exist_ok=True)
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    download_if_missing(MEDI_ICD_URL, MEDI_ICD_FILE)
    download_if_missing(MEDI_UMLS_URL, MEDI_UMLS_FILE)

    icd_rows = process_medi_icd(MEDI_ICD_FILE)
    umls_rows = process_medi_umls(MEDI_UMLS_FILE)

    out = pd.DataFrame(icd_rows + umls_rows, columns=OUTPUT_COLUMNS)

    print("=" * 100)
    print("OUTPUT SUMMARY")
    print("=" * 100)
    print(f"[INFO] Combined output rows: {len(out):,}")

    if len(out) == 0:
        raise RuntimeError("No MEDI rows were produced.")

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
        "source_record_type",
    ]

    dedup = out.drop_duplicates(subset=dedup_cols)
    dedup.to_csv(DEDUP_FILE, index=False)

    icd_subset = out[out["source_record_type"] == "ICD-coded"]
    umls_subset = out[out["source_record_type"] == "UMLS-coded"]
    hps_subset = out[out["evidence_type"] == "medi_high_precision_indication"]

    relationship_type_counts = out["relationship_type"].value_counts().to_dict()
    evidence_type_counts = out["evidence_type"].value_counts().to_dict()
    identifier_type_counts = (
        out["disease_or_condition_identifier_type"]
        .value_counts()
        .to_dict()
    )
    source_record_type_counts = out["source_record_type"].value_counts().to_dict()

    summary = {
        "timestamp_utc": now_utc(),
        "source": "MEDI-2",
        "source_page": "https://www.vumc.org/wei-lab/medi",
        "icd_url": MEDI_ICD_URL,
        "umls_url": MEDI_UMLS_URL,
        "citation_doi": "10.1136/amiajnl-2012-001431",
        "license": "CC BY-NC-SA 3.0",
        "license_warning": (
            "MEDI is non-commercial and share-alike. Inclusion of MEDI-derived "
            "rows in a redistributed merged dataset may impose NC-SA obligations "
            "on the merged work."
        ),
        "interpretation": (
            "MEDI records represent medication-indication ensemble evidence. "
            "HIGH_PRECISION_SUBSET is treated as an evidence-quality flag, not "
            "as a separate biological relationship type. Therefore all MEDI "
            "records use relationship_type=indication, while evidence_type "
            "distinguishes high-precision and ensemble records."
        ),
        "icd_file": MEDI_ICD_FILE,
        "umls_file": MEDI_UMLS_FILE,
        "icd_file_sha256": sha256_file(MEDI_ICD_FILE),
        "umls_file_sha256": sha256_file(MEDI_UMLS_FILE),
        "icd_rows_written": int(len(icd_subset)),
        "umls_rows_written": int(len(umls_subset)),
        "high_precision_rows": int(len(hps_subset)),
        "total_rows_written": int(len(out)),
        "deduplicated_rows_written": int(len(dedup)),
        "unique_drugs_rxcui": int(out["drug_identifier"].nunique()),
        "unique_diseases_combined": int(
            out[
                [
                    "disease_or_condition_identifier_type",
                    "disease_or_condition_identifier",
                ]
            ]
            .drop_duplicates()
            .shape[0]
        ),
        "source_record_type_counts": source_record_type_counts,
        "identifier_type_counts": identifier_type_counts,
        "relationship_type_counts": relationship_type_counts,
        "evidence_type_counts": evidence_type_counts,
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