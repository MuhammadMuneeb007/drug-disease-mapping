#!/usr/bin/env python3
"""
PrimeKG (Precision Medicine Knowledge Graph) -> unified drug-disease CSV

Purpose
-------
Process PrimeKG v2.1 drug-disease relationships into a publication-ready,
provenance-aware CSV using the unified multi-source drug-disease schema.

Source repository:
    https://github.com/mims-harvard/PrimeKG

Data:
    https://doi.org/10.7910/DVN/IXA7BM
    Harvard Dataverse, PrimeKG v2.1

Direct Dataverse file URL:
    https://dataverse.harvard.edu/api/access/datafile/6180620

Citation:
    Chandak P, Huang K, Zitnik M.
    Building a knowledge graph to enable precision medicine.
    Scientific Data. 2023.
    doi:10.1038/s41597-023-01960-3

Important fix
-------------
Dataverse may return HTTP 403 with urllib.urlretrieve().
This script uses requests with browser-like headers and redirect handling.

If you already downloaded the file using:

    wget https://dataverse.harvard.edu/api/access/datafile/6180620

and it saved as:

    6180620

this script will automatically detect that file and copy it to:

    data/raw/primekg/kg.csv

Important interpretation
------------------------
PrimeKG is a biomedical knowledge graph. Its drug-disease relations should
be interpreted as knowledge-graph relationships, not automatically as
regulatory approval evidence.

PrimeKG relation labels retained:
    indication
    contraindication
    off-label use

Do NOT relabel PrimeKG "indication" as "approved_indication".
"""

from __future__ import annotations

import hashlib
import json
import os
import shutil
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import requests


# ============================================================
# CONFIG
# ============================================================

PRIMEKG_KG_URL = "https://dataverse.harvard.edu/api/access/datafile/6180620"

RAW_DIR = (
    "/data/ascher02/uqmmune1/ANNOVAR/drug_disease_data/"
    "drug-disease-mapping/data/raw/primekg"
)

OUTPUT_DIR = (
    "/data/ascher02/uqmmune1/ANNOVAR/drug_disease_data/"
    "drug-disease-mapping/data/processed/primekg"
)

KG_FILE = os.path.join(RAW_DIR, "kg.csv")

OUTPUT_FILE = os.path.join(OUTPUT_DIR, "primekg_drug_disease.csv")
DEDUP_FILE = os.path.join(OUTPUT_DIR, "primekg_drug_disease_deduplicated.csv")
SUMMARY_FILE = os.path.join(OUTPUT_DIR, "primekg_summary.json")


# IMPORTANT:
# PrimeKG "indication" is kept as "indication".
# It is NOT converted to "approved_indication".
DRUG_DISEASE_RELATIONS = {
    "indication": (
        "indication",
        "primekg_indication",
    ),
    "contraindication": (
        "contraindication",
        "primekg_contraindication",
    ),
    "off-label use": (
        "off_label_use",
        "primekg_off_label_use",
    ),
}


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

    # PrimeKG-specific metadata
    "primekg_relation",
    "primekg_display_relation",
    "source_node_type",
    "target_node_type",
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


def looks_like_primekg_csv(path: str | Path) -> bool:
    """
    Quick validation that the file looks like PrimeKG kg.csv.
    """
    path = Path(path)

    if not path.exists() or path.stat().st_size == 0:
        return False

    try:
        with open(path, "r", encoding="utf-8", errors="ignore") as f:
            header = f.readline().strip().split(",")

        required = {
            "relation",
            "display_relation",
            "x_index",
            "x_id",
            "x_type",
            "x_name",
            "x_source",
            "y_index",
            "y_id",
            "y_type",
            "y_name",
            "y_source",
        }

        return required.issubset(set(header))

    except Exception:
        return False


def find_existing_manual_download(target_path: str | Path) -> Path | None:
    """
    Detect common manually downloaded Dataverse file names.

    Your wget command saved the file as:
        6180620

    This function searches current directory and project directory for that file.
    """
    target_path = Path(target_path)

    candidates = [
        Path.cwd() / "6180620",
        Path.cwd() / "kg.csv",
        target_path.parent / "6180620",
        target_path.parent / "kg.csv",
        Path("/data/ascher02/uqmmune1/ANNOVAR/drug_disease_data/drug-disease-mapping/6180620"),
    ]

    for candidate in candidates:
        if candidate.exists() and candidate.stat().st_size > 0:
            if looks_like_primekg_csv(candidate):
                return candidate

    return None


def prepare_existing_or_download(url: str, path: str, force: bool = False) -> None:
    """
    Prepare kg.csv.

    Priority:
      1. Use existing valid kg.csv.
      2. Detect manually downloaded file named 6180620.
      3. Download from Dataverse using requests.
    """
    path_obj = Path(path)
    path_obj.parent.mkdir(parents=True, exist_ok=True)

    if path_obj.exists() and path_obj.stat().st_size > 0 and not force:
        print(f"[CACHE] Already present: {path_obj}")
        print(f"        Size:   {human_size(path_obj.stat().st_size)}")
        print(f"        SHA256: {sha256_file(path_obj)}")

        if not looks_like_primekg_csv(path_obj):
            raise RuntimeError(
                f"Existing file does not look like PrimeKG kg.csv: {path_obj}"
            )

        return

    if path_obj.exists() and force:
        path_obj.unlink()

    manual = find_existing_manual_download(path_obj)

    if manual is not None:
        print(f"[INFO] Found manually downloaded PrimeKG file: {manual}")
        print(f"[INFO] Copying to expected path: {path_obj}")

        if manual.resolve() != path_obj.resolve():
            shutil.copy2(manual, path_obj)

        print(f"[READY] {path_obj}")
        print(f"        Size:   {human_size(path_obj.stat().st_size)}")
        print(f"        SHA256: {sha256_file(path_obj)}")

        return

    print(f"[DOWNLOAD] {url}")
    print("[INFO] PrimeKG kg.csv is large. This may take a few minutes.")

    headers = {
        "User-Agent": (
            "Mozilla/5.0 (X11; Linux x86_64) "
            "AppleWebKit/537.36 "
            "Chrome/120.0 Safari/537.36"
        ),
        "Accept": "text/csv,application/octet-stream,*/*",
        "Connection": "keep-alive",
    }

    tmp_path = path_obj.with_suffix(path_obj.suffix + ".part")

    with requests.get(
        url,
        headers=headers,
        stream=True,
        timeout=1800,
        allow_redirects=True,
    ) as response:
        response.raise_for_status()

        with open(tmp_path, "wb") as f:
            downloaded = 0

            for chunk in response.iter_content(chunk_size=1024 * 1024):
                if chunk:
                    f.write(chunk)
                    downloaded += len(chunk)

                    if downloaded % (100 * 1024 * 1024) < 1024 * 1024:
                        print(f"[DOWNLOAD] {human_size(downloaded)} downloaded")

    tmp_path.replace(path_obj)

    if not path_obj.exists() or path_obj.stat().st_size == 0:
        raise RuntimeError(f"Downloaded PrimeKG file is missing or empty: {path_obj}")

    if not looks_like_primekg_csv(path_obj):
        raise RuntimeError(
            f"Downloaded file does not look like PrimeKG kg.csv: {path_obj}\n"
            "If wget worked, move the downloaded file manually with:\n"
            f"mkdir -p {path_obj.parent}\n"
            f"mv 6180620 {path_obj}"
        )

    print(f"[DOWNLOADED] {path_obj}")
    print(f"             Size:   {human_size(path_obj.stat().st_size)}")
    print(f"             SHA256: {sha256_file(path_obj)}")


def infer_drug_identifier_type(source: str, identifier: str) -> str:
    source = clean_text(source)
    identifier = clean_text(identifier)

    if source:
        return source

    upper = identifier.upper()

    if upper.startswith("DB"):
        return "DrugBank"

    if upper.startswith("CHEMBL"):
        return "ChEMBL"

    return "PrimeKG_DrugID"


def infer_disease_identifier_type(source: str, identifier: str) -> str:
    source = clean_text(source)
    identifier = clean_text(identifier)

    if source:
        return source

    upper = identifier.upper()

    if upper.startswith("MONDO"):
        return "MONDO"

    if upper.startswith("DOID"):
        return "DOID"

    if upper.startswith("MESH"):
        return "MeSH"

    if upper.startswith("OMIM"):
        return "OMIM"

    return "PrimeKG_DiseaseID"


# ============================================================
# MAIN
# ============================================================

def main() -> None:
    print("=" * 100)
    print("PrimeKG DRUG-DISEASE PROCESSING")
    print("=" * 100)
    print("[SOURCE] PrimeKG v2.1")
    print("[DATA URL] https://doi.org/10.7910/DVN/IXA7BM")
    print("[KG FILE URL] https://dataverse.harvard.edu/api/access/datafile/6180620")
    print("[DOI] 10.1038/s41597-023-01960-3")
    print("=" * 100)

    os.makedirs(RAW_DIR, exist_ok=True)
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    prepare_existing_or_download(PRIMEKG_KG_URL, KG_FILE)

    print(f"[INFO] Reading: {KG_FILE}")

    df = pd.read_csv(
        KG_FILE,
        dtype=str,
        keep_default_na=False,
        low_memory=False,
    )

    print(f"[INFO] Total edges loaded: {len(df):,}")
    print(f"[INFO] Columns: {list(df.columns)}")

    required_cols = [
        "relation",
        "display_relation",
        "x_index",
        "x_id",
        "x_type",
        "x_name",
        "x_source",
        "y_index",
        "y_id",
        "y_type",
        "y_name",
        "y_source",
    ]

    missing = [c for c in required_cols if c not in df.columns]

    if missing:
        raise ValueError(f"PrimeKG kg.csv missing required columns: {missing}")

    # ------------------------------------------------------------
    # Filter drug-disease relations
    # ------------------------------------------------------------
    edges_dd_both_directions = df[
        df["relation"].isin(DRUG_DISEASE_RELATIONS.keys())
    ].copy()

    print(
        f"[INFO] Drug-disease relation edges before orientation filter: "
        f"{len(edges_dd_both_directions):,}"
    )

    both_direction_counts = (
        edges_dd_both_directions["relation"].value_counts().to_dict()
    )

    print(f"[INFO] Relation counts before orientation filter: {both_direction_counts}")

    # PrimeKG stores relevant edges in both directions.
    # Keep only drug -> disease orientation.
    mask = (
        (edges_dd_both_directions["x_type"] == "drug")
        & (edges_dd_both_directions["y_type"] == "disease")
    )

    edges_dd = edges_dd_both_directions[mask].copy()

    print(f"[INFO] After keeping drug->disease orientation: {len(edges_dd):,}")

    rel_counts = edges_dd["relation"].value_counts().to_dict()
    print(f"[INFO] Relation counts after orientation filter: {rel_counts}")

    # ------------------------------------------------------------
    # Build unified records
    # ------------------------------------------------------------
    rows = []
    skipped_missing_name = 0
    skipped_unknown_relation = 0

    for _, edge in edges_dd.iterrows():
        relation = clean_text(edge["relation"])
        display_relation = clean_text(edge["display_relation"])

        drug_id = clean_text(edge["x_id"])
        drug_name = clean_text(edge["x_name"])
        drug_source = clean_text(edge["x_source"])

        disease_id = clean_text(edge["y_id"])
        disease_name = clean_text(edge["y_name"])
        disease_source = clean_text(edge["y_source"])

        if relation not in DRUG_DISEASE_RELATIONS:
            skipped_unknown_relation += 1
            continue

        if not drug_name or not disease_name:
            skipped_missing_name += 1
            continue

        relationship_type, evidence_type = DRUG_DISEASE_RELATIONS[relation]

        rows.append(
            {
                "drug_name": drug_name,
                "drug_identifier": drug_id,
                "drug_identifier_type": infer_drug_identifier_type(
                    drug_source,
                    drug_id,
                ),
                "disease_or_condition_name": disease_name,
                "disease_or_condition_identifier": disease_id,
                "disease_or_condition_identifier_type": infer_disease_identifier_type(
                    disease_source,
                    disease_id,
                ),
                "relationship_type": relationship_type,
                "evidence_type": evidence_type,
                "source": "PrimeKG",
                "internal_source": "PrimeKG v2.1",
                "primekg_relation": relation,
                "primekg_display_relation": display_relation,
                "source_node_type": "drug",
                "target_node_type": "disease",
                "evidence_text": (
                    f"relation={relation}; "
                    f"display_relation={display_relation}; "
                    f"orientation=drug_to_disease"
                ),
            }
        )

    out = pd.DataFrame(rows, columns=OUTPUT_COLUMNS)

    print(f"[INFO] Standardised rows: {len(out):,}")
    print(f"[INFO] Skipped missing name: {skipped_missing_name:,}")
    print(f"[INFO] Skipped unknown relation: {skipped_unknown_relation:,}")

    if len(out) == 0:
        raise RuntimeError("No PrimeKG drug-disease rows were produced.")

    print()
    print("[OUTPUT PREVIEW]")
    print(out.head(10).to_string(index=False))

    out.to_csv(OUTPUT_FILE, index=False)

    # ------------------------------------------------------------
    # Deduplicate
    # ------------------------------------------------------------
    dedup_cols = [
        "drug_name",
        "drug_identifier",
        "disease_or_condition_name",
        "disease_or_condition_identifier",
        "relationship_type",
    ]

    dedup = out.drop_duplicates(subset=dedup_cols)
    dedup.to_csv(DEDUP_FILE, index=False)

    # ------------------------------------------------------------
    # Summary
    # ------------------------------------------------------------
    relationship_type_counts = out["relationship_type"].value_counts().to_dict()
    evidence_type_counts = out["evidence_type"].value_counts().to_dict()

    summary = {
        "timestamp_utc": now_utc(),
        "source": "PrimeKG v2.1",
        "source_repository": "https://github.com/mims-harvard/PrimeKG",
        "data_url": "https://doi.org/10.7910/DVN/IXA7BM",
        "kg_file_url": PRIMEKG_KG_URL,
        "citation_doi": "10.1038/s41597-023-01960-3",
        "code_license": "MIT",
        "interpretation": (
            "PrimeKG records are knowledge-graph drug-disease relationships. "
            "Relations are preserved as indication, contraindication, and "
            "off_label_use. The PrimeKG relation 'indication' is not relabelled "
            "as approved_indication."
        ),
        "kg_file": KG_FILE,
        "kg_file_size_bytes": int(Path(KG_FILE).stat().st_size),
        "kg_file_sha256": sha256_file(KG_FILE),
        "total_edges_loaded": int(len(df)),
        "drug_disease_edges_before_orientation_filter": int(
            len(edges_dd_both_directions)
        ),
        "relation_counts_before_orientation_filter": both_direction_counts,
        "drug_disease_edges_drug_to_disease": int(len(edges_dd)),
        "relation_counts_after_orientation_filter": rel_counts,
        "rows_written": int(len(out)),
        "deduplicated_rows_written": int(len(dedup)),
        "skipped_missing_name": int(skipped_missing_name),
        "skipped_unknown_relation": int(skipped_unknown_relation),
        "relationship_type_counts": relationship_type_counts,
        "evidence_type_counts": evidence_type_counts,
        "unique_drugs": int(out["drug_identifier"].nunique()),
        "unique_drug_names": int(out["drug_name"].nunique()),
        "unique_diseases": int(out["disease_or_condition_identifier"].nunique()),
        "unique_disease_names": int(out["disease_or_condition_name"].nunique()),
        "drug_identifier_types": (
            out["drug_identifier_type"].value_counts().to_dict()
        ),
        "disease_identifier_types": (
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