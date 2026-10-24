#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Create a fresh Zenodo new-version draft and upload Scientific Data-ready v2 files.

This script:
1. Creates a new version draft from an existing published Zenodo record.
2. Creates README_DATA_RELEASE.txt.
3. Creates checksums_sha256.txt.
4. Uploads all required files individually.
5. Updates metadata.
6. Does NOT publish automatically.

IMPORTANT:
- Revoke any token you pasted publicly.
- Create a new Zenodo token and set it as an environment variable:
    export ZENODO_TOKEN="YOUR_NEW_ZENODO_TOKEN_HERE"

Run:
    python Upload_Zenodo_NewVersion.py
"""

from __future__ import annotations

import gzip
import hashlib
import json
import os
import shutil
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import requests


# =============================================================================
# CONFIG
# =============================================================================

ZENODO_API = "https://zenodo.org/api"

# From DOI https://doi.org/10.5281/zenodo.18308460
# This should be the published record ID, NOT the locked draft ID.
PUBLISHED_RECORD_ID = 18308460

BASE = Path(".").resolve()
PROCESSED = BASE / "data" / "processed"

ZENODO_TOKEN = "ExvHMhDdOkJJemSunCnvIlgw0gstmF7ShY2c3xFHcCBX38x1TbhwknawUEGs"


# Candidate paths allow files either in project root or data/processed.
def first_existing(*paths: Path) -> Path:
    for path in paths:
        if path.exists() and path.stat().st_size > 0:
            return path
    # Return first path for clearer missing-file reporting.
    return paths[0]


FILES_TO_UPLOAD = [
    first_existing(
        BASE / "ALL_INCLUDED_SOURCES_drug_disease_merged.csv.gz",
        PROCESSED / "ALL_INCLUDED_SOURCES_drug_disease_merged.csv.gz",
    ),
    first_existing(
        BASE / "ALL_INCLUDED_SOURCES_drug_disease_deduplicated.csv.gz",
        PROCESSED / "ALL_INCLUDED_SOURCES_drug_disease_deduplicated.csv.gz",
    ),
    first_existing(
        BASE / "ALL_INCLUDED_SOURCES_drug_disease_source_collapsed.csv.gz",
        PROCESSED / "ALL_INCLUDED_SOURCES_drug_disease_source_collapsed.csv.gz",
    ),
    first_existing(
        BASE / "data_dictionary.csv",
        PROCESSED / "data_dictionary.csv",
    ),
    first_existing(
        BASE / "data_dictionary_summary.json",
        PROCESSED / "data_dictionary_summary.json",
    ),
    first_existing(
        BASE / "ALL_INCLUDED_SOURCES_merge_summary.json",
        PROCESSED / "ALL_INCLUDED_SOURCES_merge_summary.json",
    ),
]

README_PATH = BASE / "README_DATA_RELEASE.txt"
CHECKSUMS_PATH = BASE / "checksums_sha256.txt"

AUTO_GENERATED_FILES = [
    README_PATH,
    CHECKSUMS_PATH,
]


METADATA = {
    "metadata": {
        "title": (
            "Source-aware drug-disease relationship dataset integrating "
            "11 public biomedical resources"
        ),
        "upload_type": "dataset",
        "description": (
            "This Zenodo release contains a source-aware drug-disease relationship "
            "dataset generated from 11 public biomedical resources: "
            "AACT/ClinicalTrials.gov, ChEMBL, Comparative Toxicogenomics Database "
            "(CTD), DrugCentral, Hetionet, MEDI, Open Targets Platform, PrimeKG, "
            "SIDER 4.1, the Broad Drug Repurposing Hub, and repoDB."
            "<br><br>"
            "The release includes the raw source-level merged dataset, the "
            "source-level schema-record deduplicated dataset, a source-collapsed "
            "representation, a data dictionary, and processing summary metadata. "
            "The raw merged dataset contains 10,870,837 source-level records. "
            "Schema-record deduplication produced 4,465,457 source-level "
            "deduplicated records. The source-collapsed representation contains "
            "4,464,947 records."
            "<br><br>"
            "Records preserve source provenance, relationship type, evidence type, "
            "drug and disease/condition names, identifiers where available, and "
            "source-specific metadata. The dataset should be interpreted as a "
            "heterogeneous, source-aware relationship resource, not as a uniformly "
            "approved indication database or a predictive drug repurposing model."
            "<br><br>"
            "The accompanying code repository is available at: "
            "https://github.com/MuhammadMuneeb007/drug-disease-mapping"
            "<br><br>"
            "Users should consider the licences, attribution requirements, and "
            "reuse terms of the original upstream databases when reusing or "
            "redistributing derived records."
        ),
        "creators": [
            {
                "name": "Muneeb, Muhammad",
                "affiliation": "The University of Queensland",
            },
            {
                "name": "Ascher, David B.",
                "affiliation": (
                    "The University of Queensland; "
                    "Baker Heart and Diabetes Institute"
                ),
            },
        ],
        "version": "v2",
        "access_right": "open",
        "license": "cc-by-4.0",
        "keywords": [
            "drug-disease relationships",
            "drug repurposing",
            "biomedical data integration",
            "public biomedical databases",
            "evidence semantics",
            "source provenance",
            "clinical trial data",
            "knowledge graphs",
            "pharmacology",
            "AACT",
            "ChEMBL",
            "CTD",
            "DrugCentral",
            "Hetionet",
            "MEDI",
            "Open Targets",
            "PrimeKG",
            "SIDER",
            "Broad Drug Repurposing Hub",
            "repoDB",
        ],
        "related_identifiers": [
            {
                "identifier": "https://github.com/MuhammadMuneeb007/drug-disease-mapping",
                "relation": "isSupplementedBy",
                "scheme": "url",
            }
        ],
    }
}


# =============================================================================
# HELPERS
# =============================================================================

def now_utc() -> str:
    return datetime.now(timezone.utc).isoformat()


def params() -> dict[str, str]:
    return {"access_token": ZENODO_TOKEN}


def check_response(response: requests.Response, step: str) -> None:
    if response.ok:
        return

    print("\n" + "=" * 100)
    print(f"FAILED: {step}")
    print(f"Status code: {response.status_code}")

    try:
        print(json.dumps(response.json(), indent=2))
    except Exception:
        print(response.text)

    print("=" * 100)
    response.raise_for_status()


def human_size(num_bytes: int) -> str:
    size = float(num_bytes)
    for unit in ["B", "KB", "MB", "GB", "TB"]:
        if size < 1024:
            return f"{size:.2f} {unit}"
        size /= 1024
    return f"{size:.2f} PB"


def sha256_file(path: Path, chunk_size: int = 1024 * 1024) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(chunk_size), b""):
            h.update(chunk)
    return h.hexdigest()


def gzip_if_needed(csv_path: Path, gz_path: Path) -> None:
    if gz_path.exists() and gz_path.stat().st_size > 0:
        return

    if not csv_path.exists() or csv_path.stat().st_size == 0:
        return

    print(f"[GZIP] Creating {gz_path.name} from {csv_path.name}")

    with csv_path.open("rb") as f_in:
        with gzip.open(gz_path, "wb", compresslevel=6) as f_out:
            shutil.copyfileobj(f_in, f_out)


def prepare_gz_files() -> None:
    gzip_if_needed(
        first_existing(
            BASE / "ALL_INCLUDED_SOURCES_drug_disease_merged.csv",
            PROCESSED / "ALL_INCLUDED_SOURCES_drug_disease_merged.csv",
        ),
        BASE / "ALL_INCLUDED_SOURCES_drug_disease_merged.csv.gz",
    )

    gzip_if_needed(
        first_existing(
            BASE / "ALL_INCLUDED_SOURCES_drug_disease_deduplicated.csv",
            PROCESSED / "ALL_INCLUDED_SOURCES_drug_disease_deduplicated.csv",
        ),
        BASE / "ALL_INCLUDED_SOURCES_drug_disease_deduplicated.csv.gz",
    )

    gzip_if_needed(
        first_existing(
            BASE / "ALL_INCLUDED_SOURCES_drug_disease_source_collapsed.csv",
            PROCESSED / "ALL_INCLUDED_SOURCES_drug_disease_source_collapsed.csv",
        ),
        BASE / "ALL_INCLUDED_SOURCES_drug_disease_source_collapsed.csv.gz",
    )


def create_readme() -> None:
    text = f"""Source-aware drug-disease relationship dataset integrating 11 public biomedical resources

Created/updated UTC:
{now_utc()}

Included sources:
1. AACT / ClinicalTrials.gov
2. ChEMBL
3. Comparative Toxicogenomics Database (CTD)
4. DrugCentral
5. Hetionet
6. MEDI
7. Open Targets Platform
8. PrimeKG
9. SIDER 4.1
10. Broad Drug Repurposing Hub
11. repoDB

Not included in the locked deposited dataset:
- TTD: source-processing script retained for inspection/future extension, but not included in the final locked dataset.
- openFDA: inspection/reference script retained for inspection/future extension, but not included in the final locked dataset.

Main files:
- ALL_INCLUDED_SOURCES_drug_disease_merged.csv.gz
  Raw source-level merged table before schema-record deduplication.
  Rows: 10,870,837

- ALL_INCLUDED_SOURCES_drug_disease_deduplicated.csv.gz
  Source-level schema-record deduplicated table.
  Rows: 4,465,457

- ALL_INCLUDED_SOURCES_drug_disease_source_collapsed.csv.gz
  Source-collapsed representation aggregating records sharing the same schema-level drug and disease/condition fields across sources.
  Rows: 4,464,947

- data_dictionary.csv
  Data dictionary generated from the merged source-level dataset.

- data_dictionary_summary.json
  Summary metadata for the data dictionary generation.

- ALL_INCLUDED_SOURCES_merge_summary.json
  Merge summary, source reports, row counts, deduplication keys, and output checksums from the integration script.

- checksums_sha256.txt
  SHA256 checksums for all deposited release files.

Recommended interpretation:
This dataset is a heterogeneous, source-aware drug-disease relationship resource. Records preserve source, relationship_type, evidence_type, names, identifiers where available, and source-specific metadata. Users should not interpret all records as approved therapeutic indications. Clinical trial records represent intervention-condition co-occurrence, CTD records represent chemical-disease associations, repoDB includes approved and failed records, and knowledge-graph sources retain graph-derived relationship semantics.

Code repository:
https://github.com/MuhammadMuneeb007/drug-disease-mapping

Licence and reuse:
Users should consider the licences, attribution requirements, and reuse terms of the original upstream databases when reusing or redistributing derived records.
"""
    README_PATH.write_text(text, encoding="utf-8")
    print(f"[SAVED] {README_PATH}")


def create_checksums(files: list[Path]) -> None:
    lines = []
    for path in files:
        if path.exists() and path.stat().st_size > 0:
            lines.append(f"{sha256_file(path)}  {path.name}")
    CHECKSUMS_PATH.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"[SAVED] {CHECKSUMS_PATH}")


def check_inputs() -> list[Path]:
    if not ZENODO_TOKEN:
        raise RuntimeError(
            "ZENODO_TOKEN is missing.\n\n"
            "Create a new Zenodo token, then run:\n"
            "export ZENODO_TOKEN='YOUR_NEW_ZENODO_TOKEN_HERE'\n"
        )

    prepare_gz_files()
    create_readme()

    base_files = FILES_TO_UPLOAD
    create_checksums(base_files)

    upload_files = base_files + AUTO_GENERATED_FILES

    print("\n[CHECKING FILES]")
    missing = []

    for path in upload_files:
        if not path.exists() or path.stat().st_size == 0:
            missing.append(path)
            continue
        print(f"[OK] {path.name} ({human_size(path.stat().st_size)})")

    if missing:
        print("\n[MISSING FILES]")
        for path in missing:
            print(f" - {path}")
        raise FileNotFoundError("Some required files are missing.")

    return upload_files


def create_new_version_draft() -> tuple[int, str]:
    print("\n[1] Creating new version draft from published record...")

    url = f"{ZENODO_API}/deposit/depositions/{PUBLISHED_RECORD_ID}/actions/newversion"

    response = requests.post(url, params=params())
    check_response(response, "create new version draft")

    data = response.json()
    latest_draft_url = data["links"]["latest_draft"]

    print(f"[OK] New version action accepted.")
    print(f"[INFO] Latest draft URL: {latest_draft_url}")

    response = requests.get(latest_draft_url, params=params())
    check_response(response, "get latest draft")

    draft = response.json()
    draft_id = int(draft["id"])
    bucket_url = draft["links"]["bucket"]

    print(f"[OK] New draft ID: {draft_id}")
    print("[OK] Bucket URL found.")

    return draft_id, bucket_url


def list_existing_files(draft_id: int) -> list[dict[str, Any]]:
    print("\n[2] Listing files in new version draft...")

    url = f"{ZENODO_API}/deposit/depositions/{draft_id}/files"
    response = requests.get(url, params=params())
    check_response(response, "list files")

    files = response.json()

    if not files:
        print("[OK] No existing files copied into draft.")
    else:
        print("[INFO] Files copied into draft:")
        for f in files:
            print(f" - {f.get('filename')}")

    return files


def delete_all_existing_files(draft_id: int, files: list[dict[str, Any]]) -> None:
    print("\n[3] Removing copied old files from new version draft...")

    if not files:
        print("[OK] No files to remove.")
        return

    files_url = f"{ZENODO_API}/deposit/depositions/{draft_id}/files"

    for f in files:
        filename = f.get("filename")
        file_id = f.get("id")

        print(f"[DELETE] {filename}")

        response = requests.delete(f"{files_url}/{file_id}", params=params())

        if response.status_code == 403 and "Bucket is locked" in response.text:
            raise RuntimeError(
                "The new version draft bucket is locked. Open the Zenodo dashboard, "
                "discard/delete the stuck draft if possible, then create a fresh new "
                "version manually or contact Zenodo support."
            )

        check_response(response, f"delete {filename}")

    print("[OK] Old copied files removed.")


def upload_file(bucket_url: str, path: Path) -> None:
    print(f"\n[UPLOAD] {path.name}")
    print(f"         Size: {human_size(path.stat().st_size)}")

    upload_url = f"{bucket_url}/{path.name}"

    headers = {
        "Content-Type": "application/octet-stream",
        "Content-Length": str(path.stat().st_size),
    }

    with path.open("rb") as f:
        response = requests.put(
            upload_url,
            data=f,
            headers=headers,
            params=params(),
            timeout=None,
        )

    check_response(response, f"upload {path.name}")

    print(f"[OK] Uploaded: {path.name}")


def upload_all_files(bucket_url: str, upload_files: list[Path]) -> None:
    print("\n[4] Uploading v2 files...")

    for path in upload_files:
        upload_file(bucket_url, path)

    print("[OK] All files uploaded.")


def update_metadata(draft_id: int) -> None:
    print("\n[5] Updating metadata...")

    url = f"{ZENODO_API}/deposit/depositions/{draft_id}"

    response = requests.put(
        url,
        params=params(),
        headers={"Content-Type": "application/json"},
        data=json.dumps(METADATA),
    )

    check_response(response, "update metadata")

    print("[OK] Metadata updated.")


def print_final_file_list(draft_id: int) -> None:
    print("\n[6] Verifying uploaded draft file list...")

    files = list_existing_files(draft_id)

    print("\n[FINAL FILES IN NEW DRAFT]")
    for f in files:
        print(f" - {f.get('filename')}")


def main() -> None:
    print("=" * 100)
    print("ZENODO NEW VERSION UPLOAD FOR SCIENTIFIC DATA SUBMISSION")
    print("=" * 100)
    print(f"[BASE] {BASE}")
    print(f"[PUBLISHED RECORD ID] {PUBLISHED_RECORD_ID}")

    upload_files = check_inputs()

    draft_id, bucket_url = create_new_version_draft()

    copied_files = list_existing_files(draft_id)
    delete_all_existing_files(draft_id, copied_files)

    upload_all_files(bucket_url, upload_files)
    update_metadata(draft_id)
    print_final_file_list(draft_id)

    print("\n" + "=" * 100)
    print("[DONE] Fresh Zenodo v2 draft created and updated.")
    print(f"New draft ID: {draft_id}")
    print("Open the Zenodo dashboard, review metadata/files carefully, then publish manually.")
    print("=" * 100)


if __name__ == "__main__":
    main()