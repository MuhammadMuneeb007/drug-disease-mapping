#!/usr/bin/env python3
"""
Inspect openFDA Human Drug Label download metadata and record fields.

This does NOT extract drug-disease mappings.
It only shows:
  1. What openFDA drug/label files are available for download.
  2. What fields exist inside the JSON records.
  3. What the indication/contraindication/warning sections actually look like.
  4. Whether there are exact structured disease terms available.

Run:
  python inspect_openfda_fields.py

Optional:
  python inspect_openfda_fields.py --max-records 20
"""

from __future__ import annotations

import argparse
import json
import zipfile
from collections import Counter, defaultdict
from pathlib import Path
from urllib.parse import urlparse

import requests


OPENFDA_DOWNLOAD_JSON = "https://api.fda.gov/download.json"

RAW_DIR = Path(
    "/data/ascher02/uqmmune1/ANNOVAR/drug_disease_data/"
    "drug-disease-mapping/data/raw/openfda_drug_label"
)

INSPECT_DIR = Path(
    "/data/ascher02/uqmmune1/ANNOVAR/drug_disease_data/"
    "drug-disease-mapping/data/processed/openfda/inspection"
)


def filename_from_url(url: str, index: int) -> str:
    name = Path(urlparse(url).path).name
    if not name:
        name = f"drug-label-{index:04d}.json.zip"
    if not name.endswith(".zip"):
        name += ".zip"
    return name


def get_openfda_partitions():
    print("=" * 100)
    print("READING openFDA DOWNLOAD MANIFEST")
    print("=" * 100)

    print(f"[URL] {OPENFDA_DOWNLOAD_JSON}")

    r = requests.get(OPENFDA_DOWNLOAD_JSON, timeout=120)
    r.raise_for_status()
    metadata = r.json()

    partitions = metadata["results"]["drug"]["label"]["partitions"]

    print(f"[INFO] Number of drug/label partitions: {len(partitions)}")

    print()
    print("[FIRST 10 PARTITIONS]")
    for i, part in enumerate(partitions[:10], start=1):
        print("-" * 100)
        print(f"Partition {i}")
        print(f"  file:    {part.get('file')}")
        print(f"  records: {part.get('records')}")
        print(f"  size_mb: {part.get('size_mb')}")

    return partitions


def download_first_partition_if_missing(partitions):
    RAW_DIR.mkdir(parents=True, exist_ok=True)

    part = partitions[0]
    url = part["file"]
    zip_path = RAW_DIR / filename_from_url(url, 1)

    if zip_path.exists() and zip_path.stat().st_size > 0:
        print()
        print(f"[CACHE] Using existing first ZIP: {zip_path}")
        return zip_path

    print()
    print("=" * 100)
    print("DOWNLOADING FIRST openFDA LABEL ZIP FOR INSPECTION")
    print("=" * 100)
    print(f"[URL] {url}")
    print(f"[OUT] {zip_path}")

    with requests.get(url, stream=True, timeout=300) as r:
        r.raise_for_status()
        with open(zip_path, "wb") as f:
            for chunk in r.iter_content(chunk_size=1024 * 1024):
                if chunk:
                    f.write(chunk)

    print(f"[SAVED] {zip_path}")
    print(f"[SIZE] {zip_path.stat().st_size / 1e6:.2f} MB")

    return zip_path


def iter_records_from_zip(zip_path: Path, max_records: int):
    count = 0

    with zipfile.ZipFile(zip_path, "r") as zf:
        members = [m for m in zf.namelist() if m.endswith(".json")]

        print()
        print("=" * 100)
        print("ZIP CONTENTS")
        print("=" * 100)
        print(f"[ZIP] {zip_path}")
        print(f"[JSON MEMBERS] {members}")

        for member in members:
            with zf.open(member) as handle:
                data = json.load(handle)

            print()
            print("[TOP-LEVEL JSON KEYS]")
            print(list(data.keys()))

            print()
            print("[META]")
            print(json.dumps(data.get("meta", {}), indent=2)[:3000])

            results = data.get("results", [])

            print()
            print(f"[RESULTS COUNT IN THIS JSON] {len(results):,}")

            for record in results:
                yield record
                count += 1
                if count >= max_records:
                    return


def short_value_preview(value, max_chars=500):
    if isinstance(value, list):
        text = " ".join(str(x) for x in value[:3])
    elif isinstance(value, dict):
        text = json.dumps(value, ensure_ascii=False)[:max_chars]
    else:
        text = str(value)

    text = text.replace("\n", " ").replace("\r", " ").strip()

    if len(text) > max_chars:
        text = text[:max_chars] + "..."

    return text


def inspect_records(records):
    field_counts = Counter()
    openfda_field_counts = Counter()
    section_examples = defaultdict(list)

    all_record_keys = set()
    all_openfda_keys = set()

    sections_of_interest = [
        "indications_and_usage",
        "contraindications",
        "warnings",
        "warnings_and_precautions",
        "adverse_reactions",
        "use_in_specific_populations",
        "clinical_studies",
        "description",
        "dosage_and_administration",
    ]

    print()
    print("=" * 100)
    print("INSPECTING RECORDS")
    print("=" * 100)

    example_records = []

    for i, record in enumerate(records, start=1):
        keys = sorted(record.keys())
        all_record_keys.update(keys)

        for k in keys:
            field_counts[k] += 1

        openfda = record.get("openfda", {})
        if isinstance(openfda, dict):
            for k in openfda.keys():
                openfda_field_counts[k] += 1
                all_openfda_keys.add(k)

        for sec in sections_of_interest:
            if sec in record and len(section_examples[sec]) < 5:
                section_examples[sec].append(short_value_preview(record[sec], max_chars=700))

        if len(example_records) < 3:
            example_records.append(record)

    print()
    print("[ALL RECORD-LEVEL FIELDS FOUND]")
    for k in sorted(all_record_keys):
        print(f"  - {k}")

    print()
    print("[ALL openfda SUBFIELDS FOUND]")
    for k in sorted(all_openfda_keys):
        print(f"  - openfda.{k}")

    print()
    print("[MOST COMMON RECORD FIELDS]")
    for k, v in field_counts.most_common(50):
        print(f"  {k}: {v}")

    print()
    print("[MOST COMMON openfda FIELDS]")
    for k, v in openfda_field_counts.most_common(50):
        print(f"  openfda.{k}: {v}")

    print()
    print("=" * 100)
    print("SECTION EXAMPLES")
    print("=" * 100)

    for sec in sections_of_interest:
        print()
        print("-" * 100)
        print(f"[SECTION] {sec}")
        examples = section_examples.get(sec, [])
        if not examples:
            print("No examples found in sampled records.")
            continue

        for j, ex in enumerate(examples, start=1):
            print()
            print(f"Example {j}:")
            print(ex)

    return {
        "record_fields": sorted(all_record_keys),
        "openfda_fields": sorted(all_openfda_keys),
        "record_field_counts": dict(field_counts),
        "openfda_field_counts": dict(openfda_field_counts),
        "section_examples": dict(section_examples),
        "example_records": example_records,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-records", type=int, default=25)
    args = parser.parse_args()

    INSPECT_DIR.mkdir(parents=True, exist_ok=True)

    partitions = get_openfda_partitions()

    zip_path = download_first_partition_if_missing(partitions)

    records = list(iter_records_from_zip(zip_path, max_records=args.max_records))

    print()
    print(f"[INFO] Records sampled: {len(records)}")

    report = inspect_records(records)

    report_path = INSPECT_DIR / "openfda_field_inspection_report.json"

    with open(report_path, "w", encoding="utf-8") as f:
        json.dump(report, f, indent=2, ensure_ascii=False)

    print()
    print("=" * 100)
    print("DONE")
    print("=" * 100)
    print(f"[REPORT] {report_path}")


if __name__ == "__main__":
    main()