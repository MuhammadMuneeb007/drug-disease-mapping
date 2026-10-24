
#!/usr/bin/env python3
"""Broad Drug Repurposing Hub -> merge-compatible drug/disease CSV.

Reuses data/raw/broad_repurposing/repurposing_drugs_20200324.txt
and overwrites existing outputs in data/processed/broad_repurposing.

Does not silently download missing data.

The source annotation is not proof of regulatory approval
or treatment efficacy.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
import tempfile
import urllib.request
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd


# ============================================================
# CONFIGURATION
# ============================================================

ROOT = Path(
    os.environ.get(
        "DRUG_DISEASE_ROOT",
        Path(__file__).resolve().parent,
    )
).resolve()

RAW_DIR = ROOT / "data" / "raw" / "broad_repurposing"

OUTPUT_DIR = (
    ROOT / "data" / "processed" / "broad_repurposing"
)

RAW_NAME = "repurposing_drugs_20200324.txt"

SOURCE_URL = (
    "https://s3.amazonaws.com/data.clue.io/"
    "repurposing/downloads/" + RAW_NAME
)

SOURCE_NAME = "Broad Drug Repurposing Hub"

INTERNAL_SOURCE = (
    "Broad Repurposing Hub drug annotation file"
)

REQUIRED = [
    "pert_iname",
    "clinical_phase",
    "moa",
    "target",
    "disease_area",
    "indication",
]

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
    "clinical_phase",
    "disease_area",
    "moa",
    "target",
    "original_indication_field",
    "evidence_text",
    "broad_source_row_number",
    "indication_index",
    "n_indications_in_field",
]

COLUMNS = CORE + EXTRA

# Deduplicate identical assertions and annotations.
# Keep different targets, mechanisms and disease areas
# as distinct evidence records.
#
# The original row numbers are not part of the key.
# They are aggregated as supporting provenance.

ASSERTION_KEY = [
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
    "clinical_phase",
    "disease_area",
    "moa",
    "target",
    "original_indication_field",
    "evidence_text",
]


# ============================================================
# UTILITIES
# ============================================================

def timestamp() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256(path: Path) -> str:
    h = hashlib.sha256()

    with path.open("rb") as f:
        for chunk in iter(
            lambda: f.read(4 * 1024 * 1024),
            b"",
        ):
            h.update(chunk)

    return h.hexdigest()


def normalise_space(value: str) -> str:
    return " ".join(str(value).split())


def heading(message: str) -> None:
    print("\n" + "=" * 94)
    print(message)
    print("=" * 94)


def atomic_csv(df: pd.DataFrame, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)

    with tempfile.NamedTemporaryFile(
        mode="w",
        encoding="utf-8",
        newline="",
        suffix=".tmp",
        prefix=dest.name + ".",
        dir=dest.parent,
        delete=False,
    ) as f:
        tmp = Path(f.name)

    try:
        df.to_csv(tmp, index=False)
        os.replace(tmp, dest)

    finally:
        tmp.unlink(missing_ok=True)


def atomic_json(data: dict, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)

    with tempfile.NamedTemporaryFile(
        mode="w",
        encoding="utf-8",
        newline="",
        suffix=".tmp",
        prefix=dest.name + ".",
        dir=dest.parent,
        delete=False,
    ) as f:
        tmp = Path(f.name)

    try:
        tmp.write_text(
            json.dumps(
                data,
                indent=2,
                ensure_ascii=False,
                default=str,
            ) + "\n",
            encoding="utf-8",
        )

        os.replace(tmp, dest)

    finally:
        tmp.unlink(missing_ok=True)


# ============================================================
# RAW FILE HANDLING
# ============================================================

def obtain_raw(
    path: Path,
    allow_download: bool,
) -> None:

    if path.is_file() and path.stat().st_size > 0:
        print(f"[REUSE EXISTING RAW FILE] {path}")
        return

    if not allow_download:
        raise FileNotFoundError(
            f"Pinned raw input missing: {path}. "
            "No download was attempted. "
            "Use --download_missing only if you intend "
            "to download this snapshot."
        )

    path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    tmp = path.with_name(
        path.name + ".download"
    )

    try:
        print(
            f"[DOWNLOADING PINNED SOURCE] {SOURCE_URL}"
        )

        urllib.request.urlretrieve(
            SOURCE_URL,
            tmp,
        )

        if tmp.stat().st_size == 0:
            raise RuntimeError(
                "Downloaded input is empty"
            )

        os.replace(tmp, path)

    finally:
        tmp.unlink(missing_ok=True)


def inspect_directory(
    raw: Path,
    used: Path,
) -> list[dict]:

    heading("1. DOWNLOADED FILE INVENTORY")

    entries = []

    for path in sorted(
        p for p in raw.rglob("*") if p.is_file()
    ):
        rel = str(
            path.relative_to(raw)
        )

        active = (
            path.resolve() == used.resolve()
        )

        item = {
            "name": rel,
            "size_bytes": path.stat().st_size,
            "used_for_processing": active,
        }

        entries.append(item)

        print(
            f"  [{'INPUT' if active else 'OTHER':5}] "
            f"{rel} | "
            f"{item['size_bytes']:,} bytes"
        )

    print(
        f"Physical raw files: {len(entries):,}; "
        "files used in this processor: 1"
    )

    if not any(
        x["used_for_processing"]
        for x in entries
    ):
        raise RuntimeError(
            "Processing input does not appear in "
            "the raw-directory inventory"
        )

    return entries


# ============================================================
# STRICT SOURCE READER
# ============================================================

def read_broad(
    path: Path,
) -> tuple[list[str], list[dict], int]:
    """
    Strict tab-delimited CSV reader.

    Handles:
      - leading ! metadata lines
      - quoted tabs
      - quoted newlines
      - commas inside indication names
      - string identifiers
      - malformed record detection
    """

    with path.open(
        "r",
        encoding="utf-8-sig",
        newline="",
    ) as f:

        metadata_lines = 0

        while True:
            offset = f.tell()
            line = f.readline()

            if not line:
                raise RuntimeError(
                    "Broad input has no tabular header"
                )

            if (
                line.startswith("!")
                or not line.strip()
            ):
                metadata_lines += 1
                continue

            f.seek(offset)
            break

        reader = csv.reader(
            f,
            delimiter="\t",
            quotechar='"',
            strict=True,
        )

        header = [
            x.strip()
            for x in next(reader)
        ]

        if len(header) != len(set(header)):
            raise RuntimeError(
                "Duplicate column names in Broad header"
            )

        missing = sorted(
            set(REQUIRED) - set(header)
        )

        if missing:
            raise RuntimeError(
                f"Required columns missing: {missing}; "
                f"actual={header}"
            )

        rows = []

        for record_number, fields in enumerate(
            reader,
            start=1,
        ):

            if len(fields) != len(header):
                raise RuntimeError(
                    f"Broad raw record {record_number}: "
                    f"{len(fields)} fields, "
                    f"expected {len(header)}. "
                    "Check delimiter/quoting."
                )

            rows.append(
                {
                    "broad_source_row_number":
                        record_number,
                    **dict(zip(header, fields)),
                }
            )

    return header, rows, metadata_lines


# ============================================================
# INDICATION EXPANSION
# ============================================================

def split_indications(
    value: str,
) -> tuple[list[str], int, int]:

    raw_segments = str(value).split("|")

    clean_segments = [
        normalise_space(v)
        for v in raw_segments
    ]

    terms = [
        v for v in clean_segments
        if v
    ]

    empty = sum(
        not bool(v)
        for v in clean_segments
    )

    repeated = (
        len(terms)
        - len({
            v.casefold()
            for v in terms
        })
    )

    return terms, empty, repeated


# ============================================================
# SOURCE-LEVEL PROCESSING
# ============================================================

def process_records(
    raw_records: list[dict],
) -> tuple[pd.DataFrame, dict]:

    out = []

    missing_drug = 0
    missing_indication = 0

    missing_by_phase = Counter()

    empty_segments = 0
    repeated_terms = 0

    rows_with_multiple_indications = 0
    added_by_expansion = 0
    raw_rows_accepted = 0

    for row in raw_records:

        drug = normalise_space(
            row["pert_iname"]
        )

        phase = normalise_space(
            row["clinical_phase"]
        )

        if not drug:
            missing_drug += 1
            continue

        # Preserve the original field exactly.
        original = row["indication"]

        indications, n_empty, n_repeat = (
            split_indications(original)
        )

        empty_segments += n_empty
        repeated_terms += n_repeat

        if not indications:
            missing_indication += 1

            missing_by_phase[
                phase or "Missing"
            ] += 1

            continue

        raw_rows_accepted += 1

        if len(indications) > 1:
            rows_with_multiple_indications += 1

            added_by_expansion += (
                len(indications) - 1
            )

        area = normalise_space(
            row["disease_area"]
        )

        moa = normalise_space(
            row["moa"]
        )

        targets = normalise_space(
            row["target"]
        )

        evidence_parts = [
            f"{key}={value}"
            for key, value in (
                ("clinical_phase", phase),
                ("disease_area", area),
                ("moa", moa),
                ("target", targets),
            )
            if value
        ]

        evidence_text = "; ".join(
            evidence_parts
        )

        for index, disease in enumerate(
            indications,
            start=1,
        ):

            out.append(
                {
                    "drug_name": drug,
                    "drug_identifier": "",
                    "drug_identifier_type": "",

                    "disease_or_condition_name":
                        disease,

                    "disease_or_condition_identifier":
                        "",

                    "disease_or_condition_identifier_type":
                        "",

                    "relationship_type":
                        "drug_indication",

                    "evidence_type":
                        "broad_repurposing_hub",

                    "source": SOURCE_NAME,

                    "internal_source":
                        INTERNAL_SOURCE,

                    "clinical_phase": phase,

                    "disease_area": area,

                    "moa": moa,

                    "target": targets,

                    "original_indication_field":
                        original,

                    "evidence_text": evidence_text,

                    "broad_source_row_number":
                        str(
                            row["broad_source_row_number"]
                        ),

                    "indication_index":
                        str(index),

                    "n_indications_in_field":
                        str(len(indications)),
                }
            )

    if not out:
        raise RuntimeError(
            "No eligible drug-indication records "
            "generated. Refusing to overwrite outputs."
        )

    df = pd.DataFrame.from_records(
        out,
        columns=COLUMNS,
    ).fillna("")

    if len(df) != (
        raw_rows_accepted
        + added_by_expansion
    ):
        raise RuntimeError(
            "Indication expansion accounting error"
        )

    if len(raw_records) != (
        raw_rows_accepted
        + missing_drug
        + missing_indication
    ):
        raise RuntimeError(
            "Raw input row accounting error"
        )

    summary = {
        "rows_in_source_file":
            len(raw_records),

        "rows_skipped_no_drug":
            missing_drug,

        "rows_skipped_no_indication":
            missing_indication,

        "raw_rows_accepted":
            raw_rows_accepted,

        "rows_with_multiple_indications":
            rows_with_multiple_indications,

        "additional_rows_from_indication_expansion":
            added_by_expansion,

        "expanded_rows_from_multi_indication":
            (
                added_by_expansion
                + rows_with_multiple_indications
            ),

        "empty_pipe_segments":
            empty_segments,

        "repeated_terms_within_field":
            repeated_terms,

        "rows_missing_indication_by_clinical_phase":
            dict(missing_by_phase),

        "rows_written":
            len(df),
    }

    return df, summary


# ============================================================
# PROVENANCE-AWARE DEDUPLICATION
# ============================================================

def deduplicate_preserving_support(
    df: pd.DataFrame,
) -> pd.DataFrame:
    """
    Collapse identical assertion/annotation rows,
    retaining links to all supporting source rows.

    Differing phases, mechanisms, targets, disease
    areas and original indication fields remain distinct.
    """

    groups: dict[tuple, dict] = {}

    for item in df.to_dict(orient="records"):

        key = tuple(
            item[col]
            for col in ASSERTION_KEY
        )

        if key not in groups:
            groups[key] = {
                **item,
                "broad_supporting_record_count":
                    0,
                "broad_supporting_source_records":
                    [],
            }

        entry = groups[key]

        entry[
            "broad_supporting_record_count"
        ] += 1

        entry[
            "broad_supporting_source_records"
        ].append(
            item["broad_source_row_number"]
            + ":"
            + item["indication_index"]
        )

    rows = []

    for item in groups.values():

        item[
            "broad_supporting_source_records"
        ] = "|".join(
            item[
                "broad_supporting_source_records"
            ]
        )

        rows.append(item)

    return pd.DataFrame(
        rows,
        columns=COLUMNS + [
            "broad_supporting_record_count",
            "broad_supporting_source_records",
        ],
    )


# ============================================================
# MAIN
# ============================================================

def main() -> None:

    parser = argparse.ArgumentParser(
        description=__doc__,
    )

    parser.add_argument(
        "--input",
        default=str(
            RAW_DIR / RAW_NAME
        ),
    )

    parser.add_argument(
        "--output_dir",
        default=str(OUTPUT_DIR),
    )

    parser.add_argument(
        "--download_missing",
        action="store_true",
    )

    args = parser.parse_args()

    raw = Path(args.input).resolve()
    outdir = Path(args.output_dir).resolve()

    rawdir = raw.parent

    full_file = (
        outdir
        / "broad_repurposing_drug_disease.csv"
    )

    dedup_file = (
        outdir
        / "broad_repurposing_drug_disease_deduplicated.csv"
    )

    summary_file = (
        outdir
        / "broad_repurposing_summary.json"
    )

    # --------------------------------------------------------
    # SOURCE REUSE AND FILE INVENTORY
    # --------------------------------------------------------

    obtain_raw(
        raw,
        args.download_missing,
    )

    inventory = inspect_directory(
        rawdir,
        raw,
    )

    heading(
        "2. ACTUAL SOURCE HEADER, "
        "FORMAT AND EXAMPLE"
    )

    header, records, metadata_lines = (
        read_broad(raw)
    )

    print(f"Input: {raw}")

    print(
        "Format: TAB-SEPARATED; "
        f"skipped {metadata_lines} "
        "leading metadata/blank lines"
    )

    print(
        f"Original header ({len(header)} fields): "
        + " | ".join(header)
    )

    print(
        f"Source data records: {len(records):,}"
    )

    if records:

        print(
            "One raw example "
            "(not a generated indication row):"
        )

        for key in header:
            print(
                f"  {key:<22} "
                f"{records[0][key][:160]}"
            )

    # --------------------------------------------------------
    # EXPLAIN THE TRANSFORMATION
    # --------------------------------------------------------

    heading(
        "3. EXACT SOURCE COMBINATION "
        "AND NAMING RULES"
    )

    print(
        "No table JOIN. One source TSV is expanded "
        "vertically by the 'indication' field."
    )

    print(
        "pert_iname -> drug_name; indication split "
        "on literal '|' -> disease_or_condition_name"
    )

    print(
        "clinical_phase / moa / target / disease_area "
        "remain source annotations."
    )

    print(
        "target is never split into drug-disease rows; "
        "disease_area is never split."
    )

    print(
        "Drug and disease IDs are blank: no external "
        "ID or ontology matches are inferred."
    )

    print(
        "Source labels are kept compatible with "
        "Analysis_00_MergeDatabases.py."
    )

    print(
        "relationship_type=drug_indication; "
        "evidence_type=broad_repurposing_hub"
    )

    print(
        "These labels are NOT proof of clinical "
        "efficacy or regulatory approval."
    )

    # --------------------------------------------------------
    # PROCESS
    # --------------------------------------------------------

    full, counts = process_records(
        records
    )

    dedup = deduplicate_preserving_support(
        full
    )

    heading(
        "4. PROCESSING COUNTS AND DUPLICATES"
    )

    for key, value in counts.items():
        print(
            f"  {key:<47} {value}"
        )

    print(
        f"  {'deduplicated_assertions':<47} "
        f"{len(dedup):,}"
    )

    print(
        f"  {'redundant_source_records':<47} "
        f"{len(full) - len(dedup):,}"
    )

    print("  clinical_phase counts:")

    for phase, n in (
        full["clinical_phase"]
        .value_counts(dropna=False)
        .items()
    ):
        print(
            f"    {phase or '[empty]':<32} "
            f"{n:,}"
        )

    # --------------------------------------------------------
    # SHOW THE OUTPUT SCHEMA
    # --------------------------------------------------------

    heading(
        "5. FULL MERGE-COMPATIBLE "
        "OUTPUT SCHEMA AND SAMPLE"
    )

    print(
        f"All {len(COLUMNS)} output fields:"
    )

    for i, name in enumerate(
        COLUMNS,
        start=1,
    ):
        print(f"  {i:2d}. {name}")

    print("One actual processed row:")

    for key, value in full.iloc[0].items():
        print(
            f"  {key:<42} "
            f"{str(value)[:170]}"
        )

    # --------------------------------------------------------
    # QUALITY CHECKS
    # --------------------------------------------------------

    for col in [
        "drug_name",
        "disease_or_condition_name",
        "source",
    ]:

        if full[col].eq("").any():
            raise RuntimeError(
                f"Empty required merge values in: {col}"
            )

    if not (
        full["source"] == SOURCE_NAME
    ).all():
        raise RuntimeError(
            "Wrong controlled source label"
        )

    if len(dedup) > len(full):
        raise RuntimeError(
            "Invalid deduplication result"
        )

    # --------------------------------------------------------
    # OVERWRITE CURRENT OUTPUTS
    # --------------------------------------------------------

    heading(
        "6. OVERWRITING ORIGINAL "
        "PROCESSED OUTPUTS"
    )

    atomic_csv(
        full,
        full_file,
    )

    atomic_csv(
        dedup,
        dedup_file,
    )

    # --------------------------------------------------------
    # METADATA
    # --------------------------------------------------------

    report = {
        "status": "success",

        "timestamp_utc": timestamp(),

        "source": SOURCE_NAME,

        "source_url":
            "https://repo-hub.broadinstitute.org/repurposing",

        "data_url": SOURCE_URL,

        "dataset_snapshot": RAW_NAME,

        "citation_doi":
            "10.1038/nm.4306",

        "licence":
            "Non-commercial source data",

        "licence_warning":
            "Check source terms before redistribution "
            "of Broad-derived rows.",

        "input_file": str(raw),

        "input_size_bytes":
            raw.stat().st_size,

        "input_file_sha256":
            sha256(raw),

        "raw_file_inventory":
            inventory,

        "input_format":
            "tab-separated values, "
            "metadata lines prefixed !",

        "metadata_lines_skipped":
            metadata_lines,

        "input_header": header,

        "combination":
            "one file; split pipe-separated "
            "indication into one row per "
            "indication; no JOIN",

        "field_mapping": {
            "pert_iname": "drug_name",
            "indication":
                "disease_or_condition_name",
            "clinical_phase": "clinical_phase",
            "moa": "moa",
            "target": "target",
            "disease_area": "disease_area",
        },

        "source_semantics":
            "source drug-indication annotation, "
            "not proof of approval or efficacy",

        **counts,

        "deduplicated_rows_written":
            len(dedup),

        "deduplication_columns":
            ASSERTION_KEY,

        "deduplication_provenance":
            "source row numbers and indication "
            "indices aggregated in dedup file",

        "unique_drugs":
            int(full["drug_name"].nunique()),

        "unique_diseases":
            int(
                full[
                    "disease_or_condition_name"
                ].nunique()
            ),

        "clinical_phase_counts":
            full["clinical_phase"]
            .value_counts()
            .to_dict(),

        "relationship_type_counts":
            full["relationship_type"]
            .value_counts()
            .to_dict(),

        "evidence_type_counts":
            full["evidence_type"]
            .value_counts()
            .to_dict(),

        "standard_columns": COLUMNS,

        "output_file": str(full_file),

        "deduplicated_file": str(dedup_file),

        "summary_file": str(summary_file),

        "output_sha256":
            sha256(full_file),

        "deduplicated_sha256":
            sha256(dedup_file),
    }

    atomic_json(
        report,
        summary_file,
    )

    for output in [
        full_file,
        dedup_file,
        summary_file,
    ]:
        print(f"  [WRITTEN] {output}")

    print(
        "\nSUCCESS. Use the full non-deduplicated "
        "CSV as the source for the final merger."
    )


if __name__ == "__main__":
    main()
