#!/usr/bin/env python3
"""
Broad Drug Repurposing Hub – drug-disease processing.

Source page:       https://repo-hub.broadinstitute.org/repurposing
Drug info file:    repurposing_drugs_20200324.txt
Citation:          Corsello SM et al. Nature Medicine. 2017;23(4):405-408.
                   doi:10.1038/nm.4306

LICENCE:
    The Broad Drug Repurposing Hub data is provided for non-commercial use only.

IMPORTANT:
    Including Broad-derived rows in a redistributed merged dataset may impose
    non-commercial restrictions on the merged work. If the final Zenodo release
    is intended to be fully open/reusable, consider distributing the script and
    summary only, not the Broad-derived CSV itself.

File format:
    - Metadata lines start with "!"
    - First non-"!" line is the real header
    - Expected columns:
        pert_iname
        clinical_phase
        moa
        target
        disease_area
        indication

Processing logic:
    - drug_name comes from pert_iname
    - disease_or_condition_name comes from indication
    - indication may contain multiple values separated by "|"
    - indication is expanded into one drug-disease row per indication
    - disease_area/category is NOT expanded
    - disease_area is preserved as metadata exactly as provided
    - clinical_phase is preserved as metadata only
    - relationship_type is kept as "drug_indication"
    - evidence_type is kept as "broad_repurposing_hub"
"""

import os
import json
import hashlib
import urllib.request
from datetime import datetime, timezone

import pandas as pd


# ============================================================
# CONFIG
# ============================================================

BROAD_DRUG_URL = (
    "https://s3.amazonaws.com/data.clue.io/repurposing/downloads/"
    "repurposing_drugs_20200324.txt"
)

RAW_DIR = (
    "/data/ascher02/uqmmune1/ANNOVAR/drug_disease_data/"
    "drug-disease-mapping/data/raw/broad_repurposing"
)

OUTPUT_DIR = (
    "/data/ascher02/uqmmune1/ANNOVAR/drug_disease_data/"
    "drug-disease-mapping/data/processed/broad_repurposing"
)

DRUG_FILE = os.path.join(RAW_DIR, "repurposing_drugs_20200324.txt")

OUTPUT_FILE = os.path.join(
    OUTPUT_DIR,
    "broad_repurposing_drug_disease.csv",
)

DEDUP_FILE = os.path.join(
    OUTPUT_DIR,
    "broad_repurposing_drug_disease_deduplicated.csv",
)

SUMMARY_FILE = os.path.join(
    OUTPUT_DIR,
    "broad_repurposing_summary.json",
)


STANDARD_COLUMNS = [
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
# HELPERS
# ============================================================

def clean_text(x) -> str:
    """Convert missing values to empty string and strip whitespace."""
    if pd.isna(x):
        return ""
    return str(x).strip()


def normalise_space(x) -> str:
    """Collapse repeated whitespace."""
    x = clean_text(x)
    return " ".join(x.split())


def split_indications(indication: str) -> list[str]:
    """
    Split Broad indication field into disease/condition terms.

    Example:
        headache|fever|toothache

    becomes:
        ["headache", "fever", "toothache"]
    """
    indication = clean_text(indication)

    if not indication:
        return []

    return [
        normalise_space(x)
        for x in indication.split("|")
        if normalise_space(x)
    ]


def sha256_file(path: str) -> str:
    """Return SHA256 checksum for a file."""
    h = hashlib.sha256()

    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)

    return h.hexdigest()


def download_if_missing(url: str, path: str) -> None:
    """Download file if missing or empty."""
    if os.path.exists(path) and os.path.getsize(path) > 0:
        print(f"[INFO] Already present: {path}")
        return

    os.makedirs(os.path.dirname(path), exist_ok=True)

    print(f"[INFO] Downloading: {url}")
    urllib.request.urlretrieve(url, path)

    size_kb = os.path.getsize(path) / 1024
    print(f"[INFO] Saved to: {path} ({size_kb:.1f} KB)")


def read_broad_file(file_path: str) -> pd.DataFrame:
    """
    Read Broad Hub TSV, skipping metadata lines beginning with '!'.
    """
    skip = 0

    with open(file_path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.startswith("!"):
                skip += 1
            else:
                break

    print(f"[INFO] Skipping {skip} metadata lines beginning with '!'")

    df = pd.read_csv(
        file_path,
        sep="\t",
        skiprows=skip,
        dtype=str,
        keep_default_na=False,
    )

    df.columns = [clean_text(c) for c in df.columns]

    return df


def validate_required_columns(df: pd.DataFrame) -> None:
    """Validate that required Broad columns exist."""
    required_cols = [
        "pert_iname",
        "clinical_phase",
        "moa",
        "target",
        "disease_area",
        "indication",
    ]

    missing = [c for c in required_cols if c not in df.columns]

    if missing:
        raise ValueError(
            "Broad drug file is missing required columns: "
            + ", ".join(missing)
        )


# ============================================================
# MAIN
# ============================================================

def main() -> None:
    print("=" * 100)
    print("BROAD DRUG REPURPOSING HUB PROCESSING")
    print("=" * 100)
    print("[SOURCE]  Broad Drug Repurposing Hub")
    print("[INPUT]   repurposing_drugs_20200324.txt")
    print("[LICENCE] Non-commercial use only")
    print("[NOTE]    indication is expanded into one drug-disease row per indication.")
    print("[NOTE]    disease_area/category is preserved as metadata and is NOT expanded.")
    print("[NOTE]    clinical_phase is preserved as metadata and is NOT reclassified.")
    print("=" * 100)

    os.makedirs(RAW_DIR, exist_ok=True)
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    download_if_missing(BROAD_DRUG_URL, DRUG_FILE)

    file_size = os.path.getsize(DRUG_FILE)
    file_sha256 = sha256_file(DRUG_FILE)

    print(f"[INFO] Input file size: {file_size:,} bytes")
    print(f"[INFO] Input SHA256:   {file_sha256}")

    df = read_broad_file(DRUG_FILE)

    print(f"[INFO] Rows loaded from Broad file: {len(df):,}")
    print(f"[INFO] Columns: {list(df.columns)}")

    validate_required_columns(df)

    rows = []

    skipped_no_drug = 0
    skipped_no_indication = 0
    expanded_rows_from_multi_indication = 0

    for _, row in df.iterrows():
        drug_name = normalise_space(row["pert_iname"])
        clinical_phase = normalise_space(row["clinical_phase"])
        moa = normalise_space(row["moa"])
        target = normalise_space(row["target"])
        disease_area = normalise_space(row["disease_area"])
        indication = normalise_space(row["indication"])

        if not drug_name:
            skipped_no_drug += 1
            continue

        if not indication:
            skipped_no_indication += 1
            continue

        indications = split_indications(indication)

        if not indications:
            skipped_no_indication += 1
            continue

        if len(indications) > 1:
            expanded_rows_from_multi_indication += len(indications)

        # IMPORTANT:
        # Keep original Broad interpretation.
        # Do NOT reclassify clinical_phase into approved/clinical/preclinical/withdrawn.
        relationship_type = "drug_indication"
        evidence_type = "broad_repurposing_hub"

        for disease_name in indications:
            evidence_parts = []

            if clinical_phase:
                evidence_parts.append(f"clinical_phase={clinical_phase}")

            if disease_area:
                evidence_parts.append(f"disease_area={disease_area}")

            if moa:
                evidence_parts.append(f"moa={moa}")

            if target:
                evidence_parts.append(f"target={target}")

            rows.append(
                {
                    "drug_name": drug_name,
                    "drug_identifier": "",
                    "drug_identifier_type": "",
                    "disease_or_condition_name": disease_name,
                    "disease_or_condition_identifier": "",
                    "disease_or_condition_identifier_type": "",
                    "relationship_type": relationship_type,
                    "evidence_type": evidence_type,
                    "source": "Broad Drug Repurposing Hub",
                    "internal_source": "Broad Repurposing Hub drug annotation file",
                    "clinical_phase": clinical_phase,
                    "disease_area": disease_area,
                    "moa": moa,
                    "target": target,
                    "original_indication_field": indication,
                    "evidence_text": "; ".join(evidence_parts),
                }
            )

    out = pd.DataFrame(rows, columns=STANDARD_COLUMNS)

    print("=" * 100)
    print("[PROCESSING SUMMARY]")
    print(f"[INFO] Rows in source file:                         {len(df):,}")
    print(f"[INFO] Rows skipped, missing drug:                 {skipped_no_drug:,}")
    print(f"[INFO] Rows skipped, missing indication:           {skipped_no_indication:,}")
    print(f"[INFO] Standardised output rows:                   {len(out):,}")
    print(
        f"[INFO] Expanded rows from multi-indication fields: "
        f"{expanded_rows_from_multi_indication:,}"
    )

    if len(out) == 0:
        print("[WARN] No output rows generated.")
    else:
        print()
        print("[EXAMPLE OUTPUT ROWS]")
        print(out.head(10).to_string(index=False))

    out.to_csv(OUTPUT_FILE, index=False)

    dedup_cols = [
        "drug_name",
        "drug_identifier",
        "disease_or_condition_name",
        "disease_or_condition_identifier",
        "relationship_type",
        "evidence_type",
        "clinical_phase",
    ]

    dedup = out.drop_duplicates(subset=dedup_cols)

    dedup.to_csv(DEDUP_FILE, index=False)

    summary = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "source": "Broad Drug Repurposing Hub",
        "source_url": "https://repo-hub.broadinstitute.org/repurposing",
        "data_url": BROAD_DRUG_URL,
        "citation_doi": "10.1038/nm.4306",
        "licence": "NON-COMMERCIAL use only",
        "licence_warning": (
            "Including Broad-derived rows in a redistributed merged dataset "
            "may impose non-commercial reuse restrictions on the merged work."
        ),
        "processing_note": (
            "The indication field is split on '|'. The disease_area field is "
            "preserved as category metadata and is not expanded. The clinical_phase "
            "field is preserved as metadata and is not used to reclassify "
            "relationship_type or evidence_type."
        ),
        "input_file": DRUG_FILE,
        "input_file_size_bytes": int(file_size),
        "input_file_sha256": file_sha256,
        "rows_in_source_file": int(len(df)),
        "rows_skipped_no_drug": int(skipped_no_drug),
        "rows_skipped_no_indication": int(skipped_no_indication),
        "rows_written": int(len(out)),
        "deduplicated_rows_written": int(len(dedup)),
        "expanded_rows_from_multi_indication": int(
            expanded_rows_from_multi_indication
        ),
        "unique_drugs": int(out["drug_name"].nunique()) if len(out) else 0,
        "unique_diseases": (
            int(out["disease_or_condition_name"].nunique())
            if len(out)
            else 0
        ),
        "unique_disease_areas": (
            int(out["disease_area"].nunique())
            if len(out)
            else 0
        ),
        "clinical_phase_counts": (
            out["clinical_phase"].value_counts().to_dict()
            if len(out)
            else {}
        ),
        "relationship_type_counts": (
            out["relationship_type"].value_counts().to_dict()
            if len(out)
            else {}
        ),
        "evidence_type_counts": (
            out["evidence_type"].value_counts().to_dict()
            if len(out)
            else {}
        ),
        "output_file": OUTPUT_FILE,
        "deduplicated_file": DEDUP_FILE,
        "summary_file": SUMMARY_FILE,
        "standard_columns": STANDARD_COLUMNS,
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