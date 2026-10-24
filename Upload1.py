#!/usr/bin/env python3
"""
Create a data dictionary from the merged drug-disease dataset only.

Input:
    data/processed/ALL_INCLUDED_SOURCES_drug_disease_merged.csv
    or
    data/processed/ALL_INCLUDED_SOURCES_drug_disease_merged.csv.gz

Outputs:
    data/processed/data_dictionary.csv
    data/processed/data_dictionary.md
    data/processed/data_dictionary_summary.json

Run:
    python Make_Data_Dictionary_MergedOnly.py

Requirements:
    pip install polars
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import polars as pl


# =============================================================================
# CONFIG
# =============================================================================

BASE = Path(
    "/data/ascher02/uqmmune1/ANNOVAR/drug_disease_data/"
    "drug-disease-mapping/data/processed"
)

MERGED_CANDIDATES = [
    BASE / "ALL_INCLUDED_SOURCES_drug_disease_merged.csv",
    BASE / "ALL_INCLUDED_SOURCES_drug_disease_merged.csv.gz",
]

OUTPUT_CSV = BASE / "data_dictionary.csv"
OUTPUT_MD = BASE / "data_dictionary.md"
OUTPUT_JSON = BASE / "data_dictionary_summary.json"


# =============================================================================
# COLUMN DEFINITIONS
# =============================================================================

COLUMN_DEFINITIONS: dict[str, dict[str, str]] = {
    # -------------------------------------------------------------------------
    # Core harmonised columns
    # -------------------------------------------------------------------------
    "source": {
        "description": "Top-level source database or resource from which the record was derived.",
        "interpretation": "Use this field to filter records by source provenance.",
        "example": "CTD",
    },
    "internal_source": {
        "description": "Source-specific subcategory, processing layer, or internal evidence origin where available.",
        "interpretation": "Useful for distinguishing subcomponents within a source.",
        "example": "CTD chemicals-diseases aggregate",
    },
    "drug_name": {
        "description": "Human-readable drug, compound, chemical, medication, or clinical intervention name.",
        "interpretation": "Primary name-level drug field used for inspection and name-based overlap analyses.",
        "example": "metformin",
    },
    "drug_identifier": {
        "description": "Source-provided or harmonised drug identifier where available.",
        "interpretation": "Identifier namespace depends on drug_identifier_type.",
        "example": "CHEMBL1431",
    },
    "drug_identifier_type": {
        "description": "Identifier namespace for the drug_identifier field.",
        "interpretation": "Examples include ChEMBL, CTD Chemical ID, RxCUI, DrugBank-like IDs, PubChem, or source-specific identifiers.",
        "example": "ChEMBL",
    },
    "disease_or_condition_name": {
        "description": "Human-readable disease, condition, indication, phenotype, or clinical trial condition name.",
        "interpretation": "Primary name-level disease or condition field used for inspection and name-based overlap analyses.",
        "example": "type 2 diabetes mellitus",
    },
    "disease_or_condition_identifier": {
        "description": "Source-provided or harmonised disease or condition identifier where available.",
        "interpretation": "Identifier namespace depends on disease_or_condition_identifier_type.",
        "example": "MONDO:0005148",
    },
    "disease_or_condition_identifier_type": {
        "description": "Identifier namespace for the disease_or_condition_identifier field.",
        "interpretation": "Examples include MeSH, MONDO, EFO, UMLS, MedDRA, OMIM, ICD, or source-specific identifiers.",
        "example": "MONDO",
    },
    "relationship_type": {
        "description": "Harmonised relationship class describing the broad drug-disease relationship meaning.",
        "interpretation": "Use this field to separate indications, contraindications, off-label uses, failed indications, clinical trial co-occurrence, chemical-disease associations, and knowledge-graph relationships.",
        "example": "indication",
    },
    "evidence_type": {
        "description": "Source-specific evidence category describing how the record should be interpreted.",
        "interpretation": "Use this field together with source and relationship_type before treating a record as therapeutic evidence.",
        "example": "curated_drug_indication",
    },

    # -------------------------------------------------------------------------
    # Merge provenance columns
    # -------------------------------------------------------------------------
    "merge_source_label": {
        "description": "Controlled source label applied during the merge process.",
        "interpretation": "Used to preserve controlled source naming during integration.",
        "example": "CTD",
    },
    "merge_internal_source_label": {
        "description": "Controlled internal source label applied during the merge process.",
        "interpretation": "Used to preserve controlled internal provenance during integration.",
        "example": "CTD chemicals-diseases aggregate",
    },
    "merge_input_file": {
        "description": "Full path of the source-level input file used during merging.",
        "interpretation": "Useful for audit trails and reproducibility checks.",
        "example": "/path/to/ctd_aggregate_chemical_disease.csv",
    },
    "merge_input_filename": {
        "description": "Filename of the source-level input file used during merging.",
        "interpretation": "Useful for identifying the processed source file.",
        "example": "ctd_aggregate_chemical_disease.csv",
    },

    # -------------------------------------------------------------------------
    # Common source-specific metadata columns
    # -------------------------------------------------------------------------
    "nct_id": {
        "description": "ClinicalTrials.gov trial identifier from AACT or trial-derived records.",
        "interpretation": "Use this field to trace clinical trial intervention-condition co-occurrence records.",
        "example": "NCT00968097",
    },
    "brief_title": {
        "description": "Brief clinical trial title, usually from AACT or ClinicalTrials.gov-derived records.",
        "interpretation": "Source-specific metadata for trial-derived records.",
        "example": "A Safety Study of SEA-TGT (SGN-TGT) in Advanced Cancer",
    },
    "official_title": {
        "description": "Official clinical trial title, usually from AACT or ClinicalTrials.gov-derived records.",
        "interpretation": "Source-specific metadata for trial-derived records.",
        "example": "A Phase 2 Clinical Study",
    },
    "overall_status": {
        "description": "Overall trial status or source-specific status label.",
        "interpretation": "For trial-derived records, this describes trial status and should not be interpreted as drug approval or efficacy.",
        "example": "RECRUITING",
    },
    "study_phase": {
        "description": "Clinical trial phase from trial-derived records.",
        "interpretation": "Use only as trial metadata; it does not imply successful efficacy.",
        "example": "PHASE3",
    },
    "clinicaltrials_url": {
        "description": "ClinicalTrials.gov URL for trial-derived records.",
        "interpretation": "Use this field to trace a record back to the original trial page.",
        "example": "https://clinicaltrials.gov/study/NCT00968097",
    },
    "source_url": {
        "description": "Source URL or download location associated with the record or source-processing script.",
        "interpretation": "Useful for provenance and reproducibility.",
        "example": "https://ctdbase.org/reports/CTD_chemicals_diseases.tsv.gz",
    },
    "chemical_id": {
        "description": "Source-specific chemical identifier, usually retained from CTD or chemical-oriented sources.",
        "interpretation": "Use together with source and drug_identifier_type to interpret the identifier namespace.",
        "example": "D000086",
    },
    "chemical_name": {
        "description": "Source-specific chemical name retained from chemical-oriented sources.",
        "interpretation": "May duplicate or complement the harmonised drug_name field.",
        "example": "Aspirin",
    },
    "cas_rn": {
        "description": "CAS Registry Number where provided by the source.",
        "interpretation": "Source-specific chemical identifier metadata.",
        "example": "50-78-2",
    },
    "direct_evidence": {
        "description": "Direct evidence label from CTD chemical-disease records where available.",
        "interpretation": "Use this field to distinguish direct curated evidence from inferred chemical-disease associations.",
        "example": "therapeutic",
    },
    "inference_gene_symbol": {
        "description": "Gene symbol supporting inferred CTD chemical-disease associations where available.",
        "interpretation": "Source-specific CTD inference metadata.",
        "example": "TNF",
    },
    "inference_score": {
        "description": "Inference score associated with CTD inferred chemical-disease records where available.",
        "interpretation": "Source-specific CTD score; interpretation should follow CTD documentation.",
        "example": "12.34",
    },
    "pubmed_ids": {
        "description": "PubMed identifiers associated with a source record where available.",
        "interpretation": "Use this field to trace literature support for source-specific records.",
        "example": "12345678|23456789",
    },
    "omim_ids": {
        "description": "OMIM identifiers associated with disease records where available.",
        "interpretation": "Source-specific disease identifier metadata.",
        "example": "125853",
    },
    "mesh_id": {
        "description": "MeSH identifier retained from source-specific drug indication or disease metadata.",
        "interpretation": "Source-specific disease or indication identifier metadata.",
        "example": "D003924",
    },
    "mesh_heading": {
        "description": "MeSH heading retained from source-specific disease metadata.",
        "interpretation": "May duplicate or complement disease_or_condition_name.",
        "example": "Diabetes Mellitus, Type 2",
    },
    "max_phase": {
        "description": "Maximum clinical phase associated with a drug indication record where available.",
        "interpretation": "Usually source-specific ChEMBL or clinical indication metadata.",
        "example": "4",
    },
    "clinical_phase": {
        "description": "Clinical phase associated with a drug or indication record where available.",
        "interpretation": "Source-specific clinical phase metadata.",
        "example": "Phase 2",
    },
    "clinical_phase_raw": {
        "description": "Original raw clinical phase label retained from the source.",
        "interpretation": "Use when exact source-level clinical phase wording is needed.",
        "example": "Phase II",
    },
    "clinical_report_ids": {
        "description": "Clinical report identifiers retained from Open Targets or related clinical indication metadata where available.",
        "interpretation": "Source-specific clinical indication metadata.",
        "example": "",
    },
    "opentargets_indication_id": {
        "description": "Open Targets indication identifier where available.",
        "interpretation": "Source-specific identifier for Open Targets clinical indication records.",
        "example": "",
    },
    "drugcentral_relationship": {
        "description": "Original DrugCentral relationship label where available.",
        "interpretation": "Use to inspect original DrugCentral semantics before harmonisation.",
        "example": "indication",
    },
    "drugcentral_struct_id": {
        "description": "DrugCentral structure identifier where available.",
        "interpretation": "Source-specific DrugCentral drug identifier metadata.",
        "example": "1234",
    },
    "source_relationship_raw": {
        "description": "Original source-specific relationship label before harmonisation.",
        "interpretation": "Use this field to audit how relationship_type was derived.",
        "example": "contraindication",
    },
    "source_table": {
        "description": "Original source table from which the record was extracted where available.",
        "interpretation": "Useful for source-level audit trails.",
        "example": "drug_indication",
    },
    "source_record_type": {
        "description": "Source-specific record category retained during processing.",
        "interpretation": "Use with source and evidence_type for source-level interpretation.",
        "example": "drug_indication",
    },
    "meddra_concept_type": {
        "description": "MedDRA concept type from SIDER records where available.",
        "interpretation": "Source-specific SIDER metadata.",
        "example": "PT",
    },
    "source_concept_id": {
        "description": "Source concept identifier retained from SIDER or other source-specific files.",
        "interpretation": "Source-specific concept metadata.",
        "example": "",
    },
    "source_concept_name": {
        "description": "Source concept name retained from SIDER or other source-specific files.",
        "interpretation": "Source-specific concept metadata.",
        "example": "",
    },
    "sider_drug_id": {
        "description": "SIDER drug identifier where available.",
        "interpretation": "Source-specific SIDER drug identifier metadata.",
        "example": "",
    },
    "extraction_method": {
        "description": "Extraction method label retained from SIDER or other source-specific records.",
        "interpretation": "Use to distinguish manual, NLP-derived, or source-specific extraction classes where available.",
        "example": "",
    },
    "rxcui": {
        "description": "RxNorm concept identifier where available.",
        "interpretation": "Usually present in MEDI-derived medication-indication records.",
        "example": "860975",
    },
    "high_precision_subset": {
        "description": "Boolean or label indicating whether a MEDI record belongs to a high-precision subset.",
        "interpretation": "Source-specific MEDI evidence quality metadata.",
        "example": "true",
    },
    "number_of_resources": {
        "description": "Number of resources supporting a MEDI medication-indication record where available.",
        "interpretation": "Source-specific MEDI support-count metadata.",
        "example": "3",
    },
    "resources_mentioning": {
        "description": "Resources contributing to or mentioning a MEDI medication-indication record.",
        "interpretation": "Source-specific MEDI provenance metadata.",
        "example": "MEDLINEPLUS|SIDER",
    },
    "status": {
        "description": "Source-specific status field, commonly used for repoDB indication status.",
        "interpretation": "Use with source to distinguish approved, failed, terminated, withdrawn, or suspended records.",
        "example": "Approved",
    },
    "phase": {
        "description": "Source-specific phase field, commonly used for repoDB or clinical indication records.",
        "interpretation": "Interpret according to the source.",
        "example": "Phase 3",
    },
    "detailed_status": {
        "description": "Detailed source-specific status label where available.",
        "interpretation": "Use with status and source for fine-grained interpretation.",
        "example": "terminated",
    },
    "metaedge": {
        "description": "Hetionet metaedge label.",
        "interpretation": "Use to distinguish compound-treats-disease and compound-palliates-disease relationships.",
        "example": "CtD",
    },
    "source_node_id": {
        "description": "Source node identifier from knowledge graph records.",
        "interpretation": "Usually represents the drug or compound node in retained drug-to-disease records.",
        "example": "Compound::DB01050",
    },
    "source_node_type": {
        "description": "Source node type from knowledge graph records.",
        "interpretation": "Used to validate drug-to-disease orientation.",
        "example": "drug",
    },
    "target_node_id": {
        "description": "Target node identifier from knowledge graph records.",
        "interpretation": "Usually represents the disease node in retained drug-to-disease records.",
        "example": "Disease::DOID:9352",
    },
    "target_node_type": {
        "description": "Target node type from knowledge graph records.",
        "interpretation": "Used to validate drug-to-disease orientation.",
        "example": "disease",
    },
    "primekg_relation": {
        "description": "Original PrimeKG relation label.",
        "interpretation": "Use to inspect original PrimeKG semantics before harmonisation.",
        "example": "indication",
    },
    "primekg_display_relation": {
        "description": "PrimeKG display relation label.",
        "interpretation": "Source-specific PrimeKG relationship display metadata.",
        "example": "indication",
    },
    "moa": {
        "description": "Mechanism of action from the Broad Drug Repurposing Hub where available.",
        "interpretation": "Source-specific compound annotation metadata.",
        "example": "AMPK activator",
    },
    "target": {
        "description": "Target annotation from the Broad Drug Repurposing Hub or similar source where available.",
        "interpretation": "Source-specific target metadata.",
        "example": "AMPK",
    },
    "disease_area": {
        "description": "Disease area annotation from the Broad Drug Repurposing Hub where available.",
        "interpretation": "Source-specific broad therapeutic area metadata.",
        "example": "metabolic disease",
    },
    "original_indication_field": {
        "description": "Original unparsed indication field from the Broad Drug Repurposing Hub or similar source.",
        "interpretation": "Useful for auditing how multiple indications were expanded.",
        "example": "diabetes; obesity",
    },
    "dataset": {
        "description": "Source-specific dataset or release label where available.",
        "interpretation": "Use with source to interpret release-specific provenance.",
        "example": "",
    },
    "release": {
        "description": "Source-specific release label or release number where available.",
        "interpretation": "Use to audit which upstream release contributed a record.",
        "example": "2025",
    },
    "concept_kind": {
        "description": "Source-specific concept kind field where available.",
        "interpretation": "Use with source-specific documentation.",
        "example": "",
    },
    "drug_id": {
        "description": "Original source-specific drug identifier retained before or alongside harmonisation.",
        "interpretation": "May duplicate or complement drug_identifier.",
        "example": "",
    },
    "disease_id": {
        "description": "Original source-specific disease identifier retained before or alongside harmonisation.",
        "interpretation": "May duplicate or complement disease_or_condition_identifier.",
        "example": "",
    },
    "disease_name": {
        "description": "Original source-specific disease name retained before or alongside harmonisation.",
        "interpretation": "May duplicate or complement disease_or_condition_name.",
        "example": "",
    },
    "relationship": {
        "description": "Original source-specific relationship label retained before or alongside harmonisation.",
        "interpretation": "May duplicate or complement relationship_type.",
        "example": "",
    },
    "evidence_text": {
        "description": "Source-specific evidence text or evidence description where available.",
        "interpretation": "Use with source-specific documentation and source provenance.",
        "example": "",
    },
    "raw_col_0": {
        "description": "Unmapped raw source column retained during source parsing.",
        "interpretation": "This column should be interpreted only by checking the original source file and source-specific processing script.",
        "example": "",
    },
    "raw_col_1": {
        "description": "Unmapped raw source column retained during source parsing.",
        "interpretation": "This column should be interpreted only by checking the original source file and source-specific processing script.",
        "example": "",
    },
    "raw_col_2": {
        "description": "Unmapped raw source column retained during source parsing.",
        "interpretation": "This column should be interpreted only by checking the original source file and source-specific processing script.",
        "example": "",
    },
    "raw_col_3": {
        "description": "Unmapped raw source column retained during source parsing.",
        "interpretation": "This column should be interpreted only by checking the original source file and source-specific processing script.",
        "example": "",
    },
    "raw_col_4": {
        "description": "Unmapped raw source column retained during source parsing.",
        "interpretation": "This column should be interpreted only by checking the original source file and source-specific processing script.",
        "example": "",
    },
    "raw_col_5": {
        "description": "Unmapped raw source column retained during source parsing.",
        "interpretation": "This column should be interpreted only by checking the original source file and source-specific processing script.",
        "example": "",
    },
    "raw_col_6": {
        "description": "Unmapped raw source column retained during source parsing.",
        "interpretation": "This column should be interpreted only by checking the original source file and source-specific processing script.",
        "example": "",
    },
}


DEFAULT_DESCRIPTION = {
    "description": "Source-specific metadata column retained from one or more input resources.",
    "interpretation": "Interpretation depends on the source and should be checked against the source-specific processing scripts and metadata.",
    "example": "",
}


# =============================================================================
# HELPERS
# =============================================================================

def now_utc() -> str:
    return datetime.now(timezone.utc).isoformat()


def find_merged_file() -> Path:
    for path in MERGED_CANDIDATES:
        if path.exists() and path.stat().st_size > 0:
            return path

    candidates = "\n".join(str(p) for p in MERGED_CANDIDATES)
    raise FileNotFoundError(f"No merged file found. Checked:\n{candidates}")


def read_sample(path: Path, n_rows: int = 50000) -> pl.DataFrame:
    return pl.read_csv(
        path,
        n_rows=n_rows,
        infer_schema_length=10000,
        ignore_errors=True,
        null_values=["", "NA", "N/A", "nan", "NaN", "None", "NULL", "null"],
    )


def get_schema(path: Path) -> dict[str, str]:
    df = read_sample(path, n_rows=1000)
    return {col: str(dtype) for col, dtype in zip(df.columns, df.dtypes)}


def get_row_count(path: Path) -> int:
    return int(
        pl.scan_csv(path, ignore_errors=True)
        .select(pl.len().alias("n"))
        .collect()
        .item()
    )


def get_non_empty_count(path: Path, column: string) -> int:
    result = (
        pl.scan_csv(path, ignore_errors=True)
        .select(
            pl.col(column)
            .cast(pl.Utf8)
            .fill_null("")
            .str.strip_chars()
            .ne("")
            .sum()
            .alias("non_empty_count")
        )
        .collect()
    )

    return int(result["non_empty_count"][0])


def get_example_values(path: Path, column: str, fallback: str = "", max_examples: int = 5) -> str:
    if fallback:
        return fallback

    try:
        df = read_sample(path, n_rows=50000)

        if column not in df.columns:
            return ""

        values = (
            df.select(pl.col(column).cast(pl.Utf8).fill_null("").str.strip_chars().alias(column))
            .filter(pl.col(column) != "")
            .unique()
            .head(max_examples)
            .to_series()
            .to_list()
        )

        return " ; ".join(str(v) for v in values if v is not None)

    except Exception:
        return ""


def markdown_escape(value: Any) -> str:
    text = "" if value is None else str(value)

    # Escape Markdown table pipe characters.
    text = text.replace("|", "\\|")

    # Keep table rows on one line.
    text = text.replace("\n", " ").replace("\r", " ")

    return text


def build_dictionary(merged_path: Path) -> tuple[pl.DataFrame, dict[str, Any]]:
    schema = get_schema(merged_path)
    row_count = get_row_count(merged_path)

    rows: list[dict[str, Any]] = []

    for column in schema.keys():
        definition = COLUMN_DEFINITIONS.get(column, DEFAULT_DESCRIPTION)

        non_empty_count = get_non_empty_count(merged_path, column)
        non_empty_percent = round((non_empty_count / row_count) * 100, 3) if row_count else 0.0
        example_values = get_example_values(
            path=merged_path,
            column=column,
            fallback=definition.get("example", ""),
            max_examples=5,
        )

        rows.append(
            {
                "column_name": column,
                "description": definition["description"],
                "interpretation": definition["interpretation"],
                "present_in_file": merged_path.name,
                "detected_dtype": schema[column],
                "example_values": example_values,
                "non_empty_count": non_empty_count,
                "total_rows": row_count,
                "non_empty_percent": non_empty_percent,
            }
        )

    dictionary_df = pl.DataFrame(rows)

    summary = {
        "created_utc": now_utc(),
        "base_directory": str(BASE),
        "input_file": str(merged_path),
        "row_count": row_count,
        "n_columns_documented": len(schema),
        "output_csv": str(OUTPUT_CSV),
        "output_markdown": str(OUTPUT_MD),
        "output_json": str(OUTPUT_JSON),
        "note": (
            "This data dictionary was generated from the merged source-level "
            "dataset only. The merged file contains the full set of harmonised "
            "columns and source-specific metadata columns retained during integration."
        ),
    }

    return dictionary_df, summary


def write_markdown(dictionary_df: pl.DataFrame, path: Path) -> None:
    lines = []
    lines.append("# Data dictionary")
    lines.append("")
    lines.append(
        "This data dictionary describes columns in the merged source-aware "
        "drug--disease relationship dataset."
    )
    lines.append("")
    lines.append("| Column | Description | Interpretation | Data type | Example values | Completeness |")
    lines.append("|---|---|---|---|---|---|")

    for row in dictionary_df.to_dicts():
        completeness = (
            f"{row['non_empty_count']:,} / {row['total_rows']:,} "
            f"({row['non_empty_percent']}%)"
        )

        lines.append(
            f"| `{markdown_escape(row['column_name'])}` | "
            f"{markdown_escape(row['description'])} | "
            f"{markdown_escape(row['interpretation'])} | "
            f"{markdown_escape(row['detected_dtype'])} | "
            f"{markdown_escape(row['example_values'])} | "
            f"{markdown_escape(completeness)} |"
        )

    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    print("=" * 100)
    print("CREATE DATA DICTIONARY FROM MERGED FILE ONLY")
    print("=" * 100)

    merged_path = find_merged_file()

    print(f"[INPUT] {merged_path}")

    dictionary_df, summary = build_dictionary(merged_path)

    OUTPUT_CSV.parent.mkdir(parents=True, exist_ok=True)

    dictionary_df.write_csv(OUTPUT_CSV)
    write_markdown(dictionary_df, OUTPUT_MD)

    with open(OUTPUT_JSON, "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    print(f"[SAVED] {OUTPUT_CSV}")
    print(f"[SAVED] {OUTPUT_MD}")
    print(f"[SAVED] {OUTPUT_JSON}")

    print("\n[SUMMARY]")
    print(json.dumps(summary, indent=2))

    print("=" * 100)
    print("DONE")
    print("=" * 100)


if __name__ == "__main__":
    main()