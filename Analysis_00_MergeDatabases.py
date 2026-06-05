#!/usr/bin/env python3
"""
Merge locked/included drug-disease association sources into one unified dataset.

Included sources
----------------
1. AACT / ClinicalTrials.gov
2. Broad Drug Repurposing Hub
3. ChEMBL
4. CTD
5. DrugCentral
6. Hetionet
7. MEDI-2
8. Open Targets
9. PrimeKG
10. repoDB
11. SIDER

Excluded sources
----------------
openFDA:
    Excluded because the current openFDA script is only an inspection script
    and does not extract final drug-disease mappings.

TTD:
    Excluded because the current TTD run produced zero rows due to invalid/
    blocked 461-byte raw input files. Do not include until the source is fixed
    and produces validated non-empty output.

Outputs
-------
data/processed/ALL_INCLUDED_SOURCES_drug_disease_merged.csv
data/processed/ALL_INCLUDED_SOURCES_drug_disease_deduplicated.csv
data/processed/ALL_INCLUDED_SOURCES_drug_disease_source_collapsed.csv
data/processed/ALL_INCLUDED_SOURCES_merge_summary.json

Run
---
python mergedanalysis.py

Requirements
------------
pip install polars
"""

from __future__ import annotations

import hashlib
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

OUTPUT_MERGED = BASE / "ALL_INCLUDED_SOURCES_drug_disease_merged.csv"
OUTPUT_DEDUP = BASE / "ALL_INCLUDED_SOURCES_drug_disease_deduplicated.csv"
OUTPUT_COLLAPSED = BASE / "ALL_INCLUDED_SOURCES_drug_disease_source_collapsed.csv"
OUTPUT_SUMMARY = BASE / "ALL_INCLUDED_SOURCES_merge_summary.json"


UNIFIED_COLS = [
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


# =============================================================================
# INCLUDED SOURCES ONLY
# =============================================================================

SOURCES = [
    {
        "source": "AACT",
        "internal_source": "ClinicalTrials.gov via AACT",
        "candidates": [
            "aact/aact_drug_disease.csv",
            "aact_drug_disease.csv",
        ],
    },
    {
        "source": "Broad Drug Repurposing Hub",
        "internal_source": "Broad Repurposing Hub drug annotation file",
        "candidates": [
            "broad_repurposing/broad_repurposing_drug_disease.csv",
            "broad_repurposing_drug_disease.csv",
        ],
    },
    {
        "source": "ChEMBL",
        "internal_source": "ChEMBL drug_indication",
        "candidates": [
            "chembl/chembl_drug_disease.csv",
            "chembl_drug_disease.csv",
            "chembl.csv",
        ],
    },
    {
        "source": "CTD",
        "internal_source": "CTD chemicals-diseases aggregate",
        "candidates": [
            "ctd/ctd_aggregate_chemical_disease.csv",
            "ctd_aggregate_chemical_disease.csv",
            "ctd_drug_disease.csv",
        ],
    },
    {
        "source": "DrugCentral",
        "internal_source": "DrugCentral clinical drug-condition relationships",
        "candidates": [
            "drugcentral/drugcentral_drug_disease.csv",
            "drugcentral_drug_disease.csv",
        ],
    },
    {
        "source": "Hetionet",
        "internal_source": "Hetionet v1.0",
        "candidates": [
            "hetionet/hetionet_drug_disease.csv",
            "hetionet_drug_disease.csv",
        ],
    },
    {
        "source": "MEDI",
        "internal_source": "MEDI-2",
        "candidates": [
            "medi/medi_drug_disease.csv",
            "medi_drug_disease.csv",
        ],
    },
    {
        "source": "OpenTargets",
        "internal_source": "Open Targets clinical_indication",
        "candidates": [
            "opentargets/opentargets_drug_disease.csv",
            "opentargets_drug_disease.csv",
        ],
    },
    {
        "source": "PrimeKG",
        "internal_source": "PrimeKG v2.1",
        "candidates": [
            "primekg/primekg_drug_disease.csv",
            "primekg_drug_disease.csv",
        ],
    },
    {
        "source": "repoDB",
        "internal_source": "repoDB full database",
        "candidates": [
            "repodb/repodb_drug_disease.csv",
            "repodb_drug_disease.csv",
        ],
    },
    {
        "source": "SIDER",
        "internal_source": "SIDER 4.1",
        "candidates": [
            "sider/sider_drug_disease.csv",
            "sider_drug_disease.csv",
            "sider/sider_drug_indication.csv",
            "sider_drug_indication.csv",
        ],
    },
]


# =============================================================================
# COLUMN MAPPING
# =============================================================================

COLUMN_ALIASES = {
    "drug_name": [
        "drug_name",
        "chemical_name",
        "compound_name",
        "pert_iname",
        "intervention_name",
        "drug",
        "name",
    ],
    "drug_identifier": [
        "drug_identifier",
        "drug_id",
        "chemical_id",
        "compound_id",
        "molecule_chembl_id",
        "chembl_id",
        "rxcui",
        "RXCUI",
        "stitch_id",
    ],
    "drug_identifier_type": [
        "drug_identifier_type",
        "chemical_identifier_type",
        "drug_id_type",
    ],
    "disease_or_condition_name": [
        "disease_or_condition_name",
        "disease_name",
        "condition_name",
        "indication_name",
        "ind_name",
        "disease",
        "condition",
    ],
    "disease_or_condition_identifier": [
        "disease_or_condition_identifier",
        "disease_identifier",
        "disease_id",
        "condition_id",
        "ind_id",
        "mesh_id",
        "meddra_id",
    ],
    "disease_or_condition_identifier_type": [
        "disease_or_condition_identifier_type",
        "disease_identifier_type",
        "condition_identifier_type",
    ],
    "relationship_type": [
        "relationship_type",
        "relationship",
        "relation",
    ],
    "evidence_type": [
        "evidence_type",
        "evidence",
        "direct_evidence",
    ],
    "source": [
        "source",
    ],
    "internal_source": [
        "internal_source",
    ],
}


DEFAULT_RELATIONSHIP = "drug_disease_association"
DEFAULT_EVIDENCE = "source_specific_drug_disease_evidence"


# =============================================================================
# HELPERS
# =============================================================================

def now_utc() -> str:
    return datetime.now(timezone.utc).isoformat()


def print_section(title: str) -> None:
    print("\n" + "=" * 100)
    print(title)
    print("=" * 100)


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


def find_existing_file(candidates: list[str]) -> Path | None:
    for rel in candidates:
        path = BASE / rel

        if path.exists() and path.stat().st_size > 0:
            return path

    return None


def read_csv_safely(path: Path) -> pl.DataFrame | None:
    if not path.exists():
        return None

    if path.stat().st_size == 0:
        print(f"[SKIP] Empty file: {path}")
        return None

    try:
        df = pl.read_csv(
            path,
            infer_schema_length=10000,
            ignore_errors=True,
            null_values=[
                "",
                "NA",
                "N/A",
                "nan",
                "NaN",
                "None",
                "NULL",
                "null",
            ],
        )

    except pl.exceptions.NoDataError:
        print(f"[SKIP] No data in file: {path}")
        return None

    except Exception as e:
        print(f"[SKIP] Could not read file: {path}")
        print(f"       Reason: {type(e).__name__}: {e}")
        return None

    if df.height == 0:
        print(f"[SKIP] Header-only or zero-row file: {path}")
        return None

    return df


def get_first_existing_col(df: pl.DataFrame, candidates: list[str]) -> pl.Expr:
    for col in candidates:
        if col in df.columns:
            return pl.col(col).cast(pl.Utf8)

    return pl.lit("", dtype=pl.Utf8)


def normalise_text_expr(col_name: str) -> pl.Expr:
    return (
        pl.col(col_name)
        .cast(pl.Utf8)
        .fill_null("")
        .str.strip_chars()
        .str.replace_all(r"\s+", " ")
        .alias(col_name)
    )


def add_missing_unified_columns(
    df: pl.DataFrame,
    source: str,
    internal_source: str,
) -> pl.DataFrame:
    """
    Map source-specific columns into the unified schema.

    Important:
    ----------
    The top-level source and internal_source columns are intentionally forced
    from the controlled SOURCES configuration.

    This prevents problems such as SIDER being counted as SIDER4.1 instead of
    SIDER in source-level summaries.

    The original input source labels are still preserved in merge metadata
    columns added later.
    """
    expressions = []

    for out_col in UNIFIED_COLS:
        # Controlled source names must always be forced.
        if out_col == "source":
            expr = pl.lit(source, dtype=pl.Utf8).alias(out_col)

        elif out_col == "internal_source":
            expr = pl.lit(internal_source, dtype=pl.Utf8).alias(out_col)

        elif out_col == "relationship_type":
            raw_expr = get_first_existing_col(
                df,
                COLUMN_ALIASES.get(out_col, [out_col]),
            )

            expr = (
                pl.when(raw_expr.is_null() | (raw_expr.str.strip_chars() == ""))
                .then(pl.lit(DEFAULT_RELATIONSHIP))
                .otherwise(raw_expr)
                .alias(out_col)
            )

        elif out_col == "evidence_type":
            raw_expr = get_first_existing_col(
                df,
                COLUMN_ALIASES.get(out_col, [out_col]),
            )

            expr = (
                pl.when(raw_expr.is_null() | (raw_expr.str.strip_chars() == ""))
                .then(pl.lit(DEFAULT_EVIDENCE))
                .otherwise(raw_expr)
                .alias(out_col)
            )

        else:
            expr = get_first_existing_col(
                df,
                COLUMN_ALIASES.get(out_col, [out_col]),
            ).alias(out_col)

        expressions.append(expr)

    return df.with_columns(expressions)


def clean_text_columns(df: pl.DataFrame) -> pl.DataFrame:
    return df.with_columns([normalise_text_expr(c) for c in df.columns])


def attach_source_metadata(
    df: pl.DataFrame,
    source: str,
    internal_source: str,
    file_path: Path,
) -> pl.DataFrame:
    """
    Preserve merge-level provenance.

    These columns make it clear which controlled source label was used,
    which file was loaded, and where the row came from.
    """
    return df.with_columns(
        [
            pl.lit(source, dtype=pl.Utf8).alias("merge_source_label"),
            pl.lit(internal_source, dtype=pl.Utf8).alias("merge_internal_source_label"),
            pl.lit(str(file_path), dtype=pl.Utf8).alias("merge_input_file"),
            pl.lit(file_path.name, dtype=pl.Utf8).alias("merge_input_filename"),
        ]
    )


def load_one_source(source_cfg: dict[str, Any]) -> tuple[pl.DataFrame | None, dict[str, Any]]:
    source = source_cfg["source"]
    internal_source = source_cfg["internal_source"]
    candidates = source_cfg["candidates"]

    source_report = {
        "source": source,
        "internal_source": internal_source,
        "candidate_paths": candidates,
        "selected_file": None,
        "status": "not_loaded",
        "rows_original": 0,
        "rows_after_core_filter": 0,
        "columns_original": [],
        "file_size_bytes": 0,
        "file_sha256": "",
        "reason": "",
    }

    file_path = find_existing_file(candidates)

    if file_path is None:
        source_report["status"] = "skipped"
        source_report["reason"] = "No candidate file found."
        print(f"[SKIP] {source}: no candidate file found.")
        return None, source_report

    source_report["selected_file"] = str(file_path)
    source_report["file_size_bytes"] = int(file_path.stat().st_size)
    source_report["file_sha256"] = sha256_file(file_path)

    df = read_csv_safely(file_path)

    if df is None:
        source_report["status"] = "skipped"
        source_report["reason"] = "File unreadable, empty, or header-only."
        return None, source_report

    source_report["rows_original"] = int(df.height)
    source_report["columns_original"] = list(df.columns)

    df = add_missing_unified_columns(
        df=df,
        source=source,
        internal_source=internal_source,
    )

    df = attach_source_metadata(
        df=df,
        source=source,
        internal_source=internal_source,
        file_path=file_path,
    )

    df = clean_text_columns(df)

    before = df.height

    df = df.filter(
        (pl.col("drug_name").is_not_null())
        & (pl.col("drug_name") != "")
        & (pl.col("disease_or_condition_name").is_not_null())
        & (pl.col("disease_or_condition_name") != "")
    )

    after = df.height

    source_report["rows_after_core_filter"] = int(after)

    if after == 0:
        source_report["status"] = "skipped"
        source_report["reason"] = (
            f"All {before:,} rows dropped because drug_name or "
            "disease_or_condition_name was missing."
        )
        print(f"[SKIP] {source}: no valid rows after core filter.")
        return None, source_report

    source_report["status"] = "loaded"

    print(
        f"[LOAD] {source}: {after:,} rows "
        f"from {file_path} "
        f"({human_size(file_path.stat().st_size)})"
    )

    return df, source_report


def count_by_column(df: pl.DataFrame, column: str) -> dict[str, int]:
    if df.height == 0 or column not in df.columns:
        return {}

    counts = (
        df.group_by(column)
        .len()
        .sort("len", descending=True)
    )

    return {
        str(row[column]): int(row["len"])
        for row in counts.to_dicts()
    }


def collapse_sources(df: pl.DataFrame) -> pl.DataFrame:
    """
    Collapse identical drug-disease pairs across sources while preserving
    source support.

    This file is useful for multi-source evidence analysis.
    """
    collapse_keys = [
        "drug_name",
        "drug_identifier",
        "drug_identifier_type",
        "disease_or_condition_name",
        "disease_or_condition_identifier",
        "disease_or_condition_identifier_type",
    ]

    collapsed = (
        df.group_by(collapse_keys)
        .agg(
            [
                pl.col("source").unique().sort().str.join("|").alias("sources"),
                pl.col("relationship_type").unique().sort().str.join("|").alias("relationship_types"),
                pl.col("evidence_type").unique().sort().str.join("|").alias("evidence_types"),
                pl.col("internal_source").unique().sort().str.join("|").alias("internal_sources"),
                pl.len().alias("row_count_before_collapse"),
            ]
        )
        .sort(["drug_name", "disease_or_condition_name"])
    )

    return collapsed


def validate_no_unexpected_source_labels(
    df: pl.DataFrame,
    expected_sources: list[str],
) -> None:
    observed_sources = sorted(
        df.select(pl.col("source").unique())
        .to_series()
        .to_list()
    )

    expected_sources_sorted = sorted(expected_sources)

    unexpected = sorted(set(observed_sources) - set(expected_sources_sorted))
    missing = sorted(set(expected_sources_sorted) - set(observed_sources))

    print_section("SOURCE LABEL VALIDATION")

    print("[EXPECTED SOURCES]")
    for s in expected_sources_sorted:
        print(f"  - {s}")

    print("\n[OBSERVED SOURCES]")
    for s in observed_sources:
        print(f"  - {s}")

    if unexpected:
        print("\n[WARNING] Unexpected source labels found:")
        for s in unexpected:
            print(f"  - {s}")
    else:
        print("\n[OK] No unexpected source labels found.")

    if missing:
        print("\n[WARNING] Expected source labels missing from merged output:")
        for s in missing:
            print(f"  - {s}")
    else:
        print("[OK] All expected loaded source labels are present.")


# =============================================================================
# MAIN
# =============================================================================

def main() -> None:
    print_section("MERGING INCLUDED DRUG-DISEASE SOURCES")
    print(f"[BASE] {BASE}")

    print("\n[INCLUDED SOURCES]")
    for cfg in SOURCES:
        print(f"  - {cfg['source']}")

    print("\n[EXCLUDED SOURCES]")
    print("  - openFDA: inspection only, no final mapping output")
    print("  - TTD: excluded until validated non-empty output is produced")

    all_frames: list[pl.DataFrame] = []
    source_reports: list[dict[str, Any]] = []

    for cfg in SOURCES:
        df, report = load_one_source(cfg)
        source_reports.append(report)

        if df is not None:
            all_frames.append(df)

    if not all_frames:
        raise RuntimeError("No included source files were loaded. Nothing to merge.")

    print_section("CONCATENATING SOURCES")

    merged = pl.concat(all_frames, how="diagonal_relaxed")
    merged = clean_text_columns(merged)

    remaining_cols = [c for c in merged.columns if c not in UNIFIED_COLS]
    merged = merged.select(UNIFIED_COLS + remaining_cols)

    print(f"[MERGED ROWS BEFORE FINAL FILTER] {merged.height:,}")

    before_final_filter = merged.height

    merged = merged.filter(
        (pl.col("drug_name").is_not_null())
        & (pl.col("drug_name") != "")
        & (pl.col("disease_or_condition_name").is_not_null())
        & (pl.col("disease_or_condition_name") != "")
    )

    rows_dropped_final_filter = before_final_filter - merged.height

    print(f"[MERGED ROWS AFTER FINAL FILTER]  {merged.height:,}")
    print(f"[ROWS DROPPED FINAL FILTER]      {rows_dropped_final_filter:,}")

    loaded_sources = [
        r["source"]
        for r in source_reports
        if r["status"] == "loaded"
    ]

    validate_no_unexpected_source_labels(
        merged,
        expected_sources=loaded_sources,
    )

    print_section("DEDUPLICATING")

    dedup_keys = [
        "drug_name",
        "drug_identifier",
        "drug_identifier_type",
        "disease_or_condition_name",
        "disease_or_condition_identifier",
        "disease_or_condition_identifier_type",
        "relationship_type",
        "evidence_type",
        "source",
    ]

    dedup = merged.unique(subset=dedup_keys, keep="first")

    print(f"[DEDUP ROWS]          {dedup.height:,}")
    print(f"[DUPLICATES REMOVED] {merged.height - dedup.height:,}")

    print_section("COLLAPSING SOURCE EVIDENCE")

    collapsed = collapse_sources(dedup)

    print(f"[COLLAPSED DRUG-DISEASE PAIRS] {collapsed.height:,}")

    print_section("WRITING OUTPUTS")

    OUTPUT_MERGED.parent.mkdir(parents=True, exist_ok=True)

    merged.write_csv(OUTPUT_MERGED)
    dedup.write_csv(OUTPUT_DEDUP)
    collapsed.write_csv(OUTPUT_COLLAPSED)

    skipped_sources = {
        r["source"]: r["reason"]
        for r in source_reports
        if r["status"] != "loaded"
    }

    merged_source_counts = count_by_column(merged, "source")
    dedup_source_counts = count_by_column(dedup, "source")
    relationship_counts = count_by_column(dedup, "relationship_type")
    evidence_counts = count_by_column(dedup, "evidence_type")

    summary = {
        "timestamp_utc": now_utc(),
        "base_directory": str(BASE),
        "included_sources_requested": [cfg["source"] for cfg in SOURCES],
        "excluded_sources": {
            "openFDA": "Inspection script only; no final drug-disease mapping output.",
            "TTD": "Excluded until valid non-empty outputs are produced.",
        },
        "loaded_sources": loaded_sources,
        "skipped_sources": skipped_sources,
        "source_reports": source_reports,
        "unified_columns": UNIFIED_COLS,
        "controlled_source_label_note": (
            "The top-level source and internal_source columns are forced from "
            "the merge configuration to maintain controlled source labels. "
            "Original file provenance is retained in merge_input_file, "
            "merge_input_filename, merge_source_label, and "
            "merge_internal_source_label."
        ),
        "deduplication_keys": dedup_keys,
        "rows_merged": int(merged.height),
        "rows_deduplicated": int(dedup.height),
        "duplicates_removed": int(merged.height - dedup.height),
        "rows_collapsed_drug_disease_pairs": int(collapsed.height),
        "rows_dropped_final_filter": int(rows_dropped_final_filter),
        "per_source_rows_merged": merged_source_counts,
        "per_source_rows_deduplicated": dedup_source_counts,
        "relationship_type_counts_deduplicated": relationship_counts,
        "evidence_type_counts_deduplicated": evidence_counts,
        "outputs": {
            "merged": str(OUTPUT_MERGED),
            "deduplicated": str(OUTPUT_DEDUP),
            "source_collapsed": str(OUTPUT_COLLAPSED),
            "summary": str(OUTPUT_SUMMARY),
        },
        "output_sha256": {
            "merged": sha256_file(OUTPUT_MERGED),
            "deduplicated": sha256_file(OUTPUT_DEDUP),
            "source_collapsed": sha256_file(OUTPUT_COLLAPSED),
        },
    }

    with open(OUTPUT_SUMMARY, "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    print(f"[SAVED] {OUTPUT_MERGED}")
    print(f"        rows={merged.height:,}")
    print(f"        sha256={sha256_file(OUTPUT_MERGED)}")

    print(f"[SAVED] {OUTPUT_DEDUP}")
    print(f"        rows={dedup.height:,}")
    print(f"        sha256={sha256_file(OUTPUT_DEDUP)}")

    print(f"[SAVED] {OUTPUT_COLLAPSED}")
    print(f"        rows={collapsed.height:,}")
    print(f"        sha256={sha256_file(OUTPUT_COLLAPSED)}")

    print(f"[SAVED] {OUTPUT_SUMMARY}")

    print_section("SUMMARY")

    print("[LOADED SOURCES]")
    for source in loaded_sources:
        n = merged_source_counts.get(source, 0)
        print(f"  - {source}: {n:,} merged rows")

    if skipped_sources:
        print("\n[SKIPPED SOURCES]")
        for source, reason in skipped_sources.items():
            print(f"  - {source}: {reason}")

    print()
    print(f"[TOTAL MERGED ROWS]       {merged.height:,}")
    print(f"[TOTAL DEDUPLICATED ROWS] {dedup.height:,}")
    print(f"[COLLAPSED PAIRS]         {collapsed.height:,}")

    print("=" * 100)
    print("DONE")
    print("=" * 100)


if __name__ == "__main__":
    main()