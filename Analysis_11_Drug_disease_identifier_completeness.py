
#!/usr/bin/env python3
"""Analysis 06: Identifier completeness and vocabulary quality (locked paper-11).

Read-only with respect to locked files. Uses source-specific assertion CSV for
identifier/mapping statistics, original event ledger for genuine source-name
versus identifier-placeholder counts, and merge manifest for reconciliation.
Name normalisation matches Analyses 01/04/08: lowercase, trim, collapse space.

Run: POLARS_MAX_THREADS=16 python Analysis_06_Identifier_and_Vocabulary_QC.py --save-tables
Dependencies: polars, tqdm
"""
from __future__ import annotations

import argparse
import csv
import json
import shutil
from pathlib import Path

import polars as pl
from tqdm import tqdm

SOURCES = {
    "AACT", "Broad Drug Repurposing Hub", "ChEMBL", "CTD", "DrugCentral",
    "Hetionet", "MEDI", "OpenTargets", "PrimeKG", "repoDB", "SIDER",
}
EXPECTED = {
    "assertions": 4_534_385,
    "events": 10_895_395,
    "drug_names": 158_364,
    "disease_names": 72_563,
}
FNS = {
    "assertions": "ALL_INCLUDED_SOURCES_drug_disease_deduplicated.csv",
    "events": "ALL_INCLUDED_SOURCES_drug_disease_merged.csv",
    "manifest": "ALL_INCLUDED_SOURCES_merge_summary.json",
}
FIELDS = [
    "source",
    "drug_name",
    "drug_identifier",
    "drug_identifier_type",
    "disease_or_condition_name",
    "disease_or_condition_identifier",
    "disease_or_condition_identifier_type",
    "supporting_event_count",
]
INVALID_ID = ["", "NA", "N/A", "NAN", "NULL", "NONE"]


def checked(path: Path, required: list[str]) -> None:
    if not path.is_file() or path.stat().st_size == 0:
        raise RuntimeError(f"Missing or empty input: {path}")
    with path.open(encoding="utf-8-sig", newline="") as f:
        cols = next(csv.reader(f, strict=True), [])
    if len(cols) != len(set(cols)) or not set(required).issubset(cols):
        raise RuntimeError(
            f"Invalid schema in {path.name}; "
            f"missing {set(required) - set(cols)}"
        )


def scan(path: Path, cols: list[str]) -> pl.LazyFrame:
    return pl.scan_csv(
        str(path),
        infer_schema=False,
        ignore_errors=False,
        truncate_ragged_lines=False,
        encoding="utf8",
    ).select(cols)


def name_key(col: str) -> pl.Expr:
    return (
        pl.col(col)
        .fill_null("")
        .str.strip_chars()
        .str.replace_all(r"\s+", " ")
        .str.to_lowercase()
    )


def cleaned(col: str) -> pl.Expr:
    return pl.col(col).fill_null("").str.strip_chars()


def is_present(col: str) -> pl.Expr:
    return ~cleaned(col).str.to_uppercase().is_in(INVALID_ID)


def require(label: str, observed, expected) -> None:
    if observed != expected:
        raise RuntimeError(
            f"[FAIL] {label}: observed {observed!r}, "
            f"expected {expected!r}"
        )
    print(
        f"[PASS] {label}: {observed:,}"
        if isinstance(observed, int)
        else f"[PASS] {label}",
        flush=True,
    )


def print_tsv(title: str, df: pl.DataFrame) -> None:
    print("\n" + "=" * 112)
    print(f"{title} [{df.height:,} rows]")
    print("=" * 112)
    print("\t".join(df.columns))
    for row in df.iter_rows():
        print("\t".join(
            "" if v is None
            else f"{v:.3f}" if isinstance(v, float)
            else str(v).replace("\t", " ").replace("\n", " ").replace("\r", " ")
            for v in row
        ))


def entity_mapping_stats(
    data: pl.LazyFrame,
    entity: str,
    top: int,
):
    if entity == "drug":
        key, id_col, ns = "_drug", "_did", "_dns"
        has_id = "_has_drug_id"
    else:
        key, id_col, ns = "_disease", "_cid", "_cns"
        has_id = "_has_disease_id"

    # Each source/name/namespace/ID membership occurs at most once.
    # Keeping the namespace avoids conflating unrelated raw IDs.
    source_map = (
        data.filter(pl.col(has_id))
        .select(["source", key, ns, id_col])
        .unique()
        .collect(engine="streaming")
    )

    global_map = source_map.drop("source").unique()

    names_to_ids = (
        global_map.group_by(key).agg([
            pl.len().alias("n_qualified_ids"),
            pl.col(ns).n_unique().alias("n_namespaces"),
        ])
    )

    names_in_ns = (
        global_map.group_by([key, ns])
        .agg(pl.len().alias("n_ids_within_namespace"))
    )

    ids_to_names = (
        global_map.group_by([ns, id_col])
        .agg(pl.col(key).n_unique().alias("n_names"))
    )

    summary = pl.DataFrame([
        {
            "metric": "names_with_identifier",
            "value": names_to_ids.height,
        },
        {
            "metric": "qualified_distinct_identifiers",
            "value": ids_to_names.height,
        },
        {
            "metric": "names_with_multiple_qualified_identifiers",
            "value": names_to_ids.filter(
                pl.col("n_qualified_ids") > 1
            ).height,
        },
        {
            "metric": "names_with_multiple_namespaces",
            "value": names_to_ids.filter(
                pl.col("n_namespaces") > 1
            ).height,
        },
        {
            "metric": "name_namespace_groups_with_multiple_ids",
            "value": names_in_ns.filter(
                pl.col("n_ids_within_namespace") > 1
            ).height,
        },
        {
            "metric": "qualified_ids_with_multiple_names",
            "value": ids_to_names.filter(
                pl.col("n_names") > 1
            ).height,
        },
        {
            "metric": "max_qualified_ids_per_name",
            "value": int(
                names_to_ids["n_qualified_ids"].max() or 0
            ),
        },
        {
            "metric": "max_names_per_qualified_id",
            "value": int(
                ids_to_names["n_names"].max() or 0
            ),
        },
    ])

    # Mapping multiplicity examples are not automatically errors.
    many_ids = (
        names_in_ns.filter(
            pl.col("n_ids_within_namespace") > 1
        )
        .sort("n_ids_within_namespace", descending=True)
        .head(top)
    )

    many_names = (
        ids_to_names.filter(pl.col("n_names") > 1)
        .sort("n_names", descending=True)
        .head(top)
    )

    per_source_namespace = (
        source_map.group_by(["source", ns])
        .agg([
            pl.len().alias(
                "distinct_source_name_id_mappings"
            ),
            pl.col(key).n_unique().alias("names"),
            pl.col(id_col).n_unique().alias("identifiers"),
        ])
        .sort(
            ["source", "distinct_source_name_id_mappings"],
            descending=[False, True],
        )
    )

    return (
        summary,
        many_ids,
        many_names,
        per_source_namespace,
    )


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--processed-dir",
        type=Path,
        default=(
            Path(__file__).resolve().parent
            / "data/processed"
        ),
    )
    p.add_argument("--save-tables", action="store_true")
    p.add_argument("--top", type=int, default=15)
    p.add_argument(
        "--drug",
        default="",
        help="Optional case-insensitive drug substring",
    )
    p.add_argument(
        "--disease",
        default="",
        help="Optional case-insensitive disease substring",
    )
    p.add_argument("--keep-work", action="store_true")

    args = p.parse_args()

    if not 1 <= args.top <= 100:
        p.error("--top must be 1..100")

    base = args.processed_dir.expanduser().resolve()
    paths = {k: base / v for k, v in FNS.items()}

    checked(paths["assertions"], FIELDS)
    checked(
        paths["events"],
        [
            "source",
            "merge_drug_name_status",
            "merge_disease_name_status",
        ],
    )

    if not paths["manifest"].is_file():
        raise RuntimeError(
            f"Missing merge manifest: {paths['manifest']}"
        )

    manifest = json.loads(
        paths["manifest"].read_text(encoding="utf-8")
    )

    require(
        "locked merge status",
        manifest.get("status"),
        "success",
    )
    require(
        "locked assertion count",
        int(manifest["assertions"]),
        EXPECTED["assertions"],
    )
    require(
        "locked event count",
        int(manifest["source_events"]),
        EXPECTED["events"],
    )

    previous_sources = {
        x["source"]: int(x["observed"])
        for x in manifest.get("sources", [])
    }

    scratch = (
        base / ".analysis06_identifier_vocabulary_work"
    )

    if scratch.exists():
        raise RuntimeError(
            f"Stale scratch directory must be inspected: {scratch}"
        )

    scratch.mkdir(parents=True)
    done = False

    try:
        with tqdm(
            total=5,
            desc="Analysis 06",
            unit="stage",
        ) as bar:

            # ====================================================
            # STEP 1: NARROW ASSERTION TABLE
            # ====================================================

            print(
                "[STEP 1] Stage narrow assertions "
                "with stable name keys",
                flush=True,
            )

            parquet = scratch / "narrow_assertions.parquet"

            query = (
                scan(paths["assertions"], FIELDS)
                .with_columns([
                    name_key("drug_name").alias("_drug"),
                    name_key(
                        "disease_or_condition_name"
                    ).alias("_disease"),
                    cleaned(
                        "drug_identifier"
                    ).alias("_did"),
                    cleaned(
                        "disease_or_condition_identifier"
                    ).alias("_cid"),
                    cleaned(
                        "drug_identifier_type"
                    ).str.to_lowercase().alias("_dns"),
                    cleaned(
                        "disease_or_condition_identifier_type"
                    ).str.to_lowercase().alias("_cns"),
                    is_present(
                        "drug_identifier"
                    ).alias("_has_drug_id"),
                    is_present(
                        "disease_or_condition_identifier"
                    ).alias("_has_disease_id"),
                    pl.col("supporting_event_count")
                    .cast(pl.Int64, strict=True)
                    .alias("_weight"),
                ])
            )

            query.sink_parquet(
                parquet,
                engine="streaming",
                compression="zstd",
            )

            data = pl.scan_parquet(parquet)
            bar.update(1)

            # ====================================================
            # STEP 2: IDENTIFIER COVERAGE
            # ====================================================

            print(
                "[STEP 2] Source completeness "
                "and vocabulary coverage",
                flush=True,
            )

            data = data.with_columns([
                (
                    pl.col("_has_drug_id")
                    & pl.col("_has_disease_id")
                ).alias("_both"),

                (
                    ~pl.col("_has_drug_id")
                    & ~pl.col("_has_disease_id")
                ).alias("_neither"),

                (
                    pl.col("_has_drug_id")
                    & (pl.col("_dns") == "")
                ).alias("_drug_ns_missing"),

                (
                    pl.col("_has_disease_id")
                    & (pl.col("_cns") == "")
                ).alias("_disease_ns_missing"),
            ])

            by_source = (
                data.group_by("source")
                .agg([
                    pl.len().alias("assertions"),

                    pl.col("_weight")
                    .sum().alias("events"),

                    pl.col("_drug")
                    .n_unique().alias("drug_names"),

                    pl.col("_disease")
                    .n_unique().alias("disease_names"),

                    (pl.col("_drug") == "")
                    .sum().alias("blank_drug_keys"),

                    (pl.col("_disease") == "")
                    .sum().alias("blank_disease_keys"),

                    pl.col("_has_drug_id")
                    .sum().alias("with_drug_id"),

                    pl.col("_has_disease_id")
                    .sum().alias("with_disease_id"),

                    pl.col("_both")
                    .sum().alias("with_both_ids"),

                    pl.col("_neither")
                    .sum().alias("with_neither_id"),

                    pl.col("_drug_ns_missing")
                    .sum().alias("drug_id_without_namespace"),

                    pl.col("_disease_ns_missing")
                    .sum().alias("disease_id_without_namespace"),

                    pl.when(pl.col("_both"))
                    .then(pl.col("_weight"))
                    .otherwise(0)
                    .sum().alias("events_with_both_ids"),

                    pl.when(pl.col("_neither"))
                    .then(pl.col("_weight"))
                    .otherwise(0)
                    .sum().alias("events_with_neither_id"),
                ])
                .collect(engine="streaming")
                .with_columns([
                    (
                        100
                        * pl.col("with_drug_id")
                        / pl.col("assertions")
                    )
                    .round(3).alias("pct_assertions_drug_id"),

                    (
                        100
                        * pl.col("with_disease_id")
                        / pl.col("assertions")
                    )
                    .round(3).alias("pct_assertions_disease_id"),

                    (
                        100
                        * pl.col("with_both_ids")
                        / pl.col("assertions")
                    )
                    .round(3).alias("pct_assertions_both_ids"),
                ])
                .sort("assertions", descending=True)
            )

            require(
                "eleven named sources",
                set(by_source["source"].to_list()),
                SOURCES,
            )
            require(
                "source-specific assertions",
                int(by_source["assertions"].sum()),
                EXPECTED["assertions"],
            )
            require(
                "source-support events",
                int(by_source["events"].sum()),
                EXPECTED["events"],
            )
            require(
                "source event counts vs manifest",
                {
                    r["source"]: int(r["events"])
                    for r in by_source.to_dicts()
                },
                previous_sources,
            )
            require(
                "blank search name keys",
                int(
                    by_source["blank_drug_keys"].sum()
                    + by_source["blank_disease_keys"].sum()
                ),
                0,
            )

            core = (
                data.select([
                    pl.col("_drug")
                    .n_unique().alias("drug_names"),

                    pl.col("_disease")
                    .n_unique().alias("disease_names"),

                    (pl.col("_weight") <= 0)
                    .sum().alias("invalid_event_weights"),
                ])
                .collect(engine="streaming")
                .row(0, named=True)
            )

            require(
                "global drug search-name vocabulary",
                int(core["drug_names"]),
                EXPECTED["drug_names"],
            )
            require(
                "global disease search-name vocabulary",
                int(core["disease_names"]),
                EXPECTED["disease_names"],
            )
            require(
                "invalid event weights",
                int(core["invalid_event_weights"]),
                0,
            )

            bar.update(1)

            # ====================================================
            # STEP 3: ORIGINAL SOURCE NAMES
            # ====================================================

            print(
                "[STEP 3] Audit original names "
                "vs identifier placeholders",
                flush=True,
            )

            provenance = (
                scan(
                    paths["events"],
                    [
                        "source",
                        "merge_drug_name_status",
                        "merge_disease_name_status",
                    ],
                )
                .group_by("source")
                .agg([
                    pl.len().alias("events_recounted"),

                    (
                        pl.col("merge_drug_name_status")
                        != "source_name"
                    )
                    .sum().alias("drug_name_placeholder_events"),

                    (
                        pl.col("merge_disease_name_status")
                        != "source_name"
                    )
                    .sum().alias("disease_name_placeholder_events"),

                    (
                        pl.col("merge_drug_name_status")
                        == "source_name"
                    )
                    .sum().alias("original_drug_name_events"),

                    (
                        pl.col("merge_disease_name_status")
                        == "source_name"
                    )
                    .sum().alias("original_disease_name_events"),
                ])
                .collect(engine="streaming")
                .sort(
                    "events_recounted",
                    descending=True,
                )
            )

            require(
                "original event-ledger row count",
                int(provenance["events_recounted"].sum()),
                EXPECTED["events"],
            )
            require(
                "event counts match original-name provenance",
                {
                    r["source"]: int(r["events_recounted"])
                    for r in provenance.to_dicts()
                },
                previous_sources,
            )

            bar.update(1)

            # ====================================================
            # STEP 4: VOCABULARY MAPPING
            # ====================================================

            print(
                "[STEP 4] Namespace-aware "
                "identifier-to-name QC",
                flush=True,
            )

            (
                drug_metrics,
                drug_many_ids,
                drug_many_names,
                drug_namespace,
            ) = entity_mapping_stats(
                data,
                "drug",
                args.top,
            )

            (
                disease_metrics,
                disease_many_ids,
                disease_many_names,
                disease_namespace,
            ) = entity_mapping_stats(
                data,
                "disease",
                args.top,
            )

            status = (
                data.with_columns(
                    pl.when(pl.col("_both"))
                    .then(pl.lit("both_ids"))
                    .when(pl.col("_has_drug_id"))
                    .then(pl.lit("drug_id_only"))
                    .when(pl.col("_has_disease_id"))
                    .then(pl.lit("disease_id_only"))
                    .otherwise(pl.lit("neither_id"))
                    .alias("id_status")
                )
                .group_by(["source", "id_status"])
                .agg([
                    pl.len().alias("assertions"),
                    pl.col("_weight").sum().alias("events"),
                ])
                .collect(engine="streaming")
                .sort(["source", "id_status"])
            )

            require(
                "ID-class assertions conserved",
                int(status["assertions"].sum()),
                EXPECTED["assertions"],
            )
            require(
                "ID-class supporting events conserved",
                int(status["events"].sum()),
                EXPECTED["events"],
            )

            bar.update(1)

            # ====================================================
            # STEP 5: FINAL SUMMARY
            # ====================================================

            print(
                "[STEP 5] Final summary",
                flush=True,
            )

            overall = pl.DataFrame([
                {
                    "metric": "source_specific_assertions",
                    "value": EXPECTED["assertions"],
                },
                {
                    "metric": "supporting_evidence_events",
                    "value": EXPECTED["events"],
                },
                {
                    "metric": "normalised_drug_name_search_keys",
                    "value": int(core["drug_names"]),
                },
                {
                    "metric": "normalised_disease_name_search_keys",
                    "value": int(core["disease_names"]),
                },
                {
                    "metric": "assertions_with_drug_id",
                    "value": int(
                        by_source["with_drug_id"].sum()
                    ),
                },
                {
                    "metric": "assertions_with_disease_id",
                    "value": int(
                        by_source["with_disease_id"].sum()
                    ),
                },
                {
                    "metric": "assertions_with_both_ids",
                    "value": int(
                        by_source["with_both_ids"].sum()
                    ),
                },
                {
                    "metric": "assertions_with_neither_id",
                    "value": int(
                        by_source["with_neither_id"].sum()
                    ),
                },
                {
                    "metric": "events_with_both_ids",
                    "value": int(
                        by_source["events_with_both_ids"].sum()
                    ),
                },
                {
                    "metric": "events_with_neither_id",
                    "value": int(
                        by_source["events_with_neither_id"].sum()
                    ),
                },
                {
                    "metric": "drug_ids_without_namespace_assertions",
                    "value": int(
                        by_source["drug_id_without_namespace"].sum()
                    ),
                },
                {
                    "metric": "disease_ids_without_namespace_assertions",
                    "value": int(
                        by_source["disease_id_without_namespace"].sum()
                    ),
                },
                {
                    "metric": "original_drug_name_events",
                    "value": int(
                        provenance["original_drug_name_events"].sum()
                    ),
                },
                {
                    "metric": "drug_name_placeholder_events",
                    "value": int(
                        provenance["drug_name_placeholder_events"].sum()
                    ),
                },
                {
                    "metric": "original_disease_name_events",
                    "value": int(
                        provenance["original_disease_name_events"].sum()
                    ),
                },
                {
                    "metric": "disease_name_placeholder_events",
                    "value": int(
                        provenance["disease_name_placeholder_events"].sum()
                    ),
                },
            ])

            bar.update(1)

        # ========================================================
        # PRINT ALL RESULTS
        # ========================================================

        print_tsv(
            "1. OVERALL IDENTIFIER AND SEARCH VOCABULARY QC",
            overall,
        )
        print_tsv(
            "2. SOURCE-LEVEL IDENTIFIER COMPLETENESS",
            by_source,
        )
        print_tsv(
            "3. IDENTIFIER STATUS PER SOURCE",
            status,
        )
        print_tsv(
            "4. ORIGINAL SOURCE NAMES VS ID PLACEHOLDERS",
            provenance,
        )
        print_tsv(
            "5. DRUG NAME / IDENTIFIER MAPPING MULTIPLICITY",
            drug_metrics,
        )
        print_tsv(
            "6. DISEASE NAME / IDENTIFIER MAPPING MULTIPLICITY",
            disease_metrics,
        )
        print_tsv(
            "7. DRUG IDENTIFIER NAMESPACES BY SOURCE",
            drug_namespace,
        )
        print_tsv(
            "8. DISEASE IDENTIFIER NAMESPACES BY SOURCE",
            disease_namespace,
        )
        print_tsv(
            "9. EXAMPLES: DRUG NAME + NAMESPACE WITH MULTIPLE IDS",
            drug_many_ids,
        )
        print_tsv(
            "10. EXAMPLES: DISEASE NAME + NAMESPACE WITH MULTIPLE IDS",
            disease_many_ids,
        )
        print_tsv(
            "11. EXAMPLES: DRUG ID + NAMESPACE WITH MULTIPLE NAMES",
            drug_many_names,
        )
        print_tsv(
            "12. EXAMPLES: DISEASE ID + NAMESPACE WITH MULTIPLE NAMES",
            disease_many_names,
        )

        # ========================================================
        # OPTIONAL NAME SEARCH
        # ========================================================

        if args.drug.strip() or args.disease.strip():
            q = data

            if args.drug.strip():
                q = q.filter(
                    pl.col("_drug").str.contains(
                        args.drug.strip().lower(),
                        literal=True,
                    )
                )

            if args.disease.strip():
                q = q.filter(
                    pl.col("_disease").str.contains(
                        args.disease.strip().lower(),
                        literal=True,
                    )
                )

            examples = (
                q.select([
                    "source",
                    "drug_name",
                    "drug_identifier",
                    "drug_identifier_type",
                    "disease_or_condition_name",
                    "disease_or_condition_identifier",
                    "disease_or_condition_identifier_type",
                    "supporting_event_count",
                ])
                .limit(args.top)
                .collect(engine="streaming")
            )

            print_tsv(
                "13. OPTIONAL DRUG / DISEASE SEARCH",
                examples,
            )

        # ========================================================
        # INTERPRETATION
        # ========================================================

        print("\n[INTERPRETATION]")
        print(
            "Name vocabulary counts reflect lowercase/"
            "whitespace-normalised labels, not unique "
            "molecules or disease concepts."
        )
        print(
            "A name present in the merge may be an identifier "
            "placeholder; original-source name coverage "
            "is reported separately."
        )
        print(
            "Namespace-qualified IDs are distinguished; "
            "the same raw ID in different namespaces "
            "is not treated as identical."
        )
        print(
            "Multiple IDs/names can reflect synonyms, "
            "source conventions, ambiguous concepts or "
            "formulation distinctions—not necessarily errors."
        )
        print(
            "Missing IDs do not invalidate searchable "
            "source-name evidence; neither ID completeness "
            "nor multiple sources proves approval."
        )

        print(
            "[PASS] ANALYSIS 06 COMPLETE; "
            "locked inputs unchanged"
        )

        # ========================================================
        # OPTIONAL OUTPUT TABLES
        # ========================================================

        if args.save_tables:
            output = (
                base.parent.parent / "tables"
            )
            output.mkdir(
                parents=True,
                exist_ok=True,
            )

            tables = {
                "overall": overall,
                "source_completeness": by_source,
                "source_id_status": status,
                "source_original_name_provenance": provenance,
                "drug_mapping_multiplicity": drug_metrics,
                "disease_mapping_multiplicity": disease_metrics,
                "drug_namespaces": drug_namespace,
                "disease_namespaces": disease_namespace,
                "drug_name_namespace_multi_ids": drug_many_ids,
                "disease_name_namespace_multi_ids": disease_many_ids,
                "drug_qualified_id_multi_names": drug_many_names,
                "disease_qualified_id_multi_names": disease_many_names,
            }

            for label, frame in tables.items():
                outfile = (
                    output / f"Analysis_06_{label}.csv"
                )
                frame.write_csv(outfile)
                print("[SAVED]", outfile)

            summary = {
                row["metric"]: int(row["value"])
                for row in overall.to_dicts()
            }

            (
                output / "Analysis_06_summary.json"
            ).write_text(
                json.dumps(
                    {
                        "snapshot": "paper11-v1",
                        "polars": pl.__version__,
                        "metrics": summary,
                    },
                    indent=2,
                ) + "\n",
                encoding="utf-8",
            )

            print(
                "[SAVED]",
                output / "Analysis_06_summary.json",
            )

        done = True

    finally:
        if done or not args.keep_work:
            shutil.rmtree(
                scratch,
                ignore_errors=True,
            )
        else:
            print(
                "[SCRATCH RETAINED]",
                scratch,
            )


if __name__ == "__main__":
    main()
