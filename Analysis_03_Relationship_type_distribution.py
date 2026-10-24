
#!/usr/bin/env python3
"""Analysis 03: original relationship types and analytical evidence families.

Read-only analysis of the locked 11-source merged dataset.
Run: POLARS_MAX_THREADS=16 python Analysis_03_Relationship_type_distribution.py
Optional: --save-tables (writes only new summary tables under tables/)
Requires: polars, tqdm
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import polars as pl
from tqdm import tqdm

EXPECTED = {
    "events": 10_895_395,
    "assertions": 4_534_385,
    "pairs": 4_470_601,
    "conflicts": 87,
}

SOURCES = {
    "AACT",
    "Broad Drug Repurposing Hub",
    "ChEMBL",
    "CTD",
    "DrugCentral",
    "Hetionet",
    "MEDI",
    "OpenTargets",
    "PrimeKG",
    "repoDB",
    "SIDER",
}

FILES = {
    "events": "ALL_INCLUDED_SOURCES_drug_disease_merged.csv",
    "assertions": "ALL_INCLUDED_SOURCES_drug_disease_deduplicated.csv",
    "pairs": "ALL_INCLUDED_SOURCES_drug_disease_source_collapsed.csv",
    "manifest": "ALL_INCLUDED_SOURCES_merge_summary.json",
}


def verify_columns(path: Path, required: list[str]) -> None:
    if not path.is_file() or not path.stat().st_size:
        raise RuntimeError(f"Missing/empty input: {path}")

    with path.open(
        "r", encoding="utf-8-sig", newline=""
    ) as handle:
        names = next(csv.reader(handle, strict=True), [])

    if (
        len(names) != len(set(names))
        or not set(required).issubset(names)
    ):
        raise RuntimeError(
            f"Unexpected header in {path.name}; required={required}"
        )


def lazy_csv(path: Path, cols: list[str]) -> pl.LazyFrame:
    return pl.scan_csv(
        path,
        infer_schema=False,
        ignore_errors=False,
        truncate_ragged_lines=False,
        encoding="utf8",
    ).select(cols)


def require_equal(description: str, actual, expected) -> None:
    if actual != expected:
        raise RuntimeError(
            f"[FAIL] {description}: {actual!r} != {expected!r}"
        )

    if isinstance(actual, int):
        print(
            f"[PASS] {description}: {actual:,}",
            flush=True,
        )
    else:
        print(f"[PASS] {description}", flush=True)


def print_tsv(title: str, df: pl.DataFrame) -> None:
    """Print every row without Polars terminal truncation."""
    print("\n" + "=" * 105)
    print(title + f"  [rows={df.height:,}]")
    print("=" * 105)
    print("\t".join(df.columns))

    for row in df.iter_rows():
        values = []

        for v in row:
            if v is None:
                values.append("")
            elif isinstance(v, float):
                values.append(f"{v:.3f}")
            else:
                values.append(
                    str(v)
                    .replace("\n", " ")
                    .replace("\t", " ")
                )

        print("\t".join(values))


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__
    )

    parser.add_argument(
        "--processed-dir",
        type=Path,
        default=(
            Path(__file__).resolve().parent
            / "data"
            / "processed"
        ),
    )

    parser.add_argument(
        "--save-tables",
        action="store_true",
    )

    parser.add_argument(
        "--examples",
        type=int,
        default=10,
        help=(
            "Number of flagged potential conflict "
            "examples to print"
        ),
    )

    args = parser.parse_args()

    if not 0 <= args.examples <= 100:
        parser.error(
            "--examples must be between 0 and 100"
        )

    base = (
        args.processed_dir
        .expanduser()
        .resolve()
    )

    paths = {
        k: base / v
        for k, v in FILES.items()
    }

    verify_columns(
        paths["events"],
        [
            "source",
            "relationship_type",
            "merge_evidence_family",
        ],
    )

    verify_columns(
        paths["assertions"],
        [
            "source",
            "relationship_type",
            "supporting_event_count",
        ],
    )

    verify_columns(
        paths["pairs"],
        [
            "supporting_event_count",
            "supporting_assertion_count",
            "potential_conflict",
            "indication_like_event_count",
            "contraindication_event_count",
        ],
    )

    if not paths["manifest"].is_file():
        raise RuntimeError(
            f"Missing merge summary: {paths['manifest']}"
        )

    manifest = json.loads(
        paths["manifest"].read_text(
            encoding="utf-8"
        )
    )

    require_equal(
        "Merge manifest status",
        manifest.get("status"),
        "success",
    )

    print(
        "\nANALYSIS 03 — ORIGINAL RELATIONSHIP TYPES "
        "| LOCKED 11-SOURCE DATA"
    )

    print(
        "Original relationship_type labels "
        "are NOT approval claims."
    )

    print(
        "Event totals and assertion totals "
        "use separate denominators.\n"
    )

    with tqdm(
        total=4,
        desc="Analysis 03",
        unit="stage",
    ) as progress:

        # ----------------------------------------------------
        # STEP 1: ORIGINAL EVENT RELATIONSHIP COUNTS
        # ----------------------------------------------------

        crosswalk = (
            lazy_csv(
                paths["events"],
                [
                    "source",
                    "relationship_type",
                    "merge_evidence_family",
                ],
            )
            .group_by(
                [
                    "source",
                    "relationship_type",
                    "merge_evidence_family",
                ]
            )
            .agg(
                pl.len().alias("events")
            )
            .collect(engine="streaming")
        )

        progress.update(1)

        # ----------------------------------------------------
        # STEP 2: DISTINCT SOURCE ASSERTIONS
        # ----------------------------------------------------

        assertions = (
            lazy_csv(
                paths["assertions"],
                [
                    "source",
                    "relationship_type",
                    "supporting_event_count",
                ],
            )
            .with_columns(
                pl.col("supporting_event_count")
                .cast(pl.Int64, strict=True)
            )
            .group_by(
                [
                    "source",
                    "relationship_type",
                ]
            )
            .agg(
                [
                    pl.len().alias("assertions"),

                    pl.col("supporting_event_count")
                    .sum()
                    .alias("linked_events"),
                ]
            )
            .collect(engine="streaming")
        )

        progress.update(1)

        # ----------------------------------------------------
        # STEP 3: VALIDATE EXACT-IDENTITY PAIR COUNTS
        # ----------------------------------------------------

        pair_base = (
            lazy_csv(
                paths["pairs"],
                [
                    "supporting_event_count",
                    "supporting_assertion_count",
                    "potential_conflict",
                    "indication_like_event_count",
                    "contraindication_event_count",
                ],
            )
            .with_columns(
                [
                    pl.col(c)
                    .cast(pl.Int64, strict=True)
                    for c in [
                        "supporting_event_count",
                        "supporting_assertion_count",
                        "indication_like_event_count",
                        "contraindication_event_count",
                    ]
                ]
            )
            .with_columns(
                [
                    (
                        pl.col("potential_conflict")
                        .fill_null("")
                        .str.to_lowercase()
                        == "true"
                    ).alias("_flagged"),

                    (
                        (
                            pl.col(
                                "indication_like_event_count"
                            ) > 0
                        )
                        &
                        (
                            pl.col(
                                "contraindication_event_count"
                            ) > 0
                        )
                    ).alias("_recomputed"),
                ]
            )
        )

        pair_stats = (
            pair_base
            .select(
                [
                    pl.len().alias("pairs"),

                    pl.col("supporting_event_count")
                    .sum()
                    .alias("events"),

                    pl.col("supporting_assertion_count")
                    .sum()
                    .alias("assertions"),

                    pl.col("_flagged")
                    .sum()
                    .alias("conflicts"),

                    (
                        pl.col("_flagged")
                        != pl.col("_recomputed")
                    )
                    .sum()
                    .alias("mismatched_conflict_flags"),
                ]
            )
            .collect(engine="streaming")
            .row(0, named=True)
        )

        progress.update(1)

        # ----------------------------------------------------
        # STEP 4: RECONCILE AND GENERATE SUMMARY TABLES
        # ----------------------------------------------------

        relation_source_events = (
            crosswalk
            .group_by(
                [
                    "source",
                    "relationship_type",
                ]
            )
            .agg(
                pl.col("events")
                .sum()
                .alias("events")
            )
        )

        by_source = (
            relation_source_events
            .join(
                assertions,
                on=[
                    "source",
                    "relationship_type",
                ],
                how="full",
                coalesce=True,
            )
        )

        invalid = by_source.filter(
            pl.col("events").is_null()
            | pl.col("assertions").is_null()
            | pl.col("linked_events").is_null()
            | (
                pl.col("events")
                != pl.col("linked_events")
            )
            | (
                pl.col("assertions")
                > pl.col("events")
            )
            | (
                pl.col("assertions") <= 0
            )
        )

        if invalid.height:
            raise RuntimeError(
                "Relationship-level event/assertion "
                "support mismatch"
            )

        for c in [
            "source",
            "relationship_type",
            "merge_evidence_family",
        ]:
            bad = crosswalk.filter(
                pl.col(c).is_null()
                |
                (
                    pl.col(c)
                    .str.strip_chars()
                    == ""
                )
            )

            if bad.height:
                raise RuntimeError(
                    f"Missing original or analytical label: {c}"
                )

        sources_seen = set(
            by_source["source"].to_list()
        )

        require_equal(
            "All 11 source labels",
            sources_seen,
            SOURCES,
        )

        events_n = int(
            by_source["events"].sum()
        )

        assertions_n = int(
            by_source["assertions"].sum()
        )

        require_equal(
            "Evidence events",
            events_n,
            EXPECTED["events"],
        )

        require_equal(
            "Distinct source-specific assertions",
            assertions_n,
            EXPECTED["assertions"],
        )

        require_equal(
            "Event counts from pair file",
            int(pair_stats["events"]),
            events_n,
        )

        require_equal(
            "Assertion counts from pair file",
            int(pair_stats["assertions"]),
            assertions_n,
        )

        require_equal(
            "Exact-identity pair count",
            int(pair_stats["pairs"]),
            EXPECTED["pairs"],
        )

        require_equal(
            "Potential conflict count",
            int(pair_stats["conflicts"]),
            EXPECTED["conflicts"],
        )

        require_equal(
            "Pair conflict flags internally consistent",
            int(pair_stats["mismatched_conflict_flags"]),
            0,
        )

        require_equal(
            "Manifest event count",
            int(manifest["source_events"]),
            events_n,
        )

        require_equal(
            "Manifest assertion count",
            int(manifest["assertions"]),
            assertions_n,
        )

        require_equal(
            "Manifest pair count",
            int(manifest["exact_identity_pairs"]),
            int(pair_stats["pairs"]),
        )

        require_equal(
            "Manifest conflict count",
            int(manifest["potential_conflicting_pairs"]),
            int(pair_stats["conflicts"]),
        )

        require_equal(
            "Family crosswalk conserves events",
            int(crosswalk["events"].sum()),
            events_n,
        )

        # ----------------------------------------------------
        # PER-SOURCE COMPOSITION
        # ----------------------------------------------------

        summary = (
            by_source
            .group_by("source")
            .agg(
                [
                    pl.col("events")
                    .sum()
                    .alias("events"),

                    pl.col("assertions")
                    .sum()
                    .alias("assertions"),

                    pl.col("relationship_type")
                    .n_unique()
                    .alias("distinct_relationship_types"),
                ]
            )
            .with_columns(
                [
                    (
                        pl.col("events")
                        - pl.col("assertions")
                    ).alias("additional_support_events"),

                    (
                        100
                        * pl.col("events")
                        / events_n
                    )
                    .round(3)
                    .alias("pct_all_events"),
                ]
            )
            .sort(
                "events",
                descending=True,
            )
        )

        manifest_counts = {
            x["source"]: int(x["observed"])
            for x in manifest.get("sources", [])
        }

        require_equal(
            "Per-source event counts",
            {
                x["source"]: x["events"]
                for x in summary.to_dicts()
            },
            manifest_counts,
        )

        # ----------------------------------------------------
        # OVERALL RELATIONSHIP DISTRIBUTION
        # ----------------------------------------------------

        overall = (
            by_source
            .group_by("relationship_type")
            .agg(
                [
                    pl.col("events")
                    .sum()
                    .alias("events"),

                    pl.col("assertions")
                    .sum()
                    .alias("assertions"),

                    pl.col("source")
                    .n_unique()
                    .alias("supporting_sources"),
                ]
            )
            .with_columns(
                [
                    (
                        pl.col("events")
                        - pl.col("assertions")
                    ).alias("additional_support_events"),

                    (
                        100
                        * pl.col("events")
                        / events_n
                    )
                    .round(3)
                    .alias("pct_all_events"),

                    (
                        100
                        * pl.col("assertions")
                        / assertions_n
                    )
                    .round(3)
                    .alias("pct_all_assertions"),
                ]
            )
            .sort(
                "events",
                descending=True,
            )
        )

        # ----------------------------------------------------
        # RELATIONSHIPS WITHIN INDIVIDUAL SOURCES
        # ----------------------------------------------------

        source_totals = summary.select(
            [
                "source",

                pl.col("events")
                .alias("total_source_events"),

                pl.col("assertions")
                .alias("total_source_assertions"),
            ]
        )

        by_source = (
            by_source
            .join(
                source_totals,
                on="source",
                how="left",
            )
            .with_columns(
                [
                    (
                        pl.col("events")
                        - pl.col("assertions")
                    ).alias("additional_support_events"),

                    (
                        100
                        * pl.col("events")
                        / pl.col("total_source_events")
                    )
                    .round(3)
                    .alias("pct_within_source_events"),

                    (
                        100
                        * pl.col("assertions")
                        / pl.col("total_source_assertions")
                    )
                    .round(3)
                    .alias("pct_within_source_assertions"),
                ]
            )
            .drop(
                "total_source_events",
                "total_source_assertions",
                "linked_events",
            )
            .sort(
                ["source", "events"],
                descending=[False, True],
            )
        )

        # ----------------------------------------------------
        # RELATIONSHIP -> EVIDENCE FAMILY CROSSWALK
        # ----------------------------------------------------

        relation_family = (
            crosswalk
            .group_by(
                [
                    "relationship_type",
                    "merge_evidence_family",
                ]
            )
            .agg(
                pl.col("events")
                .sum()
                .alias("events")
            )
            .join(
                overall.select(
                    "relationship_type",
                    pl.col("events")
                    .alias("relationship_events"),
                ),
                on="relationship_type",
                how="left",
            )
            .with_columns(
                (
                    100
                    * pl.col("events")
                    / pl.col("relationship_events")
                )
                .round(3)
                .alias("pct_within_relationship")
            )
            .drop("relationship_events")
            .sort(
                [
                    "relationship_type",
                    "events",
                ],
                descending=[False, True],
            )
        )

        source_family = (
            crosswalk
            .group_by(
                [
                    "source",
                    "relationship_type",
                    "merge_evidence_family",
                ]
            )
            .agg(
                pl.col("events")
                .sum()
                .alias("events")
            )
            .sort(
                ["source", "events"],
                descending=[False, True],
            )
        )

        # ----------------------------------------------------
        # SENSITIVITY VIEW: EXCLUDE CTD
        # ----------------------------------------------------

        non_ctd = (
            by_source
            .filter(
                pl.col("source") != "CTD"
            )
            .group_by("relationship_type")
            .agg(
                [
                    pl.col("events")
                    .sum()
                    .alias("events"),

                    pl.col("assertions")
                    .sum()
                    .alias("assertions"),

                    pl.col("source")
                    .n_unique()
                    .alias("sources"),
                ]
            )
        )

        non_ctd_events = int(
            non_ctd["events"].sum()
        )

        non_ctd = (
            non_ctd
            .with_columns(
                (
                    100
                    * pl.col("events")
                    / non_ctd_events
                )
                .round(3)
                .alias("pct_non_ctd_events")
            )
            .sort(
                "events",
                descending=True,
            )
        )

        ctd_n = int(
            summary.filter(
                pl.col("source") == "CTD"
            )["events"][0]
        )

        require_equal(
            "No-CTD plus CTD equals all events",
            non_ctd_events + ctd_n,
            events_n,
        )

        require_equal(
            "Overall relationships conserve events",
            int(overall["events"].sum()),
            events_n,
        )

        require_equal(
            "Overall relationships conserve assertions",
            int(overall["assertions"].sum()),
            assertions_n,
        )

        require_equal(
            "Relationship-family crosswalk conserves events",
            int(relation_family["events"].sum()),
            events_n,
        )

        progress.update(1)

    # ========================================================
    # PRINT COMPLETE TABLES
    # ========================================================

    print_tsv(
        "1. PER-SOURCE RELATIONSHIP COMPOSITION",
        summary,
    )

    print_tsv(
        "2. OVERALL ORIGINAL RELATIONSHIP TYPES",
        overall,
    )

    print_tsv(
        "3. RELATIONSHIP TYPES WITHIN EACH SOURCE",
        by_source,
    )

    print_tsv(
        "4. RELATIONSHIP TYPE -> EVIDENCE FAMILY CROSSWALK",
        relation_family,
    )

    print_tsv(
        "5. SOURCE / RELATIONSHIP / FAMILY",
        source_family,
    )

    print_tsv(
        "6. ORIGINAL RELATIONSHIP TYPES EXCLUDING CTD",
        non_ctd,
    )

    # ========================================================
    # PRINT EXAMPLES OF POTENTIAL CONFLICTS
    # ========================================================

    if args.examples:
        example_cols = [
            "drug_name",
            "drug_identifier",
            "disease_or_condition_name",
            "disease_or_condition_identifier",
            "supporting_sources",
            "relationship_types",
            "supporting_event_count",
            "merge_pair_id",
        ]

        cases = (
            lazy_csv(
                paths["pairs"],
                ["potential_conflict"] + example_cols,
            )
            .filter(
                pl.col("potential_conflict")
                .str.to_lowercase()
                == "true"
            )
            .select(example_cols)
            .limit(args.examples)
            .collect(engine="streaming")
        )

        print_tsv(
            "7. EXAMPLES OF POTENTIALLY CONFLICTING "
            "EXACT-IDENTITY PAIRS",
            cases,
        )

    # ========================================================
    # MANUSCRIPT INTERPRETATION
    # ========================================================

    print("\n" + "=" * 105)
    print("INTERPRETATION / MANUSCRIPT NOTES")
    print("=" * 105)

    leading = overall.row(
        0,
        named=True,
    )

    print(
        f"[PASS] {events_n:,} evidence events; "
        f"{assertions_n:,} source assertions; "
        f"{int(pair_stats['pairs']):,} exact pairs; "
        f"{int(pair_stats['conflicts']):,} flagged pairs"
    )

    print(
        f"[PASS] {events_n - assertions_n:,} additional "
        "supporting evidence events preserved"
    )

    print(
        f"[RESULT] {overall.height} original relationship_type "
        f"labels across {len(SOURCES)} sources"
    )

    print(
        f"[RESULT] Most common relationship: "
        f"{leading['relationship_type']} "
        f"({int(leading['events']):,} events; "
        f"{leading['pct_all_events']:.3f}%)"
    )

    print(
        "[CAUTION] Original relationships, source evidence "
        "types and reclassified families are different "
        "concepts; 'indication' does not prove approval."
    )

    print(
        "[CAUTION] The 87 flagged pairs are exact-key "
        "label conflicts, not 87 independently validated "
        "clinical contradictions."
    )

    print(
        "[CAUTION] This analysis preserves records and "
        "does not resolve synonyms or ontologies."
    )

    # ========================================================
    # OPTIONAL OUTPUT TABLES
    # ========================================================

    if args.save_tables:
        folder = (
            base.parent.parent
            / "tables"
        )

        folder.mkdir(
            parents=True,
            exist_ok=True,
        )

        tables = {
            "source_summary": summary,
            "overall_relationships": overall,
            "relationships_by_source": by_source,
            "relationship_family_crosswalk": relation_family,
            "source_relationship_family": source_family,
            "relationships_excluding_ctd": non_ctd,
        }

        for label, frame in tables.items():
            target = (
                folder
                / f"Analysis_03_{label}.csv"
            )

            frame.write_csv(target)

            print(f"[SAVED] {target}")

    print(
        "\nANALYSIS 03 COMPLETE — "
        "locked source and merger files unchanged"
    )


if __name__ == "__main__":
    main()
