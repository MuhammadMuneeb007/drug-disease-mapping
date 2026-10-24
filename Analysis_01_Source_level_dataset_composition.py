
#!/usr/bin/env python3
"""Analysis 01: Source-level composition of the locked 11-source dataset.

Reads only:
  ALL_INCLUDED_SOURCES_drug_disease_merged.csv
  ALL_INCLUDED_SOURCES_drug_disease_deduplicated.csv
  ALL_INCLUDED_SOURCES_drug_disease_source_collapsed.csv
  ALL_INCLUDED_SOURCES_merge_summary.json

Prints:
  - Source-level event and assertion counts
  - Additional supporting evidence events
  - Unique drug and disease name coverage
  - Drug-disease name-pair coverage
  - Evidence-family distributions
  - Cross-source support distribution
  - Optional drug/disease search examples

Does not re-merge or change any existing dataset.

Requirements:
  pip install polars tqdm
"""

from __future__ import annotations

import argparse
import csv
import json

from pathlib import Path
from datetime import datetime, timezone

import polars as pl
from tqdm import tqdm


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

SOURCES = [
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
]

EXPECTED = {
    "events": 10_895_395,
    "assertions": 4_534_385,
    "pairs": 4_470_601,
}

FILENAMES = {
    "events": (
        "ALL_INCLUDED_SOURCES_drug_disease_merged.csv"
    ),
    "assertions": (
        "ALL_INCLUDED_SOURCES_drug_disease_deduplicated.csv"
    ),
    "pairs": (
        "ALL_INCLUDED_SOURCES_drug_disease_source_collapsed.csv"
    ),
    "manifest": (
        "ALL_INCLUDED_SOURCES_merge_summary.json"
    ),
}


def section(title):
    print("\n" + "=" * 104)
    print(title)
    print("=" * 104, flush=True)


def scan_csv(path):
    return pl.scan_csv(
        str(path),
        infer_schema=False,
        ignore_errors=False,
        truncate_ragged_lines=False,
        encoding="utf8",
    )


def header(path):
    if not path.is_file() or path.stat().st_size == 0:
        raise FileNotFoundError(
            f"Missing/empty merged dataset: {path}"
        )

    with path.open(
        "r",
        encoding="utf-8-sig",
        newline="",
    ) as handle:
        result = next(
            csv.reader(handle, strict=True),
            [],
        )

    if len(result) != len(set(result)):
        raise RuntimeError(
            f"Duplicate columns found in {path}"
        )

    return result


def require(path, names):
    missing = set(names) - set(header(path))

    if missing:
        raise RuntimeError(
            f"{path.name}: required columns absent: "
            f"{sorted(missing)}"
        )


def pretty(frame):
    with pl.Config(
        tbl_rows=-1,
        tbl_cols=-1,
        tbl_width_chars=180,
        fmt_str_lengths=48,
    ):
        print(frame)


def key(column):
    """Conservative label normalisation; not ontology mapping."""
    return (
        pl.col(column)
        .fill_null("")
        .str.strip_chars()
        .str.replace_all(r"\s+", " ")
        .str.to_lowercase()
    )


def assert_equal(label, actual, expected):
    if actual != expected:
        raise RuntimeError(
            f"{label} FAILED: "
            f"observed={actual}, expected={expected}"
        )

    if isinstance(actual, int):
        print(f"[OK] {label}: {actual:,}", flush=True)
    else:
        print(f"[OK] {label}: {actual}", flush=True)


def query_examples(events_path, drug, disease, n):
    if not (drug or disease):
        return

    section(
        "OPTIONAL NAME SEARCH "
        "(LITERAL SUBSTRING; NOT ONTOLOGY MATCHING)"
    )

    query = scan_csv(events_path)

    if drug:
        query = query.filter(
            pl.col("drug_name")
            .fill_null("")
            .str.to_lowercase()
            .str.contains(
                drug.casefold(),
                literal=True,
            )
        )

    if disease:
        query = query.filter(
            pl.col("disease_or_condition_name")
            .fill_null("")
            .str.to_lowercase()
            .str.contains(
                disease.casefold(),
                literal=True,
            )
        )

    selected = [
        "drug_name",
        "drug_identifier",
        "disease_or_condition_name",
        "disease_or_condition_identifier",
        "relationship_type",
        "evidence_type",
        "merge_evidence_family",
        "source",
        "merge_event_id",
    ]

    matched = query.select(selected)

    breakdown = (
        matched
        .group_by([
            "source",
            "merge_evidence_family",
        ])
        .agg(
            pl.len().alias("events")
        )
        .sort("events", descending=True)
        .collect(engine="streaming")
    )

    print(
        f"Search terms: "
        f"drug={drug!r}, disease={disease!r}"
    )

    print(
        "Matching original evidence events: "
        f"{int(breakdown['events'].sum() or 0):,}"
    )

    pretty(breakdown)

    print(
        f"FIRST {n} matching evidence events "
        "(not sorted by clinical strength):"
    )

    pretty(
        matched
        .limit(n)
        .collect(engine="streaming")
    )

    print(
        "Name search cannot automatically find "
        "synonyms, alternative spellings or "
        "differently named formulations."
    )


def main():
    parser = argparse.ArgumentParser(
        description=__doc__
    )

    parser.add_argument(
        "--processed-dir",
        type=Path,
        default=(
            Path(__file__).resolve().parent
            / "data/processed"
        ),
    )

    parser.add_argument(
        "--drug",
        default="",
        help="Optional case-insensitive drug-name substring",
    )

    parser.add_argument(
        "--disease",
        default="",
        help="Optional disease-name substring",
    )

    parser.add_argument(
        "--examples",
        type=int,
        default=12,
    )

    parser.add_argument(
        "--save-tables",
        action="store_true",
        help=(
            "Write reproducible Analysis 01 "
            "CSV and JSON tables to tables/"
        ),
    )

    args = parser.parse_args()

    if not 1 <= args.examples <= 100:
        parser.error(
            "--examples must be 1..100"
        )

    base = (
        args.processed_dir
        .expanduser()
        .resolve()
    )

    files = {
        name: base / filename
        for name, filename in FILENAMES.items()
    }

    if not files["manifest"].is_file():
        raise FileNotFoundError(
            f"Missing locked merge manifest: "
            f"{files['manifest']}"
        )

    manifest = json.loads(
        files["manifest"].read_text(
            encoding="utf-8"
        )
    )

    if manifest.get("status") != "success":
        raise RuntimeError(
            "Merge manifest is not marked success"
        )

    require(
        files["events"],
        CORE + [
            "merge_evidence_family",
            "merge_pair_id",
        ],
    )

    require(
        files["assertions"],
        CORE + [
            "merge_pair_id",
            "supporting_event_count",
        ],
    )

    require(
        files["pairs"],
        CORE[:6] + [
            "supporting_source_count",
            "supporting_event_count",
            "supporting_assertion_count",
            "potential_conflict",
        ],
    )

    section(
        "ANALYSIS 01 | SOURCE-LEVEL COMPOSITION "
        "| LOCKED PAPER-11 DATA"
    )

    print(f"Directory: {base}")

    print(
        "Source events are preserved; "
        "assertion summarisation did NOT "
        "delete supporting events."
    )

    print(
        "Unique name coverage is label-based "
        "(case/whitespace normalisation only)."
    )

    progress = tqdm(
        total=4,
        desc="Analysis 01",
        unit="stage",
    )

    # --------------------------------------------------------
    # STAGE 1: SOURCE EVENT COUNTS AND EVIDENCE FAMILIES
    # --------------------------------------------------------

    event_source_families = (
        scan_csv(files["events"])
        .select([
            "source",
            "merge_evidence_family",
        ])
        .group_by([
            "source",
            "merge_evidence_family",
        ])
        .agg(
            pl.len().alias("events")
        )
        .collect(engine="streaming")
    )

    event_sources = (
        event_source_families
        .group_by("source")
        .agg(
            pl.col("events")
            .sum()
            .alias("source_events")
        )
    )

    progress.update(1)

    # --------------------------------------------------------
    # STAGE 2: ASSERTIONS AND NAME-BASED COVERAGE
    # --------------------------------------------------------

    a = (
        scan_csv(files["assertions"])
        .select([
            "source",
            "drug_name",
            "disease_or_condition_name",
            "merge_pair_id",
            "supporting_event_count",
        ])
        .with_columns([
            key("drug_name").alias("_drug_key"),

            key(
                "disease_or_condition_name"
            ).alias("_disease_key"),

            pl.col("supporting_event_count")
            .cast(pl.Int64),
        ])
    )

    # Count distinct source-level names and pairs.
    # A struct avoids ambiguity from joining labels
    # using a string separator.

    assertion_stats = (
        a.group_by("source")
        .agg([
            pl.len().alias(
                "source_assertions"
            ),

            pl.col("supporting_event_count")
            .sum()
            .alias("reconstructed_events"),

            pl.col("merge_pair_id")
            .n_unique()
            .alias("source_exact_identity_pairs"),

            pl.col("_drug_key")
            .n_unique()
            .alias("unique_drug_name_labels"),

            pl.col("_disease_key")
            .n_unique()
            .alias("unique_disease_name_labels"),

            pl.struct([
                "_drug_key",
                "_disease_key",
            ])
            .n_unique()
            .alias(
                "unique_drug_disease_name_pairs"
            ),
        ])
        .collect(engine="streaming")
    )

    progress.update(1)

    # --------------------------------------------------------
    # STAGE 3: GLOBAL EXACT-IDENTITY PAIR SUMMARY
    # --------------------------------------------------------

    pairs_lf = (
        scan_csv(files["pairs"])
        .select([
            pl.col("supporting_event_count")
            .cast(pl.Int64),

            pl.col("supporting_assertion_count")
            .cast(pl.Int64),

            pl.col("supporting_source_count")
            .cast(pl.Int64),

            pl.col("potential_conflict"),
        ])
    )

    pairs_totals = (
        pairs_lf
        .select([
            pl.len().alias("pairs"),

            pl.col("supporting_event_count")
            .sum()
            .alias("linked_events"),

            pl.col("supporting_assertion_count")
            .sum()
            .alias("linked_assertions"),

            (
                pl.col("potential_conflict")
                .str.to_lowercase()
                == "true"
            )
            .sum()
            .alias("conflicts"),
        ])
        .collect(engine="streaming")
        .row(0, named=True)
    )

    source_overlap = (
        pairs_lf
        .group_by("supporting_source_count")
        .agg(
            pl.len().alias(
                "exact_identity_pairs"
            )
        )
        .sort("supporting_source_count")
        .collect(engine="streaming")
    )

    progress.update(1)

    # --------------------------------------------------------
    # STAGE 4: JOIN THE SMALL SOURCE SUMMARIES
    # --------------------------------------------------------

    composition = (
        event_sources
        .join(
            assertion_stats,
            on="source",
            how="full",
            coalesce=True,
        )
        .with_columns([
            (
                pl.col("source_events")
                - pl.col("source_assertions")
            ).alias(
                "additional_supporting_events"
            ),

            (
                100
                * pl.col("source_assertions")
                / pl.col("source_events")
            )
            .round(2)
            .alias(
                "assertions_as_pct_of_events"
            ),
        ])
        .sort(
            "source_events",
            descending=True,
        )
    )

    events_n = int(
        composition["source_events"].sum()
    )

    assertions_n = int(
        composition["source_assertions"].sum()
    )

    pairs_n = int(
        pairs_totals["pairs"]
    )

    # --------------------------------------------------------
    # VALIDATE AGAINST THE LOCKED MERGE
    # --------------------------------------------------------

    assert_equal(
        "11 source labels",
        sorted(composition["source"].to_list()),
        sorted(SOURCES),
    )

    assert_equal(
        "event records",
        events_n,
        EXPECTED["events"],
    )

    assert_equal(
        "source-specific assertions",
        assertions_n,
        EXPECTED["assertions"],
    )

    assert_equal(
        "exact-identity pairs",
        pairs_n,
        EXPECTED["pairs"],
    )

    assert_equal(
        "assertion support -> events",
        int(
            composition[
                "reconstructed_events"
            ].sum()
        ),
        events_n,
    )

    assert_equal(
        "pair support -> events",
        int(
            pairs_totals["linked_events"]
        ),
        events_n,
    )

    assert_equal(
        "pair support -> assertions",
        int(
            pairs_totals["linked_assertions"]
        ),
        assertions_n,
    )

    assert_equal(
        "event family count -> events",
        int(
            event_source_families["events"].sum()
        ),
        events_n,
    )

    for field, value in [
        ("source_events", events_n),
        ("assertions", assertions_n),
        ("exact_identity_pairs", pairs_n),
    ]:
        assert_equal(
            f"manifest {field}",
            int(manifest[field]),
            value,
        )

    for item in manifest.get(
        "sources",
        [],
    ):
        row = composition.filter(
            pl.col("source") == item["source"]
        )

        assert_equal(
            f"manifest source {item['source']}",
            int(
                row["source_events"][0]
            ),
            int(item["observed"]),
        )

    progress.update(1)
    progress.close()

    # --------------------------------------------------------
    # DISPLAY 1: SOURCE-LEVEL DATASET COMPOSITION
    # --------------------------------------------------------

    section(
        "PER-SOURCE DATASET COMPOSITION"
    )

    pretty(composition)

    # --------------------------------------------------------
    # DISPLAY 2: ANALYTICAL EVIDENCE FAMILY DISTRIBUTION
    # --------------------------------------------------------

    section(
        "ANALYTICAL EVIDENCE FAMILIES"
    )

    family_totals = (
        event_source_families
        .group_by("merge_evidence_family")
        .agg(
            pl.col("events")
            .sum()
            .alias("events")
        )
        .with_columns(
            (
                pl.col("events")
                * 100
                / events_n
            )
            .round(3)
            .alias("pct_all_events")
        )
        .sort(
            "events",
            descending=True,
        )
    )

    pretty(family_totals)

    # --------------------------------------------------------
    # DISPLAY 3: CROSS-SOURCE EXACT-IDENTITY SUPPORT
    # --------------------------------------------------------

    section(
        "EXACT-IDENTITY CROSS-SOURCE SUPPORT"
    )

    pretty(source_overlap)

    # --------------------------------------------------------
    # DISPLAY 4: MANUSCRIPT-READY TAB-SEPARATED TABLE
    # --------------------------------------------------------

    section(
        "MANUSCRIPT-READY SOURCE SUMMARY"
    )

    print(
        "source\t"
        "events\t"
        "assertions\t"
        "additional_supporting_events\t"
        "exact_identity_pairs\t"
        "unique_drug_name_labels\t"
        "unique_disease_name_labels\t"
        "unique_name_pairs"
    )

    columns = [
        "source",
        "source_events",
        "source_assertions",
        "additional_supporting_events",
        "source_exact_identity_pairs",
        "unique_drug_name_labels",
        "unique_disease_name_labels",
        "unique_drug_disease_name_pairs",
    ]

    for row in composition.iter_rows(named=True):
        print(
            "\t".join(
                str(row[c])
                for c in columns
            )
        )

    # --------------------------------------------------------
    # DISPLAY 5: KEY RESULTS AND SCIENTIFIC INTERPRETATION
    # --------------------------------------------------------

    section(
        "KEY TOTALS AND INTERPRETATION"
    )

    ctd_events = int(
        composition.filter(
            pl.col("source") == "CTD"
        )["source_events"][0]
    )

    print(
        f"Full source evidence events      : "
        f"{events_n:,}"
    )

    print(
        f"Distinct source assertions       : "
        f"{assertions_n:,}"
    )

    print(
        f"Additional supporting events     : "
        f"{events_n - assertions_n:,} "
        "(NOT deleted)"
    )

    print(
        f"Exact-identity drug-disease pairs: "
        f"{pairs_n:,}"
    )

    print(
        f"Potential conflicting pairs      : "
        f"{int(pairs_totals['conflicts']):,}"
    )

    print(
        f"CTD share of events              : "
        f"{100 * ctd_events / events_n:.2f}%"
    )

    print(
        "Per-source exact-identity pair counts "
        "must NOT be added to infer a global unique count."
    )

    print(
        "Normalised name counts are for search/coverage, "
        "not proof of chemical/ontology equivalence."
    )

    print(
        "Evidence includes inferred chemical links, "
        "trials, contraindications and text mentions."
    )

    # --------------------------------------------------------
    # OPTIONAL: SAVE SUMMARY TABLES
    # --------------------------------------------------------

    if args.save_tables:
        output = (
            Path(__file__).resolve().parent
            / "tables"
        )

        output.mkdir(
            parents=True,
            exist_ok=True,
        )

        composition.write_csv(
            output
            / "Analysis_01_source_composition.csv"
        )

        family_totals.write_csv(
            output
            / "Analysis_01_evidence_family_composition.csv"
        )

        source_overlap.write_csv(
            output
            / "Analysis_01_exact_identity_source_overlap.csv"
        )

        summary = {
            "timestamp_utc": (
                datetime.now(timezone.utc).isoformat()
            ),
            "data_snapshot": "paper11-locked-v1",
            "polars_version": pl.__version__,
            "events": events_n,
            "assertions": assertions_n,
            "exact_identity_pairs": pairs_n,
            "potential_conflicting_pairs": int(
                pairs_totals["conflicts"]
            ),
            "ctd_events": ctd_events,
            "note": (
                "Label-based name coverage; "
                "source assertions retain evidence-event links."
            ),
        }

        (
            output / "Analysis_01_summary.json"
        ).write_text(
            json.dumps(summary, indent=2) + "\n",
            encoding="utf-8",
        )

        print(
            "\n[SAVED] Source composition, "
            "evidence families, overlap and JSON summary: "
            f"{output}"
        )

    # --------------------------------------------------------
    # OPTIONAL: DRUG/DISEASE SEARCH DEMONSTRATION
    # --------------------------------------------------------

    query_examples(
        files["events"],
        args.drug.strip(),
        args.disease.strip(),
        args.examples,
    )

    section(
        "ANALYSIS 01 COMPLETE | INPUTS UNCHANGED"
    )


if __name__ == "__main__":
    main()
