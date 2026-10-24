
#!/usr/bin/env python3
"""Analysis 08: name-level source overlap; read-only with respect to locked data.

Combines source overlap, multi-source support and 55 pairwise source comparisons.
Uses the exact Analysis 01/04 name key: lowercase, trim and collapse whitespace.
Does NOT merge ontologies, deduplicate source events, or infer approval.

Requires: polars, tqdm.
Run: POLARS_MAX_THREADS=16 python Analysis_08_Cross_source_drug_disease_pair_overlap.py --save-tables
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import shutil
from collections import Counter, defaultdict
from itertools import combinations
from pathlib import Path

import polars as pl
from tqdm import tqdm

SOURCES = sorted([
    "AACT", "Broad Drug Repurposing Hub", "ChEMBL", "CTD", "DrugCentral",
    "Hetionet", "MEDI", "OpenTargets", "PrimeKG", "repoDB", "SIDER",
])
BITS = {name: 1 << i for i, name in enumerate(SOURCES)}
MULTI_MASKS = [m for m in range(1 << len(SOURCES)) if m.bit_count() >= 2]

EXPECTED = {
    "assertions": 4_534_385,
    "events": 10_895_395,
    "strict_pairs": 4_470_601,
    "drug_names": 158_364,
    "disease_names": 72_563,
    "name_pairs": 4_313_189,
    "multi_source_name_pairs": 83_878,
}

ASSERT_FILE = "ALL_INCLUDED_SOURCES_drug_disease_deduplicated.csv"
PAIR_FILE = "ALL_INCLUDED_SOURCES_drug_disease_source_collapsed.csv"
MANIFEST = "ALL_INCLUDED_SOURCES_merge_summary.json"


def check(path, columns):
    if not path.is_file() or path.stat().st_size == 0:
        raise RuntimeError(f"Missing or empty input: {path}")

    with path.open("r", newline="", encoding="utf-8-sig") as f:
        head = next(csv.reader(f, strict=True), [])

    if len(head) != len(set(head)) or not set(columns).issubset(head):
        raise RuntimeError(
            f"Invalid input schema: {path.name}; need {columns}"
        )


def scan(path, fields):
    return pl.scan_csv(
        str(path),
        infer_schema=False,
        ignore_errors=False,
        truncate_ragged_lines=False,
        encoding="utf8",
    ).select(fields)


def norm(name):
    return (
        pl.col(name)
        .fill_null("")
        .str.strip_chars()
        .str.replace_all(r"\s+", " ")
        .str.to_lowercase()
    )


def source_bit():
    expr = pl.lit(0, dtype=pl.Int32)

    for name, flag in BITS.items():
        expr = (
            pl.when(pl.col("source") == name)
            .then(pl.lit(flag))
            .otherwise(expr)
        )

    return expr.alias("_bit")


def equal(label, actual, expected):
    if actual != expected:
        raise RuntimeError(
            f"[FAIL] {label}: {actual!r} != {expected!r}"
        )

    print(
        f"[PASS] {label}: {actual:,}"
        if isinstance(actual, int)
        else f"[PASS] {label}",
        flush=True,
    )


def print_table(title, df):
    print("\n" + "=" * 100)
    print(f"{title} [{df.height:,} rows]")
    print("=" * 100)
    print("\t".join(df.columns))

    for record in df.iter_rows():
        print("\t".join(
            "" if v is None
            else f"{v:.6f}" if isinstance(v, float)
            else str(v).replace("\t", " ").replace("\n", " ")
            for v in record
        ))


def masks_from_grouped(grouped):
    return Counter({
        int(k): int(n)
        for k, n in (
            grouped.group_by("mask")
            .len()
            .select("mask", "len")
            .iter_rows()
        )
    })


def summary_from_masks(label, masks):
    """Count each name key at most once per source."""
    total = sum(masks.values())

    source_n = {
        s: sum(
            n for mask, n in masks.items()
            if mask & bit
        )
        for s, bit in BITS.items()
    }

    dist = Counter()

    for mask, count in masks.items():
        if mask == 0:
            raise RuntimeError(
                f"{label}: encountered unknown source mask"
            )

        dist[mask.bit_count()] += count

    multi = sum(
        n for num, n in dist.items()
        if num >= 2
    )

    return total, source_n, dist, multi


def overlap_tables(label, masks):
    total, source_n, dist, multi = summary_from_masks(
        label, masks
    )

    long = []

    for a, b in combinations(SOURCES, 2):
        na, nb = source_n[a], source_n[b]

        shared = sum(
            n for mask, n in masks.items()
            if mask & BITS[a] and mask & BITS[b]
        )

        union = na + nb - shared

        long.append({
            "source_a": a,
            "source_b": b,
            "a_unique_keys": na,
            "b_unique_keys": nb,
            "intersection": shared,
            "union": union,
            "jaccard": shared / union if union else 0.0,
            "overlap_coefficient":
                shared / min(na, nb)
                if min(na, nb) else 0.0,
            "pct_source_a":
                100 * shared / na if na else 0.0,
            "pct_source_b":
                100 * shared / nb if nb else 0.0,
        })

    pairs = pl.DataFrame(long).sort(
        ["intersection", "jaccard"],
        descending=True,
    )

    by_source = pl.DataFrame({
        "source": SOURCES,
        "unique_name_keys": [
            source_n[s] for s in SOURCES
        ],
    }).sort(
        "unique_name_keys",
        descending=True,
    )

    distribution = pl.DataFrame({
        "n_sources": sorted(dist),
        "unique_keys": [
            dist[i] for i in sorted(dist)
        ],
    }).with_columns(
        (
            100 * pl.col("unique_keys") / total
        ).round(4).alias("percent_of_keys")
    )

    matrices = {}

    for metric in (
        "intersection",
        "jaccard",
        "overlap_coefficient",
    ):
        rows = []

        for a in SOURCES:
            row = {"source": a}

            for b in SOURCES:
                if a == b:
                    row[b] = (
                        source_n[a]
                        if metric == "intersection"
                        else 1.0
                    )
                else:
                    shared = sum(
                        n for mask, n in masks.items()
                        if mask & BITS[a] and mask & BITS[b]
                    )

                    union = (
                        source_n[a]
                        + source_n[b]
                        - shared
                    )

                    if metric == "intersection":
                        row[b] = shared

                    elif metric == "jaccard":
                        row[b] = (
                            shared / union
                            if union else 0.0
                        )

                    else:
                        denom = min(
                            source_n[a], source_n[b]
                        )

                        row[b] = (
                            shared / denom
                            if denom else 0.0
                        )

            rows.append(row)

        matrices[metric] = pl.DataFrame(rows)

    top = pairs.row(0, named=True)

    description = {
        "level": label,
        "unique_name_keys": total,
        "single_source_keys": dist.get(1, 0),
        "multi_source_keys": multi,
        "pct_multi_source": round(
            100 * multi / total, 4
        ),
        "max_sources_per_key": max(dist),
        "nonzero_source_pairs": int(
            pairs.filter(
                pl.col("intersection") > 0
            ).height
        ),
        "largest_overlap_sources":
            f"{top['source_a']} | {top['source_b']}",
        "largest_overlap_keys":
            int(top["intersection"]),
        "largest_overlap_jaccard":
            round(top["jaccard"], 6),
    }

    equal(
        f"{label} source support histogram",
        int(distribution["unique_keys"].sum()),
        total,
    )

    equal(
        f"{label} number of pairwise source comparisons",
        pairs.height,
        55,
    )

    return {
        "summary": description,
        "per_source": by_source,
        "support": distribution,
        "pairwise": pairs,
        "matrices": matrices,
    }


def main():
    ap = argparse.ArgumentParser(
        description=__doc__
    )

    ap.add_argument(
        "--processed-dir",
        type=Path,
        default=(
            Path(__file__).resolve().parent
            / "data" / "processed"
        ),
    )

    ap.add_argument(
        "--buckets",
        type=int,
        default=128,
    )

    ap.add_argument(
        "--top",
        type=int,
        default=30,
    )

    ap.add_argument(
        "--save-tables",
        action="store_true",
    )

    ap.add_argument(
        "--keep-work",
        action="store_true",
    )

    args = ap.parse_args()

    if (
        not 8 <= args.buckets <= 1024
        or not 1 <= args.top <= 100
    ):
        ap.error(
            "--buckets must be 8..1024 and --top 1..100"
        )

    base = args.processed_dir.expanduser().resolve()

    af = base / ASSERT_FILE
    pf = base / PAIR_FILE
    mf = base / MANIFEST

    required = [
        "source",
        "drug_name",
        "disease_or_condition_name",
        "relationship_type",
        "evidence_type",
        "supporting_event_count",
    ]

    check(af, required)

    check(
        pf,
        ["supporting_source_count", "merge_pair_id"],
    )

    if not mf.is_file():
        raise RuntimeError(
            f"Missing merge manifest: {mf}"
        )

    manifest = json.loads(
        mf.read_text(encoding="utf8")
    )

    equal(
        "merge status",
        manifest.get("status"),
        "success",
    )

    equal(
        "manifest assertions",
        int(manifest["assertions"]),
        EXPECTED["assertions"],
    )

    equal(
        "manifest events",
        int(manifest["source_events"]),
        EXPECTED["events"],
    )

    equal(
        "manifest exact-ID pairs",
        int(manifest["exact_identity_pairs"]),
        EXPECTED["strict_pairs"],
    )

    work = base / ".analysis08_polars_work"

    if work.exists():
        raise RuntimeError(
            f"Existing scratch; inspect or relocate "
            f"before rerun: {work}"
        )

    work.mkdir()

    completed = False

    output_dir = (
        base
        / "analysis_08_cross_source_overlap_namekeys_v2"
    )

    writer = None
    full_handle = None
    full_staged = None

    try:
        print(
            "[STEP 1] Stage narrow assertion records as Parquet",
            flush=True,
        )

        narrow_file = work / "assertion_keys.parquet"

        lf = scan(af, required).with_columns([
            norm("drug_name").alias("_drug"),
            norm(
                "disease_or_condition_name"
            ).alias("_disease"),
            pl.col("supporting_event_count")
            .cast(pl.Int64, strict=True),
            source_bit(),
        ])

        lf.sink_parquet(
            narrow_file,
            engine="streaming",
            compression="zstd",
        )

        data = pl.scan_parquet(narrow_file)

        stats = data.select([
            pl.len().alias("assertions"),
            pl.col("supporting_event_count")
            .sum().alias("events"),
            (pl.col("_bit") == 0)
            .sum().alias("invalid_sources"),
            (
                (pl.col("_drug") == "")
                | (pl.col("_disease") == "")
            ).sum().alias("empty_keys"),
        ]).collect(
            engine="streaming"
        ).row(0, named=True)

        equal(
            "source assertions",
            int(stats["assertions"]),
            EXPECTED["assertions"],
        )

        equal(
            "supporting source events",
            int(stats["events"]),
            EXPECTED["events"],
        )

        equal(
            "unrecognised source labels",
            int(stats["invalid_sources"]),
            0,
        )

        equal(
            "blank name keys",
            int(stats["empty_keys"]),
            0,
        )

        # --------------------------------------------------
        # Drug and disease coverage
        # --------------------------------------------------

        print(
            "[STEP 2] Calculate distinct drug and disease "
            "name source membership",
            flush=True,
        )

        masks = {}

        for name, key in (
            ("drug", "_drug"),
            ("disease", "_disease"),
        ):
            grouped = (
                data.group_by(key)
                .agg(
                    pl.col("_bit")
                    .unique()
                    .sum()
                    .alias("mask")
                )
                .collect(engine="streaming")
            )

            masks[name] = masks_from_grouped(grouped)

            equal(
                f"unique {name} names",
                grouped.height,
                EXPECTED[
                    "drug_names"
                    if name == "drug"
                    else "disease_names"
                ],
            )

        del grouped

        # --------------------------------------------------
        # Pair partitions
        # --------------------------------------------------

        print(
            "[STEP 3] Partition name pairs "
            "for bounded-memory processing",
            flush=True,
        )

        folder = work / "name_pair_buckets"

        data.with_columns(
            (
                pl.struct(
                    ["_drug", "_disease"]
                ).hash(seed=202604)
                % args.buckets
            )
            .cast(pl.Int32)
            .alias("_bucket")
        ).sink_parquet(
            pl.PartitionBy(
                folder,
                key="_bucket",
                include_key=True,
            ),
            engine="streaming",
            compression="zstd",
            mkdir=True,
        )

        files_by_bucket = defaultdict(list)

        for part in folder.rglob("*.parquet"):
            b = int(
                pl.read_parquet(
                    part,
                    columns=["_bucket"],
                    n_rows=1,
                )["_bucket"][0]
            )

            files_by_bucket[b].append(part)

        if not files_by_bucket:
            raise RuntimeError(
                "No name-pair partitions created"
            )

        # Save complete multi-source pairs only when requested.
        if args.save_tables:
            output_dir.mkdir(
                parents=True,
                exist_ok=True,
            )

            full_staged = (
                output_dir
                / "all_multisource_name_pairs.csv.merging.tmp"
            )

            if full_staged.exists():
                raise RuntimeError(
                    "Previous staged output needs inspection: "
                    f"{full_staged}"
                )

            full_handle = full_staged.open(
                "w",
                newline="",
                encoding="utf-8",
            )

            cols = [
                "drug_name_key",
                "disease_name_key",
                "n_sources",
                "assertions",
                "evidence_events",
                "sources_present",
                "n_relationship_types",
                "relationship_types_present",
                "n_evidence_types",
                "evidence_types_present",
            ]

            writer = csv.DictWriter(
                full_handle,
                fieldnames=cols,
            )

            writer.writeheader()

        # --------------------------------------------------
        # Aggregate each pair partition
        # --------------------------------------------------

        masks["pair"] = Counter()

        seen_assertions = 0
        events_in_buckets = 0

        source_multi_pairs = Counter()
        relation_multi_pairs = Counter()
        evidence_multi_pairs = Counter()

        leading_pairs = []

        print(
            "[STEP 4] Aggregate pair source masks "
            "and evidence details",
            flush=True,
        )

        for bucket, part_paths in tqdm(
            sorted(files_by_bucket.items()),
            desc="Name-pair buckets",
            unit="bucket",
        ):
            block = pl.read_parquet(
                part_paths
            ).drop("_bucket")

            seen_assertions += block.height

            events_in_buckets += int(
                block["supporting_event_count"].sum()
            )

            pair_groups = block.group_by(
                ["_drug", "_disease"]
            ).agg([
                pl.col("_bit")
                .unique()
                .sum()
                .alias("mask"),

                pl.len().alias("assertions"),

                pl.col("supporting_event_count")
                .sum()
                .alias("evidence_events"),
            ])

            masks["pair"].update(
                masks_from_grouped(pair_groups)
            )

            multi_keys = pair_groups.filter(
                pl.col("mask").is_in(MULTI_MASKS)
            ).select(
                ["_drug", "_disease"]
            )

            if not multi_keys.height:
                continue

            matches = block.join(
                multi_keys,
                on=["_drug", "_disease"],
                how="inner",
            )

            # Preserve every relationship/evidence annotation.
            pair_details = matches.group_by(
                ["_drug", "_disease"]
            ).agg([
                pl.col("source")
                .n_unique()
                .alias("n_sources"),

                pl.len().alias("assertions"),

                pl.col("supporting_event_count")
                .sum()
                .alias("evidence_events"),

                pl.col("source")
                .unique().sort()
                .str.join(" | ")
                .alias("sources_present"),

                pl.col("relationship_type")
                .n_unique()
                .alias("n_relationship_types"),

                pl.col("relationship_type")
                .unique().sort()
                .str.join(" | ")
                .alias("relationship_types_present"),

                pl.col("evidence_type")
                .n_unique()
                .alias("n_evidence_types"),

                pl.col("evidence_type")
                .unique().sort()
                .str.join(" | ")
                .alias("evidence_types_present"),
            ]).rename({
                "_drug": "drug_name_key",
                "_disease": "disease_name_key",
            })

            for row in pair_details.select(
                ["sources_present"]
            ).iter_rows():
                for source in row[0].split(" | "):
                    source_multi_pairs[source] += 1

            # Count distinct PAIRS carrying each relationship.
            for rel, n in (
                matches.select([
                    "_drug",
                    "_disease",
                    "relationship_type",
                ])
                .unique()
                .group_by("relationship_type")
                .len()
                .iter_rows()
            ):
                relation_multi_pairs[rel] += int(n)

            # Count distinct PAIRS carrying each evidence type.
            for evidence, n in (
                matches.select([
                    "_drug",
                    "_disease",
                    "evidence_type",
                ])
                .unique()
                .group_by("evidence_type")
                .len()
                .iter_rows()
            ):
                evidence_multi_pairs[evidence] += int(n)

            rank = (
                pair_details.sort(
                    [
                        "n_sources",
                        "assertions",
                        "drug_name_key",
                    ],
                    descending=[
                        True,
                        True,
                        False,
                    ],
                )
                .head(args.top)
                .to_dicts()
            )

            leading_pairs.extend(rank)

            if writer:
                writer.writerows(
                    pair_details.to_dicts()
                )

        if full_handle:
            full_handle.close()
            full_handle = None

        # --------------------------------------------------
        # Reconciliation
        # --------------------------------------------------

        equal(
            "bucket assertion conservation",
            seen_assertions,
            EXPECTED["assertions"],
        )

        equal(
            "bucket event conservation",
            events_in_buckets,
            EXPECTED["events"],
        )

        equal(
            "global name pairs",
            sum(masks["pair"].values()),
            EXPECTED["name_pairs"],
        )

        pair_multi = sum(
            n for mask, n in masks["pair"].items()
            if mask.bit_count() >= 2
        )

        equal(
            "multi-source name pairs",
            pair_multi,
            EXPECTED["multi_source_name_pairs"],
        )

        equal(
            "source memberships of multi-source name pairs",
            sum(source_multi_pairs.values()),
            sum(
                mask.bit_count() * n
                for mask, n in masks["pair"].items()
                if mask.bit_count() >= 2
            ),
        )

        # --------------------------------------------------
        # Source overlap calculations
        # --------------------------------------------------

        print(
            "[STEP 5] Build pairwise matrices "
            "and overlap summaries",
            flush=True,
        )

        result = {
            level: overlap_tables(
                level, masks[level]
            )
            for level in [
                "drug",
                "disease",
                "pair",
            ]
        }

        combined = pl.DataFrame([
            result[level]["summary"]
            for level in [
                "drug",
                "disease",
                "pair",
            ]
        ])

        strict = scan(
            pf,
            [
                "supporting_source_count",
                "merge_pair_id",
            ],
        )

        strict_stats = (
            strict.with_columns(
                pl.col("supporting_source_count")
                .cast(pl.Int64)
            )
            .group_by("supporting_source_count")
            .agg(
                pl.len().alias(
                    "exact_identity_pairs"
                )
            )
            .sort("supporting_source_count")
            .collect(engine="streaming")
        )

        equal(
            "strict source-collapsed pairs",
            int(
                strict_stats[
                    "exact_identity_pairs"
                ].sum()
            ),
            EXPECTED["strict_pairs"],
        )

        equal(
            "strict multi-source pairs",
            int(
                strict_stats.filter(
                    pl.col("supporting_source_count") > 1
                )["exact_identity_pairs"].sum() or 0
            ),
            9,
        )

        support_tiers = (
            result["pair"]["support"]
            .with_columns(
                pl.when(
                    pl.col("n_sources") == 1
                ).then(pl.lit("1 source"))
                .when(
                    pl.col("n_sources") == 2
                ).then(pl.lit("2 sources"))
                .when(
                    pl.col("n_sources") == 3
                ).then(pl.lit("3 sources"))
                .when(
                    pl.col("n_sources") == 4
                ).then(pl.lit("4 sources"))
                .otherwise(pl.lit("5+ sources"))
                .alias("support_tier")
            )
            .group_by("support_tier")
            .agg(
                pl.col("unique_keys")
                .sum()
                .alias("name_pairs")
            )
        )

        source_contrib = pl.DataFrame({
            "source": SOURCES,
            "multi_source_name_pairs": [
                source_multi_pairs[s]
                for s in SOURCES
            ],
        }).with_columns(
            (
                100
                * pl.col("multi_source_name_pairs")
                / pair_multi
            ).round(3)
            .alias("pct_all_multisource_pairs")
        ).sort(
            "multi_source_name_pairs",
            descending=True,
        )

        rel_multi = pl.DataFrame({
            "relationship_type": sorted(
                relation_multi_pairs
            ),
            "multi_source_name_pairs_with_relationship": [
                relation_multi_pairs[k]
                for k in sorted(relation_multi_pairs)
            ],
        }).with_columns(
            (
                100
                * pl.col(
                    "multi_source_name_pairs_with_relationship"
                )
                / pair_multi
            ).round(3)
            .alias("pct_of_multisource_pairs_nonexclusive")
        ).sort(
            "multi_source_name_pairs_with_relationship",
            descending=True,
        )

        evid_multi = pl.DataFrame({
            "evidence_type": sorted(
                evidence_multi_pairs
            ),
            "multi_source_name_pairs_with_evidence": [
                evidence_multi_pairs[k]
                for k in sorted(evidence_multi_pairs)
            ],
        }).with_columns(
            (
                100
                * pl.col(
                    "multi_source_name_pairs_with_evidence"
                )
                / pair_multi
            ).round(3)
            .alias("pct_of_multisource_pairs_nonexclusive")
        ).sort(
            "multi_source_name_pairs_with_evidence",
            descending=True,
        )

        leading_pairs.sort(
            key=lambda v: (
                -v["n_sources"],
                -v["assertions"],
                v["drug_name_key"],
                v["disease_name_key"],
            )
        )

        top_df = pl.DataFrame(
            leading_pairs[:args.top]
        )

        high_rows = [
            r for r in leading_pairs
            if r["n_sources"] >= 3
        ][:args.top]

        high_df = (
            pl.DataFrame(high_rows)
            if high_rows
            else top_df.head(0)
        )

        # --------------------------------------------------
        # Print results
        # --------------------------------------------------

        print_table(
            "COMBINED DRUG / DISEASE / PAIR OVERLAP",
            combined,
        )

        for level in (
            "drug",
            "disease",
            "pair",
        ):
            print_table(
                f"{level.upper()}: SOURCE COUNTS",
                result[level]["per_source"],
            )

            print_table(
                f"{level.upper()}: SOURCE SUPPORT DISTRIBUTION",
                result[level]["support"],
            )

            print_table(
                f"{level.upper()}: PAIRWISE SOURCE OVERLAP "
                "(55 COMPARISONS)",
                result[level]["pairwise"],
            )

            print_table(
                f"{level.upper()}: INTERSECTION MATRIX",
                result[level]["matrices"]["intersection"],
            )

            print_table(
                f"{level.upper()}: JACCARD MATRIX",
                result[level]["matrices"]["jaccard"],
            )

            print_table(
                f"{level.upper()}: OVERLAP COEFFICIENT MATRIX",
                result[level]["matrices"]["overlap_coefficient"],
            )

        print_table(
            "PAIR SUPPORT TIERS",
            support_tiers,
        )

        print_table(
            "SOURCE CONTRIBUTION TO MULTI-SOURCE PAIRS",
            source_contrib,
        )

        print_table(
            "RELATIONSHIP TYPES SEEN IN MULTI-SOURCE PAIRS",
            rel_multi,
        )

        print_table(
            "EVIDENCE TYPES SEEN IN MULTI-SOURCE PAIRS",
            evid_multi,
        )

        print_table(
            "TOP MULTI-SOURCE NAME PAIRS",
            top_df,
        )

        print_table(
            "TOP PAIRS FOUND IN THREE OR MORE SOURCES",
            high_df,
        )

        print_table(
            "STRICT NAME+IDENTIFIER OVERLAP "
            "(COMPARISON ONLY)",
            strict_stats,
        )

        # --------------------------------------------------
        # Optional supplementary-material tables
        # --------------------------------------------------

        if args.save_tables:
            for level in [
                "drug",
                "disease",
                "pair",
            ]:
                for name in [
                    "per_source",
                    "support",
                    "pairwise",
                ]:
                    result[level][name].write_csv(
                        output_dir / f"{level}_{name}.tsv",
                        separator="\t",
                    )

                for name, df in (
                    result[level]["matrices"].items()
                ):
                    df.write_csv(
                        output_dir
                        / f"{level}_{name}_matrix.tsv",
                        separator="\t",
                    )

            tables = [
                ("combined_overlap_summary", combined),
                ("pair_support_tiers", support_tiers),
                (
                    "source_contribution_to_multisource_pairs",
                    source_contrib,
                ),
                (
                    "relationship_types_in_multisource_pairs",
                    rel_multi,
                ),
                (
                    "evidence_types_in_multisource_pairs",
                    evid_multi,
                ),
                ("top_multisource_pairs", top_df),
                (
                    "high_recurrence_pairs_3plus",
                    high_df,
                ),
                (
                    "strict_identity_source_support",
                    strict_stats,
                ),
            ]

            for name, frame in tables:
                frame.write_csv(
                    output_dir / f"{name}.tsv",
                    separator="\t",
                )

            if full_staged is not None:
                os.replace(
                    full_staged,
                    output_dir
                    / "all_multisource_name_pairs.csv",
                )

            print(
                "[SAVED] All Analysis 08 tables:",
                output_dir,
            )

        print("\n[INTERPRETATION]")

        print(
            f"{result['pair']['summary']['multi_source_keys']:,} "
            f"of {result['pair']['summary']['unique_name_keys']:,} "
            "name pairs occur in at least two databases."
        )

        print(
            "Matching search labels may refer to "
            "different formulations or concepts."
        )

        print(
            "Multiple sources do not necessarily "
            "provide independent evidence."
        )

        print(
            "Relationship and evidence labels are "
            "preserved; overlaps are not approval claims."
        )

        print(
            "ANALYSIS 08 PASS; LOCKED INPUT FILES UNCHANGED"
        )

        completed = True

    finally:
        if full_handle:
            full_handle.close()

        if (
            not completed
            and full_staged is not None
            and not args.keep_work
        ):
            full_staged.unlink(missing_ok=True)

        if completed or not args.keep_work:
            shutil.rmtree(
                work,
                ignore_errors=True,
            )
        else:
            print(
                "[SCRATCH RETAINED]",
                work,
            )


if __name__ == "__main__":
    main()
