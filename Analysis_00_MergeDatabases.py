
#!/usr/bin/env python3
"""Polars streaming merger, 11 manuscript sources only.

Dependencies: pip install -U polars tqdm
Run: python mergedanalysis.py --recount-merged | tee mergedanalysis.log

Full events are never deduplicated. Assertion and exact-identity pair
summaries link to events using two-seed Polars hash IDs, which are
version-dependent; the manifest records the Polars version. No ontology
crosswalk, regulatory adjudication, or independent-evidence claim is made.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import shutil
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import polars as pl
from tqdm import tqdm

CORE = [
    "drug_name", "drug_identifier", "drug_identifier_type",
    "disease_or_condition_name", "disease_or_condition_identifier",
    "disease_or_condition_identifier_type", "relationship_type",
    "evidence_type", "source", "internal_source",
]
IDENTITY = CORE[:6]
ASSERTION = CORE[:9]
MISSING = ["", "NA", "N/A", "NAN", "NULL", "NONE"]

SOURCES = [
    ("AACT", "aact/aact_drug_disease.csv", 688745),
    (
        "Broad Drug Repurposing Hub",
        "broad_repurposing/broad_repurposing_drug_disease.csv",
        3892,
    ),
    ("ChEMBL", "chembl_drug_disease.csv", 60055),
    ("CTD", "ctd_aggregate_chemical_disease.csv", 9735698),
    ("DrugCentral", "drugcentral_drug_disease.csv", 42515),
    ("Hetionet", "hetionet/hetionet_drug_disease.csv", 1145),
    ("MEDI", "medi/medi_drug_disease.csv", 222412),
    ("OpenTargets", "opentargets_drug_disease.csv", 53950),
    ("PrimeKG", "primekg/primekg_drug_disease.csv", 42631),
    ("repoDB", "repodb/repodb_drug_disease.csv", 13558),
    ("SIDER", "sider/sider_drug_disease.csv", 30794),
]

FILENAMES = {
    "merged": "ALL_INCLUDED_SOURCES_drug_disease_merged.csv",
    "dedup": "ALL_INCLUDED_SOURCES_drug_disease_deduplicated.csv",
    "pairs": "ALL_INCLUDED_SOURCES_drug_disease_source_collapsed.csv",
    "summary": "ALL_INCLUDED_SOURCES_merge_summary.json",
}

ADDED = [
    "merge_event_id",
    "merge_assertion_id",
    "merge_pair_id",
    "merge_source_row_number",
    "merge_input_file",
    "merge_original_source",
    "merge_original_internal_source",
    "merge_original_drug_name",
    "merge_original_disease_name",
    "merge_drug_name_status",
    "merge_disease_name_status",
    "merge_evidence_family",
    "merge_bucket",
]


def now():
    return datetime.now(timezone.utc).isoformat()


def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def csv_header(path):
    if not path.is_file() or not path.stat().st_size:
        raise RuntimeError(f"Missing manuscript input: {path}")

    with path.open(
        "r", encoding="utf-8-sig", newline=""
    ) as f:
        head = next(csv.reader(f, strict=True), None)

    if (
        not head
        or len(head) != len(set(head))
        or head[:10] != CORE
    ):
        raise RuntimeError(
            f"Nonconforming 10-column core "
            f"or duplicate header: {path}: {head}"
        )

    if set(head).intersection(ADDED):
        raise RuntimeError(
            f"Reserved merge column in {path}"
        )

    return head


def blank(column):
    return (
        pl.col(column)
        .fill_null("")
        .str.strip_chars()
        .str.to_uppercase()
        .is_in(MISSING)
    )


def normalise_evidence(source, columns):
    """Source-aware interpretation, not independent approval."""

    r = pl.col("relationship_type").str.to_lowercase()
    e = pl.col("evidence_type").str.to_lowercase()

    if "clinical_stage_is_explicit_approval" in columns:
        stage = (
            pl.col("clinical_stage_is_explicit_approval")
            .fill_null("")
            .str.to_lowercase() == "true"
        )
    else:
        stage = pl.lit(False)

    return (
        pl.when(r == "contraindication")
        .then(pl.lit("contraindication"))

        .when(r == "discontinued_development")
        .then(pl.lit("halted_development"))

        .when(source == "AACT")
        .then(pl.lit("trial_cooccurrence"))

        .when(source == "CTD")
        .then(
            pl.when(r == "ctd_direct_therapeutic")
            .then(pl.lit("ctd_direct_therapeutic"))
            .when(r.str.contains("marker|mechanism"))
            .then(pl.lit("ctd_marker_mechanism"))
            .otherwise(pl.lit("ctd_inferred_or_other"))
        )

        .when(source == "SIDER")
        .then(
            pl.when(r == "label_nlp_indication")
            .then(pl.lit("label_nlp_indication"))
            .when(r == "label_precondition")
            .then(pl.lit("label_precondition"))
            .otherwise(
                pl.lit("label_text_mention_or_other")
            )
        )

        .when(
            (source == "repoDB")
            & (e == "repodb_approved_indication")
        )
        .then(pl.lit("source_reported_approval"))

        .when(
            (source == "OpenTargets") & stage
        )
        .then(pl.lit("source_reported_approval"))

        .when(r == "off_label_use")
        .then(pl.lit("off_label_use"))

        .when(r == "compound_palliates_disease")
        .then(pl.lit("palliative_relation"))

        .when(
            r.is_in([
                "indication",
                "drug_indication",
                "clinical_indication",
                "compound_treats_disease",
            ])
        )
        .then(pl.lit("source_reported_indication"))

        .otherwise(pl.lit("other_or_unclassified"))
        .alias("merge_evidence_family")
    )


def hashed_id(cols, prefix):
    """Version-pinned Polars hash IDs; not ontology crosswalks."""

    x = pl.struct([pl.col(c) for c in cols])

    return pl.concat_str(
        [
            pl.lit(prefix),
            x.hash(seed=12345).cast(pl.String),
            pl.lit("_"),
            x.hash(seed=90817).cast(pl.String),
        ],
        separator="",
    )


def make_source_lf(name, number, rel, path, header, buckets):
    lf = pl.scan_csv(
        path,
        infer_schema=False,
        row_index_name="_source_row_zero",
        ignore_errors=False,
        truncate_ragged_lines=False,
        encoding="utf8",
    )

    lf = lf.with_columns([
        pl.col(c).fill_null("")
        for c in header
    ])

    lf = lf.with_columns([
        pl.col("source").alias("merge_original_source"),
        pl.col("internal_source").alias(
            "merge_original_internal_source"
        ),
        pl.col("drug_name").alias(
            "merge_original_drug_name"
        ),
        pl.col("disease_or_condition_name").alias(
            "merge_original_disease_name"
        ),
    ])

    lf = lf.with_columns([
        pl.col(c).str.strip_chars().alias(c)
        for c in CORE
    ])

    lf = lf.with_columns([
        pl.when(
            blank("drug_name")
            & ~blank("drug_identifier")
        )
        .then(pl.col("drug_identifier"))
        .otherwise(pl.col("drug_name"))
        .alias("drug_name"),

        pl.when(
            blank("disease_or_condition_name")
            & ~blank("disease_or_condition_identifier")
        )
        .then(pl.col("disease_or_condition_identifier"))
        .otherwise(pl.col("disease_or_condition_name"))
        .alias("disease_or_condition_name"),

        pl.when(blank("drug_name"))
        .then(pl.lit(
            "identifier_placeholder_or_unresolved"
        ))
        .otherwise(pl.lit("source_name"))
        .alias("merge_drug_name_status"),

        pl.when(blank("disease_or_condition_name"))
        .then(pl.lit(
            "identifier_placeholder_or_unresolved"
        ))
        .otherwise(pl.lit("source_name"))
        .alias("merge_disease_name_status"),

        pl.lit(name).alias("source"),
        pl.lit(rel).alias("merge_input_file"),

        (pl.col("_source_row_zero") + 1)
        .alias("_source_row_one"),
    ])

    lf = lf.with_columns([
        pl.concat_str([
            pl.lit(f"{number:02d}:"),
            pl.col("_source_row_one").cast(pl.String),
        ]).alias("merge_event_id"),

        pl.col("_source_row_one")
        .cast(pl.String)
        .alias("merge_source_row_number"),

        normalise_evidence(name, header),
    ])

    lf = lf.with_columns([
        hashed_id(IDENTITY, "P").alias(
            "merge_pair_id"
        ),

        hashed_id(ASSERTION, "A").alias(
            "merge_assertion_id"
        ),

        (
            pl.struct([
                pl.col(c) for c in IDENTITY
            ]).hash(seed=771) % buckets
        )
        .cast(pl.Int32)
        .alias("merge_bucket"),
    ])

    return lf.drop(
        "_source_row_zero",
        "_source_row_one",
    )


def count_lf(lf):
    return int(
        lf.select(pl.len())
        .collect(engine="streaming")
        .item()
    )


def write_grouped_buckets(outputs, bucket_groups):
    dedup_cols = CORE + [
        "merge_assertion_id",
        "merge_pair_id",
        "supporting_event_count",
        "representative_event_id",
    ]

    pair_cols = IDENTITY + [
        "merge_pair_id",
        "supporting_event_count",
        "supporting_assertion_count",
        "supporting_source_count",
        "supporting_sources",
        "relationship_types",
        "evidence_types",
        "evidence_families",
        "indication_like_event_count",
        "contraindication_event_count",
        "potential_conflict",
        "identity_note",
    ]

    totals = {
        "assertions": 0,
        "pairs": 0,
        "events": 0,
        "potential_conflicts": 0,
    }

    with (
        outputs["dedup"].open(
            "w", encoding="utf8", newline=""
        ) as da,
        outputs["pairs"].open(
            "w", encoding="utf8", newline=""
        ) as pa,
    ):
        csv.writer(da).writerow(dedup_cols)
        csv.writer(pa).writerow(pair_cols)

        for part_files in tqdm(
            bucket_groups,
            desc="Aggregating identity buckets",
            unit="bucket",
        ):
            # Read all fragments belonging to the same
            # hash bucket together to avoid splitting groups.

            event = pl.read_parquet(part_files)

            if event.height == 0:
                continue

            totals["events"] += event.height

            # Distinct source-level assertions.
            # Each assertion retains an exact support count.
            a = event.group_by(ASSERTION).agg([
                pl.len().alias(
                    "supporting_event_count"
                ),
                pl.col("merge_event_id")
                .min()
                .alias("representative_event_id"),
                pl.col("merge_assertion_id").first(),
                pl.col("merge_pair_id").first(),
                pl.col("internal_source").first(),
            ])

            if (
                a.select(
                    pl.col(
                        "merge_assertion_id"
                    ).n_unique()
                ).item()
                != a.height
            ):
                raise RuntimeError(
                    "Assertion hash collision detected"
                )

            a = a.select(dedup_cols)
            a.write_csv(
                da, include_header=False
            )

            totals["assertions"] += a.height

            # Exact six-column identity groups.
            # These are not cross-ontology mappings.
            b = event.group_by(IDENTITY).agg([
                pl.len().alias(
                    "supporting_event_count"
                ),

                pl.col("merge_assertion_id")
                .n_unique()
                .alias("supporting_assertion_count"),

                pl.col("source")
                .n_unique()
                .alias("supporting_source_count"),

                pl.col("source")
                .unique()
                .sort()
                .alias("_sources"),

                pl.col("relationship_type")
                .unique()
                .sort()
                .alias("_rels"),

                pl.col("evidence_type")
                .unique()
                .sort()
                .alias("_evidence"),

                pl.col("merge_evidence_family")
                .unique()
                .sort()
                .alias("_families"),

                pl.col("merge_pair_id").first(),

                (
                    pl.col("merge_evidence_family")
                    .is_in([
                        "source_reported_approval",
                        "source_reported_indication",
                        "ctd_direct_therapeutic",
                        "label_nlp_indication",
                        "palliative_relation",
                    ])
                )
                .sum()
                .alias("indication_like_event_count"),

                (
                    pl.col("merge_evidence_family")
                    == "contraindication"
                )
                .sum()
                .alias("contraindication_event_count"),
            ])

            if (
                b.select(
                    pl.col("merge_pair_id").n_unique()
                ).item()
                != b.height
            ):
                raise RuntimeError(
                    "Pair hash collision detected"
                )

            b = b.with_columns([
                pl.col("_sources")
                .list.join("|")
                .alias("supporting_sources"),

                pl.col("_rels")
                .list.join("|")
                .alias("relationship_types"),

                pl.col("_evidence")
                .list.join("|")
                .alias("evidence_types"),

                pl.col("_families")
                .list.join("|")
                .alias("evidence_families"),

                (
                    (pl.col(
                        "indication_like_event_count"
                    ) > 0)
                    &
                    (pl.col(
                        "contraindication_event_count"
                    ) > 0)
                ).alias("potential_conflict"),

                pl.lit(
                    "exact six identity columns only; "
                    "NOT ontology-resolved"
                ).alias("identity_note"),
            ]).select(pair_cols)

            b.write_csv(
                pa, include_header=False
            )

            totals["pairs"] += b.height

            totals["potential_conflicts"] += int(
                b["potential_conflict"].sum()
            )

            # Every event must be represented
            # once in each summarisation scheme.
            if (
                int(
                    a["supporting_event_count"].sum()
                ) != event.height
                or
                int(
                    b["supporting_event_count"].sum()
                ) != event.height
            ):
                raise RuntimeError(
                    "Source-event conservation failed "
                    "inside identity bucket"
                )

            if (
                b["supporting_assertion_count"].sum()
                != a.height
            ):
                raise RuntimeError(
                    "Assertion-to-pair conservation failed"
                )

    return totals, dedup_cols, pair_cols


def main():
    ap = argparse.ArgumentParser(
        description=__doc__
    )

    root = Path(
        os.environ.get(
            "DRUG_DISEASE_ROOT",
            Path(__file__).resolve().parent,
        )
    )

    ap.add_argument(
        "--processed-dir",
        type=Path,
        default=root / "data/processed",
    )

    ap.add_argument(
        "--buckets",
        type=int,
        default=128,
    )

    ap.add_argument(
        "--allow-count-drift",
        action="store_true",
    )

    ap.add_argument(
        "--recount-merged",
        action="store_true",
    )

    ap.add_argument(
        "--keep-work",
        action="store_true",
    )

    args = ap.parse_args()

    if args.buckets < 8 or args.buckets > 1024:
        ap.error(
            "--buckets must be 8..1024"
        )

    base = (
        args.processed_dir
        .expanduser()
        .resolve()
    )

    if not base.is_dir():
        ap.error(
            f"Processed source directory missing: {base}"
        )

    work = base / ".polars_merge_work"

    if work.exists():
        raise RuntimeError(
            f"Existing scratch folder: {work}; "
            "inspect/remove it explicitly"
        )

    paths = {
        k: base / v
        for k, v in FILENAMES.items()
    }

    staged = {
        k: base / (v + ".merging.tmp")
        for k, v in FILENAMES.items()
    }

    if any(
        p.exists() for p in staged.values()
    ):
        raise RuntimeError(
            "A staged merger output already exists; "
            "inspect/remove *.merging.tmp"
        )

    # Validate all eleven sources before
    # starting any processing.
    specs = []

    for j, (
        source,
        rel,
        expected,
    ) in enumerate(SOURCES, 1):
        path = base / rel

        specs.append((
            j,
            source,
            rel,
            expected,
            path,
            csv_header(path),
        ))

    work.mkdir(parents=True)
    completed = False

    try:
        print(
            "[STEP 1] Polars streaming source "
            "validation and Parquet staging",
            flush=True,
        )

        manifests = []
        individual = []

        for (
            j,
            source,
            rel,
            expected,
            path,
            header,
        ) in tqdm(
            specs,
            desc="Eleven sources",
            unit="source",
        ):
            lf = make_source_lf(
                source,
                j,
                rel,
                path,
                header,
                args.buckets,
            )

            temp = (
                work / f"source_{j:02d}.parquet"
            )

            lf.sink_parquet(
                temp,
                engine="streaming",
                compression="zstd",
            )

            narrow = pl.scan_parquet(temp)
            n = count_lf(narrow)

            invalid = count_lf(
                narrow.filter(
                    blank("drug_name")
                    | blank("disease_or_condition_name")
                    | blank("relationship_type")
                    | blank("evidence_type")
                    | blank("merge_original_source")
                )
            )

            if invalid:
                raise RuntimeError(
                    f"{source}: {invalid} ineligible "
                    "or missing-evidence rows; "
                    "no silent drop"
                )

            if (
                not args.allow_count_drift
                and n != expected
            ):
                raise RuntimeError(
                    f"{source}: expected {expected:,}, "
                    f"observed {n:,}; stop"
                )

            manifests.append({
                "source": source,
                "path": rel,
                "expected": expected,
                "observed": n,
                "sha256": digest(path),
                "source_columns": header,
            })

            individual.append(temp)

            tqdm.write(
                f"[OK] {source}: {n:,} source events"
            )

        expected_sum = sum(
            m["observed"]
            for m in manifests
        )

        print(
            f"[STEP 2] Stream the COMPLETE "
            f"event ledger "
            f"({expected_sum:,} source events)",
            flush=True,
        )

        # Diagonal concatenation retains
        # all source-specific columns.
        full = pl.concat(
            [
                pl.scan_parquet(p)
                for p in individual
            ],
            how="diagonal_relaxed",
        )

        event_fields = [
            c for c in full.collect_schema().names()
            if c not in CORE
        ]

        full = full.select(
            CORE + event_fields
        )

        full.sink_csv(
            staged["merged"],
            engine="streaming",
            maintain_order=True,
        )

        print(
            "[STEP 3] Bucket narrow fields for "
            "bounded-memory grouping",
            flush=True,
        )

        narrow = full.select(
            CORE + [
                "merge_event_id",
                "merge_assertion_id",
                "merge_pair_id",
                "merge_evidence_family",
                "merge_bucket",
            ]
        )

        bucket_root = work / "buckets"

        narrow.sink_parquet(
            pl.PartitionBy(
                bucket_root,
                key="merge_bucket",
                include_key=True,
            ),
            engine="streaming",
            compression="zstd",
            mkdir=True,
        )

        bucket_files = sorted(
            bucket_root.rglob("*.parquet")
        )

        if not bucket_files:
            raise RuntimeError(
                "No bucket Parquet files written"
            )

        # All fragments of one bucket must
        # be read together before grouping.
        grouped = defaultdict(list)

        for part_file in bucket_files:
            grouped[part_file.parent].append(
                part_file
            )

        bucket_groups = [
            grouped[p]
            for p in sorted(grouped)
        ]

        print(
            f"[STEP 4] Aggregate "
            f"{len(bucket_groups)} buckets "
            f"across {len(bucket_files)} fragments",
            flush=True,
        )

        (
            counts,
            assertion_cols,
            pair_cols,
        ) = write_grouped_buckets(
            staged,
            bucket_groups,
        )

        if (
            counts["events"] != expected_sum
            or counts["assertions"] > counts["events"]
            or counts["pairs"] > counts["assertions"]
        ):
            raise RuntimeError(
                "Global reconciliation failed: "
                f"{counts}; inputs={expected_sum}"
            )

        if args.recount_merged:
            print(
                "[STEP 5] Recount merged CSV "
                "(extra scan)",
                flush=True,
            )

            actual = count_lf(
                pl.scan_csv(
                    staged["merged"],
                    infer_schema=False,
                )
            )

            if actual != expected_sum:
                raise RuntimeError(
                    f"Merged CSV readback "
                    f"{actual:,} != {expected_sum:,}"
                )

        report = {
            "status": "success",
            "completed_at_utc": now(),
            "polars_version": pl.__version__,
            "scope": "11 manuscript sources only",
            "pinned_counts_enforced":
                not args.allow_count_drift,
            "sources": manifests,
            "source_events": expected_sum,
            "assertions": counts["assertions"],
            "exact_identity_pairs": counts["pairs"],
            "potential_conflicting_pairs":
                counts["potential_conflicts"],
            "assertion_key_columns": ASSERTION,
            "pair_key_columns": IDENTITY,
            "note": (
                "Exact lexical/namespace equality only. "
                "Distinct evidence classes preserved. "
                "Polars hash keys depend on installed version."
            ),
            "outputs": {
                k: str(v)
                for k, v in paths.items()
            },
            "output_sha256": {
                k: digest(staged[k])
                for k in (
                    "merged",
                    "dedup",
                    "pairs",
                )
            },
            "schema": {
                "full": CORE + event_fields,
                "assertions": assertion_cols,
                "pairs": pair_cols,
            },
        }

        staged["summary"].write_text(
            json.dumps(
                report,
                indent=2,
                ensure_ascii=False,
            ) + "\n",
            encoding="utf8",
        )

        print(
            "[STEP 6] Publish validated outputs "
            "(summary last)",
            flush=True,
        )

        # Each replacement is atomic individually,
        # but the four-file set is not transactional.
        for k in (
            "merged",
            "dedup",
            "pairs",
            "summary",
        ):
            os.replace(
                staged[k],
                paths[k],
            )

            print(
                f"[PUBLISHED] {paths[k]}",
                flush=True,
            )

        completed = True

        print(
            f"SUCCESS: {expected_sum:,} events; "
            f"{counts['assertions']:,} assertions; "
            f"{counts['pairs']:,} exact-identity pairs; "
            f"{counts['potential_conflicts']:,} "
            f"potential conflicts",
            flush=True,
        )

    finally:
        if completed or not args.keep_work:
            shutil.rmtree(
                work,
                ignore_errors=True,
            )

        if not completed:
            for path in staged.values():
                path.unlink(missing_ok=True)


if __name__ == "__main__":
    main()
