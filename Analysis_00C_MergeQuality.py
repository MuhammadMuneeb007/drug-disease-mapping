
#!/usr/bin/env python3
"""Independent merge-quality audit for the eleven manuscript databases.

Inputs:
  ALL_INCLUDED_SOURCES_drug_disease_merged.csv
  ALL_INCLUDED_SOURCES_drug_disease_deduplicated.csv
  ALL_INCLUDED_SOURCES_drug_disease_source_collapsed.csv
  ALL_INCLUDED_SOURCES_merge_summary.json

Outputs:
  Analysis_00C_merge_quality_report.json

No downloading, source reprocessing or alteration of merged CSVs.

Run with the SAME Polars version used by the merger.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import shutil
import sys
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

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

IDENTITY = CORE[:6]
ASSERTION = CORE[:9]

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

INDICATION_LIKE = [
    "source_reported_approval",
    "source_reported_indication",
    "ctd_direct_therapeutic",
    "label_nlp_indication",
    "palliative_relation",
]

FILES = {
    "events":
        "ALL_INCLUDED_SOURCES_drug_disease_merged.csv",
    "assertions":
        "ALL_INCLUDED_SOURCES_drug_disease_deduplicated.csv",
    "pairs":
        "ALL_INCLUDED_SOURCES_drug_disease_source_collapsed.csv",
    "manifest":
        "ALL_INCLUDED_SOURCES_merge_summary.json",
}

EVENT_EXTRA = [
    "merge_event_id",
    "merge_assertion_id",
    "merge_pair_id",
    "merge_source_row_number",
    "merge_bucket",
    "merge_evidence_family",
    "merge_original_source",
    "merge_input_file",
]

ASSERT_EXTRA = [
    "merge_assertion_id",
    "merge_pair_id",
    "supporting_event_count",
    "representative_event_id",
]

PAIR_EXTRA = [
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
]


def timestamp():
    return datetime.now(timezone.utc).isoformat()


def checksum(path):
    h = hashlib.sha256()

    with path.open("rb") as fh:
        for chunk in iter(
            lambda: fh.read(8 * 1024 * 1024),
            b"",
        ):
            h.update(chunk)

    return h.hexdigest()


def original_header(path):
    if not path.is_file() or path.stat().st_size == 0:
        raise RuntimeError(
            f"Required database file missing or empty: {path}"
        )

    with path.open(
        "r",
        encoding="utf-8-sig",
        newline="",
    ) as fh:
        fields = next(
            csv.reader(fh, strict=True),
            None,
        )

    if not fields or len(fields) != len(set(fields)):
        raise RuntimeError(
            f"Missing/duplicate columns: {path}"
        )

    return fields


def scan(path):
    return pl.scan_csv(
        path,
        infer_schema=False,
        ignore_errors=False,
        truncate_ragged_lines=False,
        encoding="utf8",
    )


def twohash(cols, prefix):
    """Reproduce the original merge fingerprint."""
    s = pl.struct([
        pl.col(c) for c in cols
    ])

    return pl.concat_str(
        [
            pl.lit(prefix),
            s.hash(seed=12345).cast(pl.String),
            pl.lit("_"),
            s.hash(seed=90817).cast(pl.String),
        ],
        separator="",
    )


def audit_bucket_expr(buckets):
    return (
        (
            pl.struct([
                pl.col(c) for c in IDENTITY
            ]).hash(seed=771) % buckets
        )
        .cast(pl.Int32)
        .alias("audit_bucket")
    )


def partition(lf, folder):
    lf.sink_parquet(
        pl.PartitionBy(
            folder,
            key="audit_bucket",
            include_key=True,
        ),
        engine="streaming",
        compression="zstd",
        mkdir=True,
    )

    found = defaultdict(list)

    for path in folder.rglob("*.parquet"):
        bucket = int(
            pl.scan_parquet(path)
            .select("audit_bucket")
            .head(1)
            .collect(engine="streaming")
            .item()
        )

        found[bucket].append(path)

    if not found:
        raise RuntimeError(
            f"Partitioning yielded no files: {folder}"
        )

    return found


def n_bad(frame, predicate):
    return int(
        frame.select(
            predicate.fill_null(True).sum()
        ).item()
    )


def expected_assertions(events):
    """Reconstruct assertions from the actual event fields."""
    return events.group_by(ASSERTION).agg(
        [
            pl.len().alias("calc_events"),

            pl.col("merge_assertion_id")
            .n_unique()
            .alias("number_assertion_ids"),

            pl.col("merge_pair_id")
            .n_unique()
            .alias("number_pair_ids"),

            pl.col("merge_assertion_id")
            .first()
            .alias("calc_aid"),

            pl.col("merge_pair_id")
            .first()
            .alias("calc_pid"),
        ]
    )


def expected_pairs(events):
    """Independently reconstruct exact-identity pairs."""
    raw = events.group_by(IDENTITY).agg(
        [
            pl.len().alias("calc_events"),

            pl.col("merge_assertion_id")
            .n_unique()
            .alias("calc_assertions"),

            pl.col("source")
            .n_unique()
            .alias("calc_sources"),

            pl.col("source")
            .unique()
            .sort()
            .alias("_sources"),

            pl.col("relationship_type")
            .unique()
            .sort()
            .alias("_relationships"),

            pl.col("evidence_type")
            .unique()
            .sort()
            .alias("_evidences"),

            pl.col("merge_evidence_family")
            .unique()
            .sort()
            .alias("_families"),

            pl.col("merge_pair_id")
            .n_unique()
            .alias("number_pair_ids"),

            pl.col("merge_pair_id")
            .first()
            .alias("calc_pid"),

            pl.col("merge_evidence_family")
            .is_in(INDICATION_LIKE)
            .sum()
            .alias("calc_indication"),

            (
                pl.col("merge_evidence_family")
                == "contraindication"
            )
            .sum()
            .alias("calc_contra"),
        ]
    )

    return raw.with_columns(
        [
            pl.col("_sources")
            .list.join("|")
            .alias("calc_source_names"),

            pl.col("_relationships")
            .list.join("|")
            .alias("calc_relation_names"),

            pl.col("_evidences")
            .list.join("|")
            .alias("calc_evidence_names"),

            pl.col("_families")
            .list.join("|")
            .alias("calc_family_names"),

            (
                (pl.col("calc_indication") > 0)
                & (pl.col("calc_contra") > 0)
            ).alias("calc_conflict"),
        ]
    ).drop(
        "_sources",
        "_relationships",
        "_evidences",
        "_families",
    )


def query_pairs(
    pair_file,
    term_drug,
    term_disease,
    limit,
):
    """Run case-insensitive example queries on the pair database."""
    lf = scan(pair_file)

    if term_drug:
        lf = lf.filter(
            pl.col("drug_name")
            .str.to_lowercase()
            .str.contains(
                term_drug.casefold(),
                literal=True,
            )
        )

    if term_disease:
        lf = lf.filter(
            pl.col("disease_or_condition_name")
            .str.to_lowercase()
            .str.contains(
                term_disease.casefold(),
                literal=True,
            )
        )

    if not (term_drug or term_disease):
        return None

    columns = [
        "drug_name",
        "drug_identifier",
        "disease_or_condition_name",
        "disease_or_condition_identifier",
        "supporting_event_count",
        "supporting_source_count",
        "supporting_sources",
        "relationship_types",
        "potential_conflict",
        "merge_pair_id",
    ]

    n = int(
        lf.select(pl.len())
        .collect(engine="streaming")
        .item()
    )

    data = (
        lf.with_columns(
            pl.col("supporting_event_count")
            .cast(pl.Int64)
            .alias("_count")
        )
        .sort("_count", descending=True)
        .select(columns)
        .head(limit)
        .collect(engine="streaming")
        .to_dicts()
    )

    return {
        "matching_exact_identity_pairs": n,
        "examples": data,
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
            / "data/processed"
        ),
    )

    ap.add_argument(
        "--buckets",
        type=int,
        default=128,
        help="Must match the bucket count used by the merger",
    )

    ap.add_argument(
        "--verify-hashes",
        action="store_true",
    )

    ap.add_argument(
        "--verify-source-hashes",
        action="store_true",
    )

    ap.add_argument(
        "--drug",
        default="",
    )

    ap.add_argument(
        "--disease",
        default="",
    )

    ap.add_argument(
        "--examples",
        type=int,
        default=15,
    )

    ap.add_argument(
        "--keep-work",
        action="store_true",
    )

    args = ap.parse_args()

    if not 8 <= args.buckets <= 1024:
        ap.error(
            "--buckets must be between 8 and 1024"
        )

    if not 1 <= args.examples <= 100:
        ap.error(
            "--examples must be between 1 and 100"
        )

    base = (
        args.processed_dir
        .expanduser()
        .resolve()
    )

    inputs = {
        key: base / filename
        for key, filename in FILES.items()
    }

    if not inputs["manifest"].is_file():
        raise RuntimeError(
            f"Missing merge manifest: {inputs['manifest']}"
        )

    manifest = json.loads(
        inputs["manifest"].read_text(
            encoding="utf-8"
        )
    )

    if manifest.get("status") != "success":
        raise RuntimeError(
            "Merge manifest status is not success"
        )

    headers = {
        key: original_header(inputs[key])
        for key in (
            "events",
            "assertions",
            "pairs",
        )
    }

    if (
        headers["events"][:10] != CORE
        or headers["assertions"][:10] != CORE
    ):
        raise RuntimeError(
            "Invalid common-field header order"
        )

    manifest_schema = manifest.get(
        "schema",
        {},
    )

    for file_key, schema_key in (
        ("events", "full"),
        ("assertions", "assertions"),
        ("pairs", "pairs"),
    ):
        if (
            manifest_schema.get(schema_key)
            and headers[file_key]
            != manifest_schema[schema_key]
        ):
            raise RuntimeError(
                f"{file_key}: CSV header differs "
                "from merge manifest"
            )

    for key, required in (
        ("events", CORE + EVENT_EXTRA),
        ("assertions", CORE + ASSERT_EXTRA),
        ("pairs", IDENTITY + PAIR_EXTRA),
    ):
        absent = sorted(
            set(required) - set(headers[key])
        )

        if absent:
            raise RuntimeError(
                f"{key} missing columns: {absent}"
            )

    merge_version = manifest.get(
        "polars_version",
        "unknown",
    )

    if merge_version != pl.__version__:
        raise RuntimeError(
            "Polars version mismatch: "
            f"merge={merge_version}, "
            f"QC={pl.__version__}. "
            "Run QC using the same Polars environment "
            "as the original merge."
        )

    out = (
        base
        / "Analysis_00C_merge_quality_report.json"
    )

    scratch = (
        base
        / ".merge_quality_work"
    )

    if scratch.exists():
        raise RuntimeError(
            f"Scratch directory already exists: {scratch}"
        )

    scratch.mkdir(parents=True)

    checks = []
    warnings = []

    per_source = Counter()
    per_source_row_sums = Counter()
    per_source_row_min = {}
    per_source_row_max = {}

    per_family = Counter()
    per_source_family = Counter()

    namespace_drug = Counter()
    namespace_disease = Counter()

    relationship_counts = Counter()
    source_overlap = Counter()

    total_events = 0
    total_assertions = 0
    total_pairs = 0
    total_conflicts = 0

    total_empty_drug_id = 0
    total_empty_disease_id = 0
    total_invalid_reclassification = 0

    tops = []
    conflicts_examples = []

    def check(name, observed, expected, details=None):
        passed = observed == expected

        checks.append(
            {
                "name": name,
                "pass": passed,
                "observed": observed,
                "expected": expected,
                "details": details,
            }
        )

        if not passed:
            print(
                f"[FAIL] {name}: "
                f"{observed} != {expected}",
                flush=True,
            )

    try:
        print(
            "[1/5] Partitioning evidence and summaries",
            flush=True,
        )

        extra_source_fields = [
            c for c in (
                "clinical_stage_is_explicit_approval",
                "ctd_evidence_class",
                "medi_confidence_tier",
                "extraction_method",
                "status",
            )
            if c in headers["events"]
        ]

        evt = scan(inputs["events"]).select(
            CORE + EVENT_EXTRA + extra_source_fields
        )

        assertions = scan(
            inputs["assertions"]
        ).select(
            CORE + ASSERT_EXTRA
        )

        pairs = scan(
            inputs["pairs"]
        ).select(
            IDENTITY + PAIR_EXTRA
        )

        # Preserve the original merge's empty-string
        # representation when recomputing hash IDs.

        evt = evt.with_columns(
            [
                pl.col(c).fill_null("")
                for c in (
                    CORE + EVENT_EXTRA
                    + extra_source_fields
                )
            ]
        )

        assertions = assertions.with_columns(
            [
                pl.col(c).fill_null("")
                for c in CORE + ASSERT_EXTRA
            ]
        )

        pairs = pairs.with_columns(
            [
                pl.col(c).fill_null("")
                for c in IDENTITY + PAIR_EXTRA
            ]
        )

        evt = evt.with_columns(
            audit_bucket_expr(args.buckets)
        )

        assertions = assertions.with_columns(
            audit_bucket_expr(args.buckets)
        )

        pairs = pairs.with_columns(
            audit_bucket_expr(args.buckets)
        )

        parts = {}

        for label, query in tqdm(
            [
                ("events", evt),
                ("assertions", assertions),
                ("pairs", pairs),
            ],
            desc="Streaming partition steps",
            unit="file",
        ):
            parts[label] = partition(
                query,
                scratch / label,
            )

        all_buckets = sorted(
            set().union(
                *[
                    set(x)
                    for x in parts.values()
                ]
            )
        )

        print(
            f"[2/5] Validating {len(all_buckets)} buckets",
            flush=True,
        )

        for bucket in tqdm(
            all_buckets,
            desc="Cross-checking evidence",
            unit="bucket",
        ):
            group = {}

            for label in (
                "events",
                "assertions",
                "pairs",
            ):
                paths = parts[label].get(
                    bucket,
                    [],
                )

                if not paths:
                    raise RuntimeError(
                        f"Bucket {bucket}: missing {label}"
                    )

                group[label] = pl.read_parquet(
                    paths
                )

                observed = (
                    group[label]["audit_bucket"]
                    .unique()
                    .to_list()
                )

                check(
                    f"bucket {bucket} {label} partition",
                    observed,
                    [bucket],
                )

            e = group["events"]
            a = group["assertions"]
            p = group["pairs"]

            total_events += e.height
            total_assertions += a.height
            total_pairs += p.height

            total_conflicts += int(
                p.select(
                    (
                        pl.col("potential_conflict")
                        .str.to_lowercase()
                        == "true"
                    ).sum()
                ).item()
            )

            # Recompute pair and assertion IDs.

            bad_pair_hash = n_bad(
                e,
                pl.col("merge_pair_id")
                != twohash(IDENTITY, "P"),
            )

            bad_assert_hash = n_bad(
                e,
                pl.col("merge_assertion_id")
                != twohash(ASSERTION, "A"),
            )

            bad_partition = n_bad(
                e,
                pl.col("merge_bucket")
                .cast(pl.Int32)
                != pl.col("audit_bucket"),
            )

            check(
                f"bucket {bucket} pair fingerprint",
                bad_pair_hash,
                0,
            )

            check(
                f"bucket {bucket} assertion fingerprint",
                bad_assert_hash,
                0,
            )

            check(
                f"bucket {bucket} partition consistency",
                bad_partition,
                0,
            )

            # Verify mandatory values.

            for field in (
                "drug_name",
                "disease_or_condition_name",
                "relationship_type",
                "evidence_type",
            ):
                bad = n_bad(
                    e,
                    pl.col(field)
                    .str.strip_chars()
                    == "",
                )

                check(
                    f"bucket {bucket} blank {field}",
                    bad,
                    0,
                )

            check(
                f"bucket {bucket} event IDs unique",
                e["merge_event_id"].n_unique(),
                e.height,
            )

            # Independently verify source event IDs.

            source_code_expr = pl.lit("?")

            for source_idx, source_name in enumerate(
                SOURCES,
                1,
            ):
                source_code_expr = (
                    pl.when(
                        pl.col("source")
                        == source_name
                    )
                    .then(
                        pl.lit(f"{source_idx:02d}")
                    )
                    .otherwise(
                        source_code_expr
                    )
                )

            expected_event_id = pl.concat_str(
                [
                    source_code_expr,
                    pl.lit(":"),
                    pl.col("merge_source_row_number"),
                ],
                separator="",
            )

            check(
                f"bucket {bucket} event ID linkage",
                n_bad(
                    e,
                    pl.col("merge_event_id")
                    != expected_event_id,
                ),
                0,
            )

            check(
                f"bucket {bucket} source labels",
                e.filter(
                    ~pl.col("source").is_in(SOURCES)
                ).height,
                0,
            )

            check(
                f"bucket {bucket} original source present",
                e.filter(
                    pl.col("merge_original_source")
                    .str.strip_chars()
                    == ""
                ).height,
                0,
            )

            # Independent evidence-reclassification checks.

            relation = (
                pl.col("relationship_type")
                .str.to_lowercase()
            )

            family = pl.col(
                "merge_evidence_family"
            )

            bad_class = (
                (
                    (relation == "contraindication")
                    & (family != "contraindication")
                )
                |
                (
                    (pl.col("source") == "AACT")
                    & (family != "trial_cooccurrence")
                )
                |
                (
                    (pl.col("source") == "CTD")
                    & (relation == "ctd_direct_therapeutic")
                    & (family != "ctd_direct_therapeutic")
                )
                |
                (
                    (pl.col("source") == "SIDER")
                    & (relation == "label_nlp_indication")
                    & (family != "label_nlp_indication")
                )
                |
                (
                    (pl.col("source") == "SIDER")
                    & (relation == "label_precondition")
                    & (family != "label_precondition")
                )
                |
                (
                    (pl.col("source") == "repoDB")
                    & (
                        pl.col("evidence_type")
                        == "repodb_approved_indication"
                    )
                    & (
                        family != "source_reported_approval"
                    )
                )
            )

            if (
                "clinical_stage_is_explicit_approval"
                in e.columns
            ):
                bad_class = bad_class | (
                    (
                        pl.col("source")
                        == "OpenTargets"
                    )
                    & (
                        pl.col(
                            "clinical_stage_is_explicit_approval"
                        )
                        .str.to_lowercase()
                        == "true"
                    )
                    & (
                        family != "source_reported_approval"
                    )
                )

            classification_bad = n_bad(
                e,
                bad_class,
            )

            total_invalid_reclassification += (
                classification_bad
            )

            check(
                f"bucket {bucket} evidence reclassification",
                classification_bad,
                0,
            )

            total_empty_drug_id += e.filter(
                pl.col("drug_identifier")
                .str.strip_chars()
                == ""
            ).height

            total_empty_disease_id += e.filter(
                pl.col("disease_or_condition_identifier")
                .str.strip_chars()
                == ""
            ).height

            # Build independent source summaries.

            source_stats = (
                e.with_columns(
                    pl.col("merge_source_row_number")
                    .cast(pl.Int64)
                    .alias("_row")
                )
                .group_by("source")
                .agg(
                    pl.len().alias("n"),
                    pl.col("_row")
                    .sum()
                    .alias("row_sum"),
                    pl.col("_row")
                    .min()
                    .alias("row_min"),
                    pl.col("_row")
                    .max()
                    .alias("row_max"),
                )
                .to_dicts()
            )

            for item in source_stats:
                src = item["source"]

                per_source[src] += int(
                    item["n"]
                )

                per_source_row_sums[src] += int(
                    item["row_sum"]
                )

                per_source_row_min[src] = min(
                    per_source_row_min.get(
                        src,
                        item["row_min"],
                    ),
                    item["row_min"],
                )

                per_source_row_max[src] = max(
                    per_source_row_max.get(
                        src,
                        item["row_max"],
                    ),
                    item["row_max"],
                )

            for item in (
                e.group_by(
                    "merge_evidence_family"
                ).len().to_dicts()
            ):
                per_family[
                    item["merge_evidence_family"]
                ] += int(item["len"])

            for item in (
                e.group_by(
                    [
                        "source",
                        "merge_evidence_family",
                    ]
                ).len().to_dicts()
            ):
                per_source_family[
                    (
                        item["source"],
                        item["merge_evidence_family"],
                    )
                ] += int(item["len"])

            for item in (
                e.group_by(
                    "drug_identifier_type"
                ).len().to_dicts()
            ):
                namespace_drug[
                    item["drug_identifier_type"]
                ] += int(item["len"])

            for item in (
                e.group_by(
                    "disease_or_condition_identifier_type"
                ).len().to_dicts()
            ):
                namespace_disease[
                    item[
                        "disease_or_condition_identifier_type"
                    ]
                ] += int(item["len"])

            for item in (
                e.group_by(
                    "relationship_type"
                ).len().to_dicts()
            ):
                relationship_counts[
                    item["relationship_type"]
                ] += int(item["len"])

            # Recalculate assertions from the events.

            ea = expected_assertions(e)

            check(
                f"bucket {bucket} assertion ID consistency",
                ea.filter(
                    (
                        pl.col("number_assertion_ids")
                        != 1
                    )
                    |
                    (
                        pl.col("number_pair_ids")
                        != 1
                    )
                ).height,
                0,
            )

            check(
                f"bucket {bucket} distinct assertion IDs",
                ea["calc_aid"].n_unique(),
                ea.height,
            )

            check(
                f"bucket {bucket} published assertion count",
                a.height,
                ea.height,
            )

            matches_a = a.join(
                ea,
                on=ASSERTION,
                how="full",
                coalesce=True,
            )

            mism_a = matches_a.filter(
                pl.col("calc_events").is_null()
                | pl.col("supporting_event_count").is_null()
                | (
                    pl.col("supporting_event_count")
                    .cast(pl.Int64)
                    != pl.col("calc_events")
                )
                | (
                    pl.col("merge_assertion_id")
                    != pl.col("calc_aid")
                )
                | (
                    pl.col("merge_pair_id")
                    != pl.col("calc_pid")
                )
            )

            check(
                f"bucket {bucket} assertion evidence support",
                mism_a.height,
                0,
                (
                    mism_a.head(2).to_dicts()
                    if mism_a.height else None
                ),
            )

            check(
                f"bucket {bucket} assertion event conservation",
                int(ea["calc_events"].sum()),
                e.height,
            )

            # Check representative event IDs.

            reps = a.select(
                "merge_assertion_id",
                "representative_event_id",
            ).join(
                e.select(
                    "merge_assertion_id",
                    "merge_event_id",
                ).unique(),
                left_on=[
                    "merge_assertion_id",
                    "representative_event_id",
                ],
                right_on=[
                    "merge_assertion_id",
                    "merge_event_id",
                ],
                how="anti",
            )

            check(
                f"bucket {bucket} representative event links",
                reps.height,
                0,
            )

            # Recalculate six-field pair summaries.

            ep = expected_pairs(e)

            check(
                f"bucket {bucket} pair ID consistency",
                ep.filter(
                    pl.col("number_pair_ids")
                    != 1
                ).height,
                0,
            )

            check(
                f"bucket {bucket} distinct pair IDs",
                ep["calc_pid"].n_unique(),
                ep.height,
            )

            check(
                f"bucket {bucket} published pair count",
                p.height,
                ep.height,
            )

            matches_p = p.join(
                ep,
                on=IDENTITY,
                how="full",
                coalesce=True,
            )

            mism_p = matches_p.filter(
                pl.col("calc_events").is_null()
                | pl.col("supporting_event_count").is_null()
                | (
                    pl.col("supporting_event_count")
                    .cast(pl.Int64)
                    != pl.col("calc_events")
                )
                | (
                    pl.col("supporting_assertion_count")
                    .cast(pl.Int64)
                    != pl.col("calc_assertions")
                )
                | (
                    pl.col("supporting_source_count")
                    .cast(pl.Int64)
                    != pl.col("calc_sources")
                )
                | (
                    pl.col("merge_pair_id")
                    != pl.col("calc_pid")
                )
                | (
                    pl.col("supporting_sources")
                    != pl.col("calc_source_names")
                )
                | (
                    pl.col("relationship_types")
                    != pl.col("calc_relation_names")
                )
                | (
                    pl.col("evidence_types")
                    != pl.col("calc_evidence_names")
                )
                | (
                    pl.col("evidence_families")
                    != pl.col("calc_family_names")
                )
                | (
                    pl.col("indication_like_event_count")
                    .cast(pl.Int64)
                    != pl.col("calc_indication")
                )
                | (
                    pl.col("contraindication_event_count")
                    .cast(pl.Int64)
                    != pl.col("calc_contra")
                )
                | (
                    (
                        pl.col("potential_conflict")
                        .str.to_lowercase()
                        == "true"
                    )
                    != pl.col("calc_conflict")
                )
            )

            check(
                f"bucket {bucket} pair evidence support",
                mism_p.height,
                0,
                (
                    mism_p.head(2).to_dicts()
                    if mism_p.height else None
                ),
            )

            check(
                f"bucket {bucket} pair event conservation",
                int(ep["calc_events"].sum()),
                e.height,
            )

            check(
                f"bucket {bucket} assertion-to-pair conservation",
                int(ep["calc_assertions"].sum()),
                a.height,
            )

            # Distribution of one-source and multi-source pairs.

            for item in (
                p.group_by(
                    "supporting_source_count"
                ).len().to_dicts()
            ):
                source_overlap[
                    int(item["supporting_source_count"])
                ] += int(item["len"])

            # Retain representative high-support examples.

            tops.extend(
                p.with_columns(
                    pl.col("supporting_event_count")
                    .cast(pl.Int64)
                    .alias("_support")
                )
                .sort(
                    "_support",
                    descending=True,
                )
                .select(
                    IDENTITY + [
                        "merge_pair_id",
                        "supporting_event_count",
                        "supporting_sources",
                        "relationship_types",
                    ]
                )
                .head(4)
                .to_dicts()
            )

            # Inspect potential conflicting relationship labels.

            if len(conflicts_examples) < args.examples:
                subset = p.filter(
                    pl.col("potential_conflict")
                    .str.to_lowercase()
                    == "true"
                )

                conflicts_examples.extend(
                    subset.select(
                        IDENTITY + [
                            "merge_pair_id",
                            "supporting_sources",
                            "relationship_types",
                            "supporting_event_count",
                        ]
                    )
                    .head(
                        args.examples
                        - len(conflicts_examples)
                    )
                    .to_dicts()
                )

        print(
            "[3/5] Global conservation and hash checks",
            flush=True,
        )

        check(
            "global events",
            total_events,
            int(manifest["source_events"]),
        )

        check(
            "global assertions",
            total_assertions,
            int(manifest["assertions"]),
        )

        check(
            "global exact pairs",
            total_pairs,
            int(manifest["exact_identity_pairs"]),
        )

        check(
            "global potential conflicts",
            total_conflicts,
            int(
                manifest["potential_conflicting_pairs"]
            ),
        )

        check(
            "all eleven sources represented",
            sorted(per_source),
            sorted(SOURCES),
        )

        declared = {
            x["source"]: int(x["observed"])
            for x in manifest.get("sources", [])
        }

        check(
            "per-source counts agree with manifest",
            dict(per_source),
            declared,
        )

        for src, n in per_source.items():
            check(
                f"{src} row indices begin at one",
                per_source_row_min[src],
                1,
            )

            check(
                f"{src} row indices reach source size",
                per_source_row_max[src],
                n,
            )

            check(
                f"{src} arithmetic row-index conservation",
                per_source_row_sums[src],
                n * (n + 1) // 2,
            )

        check(
            "evidence family counts conserve rows",
            sum(per_family.values()),
            total_events,
        )

        check(
            "source-family counts conserve rows",
            sum(per_source_family.values()),
            total_events,
        )

        check(
            "pair groups conserve rows",
            sum(source_overlap.values()),
            total_pairs,
        )

        check(
            "tested evidence reclassifications",
            total_invalid_reclassification,
            0,
        )

        if args.verify_hashes:
            for key, filename_key in (
                ("events", "merged"),
                ("assertions", "dedup"),
                ("pairs", "pairs"),
            ):
                expected = (
                    manifest.get(
                        "output_sha256",
                        {},
                    ).get(filename_key)
                )

                if not expected:
                    check(
                        f"SHA256 {key} present",
                        "missing",
                        "available",
                    )
                else:
                    check(
                        f"SHA256 {key}",
                        checksum(inputs[key]),
                        expected,
                    )

        if args.verify_source_hashes:
            for spec in manifest.get(
                "sources",
                [],
            ):
                path = base / spec["path"]

                check(
                    f"source checksum {spec['source']}",
                    (
                        checksum(path)
                        if path.is_file()
                        else "missing"
                    ),
                    spec["sha256"],
                )

        # Scientifically important warnings.

        if (
            source_overlap.get(2, 0)
            + sum(
                count
                for n, count in source_overlap.items()
                if n >= 3
            )
        ) == 0:
            warnings.append(
                "No exact cross-source pairs: "
                "investigate name and namespace fragmentation"
            )

        if total_conflicts != 87:
            warnings.append(
                "Conflict count differs from the supplied "
                "87-conflict merger run"
            )

        if per_family.get(
            "other_or_unclassified",
            0,
        ):
            warnings.append(
                "Unclassified evidence families need review"
            )

        if (
            total_empty_drug_id
            or total_empty_disease_id
        ):
            warnings.append(
                "Some original external identifiers are blank; "
                "do not claim complete ontology mapping"
            )

        print(
            "[4/5] Querying selected drug/disease",
            flush=True,
        )

        query = query_pairs(
            inputs["pairs"],
            args.drug,
            args.disease,
            args.examples,
        )

        if query:
            print(
                json.dumps(
                    query,
                    indent=2,
                    ensure_ascii=False,
                )[:16000]
            )

        passed = all(
            result["pass"]
            for result in checks
        )

        report = {
            "generated_at_utc": timestamp(),

            "verdict": (
                "PASS_WITH_WARNINGS"
                if passed and warnings
                else "PASS"
                if passed
                else "FAIL"
            ),

            "scientific_lock_recommendation":
                "hold_for_scientific_review",

            "notes": [
                "Technical consistency does not prove clinical validity.",
                "Exact identities do not resolve heterogeneous ontologies.",
                "Reported approval is not independent regulatory verification.",
                "Conflicts are exact-key label conflicts only.",
                "Database recurrence is not necessarily independent evidence.",
            ],

            "polars_version": pl.__version__,
            "merge_polars_version": merge_version,

            "events": total_events,
            "assertions": total_assertions,
            "exact_identity_pairs": total_pairs,
            "potential_conflicts": total_conflicts,

            "known_reclassification_mismatches":
                total_invalid_reclassification,

            "empty_drug_ids":
                total_empty_drug_id,

            "empty_disease_ids":
                total_empty_disease_id,

            "source_counts":
                dict(per_source),

            "evidence_family_counts":
                dict(per_family),

            "source_evidence_family_counts": {
                f"{src}::{family}": n
                for (src, family), n
                in per_source_family.items()
            },

            "relationship_counts":
                dict(relationship_counts),

            "drug_identifier_namespace_counts":
                dict(namespace_drug),

            "disease_identifier_namespace_counts":
                dict(namespace_disease),

            "exact_pairs_by_number_of_supporting_sources": {
                str(n): count
                for n, count
                in sorted(source_overlap.items())
            },

            "top_supported_exact_pairs": sorted(
                tops,
                key=lambda x: int(
                    x["supporting_event_count"]
                ),
                reverse=True,
            )[:args.examples],

            "potential_conflict_examples":
                conflicts_examples,

            "query": query,

            "warnings": warnings,

            "checks_run": len(checks),

            "checks_failed": sum(
                not c["pass"]
                for c in checks
            ),

            "failed_checks": [
                c for c in checks
                if not c["pass"]
            ][:100],
        }

        out.write_text(
            json.dumps(
                report,
                indent=2,
                ensure_ascii=False,
            ) + "\n",
            encoding="utf-8",
        )

        print(
            "[5/5] Final verdict:",
            report["verdict"],
        )

        print(
            "[CHECKS]",
            report["checks_run"],
            "failed:",
            report["checks_failed"],
        )

        print(
            "[EVENTS]",
            total_events,
            "[ASSERTIONS]",
            total_assertions,
            "[PAIRS]",
            total_pairs,
            "[CONFLICTS]",
            total_conflicts,
        )

        print(
            "[SOURCE COUNTS]",
            json.dumps(
                dict(per_source),
                indent=2,
            ),
        )

        print(
            "[EVIDENCE FAMILIES]",
            json.dumps(
                dict(per_family),
                indent=2,
            ),
        )

        print(
            "[CROSS-SOURCE SUPPORT]",
            json.dumps(
                dict(source_overlap),
                indent=2,
            ),
        )

        print(
            "[WARNINGS]",
            json.dumps(
                warnings,
                indent=2,
            ),
        )

        print(
            "[TOP EXACT PAIRS]",
            json.dumps(
                report["top_supported_exact_pairs"][:8],
                indent=2,
                ensure_ascii=False,
            ),
        )

        print(
            "[EXAMPLE POTENTIAL CONFLICTS]",
            json.dumps(
                conflicts_examples[:8],
                indent=2,
                ensure_ascii=False,
            ),
        )

        print("[REPORT]", out)

        return 0 if passed else 2

    finally:
        if args.keep_work:
            print(
                "[SCRATCH RETAINED]",
                scratch,
            )
        else:
            shutil.rmtree(
                scratch,
                ignore_errors=True,
            )


if __name__ == "__main__":
    try:
        sys.exit(main())

    except Exception as exc:
        print(
            f"[QC ERROR] "
            f"{type(exc).__name__}: {exc}",
            file=sys.stderr,
        )
        raise
