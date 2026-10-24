#!/usr/bin/env python3
"""Analysis 09: manuscript-ready pairwise drug--disease name-pair overlap heatmap.

Uses the same lowercase+whitespace-only matching as Analyses 01, 04 and 08.

Reads the Analysis 08 pairwise TSV if available. Otherwise independently
reconstructs all 55 overlaps from the locked assertion CSV using streamed,
hash-partitioned Polars processing.

Never alters the locked datasets.

Outputs:
    figures/pairwise_drug_disease_pair_overlap_heatmap.png
    figures/pairwise_drug_disease_pair_overlap_heatmap.pdf
    tables/pairwise_drug_disease_pair_overlap_heatmap.csv

Requirements:
    pip install polars numpy matplotlib tqdm

Run:
    POLARS_MAX_THREADS=16 python Analysis_09_Pairwise_drug_disease_pair_overlap_heatmap.py
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import shutil

from collections import Counter, defaultdict
from itertools import combinations
from pathlib import Path

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import polars as pl

from tqdm import tqdm


# ============================================================
# CONFIGURATION
# ============================================================

SOURCES = [
    "AACT",
    "Broad Drug Repurposing Hub",
    "CTD",
    "ChEMBL",
    "DrugCentral",
    "Hetionet",
    "MEDI",
    "OpenTargets",
    "PrimeKG",
    "SIDER",
    "repoDB",
]

SHORT = {
    "AACT": "AACT",
    "Broad Drug Repurposing Hub": "Broad Hub",
    "CTD": "CTD",
    "ChEMBL": "ChEMBL",
    "DrugCentral": "DrugCentral",
    "Hetionet": "Hetionet",
    "MEDI": "MEDI",
    "OpenTargets": "Open Targets",
    "PrimeKG": "PrimeKG",
    "SIDER": "SIDER",
    "repoDB": "repoDB",
}

BITS = {
    source: 1 << i
    for i, source in enumerate(SOURCES)
}

# Locked results from Analysis 08
EXPECTED_SOURCE_PAIRS = {
    "AACT": 466990,
    "Broad Drug Repurposing Hub": 3885,
    "CTD": 3529952,
    "ChEMBL": 59681,
    "DrugCentral": 42476,
    "Hetionet": 1145,
    "MEDI": 206929,
    "OpenTargets": 53942,
    "PrimeKG": 42383,
    "SIDER": 16875,
    "repoDB": 10950,
}

EXPECTED_TOTAL_PAIRS = 4_313_189
EXPECTED_MULTI_SOURCE_PAIRS = 83_878

EXPECTED_ASSERTIONS = 4_534_385
EXPECTED_EVENTS = 10_895_395

EXPECTED_PAIRWISE_INTERSECTION = 43_401


# ============================================================
# INPUT VALIDATION
# ============================================================

def ensure_schema(path, needed):

    if not path.is_file() or path.stat().st_size == 0:
        raise RuntimeError(
            f"Required file is missing or empty: {path}"
        )

    with path.open(
        "r",
        encoding="utf-8-sig",
        newline="",
    ) as handle:
        cols = next(
            csv.reader(handle, strict=True),
            [],
        )

    if (
        len(cols) != len(set(cols))
        or not set(needed).issubset(cols)
    ):
        raise RuntimeError(
            f"Unexpected CSV schema in {path}"
        )


# ============================================================
# NAME NORMALISATION
# ============================================================

def normalized_expr(col):

    return (
        pl.col(col)
        .fill_null("")
        .str.strip_chars()
        .str.replace_all(r"\s+", " ")
        .str.to_lowercase()
    )


def source_bit_expr():

    expr = pl.lit(0, dtype=pl.Int32)

    for source, bit in BITS.items():
        expr = (
            pl.when(pl.col("source") == source)
            .then(pl.lit(bit))
            .otherwise(expr)
        )

    return expr.alias("_bit")


# ============================================================
# COMPUTE OVERLAPS FROM SOURCE MEMBERSHIP
# ============================================================

def overlaps_from_masks(histogram):

    source_size = {
        source: sum(
            n
            for mask, n in histogram.items()
            if mask & bit
        )
        for source, bit in BITS.items()
    }

    n_pairs = sum(histogram.values())

    n_multi = sum(
        n
        for mask, n in histogram.items()
        if mask.bit_count() > 1
    )

    if (
        n_pairs != EXPECTED_TOTAL_PAIRS
        or n_multi != EXPECTED_MULTI_SOURCE_PAIRS
    ):
        raise RuntimeError(
            "Name-pair regression mismatch: "
            f"{n_pairs:,}; multi={n_multi:,}"
        )

    if source_size != EXPECTED_SOURCE_PAIRS:
        raise RuntimeError(
            "Per-source coverage does not match "
            "the locked Analysis 08 results"
        )

    rows = []

    for a, b in combinations(SOURCES, 2):

        na = source_size[a]
        nb = source_size[b]

        common = sum(
            n
            for mask, n in histogram.items()
            if (
                (mask & BITS[a])
                and (mask & BITS[b])
            )
        )

        union = na + nb - common

        rows.append({
            "source_a": a,
            "source_b": b,
            "a_unique_keys": na,
            "b_unique_keys": nb,
            "intersection": common,
            "union": union,
            "jaccard": (
                common / union
                if union else 0.0
            ),
            "overlap_coefficient": (
                common / min(na, nb)
                if min(na, nb) else 0.0
            ),
            "pct_source_a": (
                100 * common / na
                if na else 0.0
            ),
            "pct_source_b": (
                100 * common / nb
                if nb else 0.0
            ),
        })

    return rows


# ============================================================
# POLARS RECOMPUTATION
# ============================================================

def recompute_pairwise(
    assertions,
    manifest,
    work,
    buckets,
):

    ensure_schema(
        assertions,
        [
            "source",
            "drug_name",
            "disease_or_condition_name",
            "supporting_event_count",
        ],
    )

    work.mkdir(
        parents=True,
        exist_ok=False,
    )

    print(
        "[STEP 1] Streaming assertion records",
        flush=True,
    )

    narrow_file = work / "assertions.parquet"

    lf = (
        pl.scan_csv(
            assertions,
            infer_schema=False,
            ignore_errors=False,
            truncate_ragged_lines=False,
            encoding="utf8",
        )
        .select([
            "source",
            "drug_name",
            "disease_or_condition_name",
            "supporting_event_count",
        ])
        .with_columns([
            normalized_expr("drug_name").alias(
                "_drug"
            ),

            normalized_expr(
                "disease_or_condition_name"
            ).alias("_disease"),

            pl.col("supporting_event_count")
            .cast(pl.Int64, strict=True),

            source_bit_expr(),
        ])
        .select([
            "_drug",
            "_disease",
            "_bit",
            "supporting_event_count",
        ])
    )

    lf.sink_parquet(
        narrow_file,
        engine="streaming",
        compression="zstd",
    )

    data = pl.scan_parquet(narrow_file)

    stats = (
        data.select([
            pl.len().alias("assertions"),

            pl.col("supporting_event_count")
            .sum()
            .alias("events"),

            (pl.col("_bit") == 0)
            .sum()
            .alias("unknown_source"),

            (
                (pl.col("_drug") == "")
                |
                (pl.col("_disease") == "")
            )
            .sum()
            .alias("blank_key"),
        ])
        .collect(engine="streaming")
        .row(0, named=True)
    )

    if (
        stats["assertions"] != EXPECTED_ASSERTIONS
        or stats["events"] != EXPECTED_EVENTS
        or stats["unknown_source"]
        or stats["blank_key"]
    ):
        raise RuntimeError(
            f"Assertion conservation failed: {stats}"
        )

    if (
        int(manifest["assertions"])
        != stats["assertions"]
        or int(manifest["source_events"])
        != stats["events"]
    ):
        raise RuntimeError(
            "Merge manifest disagrees with assertions"
        )

    print(
        "[STEP 2] Partitioning normalised pairs",
        flush=True,
    )

    bucket_dir = work / "buckets"

    (
        data.with_columns(
            (
                pl.struct([
                    "_drug",
                    "_disease",
                ]).hash(seed=202604)
                % buckets
            )
            .cast(pl.Int32)
            .alias("_bucket")
        )
        .sink_parquet(
            pl.PartitionBy(
                bucket_dir,
                key="_bucket",
                include_key=True,
            ),
            engine="streaming",
            compression="zstd",
            mkdir=True,
        )
    )

    paths = defaultdict(list)

    for file in bucket_dir.rglob("*.parquet"):

        bucket = int(
            pl.read_parquet(
                file,
                columns=["_bucket"],
                n_rows=1,
            )["_bucket"][0]
        )

        paths[bucket].append(file)

    if not paths:
        raise RuntimeError(
            "Polars partitioning produced no files"
        )

    histogram = Counter()

    assertions_recounted = 0
    events_recounted = 0

    print(
        "[STEP 3] Aggregate source overlap",
        flush=True,
    )

    for bucket, group_files in tqdm(
        sorted(paths.items()),
        desc="Overlap partitions",
        unit="bucket",
    ):

        block = pl.read_parquet(group_files)

        assertions_recounted += block.height

        events_recounted += int(
            block["supporting_event_count"].sum()
        )

        grouped = (
            block.group_by([
                "_drug",
                "_disease",
            ])
            .agg(
                pl.col("_bit")
                .unique()
                .sum()
                .alias("mask")
            )
        )

        histogram.update({
            int(mask): int(count)
            for mask, count in (
                grouped.group_by("mask")
                .len()
                .iter_rows()
            )
        })

    if (
        assertions_recounted != EXPECTED_ASSERTIONS
        or events_recounted != EXPECTED_EVENTS
    ):
        raise RuntimeError(
            "Partition-level conservation failed"
        )

    return overlaps_from_masks(histogram)


# ============================================================
# LOAD PRECOMPUTED ANALYSIS 08 RESULTS
# ============================================================

def load_saved_pairs(path):

    fields = [
        "source_a",
        "source_b",
        "a_unique_keys",
        "b_unique_keys",
        "intersection",
        "union",
        "jaccard",
        "overlap_coefficient",
        "pct_source_a",
        "pct_source_b",
    ]

    if not path.is_file():
        return None

    print(
        f"[REUSE] Analysis 08 table: {path}",
        flush=True,
    )

    rows = []

    with path.open(
        "r",
        encoding="utf-8-sig",
        newline="",
    ) as handle:

        reader = csv.DictReader(
            handle,
            delimiter="\t",
        )

        if (
            not reader.fieldnames
            or not set(fields).issubset(
                reader.fieldnames
            )
        ):
            raise RuntimeError(
                "Analysis 08 TSV has incompatible columns"
            )

        for r in reader:

            row = {}

            for c in fields:

                if c in ("source_a", "source_b"):
                    row[c] = r[c]

                elif c in (
                    "a_unique_keys",
                    "b_unique_keys",
                    "intersection",
                    "union",
                ):
                    row[c] = int(r[c])

                else:
                    row[c] = float(r[c])

            rows.append(row)

    return rows


# ============================================================
# INDEPENDENT VALIDATION
# ============================================================

def verify_pairwise(rows):

    if len(rows) != 55:
        raise RuntimeError(
            f"Expected 55 source comparisons; "
            f"found {len(rows)}"
        )

    expected_pairs = {
        frozenset(x)
        for x in combinations(SOURCES, 2)
    }

    observed = {
        frozenset((
            r["source_a"],
            r["source_b"],
        ))
        for r in rows
    }

    if (
        observed != expected_pairs
        or len(observed) != len(rows)
    ):
        raise RuntimeError(
            "Missing or duplicate source comparisons"
        )

    for r in rows:

        a = r["source_a"]
        b = r["source_b"]

        na = r["a_unique_keys"]
        nb = r["b_unique_keys"]
        k = r["intersection"]

        if (
            na != EXPECTED_SOURCE_PAIRS[a]
            or nb != EXPECTED_SOURCE_PAIRS[b]
        ):
            raise RuntimeError(
                f"Source count mismatch: {a} vs {b}"
            )

        if (
            not 0 <= k <= min(na, nb)
            or r["union"] != na + nb - k
        ):
            raise RuntimeError(
                f"Invalid intersection: {a} vs {b}"
            )

        expected_jaccard = k / (na + nb - k)

        if not math.isclose(
            float(r["jaccard"]),
            expected_jaccard,
            abs_tol=1e-6,
        ):
            raise RuntimeError(
                f"Jaccard mismatch: {a} vs {b}"
            )

    selected = next(
        r for r in rows
        if {
            r["source_a"],
            r["source_b"],
        } == {
            "ChEMBL",
            "OpenTargets",
        }
    )

    if (
        selected["intersection"]
        != EXPECTED_PAIRWISE_INTERSECTION
    ):
        raise RuntimeError(
            "ChEMBL–OpenTargets overlap mismatch"
        )

    print(
        "[PASS] 11 sources and 55 comparisons",
        flush=True,
    )

    print(
        "[PASS] ChEMBL–OpenTargets: "
        f"{selected['intersection']:,} shared pairs; "
        f"Jaccard={selected['jaccard']:.6f}",
        flush=True,
    )


# ============================================================
# PUBLICATION-QUALITY HEATMAP
# ============================================================

def plot_heatmap(
    rows,
    png_file,
    pdf_file,
):

    n = len(SOURCES)

    counts = np.zeros(
        (n, n),
        dtype=np.int64,
    )

    jaccard_pct = np.full(
        (n, n),
        np.nan,
        dtype=float,
    )

    for i, name in enumerate(SOURCES):
        counts[i, i] = EXPECTED_SOURCE_PAIRS[name]

    indices = {
        source: i
        for i, source in enumerate(SOURCES)
    }

    for r in rows:

        i = indices[r["source_a"]]
        j = indices[r["source_b"]]

        counts[i, j] = r["intersection"]
        counts[j, i] = r["intersection"]

        value = 100 * r["jaccard"]

        jaccard_pct[i, j] = value
        jaccard_pct[j, i] = value

    vmax = float(
        np.nanmax(jaccard_pct)
    )

    cmap = plt.get_cmap("YlOrRd").copy()

    cmap.set_bad("#344054")

    fig, ax = plt.subplots(
        figsize=(14.8, 12.2)
    )

    img = ax.imshow(
        jaccard_pct,
        vmin=0,
        vmax=vmax,
        cmap=cmap,
        interpolation="nearest",
        aspect="equal",
    )

    # --------------------------------------------------------
    # AXIS LABELS
    # --------------------------------------------------------

    ax.set_xticks(
        range(n),
        labels=[
            SHORT[s]
            for s in SOURCES
        ],
        rotation=44,
        ha="right",
    )

    ax.set_yticks(
        range(n),
        labels=[
            SHORT[s]
            for s in SOURCES
        ],
    )

    ax.tick_params(
        axis="both",
        which="major",
        length=0,
        labelsize=10.5,
    )

    # --------------------------------------------------------
    # CELL BORDERS
    # --------------------------------------------------------

    ax.set_xticks(
        np.arange(-0.5, n, 1),
        minor=True,
    )

    ax.set_yticks(
        np.arange(-0.5, n, 1),
        minor=True,
    )

    ax.grid(
        which="minor",
        color="white",
        linewidth=1.35,
    )

    ax.tick_params(
        which="minor",
        bottom=False,
        left=False,
    )

    for spine in ax.spines.values():
        spine.set_visible(False)

    # --------------------------------------------------------
    # CELL ANNOTATIONS
    # --------------------------------------------------------

    for i in range(n):
        for j in range(n):

            count = f"{counts[i, j]:,}"

            if i == j:

                color = "white"
                sub = "source size"

            else:

                percentage = jaccard_pct[i, j]

                color = (
                    "white"
                    if percentage / vmax >= 0.62
                    else "#253043"
                )

                sub = f"{percentage:.2f}%"

            ax.text(
                j,
                i - 0.13,
                count,
                fontsize=8.0,
                weight="semibold",
                color=color,
                ha="center",
                va="center",
            )

            ax.text(
                j,
                i + 0.17,
                sub,
                fontsize=7.7,
                color=color,
                ha="center",
                va="center",
            )

    # --------------------------------------------------------
    # TITLE
    # --------------------------------------------------------

    ax.set_title(
        "Pairwise overlap of drug–disease name pairs "
        "across 11 databases",
        fontweight="bold",
        fontsize=16,
        pad=23,
    )

    # --------------------------------------------------------
    # COLOUR BAR
    # --------------------------------------------------------

    cbar = fig.colorbar(
        img,
        ax=ax,
        fraction=0.034,
        pad=0.03,
        shrink=0.86,
    )

    cbar.set_label(
        "Jaccard similarity (%)",
        fontsize=11,
    )

    cbar.ax.tick_params(
        labelsize=9
    )

    # --------------------------------------------------------
    # FIGURE FOOTNOTE
    # --------------------------------------------------------

    fig.text(
        0.5,
        0.035,
        "Each off-diagonal cell: shared name pairs (top) "
        "and Jaccard similarity (bottom). "
        "Diagonal: source pair count.\n"
        "Names are matched case-insensitively after "
        "whitespace normalisation; overlap does not "
        "imply independent evidence or therapeutic approval.",
        ha="center",
        va="center",
        fontsize=9,
        color="#475467",
    )

    fig.subplots_adjust(
        left=0.13,
        right=0.88,
        bottom=0.16,
        top=0.90,
    )

    # --------------------------------------------------------
    # SAVE FIGURES
    # --------------------------------------------------------

    fig.savefig(
        png_file,
        format="png",
        dpi=300,
        bbox_inches="tight",
        facecolor="white",
    )

    fig.savefig(
        pdf_file,
        format="pdf",
        bbox_inches="tight",
        facecolor="white",
    )

    plt.close(fig)

    return vmax


# ============================================================
# MAIN
# ============================================================

def main():

    ap = argparse.ArgumentParser(
        description=__doc__
    )

    root = Path(__file__).resolve().parent

    ap.add_argument(
        "--processed-dir",
        type=Path,
        default=root / "data" / "processed",
    )

    ap.add_argument(
        "--output-dir",
        type=Path,
        default=root,
    )

    ap.add_argument(
        "--buckets",
        type=int,
        default=128,
    )

    ap.add_argument(
        "--recompute",
        action="store_true",
        help=(
            "Recalculate overlaps from locked "
            "assertions instead of reading Analysis 08 TSV"
        ),
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

    base = (
        args.processed_dir
        .expanduser()
        .resolve()
    )

    root = (
        args.output_dir
        .expanduser()
        .resolve()
    )

    manifest_file = (
        base
        / "ALL_INCLUDED_SOURCES_merge_summary.json"
    )

    if not manifest_file.is_file():
        raise RuntimeError(
            f"Missing locked manifest: {manifest_file}"
        )

    manifest = json.loads(
        manifest_file.read_text(
            encoding="utf8"
        )
    )

    if (
        manifest.get("status") != "success"
        or int(manifest.get("assertions", 0))
        != EXPECTED_ASSERTIONS
        or int(manifest.get("source_events", 0))
        != EXPECTED_EVENTS
    ):
        raise RuntimeError(
            "Locked merge manifest does not reconcile"
        )

    work = (
        base
        / ".analysis09_heatmap_work"
    )

    saved_tsv = (
        base
        / "analysis_08_cross_source_overlap_namekeys_v2"
        / "pair_pairwise.tsv"
    )

    if work.exists():
        raise RuntimeError(
            f"Old scratch directory exists: {work}"
        )

    try:

        # ----------------------------------------------------
        # LOAD OR RECOMPUTE OVERLAPS
        # ----------------------------------------------------

        rows = (
            None
            if args.recompute
            else load_saved_pairs(saved_tsv)
        )

        if rows is None:

            print(
                "[INFO] Analysis 08 TSV unavailable; "
                "recomputing from locked assertions",
                flush=True,
            )

            rows = recompute_pairwise(
                base
                / "ALL_INCLUDED_SOURCES_drug_disease_deduplicated.csv",
                manifest,
                work,
                args.buckets,
            )

        verify_pairwise(rows)

        # ----------------------------------------------------
        # OUTPUT DIRECTORIES
        # ----------------------------------------------------

        figdir = root / "figures"
        tabdir = root / "tables"

        figdir.mkdir(
            parents=True,
            exist_ok=True,
        )

        tabdir.mkdir(
            parents=True,
            exist_ok=True,
        )

        png = (
            figdir
            / "pairwise_drug_disease_pair_overlap_heatmap.png"
        )

        pdf = (
            figdir
            / "pairwise_drug_disease_pair_overlap_heatmap.pdf"
        )

        csv_file = (
            tabdir
            / "pairwise_drug_disease_pair_overlap_heatmap.csv"
        )

        png_tmp = png.with_name(
            png.name + ".tmp"
        )

        pdf_tmp = pdf.with_name(
            pdf.name + ".tmp"
        )

        csv_tmp = csv_file.with_name(
            csv_file.name + ".tmp"
        )

        # ----------------------------------------------------
        # WRITE PAIRWISE SUMMARY
        # ----------------------------------------------------

        fields = [
            "source_a",
            "source_b",
            "a_unique_keys",
            "b_unique_keys",
            "intersection",
            "union",
            "jaccard",
            "overlap_coefficient",
            "pct_source_a",
            "pct_source_b",
        ]

        with csv_tmp.open(
            "w",
            encoding="utf8",
            newline="",
        ) as output:

            writer = csv.DictWriter(
                output,
                fieldnames=fields,
            )

            writer.writeheader()

            writer.writerows(
                sorted(
                    rows,
                    key=lambda r: (
                        r["source_a"],
                        r["source_b"],
                    ),
                )
            )

        # ----------------------------------------------------
        # GENERATE FIGURE
        # ----------------------------------------------------

        vmax = plot_heatmap(
            rows,
            png_tmp,
            pdf_tmp,
        )

        # ----------------------------------------------------
        # PUBLISH RESULTS
        # ----------------------------------------------------

        for staged, final in (
            (png_tmp, png),
            (pdf_tmp, pdf),
            (csv_tmp, csv_file),
        ):
            os.replace(staged, final)

        # ----------------------------------------------------
        # PRINT ALL 55 COMPARISONS
        # ----------------------------------------------------

        print(
            "\nPAIRWISE NAME-PAIR OVERLAP: "
            "ALL 55 SOURCE COMPARISONS"
        )

        print(
            "source_a\tsource_b\tshared_pairs\t"
            "jaccard_pct\ta_pairs\tb_pairs"
        )

        for r in sorted(
            rows,
            key=lambda x: x["intersection"],
            reverse=True,
        ):

            print(
                f"{r['source_a']}\t"
                f"{r['source_b']}\t"
                f"{r['intersection']:,}\t"
                f"{100*r['jaccard']:.3f}\t"
                f"{r['a_unique_keys']:,}\t"
                f"{r['b_unique_keys']:,}"
            )

        print(
            f"\n[RESULT] Highest off-diagonal Jaccard: "
            f"{vmax:.3f}%"
        )

        print(f"[SAVED] {png}")
        print(f"[SAVED] {pdf}")
        print(f"[SAVED] {csv_file}")

        print(
            "[PASS] Analysis 09 complete: "
            "11 databases; 55 verified comparisons"
        )

        print(
            "[NOTE] Shared name pairs do not necessarily "
            "represent independent therapeutic evidence."
        )

    finally:

        if work.exists():

            if args.keep_work:
                print(
                    f"[SCRATCH RETAINED] {work}"
                )

            else:
                shutil.rmtree(work)


if __name__ == "__main__":
    main()