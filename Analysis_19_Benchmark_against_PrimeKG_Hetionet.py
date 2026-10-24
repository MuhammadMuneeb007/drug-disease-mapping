#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Analysis 19: Benchmark locked 11-source drug-disease data vs PrimeKG/Hetionet.

Computes name-key coverage and increments from the LOCKED source-assertion data,
plus a separate indication-oriented sensitivity check using original event-level
merge_evidence_family. Does not modify locked data or interpret recurrence as
independent therapeutic validation. Optional output: tables/Analysis_19_*.csv.

Requires: polars, tqdm
Run: POLARS_MAX_THREADS=16 python Analysis_19_Benchmark_against_PrimeKG_Hetionet.py --save-tables
"""
from __future__ import annotations

import argparse
import csv
import json
import shutil
from collections import Counter, defaultdict
from pathlib import Path

import polars as pl
from tqdm import tqdm

SOURCES = [
    "AACT", "Broad Drug Repurposing Hub", "CTD", "ChEMBL", "DrugCentral",
    "Hetionet", "MEDI", "OpenTargets", "PrimeKG", "SIDER", "repoDB",
]
BITS = {name: 1 << i for i, name in enumerate(SOURCES)}
ALL_MASK = sum(BITS.values())
COMPARE = ("PrimeKG", "Hetionet")
COMPARATOR_MASK = BITS["PrimeKG"] | BITS["Hetionet"]
OTHER_MASK = ALL_MASK ^ COMPARATOR_MASK
NON_CTD_MASK = ALL_MASK ^ BITS["CTD"]

EXPECTED = {
    "assertions": 4_534_385,
    "events": 10_895_395,
    "drug_names": 158_364,
    "disease_names": 72_563,
    "pairs": 4_313_189,
    "multi_source_pairs": 83_878,
    "primekg_pairs": 42_383,
    "hetionet_pairs": 1_145,
    "prime_het_overlap": 331,
}
SOURCE_PAIRS_08 = {
    "AACT": 466_990, "Broad Drug Repurposing Hub": 3_885,
    "CTD": 3_529_952, "ChEMBL": 59_681, "DrugCentral": 42_476,
    "Hetionet": 1_145, "MEDI": 206_929, "OpenTargets": 53_942,
    "PrimeKG": 42_383, "SIDER": 16_875, "repoDB": 10_950,
}
FAMILIES = [
    "source_reported_approval",
    "source_reported_indication",
    "ctd_direct_therapeutic",
    "label_nlp_indication",
]
ASSERT_FILE = "ALL_INCLUDED_SOURCES_drug_disease_deduplicated.csv"
EVENT_FILE = "ALL_INCLUDED_SOURCES_drug_disease_merged.csv"
MANIFEST_FILE = "ALL_INCLUDED_SOURCES_merge_summary.json"


def check_file(path: Path, required: list[str]) -> None:
    if not path.is_file() or not path.stat().st_size:
        raise RuntimeError(f"Missing or empty input: {path}")
    with path.open("r", encoding="utf-8-sig", newline="") as f:
        columns = next(csv.reader(f, strict=True), [])
    if len(columns) != len(set(columns)) or not set(required).issubset(columns):
        raise RuntimeError(
            f"Unexpected schema: {path}; missing={sorted(set(required) - set(columns))}"
        )


def scan(path: Path, columns: list[str]) -> pl.LazyFrame:
    return pl.scan_csv(
        str(path), infer_schema=False, ignore_errors=False,
        truncate_ragged_lines=False, encoding="utf8",
    ).select(columns)


def key(column: str) -> pl.Expr:
    # Precisely the name-key definition used in Analyses 01, 04, 08 and 09.
    return (pl.col(column).fill_null("").str.strip_chars()
            .str.replace_all(r"\s+", " ").str.to_lowercase())


def source_bit() -> pl.Expr:
    expr = pl.lit(0, dtype=pl.Int32)
    for name, bit in BITS.items():
        expr = pl.when(pl.col("source") == name).then(pl.lit(bit)).otherwise(expr)
    return expr.alias("_bit")


def eq(label, observed, expected) -> None:
    if observed != expected:
        raise RuntimeError(f"[FAIL] {label}: observed={observed!r}, expected={expected!r}")
    detail = f"{observed:,}" if isinstance(observed, int) else "OK"
    print(f"[PASS] {label}: {detail}", flush=True)


def show(title: str, frame: pl.DataFrame) -> None:
    print("\n" + "=" * 106)
    print(f"{title} [{frame.height:,} rows]")
    print("=" * 106)
    if not frame.height:
        print("[EMPTY]")
        return
    print("\t".join(frame.columns))
    for record in frame.iter_rows():
        print("\t".join(
            "" if v is None else f"{v:.3f}" if isinstance(v, float)
            else str(v).replace("\n", " ").replace("\r", " ").replace("\t", " ")
            for v in record
        ))


def count_where(hist: Counter, predicate) -> int:
    return sum(n for mask, n in hist.items() if predicate(mask))


def pair_totals(hist: Counter) -> dict[str, int]:
    pk, ht = BITS["PrimeKG"], BITS["Hetionet"]
    return {
        "integrated": sum(hist.values()),
        "primekg": count_where(hist, lambda m: bool(m & pk)),
        "hetionet": count_where(hist, lambda m: bool(m & ht)),
        "both": count_where(hist, lambda m: (m & COMPARATOR_MASK) == COMPARATOR_MASK),
        "comparator_union": count_where(hist, lambda m: bool(m & COMPARATOR_MASK)),
        "outside_both": count_where(hist, lambda m: not (m & COMPARATOR_MASK)),
        "noncomparator_union": count_where(hist, lambda m: bool(m & OTHER_MASK)),
        "non_ctd_union": count_where(hist, lambda m: bool(m & NON_CTD_MASK)),
        "outside_both_non_ctd": count_where(
            hist, lambda m: not (m & COMPARATOR_MASK) and bool(m & NON_CTD_MASK)
        ),
        "unique_single_source": count_where(hist, lambda m: m.bit_count() == 1),
        "multi_source": count_where(hist, lambda m: m.bit_count() >= 2),
    }


def scope_row(name: str, data: pl.LazyFrame, mask_filter, hist: Counter) -> dict:
    lf = data.filter(mask_filter) if mask_filter is not None else data
    stats = lf.select([
        pl.len().alias("assertions"),
        pl.col("_weight").sum().alias("events"),
        pl.col("_drug").n_unique().alias("drug_name_keys"),
        pl.col("_disease").n_unique().alias("disease_name_keys"),
        pl.col("source").n_unique().alias("sources"),
        pl.col("relationship_type").n_unique().alias("relationship_types"),
        pl.col("evidence_type").n_unique().alias("evidence_types"),
    ]).collect(engine="streaming").row(0, named=True)
    totals = pair_totals(hist)
    pairs = {
        "Integrated": totals["integrated"],
        "PrimeKG only": totals["primekg"],
        "Hetionet only": totals["hetionet"],
        "PrimeKG + Hetionet union": totals["comparator_union"],
        "Other 9 sources union": totals["noncomparator_union"],
        "Integrated without CTD": totals["non_ctd_union"],
    }[name]
    return {"resource": name, **{k: int(v or 0) for k, v in stats.items()},
            "unique_name_pairs": pairs,
            "pct_integrated_pairs": round(100 * pairs / totals["integrated"], 3)}


def sem_table(data: pl.LazyFrame, column: str, selected: str) -> pl.DataFrame:
    group_defs = [
        ("Integrated", None),
        ("PrimeKG only", pl.col("source") == "PrimeKG"),
        ("Hetionet only", pl.col("source") == "Hetionet"),
        ("PrimeKG + Hetionet", pl.col("source").is_in(COMPARE)),
        ("Other 9 sources", ~pl.col("source").is_in(COMPARE)),
        ("Without CTD", pl.col("source") != "CTD"),
    ]
    out = []
    for name, criterion in group_defs:
        q = data if criterion is None else data.filter(criterion)
        group = q.group_by(column).agg([
            pl.len().alias("assertions"),
            pl.col("_weight").sum().alias("events"),
        ]).collect(engine="streaming")
        for row in group.iter_rows(named=True):
            out.append({"resource": name, selected: row[column],
                        "assertions": int(row["assertions"]),
                        "events": int(row["events"])})
    return pl.DataFrame(out).sort(["resource", "events"], descending=[False, True])


def indication_sensitivity(event_file: Path):
    """Only source-defined indication-oriented families from event ledger.

    This is a retrieval sensitivity analysis, NOT a regulator-verified treatment
    list; families are broad classifications used in Analysis 18.
    """
    print("[STEP 5] Event-level indication-oriented sensitivity", flush=True)
    cols = ["source", "drug_name", "disease_or_condition_name", "merge_evidence_family"]
    check_file(event_file, cols)
    selected = (
        scan(event_file, cols)
        .filter(pl.col("merge_evidence_family").is_in(FAMILIES))
        .with_columns([
            key("drug_name").alias("_drug"),
            key("disease_or_condition_name").alias("_disease"),
            source_bit(),
        ])
        .select(["source", "_drug", "_disease", "_bit", "merge_evidence_family"])
        .collect(engine="streaming")
    )
    if selected.is_empty():
        raise RuntimeError("No indication-oriented family events were found; inspect merge classifications")
    invalid = selected.filter(
        (pl.col("_drug") == "") | (pl.col("_disease") == "") | (pl.col("_bit") == 0)
    )
    eq("indication-oriented events with invalid keys", invalid.height, 0)
    grouped = selected.group_by(["_drug", "_disease"]).agg(
        pl.col("_bit").unique().sum().alias("mask")
    )
    hist = Counter(dict(
        (int(mask), int(n))
        for mask, n in grouped.group_by("mask").len().iter_rows()
    ))
    results = pair_totals(hist)
    eq("indication-family name-pair conservation", sum(hist.values()), grouped.height)
    source_families = (
        selected.group_by(["source", "merge_evidence_family"])
        .agg(pl.len().alias("events"))
        .sort(["source", "events"], descending=[False, True])
    )
    return selected.height, results, source_families


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--processed-dir", type=Path,
                        default=Path(__file__).resolve().parent / "data" / "processed")
    parser.add_argument("--buckets", type=int, default=128)
    parser.add_argument("--top", type=int, default=20)
    parser.add_argument("--save-tables", action="store_true")
    parser.add_argument("--keep-work", action="store_true")
    args = parser.parse_args()
    if not (8 <= args.buckets <= 1024 and 1 <= args.top <= 100):
        parser.error("--buckets 8..1024 and --top 1..100")

    base = args.processed_dir.expanduser().resolve()
    assertions = base / ASSERT_FILE
    events = base / EVENT_FILE
    manifest_file = base / MANIFEST_FILE
    required = ["source", "drug_name", "disease_or_condition_name",
                "relationship_type", "evidence_type", "supporting_event_count"]
    check_file(assertions, required)
    if not manifest_file.is_file():
        raise RuntimeError(f"Missing locked merge manifest: {manifest_file}")
    manifest = json.loads(manifest_file.read_text(encoding="utf-8"))
    eq("merge manifest status", manifest.get("status"), "success")
    eq("locked source assertions", int(manifest["assertions"]), EXPECTED["assertions"])
    eq("locked source events", int(manifest["source_events"]), EXPECTED["events"])
    eq("manifest sources", {x["source"] for x in manifest["sources"]}, set(SOURCES))

    work = base / ".analysis19_primekg_hetionet_work"
    if work.exists():
        raise RuntimeError(f"Previous scratch exists; inspect/rename: {work}")
    work.mkdir()
    complete = False
    try:
        print("\nANALYSIS 19 | INTEGRATED DATA VS PRIMEKG / HETIONET", flush=True)
        print("[NOTE] Drug/disease name keys: lowercase + trim + collapse whitespace")
        print("[NOTE] Comparison is to included PrimeKG/Hetionet source subsets, not all edges in their full graphs")
        print("[NOTE] More pairs does not imply more independent clinical evidence\n")
        print("[STEP 1] Stage projected assertion columns with conservative name keys", flush=True)
        parquet = work / "assertion_keys.parquet"
        (
            scan(assertions, required)
            .with_columns([
                key("drug_name").alias("_drug"),
                key("disease_or_condition_name").alias("_disease"),
                pl.col("supporting_event_count").cast(pl.Int64, strict=True).alias("_weight"),
                source_bit(),
            ])
            .select(["source", "_drug", "_disease", "_bit",
                     "relationship_type", "evidence_type", "_weight"])
            .sink_parquet(parquet, compression="zstd", engine="streaming")
        )
        data = pl.scan_parquet(parquet)
        stats = data.select([
            pl.len().alias("assertions"),
            pl.col("_weight").sum().alias("events"),
            pl.col("_drug").n_unique().alias("drugs"),
            pl.col("_disease").n_unique().alias("diseases"),
            (pl.col("_bit") == 0).sum().alias("unknown_sources"),
            ((pl.col("_drug") == "") | (pl.col("_disease") == "")).sum().alias("blank_keys"),
        ]).collect(engine="streaming").row(0, named=True)
        for name, field in [("assertions", "assertions"), ("events", "events"),
                            ("drug names", "drugs"), ("disease names", "diseases")]:
            expected = EXPECTED["drug_names"] if name == "drug names" else (
                EXPECTED["disease_names"] if name == "disease names" else EXPECTED[field])
            eq(name, int(stats[field]), expected)
        eq("invalid source labels", int(stats["unknown_sources"]), 0)
        eq("empty name keys", int(stats["blank_keys"]), 0)

        by_source = (data.group_by("source").agg([
            pl.len().alias("assertions"),
            pl.col("_weight").sum().alias("events"),
            pl.col("_drug").n_unique().alias("unique_drug_names"),
            pl.col("_disease").n_unique().alias("unique_disease_names"),
            pl.col("relationship_type").n_unique().alias("relationship_types"),
            pl.col("evidence_type").n_unique().alias("evidence_types"),
        ]).collect(engine="streaming"))
        eq("source-level events vs manifest",
           {r["source"]: r["events"] for r in by_source.to_dicts()},
           {r["source"]: int(r["observed"]) for r in manifest["sources"]})

        print("[STEP 2] Partition drug-disease name pairs into bounded-memory buckets", flush=True)
        bucket_root = work / "buckets"
        data.select(["_drug", "_disease", "_bit", "_weight"]).with_columns(
            (pl.struct(["_drug", "_disease"]).hash(seed=202604) % args.buckets)
            .cast(pl.Int32).alias("_bucket")
        ).sink_parquet(
            pl.PartitionBy(bucket_root, key="_bucket", include_key=True),
            compression="zstd", engine="streaming", mkdir=True,
        )
        parts = defaultdict(list)
        for part in bucket_root.rglob("*.parquet"):
            bucket = int(pl.read_parquet(part, columns=["_bucket"], n_rows=1)["_bucket"][0])
            parts[bucket].append(part)
        if not parts:
            raise RuntimeError("No pair buckets written")

        print("[STEP 3] Count comparator pair membership; examples and source contributions", flush=True)
        hist = Counter()
        seen_assertions = 0
        seen_events = 0
        top_new = []
        top_new_nonctd = []
        eligible_masks = [m for m in range(1, 1 << len(SOURCES))
                          if (m & COMPARATOR_MASK) == 0]
        nonctd_masks = [m for m in eligible_masks if m & NON_CTD_MASK]
        for _, paths in tqdm(sorted(parts.items()), desc="Analysis 19 partitions", unit="bucket"):
            block = pl.read_parquet(paths)
            seen_assertions += block.height
            seen_events += int(block["_weight"].sum())
            pairs = block.group_by(["_drug", "_disease"]).agg([
                pl.col("_bit").unique().sum().alias("mask"),
                pl.len().alias("assertions"),
                pl.col("_weight").sum().alias("events"),
            ])
            hist.update({int(mask): int(n) for mask, n in
                         pairs.group_by("mask").len().iter_rows()})
            for allowed, destination in [(eligible_masks, top_new),
                                         (nonctd_masks, top_new_nonctd)]:
                records = (pairs.filter(pl.col("mask").is_in(allowed))
                           .sort(["assertions", "events"], descending=True)
                           .head(args.top).to_dicts())
                destination.extend(records)
        eq("partitioned assertion conservation", seen_assertions, EXPECTED["assertions"])
        eq("partitioned evidence-event conservation", seen_events, EXPECTED["events"])
        totals = pair_totals(hist)
        eq("integrated name pairs", totals["integrated"], EXPECTED["pairs"])
        eq("multi-source name pairs", totals["multi_source"], EXPECTED["multi_source_pairs"])
        eq("PrimeKG name pairs", totals["primekg"], EXPECTED["primekg_pairs"])
        eq("Hetionet name pairs", totals["hetionet"], EXPECTED["hetionet_pairs"])
        eq("PrimeKG-Hetionet shared name pairs", totals["both"], EXPECTED["prime_het_overlap"])
        eq("comparator union plus outside = integrated", totals["comparator_union"] + totals["outside_both"], totals["integrated"])
        eq("comparator-union inclusion-exclusion", totals["comparator_union"],
           totals["primekg"] + totals["hetionet"] - totals["both"])

        print("[STEP 4] Compile resource scope, gain, Venn regions and source contributions", flush=True)
        selectors = [
            ("Integrated", None),
            ("PrimeKG only", pl.col("source") == "PrimeKG"),
            ("Hetionet only", pl.col("source") == "Hetionet"),
            ("PrimeKG + Hetionet union", pl.col("source").is_in(COMPARE)),
            ("Other 9 sources union", ~pl.col("source").is_in(COMPARE)),
            ("Integrated without CTD", pl.col("source") != "CTD"),
        ]
        scope = pl.DataFrame([scope_row(label, data, filt, hist) for label, filt in selectors])
        size = totals["integrated"]
        pk, ht = BITS["PrimeKG"], BITS["Hetionet"]
        venn = pl.DataFrame([
            {"region": "Not in PrimeKG or Hetionet", "name_pairs": totals["outside_both"]},
            {"region": "PrimeKG, not Hetionet", "name_pairs": totals["primekg"] - totals["both"]},
            {"region": "Hetionet, not PrimeKG", "name_pairs": totals["hetionet"] - totals["both"]},
            {"region": "In both PrimeKG and Hetionet", "name_pairs": totals["both"]},
        ]).with_columns((100 * pl.col("name_pairs") / size).round(4).alias("pct_integrated_pairs"))
        eq("mutually exclusive comparator regions", int(venn["name_pairs"].sum()), size)
        gain = pl.DataFrame([
            {"comparison": "Integrated pairs", "name_pairs": size},
            {"comparison": "PrimeKG pairs", "name_pairs": totals["primekg"]},
            {"comparison": "Hetionet pairs", "name_pairs": totals["hetionet"]},
            {"comparison": "PrimeKG + Hetionet union", "name_pairs": totals["comparator_union"]},
            {"comparison": "Additional beyond PrimeKG", "name_pairs": size - totals["primekg"]},
            {"comparison": "Additional beyond Hetionet", "name_pairs": size - totals["hetionet"]},
            {"comparison": "Additional beyond both", "name_pairs": totals["outside_both"]},
            {"comparison": "Additional beyond both (also in non-CTD source)",
             "name_pairs": totals["outside_both_non_ctd"]},
            {"comparison": "Integrated without CTD", "name_pairs": totals["non_ctd_union"]},
        ]).with_columns((100 * pl.col("name_pairs") / size).round(4).alias("pct_integrated_pairs"))

        source_pair_rows = []
        for source in SOURCES:
            bit = BITS[source]
            in_source = count_where(hist, lambda m, bit=bit: bool(m & bit))
            outside_comp = count_where(
                hist, lambda m, bit=bit: bool(m & bit) and not (m & COMPARATOR_MASK)
            )
            only_this = int(hist.get(bit, 0))
            source_pair_rows.append({
                "source": source, "unique_name_pairs": in_source,
                "pairs_not_in_either_comparator": outside_comp,
                "pairs_unique_to_this_source": only_this,
                "pairs_shared_with_either_comparator": in_source - outside_comp,
                "pct_source_pairs_not_in_comparators": round(100 * outside_comp / in_source, 3),
            })
            eq(f"Analysis 08 source pairs: {source}", in_source, SOURCE_PAIRS_08[source])
        per_source_pairs = pl.DataFrame(source_pair_rows)
        contribution = (by_source.join(per_source_pairs, on="source", how="left")
                        .sort("unique_name_pairs", descending=True))
        new_only = contribution.filter(~pl.col("source").is_in(COMPARE))
        relationship_scope = sem_table(data, "relationship_type", "relationship_type")
        evidence_scope = sem_table(data, "evidence_type", "evidence_type")

        def top_frame(records):
            records = sorted(records, key=lambda x: (-x["assertions"], -x["events"],
                                                      x["_drug"], x["_disease"]))[:args.top]
            return pl.DataFrame([
                {"drug_name_key": r["_drug"], "disease_name_key": r["_disease"],
                 "assertions": r["assertions"], "events": r["events"],
                 "n_sources": int(r["mask"]).bit_count(),
                 "sources": " | ".join(s for s in SOURCES if r["mask"] & BITS[s])}
                for r in records
            ]) if records else pl.DataFrame({
                "drug_name_key": [], "disease_name_key": [], "assertions": [],
                "events": [], "n_sources": [], "sources": [],
            })
        top_new_frame = top_frame(top_new)
        top_new_nonctd_frame = top_frame(top_new_nonctd)

        # Unlike raw coverage, this sensitivity compares only selected
        # indication-oriented source-reported evidence categories.
        ind_events, ind, ind_source_families = indication_sensitivity(events)
        indication_table = pl.DataFrame([
            {"comparison": "Indication-oriented source events", "value": ind_events},
            {"comparison": "Indication-oriented unique name pairs", "value": ind["integrated"]},
            {"comparison": "Indication-oriented pairs in PrimeKG", "value": ind["primekg"]},
            {"comparison": "Indication-oriented pairs in Hetionet", "value": ind["hetionet"]},
            {"comparison": "Indication-oriented PrimeKG + Hetionet union", "value": ind["comparator_union"]},
            {"comparison": "Indication-oriented pairs outside both", "value": ind["outside_both"]},
            {"comparison": "Indication-oriented pairs outside both, non-CTD supported",
             "value": ind["outside_both_non_ctd"]},
        ])

        tables = {
            "resource_scope": scope,
            "incremental_coverage": gain,
            "mutually_exclusive_comparator_regions": venn,
            "source_contribution": contribution,
            "non_comparator_source_contribution": new_only,
            "original_relationship_scope": relationship_scope,
            "original_evidence_scope": evidence_scope,
            "top_additional_pair_examples": top_new_frame,
            "top_additional_nonctd_pair_examples": top_new_nonctd_frame,
            "indication_oriented_sensitivity": indication_table,
            "indication_family_by_source": ind_source_families,
        }
        for i, (label, frame) in enumerate(tables.items(), 1):
            show(f"{i}. {label.upper().replace('_', ' ')}", frame)

        print("\n[MANUSCRIPT-READY SUMMARY]")
        print(
            f"The integrated 11-source resource contained {size:,} distinct "
            f"drug-disease name pairs. PrimeKG contained {totals['primekg']:,} pairs "
            f"and Hetionet contained {totals['hetionet']:,}, with "
            f"{totals['both']:,} shared by the two sources. Their combined coverage "
            f"was {totals['comparator_union']:,} pairs; "
            f"{totals['outside_both']:,} ({100*totals['outside_both']/size:.3f}%) "
            "integrated pairs were absent from both source-specific subsets. "
            f"Of these, {totals['outside_both_non_ctd']:,} were observed in at "
            "least one non-CTD source. These are name-based coverage comparisons, "
            "not evidence of regulatory approval, clinical benefit or novel discovery."
        )
        print(
            "[CAUTION] The integrated dataset necessarily contains both comparator "
            "subsets, so their 'outside integrated' count is zero by construction."
        )
        print(
            "[CAUTION] CTD inferred associations inflate broad coverage. "
            "Indication-oriented figures use source-reported classifications, "
            "not externally validated approval status."
        )
        print(
            "[CAUTION] Full PrimeKG and Hetionet include broader graph-edge types "
            "that are not evaluated by their extracted drug-disease subsets here."
        )
        if args.save_tables:
            outdir = base.parent.parent / "tables"
            outdir.mkdir(parents=True, exist_ok=True)
            for name, frame in tables.items():
                output = outdir / f"Analysis_19_{name}.csv"
                frame.write_csv(output)
                print("[SAVED]", output)
            report = {
                "snapshot": "paper11-locked-v1",
                "polars_version": pl.__version__,
                "expected_assertions": EXPECTED["assertions"],
                "expected_events": EXPECTED["events"],
                "all_pairs": totals,
                "indication_families": FAMILIES,
                "indication_event_count": ind_events,
                "indication_pairs": ind,
            }
            outfile = outdir / "Analysis_19_summary.json"
            outfile.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
            print("[SAVED]", outfile)
        complete = True
        print("\n[PASS] ANALYSIS 19 COMPLETE; LOCKED INPUTS UNCHANGED")
    finally:
        if complete or not args.keep_work:
            shutil.rmtree(work, ignore_errors=True)
        else:
            print("[SCRATCH RETAINED]", work)


if __name__ == "__main__":
    main()
