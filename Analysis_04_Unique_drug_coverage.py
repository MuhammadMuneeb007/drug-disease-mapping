
#!/usr/bin/env python3
"""Analysis 04: name-search coverage of the locked 11-source drug--disease database.

Exact name-key coverage follows Analysis 01: case-folding and whitespace collapsing
only. This is search-label overlap, not synonym/ontology resolution or efficacy.

Inputs: locked assertion CSV, pair CSV, merge manifest. Nothing is modified.
Default: print all summary tables. Optional: --save-tables.
Requires: polars, tqdm. Keeps RAM bounded by pair-hash partitioning.
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

SOURCES = {
    "AACT", "Broad Drug Repurposing Hub", "ChEMBL", "CTD", "DrugCentral",
    "Hetionet", "MEDI", "OpenTargets", "PrimeKG", "repoDB", "SIDER",
}
FILES = {
    "assertions": "ALL_INCLUDED_SOURCES_drug_disease_deduplicated.csv",
    "pairs": "ALL_INCLUDED_SOURCES_drug_disease_source_collapsed.csv",
    "manifest": "ALL_INCLUDED_SOURCES_merge_summary.json",
}
EXPECTED_ASSERTIONS = 4_534_385
EXPECTED_PAIRS = 4_470_601
ID_COLS = ["drug_identifier", "disease_or_condition_identifier"]


def check_file(path, required):
    if not path.is_file() or path.stat().st_size == 0:
        raise RuntimeError(f"Missing/empty input: {path}")
    with path.open("r", encoding="utf-8-sig", newline="") as f:
        header = next(csv.reader(f, strict=True), [])
    if len(header) != len(set(header)) or not set(required).issubset(header):
        raise RuntimeError(f"Invalid CSV schema: {path.name}; expected {required}")


def scan(path, cols):
    return pl.scan_csv(
        str(path), infer_schema=False, ignore_errors=False,
        truncate_ragged_lines=False, encoding="utf8",
    ).select(cols)


def name_key(col):
    return (
        pl.col(col).fill_null("").str.strip_chars()
        .str.replace_all(r"\s+", " ").str.to_lowercase()
    )


def show(title, frame):
    print("\n" + "=" * 110)
    print(f"{title}  [rows={frame.height:,}]")
    print("=" * 110)
    print("\t".join(frame.columns))
    for row in frame.iter_rows():
        print("\t".join(
            "" if value is None else str(value).replace("\n", " ").replace("\t", " ")
            for value in row
        ))


def must_equal(name, actual, expected):
    if actual != expected:
        raise RuntimeError(f"[FAIL] {name}: observed {actual!r}; expected {expected!r}")
    print(f"[PASS] {name}: {actual:,}" if isinstance(actual, int) else f"[PASS] {name}")


def source_support(frame, key_col, label, top):
    """Source-key rows are distinct: one source contributes at most once per key."""
    per_source = frame.group_by("source").len().rename({"len": f"unique_{label}_names"})
    global_support = frame.group_by(key_col).agg([
        pl.col("source").n_unique().alias("n_sources"),
        pl.col("source").unique().sort().str.join(" | ").alias("sources_present"),
    ])
    distribution = (
        global_support.group_by("n_sources")
        .len().rename({"len": f"unique_{label}_names"})
        .sort("n_sources")
    )
    leaders = (
        global_support.filter(pl.col("n_sources") > 1)
        .sort(["n_sources", key_col], descending=[True, False])
        .head(top)
    )
    return per_source, global_support.height, distribution, leaders


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--processed-dir", type=Path,
                   default=Path(__file__).resolve().parent / "data" / "processed")
    p.add_argument("--buckets", type=int, default=128)
    p.add_argument("--top", type=int, default=25)
    p.add_argument("--save-tables", action="store_true")
    p.add_argument("--drug", default="", help="Optional literal drug-name search")
    p.add_argument("--disease", default="", help="Optional literal disease-name search")
    p.add_argument("--examples", type=int, default=15)
    p.add_argument("--keep-work", action="store_true")
    args = p.parse_args()
    if not (8 <= args.buckets <= 1024 and 1 <= args.top <= 100 and 1 <= args.examples <= 100):
        p.error("--buckets 8..1024, --top 1..100, and --examples 1..100 required")

    base = args.processed_dir.expanduser().resolve()
    files = {k: base / v for k, v in FILES.items()}
    if not files["manifest"].is_file():
        raise RuntimeError(f"Missing locked manifest: {files['manifest']}")
    manifest = json.loads(files["manifest"].read_text(encoding="utf-8"))
    must_equal("merge manifest status", manifest.get("status"), "success")
    check_file(files["assertions"], [
        "source", "drug_name", "disease_or_condition_name", *ID_COLS,
        "relationship_type", "evidence_type", "supporting_event_count",
    ])
    check_file(files["pairs"], ["supporting_source_count", "merge_pair_id"])

    work = base / ".analysis_04_name_coverage_work"
    if work.exists():
        raise RuntimeError(f"Scratch already exists; inspect manually: {work}")
    work.mkdir(parents=True)
    completed = False
    try:
        print("[STEP 1] Streaming narrow, normalised-name assertion table", flush=True)
        narrow = work / "name_keys.parquet"
        lf = scan(files["assertions"], [
            "source", "drug_name", "disease_or_condition_name", *ID_COLS,
            "relationship_type", "evidence_type", "supporting_event_count",
        ]).with_columns([
            name_key("drug_name").alias("_drug"),
            name_key("disease_or_condition_name").alias("_disease"),
        ]).select([
            "source", "drug_name", "disease_or_condition_name",
            "_drug", "_disease", *ID_COLS,
            "relationship_type", "evidence_type", "supporting_event_count",
        ])
        lf.sink_parquet(narrow, engine="streaming", compression="zstd")
        data = pl.scan_parquet(narrow)

        print("[STEP 2] Validate totals, identifier coverage and names", flush=True)
        per_source = (
            data.group_by("source").agg([
                pl.len().alias("assertions"),
                pl.col("supporting_event_count")
                .cast(pl.Int64, strict=True).sum().alias("linked_events"),
                (pl.col("drug_identifier").fill_null("").str.strip_chars() != "")
                .sum().alias("assertions_with_drug_id"),
                (pl.col("disease_or_condition_identifier").fill_null("").str.strip_chars() != "")
                .sum().alias("assertions_with_disease_id"),
                ((pl.col("drug_identifier").fill_null("").str.strip_chars() != "")
                 & (pl.col("disease_or_condition_identifier").fill_null("").str.strip_chars() != ""))
                .sum().alias("assertions_with_both_ids"),
                ((pl.col("_drug") == "") | (pl.col("_disease") == ""))
                .sum().alias("blank_name_key_assertions"),
            ])
            .collect(engine="streaming")
        )
        assertion_total = int(per_source["assertions"].sum())
        must_equal("assertion rows", assertion_total, EXPECTED_ASSERTIONS)
        must_equal("manifest assertions", assertion_total, int(manifest["assertions"]))
        must_equal("assertion support -> events", int(per_source["linked_events"].sum()),
                   int(manifest["source_events"]))
        must_equal("11 source labels", set(per_source["source"].to_list()), SOURCES)
        must_equal("blank drug/disease name keys", int(per_source["blank_name_key_assertions"].sum()), 0)
        source_manifest = {row["source"]: int(row["observed"])
                           for row in manifest.get("sources", [])}
        must_equal("source event counts", {
            r["source"]: r["linked_events"] for r in per_source.to_dicts()
        }, source_manifest)

        source_drugs = data.select("source", "_drug").unique().collect(engine="streaming")
        source_diseases = data.select("source", "_disease").unique().collect(engine="streaming")
        drug_per, n_drugs, drug_dist, drug_top = source_support(
            source_drugs, "_drug", "drug", args.top)
        disease_per, n_diseases, disease_dist, disease_top = source_support(
            source_diseases, "_disease", "disease", args.top)

        print("[STEP 3] Partition name pairs for bounded-memory overlap analyses", flush=True)
        bucket_root = work / "name_pair_buckets"
        data.select([
            "source", "_drug", "_disease", "relationship_type", "evidence_type",
        ]).with_columns(
            (pl.struct(["_drug", "_disease"]).hash(seed=202604) % args.buckets)
            .cast(pl.Int32).alias("_bucket")
        ).sink_parquet(
            pl.PartitionBy(bucket_root, key="_bucket", include_key=True),
            engine="streaming", compression="zstd", mkdir=True,
        )
        parts = defaultdict(list)
        for path in bucket_root.rglob("*.parquet"):
            b = int(pl.read_parquet(path, columns=["_bucket"], n_rows=1)["_bucket"][0])
            parts[b].append(path)
        if not parts:
            raise RuntimeError("No pair partitions written")

        total_pairs = 0
        source_pairs = Counter()
        pair_distribution = Counter()
        diversity = Counter()
        multi_pair_total = 0
        top_multi = []
        bucket_assertions = 0
        print("[STEP 4] Count pairs, overlapping sources, relationship diversity", flush=True)
        for bucket, paths in tqdm(sorted(parts.items()), desc="Name-pair buckets", unit="bucket"):
            block = pl.read_parquet(paths).drop("_bucket")
            bucket_assertions += block.height
            source_pair_frame = block.select("source", "_drug", "_disease").unique()
            for rec in source_pair_frame.group_by("source").len().iter_rows(named=True):
                source_pairs[rec["source"]] += int(rec["len"])

            pairs = block.group_by(["_drug", "_disease"]).agg([
                pl.len().alias("assertions"),
                pl.col("source").n_unique().alias("n_sources"),
                pl.col("relationship_type").n_unique().alias("n_relationship_types"),
                pl.col("evidence_type").n_unique().alias("n_evidence_types"),
            ])
            total_pairs += pairs.height
            for rec in pairs.group_by("n_sources").len().iter_rows(named=True):
                pair_distribution[int(rec["n_sources"])] += int(rec["len"])
            for rec in pairs.group_by([
                "n_sources", "n_relationship_types", "n_evidence_types"
            ]).len().iter_rows(named=True):
                diversity[(int(rec["n_sources"]), int(rec["n_relationship_types"]),
                           int(rec["n_evidence_types"]))] += int(rec["len"])

            multis = pairs.filter(pl.col("n_sources") > 1)
            multi_pair_total += multis.height
            if multis.height:
                src = source_pair_frame.group_by(["_drug", "_disease"]).agg(
                    pl.col("source").sort().str.join(" | ").alias("sources_present")
                )
                leaders = (multis.sort(["n_sources", "assertions"], descending=True)
                           .head(args.top).join(src, on=["_drug", "_disease"], how="left"))
                top_multi.extend(leaders.to_dicts())

        must_equal("pair-bucket assertion conservation", bucket_assertions, assertion_total)
        must_equal("unique name-pair support distribution", sum(pair_distribution.values()), total_pairs)
        must_equal("name-pair semantic diversity conservation", sum(diversity.values()), total_pairs)
        must_equal("multi-source name-pair conservation", sum(
            count for n, count in pair_distribution.items() if n > 1), multi_pair_total)
        top_multi.sort(key=lambda r: (-r["n_sources"], -r["assertions"], r["_drug"], r["_disease"]))
        top_multi = top_multi[:args.top]

        summary = (per_source.join(drug_per, on="source", how="left")
                   .join(disease_per, on="source", how="left")
                   .join(pl.DataFrame({"source": list(source_pairs),
                                       "unique_name_pairs": list(source_pairs.values())}),
                         on="source", how="left")
                   .sort("assertions", descending=True))
        must_equal("all source name-pair coverage measured", summary.height, 11)

        pair_counts = (scan(files["pairs"], ["merge_pair_id", "supporting_source_count"])
                       .with_columns(pl.col("supporting_source_count").cast(pl.Int64))
                       .group_by("supporting_source_count")
                       .agg(pl.len().alias("exact_identity_pairs"))
                       .sort("supporting_source_count")
                       .collect(engine="streaming"))
        must_equal("locked exact identity pairs", int(pair_counts["exact_identity_pairs"].sum()),
                   EXPECTED_PAIRS)
        must_equal("manifest exact identity pairs", int(manifest["exact_identity_pairs"]),
                   EXPECTED_PAIRS)

        def distribution(counter, keyname, countname):
            return pl.DataFrame({keyname: sorted(counter),
                                 countname: [counter[n] for n in sorted(counter)]})

        pair_dist = distribution(pair_distribution, "n_sources", "unique_name_pairs")
        semantics = pl.DataFrame([
            {"n_sources": k[0], "n_relationship_types": k[1],
             "n_evidence_types": k[2], "unique_name_pairs": n}
            for k, n in sorted(diversity.items())
        ])
        top_pair_df = pl.DataFrame(top_multi) if top_multi else pl.DataFrame({
            "_drug": [], "_disease": [], "n_sources": [], "assertions": [],
            "n_relationship_types": [], "n_evidence_types": [], "sources_present": [],
        })

        show("1. OVERALL NAME-BASED COVERAGE", pl.DataFrame({
            "metric": ["source_assertions", "original_events_linked", "sources",
                       "unique_drug_name_keys", "unique_disease_name_keys",
                       "unique_name_pairs", "name_pairs_supported_by_multiple_sources",
                       "locked_exact_identity_pairs"],
            "value": [assertion_total, int(per_source["linked_events"].sum()), 11,
                      n_drugs, n_diseases, total_pairs, multi_pair_total,
                      EXPECTED_PAIRS],
        }))
        show("2. SOURCE COVERAGE AND ASSERTION-LEVEL ID COMPLETENESS", summary)
        show("3. DRUG NAME: NUMBER OF SOURCES", drug_dist)
        show("4. DISEASE NAME: NUMBER OF SOURCES", disease_dist)
        show("5. NAME PAIR: NUMBER OF SOURCES", pair_dist)
        show("6. TOP MULTI-SOURCE DRUG NAMES", drug_top)
        show("7. TOP MULTI-SOURCE DISEASE NAMES", disease_top)
        show("8. TOP MULTI-SOURCE DRUG-DISEASE NAME PAIRS", top_pair_df)
        show("9. PAIR RELATIONSHIP/EVIDENCE DIVERSITY", semantics)
        show("10. STRICT MERGE IDENTITY OVERLAP FOR COMPARISON", pair_counts)

        print("\n[INTERPRETATION]")
        print("Source assertions (4,534,385) are not original evidence events (10,895,395).")
        print("Name keys use lowercase and collapsed whitespace, matching Analysis 01.")
        print("Same name across sources = name-level retrieval overlap, NOT biological equivalence.")
        print("Overlap is not independent corroboration and does not imply approved use.")
        print("Summing per-source unique pair counts does not equal global unique pairs.")
        print(f"[PASS] Analysis 04: {n_drugs:,} drug keys; {n_diseases:,} disease keys; "
              f"{total_pairs:,} name pairs; {multi_pair_total:,} multi-source name pairs")

        if args.drug.strip() or args.disease.strip():
            print("\n[OPTIONAL SEARCH] Literal substring of normalised names (assertion rows)")
            query = data
            if args.drug.strip():
                query = query.filter(pl.col("_drug").str.contains(
                    args.drug.casefold().strip(), literal=True))
            if args.disease.strip():
                query = query.filter(pl.col("_disease").str.contains(
                    args.disease.casefold().strip(), literal=True))
            query_source = (query.group_by("source").agg([
                pl.len().alias("assertions"),
                pl.col("supporting_event_count").cast(pl.Int64).sum().alias("events"),
            ]).collect(engine="streaming").sort("events", descending=True))
            show("SEARCH MATCHES BY SOURCE", query_source)
            show("FIRST MATCHING ASSERTIONS", query.select([
                "drug_name", "disease_or_condition_name", "source",
                "relationship_type", "evidence_type", "supporting_event_count",
            ]).head(args.examples).collect(engine="streaming"))

        if args.save_tables:
            target = base.parent.parent / "tables"
            target.mkdir(parents=True, exist_ok=True)
            tables = {
                "source_coverage": summary,
                "drug_source_support": drug_dist,
                "disease_source_support": disease_dist,
                "name_pair_source_support": pair_dist,
                "top_multi_source_drugs": drug_top,
                "top_multi_source_diseases": disease_top,
                "top_multi_source_name_pairs": top_pair_df,
                "pair_semantic_diversity": semantics,
                "exact_identity_pair_support": pair_counts,
            }
            for key_name, frame in tables.items():
                out = target / f"Analysis_04_{key_name}.csv"
                frame.write_csv(out)
                print("[SAVED]", out)
        completed = True
        print("ANALYSIS 04 COMPLETE; LOCKED DATA UNCHANGED")
    finally:
        if completed or not args.keep_work:
            shutil.rmtree(work, ignore_errors=True)
        else:
            print("[SCRATCH RETAINED]", work)


if __name__ == "__main__":
    main()
