#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Analysis 18: Type 2 diabetes worked-example retrieval from locked paper-11 data.

Queries preserved drug/disease labels, original relationships/evidence labels,
analytical evidence families, source identifiers, and event support. Does not
modify locked inputs. A name match is not drug approval or independent evidence.

Run:
    POLARS_MAX_THREADS=16 python Analysis_18_Worked_example_type2_diabetes_indication_use_case.py --save-tables
Requires: polars, tqdm
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import polars as pl
from tqdm import tqdm

SOURCES = {
    "AACT", "Broad Drug Repurposing Hub", "ChEMBL", "CTD", "DrugCentral",
    "Hetionet", "MEDI", "OpenTargets", "PrimeKG", "repoDB", "SIDER",
}
EXPECTED_EVENTS = 10_895_395
EXPECTED_ASSERTIONS = 4_534_385

# Grouped evidence is source-reported, not independent therapeutic verification.
INDICATION_FAMILIES = [
    "source_reported_approval", "source_reported_indication",
    "ctd_direct_therapeutic", "label_nlp_indication",
]

# Disease-query aliases are purposefully constrained. Generic diabetes, Type 1,
# prediabetes, gestational diabetes and ambiguous adult-onset diabetes do not
# match by default. The query key is not used to merge ontology identities.
DISEASE_RE = (
    r"\b(?:type\s*(?:2|ii)\s+diabetes|diabetes\s+(?:mellitus\s+)?"
    r"type\s*(?:2|ii)|non\s+insulin\s+dependent\s+diabetes|t2dm|t2d|niddm)\b"
)

ASSERT_FIELDS = [
    "merge_assertion_id", "source", "drug_name", "drug_identifier",
    "drug_identifier_type", "disease_or_condition_name",
    "disease_or_condition_identifier", "disease_or_condition_identifier_type",
    "relationship_type", "evidence_type", "supporting_event_count",
]
EVENT_FIELDS = [
    "merge_assertion_id", "merge_event_id", "source", "drug_name",
    "drug_identifier", "drug_identifier_type", "disease_or_condition_name",
    "disease_or_condition_identifier", "relationship_type", "evidence_type",
    "merge_evidence_family",
]


def header_checked(path: Path, fields: list[str]) -> None:
    if not path.is_file() or not path.stat().st_size:
        raise RuntimeError(f"Missing/empty locked input: {path}")
    with path.open("r", encoding="utf-8-sig", newline="") as fh:
        header = next(csv.reader(fh, strict=True), [])
    if len(header) != len(set(header)) or not set(fields).issubset(header):
        raise RuntimeError(f"Schema mismatch: {path.name}; missing={sorted(set(fields)-set(header))}")


def scan(path: Path, fields: list[str]) -> pl.LazyFrame:
    return pl.scan_csv(
        str(path), infer_schema=False, ignore_errors=False,
        truncate_ragged_lines=False, encoding="utf8",
    ).select(fields)


def name_key(field: str) -> pl.Expr:
    # Exactly the same conservative name key as Analyses 01, 04, 08.
    return (pl.col(field).fill_null("").str.strip_chars()
            .str.replace_all(r"\s+", " ").str.to_lowercase())


def disease_query_key(field: str) -> pl.Expr:
    # Applies ONLY to the query; does not rewrite stored disease names.
    return (name_key(field).str.replace_all(r"[-_/(),;:]+", " ")
            .str.replace_all(r"\s+", " ").str.strip_chars())


def require(label: str, actual, expected) -> None:
    if actual != expected:
        raise RuntimeError(f"[FAIL] {label}: got {actual!r}, expected {expected!r}")
    print(f"[PASS] {label}: {actual:,}" if isinstance(actual, int)
          else f"[PASS] {label}", flush=True)


def print_tsv(title: str, frame: pl.DataFrame) -> None:
    print("\n" + "=" * 105)
    print(f"{title} [{frame.height:,} rows]")
    print("=" * 105)
    if not frame.height:
        print("[NO MATCHES]")
        return
    print("\t".join(frame.columns))
    for record in frame.iter_rows():
        print("\t".join(
            "" if value is None else str(value).replace("\t", " ")
            .replace("\r", " ").replace("\n", " ")
            for value in record
        ))


def add_keys(lf: pl.LazyFrame) -> pl.LazyFrame:
    return lf.with_columns([
        name_key("drug_name").alias("_drug_key"),
        name_key("disease_or_condition_name").alias("_disease_key"),
        disease_query_key("disease_or_condition_name").alias("_query_disease"),
    ]).filter(pl.col("_query_disease").str.contains(DISEASE_RE))


def stats(frame: pl.DataFrame, label: str) -> dict:
    if not frame.height:
        return {
            "group": label, "rows": 0, "drug_names": 0, "disease_labels": 0,
            "name_pairs": 0, "sources": 0,
        }
    pairs = frame.select(["_drug_key", "_disease_key"]).unique()
    return {
        "group": label,
        "rows": frame.height,
        "drug_names": frame["_drug_key"].n_unique(),
        "disease_labels": frame["_disease_key"].n_unique(),
        "name_pairs": pairs.height,
        "sources": frame["source"].n_unique(),
    }


def summary_table(frame: pl.DataFrame, group_cols: list[str], label: str) -> pl.DataFrame:
    return (
        frame.group_by(group_cols)
        .agg([
            pl.len().alias(label + "_events"),
            pl.col("_drug_key").n_unique().alias("unique_drug_names"),
            pl.struct(["_drug_key", "_disease_key"]).n_unique().alias("unique_name_pairs"),
        ])
        .sort(label + "_events", descending=True)
    )


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--processed-dir", type=Path,
                   default=Path(__file__).resolve().parent / "data" / "processed")
    p.add_argument("--top", type=int, default=50)
    p.add_argument("--examples", type=int, default=20)
    p.add_argument("--save-tables", action="store_true")
    args = p.parse_args()
    if not 1 <= args.top <= 200 or not 1 <= args.examples <= 100:
        p.error("--top must be 1..200 and --examples 1..100")

    base = args.processed_dir.expanduser().resolve()
    assertion_file = base / "ALL_INCLUDED_SOURCES_drug_disease_deduplicated.csv"
    event_file = base / "ALL_INCLUDED_SOURCES_drug_disease_merged.csv"
    manifest_file = base / "ALL_INCLUDED_SOURCES_merge_summary.json"
    header_checked(assertion_file, ASSERT_FIELDS)
    header_checked(event_file, EVENT_FIELDS)
    if not manifest_file.is_file():
        raise RuntimeError(f"Missing merge manifest: {manifest_file}")
    manifest = json.loads(manifest_file.read_text(encoding="utf-8"))
    require("merge status", manifest.get("status"), "success")
    require("locked full source events", int(manifest["source_events"]), EXPECTED_EVENTS)
    require("locked assertions", int(manifest["assertions"]), EXPECTED_ASSERTIONS)

    print("\nANALYSIS 18 | TYPE 2 DIABETES SEARCH AND INDICATION-EVIDENCE REVIEW")
    print("[QUERY] Explicit Type 2 diabetes synonyms/abbreviations only")
    print("[RULE] Search labels are not ontology harmonisation or clinical validation")

    with tqdm(total=3, desc="Analysis 18", unit="step") as bar:
        print("[STEP 1] Find matching source assertions", flush=True)
        assertions = (
            add_keys(scan(assertion_file, ASSERT_FIELDS))
            .with_columns(pl.col("supporting_event_count").cast(pl.Int64, strict=True))
            .collect(engine="streaming")
        )
        bar.update(1)

        print("[STEP 2] Retrieve matching original evidence events", flush=True)
        events = add_keys(scan(event_file, EVENT_FIELDS)).collect(engine="streaming")
        bar.update(1)

        if not assertions.height or not events.height:
            raise RuntimeError(
                "Type 2 diabetes query returned no matching assertions/events. "
                "Review matched disease terminology before interpreting zero as no evidence."
            )
        if not set(assertions["source"].to_list()).issubset(SOURCES):
            raise RuntimeError("Unexpected source in selected assertions")
        if not set(events["source"].to_list()).issubset(SOURCES):
            raise RuntimeError("Unexpected source in selected events")

        require("unique retrieved assertion IDs",
                assertions["merge_assertion_id"].n_unique(), assertions.height)
        require("matched evidence-event preservation",
                int(assertions["supporting_event_count"].sum()), events.height)

        by_aid = (events.group_by("merge_assertion_id")
                  .agg(pl.len().alias("recounted_events")))
        audit = assertions.select([
            "merge_assertion_id", "supporting_event_count",
        ]).join(by_aid, on="merge_assertion_id", how="full", coalesce=True)
        mismatch = audit.filter(
            pl.col("supporting_event_count").is_null()
            | pl.col("recounted_events").is_null()
            | (pl.col("supporting_event_count") != pl.col("recounted_events"))
        )
        require("missing/unreconciled assertion support rows", mismatch.height, 0)
        require("unique matching original event IDs",
                events["merge_event_id"].n_unique(), events.height)

        # Every queried evidence event stays in the broad search output;
        # the indication-oriented view is a distinct evidence-based filter.
        selected = events.filter(
            pl.col("merge_evidence_family").is_in(INDICATION_FAMILIES)
        )
        excluded = events.filter(
            ~pl.col("merge_evidence_family").is_in(INDICATION_FAMILIES)
        )
        require("indication/other partition conserves events",
                selected.height + excluded.height, events.height)

        # Tables preserve source-specific evidence semantics.
        source_all = summary_table(events, ["source"], "matched")
        relationship_all = summary_table(events, ["relationship_type"], "matched")
        evidence_all = summary_table(events, ["evidence_type"], "matched")
        families = summary_table(events, ["merge_evidence_family"], "matched")
        source_family = summary_table(
            events, ["source", "merge_evidence_family"], "matched")
        source_indication = summary_table(selected, ["source"], "indication_oriented")
        relationship_indication = summary_table(
            selected, ["relationship_type"], "indication_oriented")
        evidence_indication = summary_table(
            selected, ["evidence_type"], "indication_oriented")
        matched_diseases = (
            events.group_by([
                "disease_or_condition_name", "disease_or_condition_identifier",
            ]).agg([
                pl.len().alias("events"),
                pl.col("source").n_unique().alias("sources"),
                pl.col("_drug_key").n_unique().alias("drug_names"),
            ]).sort("events", descending=True)
        )

        def drugs_overview(df: pl.DataFrame) -> pl.DataFrame:
            if not df.height:
                return pl.DataFrame({
                    "drug_name_key": [], "source_count": [], "source_names": [],
                    "event_count": [], "assertion_count": [], "evidence_families": [],
                    "drug_names": [], "drug_identifiers": [],
                    "disease_labels": [], "name_pairs": [],
                })
            return (
                df.group_by("_drug_key").agg([
                    pl.len().alias("event_count"),
                    pl.col("merge_assertion_id").n_unique().alias("assertion_count"),
                    pl.col("source").n_unique().alias("source_count"),
                    pl.col("source").unique().sort().str.join(" | ").alias("source_names"),
                    pl.col("merge_evidence_family").unique().sort().str.join(" | ")
                    .alias("evidence_families"),
                    pl.col("drug_name").unique().sort().str.join(" | ").alias("drug_names"),
                    pl.concat_str(["drug_identifier_type", "drug_identifier"], separator=":")
                    .filter(pl.col("drug_identifier").str.strip_chars() != "")
                    .unique().sort().str.join(" | ").alias("drug_identifiers"),
                    pl.col("disease_or_condition_name").unique().sort().str.join(" | ")
                    .alias("disease_labels"),
                    pl.col("_disease_key").n_unique().alias("name_pairs"),
                ]).rename({"_drug_key": "drug_name_key"})
                .sort(["source_count", "event_count", "drug_name_key"],
                      descending=[True, True, False])
            )

        broad_drugs = drugs_overview(events)
        indication_drugs = drugs_overview(selected)
        matching_assertions_by_source = (
            assertions.group_by("source").agg([
                pl.len().alias("assertions"),
                pl.col("supporting_event_count").sum().alias("events"),
            ]).sort("events", descending=True)
        )
        bar.update(1)

    overview = pl.DataFrame([
        stats(assertions, "all_matched_assertions"),
        stats(events, "all_matched_evidence_events"),
        stats(selected, "indication_oriented_evidence_events"),
        stats(excluded, "other_evidence_events"),
    ])

    tables = {
        "overall": overview,
        "matched_assertions_by_source": matching_assertions_by_source,
        "all_evidence_by_source": source_all,
        "original_relationships": relationship_all,
        "original_evidence_types": evidence_all,
        "analytical_evidence_families": families,
        "source_evidence_families": source_family,
        "indication_oriented_by_source": source_indication,
        "indication_oriented_relationships": relationship_indication,
        "indication_oriented_evidence_types": evidence_indication,
        "matched_disease_labels": matched_diseases,
        "top_broad_search_drugs": broad_drugs.head(args.top),
        "top_indication_oriented_drugs": indication_drugs.head(args.top),
        "example_original_evidence": events.select([
            "drug_name", "drug_identifier", "disease_or_condition_name",
            "disease_or_condition_identifier", "source", "relationship_type",
            "evidence_type", "merge_evidence_family", "merge_event_id",
        ]).head(args.examples),
    }
    for i, (name, frame) in enumerate(tables.items(), 1):
        print_tsv(f"{i}. {name.replace('_', ' ').upper()}", frame)

    print("\n[MANUSCRIPT-READY WORKED EXAMPLE]")
    print(
        "A Type 2 diabetes name-query retrieved "
        f"{assertions.height:,} source-specific assertions linked to "
        f"{events.height:,} preserved evidence events from "
        f"{events['source'].n_unique()} sources. "
        f"These comprised {events['_drug_key'].n_unique():,} normalised "
        "drug-name keys. Evidence-family selection identified "
        f"{selected.height:,} indication-oriented source events representing "
        f"{selected['_drug_key'].n_unique() if selected.height else 0:,} "
        "drug-name keys; other retrieved events were retained separately. "
        "This example demonstrates source-aware retrieval rather than verified "
        "regulatory approval or treatment efficacy."
    )
    print("[NOTE] The indication-oriented list is an evidence retrieval result, not a treatment recommendation.")
    print("[NOTE] Generic diabetes/Type 1/gestational terms are not automatically pooled.")
    print("[NOTE] Original drug names may include source-supplied identifier-looking labels.")
    if not selected.height:
        print("[REVIEW] No indication-oriented evidence families found. Check family labels and disease terms.")
    else:
        print("[PASS] Type 2 diabetes indication-oriented retrieval and support reconciliation")

    if args.save_tables:
        folder = base.parent.parent / "tables"
        folder.mkdir(parents=True, exist_ok=True)
        for name, frame in tables.items():
            path = folder / f"Analysis_18_{name}.csv"
            frame.write_csv(path)
            print("[SAVED]", path)
        report = {
            "query": "Type 2 diabetes controlled synonym expression",
            "indication_families": INDICATION_FAMILIES,
            "manifest_events": EXPECTED_EVENTS,
            "manifest_assertions": EXPECTED_ASSERTIONS,
            "matched_assertions": assertions.height,
            "matched_events": events.height,
            "indication_oriented_events": selected.height,
            "other_evidence_events": excluded.height,
            "matched_source_count": events["source"].n_unique(),
            "polars_version": pl.__version__,
        }
        path = folder / "Analysis_18_summary.json"
        path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
        print("[SAVED]", path)
    print("ANALYSIS 18 COMPLETE | LOCKED INPUTS UNCHANGED")


if __name__ == "__main__":
    main()