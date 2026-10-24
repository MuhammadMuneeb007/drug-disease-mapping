#!/usr/bin/env python3
"""Analysis 02: evidence types in locked 11-source drug--disease dataset.

Reads full evidence events, deduplicated source assertions, and merger manifest.
Original evidence_type versus analytical merge_evidence_family are separate.
Nothing is deduplicated or changed. No outputs unless --save-tables.

Install: python -m pip install polars tqdm rich
Run: POLARS_MAX_THREADS=16 python Analysis_02_Evidence_type_distribution.py
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import polars as pl
from rich.console import Console
from rich.table import Table
from tqdm import tqdm

SOURCES = {
    "AACT", "Broad Drug Repurposing Hub", "ChEMBL", "CTD",
    "DrugCentral", "Hetionet", "MEDI", "OpenTargets",
    "PrimeKG", "repoDB", "SIDER",
}
EXPECTED_EVENTS = 10_895_395
EXPECTED_ASSERTIONS = 4_534_385
EVENT_NAME = "ALL_INCLUDED_SOURCES_drug_disease_merged.csv"
ASSERT_NAME = "ALL_INCLUDED_SOURCES_drug_disease_deduplicated.csv"
MANIFEST_NAME = "ALL_INCLUDED_SOURCES_merge_summary.json"
console = Console()


def fail(message):
    raise RuntimeError(message)


def check_file(path, required):
    if not path.is_file() or path.stat().st_size == 0:
        fail(f"Input file missing or empty: {path}")
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        header = next(csv.reader(handle, strict=True), [])
    if len(header) != len(set(header)) or not set(required).issubset(header):
        fail(f"Invalid header in {path.name}; required={required}; actual={header}")


def scan(path, fields):
    return pl.scan_csv(
        path,
        infer_schema=False,
        ignore_errors=False,
        truncate_ragged_lines=False,
        encoding="utf8",
    ).select(fields)


def display(title, frame):
    console.rule(f"[bold cyan]{title}[/bold cyan]")
    table = Table(show_lines=False, header_style="bold", expand=False)
    for field, dtype in zip(frame.columns, frame.dtypes):
        table.add_column(
            field.replace("_", " "),
            justify="right" if dtype.is_numeric() else "left",
            no_wrap=True,
        )
    for item in frame.iter_rows():
        rendered = []
        for value in item:
            if isinstance(value, int):
                rendered.append(f"{value:,}")
            elif isinstance(value, float):
                rendered.append(f"{value:.3f}")
            else:
                rendered.append("" if value is None else str(value))
        table.add_row(*rendered)
    console.print(table)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--processed-dir", type=Path,
        default=Path(__file__).resolve().parent / "data" / "processed",
    )
    parser.add_argument(
        "--save-tables", action="store_true",
        help="Save Analysis 02 summary tables under repository tables/",
    )
    args = parser.parse_args()
    base = args.processed_dir.expanduser().resolve()
    event_file = base / EVENT_NAME
    assertion_file = base / ASSERT_NAME
    manifest_file = base / MANIFEST_NAME
    check_file(event_file, ["source", "evidence_type", "merge_evidence_family"])
    check_file(assertion_file, ["source", "evidence_type", "supporting_event_count"])
    if not manifest_file.is_file():
        fail(f"Missing merger manifest: {manifest_file}")
    manifest = json.loads(manifest_file.read_text(encoding="utf-8"))
    if manifest.get("status") != "success":
        fail("Manifest does not report successful merge")

    console.rule("[bold]ANALYSIS 02 | EVIDENCE TYPE DISTRIBUTION[/bold]")
    console.print(f"Source events: [cyan]{event_file}[/cyan]")
    console.print(f"Assertion summary: [cyan]{assertion_file}[/cyan]")
    console.print(
        "Event and assertion counts are different denominators; "
        "neither is an approved-indication count.\n"
    )

    with tqdm(total=4, desc="Analysis 02", unit="stage") as progress:
        # Full events: source evidence labels AND merger's analytical classes.
        # Grouping only these three columns avoids collecting 200-wide rows.
        crosswalk = (
            scan(event_file, ["source", "evidence_type", "merge_evidence_family"])
            .group_by(["source", "evidence_type", "merge_evidence_family"])
            .agg(pl.len().alias("events"))
            .collect(engine="streaming")
        )
        progress.update(1)

        # One assertion row can represent multiple original evidence events.
        assertion_rows = (
            scan(assertion_file, ["source", "evidence_type", "supporting_event_count"])
            .with_columns(
                pl.col("supporting_event_count").cast(pl.Int64, strict=True)
            )
            .group_by(["source", "evidence_type"])
            .agg([
                pl.len().alias("assertions"),
                pl.col("supporting_event_count").sum().alias("linked_events"),
            ])
            .collect(engine="streaming")
        )
        progress.update(1)

        events_by_label = (
            crosswalk.group_by(["source", "evidence_type"])
            .agg(pl.col("events").sum().alias("events"))
        )
        by_source_type = (
            events_by_label.join(
                assertion_rows, on=["source", "evidence_type"],
                how="full", coalesce=True,
            )
            .with_columns([
                pl.col("events").fill_null(0),
                pl.col("assertions").fill_null(0),
                pl.col("linked_events").fill_null(0),
            ])
        )
        sources_seen = set(by_source_type["source"].to_list())
        if sources_seen != SOURCES:
            fail(
                f"Source mismatch: missing={SOURCES-sources_seen}, "
                f"unexpected={sources_seen-SOURCES}"
            )
        if by_source_type.filter(
            (pl.col("events") != pl.col("linked_events"))
            | (pl.col("assertions") > pl.col("events"))
            | (pl.col("assertions") == 0)
        ).height:
            fail("Evidence-type conservation mismatch between events and assertions")

        event_total = int(by_source_type["events"].sum())
        assertion_total = int(by_source_type["assertions"].sum())
        if event_total != EXPECTED_EVENTS or assertion_total != EXPECTED_ASSERTIONS:
            fail(
                "Unexpected locked snapshot totals: "
                f"events={event_total}, assertions={assertion_total}"
            )
        if (
            event_total != int(manifest["source_events"])
            or assertion_total != int(manifest["assertions"])
        ):
            fail("Event/assertion counts do not match merge manifest")

        manifest_source_counts = {
            x["source"]: int(x["observed"])
            for x in manifest.get("sources", [])
        }

        sources = (
            by_source_type.group_by("source")
            .agg([
                pl.col("events").sum().alias("events"),
                pl.col("assertions").sum().alias("assertions"),
                pl.col("evidence_type").n_unique().alias("evidence_types"),
            ])
            .with_columns([
                (pl.col("events") - pl.col("assertions"))
                .alias("additional_support_events"),
                (100 * pl.col("events") / event_total)
                .round(3).alias("pct_all_events"),
            ])
            .sort("events", descending=True)
        )
        if {
            r["source"]: r["events"] for r in sources.to_dicts()
        } != manifest_source_counts:
            fail("Source row counts do not reconcile with merge manifest")

        overall = (
            by_source_type.group_by("evidence_type")
            .agg([
                pl.col("events").sum().alias("events"),
                pl.col("assertions").sum().alias("assertions"),
                pl.col("source").n_unique().alias("sources"),
            ])
            .with_columns([
                (pl.col("events") - pl.col("assertions"))
                .alias("additional_support_events"),
                (100 * pl.col("events") / event_total)
                .round(3).alias("pct_all_events"),
                (100 * pl.col("assertions") / assertion_total)
                .round(3).alias("pct_all_assertions"),
            ])
            .sort("events", descending=True)
        )

        by_source_type = (
            by_source_type.join(
                sources.select([
                    "source",
                    pl.col("events").alias("source_total_events"),
                    pl.col("assertions").alias("source_total_assertions"),
                ]),
                on="source", how="left",
            )
            .with_columns([
                (pl.col("events") - pl.col("assertions"))
                .alias("additional_support_events"),
                (100 * pl.col("events") / pl.col("source_total_events"))
                .round(3).alias("pct_source_events"),
                (100 * pl.col("assertions") / pl.col("source_total_assertions"))
                .round(3).alias("pct_source_assertions"),
            ])
            .drop("source_total_events", "source_total_assertions", "linked_events")
            .sort(["source", "events"], descending=[False, True])
        )

        family = (
            crosswalk.group_by("merge_evidence_family")
            .agg(pl.col("events").sum().alias("events"))
            .with_columns(
                (100 * pl.col("events") / event_total)
                .round(3).alias("pct_all_events")
            )
            .sort("events", descending=True)
        )

        family_by_source = (
            crosswalk.group_by(["source", "merge_evidence_family"])
            .agg(pl.col("events").sum().alias("events"))
            .sort(["source", "events"], descending=[False, True])
        )

        families_without_ctd = (
            family_by_source.filter(pl.col("source") != "CTD")
            .group_by("merge_evidence_family")
            .agg(pl.col("events").sum().alias("events"))
        )
        without_ctd_n = int(families_without_ctd["events"].sum())
        families_without_ctd = (
            families_without_ctd.with_columns(
                (100 * pl.col("events") / without_ctd_n)
                .round(3).alias("pct_non_ctd_events")
            )
            .sort("events", descending=True)
        )

        ambiguous = (
            crosswalk.group_by(["source", "evidence_type"])
            .agg([
                pl.col("merge_evidence_family").n_unique().alias("families"),
                pl.col("events").sum().alias("events"),
            ])
            .filter(pl.col("families") > 1)
            .sort("events", descending=True)
        )
        progress.update(1)

        # Hard invariants: neither data product may lose records.
        if int(overall["events"].sum()) != event_total:
            fail("Evidence-type event totals are not conserved")
        if int(overall["assertions"].sum()) != assertion_total:
            fail("Evidence-type assertion totals are not conserved")
        if int(family["events"].sum()) != event_total:
            fail("Evidence-family event totals are not conserved")
        if int(family_by_source["events"].sum()) != event_total:
            fail("Source-family event totals are not conserved")
        if without_ctd_n + manifest_source_counts["CTD"] != event_total:
            fail("Non-CTD sensitivity totals do not reconcile")

        if crosswalk.filter(
            pl.col("source").is_null()
            | (pl.col("source").str.strip_chars() == "")
            | pl.col("evidence_type").is_null()
            | (pl.col("evidence_type").str.strip_chars() == "")
            | pl.col("merge_evidence_family").is_null()
            | (pl.col("merge_evidence_family").str.strip_chars() == "")
        ).height:
            fail("Missing source/evidence/family labels detected")
        progress.update(1)

    display("1. SOURCE-LEVEL EVIDENCE COMPOSITION", sources)
    display("2. OVERALL ORIGINAL evidence_type DISTRIBUTION", overall)
    display("3. ORIGINAL evidence_type WITHIN EACH SOURCE", by_source_type)
    display("4. RECLASSIFIED EVIDENCE FAMILIES (EVENTS)", family)
    display("5. EVIDENCE FAMILIES BY SOURCE", family_by_source)
    display("6. EVIDENCE FAMILIES EXCLUDING CTD", families_without_ctd)
    display(
        "7. ORIGINAL -> RECLASSIFIED LABEL CROSSWALK",
        crosswalk.sort(["source", "events"], descending=[False, True]),
    )

    if ambiguous.height:
        display(
            "CHECK: ONE ORIGINAL LABEL WITH MULTIPLE ANALYTICAL FAMILIES",
            ambiguous,
        )
    else:
        console.print(
            "[green]Each (source, evidence_type) maps to one "
            "analytical family in this snapshot.[/green]"
        )

    console.rule("[bold]RECONCILIATION / MANUSCRIPT NOTES[/bold]")
    console.print(
        f"[green]PASS[/green] 11 sources; "
        f"{event_total:,} preserved evidence events"
    )
    console.print(
        f"[green]PASS[/green] "
        f"{assertion_total:,} distinct source-specific assertions"
    )
    console.print(
        f"[green]PASS[/green] "
        f"{event_total - assertion_total:,} "
        "additional supporting events (not deleted)"
    )
    console.print(
        "[green]PASS[/green] original evidence labels "
        "and analytical families conserved"
    )
    console.print(
        f"CTD: {manifest_source_counts['CTD']:,} events "
        f"({100 * manifest_source_counts['CTD'] / event_total:.2f}% "
        "of all events)"
    )

    other = family.filter(
        pl.col("merge_evidence_family") == "other_or_unclassified"
    )
    if other.height:
        console.print(
            f"[yellow]REVIEW[/yellow] {int(other['events'][0]):,} "
            "other/unclassified events; see crosswalk."
        )

    console.print(
        "Interpretation: evidence events are not unique drugs, "
        "pairs, independent studies or approved indications."
    )
    console.print(
        "Assertion percentages use 4,534,385 as denominator; "
        "event percentages use 10,895,395."
    )

    first = overall.row(0, named=True)
    console.print(
        f"Manuscript-ready: Across 11 sources, {event_total:,} "
        f"source evidence events were summarised into "
        f"{assertion_total:,} distinct source-specific assertions. "
        f"The most frequent original evidence type was "
        f"'{first['evidence_type']}' "
        f"({int(first['events']):,} events; "
        f"{float(first['pct_all_events']):.3f}%). "
        "The original evidence labels were preserved and separately "
        "grouped into analytical evidence families without equating "
        "trial, inference, label, contraindication or approval evidence."
    )

    if args.save_tables:
        target = base.parent.parent / "tables"
        target.mkdir(parents=True, exist_ok=True)

        results = {
            "sources": sources,
            "evidence_types": overall,
            "evidence_types_by_source": by_source_type,
            "evidence_families": family,
            "families_by_source": family_by_source,
            "families_without_ctd": families_without_ctd,
            "type_to_family_crosswalk": crosswalk,
            "ambiguous_type_families": ambiguous,
        }

        for name, frame in results.items():
            output = target / f"Analysis_02_{name}.csv"
            frame.write_csv(output)
            console.print(f"[green]SAVED[/green] {output}")

    console.rule(
        "[green]ANALYSIS 02 COMPLETE; LOCKED INPUTS UNCHANGED[/green]"
    )


if __name__ == "__main__":
    main()