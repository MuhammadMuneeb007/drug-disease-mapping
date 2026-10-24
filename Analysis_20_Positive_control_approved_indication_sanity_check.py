#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Analysis 20: source-aware positive-control indication retrieval.

Tests the original 20 drug-disease controls without updating the locked merge.
Separates strict synonym/identifier + disease-label evidence from contextual
name mentions and non-indication associations. Does NOT validate approvals.

Run:
 POLARS_MAX_THREADS=16 python Analysis_20_Positive_control_approved_indication_sanity_check.py --save-tables
"""
from __future__ import annotations

import argparse
import csv
import json
import re
from pathlib import Path

import polars as pl
from tqdm import tqdm

EXPECTED_ASSERTIONS = 4_534_385
EXPECTED_EVENTS = 10_895_395
SOURCES = {"AACT", "Broad Drug Repurposing Hub", "ChEMBL", "CTD", "DrugCentral",
           "Hetionet", "MEDI", "OpenTargets", "PrimeKG", "repoDB", "SIDER"}
INDICATION_FAMILIES = {"source_reported_indication", "source_reported_approval",
                       "ctd_direct_therapeutic", "label_nlp_indication"}
CAUTION_FAMILIES = {"contraindication", "off_label_use", "discontinued_development"}

# User-provided control definitions preserved from the original Analysis 20.
# Each row: ID | drug | expected condition | drug aliases | drug IDs | disease aliases
# Aliases and identifiers are declared test inputs, NOT independently verified here.
CONTROL_DATA = """PC01|Metformin|Type 2 diabetes mellitus|metformin;metformin hydrochloride|DB00331;CHEMBL1431|type 2 diabetes;type 2 diabetes mellitus;type ii diabetes;type ii diabetes mellitus;diabetes mellitus type 2;diabetes mellitus, type 2;non-insulin-dependent diabetes;non insulin dependent diabetes;t2d;t2dm
PC02|Lisinopril|Hypertension|lisinopril|DB00722;CHEMBL1237|hypertension;high blood pressure;essential hypertension
PC03|Sildenafil|Erectile dysfunction|sildenafil;sildenafil citrate;viagra|DB00203;CHEMBL192|erectile dysfunction;male erectile dysfunction;impotence
PC04|Imatinib|Chronic myeloid leukemia|imatinib;imatinib mesylate;gleevec;glivec|DB00619;CHEMBL941|chronic myeloid leukemia;chronic myelogenous leukemia;cml;philadelphia chromosome positive chronic myeloid leukemia;ph cml
PC05|Atorvastatin|Hypercholesterolemia|atorvastatin;atorvastatin calcium;lipitor|DB01076;CHEMBL1487|hypercholesterolemia;hypercholesterolaemia;high cholesterol;hyperlipidemia;dyslipidemia
PC06|Warfarin|Thrombosis|warfarin;warfarin sodium;coumadin|DB00682;CHEMBL1464|thrombosis;venous thrombosis;deep vein thrombosis;pulmonary embolism;thromboembolism;embolism and thrombosis
PC07|Albuterol|Asthma|albuterol;salbutamol;albuterol sulfate|DB01001;CHEMBL714|asthma;bronchial asthma;bronchospasm;reversible obstructive airway disease
PC08|Omeprazole|Gastroesophageal reflux disease|omeprazole;prilosec|DB00338;CHEMBL1503|gastroesophageal reflux disease;gastro-oesophageal reflux disease;gastroesophageal reflux;gastro-oesophageal reflux;gerd;reflux esophagitis;reflux oesophagitis
PC09|Levothyroxine|Hypothyroidism|levothyroxine;levothyroxine sodium;thyroxine|DB00451;CHEMBL1624|hypothyroidism;thyroid hormone deficiency;underactive thyroid
PC10|Sertraline|Major depressive disorder|sertraline;sertraline hydrochloride;zoloft|DB01104;CHEMBL809|major depressive disorder;major depression;depression;depressive disorder
PC11|Fluoxetine|Major depressive disorder|fluoxetine;fluoxetine hydrochloride;prozac|DB00472;CHEMBL41|major depressive disorder;major depression;depression;depressive disorder
PC12|Amlodipine|Hypertension|amlodipine;amlodipine besylate;amlodipine besilate;norvasc|DB00381;CHEMBL1491|hypertension;high blood pressure;essential hypertension
PC13|Losartan|Hypertension|losartan;losartan potassium;cozaar|DB00678;CHEMBL191|hypertension;high blood pressure;essential hypertension
PC14|Ciprofloxacin|Urinary tract infection|ciprofloxacin;ciprofloxacin hydrochloride;cipro|DB00537;CHEMBL8|urinary tract infection;urinary tract infections;uti;cystitis;pyelonephritis
PC15|Oseltamivir|Influenza|oseltamivir;oseltamivir phosphate;tamiflu|DB00198;CHEMBL1229|influenza;flu;influenza a;influenza b
PC16|Aciclovir|Herpes simplex virus infection|aciclovir;acyclovir;acyclovir sodium;zovirax|DB00787;CHEMBL184|herpes simplex;herpes simplex infection;herpes simplex virus infection;genital herpes;cold sores
PC17|Methotrexate|Rheumatoid arthritis|methotrexate;methotrexate sodium|DB00563;CHEMBL34259|rheumatoid arthritis;ra
PC18|Adalimumab|Rheumatoid arthritis|adalimumab;humira|DB00051;CHEMBL1201580|rheumatoid arthritis;ra
PC19|Tamoxifen|Breast cancer|tamoxifen;tamoxifen citrate;nolvadex|DB00675;CHEMBL83|breast cancer;breast carcinoma;malignant neoplasm of breast
PC20|Trastuzumab|HER2-positive breast cancer|trastuzumab;herceptin|DB00072;CHEMBL1201585|her2 positive breast cancer;erbb2 positive breast cancer;breast cancer;breast carcinoma
"""

# Broad or ambiguous disease aliases are permitted only as contextual matches.
# Contextual matches cannot generate a strict positive-control pass.
CONTEXT_ONLY = {
    "PC03": {"impotence"},
    "PC05": {"hyperlipidemia", "dyslipidemia"},
    "PC06": {"embolism and thrombosis", "embolism"},
    "PC07": {"bronchospasm", "reversible obstructive airway disease"},
    "PC08": {"reflux esophagitis", "reflux oesophagitis"},
    "PC10": {"depression", "depressive disorder"},
    "PC11": {"depression", "depressive disorder"},
    "PC14": {"cystitis", "pyelonephritis"},
    "PC16": {"genital herpes", "cold sores"},
    "PC17": {"ra"},
    "PC18": {"ra"},
    "PC20": {"breast cancer", "breast carcinoma"},
}

ASSERT_FILE = "ALL_INCLUDED_SOURCES_drug_disease_deduplicated.csv"
EVENT_FILE = "ALL_INCLUDED_SOURCES_drug_disease_merged.csv"
MANIFEST_FILE = "ALL_INCLUDED_SOURCES_merge_summary.json"
ASSERT_FIELDS = [
    "merge_assertion_id", "source", "drug_name", "drug_identifier",
    "drug_identifier_type", "disease_or_condition_name",
    "disease_or_condition_identifier", "relationship_type", "evidence_type",
    "supporting_event_count",
]
EVENT_FIELDS = [
    "merge_assertion_id", "source", "relationship_type", "evidence_type",
    "merge_evidence_family", "merge_event_id",
]


def name_key(value):
    """Exact Analysis 01/04/08 name-key rule; not an ontology mapping."""
    return re.sub(r"\s+", " ", str(value or "").strip()).lower()


def expr_name_key(name):
    return (pl.col(name).fill_null("").str.strip_chars()
            .str.replace_all(r"\s+", " ").str.to_lowercase())


def normalized_identifier(name):
    return (pl.col(name).fill_null("").str.strip_chars()
            .str.replace(r"(?i)^(?:drugbank|chembl)\s*:\s*", "")
            .str.to_uppercase())


def normalized_namespace(name):
    return (pl.col(name).fill_null("").str.to_lowercase()
            .str.replace_all(r"[^a-z0-9]+", ""))


def qualified_ids(control):
    return {
        "drugbank": [x.upper() for x in control["drug_ids"] if x.upper().startswith("DB")],
        "chembl": [x.upper() for x in control["drug_ids"] if x.upper().startswith("CHEMBL")],
    }


def id_match_expr(control):
    ids = qualified_ids(control)
    return ((pl.col("_ns") == "drugbank") & pl.col("_id").is_in(ids["drugbank"])) | (
        (pl.col("_ns") == "chembl") & pl.col("_id").is_in(ids["chembl"])
    )


def parse_controls():
    result = []
    for line in CONTROL_DATA.splitlines():
        fields = line.strip().split("|")
        if len(fields) != 6:
            raise RuntimeError(f"Invalid positive-control specification: {line}")
        cid, drug, disease, drug_aliases, ids, disease_aliases = fields
        da = {name_key(x) for x in drug_aliases.split(";") if x.strip()}
        dd = {name_key(x) for x in disease_aliases.split(";") if x.strip()}
        strict = dd - {name_key(x) for x in CONTEXT_ONLY.get(cid, set())}
        if cid == "PC20":
            strict |= {"her2-positive breast cancer", "her2 positive breast cancer",
                       "erbb2-positive breast cancer", "erbb2 positive breast cancer",
                       "her2-positive breast carcinoma", "her2 positive breast carcinoma"}
            dd |= strict
        if not da or not dd or not strict:
            raise RuntimeError(f"Empty alias set for {cid}")
        result.append({"control_id": cid, "drug": drug, "disease": disease,
                       "drug_aliases": da, "drug_ids": ids.split(";"),
                       "disease_aliases": dd, "strict_disease_aliases": strict})
    if len(result) != 20 or len({c["control_id"] for c in result}) != 20:
        raise RuntimeError("Expected 20 distinct positive-control definitions")
    return result


def validate_header(path, required):
    if not path.is_file() or path.stat().st_size == 0:
        raise RuntimeError(f"Missing or empty input: {path}")
    with path.open(newline="", encoding="utf-8-sig") as f:
        header = next(csv.reader(f, strict=True), [])
    if len(header) != len(set(header)) or not set(required).issubset(header):
        raise RuntimeError(f"Schema mismatch: {path.name}; missing={set(required)-set(header)}")


def scan_csv(path, cols):
    return pl.scan_csv(str(path), infer_schema=False, ignore_errors=False,
                       truncate_ragged_lines=False, encoding="utf8").select(cols)


def require(label, observed, expected):
    if observed != expected:
        raise RuntimeError(f"[FAIL] {label}: observed={observed!r}; expected={expected!r}")
    print(f"[PASS] {label}: {observed:,}" if isinstance(observed, int)
          else f"[PASS] {label}", flush=True)


def print_table(title, frame):
    print("\n" + "=" * 105)
    print(f"{title} [{frame.height:,} rows]")
    print("=" * 105)
    if not frame.height:
        print("[EMPTY]")
        return
    print("\t".join(frame.columns))
    for row in frame.iter_rows():
        print("\t".join(
            ("" if x is None else f"{x:.3f}" if isinstance(x, float) else str(x))
            .replace("\t", " ").replace("\r", " ").replace("\n", " ")
            for x in row
        ))


def contextual_match_regex(aliases):
    # Token boundaries prevent 'flu' matching influenza, fluorine, etc.
    # Abbreviations <=3 characters require full-label equality.
    values = [x for x in aliases if len(x) > 3]
    if not values:
        return None
    parts = sorted((re.escape(x) for x in values), key=len, reverse=True)
    return re.compile(r"(?<![a-z0-9])(?:" + "|".join(parts) + r")(?![a-z0-9])")


def matched_disease_tiers(values, control):
    strict = control["strict_disease_aliases"]
    aliases = control["disease_aliases"]
    contextual = contextual_match_regex(aliases)
    tiers = {}
    for term in values:
        if term in strict:
            tiers[term] = "strict"
        elif term in aliases or (contextual is not None and contextual.search(term)):
            tiers[term] = "contextual"
    return tiers


def rows_summary(frame, by, count_field="events"):
    if not frame.height:
        return pl.DataFrame({**{k: [] for k in by}, count_field: []})
    return (frame.group_by(by).agg(pl.col("events").sum().alias(count_field))
            .sort(count_field, descending=True))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--processed-dir", type=Path,
                    default=Path(__file__).resolve().parent / "data" / "processed")
    ap.add_argument("--save-tables", action="store_true")
    ap.add_argument("--examples", type=int, default=15)
    ap.add_argument("--min-strict-passes", type=int, default=0,
                    help="Optional scientific gate (default 0, report rather than fail)")
    args = ap.parse_args()
    if not 1 <= args.examples <= 100 or not 0 <= args.min_strict_passes <= 20:
        ap.error("--examples 1..100; --min-strict-passes 0..20")

    base = args.processed_dir.expanduser().resolve()
    af = base / ASSERT_FILE
    ef = base / EVENT_FILE
    mf = base / MANIFEST_FILE
    validate_header(af, ASSERT_FIELDS)
    validate_header(ef, EVENT_FIELDS)
    if not mf.is_file():
        raise RuntimeError(f"Missing locked manifest: {mf}")
    manifest = json.loads(mf.read_text(encoding="utf-8"))
    require("merge status", manifest.get("status"), "success")
    require("manifest assertions", int(manifest["assertions"]), EXPECTED_ASSERTIONS)
    require("manifest events", int(manifest["source_events"]), EXPECTED_EVENTS)
    controls = parse_controls()

    print("\nANALYSIS 20 | POSITIVE-CONTROL SOURCE-INDICATION SANITY CHECK")
    print("[POLICY] Strict disease label matches are distinguished from contextual matches.")
    print("[POLICY] Source-reported indication does not mean approval or efficacy.")

    with tqdm(total=4, desc="Analysis 20", unit="stage") as progress:
        print("[STEP 1] Validate assertion source totals", flush=True)
        narrow = scan_csv(af, ASSERT_FIELDS)
        counts = (narrow.select([
            pl.len().alias("assertions"),
            pl.col("supporting_event_count").cast(pl.Int64, strict=True)
            .sum().alias("events"),
        ]).collect(engine="streaming").row(0, named=True))
        require("locked assertion row count", int(counts["assertions"]), EXPECTED_ASSERTIONS)
        require("assertion-linked events", int(counts["events"]), EXPECTED_EVENTS)
        progress.update(1)

        print("[STEP 2] Stream selected drug candidates once", flush=True)
        name_aliases = sorted(set().union(*(c["drug_aliases"] for c in controls)))
        db_ids = sorted({x for c in controls for x in qualified_ids(c)["drugbank"]})
        chembl_ids = sorted({x for c in controls for x in qualified_ids(c)["chembl"]})
        names = (narrow.with_columns([
            expr_name_key("drug_name").alias("_drug"),
            expr_name_key("disease_or_condition_name").alias("_disease"),
            normalized_identifier("drug_identifier").alias("_id"),
            normalized_namespace("drug_identifier_type").alias("_ns"),
            pl.col("supporting_event_count").cast(pl.Int64, strict=True),
        ]))
        candidates = names.filter(
            pl.col("_drug").is_in(name_aliases)
            | ((pl.col("_ns") == "drugbank") & pl.col("_id").is_in(db_ids))
            | ((pl.col("_ns") == "chembl") & pl.col("_id").is_in(chembl_ids))
        ).collect(engine="streaming")
        if candidates.filter(~pl.col("source").is_in(sorted(SOURCES))).height:
            raise RuntimeError("Unknown source label in drug candidate assertions")
        if candidates["merge_assertion_id"].n_unique() != candidates.height:
            raise RuntimeError("Duplicate assertion IDs in candidate subset")
        print(f"[INFO] Drug candidate assertions: {candidates.height:,}", flush=True)

        matched_parts = []
        for control in controls:
            drugs = candidates.filter(
                pl.col("_drug").is_in(sorted(control["drug_aliases"]))
                | id_match_expr(control)
            )
            if not drugs.height:
                continue
            tiers = matched_disease_tiers(drugs["_disease"].unique().to_list(), control)
            if not tiers:
                continue
            part = (drugs.filter(pl.col("_disease").is_in(list(tiers)))
                    .with_columns([
                        pl.lit(control["control_id"]).alias("control_id"),
                        pl.when(pl.col("_disease").is_in(sorted(control["strict_disease_aliases"])))
                        .then(pl.lit("strict")).otherwise(pl.lit("contextual"))
                        .alias("disease_match_tier"),
                    ]))
            matched_parts.append(part)

        matches = (pl.concat(matched_parts, how="vertical") if matched_parts
                   else candidates.head(0).with_columns([
                       pl.lit("").alias("control_id"),
                       pl.lit("").alias("disease_match_tier")]))
        print(f"[INFO] Matched (control, assertion) rows: {matches.height:,}", flush=True)
        progress.update(1)

        print("[STEP 3] Reconcile original evidence-event support", flush=True)
        ids = matches["merge_assertion_id"].unique().to_list()
        if ids:
            grouped_events = (scan_csv(ef, EVENT_FIELDS)
                .filter(pl.col("merge_assertion_id").is_in(ids))
                .group_by(["merge_assertion_id", "source", "relationship_type",
                           "evidence_type", "merge_evidence_family"])
                .agg([
                    pl.len().alias("events"),
                    pl.col("merge_event_id").n_unique().alias("distinct_event_ids"),
                ]).collect(engine="streaming"))
            if grouped_events.filter(pl.col("events") != pl.col("distinct_event_ids")).height:
                raise RuntimeError("Duplicate event IDs within evidence family groups")
            summed = (grouped_events.group_by("merge_assertion_id")
                      .agg(pl.col("events").sum().alias("actual_events")))
            uniq_assertions = matches.select([
                "merge_assertion_id", "supporting_event_count"
            ]).unique()
            audit = uniq_assertions.join(summed, on="merge_assertion_id", how="full", coalesce=True)
            bad = audit.filter(
                pl.col("actual_events").is_null()
                | pl.col("supporting_event_count").is_null()
                | (pl.col("supporting_event_count") != pl.col("actual_events"))
            )
            require("unreconciled matched assertions", bad.height, 0)
            require("matched event weight conservation",
                    int(uniq_assertions["supporting_event_count"].sum()),
                    int(grouped_events["events"].sum()))
            joined = matches.join(grouped_events, on="merge_assertion_id", how="inner",
                                  suffix="_original")
            if joined.filter(pl.col("source") != pl.col("source_original")).height:
                raise RuntimeError("Event/assertion source inconsistency")
        else:
            joined = pl.DataFrame({
                "control_id": [], "disease_match_tier": [], "source": [],
                "merge_evidence_family": [], "relationship_type_original": [],
                "evidence_type_original": [], "events": [],
            })
        progress.update(1)

        print("[STEP 4] Summarise strict and contextual indications", flush=True)
        results = []
        for control in controls:
            cid = control["control_id"]
            m = matches.filter(pl.col("control_id") == cid)
            e = joined.filter(pl.col("control_id") == cid)
            strict = m.filter(pl.col("disease_match_tier") == "strict")
            strict_e = e.filter(pl.col("disease_match_tier") == "strict")
            indicated = e.filter(pl.col("merge_evidence_family").is_in(sorted(INDICATION_FAMILIES)))
            strict_indicated = strict_e.filter(
                pl.col("merge_evidence_family").is_in(sorted(INDICATION_FAMILIES)))
            caution = e.filter(pl.col("merge_evidence_family").is_in(sorted(CAUTION_FAMILIES)))
            approval = strict_e.filter(pl.col("merge_evidence_family") == "source_reported_approval")
            if strict_indicated.height:
                status = "PASS_STRICT_INDICATION_EVIDENCE"
            elif indicated.height:
                status = "PARTIAL_CONTEXTUAL_INDICATION_ONLY"
            elif m.height:
                status = "PARTIAL_ASSOCIATION_ONLY"
            else:
                status = "FAIL_NOT_RETRIEVED"
            results.append({
                "control_id": cid,
                "expected_drug": control["drug"],
                "expected_disease": control["disease"],
                "status": status,
                "matched_assertions": m.height,
                "strict_disease_assertions": strict.height,
                "retrieved_evidence_events": int(m["supporting_event_count"].sum() or 0),
                "strict_indication_events": int(strict_indicated["events"].sum() or 0),
                "contextual_indication_events": int(indicated["events"].sum() or 0) - int(strict_indicated["events"].sum() or 0),
                "source_reported_approval_events": int(approval["events"].sum() or 0),
                "contraindication_events": int(e.filter(pl.col("merge_evidence_family") == "contraindication")["events"].sum() or 0),
                "caution_family_events": int(caution["events"].sum() or 0),
                "retrieving_sources": m["source"].n_unique() if m.height else 0,
                "strict_indication_sources": strict_indicated["source"].n_unique() if strict_indicated.height else 0,
                "sources_present": " | ".join(sorted(m["source"].unique().to_list())),
                "strict_indication_source_names": " | ".join(sorted(strict_indicated["source"].unique().to_list())),
            })
        table = pl.DataFrame(results)
        require("positive controls evaluated", table.height, 20)
        require("status classification conservation", table.group_by("status").len()["len"].sum(), 20)
        require("joined evidence weights by matched control",
                int(table["retrieved_evidence_events"].sum()),
                int(joined["events"].sum()) if joined.height else 0)

        distribution = (table.group_by("status").len().rename({"len": "controls"})
                        .with_columns((100 * pl.col("controls") / 20).round(1).alias("pct_controls"))
                        .sort("controls", descending=True))
        per_source = rows_summary(joined, ["control_id", "source"])
        family = rows_summary(joined, ["control_id", "disease_match_tier", "merge_evidence_family"])
        relationship = rows_summary(joined, ["control_id", "disease_match_tier", "relationship_type_original"])
        evidence = rows_summary(joined, ["control_id", "disease_match_tier", "evidence_type_original"])
        matched_diseases = (matches.group_by([
            "control_id", "disease_match_tier", "disease_or_condition_name",
        ]).agg([
            pl.len().alias("assertions"),
            pl.col("supporting_event_count").sum().alias("events"),
        ]).sort(["control_id", "events"], descending=[False, True])
            if matches.height else pl.DataFrame({"control_id": [], "disease_match_tier": [],
                                                 "disease_or_condition_name": [], "assertions": [], "events": []}))
        examples = (matches.sort(["control_id", "supporting_event_count"],
                                 descending=[False, True]).group_by("control_id", maintain_order=True)
                    .head(args.examples).select([
                        "control_id", "disease_match_tier", "drug_name", "drug_identifier",
                        "drug_identifier_type", "disease_or_condition_name", "source",
                        "relationship_type", "evidence_type", "supporting_event_count",
                    ]) if matches.height else matches.head(0))
        progress.update(1)

    outputs = {
        "controls": table,
        "status_distribution": distribution,
        "source_support_by_control": per_source,
        "family_by_control_and_match_tier": family,
        "relationship_by_control_and_match_tier": relationship,
        "evidence_by_control_and_match_tier": evidence,
        "matched_disease_labels": matched_diseases,
        "example_assertions": examples,
    }
    for index, (name, frame) in enumerate(outputs.items(), 1):
        print_table(f"{index}. {name.replace('_', ' ').upper()}", frame)

    n_pass = table.filter(pl.col("status") == "PASS_STRICT_INDICATION_EVIDENCE").height
    n_context = table.filter(pl.col("status") == "PARTIAL_CONTEXTUAL_INDICATION_ONLY").height
    n_assoc = table.filter(pl.col("status") == "PARTIAL_ASSOCIATION_ONLY").height
    n_absent = table.filter(pl.col("status") == "FAIL_NOT_RETRIEVED").height
    print("\n[MANUSCRIPT SUMMARY]")
    print(f"Of 20 predefined positive controls, {n_pass} were retrieved with strict "
          f"drug/condition matching and source-reported indication-oriented evidence; "
          f"{n_context} had such evidence only for broader contextual condition labels; "
          f"{n_assoc} were retrieved without indication-oriented evidence; "
          f"and {n_absent} were not retrieved under the declared matching rules. "
          "The outcomes describe evidence retrieval, not independent verification "
          "of marketing approval, regulatory indication or treatment efficacy.")
    print("[CAUTION] Broad context terms never confer strict PASS (including generic breast cancer for PC20).")
    print("[CAUTION] Multiple databases may reuse upstream evidence; no independence assumed.")
    print("[PASS] ANALYSIS 20 COMPUTATION COMPLETE; LOCKED INPUTS UNCHANGED")

    if args.save_tables:
        folder = base.parent.parent / "tables"
        folder.mkdir(parents=True, exist_ok=True)
        for key, frame in outputs.items():
            target = folder / f"Analysis_20_{key}.csv"
            frame.write_csv(target)
            print("[SAVED]", target)
        report = {
            "snapshot": "paper11-locked-v1", "manifest_assertions": EXPECTED_ASSERTIONS,
            "manifest_events": EXPECTED_EVENTS, "controls": 20,
            "strict_passes": n_pass, "context_only": n_context,
            "association_only": n_assoc, "not_found": n_absent,
            "indication_families": sorted(INDICATION_FAMILIES),
            "polars_version": pl.__version__,
            "method": "Conservative drug name/typed-ID matching, strict vs contextual disease labels",
        }
        target = folder / "Analysis_20_summary.json"
        target.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
        print("[SAVED]", target)

    if n_pass < args.min_strict_passes:
        raise RuntimeError(
            f"Strict pass gate failed: {n_pass}/{20} < {args.min_strict_passes}. "
            "Review control statuses instead of changing the locked dataset."
        )


if __name__ == "__main__":
    main()
