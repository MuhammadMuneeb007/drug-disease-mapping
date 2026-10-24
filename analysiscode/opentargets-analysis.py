#!/usr/bin/env python3
"""
Open Targets Analysis: Exploratory analysis of Open Targets drug-disease output.

Input:
  data/processed/opentargets_drug_disease.csv

Expected columns:
  drug_id
  drug_name
  disease_id
  disease_name
  max_phase
  clinical_phase_raw
  clinical_report_ids
  opentargets_indication_id
  source
  dataset
  release
  evidence_type

Purpose:
  1. Validate final processed Open Targets output.
  2. Summarise drug-disease clinical indication coverage.
  3. Identify drugs linked to the most diseases.
  4. Identify diseases linked to the most drugs.
  5. Analyse maximum clinical stage / phase distribution.
  6. Analyse clinical report provenance support.
  7. Identify broad disease-category landscape.
  8. Export tables, figures, report, and metadata.

Important interpretation:
  Open Targets clinical_indication records represent drug-disease clinical
  indication/clinical precedence-style relationships. They are stronger than
  broad chemical-disease association resources, but they should still be
  interpreted with clinical stage, report provenance, and disease ontology
  mapping preserved.

Requirements:
  pip install pandas matplotlib tabulate
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd
import matplotlib.pyplot as plt


REQUIRED_COLUMNS = [
    "drug_id",
    "drug_name",
    "disease_id",
    "disease_name",
    "max_phase",
    "clinical_phase_raw",
    "clinical_report_ids",
    "opentargets_indication_id",
    "source",
    "dataset",
    "release",
    "evidence_type",
]


# =============================================================================
# Helpers
# =============================================================================

def now_utc() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256_file(path: Path, chunk_size: int = 1024 * 1024) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(chunk_size), b""):
            h.update(chunk)
    return h.hexdigest()


def human_size(num_bytes: int) -> str:
    size = float(num_bytes)
    for unit in ["B", "KB", "MB", "GB", "TB"]:
        if size < 1024:
            return f"{size:.2f} {unit}"
        size /= 1024
    return f"{size:.2f} PB"


def print_section(title: str) -> None:
    print("\n" + "=" * 100)
    print(title)
    print("=" * 100)


def clean_text(x: object) -> str:
    if x is None:
        return ""
    try:
        if pd.isna(x):
            return ""
    except Exception:
        pass
    return str(x).strip()


def normalise_text(x: object) -> str:
    x = clean_text(x).lower()
    x = re.sub(r"\s+", " ", x)
    return x.strip()


def count_semicolon_items(x: object) -> int:
    x = clean_text(x)
    if not x:
        return 0
    return len([v for v in x.split(";") if v.strip()])


def save_table(df: pd.DataFrame, out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_path, index=False)
    print(f"[SAVED] {out_path}")


def barplot(
    df: pd.DataFrame,
    x_col: str,
    y_col: str,
    title: str,
    out_path: Path,
    top_n: int = 20,
) -> None:
    if df.empty:
        print(f"[SKIPPED FIGURE] Empty table: {out_path}")
        return

    plot_df = df.head(top_n).copy()
    plot_df = plot_df.iloc[::-1]

    plt.figure(figsize=(10, max(5, 0.35 * len(plot_df))))
    plt.barh(plot_df[x_col].astype(str), plot_df[y_col])
    plt.xlabel(y_col)
    plt.ylabel(x_col)
    plt.title(title)
    plt.tight_layout()

    out_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_path, dpi=300)
    plt.close()

    print(f"[FIGURE] {out_path}")


# =============================================================================
# Validation and cleaning
# =============================================================================

def validate_input(df: pd.DataFrame, input_path: Path) -> dict:
    print_section("1. INPUT VALIDATION")

    print(f"[INPUT FILE] {input_path}")
    print(f"[SIZE] {human_size(input_path.stat().st_size)}")
    print(f"[SHA256] {sha256_file(input_path)}")
    print(f"[ROWS, COLUMNS] {df.shape}")

    print("\n[COLUMNS]")
    print(list(df.columns))

    missing_required = [c for c in REQUIRED_COLUMNS if c not in df.columns]
    if missing_required:
        raise ValueError(f"Missing required columns: {missing_required}")

    print("\n[HEAD]")
    preview_cols = [
        "drug_id",
        "drug_name",
        "disease_id",
        "disease_name",
        "max_phase",
        "clinical_phase_raw",
        "source",
        "release",
    ]
    print(df[preview_cols].head(10).to_string(index=False))

    missing_df = (
        df.isna()
        .sum()
        .reset_index()
        .rename(columns={"index": "column", 0: "missing_count"})
    )
    missing_df["missing_percent"] = (missing_df["missing_count"] / len(df) * 100).round(3)

    empty_df = (
        df.astype(str)
        .apply(lambda s: s.str.strip().eq("").sum())
        .reset_index()
        .rename(columns={"index": "column", 0: "empty_string_count"})
    )
    empty_df["empty_string_percent"] = (empty_df["empty_string_count"] / len(df) * 100).round(3)

    print("\n[MISSING VALUES]")
    print(missing_df.to_string(index=False))

    print("\n[EMPTY STRINGS]")
    print(empty_df.to_string(index=False))

    return {
        "input_file": str(input_path),
        "input_size_bytes": input_path.stat().st_size,
        "input_sha256": sha256_file(input_path),
        "n_rows": int(df.shape[0]),
        "n_columns": int(df.shape[1]),
        "columns": list(df.columns),
        "missing_values": missing_df.to_dict(orient="records"),
        "empty_strings": empty_df.to_dict(orient="records"),
    }


def add_normalised_columns(df: pd.DataFrame) -> pd.DataFrame:
    print_section("2. CLEANING / NORMALISED ANALYSIS COLUMNS")

    df = df.copy()

    for col in REQUIRED_COLUMNS:
        if col in df.columns and col != "max_phase":
            df[col] = df[col].fillna("").astype(str).str.strip()

    df["drug_id_norm"] = df["drug_id"].map(normalise_text)
    df["drug_name_norm"] = df["drug_name"].map(normalise_text)
    df["disease_id_norm"] = df["disease_id"].map(normalise_text)
    df["disease_name_norm"] = df["disease_name"].map(normalise_text)

    df["max_phase_num"] = pd.to_numeric(df["max_phase"], errors="coerce").fillna(0).astype(int)
    df["clinical_phase_raw_clean"] = df["clinical_phase_raw"].replace("", "Missing")

    df["clinical_report_count"] = df["clinical_report_ids"].map(count_semicolon_items)
    df["has_clinical_report_ids"] = df["clinical_report_count"].gt(0)

    df["is_approved_or_phase4"] = df["max_phase_num"].ge(4)
    df["is_phase3_or_higher"] = df["max_phase_num"].ge(3)
    df["is_phase2_or_higher"] = df["max_phase_num"].ge(2)

    df["has_drug_name"] = df["drug_name"].ne("") & df["drug_name"].ne(df["drug_id"])
    df["has_disease_name"] = df["disease_name"].ne("") & df["disease_name"].ne(df["disease_id"])

    print("[HEAD AFTER NORMALISATION]")
    show_cols = [
        "drug_id",
        "drug_name",
        "disease_id",
        "disease_name",
        "max_phase_num",
        "clinical_phase_raw_clean",
        "clinical_report_count",
    ]
    print(df[show_cols].head(10).to_string(index=False))

    return df


# =============================================================================
# Analysis tables
# =============================================================================

def overview_stats(df: pd.DataFrame) -> pd.DataFrame:
    print_section("3. OVERVIEW STATISTICS")

    stats = {
        "total_rows": len(df),
        "unique_drug_ids": df["drug_id_norm"].replace("", pd.NA).nunique(dropna=True),
        "unique_drug_names": df["drug_name_norm"].replace("", pd.NA).nunique(dropna=True),
        "unique_disease_ids": df["disease_id_norm"].replace("", pd.NA).nunique(dropna=True),
        "unique_disease_names": df["disease_name_norm"].replace("", pd.NA).nunique(dropna=True),
        "unique_drug_disease_pairs": df[["drug_id_norm", "disease_id_norm"]].drop_duplicates().shape[0],
        "approved_or_phase4_records": int(df["is_approved_or_phase4"].sum()),
        "phase3_or_higher_records": int(df["is_phase3_or_higher"].sum()),
        "phase2_or_higher_records": int(df["is_phase2_or_higher"].sum()),
        "records_with_clinical_report_ids": int(df["has_clinical_report_ids"].sum()),
        "mean_clinical_reports_per_record": round(float(df["clinical_report_count"].mean()), 3),
        "median_clinical_reports_per_record": round(float(df["clinical_report_count"].median()), 3),
        "max_clinical_reports_per_record": int(df["clinical_report_count"].max()),
        "records_with_drug_names": int(df["has_drug_name"].sum()),
        "records_with_disease_names": int(df["has_disease_name"].sum()),
        "exact_duplicate_rows": int(df.duplicated().sum()),
    }

    out = pd.DataFrame(list(stats.items()), columns=["metric", "value"])
    print(out.to_string(index=False))
    return out


def top_drugs_by_unique_diseases(df: pd.DataFrame, top_n: int = 100) -> pd.DataFrame:
    print_section("4. DRUGS LINKED TO THE LARGEST NUMBER OF DISTINCT DISEASES")

    out = (
        df.groupby(["drug_id_norm", "drug_name"])
        .agg(
            unique_disease_ids=("disease_id_norm", "nunique"),
            unique_disease_names=("disease_name_norm", "nunique"),
            total_records=("drug_id_norm", "size"),
            approved_or_phase4_records=("is_approved_or_phase4", "sum"),
            phase3_or_higher_records=("is_phase3_or_higher", "sum"),
            records_with_clinical_reports=("has_clinical_report_ids", "sum"),
            total_clinical_report_ids=("clinical_report_count", "sum"),
            mean_clinical_report_ids=("clinical_report_count", "mean"),
            max_phase=("max_phase_num", "max"),
        )
        .reset_index()
        .rename(columns={"drug_id_norm": "drug_id"})
        .sort_values(["unique_disease_ids", "total_records"], ascending=False)
        .head(top_n)
    )

    out["mean_clinical_report_ids"] = out["mean_clinical_report_ids"].round(3)

    print(out.head(25).to_string(index=False))
    return out


def top_diseases_by_unique_drugs(df: pd.DataFrame, top_n: int = 100) -> pd.DataFrame:
    print_section("5. DISEASES LINKED TO THE LARGEST NUMBER OF DISTINCT DRUGS")

    out = (
        df.groupby(["disease_id_norm", "disease_name"])
        .agg(
            unique_drug_ids=("drug_id_norm", "nunique"),
            unique_drug_names=("drug_name_norm", "nunique"),
            total_records=("disease_id_norm", "size"),
            approved_or_phase4_records=("is_approved_or_phase4", "sum"),
            phase3_or_higher_records=("is_phase3_or_higher", "sum"),
            records_with_clinical_reports=("has_clinical_report_ids", "sum"),
            total_clinical_report_ids=("clinical_report_count", "sum"),
            mean_clinical_report_ids=("clinical_report_count", "mean"),
            max_phase=("max_phase_num", "max"),
        )
        .reset_index()
        .rename(columns={"disease_id_norm": "disease_id"})
        .sort_values(["unique_drug_ids", "total_records"], ascending=False)
        .head(top_n)
    )

    out["mean_clinical_report_ids"] = out["mean_clinical_report_ids"].round(3)

    print(out.head(25).to_string(index=False))
    return out


def most_supported_drug_disease_pairs(df: pd.DataFrame, top_n: int = 100) -> pd.DataFrame:
    print_section("6. MOST SUPPORTED DRUG-DISEASE PAIRS BY CLINICAL REPORT COUNT")

    out = (
        df.groupby(["drug_id", "drug_name", "disease_id", "disease_name"])
        .agg(
            total_records=("drug_id", "size"),
            max_phase=("max_phase_num", "max"),
            clinical_phase_values=("clinical_phase_raw_clean", lambda x: "; ".join(sorted(set(x.astype(str)))[:20])),
            total_clinical_report_ids=("clinical_report_count", "sum"),
            max_clinical_report_ids_single_record=("clinical_report_count", "max"),
            opentargets_indication_ids=("opentargets_indication_id", lambda x: ";".join(sorted(set(x.astype(str)))[:20])),
        )
        .reset_index()
        .sort_values(["total_clinical_report_ids", "max_phase", "total_records"], ascending=False)
        .head(top_n)
    )

    print(out.head(30).to_string(index=False))
    return out


def phase_distribution(df: pd.DataFrame) -> pd.DataFrame:
    print_section("7. CLINICAL PHASE / STAGE DISTRIBUTION")

    out = (
        df.groupby(["max_phase_num", "clinical_phase_raw_clean"])
        .agg(
            total_records=("drug_id_norm", "size"),
            unique_drugs=("drug_id_norm", "nunique"),
            unique_diseases=("disease_id_norm", "nunique"),
            unique_pairs=("drug_id_norm", lambda x: None),
            records_with_clinical_reports=("has_clinical_report_ids", "sum"),
        )
        .reset_index()
        .drop(columns=["unique_pairs"])
        .sort_values(["max_phase_num", "total_records"], ascending=[False, False])
    )

    pair_counts = (
        df.groupby(["max_phase_num", "clinical_phase_raw_clean"])
        .apply(
            lambda x: x[["drug_id_norm", "disease_id_norm"]].drop_duplicates().shape[0],
            include_groups=False,
        )
        .reset_index(name="unique_pairs")
    )

    out = out.merge(pair_counts, on=["max_phase_num", "clinical_phase_raw_clean"], how="left")

    print(out.to_string(index=False))
    return out


def clinical_report_support_summary(df: pd.DataFrame) -> pd.DataFrame:
    print_section("8. CLINICAL REPORT PROVENANCE SUPPORT")

    stats = {
        "records_with_clinical_report_ids": int(df["has_clinical_report_ids"].sum()),
        "records_without_clinical_report_ids": int((~df["has_clinical_report_ids"]).sum()),
        "mean_clinical_report_ids_per_record": df["clinical_report_count"].mean(),
        "median_clinical_report_ids_per_record": df["clinical_report_count"].median(),
        "min_clinical_report_ids_per_record": df["clinical_report_count"].min(),
        "max_clinical_report_ids_per_record": df["clinical_report_count"].max(),
        "q75_clinical_report_ids_per_record": df["clinical_report_count"].quantile(0.75),
        "q90_clinical_report_ids_per_record": df["clinical_report_count"].quantile(0.90),
        "q95_clinical_report_ids_per_record": df["clinical_report_count"].quantile(0.95),
        "q99_clinical_report_ids_per_record": df["clinical_report_count"].quantile(0.99),
    }

    out = pd.DataFrame(
        [{"metric": k, "value": round(float(v), 4)} for k, v in stats.items()]
    )

    print(out.to_string(index=False))
    return out


def top_pairs_by_phase_and_reports(df: pd.DataFrame, top_n: int = 100) -> pd.DataFrame:
    print_section("9. HIGH-CONFIDENCE-LIKE PAIRS: PHASE 4 / APPROVAL WITH MANY REPORTS")

    sub = df[df["is_approved_or_phase4"]].copy()

    out = (
        sub.groupby(["drug_id", "drug_name", "disease_id", "disease_name"])
        .agg(
            max_phase=("max_phase_num", "max"),
            clinical_phase_values=("clinical_phase_raw_clean", lambda x: "; ".join(sorted(set(x.astype(str)))[:20])),
            total_clinical_report_ids=("clinical_report_count", "sum"),
            max_clinical_report_ids_single_record=("clinical_report_count", "max"),
            total_records=("drug_id", "size"),
        )
        .reset_index()
        .sort_values(["total_clinical_report_ids", "total_records"], ascending=False)
        .head(top_n)
    )

    print(out.head(30).to_string(index=False))
    return out


def disease_keyword_landscape(df: pd.DataFrame) -> pd.DataFrame:
    print_section("10. BROAD DISEASE KEYWORD LANDSCAPE")

    keyword_patterns = {
        "cancer_or_tumour": r"cancer|tumou?r|neoplasm|carcinoma|lymphoma|leukemia|melanoma|glioblastoma|sarcoma|myeloma",
        "infection": r"infection|infectious|hiv|hepatitis|covid|tuberculosis|pneumonia|influenza|malaria|viral|bacterial|sepsis",
        "cardiovascular": r"heart|cardiac|cardio|hypertension|coronary|stroke|atrial|vascular|myocardial|infarction",
        "neurological": r"migraine|alzheimer|parkinson|epilepsy|seizure|sclerosis|neuropath|dementia|neuro",
        "psychiatric": r"depression|anxiety|schizophrenia|bipolar|autism|adhd|ptsd|psychiatric",
        "metabolic": r"diabetes|obesity|metabolic|hyperlipidemia|cholesterol|insulin|dyslipidemia",
        "autoimmune_inflammatory": r"arthritis|lupus|psoriasis|crohn|colitis|inflammatory|autoimmune",
        "respiratory": r"asthma|copd|pulmonary|respiratory|lung|pneumonia",
        "renal": r"kidney|renal|nephro",
        "pain": r"pain|analgesia|analgesic",
    }

    rows = []

    for category, pattern in keyword_patterns.items():
        mask = df["disease_name_norm"].str.contains(pattern, regex=True, na=False)
        sub = df[mask]

        rows.append(
            {
                "category": category,
                "records": len(sub),
                "unique_drugs": sub["drug_id_norm"].nunique(),
                "unique_diseases": sub["disease_id_norm"].nunique(),
                "unique_pairs": sub[["drug_id_norm", "disease_id_norm"]].drop_duplicates().shape[0],
                "approved_or_phase4_records": int(sub["is_approved_or_phase4"].sum()),
                "phase3_or_higher_records": int(sub["is_phase3_or_higher"].sum()),
                "records_with_clinical_reports": int(sub["has_clinical_report_ids"].sum()),
                "total_clinical_report_ids": int(sub["clinical_report_count"].sum()),
            }
        )

    out = pd.DataFrame(rows).sort_values("records", ascending=False)
    print(out.to_string(index=False))
    return out


def id_prefix_distribution(df: pd.DataFrame) -> pd.DataFrame:
    print_section("11. DISEASE ID PREFIX DISTRIBUTION")

    tmp = df.copy()
    tmp["disease_id_prefix"] = tmp["disease_id"].map(
        lambda x: clean_text(x).split("_")[0] if "_" in clean_text(x) else clean_text(x).split(":")[0]
    )
    tmp["disease_id_prefix"] = tmp["disease_id_prefix"].replace("", "Missing")

    out = (
        tmp.groupby("disease_id_prefix")
        .agg(
            total_records=("disease_id", "size"),
            unique_drugs=("drug_id_norm", "nunique"),
            unique_diseases=("disease_id_norm", "nunique"),
            approved_or_phase4_records=("is_approved_or_phase4", "sum"),
        )
        .reset_index()
        .sort_values("total_records", ascending=False)
    )

    print(out.to_string(index=False))
    return out


def quality_flags(df: pd.DataFrame) -> pd.DataFrame:
    print_section("12. QUALITY AND INTERPRETATION FLAGS")

    flags = {
        "missing_drug_id": df["drug_id"].eq(""),
        "missing_drug_name": df["drug_name"].eq(""),
        "drug_name_same_as_id": df["drug_name"].eq(df["drug_id"]),
        "missing_disease_id": df["disease_id"].eq(""),
        "missing_disease_name": df["disease_name"].eq(""),
        "disease_name_same_as_id": df["disease_name"].eq(df["disease_id"]),
        "missing_clinical_phase_raw": df["clinical_phase_raw_clean"].eq("Missing"),
        "max_phase_zero_or_unknown": df["max_phase_num"].eq(0),
        "missing_clinical_report_ids": ~df["has_clinical_report_ids"],
        "non_chembl_drug_id": ~df["drug_id_norm"].str.startswith("chembl", na=False),
        "non_efo_mondo_hp_disease_id": ~(
            df["disease_id_norm"].str.startswith("efo", na=False)
            | df["disease_id_norm"].str.startswith("mondo", na=False)
            | df["disease_id_norm"].str.startswith("hp", na=False)
        ),
    }

    rows = []

    for flag_name, mask in flags.items():
        sub = df[mask]
        rows.append(
            {
                "flag": flag_name,
                "records": len(sub),
                "percent_records": round(len(sub) / len(df) * 100, 3) if len(df) else 0,
                "unique_drugs": sub["drug_id_norm"].nunique(),
                "unique_diseases": sub["disease_id_norm"].nunique(),
                "unique_pairs": sub[["drug_id_norm", "disease_id_norm"]].drop_duplicates().shape[0],
            }
        )

    out = pd.DataFrame(rows).sort_values("records", ascending=False)
    print(out.to_string(index=False))
    return out


# =============================================================================
# Figures and report
# =============================================================================

def make_figures(
    outdir: Path,
    top_drugs: pd.DataFrame,
    top_diseases: pd.DataFrame,
    phase_dist: pd.DataFrame,
    keyword_landscape: pd.DataFrame,
    id_prefixes: pd.DataFrame,
) -> None:
    print_section("13. GENERATING FIGURES")

    figdir = outdir / "figures"
    figdir.mkdir(parents=True, exist_ok=True)

    barplot(
        top_drugs,
        x_col="drug_name",
        y_col="unique_disease_ids",
        title="Top Open Targets drugs by number of distinct diseases",
        out_path=figdir / "top_drugs_by_unique_diseases.png",
    )

    barplot(
        top_diseases,
        x_col="disease_name",
        y_col="unique_drug_ids",
        title="Top Open Targets diseases by number of distinct drugs",
        out_path=figdir / "top_diseases_by_unique_drugs.png",
    )

    phase_plot = phase_dist.copy()
    phase_plot["phase_label"] = (
        phase_plot["max_phase_num"].astype(str) + " / " + phase_plot["clinical_phase_raw_clean"].astype(str)
    )

    barplot(
        phase_plot,
        x_col="phase_label",
        y_col="total_records",
        title="Open Targets clinical phase distribution",
        out_path=figdir / "clinical_phase_distribution.png",
        top_n=30,
    )

    barplot(
        keyword_landscape,
        x_col="category",
        y_col="records",
        title="Broad Open Targets disease keyword landscape",
        out_path=figdir / "disease_keyword_landscape.png",
    )

    barplot(
        id_prefixes,
        x_col="disease_id_prefix",
        y_col="total_records",
        title="Disease ID prefix distribution",
        out_path=figdir / "disease_id_prefix_distribution.png",
    )


def write_markdown_report(
    outdir: Path,
    overview: pd.DataFrame,
    top_drugs: pd.DataFrame,
    top_diseases: pd.DataFrame,
    phase_dist: pd.DataFrame,
    report_support: pd.DataFrame,
    flags: pd.DataFrame,
) -> None:
    report_path = outdir / "OpenTargets_Analysis_Report.md"

    with open(report_path, "w", encoding="utf-8") as f:
        f.write("# Open Targets Drug-Disease Analysis Report\n\n")
        f.write(
            "This report summarises exploratory analysis of Open Targets "
            "`clinical_indication` drug-disease records.\n\n"
        )
        f.write(
            "Important interpretation: Open Targets clinical indication records "
            "represent clinical indication or clinical precedence-style "
            "drug-disease relationships. The maximum clinical stage and clinical "
            "report identifiers should be retained as provenance fields.\n\n"
        )

        f.write("## Overview statistics\n\n")
        f.write(overview.to_markdown(index=False))
        f.write("\n\n")

        f.write("## Top drugs by number of distinct diseases\n\n")
        f.write(top_drugs.head(20).to_markdown(index=False))
        f.write("\n\n")

        f.write("## Top diseases by number of distinct drugs\n\n")
        f.write(top_diseases.head(20).to_markdown(index=False))
        f.write("\n\n")

        f.write("## Clinical phase distribution\n\n")
        f.write(phase_dist.to_markdown(index=False))
        f.write("\n\n")

        f.write("## Clinical report support summary\n\n")
        f.write(report_support.to_markdown(index=False))
        f.write("\n\n")

        f.write("## Quality and interpretation flags\n\n")
        f.write(flags.to_markdown(index=False))
        f.write("\n\n")

    print(f"[REPORT] {report_path}")


# =============================================================================
# Main
# =============================================================================

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Exploratory analysis of Open Targets drug-disease output."
    )

    parser.add_argument(
        "--input",
        default="./data/processed/opentargets_drug_disease.csv",
        help="Path to final Open Targets processed CSV file.",
    )

    parser.add_argument(
        "--outdir",
        default="./analysis/opentargets_analysis",
        help="Output directory for Open Targets analysis tables and figures.",
    )

    parser.add_argument(
        "--top_n",
        type=int,
        default=100,
        help="Number of rows to keep in top summary tables.",
    )

    args = parser.parse_args()

    input_path = Path(args.input).resolve()
    outdir = Path(args.outdir).resolve()
    tables_dir = outdir / "tables"

    outdir.mkdir(parents=True, exist_ok=True)
    tables_dir.mkdir(parents=True, exist_ok=True)

    if not input_path.exists():
        raise FileNotFoundError(f"Input file not found: {input_path}")

    print_section("OPEN TARGETS ANALYSIS")
    print(f"[STARTED] {now_utc()}")

    df = pd.read_csv(input_path, dtype=str, low_memory=False)

    metadata = {
        "analysis_name": "Open Targets Drug-Disease Analysis",
        "started_at_utc": now_utc(),
        "interpretation": (
            "Open Targets clinical_indication records represent drug-disease "
            "clinical indication relationships. Clinical phase and report "
            "provenance should be retained."
        ),
    }

    metadata["input_validation"] = validate_input(df, input_path)

    df = add_normalised_columns(df)

    overview = overview_stats(df)
    top_drugs = top_drugs_by_unique_diseases(df, top_n=args.top_n)
    top_diseases = top_diseases_by_unique_drugs(df, top_n=args.top_n)
    supported_pairs = most_supported_drug_disease_pairs(df, top_n=args.top_n)
    phase_dist = phase_distribution(df)
    report_support = clinical_report_support_summary(df)
    high_confidence_pairs = top_pairs_by_phase_and_reports(df, top_n=args.top_n)
    keyword_landscape = disease_keyword_landscape(df)
    id_prefixes = id_prefix_distribution(df)
    flags = quality_flags(df)

    save_table(overview, tables_dir / "overview_statistics.csv")
    save_table(top_drugs, tables_dir / "top_drugs_by_unique_diseases.csv")
    save_table(top_diseases, tables_dir / "top_diseases_by_unique_drugs.csv")
    save_table(supported_pairs, tables_dir / "most_supported_drug_disease_pairs.csv")
    save_table(phase_dist, tables_dir / "clinical_phase_distribution.csv")
    save_table(report_support, tables_dir / "clinical_report_support_summary.csv")
    save_table(high_confidence_pairs, tables_dir / "phase4_pairs_with_many_reports.csv")
    save_table(keyword_landscape, tables_dir / "broad_disease_keyword_landscape.csv")
    save_table(id_prefixes, tables_dir / "disease_id_prefix_distribution.csv")
    save_table(flags, tables_dir / "quality_interpretation_flags.csv")

    make_figures(
        outdir=outdir,
        top_drugs=top_drugs,
        top_diseases=top_diseases,
        phase_dist=phase_dist,
        keyword_landscape=keyword_landscape,
        id_prefixes=id_prefixes,
    )

    write_markdown_report(
        outdir=outdir,
        overview=overview,
        top_drugs=top_drugs,
        top_diseases=top_diseases,
        phase_dist=phase_dist,
        report_support=report_support,
        flags=flags,
    )

    metadata["finished_at_utc"] = now_utc()
    metadata["status"] = "success"

    metadata["outputs"] = {
        "outdir": str(outdir),
        "tables_dir": str(tables_dir),
        "report": str(outdir / "OpenTargets_Analysis_Report.md"),
        "figures_dir": str(outdir / "figures"),
    }

    metadata_path = outdir / "opentargets_analysis_metadata.json"

    with open(metadata_path, "w", encoding="utf-8") as f:
        json.dump(metadata, f, indent=2, default=str)

    print_section("DONE")
    print(f"[OUTPUT DIRECTORY] {outdir}")
    print(f"[TABLES] {tables_dir}")
    print(f"[REPORT] {outdir / 'OpenTargets_Analysis_Report.md'}")
    print(f"[METADATA] {metadata_path}")


if __name__ == "__main__":
    main()