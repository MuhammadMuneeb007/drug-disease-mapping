#!/usr/bin/env python3
"""
AACT Analysis: Exploratory analysis of AACT drug intervention-condition output.

Input:
  data/processed/aact_drug_disease.csv

Expected columns:
  nct_id
  drug_name
  disease_name
  intervention_type
  study_phase
  overall_status
  source
  internal_source
  evidence_type

Purpose:
  1. Validate final processed AACT output.
  2. Print useful statistics.
  3. Identify most common drugs, diseases, trial statuses, and phases.
  4. Export analysis-ready summary tables.
  5. Generate manuscript-friendly exploratory summaries.

Interpretation:
  AACT/ClinicalTrials.gov records represent trial intervention-condition
  co-occurrence. They do not imply drug approval, efficacy, or positive trial outcome.

Requirements:
  pip install pandas matplotlib
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt


REQUIRED_COLUMNS = [
    "nct_id",
    "drug_name",
    "disease_name",
    "intervention_type",
    "study_phase",
    "overall_status",
]


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


def clean_text(x: object) -> str:
    if pd.isna(x):
        return ""
    return str(x).strip()


def normalise_text(x: object) -> str:
    x = clean_text(x).lower()
    x = re.sub(r"\s+", " ", x)
    return x.strip()


def print_section(title: str) -> None:
    print("\n" + "=" * 100)
    print(title)
    print("=" * 100)


def save_table(df: pd.DataFrame, out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_path, index=False)
    print(f"[SAVED] {out_path}")


def barplot(df: pd.DataFrame, x_col: str, y_col: str, title: str, out_path: Path, top_n: int = 20) -> None:
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


def validate_input(df: pd.DataFrame, input_path: Path) -> dict:
    print_section("1. INPUT VALIDATION")

    print(f"[INPUT FILE] {input_path}")
    print(f"[SIZE] {human_size(input_path.stat().st_size)}")
    print(f"[SHA256] {sha256_file(input_path)}")
    print(f"[ROWS, COLUMNS] {df.shape}")

    print("\n[COLUMNS]")
    print(list(df.columns))

    missing = [c for c in REQUIRED_COLUMNS if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    print("\n[HEAD]")
    print(df.head(10).to_string(index=False))

    print("\n[MISSING VALUES]")
    missing_df = (
        df.isna()
        .sum()
        .reset_index()
        .rename(columns={"index": "column", 0: "missing_count"})
    )
    missing_df["missing_percent"] = (missing_df["missing_count"] / len(df) * 100).round(3)
    print(missing_df.to_string(index=False))

    return {
        "input_file": str(input_path),
        "input_size_bytes": input_path.stat().st_size,
        "input_sha256": sha256_file(input_path),
        "n_rows": int(df.shape[0]),
        "n_columns": int(df.shape[1]),
        "columns": list(df.columns),
        "missing_values": missing_df.to_dict(orient="records"),
    }


def add_normalised_columns(df: pd.DataFrame) -> pd.DataFrame:
    print_section("2. CLEANING / NORMALISED ANALYSIS COLUMNS")

    df = df.copy()

    df["drug_name_clean"] = df["drug_name"].map(clean_text)
    df["disease_name_clean"] = df["disease_name"].map(clean_text)
    df["drug_name_norm"] = df["drug_name"].map(normalise_text)
    df["disease_name_norm"] = df["disease_name"].map(normalise_text)

    df["study_phase_clean"] = df["study_phase"].fillna("Missing").astype(str).str.strip()
    df.loc[df["study_phase_clean"].isin(["", "nan", "NaN", "None"]), "study_phase_clean"] = "Missing"

    df["overall_status_clean"] = df["overall_status"].fillna("Missing").astype(str).str.strip()
    df.loc[df["overall_status_clean"].isin(["", "nan", "NaN", "None"]), "overall_status_clean"] = "Missing"

    print("[HEAD AFTER ADDING CLEAN/NORMALISED COLUMNS]")
    cols = [
        "nct_id",
        "drug_name",
        "drug_name_norm",
        "disease_name",
        "disease_name_norm",
        "study_phase_clean",
        "overall_status_clean",
    ]
    print(df[cols].head(10).to_string(index=False))

    return df


def overview_stats(df: pd.DataFrame) -> pd.DataFrame:
    print_section("3. OVERVIEW STATISTICS")

    stats = {
        "total_rows": len(df),
        "unique_trials": df["nct_id"].nunique(),
        "unique_drug_names_raw": df["drug_name"].nunique(dropna=True),
        "unique_disease_names_raw": df["disease_name"].nunique(dropna=True),
        "unique_drug_names_normalised": df["drug_name_norm"].nunique(),
        "unique_disease_names_normalised": df["disease_name_norm"].nunique(),
        "unique_drug_disease_pairs_normalised": df[["drug_name_norm", "disease_name_norm"]].drop_duplicates().shape[0],
        "unique_trial_drug_pairs": df[["nct_id", "drug_name_norm"]].drop_duplicates().shape[0],
        "unique_trial_disease_pairs": df[["nct_id", "disease_name_norm"]].drop_duplicates().shape[0],
    }

    stats_df = pd.DataFrame(list(stats.items()), columns=["metric", "value"])
    print(stats_df.to_string(index=False))
    return stats_df


def top_drugs_by_number_of_diseases(df: pd.DataFrame, top_n: int = 50) -> pd.DataFrame:
    print_section("4. DRUGS ASSOCIATED WITH THE LARGEST NUMBER OF DISTINCT CONDITIONS")

    out = (
        df.groupby("drug_name_norm")
        .agg(
            unique_diseases=("disease_name_norm", "nunique"),
            unique_trials=("nct_id", "nunique"),
            total_records=("nct_id", "size"),
        )
        .reset_index()
        .sort_values(["unique_diseases", "unique_trials", "total_records"], ascending=False)
        .head(top_n)
    )

    print(out.head(25).to_string(index=False))
    return out


def top_diseases_by_number_of_drugs(df: pd.DataFrame, top_n: int = 50) -> pd.DataFrame:
    print_section("5. CONDITIONS STUDIED WITH THE LARGEST NUMBER OF DISTINCT DRUGS")

    out = (
        df.groupby("disease_name_norm")
        .agg(
            unique_drugs=("drug_name_norm", "nunique"),
            unique_trials=("nct_id", "nunique"),
            total_records=("nct_id", "size"),
        )
        .reset_index()
        .sort_values(["unique_drugs", "unique_trials", "total_records"], ascending=False)
        .head(top_n)
    )

    print(out.head(25).to_string(index=False))
    return out


def most_common_drugs(df: pd.DataFrame, top_n: int = 50) -> pd.DataFrame:
    print_section("6. MOST COMMONLY OCCURRING DRUG NAMES")

    out = (
        df.groupby("drug_name_norm")
        .agg(
            total_records=("nct_id", "size"),
            unique_trials=("nct_id", "nunique"),
            unique_diseases=("disease_name_norm", "nunique"),
        )
        .reset_index()
        .sort_values(["total_records", "unique_trials"], ascending=False)
        .head(top_n)
    )

    print(out.head(25).to_string(index=False))
    return out


def most_common_diseases(df: pd.DataFrame, top_n: int = 50) -> pd.DataFrame:
    print_section("7. MOST COMMONLY OCCURRING CONDITIONS")

    out = (
        df.groupby("disease_name_norm")
        .agg(
            total_records=("nct_id", "size"),
            unique_trials=("nct_id", "nunique"),
            unique_drugs=("drug_name_norm", "nunique"),
        )
        .reset_index()
        .sort_values(["total_records", "unique_trials"], ascending=False)
        .head(top_n)
    )

    print(out.head(25).to_string(index=False))
    return out


def most_common_drug_disease_pairs(df: pd.DataFrame, top_n: int = 100) -> pd.DataFrame:
    print_section("8. MOST COMMON DRUG–CONDITION PAIRS")

    out = (
        df.groupby(["drug_name_norm", "disease_name_norm"])
        .agg(
            total_records=("nct_id", "size"),
            unique_trials=("nct_id", "nunique"),
            phases=("study_phase_clean", lambda x: "; ".join(sorted(set(x.dropna().astype(str)))[:10])),
            statuses=("overall_status_clean", lambda x: "; ".join(sorted(set(x.dropna().astype(str)))[:10])),
        )
        .reset_index()
        .sort_values(["unique_trials", "total_records"], ascending=False)
        .head(top_n)
    )

    print(out.head(30).to_string(index=False))
    return out


def phase_status_analysis(df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    print_section("9. TRIAL PHASE AND STATUS DISTRIBUTION")

    phase = (
        df.groupby("study_phase_clean")
        .agg(
            total_records=("nct_id", "size"),
            unique_trials=("nct_id", "nunique"),
            unique_drugs=("drug_name_norm", "nunique"),
            unique_diseases=("disease_name_norm", "nunique"),
        )
        .reset_index()
        .sort_values("unique_trials", ascending=False)
    )

    status = (
        df.groupby("overall_status_clean")
        .agg(
            total_records=("nct_id", "size"),
            unique_trials=("nct_id", "nunique"),
            unique_drugs=("drug_name_norm", "nunique"),
            unique_diseases=("disease_name_norm", "nunique"),
        )
        .reset_index()
        .sort_values("unique_trials", ascending=False)
    )

    phase_status = (
        df.groupby(["study_phase_clean", "overall_status_clean"])
        .agg(
            total_records=("nct_id", "size"),
            unique_trials=("nct_id", "nunique"),
        )
        .reset_index()
        .sort_values(["unique_trials", "total_records"], ascending=False)
    )

    print("\n[PHASE DISTRIBUTION]")
    print(phase.to_string(index=False))

    print("\n[STATUS DISTRIBUTION]")
    print(status.to_string(index=False))

    print("\n[PHASE x STATUS TOP 30]")
    print(phase_status.head(30).to_string(index=False))

    return phase, status, phase_status


def broad_keyword_landscape(df: pd.DataFrame) -> pd.DataFrame:
    print_section("10. BROAD CONDITION KEYWORD LANDSCAPE")

    keyword_patterns = {
        "cancer_or_tumour": r"cancer|tumou?r|neoplasm|carcinoma|lymphoma|leukemia|melanoma|glioblastoma|sarcoma",
        "infection": r"infection|infectious|hiv|hepatitis|covid|tuberculosis|pneumonia|influenza|malaria",
        "cardiovascular": r"heart|cardiac|cardio|hypertension|coronary|stroke|atrial|vascular|myocardial",
        "neurological": r"migraine|alzheimer|parkinson|epilepsy|seizure|sclerosis|neuropath|stroke|dementia",
        "psychiatric": r"depression|anxiety|schizophrenia|bipolar|autism|adhd|ptsd",
        "metabolic": r"diabetes|obesity|metabolic|hyperlipidemia|cholesterol|insulin",
        "autoimmune_inflammatory": r"arthritis|lupus|psoriasis|crohn|colitis|inflammatory|autoimmune",
        "respiratory": r"asthma|copd|pulmonary|respiratory|lung",
        "healthy_volunteers": r"healthy|healthy volunteer|healthy subjects",
        "bioequivalence": r"bioequivalence|bioequivalent",
    }

    rows = []
    for label, pattern in keyword_patterns.items():
        mask = df["disease_name_norm"].str.contains(pattern, regex=True, na=False)
        sub = df[mask]
        rows.append(
            {
                "category": label,
                "records": len(sub),
                "unique_trials": sub["nct_id"].nunique(),
                "unique_drugs": sub["drug_name_norm"].nunique(),
                "unique_diseases": sub["disease_name_norm"].nunique(),
            }
        )

    out = pd.DataFrame(rows).sort_values("unique_trials", ascending=False)
    print(out.to_string(index=False))
    return out


def possible_quality_flags(df: pd.DataFrame) -> pd.DataFrame:
    print_section("11. POTENTIAL DATA QUALITY / INTERPRETATION FLAGS")

    flags = {
        "drug_name_contains_placebo": df["drug_name_norm"].str.contains(r"\bplacebo\b", regex=True, na=False),
        "condition_healthy_volunteers": df["disease_name_norm"].str.contains(r"healthy|healthy volunteer|healthy subjects", regex=True, na=False),
        "condition_bioequivalence": df["disease_name_norm"].str.contains(r"bioequivalence|bioequivalent", regex=True, na=False),
        "drug_name_contains_plus_or_slash": df["drug_name_clean"].str.contains(r"\+|/", regex=True, na=False),
        "missing_phase": df["study_phase_clean"].eq("Missing"),
        "withdrawn_or_terminated": df["overall_status_clean"].isin(["WITHDRAWN", "TERMINATED"]),
    }

    rows = []
    for flag_name, mask in flags.items():
        sub = df[mask]
        rows.append(
            {
                "flag": flag_name,
                "records": len(sub),
                "percent_records": round(len(sub) / len(df) * 100, 3),
                "unique_trials": sub["nct_id"].nunique(),
                "unique_drugs": sub["drug_name_norm"].nunique(),
                "unique_diseases": sub["disease_name_norm"].nunique(),
            }
        )

    out = pd.DataFrame(rows).sort_values("records", ascending=False)
    print(out.to_string(index=False))
    return out


def make_figures(
    outdir: Path,
    top_drugs_diseases: pd.DataFrame,
    top_diseases_drugs: pd.DataFrame,
    common_drugs: pd.DataFrame,
    common_diseases: pd.DataFrame,
    phase: pd.DataFrame,
    status: pd.DataFrame,
    keyword_landscape: pd.DataFrame,
) -> None:
    print_section("12. GENERATING FIGURES")

    figdir = outdir / "figures"
    figdir.mkdir(parents=True, exist_ok=True)

    barplot(
        top_drugs_diseases,
        x_col="drug_name_norm",
        y_col="unique_diseases",
        title="Top drugs by number of distinct AACT conditions",
        out_path=figdir / "top_drugs_by_unique_diseases.png",
    )

    barplot(
        top_diseases_drugs,
        x_col="disease_name_norm",
        y_col="unique_drugs",
        title="Top conditions by number of distinct AACT drug interventions",
        out_path=figdir / "top_conditions_by_unique_drugs.png",
    )

    barplot(
        common_drugs,
        x_col="drug_name_norm",
        y_col="unique_trials",
        title="Most common drug interventions by unique trials",
        out_path=figdir / "most_common_drugs_by_trials.png",
    )

    barplot(
        common_diseases,
        x_col="disease_name_norm",
        y_col="unique_trials",
        title="Most common conditions by unique trials",
        out_path=figdir / "most_common_conditions_by_trials.png",
    )

    barplot(
        phase,
        x_col="study_phase_clean",
        y_col="unique_trials",
        title="AACT drug-condition records by trial phase",
        out_path=figdir / "phase_distribution.png",
    )

    barplot(
        status,
        x_col="overall_status_clean",
        y_col="unique_trials",
        title="AACT drug-condition records by trial status",
        out_path=figdir / "status_distribution.png",
    )

    barplot(
        keyword_landscape,
        x_col="category",
        y_col="unique_trials",
        title="Broad condition keyword landscape",
        out_path=figdir / "broad_condition_keyword_landscape.png",
    )


def write_markdown_report(
    outdir: Path,
    overview: pd.DataFrame,
    top_drugs_diseases: pd.DataFrame,
    top_diseases_drugs: pd.DataFrame,
    common_pairs: pd.DataFrame,
    quality_flags: pd.DataFrame,
) -> None:
    report_path = outdir / "AACT_Analysis_Report.md"

    with open(report_path, "w", encoding="utf-8") as f:
        f.write("# AACT Analysis Report\n\n")
        f.write("This report summarises exploratory analysis of AACT-derived drug intervention-condition records.\n\n")
        f.write("Important interpretation: AACT records represent clinical trial intervention-condition co-occurrence. They do not imply approval, efficacy, or positive trial outcome.\n\n")

        f.write("## Overview statistics\n\n")
        f.write(overview.to_markdown(index=False))
        f.write("\n\n")

        f.write("## Top drugs by number of distinct conditions\n\n")
        f.write(top_drugs_diseases.head(20).to_markdown(index=False))
        f.write("\n\n")

        f.write("## Top conditions by number of distinct drug interventions\n\n")
        f.write(top_diseases_drugs.head(20).to_markdown(index=False))
        f.write("\n\n")

        f.write("## Most common drug-condition pairs\n\n")
        f.write(common_pairs.head(30).to_markdown(index=False))
        f.write("\n\n")

        f.write("## Potential quality and interpretation flags\n\n")
        f.write(quality_flags.to_markdown(index=False))
        f.write("\n\n")

    print(f"[REPORT] {report_path}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Exploratory analysis of AACT drug-condition output.")
    parser.add_argument(
        "--input",
        default="./data/processed/aact_drug_disease.csv",
        help="Path to final AACT processed CSV file.",
    )
    parser.add_argument(
        "--outdir",
        default="./analysis/aact_analysis",
        help="Output directory for analysis tables and figures.",
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

    print_section("AACT ANALYSIS")
    print(f"[STARTED] {now_utc()}")

    df = pd.read_csv(input_path, dtype=str, low_memory=False)

    metadata = {
        "analysis_name": "AACT Analysis",
        "started_at_utc": now_utc(),
    }

    metadata["input_validation"] = validate_input(df, input_path)

    df = add_normalised_columns(df)

    overview = overview_stats(df)
    top_drugs_diseases = top_drugs_by_number_of_diseases(df, top_n=args.top_n)
    top_diseases_drugs = top_diseases_by_number_of_drugs(df, top_n=args.top_n)
    common_drugs = most_common_drugs(df, top_n=args.top_n)
    common_diseases = most_common_diseases(df, top_n=args.top_n)
    common_pairs = most_common_drug_disease_pairs(df, top_n=args.top_n)
    phase, status, phase_status = phase_status_analysis(df)
    keyword_landscape = broad_keyword_landscape(df)
    quality_flags = possible_quality_flags(df)

    save_table(overview, tables_dir / "overview_statistics.csv")
    save_table(top_drugs_diseases, tables_dir / "top_drugs_by_unique_diseases.csv")
    save_table(top_diseases_drugs, tables_dir / "top_conditions_by_unique_drugs.csv")
    save_table(common_drugs, tables_dir / "most_common_drugs.csv")
    save_table(common_diseases, tables_dir / "most_common_conditions.csv")
    save_table(common_pairs, tables_dir / "most_common_drug_condition_pairs.csv")
    save_table(phase, tables_dir / "phase_distribution.csv")
    save_table(status, tables_dir / "status_distribution.csv")
    save_table(phase_status, tables_dir / "phase_status_distribution.csv")
    save_table(keyword_landscape, tables_dir / "broad_condition_keyword_landscape.csv")
    save_table(quality_flags, tables_dir / "quality_interpretation_flags.csv")

    make_figures(
        outdir=outdir,
        top_drugs_diseases=top_drugs_diseases,
        top_diseases_drugs=top_diseases_drugs,
        common_drugs=common_drugs,
        common_diseases=common_diseases,
        phase=phase,
        status=status,
        keyword_landscape=keyword_landscape,
    )

    write_markdown_report(
        outdir=outdir,
        overview=overview,
        top_drugs_diseases=top_drugs_diseases,
        top_diseases_drugs=top_diseases_drugs,
        common_pairs=common_pairs,
        quality_flags=quality_flags,
    )

    metadata["finished_at_utc"] = now_utc()
    metadata["status"] = "success"
    metadata_path = outdir / "aact_analysis_metadata.json"
    with open(metadata_path, "w", encoding="utf-8") as f:
        json.dump(metadata, f, indent=2)

    print_section("DONE")
    print(f"[OUTPUT DIRECTORY] {outdir}")
    print(f"[METADATA] {metadata_path}")


if __name__ == "__main__":
    main()