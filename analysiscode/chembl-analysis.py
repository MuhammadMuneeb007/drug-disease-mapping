#!/usr/bin/env python3
"""
ChEMBL Analysis: Exploratory analysis of ChEMBL drug-disease indication output.

Input:
  data/processed/chembl_drug_disease.csv

Expected columns:
  drug_id
  drug_name
  disease_id
  disease_name
  mesh_id
  mesh_heading
  max_phase
  source
  internal_source
  evidence_type

Purpose:
  1. Validate final processed ChEMBL output.
  2. Print useful statistics.
  3. Identify drugs linked to the most diseases.
  4. Identify diseases linked to the most drugs.
  5. Summarise clinical phase distribution.
  6. Detect missing identifiers and possible interpretation issues.
  7. Export analysis-ready summary tables and figures.

Interpretation:
  ChEMBL records represent curated drug indication and clinical-phase evidence.
  They should not be interpreted as clinical trial co-occurrence records.
  They should also not be treated as proof of current approval without checking
  the underlying source and clinical phase.

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

import pandas as pd
import matplotlib.pyplot as plt


REQUIRED_COLUMNS = [
    "drug_id",
    "drug_name",
    "disease_id",
    "disease_name",
    "mesh_id",
    "mesh_heading",
    "max_phase",
    "source",
    "internal_source",
    "evidence_type",
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


def print_section(title: str) -> None:
    print("\n" + "=" * 100)
    print(title)
    print("=" * 100)


def clean_text(x: object) -> str:
    if pd.isna(x):
        return ""
    return str(x).strip()


def normalise_text(x: object) -> str:
    x = clean_text(x).lower()
    x = re.sub(r"\s+", " ", x)
    return x.strip()


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

    missing_df = (
        df.isna()
        .sum()
        .reset_index()
        .rename(columns={"index": "column", 0: "missing_count"})
    )
    missing_df["missing_percent"] = (missing_df["missing_count"] / len(df) * 100).round(3)

    print("\n[MISSING VALUES]")
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

    for col in ["drug_id", "drug_name", "disease_id", "disease_name", "mesh_id", "mesh_heading"]:
        df[col] = df[col].fillna("").astype(str).str.strip()

    df["drug_id_norm"] = df["drug_id"].map(normalise_text)
    df["drug_name_norm"] = df["drug_name"].map(normalise_text)
    df["disease_id_norm"] = df["disease_id"].map(normalise_text)
    df["disease_name_norm"] = df["disease_name"].map(normalise_text)
    df["mesh_id_norm"] = df["mesh_id"].map(normalise_text)
    df["mesh_heading_norm"] = df["mesh_heading"].map(normalise_text)

    df["max_phase_num"] = pd.to_numeric(df["max_phase"], errors="coerce").fillna(0).astype(int)
    df["max_phase_label"] = df["max_phase_num"].map(
        {
            0: "Phase 0 / unknown",
            1: "Phase 1",
            2: "Phase 2",
            3: "Phase 3",
            4: "Phase 4 / approved",
        }
    ).fillna("Other / unmapped")

    print("[HEAD AFTER NORMALISATION]")
    show_cols = [
        "drug_id",
        "drug_name",
        "drug_name_norm",
        "disease_id",
        "disease_name",
        "disease_name_norm",
        "mesh_id",
        "mesh_heading",
        "max_phase",
        "max_phase_label",
    ]
    print(df[show_cols].head(10).to_string(index=False))

    return df


def overview_stats(df: pd.DataFrame) -> pd.DataFrame:
    print_section("3. OVERVIEW STATISTICS")

    stats = {
        "total_rows": len(df),
        "unique_drug_ids": df["drug_id_norm"].replace("", pd.NA).nunique(dropna=True),
        "unique_drug_names": df["drug_name_norm"].replace("", pd.NA).nunique(dropna=True),
        "unique_disease_ids": df["disease_id_norm"].replace("", pd.NA).nunique(dropna=True),
        "unique_disease_names": df["disease_name_norm"].replace("", pd.NA).nunique(dropna=True),
        "unique_mesh_ids": df["mesh_id_norm"].replace("", pd.NA).nunique(dropna=True),
        "unique_mesh_headings": df["mesh_heading_norm"].replace("", pd.NA).nunique(dropna=True),
        "unique_drug_disease_pairs_by_id": df[["drug_id_norm", "disease_id_norm"]]
        .replace("", pd.NA)
        .dropna()
        .drop_duplicates()
        .shape[0],
        "unique_drug_disease_pairs_by_name": df[["drug_name_norm", "disease_name_norm"]]
        .replace("", pd.NA)
        .dropna()
        .drop_duplicates()
        .shape[0],
        "exact_duplicate_rows": int(df.duplicated().sum()),
    }

    out = pd.DataFrame(list(stats.items()), columns=["metric", "value"])
    print(out.to_string(index=False))
    return out


def top_drugs_by_unique_diseases(df: pd.DataFrame, top_n: int = 50) -> pd.DataFrame:
    print_section("4. DRUGS LINKED TO THE LARGEST NUMBER OF DISTINCT DISEASES")

    valid = df[(df["drug_name_norm"] != "") & (df["disease_name_norm"] != "")].copy()

    out = (
        valid.groupby(["drug_id_norm", "drug_name_norm"])
        .agg(
            unique_disease_names=("disease_name_norm", "nunique"),
            unique_disease_ids=("disease_id_norm", "nunique"),
            unique_mesh_ids=("mesh_id_norm", "nunique"),
            total_records=("drug_id_norm", "size"),
            max_phase_highest=("max_phase_num", "max"),
        )
        .reset_index()
        .sort_values(
            ["unique_disease_names", "unique_disease_ids", "total_records"],
            ascending=False,
        )
        .head(top_n)
    )

    print(out.head(25).to_string(index=False))
    return out


def top_diseases_by_unique_drugs(df: pd.DataFrame, top_n: int = 50) -> pd.DataFrame:
    print_section("5. DISEASES LINKED TO THE LARGEST NUMBER OF DISTINCT DRUGS")

    valid = df[(df["drug_name_norm"] != "") & (df["disease_name_norm"] != "")].copy()

    out = (
        valid.groupby(["disease_id_norm", "disease_name_norm"])
        .agg(
            unique_drug_names=("drug_name_norm", "nunique"),
            unique_drug_ids=("drug_id_norm", "nunique"),
            total_records=("disease_name_norm", "size"),
            max_phase_highest=("max_phase_num", "max"),
        )
        .reset_index()
        .sort_values(
            ["unique_drug_names", "unique_drug_ids", "total_records"],
            ascending=False,
        )
        .head(top_n)
    )

    print(out.head(25).to_string(index=False))
    return out


def most_common_drugs(df: pd.DataFrame, top_n: int = 50) -> pd.DataFrame:
    print_section("6. MOST COMMON DRUGS IN ChEMBL INDICATION RECORDS")

    valid = df[df["drug_name_norm"] != ""].copy()

    out = (
        valid.groupby(["drug_id_norm", "drug_name_norm"])
        .agg(
            total_records=("drug_name_norm", "size"),
            unique_disease_names=("disease_name_norm", "nunique"),
            max_phase_highest=("max_phase_num", "max"),
        )
        .reset_index()
        .sort_values(["total_records", "unique_disease_names"], ascending=False)
        .head(top_n)
    )

    print(out.head(25).to_string(index=False))
    return out


def most_common_diseases(df: pd.DataFrame, top_n: int = 50) -> pd.DataFrame:
    print_section("7. MOST COMMON DISEASES IN ChEMBL INDICATION RECORDS")

    valid = df[df["disease_name_norm"] != ""].copy()

    out = (
        valid.groupby(["disease_id_norm", "disease_name_norm"])
        .agg(
            total_records=("disease_name_norm", "size"),
            unique_drug_names=("drug_name_norm", "nunique"),
            max_phase_highest=("max_phase_num", "max"),
        )
        .reset_index()
        .sort_values(["total_records", "unique_drug_names"], ascending=False)
        .head(top_n)
    )

    print(out.head(25).to_string(index=False))
    return out


def most_common_drug_disease_pairs(df: pd.DataFrame, top_n: int = 100) -> pd.DataFrame:
    print_section("8. MOST COMMON DRUG–DISEASE PAIRS")

    valid = df[(df["drug_name_norm"] != "") & (df["disease_name_norm"] != "")].copy()

    out = (
        valid.groupby(["drug_id_norm", "drug_name_norm", "disease_id_norm", "disease_name_norm"])
        .agg(
            total_records=("drug_name_norm", "size"),
            max_phase_highest=("max_phase_num", "max"),
            max_phase_labels=("max_phase_label", lambda x: "; ".join(sorted(set(x.dropna().astype(str))))),
            mesh_ids=("mesh_id_norm", lambda x: "; ".join(sorted(set([v for v in x.dropna().astype(str) if v]))[:10])),
            mesh_headings=("mesh_heading_norm", lambda x: "; ".join(sorted(set([v for v in x.dropna().astype(str) if v]))[:10])),
        )
        .reset_index()
        .sort_values(["total_records", "max_phase_highest"], ascending=False)
        .head(top_n)
    )

    print(out.head(30).to_string(index=False))
    return out


def phase_distribution(df: pd.DataFrame) -> pd.DataFrame:
    print_section("9. MAXIMUM CLINICAL PHASE DISTRIBUTION")

    out = (
        df.groupby(["max_phase_num", "max_phase_label"])
        .agg(
            total_records=("drug_id_norm", "size"),
            unique_drug_names=("drug_name_norm", "nunique"),
            unique_disease_names=("disease_name_norm", "nunique"),
            unique_drug_disease_pairs=("drug_id_norm", lambda x: None),
        )
        .reset_index()
        .drop(columns=["unique_drug_disease_pairs"])
        .sort_values("max_phase_num", ascending=True)
    )

    # Add pair counts separately to avoid awkward groupby lambda on multiple columns.
    pair_counts = (
        df.groupby(["max_phase_num", "max_phase_label"])
        .apply(lambda x: x[["drug_id_norm", "disease_name_norm"]].drop_duplicates().shape[0], include_groups=False)
        .reset_index(name="unique_drug_disease_pairs")
    )

    out = out.merge(pair_counts, on=["max_phase_num", "max_phase_label"], how="left")

    print(out.to_string(index=False))
    return out


def disease_keyword_landscape(df: pd.DataFrame) -> pd.DataFrame:
    print_section("10. BROAD DISEASE KEYWORD LANDSCAPE")

    keyword_patterns = {
        "cancer_or_tumour": r"cancer|tumou?r|neoplasm|carcinoma|lymphoma|leukemia|melanoma|glioblastoma|sarcoma|myeloma",
        "infection": r"infection|infectious|hiv|hepatitis|covid|tuberculosis|pneumonia|influenza|malaria|viral|bacterial",
        "cardiovascular": r"heart|cardiac|cardio|hypertension|coronary|stroke|atrial|vascular|myocardial",
        "neurological": r"migraine|alzheimer|parkinson|epilepsy|seizure|sclerosis|neuropath|dementia|neuro",
        "psychiatric": r"depression|anxiety|schizophrenia|bipolar|autism|adhd|ptsd|psychiatric",
        "metabolic": r"diabetes|obesity|metabolic|hyperlipidemia|cholesterol|insulin|dyslipidemia",
        "autoimmune_inflammatory": r"arthritis|lupus|psoriasis|crohn|colitis|inflammatory|autoimmune",
        "respiratory": r"asthma|copd|pulmonary|respiratory|lung",
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
                "unique_drug_names": sub["drug_name_norm"].nunique(),
                "unique_disease_names": sub["disease_name_norm"].nunique(),
                "unique_drug_disease_pairs": sub[["drug_id_norm", "disease_name_norm"]].drop_duplicates().shape[0],
                "highest_max_phase": sub["max_phase_num"].max() if len(sub) else 0,
            }
        )

    out = pd.DataFrame(rows).sort_values("records", ascending=False)
    print(out.to_string(index=False))
    return out


def missingness_and_quality_flags(df: pd.DataFrame) -> pd.DataFrame:
    print_section("11. MISSINGNESS AND INTERPRETATION FLAGS")

    flags = {
        "missing_drug_name": df["drug_name_norm"].eq(""),
        "missing_disease_name": df["disease_name_norm"].eq(""),
        "missing_disease_id": df["disease_id_norm"].eq(""),
        "missing_mesh_id": df["mesh_id_norm"].eq(""),
        "missing_mesh_heading": df["mesh_heading_norm"].eq(""),
        "max_phase_zero_or_unknown": df["max_phase_num"].eq(0),
        "phase_4_or_approved": df["max_phase_num"].ge(4),
        "drug_name_contains_combination_marker": df["drug_name"].astype(str).str.contains(r"\+|/| and | with ", case=False, regex=True, na=False),
    }

    rows = []
    for flag_name, mask in flags.items():
        sub = df[mask]
        rows.append(
            {
                "flag": flag_name,
                "records": len(sub),
                "percent_records": round(len(sub) / len(df) * 100, 3) if len(df) else 0,
                "unique_drug_names": sub["drug_name_norm"].nunique(),
                "unique_disease_names": sub["disease_name_norm"].nunique(),
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
    keyword_landscape: pd.DataFrame,
) -> None:
    print_section("12. GENERATING FIGURES")

    figdir = outdir / "figures"
    figdir.mkdir(parents=True, exist_ok=True)

    barplot(
        top_drugs_diseases,
        x_col="drug_name_norm",
        y_col="unique_disease_names",
        title="Top ChEMBL drugs by number of distinct disease indications",
        out_path=figdir / "top_drugs_by_unique_diseases.png",
    )

    barplot(
        top_diseases_drugs,
        x_col="disease_name_norm",
        y_col="unique_drug_names",
        title="Top ChEMBL diseases by number of distinct drugs",
        out_path=figdir / "top_diseases_by_unique_drugs.png",
    )

    barplot(
        common_drugs,
        x_col="drug_name_norm",
        y_col="total_records",
        title="Most common ChEMBL drug indication records",
        out_path=figdir / "most_common_drugs.png",
    )

    barplot(
        common_diseases,
        x_col="disease_name_norm",
        y_col="total_records",
        title="Most common ChEMBL disease indication records",
        out_path=figdir / "most_common_diseases.png",
    )

    barplot(
        phase,
        x_col="max_phase_label",
        y_col="total_records",
        title="ChEMBL indication records by maximum clinical phase",
        out_path=figdir / "max_phase_distribution.png",
    )

    barplot(
        keyword_landscape,
        x_col="category",
        y_col="records",
        title="Broad ChEMBL disease keyword landscape",
        out_path=figdir / "disease_keyword_landscape.png",
    )


def write_markdown_report(
    outdir: Path,
    overview: pd.DataFrame,
    top_drugs_diseases: pd.DataFrame,
    top_diseases_drugs: pd.DataFrame,
    common_pairs: pd.DataFrame,
    phase: pd.DataFrame,
    flags: pd.DataFrame,
) -> None:
    report_path = outdir / "ChEMBL_Analysis_Report.md"

    with open(report_path, "w", encoding="utf-8") as f:
        f.write("# ChEMBL Analysis Report\n\n")
        f.write(
            "This report summarises exploratory analysis of ChEMBL-derived "
            "drug-disease indication records.\n\n"
        )
        f.write(
            "Important interpretation: ChEMBL records represent curated drug "
            "indication and clinical-phase evidence. They should not be treated "
            "as identical to AACT/ClinicalTrials.gov intervention-condition "
            "co-occurrence records.\n\n"
        )

        f.write("## Overview statistics\n\n")
        f.write(overview.to_markdown(index=False))
        f.write("\n\n")

        f.write("## Top drugs by number of distinct diseases\n\n")
        f.write(top_drugs_diseases.head(20).to_markdown(index=False))
        f.write("\n\n")

        f.write("## Top diseases by number of distinct drugs\n\n")
        f.write(top_diseases_drugs.head(20).to_markdown(index=False))
        f.write("\n\n")

        f.write("## Most common drug-disease pairs\n\n")
        f.write(common_pairs.head(30).to_markdown(index=False))
        f.write("\n\n")

        f.write("## Maximum clinical phase distribution\n\n")
        f.write(phase.to_markdown(index=False))
        f.write("\n\n")

        f.write("## Missingness and interpretation flags\n\n")
        f.write(flags.to_markdown(index=False))
        f.write("\n\n")

    print(f"[REPORT] {report_path}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Exploratory analysis of ChEMBL drug-disease output.")
    parser.add_argument(
        "--input",
        default="./data/processed/chembl_drug_disease.csv",
        help="Path to final ChEMBL processed CSV file.",
    )
    parser.add_argument(
        "--outdir",
        default="./analysis/chembl_analysis",
        help="Output directory for ChEMBL analysis tables and figures.",
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

    print_section("ChEMBL ANALYSIS")
    print(f"[STARTED] {now_utc()}")

    df = pd.read_csv(input_path, dtype=str, low_memory=False)

    metadata = {
        "analysis_name": "ChEMBL Analysis",
        "started_at_utc": now_utc(),
    }

    metadata["input_validation"] = validate_input(df, input_path)

    df = add_normalised_columns(df)

    overview = overview_stats(df)
    top_drugs_diseases = top_drugs_by_unique_diseases(df, top_n=args.top_n)
    top_diseases_drugs = top_diseases_by_unique_drugs(df, top_n=args.top_n)
    common_drugs = most_common_drugs(df, top_n=args.top_n)
    common_diseases = most_common_diseases(df, top_n=args.top_n)
    common_pairs = most_common_drug_disease_pairs(df, top_n=args.top_n)
    phase = phase_distribution(df)
    keyword_landscape = disease_keyword_landscape(df)
    flags = missingness_and_quality_flags(df)

    save_table(overview, tables_dir / "overview_statistics.csv")
    save_table(top_drugs_diseases, tables_dir / "top_drugs_by_unique_diseases.csv")
    save_table(top_diseases_drugs, tables_dir / "top_diseases_by_unique_drugs.csv")
    save_table(common_drugs, tables_dir / "most_common_drugs.csv")
    save_table(common_diseases, tables_dir / "most_common_diseases.csv")
    save_table(common_pairs, tables_dir / "most_common_drug_disease_pairs.csv")
    save_table(phase, tables_dir / "max_phase_distribution.csv")
    save_table(keyword_landscape, tables_dir / "broad_disease_keyword_landscape.csv")
    save_table(flags, tables_dir / "missingness_and_quality_flags.csv")

    make_figures(
        outdir=outdir,
        top_drugs_diseases=top_drugs_diseases,
        top_diseases_drugs=top_diseases_drugs,
        common_drugs=common_drugs,
        common_diseases=common_diseases,
        phase=phase,
        keyword_landscape=keyword_landscape,
    )

    write_markdown_report(
        outdir=outdir,
        overview=overview,
        top_drugs_diseases=top_drugs_diseases,
        top_diseases_drugs=top_diseases_drugs,
        common_pairs=common_pairs,
        phase=phase,
        flags=flags,
    )

    metadata["finished_at_utc"] = now_utc()
    metadata["status"] = "success"

    metadata_path = outdir / "chembl_analysis_metadata.json"
    with open(metadata_path, "w", encoding="utf-8") as f:
        json.dump(metadata, f, indent=2)

    print_section("DONE")
    print(f"[OUTPUT DIRECTORY] {outdir}")
    print(f"[METADATA] {metadata_path}")


if __name__ == "__main__":
    main()