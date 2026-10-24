#!/usr/bin/env python3
"""
CTD Analysis: Exploratory analysis of CTD chemical-disease association output.

Input:
  data/processed/ctd_aggregate_chemical_disease.csv

Expected columns:
  chemical_name
  chemical_id
  cas_rn
  disease_name
  disease_id
  direct_evidence
  inference_gene_symbol
  inference_score
  omim_ids
  pubmed_ids
  drug_name
  drug_id
  source
  internal_source
  evidence_type

Purpose:
  1. Validate final processed CTD output.
  2. Summarise chemical-disease association coverage.
  3. Identify chemicals linked to the most diseases.
  4. Identify diseases linked to the most chemicals.
  5. Analyse direct evidence and inferred evidence.
  6. Analyse inference-score distribution.
  7. Analyse PubMed support.
  8. Export tables, figures, report, and metadata.

Important interpretation:
  CTD records represent chemical-disease associations, not necessarily
  therapeutic drug indications. Aggregate CTD includes many inferred
  associations, so these should be kept distinct from ChEMBL curated
  indication records and AACT/ClinicalTrials.gov trial co-occurrence records.

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
    "chemical_name",
    "chemical_id",
    "disease_name",
    "disease_id",
    "direct_evidence",
    "inference_gene_symbol",
    "inference_score",
    "pubmed_ids",
    "drug_name",
    "drug_id",
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


def count_pipe_items(x: object) -> int:
    x = clean_text(x)
    if not x:
        return 0
    return len([v for v in x.split("|") if v.strip()])


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

    missing_required = [c for c in REQUIRED_COLUMNS if c not in df.columns]
    if missing_required:
        raise ValueError(f"Missing required columns: {missing_required}")

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

    for col in [
        "chemical_name",
        "chemical_id",
        "disease_name",
        "disease_id",
        "direct_evidence",
        "inference_gene_symbol",
        "inference_score",
        "omim_ids",
        "pubmed_ids",
        "drug_name",
        "drug_id",
    ]:
        if col in df.columns:
            df[col] = df[col].fillna("").astype(str).str.strip()

    df["chemical_name_norm"] = df["chemical_name"].map(normalise_text)
    df["chemical_id_norm"] = df["chemical_id"].map(normalise_text)
    df["disease_name_norm"] = df["disease_name"].map(normalise_text)
    df["disease_id_norm"] = df["disease_id"].map(normalise_text)

    df["direct_evidence_clean"] = df["direct_evidence"].replace("", "Missing")
    df["inference_gene_symbol_clean"] = df["inference_gene_symbol"].replace("", "Missing")
    df["inference_score_num"] = pd.to_numeric(df["inference_score"], errors="coerce")

    df["pubmed_count"] = df["pubmed_ids"].map(count_pipe_items)
    df["omim_count"] = df["omim_ids"].map(count_pipe_items)

    df["has_direct_evidence"] = df["direct_evidence_clean"].ne("Missing")
    df["has_inference_score"] = df["inference_score_num"].notna()
    df["has_pubmed"] = df["pubmed_count"].gt(0)

    print("[HEAD AFTER NORMALISATION]")
    show_cols = [
        "chemical_name",
        "chemical_id",
        "disease_name",
        "disease_id",
        "direct_evidence_clean",
        "inference_gene_symbol_clean",
        "inference_score_num",
        "pubmed_count",
    ]
    print(df[show_cols].head(10).to_string(index=False))

    return df


def overview_stats(df: pd.DataFrame) -> pd.DataFrame:
    print_section("3. OVERVIEW STATISTICS")

    stats = {
        "total_rows": len(df),
        "unique_chemical_ids": df["chemical_id_norm"].replace("", pd.NA).nunique(dropna=True),
        "unique_chemical_names": df["chemical_name_norm"].replace("", pd.NA).nunique(dropna=True),
        "unique_disease_ids": df["disease_id_norm"].replace("", pd.NA).nunique(dropna=True),
        "unique_disease_names": df["disease_name_norm"].replace("", pd.NA).nunique(dropna=True),
        "unique_chemical_disease_pairs": df[["chemical_id_norm", "disease_id_norm"]].drop_duplicates().shape[0],
        "direct_evidence_records": int(df["has_direct_evidence"].sum()),
        "inferred_records_with_score": int(df["has_inference_score"].sum()),
        "records_with_pubmed_ids": int(df["has_pubmed"].sum()),
        "exact_duplicate_rows": int(df.duplicated().sum()),
    }

    out = pd.DataFrame(list(stats.items()), columns=["metric", "value"])
    print(out.to_string(index=False))
    return out


def top_chemicals_by_unique_diseases(df: pd.DataFrame, top_n: int = 50) -> pd.DataFrame:
    print_section("4. CHEMICALS LINKED TO THE LARGEST NUMBER OF DISTINCT DISEASES")

    out = (
        df.groupby(["chemical_id_norm", "chemical_name_norm"])
        .agg(
            unique_disease_ids=("disease_id_norm", "nunique"),
            unique_disease_names=("disease_name_norm", "nunique"),
            total_records=("chemical_id_norm", "size"),
            direct_evidence_records=("has_direct_evidence", "sum"),
            inferred_records_with_score=("has_inference_score", "sum"),
            records_with_pubmed=("has_pubmed", "sum"),
            mean_inference_score=("inference_score_num", "mean"),
            max_inference_score=("inference_score_num", "max"),
        )
        .reset_index()
        .sort_values(["unique_disease_ids", "total_records"], ascending=False)
        .head(top_n)
    )

    out["mean_inference_score"] = out["mean_inference_score"].round(3)
    out["max_inference_score"] = out["max_inference_score"].round(3)

    print(out.head(25).to_string(index=False))
    return out


def top_diseases_by_unique_chemicals(df: pd.DataFrame, top_n: int = 50) -> pd.DataFrame:
    print_section("5. DISEASES LINKED TO THE LARGEST NUMBER OF DISTINCT CHEMICALS")

    out = (
        df.groupby(["disease_id_norm", "disease_name_norm"])
        .agg(
            unique_chemical_ids=("chemical_id_norm", "nunique"),
            unique_chemical_names=("chemical_name_norm", "nunique"),
            total_records=("disease_id_norm", "size"),
            direct_evidence_records=("has_direct_evidence", "sum"),
            inferred_records_with_score=("has_inference_score", "sum"),
            records_with_pubmed=("has_pubmed", "sum"),
            mean_inference_score=("inference_score_num", "mean"),
            max_inference_score=("inference_score_num", "max"),
        )
        .reset_index()
        .sort_values(["unique_chemical_ids", "total_records"], ascending=False)
        .head(top_n)
    )

    out["mean_inference_score"] = out["mean_inference_score"].round(3)
    out["max_inference_score"] = out["max_inference_score"].round(3)

    print(out.head(25).to_string(index=False))
    return out


def most_common_chemical_disease_pairs(df: pd.DataFrame, top_n: int = 100) -> pd.DataFrame:
    print_section("6. MOST COMMON CHEMICAL–DISEASE PAIRS")

    out = (
        df.groupby(["chemical_id_norm", "chemical_name_norm", "disease_id_norm", "disease_name_norm"])
        .agg(
            total_records=("chemical_id_norm", "size"),
            direct_evidence_types=("direct_evidence_clean", lambda x: "; ".join(sorted(set(x.astype(str)))[:20])),
            inference_genes=("inference_gene_symbol_clean", lambda x: "; ".join(sorted(set(x.astype(str)))[:20])),
            mean_inference_score=("inference_score_num", "mean"),
            max_inference_score=("inference_score_num", "max"),
            total_pubmed_ids=("pubmed_count", "sum"),
        )
        .reset_index()
        .sort_values(["total_records", "total_pubmed_ids"], ascending=False)
        .head(top_n)
    )

    out["mean_inference_score"] = out["mean_inference_score"].round(3)
    out["max_inference_score"] = out["max_inference_score"].round(3)

    print(out.head(30).to_string(index=False))
    return out


def direct_evidence_distribution(df: pd.DataFrame) -> pd.DataFrame:
    print_section("7. DIRECT EVIDENCE DISTRIBUTION")

    out = (
        df.groupby("direct_evidence_clean")
        .agg(
            total_records=("chemical_id_norm", "size"),
            unique_chemicals=("chemical_id_norm", "nunique"),
            unique_diseases=("disease_id_norm", "nunique"),
            unique_pairs=("chemical_id_norm", lambda x: None),
        )
        .reset_index()
        .drop(columns=["unique_pairs"])
        .sort_values("total_records", ascending=False)
    )

    pair_counts = (
        df.groupby("direct_evidence_clean")
        .apply(lambda x: x[["chemical_id_norm", "disease_id_norm"]].drop_duplicates().shape[0], include_groups=False)
        .reset_index(name="unique_pairs")
    )

    out = out.merge(pair_counts, on="direct_evidence_clean", how="left")
    print(out.to_string(index=False))
    return out


def inference_score_distribution(df: pd.DataFrame) -> pd.DataFrame:
    print_section("8. INFERENCE SCORE DISTRIBUTION")

    scored = df[df["inference_score_num"].notna()].copy()

    if scored.empty:
        out = pd.DataFrame(columns=["metric", "value"])
        print("[INFO] No inference scores found.")
        return out

    stats = {
        "records_with_inference_score": len(scored),
        "mean": scored["inference_score_num"].mean(),
        "median": scored["inference_score_num"].median(),
        "min": scored["inference_score_num"].min(),
        "max": scored["inference_score_num"].max(),
        "q25": scored["inference_score_num"].quantile(0.25),
        "q75": scored["inference_score_num"].quantile(0.75),
        "q90": scored["inference_score_num"].quantile(0.90),
        "q95": scored["inference_score_num"].quantile(0.95),
        "q99": scored["inference_score_num"].quantile(0.99),
    }

    out = pd.DataFrame(
        [{"metric": k, "value": round(float(v), 4)} for k, v in stats.items()]
    )

    print(out.to_string(index=False))
    return out


def top_inferred_gene_symbols(df: pd.DataFrame, top_n: int = 50) -> pd.DataFrame:
    print_section("9. TOP INFERENCE GENE SYMBOLS")

    inferred = df[df["inference_gene_symbol_clean"].ne("Missing")].copy()

    out = (
        inferred.groupby("inference_gene_symbol_clean")
        .agg(
            total_records=("chemical_id_norm", "size"),
            unique_chemicals=("chemical_id_norm", "nunique"),
            unique_diseases=("disease_id_norm", "nunique"),
            mean_inference_score=("inference_score_num", "mean"),
            max_inference_score=("inference_score_num", "max"),
        )
        .reset_index()
        .sort_values(["total_records", "unique_diseases"], ascending=False)
        .head(top_n)
    )

    out["mean_inference_score"] = out["mean_inference_score"].round(3)
    out["max_inference_score"] = out["max_inference_score"].round(3)

    print(out.head(25).to_string(index=False))
    return out


def pubmed_support_analysis(df: pd.DataFrame) -> pd.DataFrame:
    print_section("10. PUBMED SUPPORT SUMMARY")

    stats = {
        "records_with_pubmed": int(df["has_pubmed"].sum()),
        "records_without_pubmed": int((~df["has_pubmed"]).sum()),
        "mean_pubmed_ids_per_record": df["pubmed_count"].mean(),
        "median_pubmed_ids_per_record": df["pubmed_count"].median(),
        "max_pubmed_ids_per_record": df["pubmed_count"].max(),
    }

    out = pd.DataFrame(
        [{"metric": k, "value": round(float(v), 4)} for k, v in stats.items()]
    )

    print(out.to_string(index=False))
    return out


def disease_keyword_landscape(df: pd.DataFrame) -> pd.DataFrame:
    print_section("11. BROAD DISEASE KEYWORD LANDSCAPE")

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
                "unique_chemicals": sub["chemical_id_norm"].nunique(),
                "unique_diseases": sub["disease_id_norm"].nunique(),
                "unique_pairs": sub[["chemical_id_norm", "disease_id_norm"]].drop_duplicates().shape[0],
                "direct_evidence_records": int(sub["has_direct_evidence"].sum()),
                "inferred_records_with_score": int(sub["has_inference_score"].sum()),
            }
        )

    out = pd.DataFrame(rows).sort_values("records", ascending=False)
    print(out.to_string(index=False))
    return out


def quality_flags(df: pd.DataFrame) -> pd.DataFrame:
    print_section("12. QUALITY AND INTERPRETATION FLAGS")

    flags = {
        "missing_direct_evidence": df["direct_evidence_clean"].eq("Missing"),
        "has_direct_evidence": df["has_direct_evidence"],
        "has_inference_score": df["has_inference_score"],
        "missing_inference_gene_symbol": df["inference_gene_symbol_clean"].eq("Missing"),
        "missing_pubmed_ids": ~df["has_pubmed"],
        "chemical_id_not_mesh_prefixed": ~df["chemical_id_norm"].str.startswith("mesh:", na=False),
        "disease_id_not_mesh_or_omim": ~(
            df["disease_id_norm"].str.startswith("mesh:", na=False)
            | df["disease_id_norm"].str.startswith("omim:", na=False)
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
                "unique_chemicals": sub["chemical_id_norm"].nunique(),
                "unique_diseases": sub["disease_id_norm"].nunique(),
                "unique_pairs": sub[["chemical_id_norm", "disease_id_norm"]].drop_duplicates().shape[0],
            }
        )

    out = pd.DataFrame(rows).sort_values("records", ascending=False)
    print(out.to_string(index=False))
    return out


def make_figures(
    outdir: Path,
    top_chemicals: pd.DataFrame,
    top_diseases: pd.DataFrame,
    direct_evidence: pd.DataFrame,
    keyword_landscape: pd.DataFrame,
) -> None:
    print_section("13. GENERATING FIGURES")

    figdir = outdir / "figures"
    figdir.mkdir(parents=True, exist_ok=True)

    barplot(
        top_chemicals,
        x_col="chemical_name_norm",
        y_col="unique_disease_ids",
        title="Top CTD chemicals by number of distinct diseases",
        out_path=figdir / "top_chemicals_by_unique_diseases.png",
    )

    barplot(
        top_diseases,
        x_col="disease_name_norm",
        y_col="unique_chemical_ids",
        title="Top CTD diseases by number of distinct chemicals",
        out_path=figdir / "top_diseases_by_unique_chemicals.png",
    )

    barplot(
        direct_evidence,
        x_col="direct_evidence_clean",
        y_col="total_records",
        title="CTD direct evidence distribution",
        out_path=figdir / "direct_evidence_distribution.png",
    )

    barplot(
        keyword_landscape,
        x_col="category",
        y_col="records",
        title="Broad CTD disease keyword landscape",
        out_path=figdir / "disease_keyword_landscape.png",
    )


def write_markdown_report(
    outdir: Path,
    overview: pd.DataFrame,
    top_chemicals: pd.DataFrame,
    top_diseases: pd.DataFrame,
    direct_evidence: pd.DataFrame,
    inference_score: pd.DataFrame,
    flags: pd.DataFrame,
) -> None:
    report_path = outdir / "CTD_Analysis_Report.md"

    with open(report_path, "w", encoding="utf-8") as f:
        f.write("# CTD Analysis Report\n\n")
        f.write(
            "This report summarises exploratory analysis of CTD-derived "
            "chemical-disease association records.\n\n"
        )
        f.write(
            "Important interpretation: CTD records represent chemical-disease "
            "associations and are not necessarily therapeutic drug indications. "
            "Aggregate CTD includes many inferred associations, so provenance "
            "and evidence type should be retained in downstream analyses.\n\n"
        )

        f.write("## Overview statistics\n\n")
        f.write(overview.to_markdown(index=False))
        f.write("\n\n")

        f.write("## Top chemicals by number of distinct diseases\n\n")
        f.write(top_chemicals.head(20).to_markdown(index=False))
        f.write("\n\n")

        f.write("## Top diseases by number of distinct chemicals\n\n")
        f.write(top_diseases.head(20).to_markdown(index=False))
        f.write("\n\n")

        f.write("## Direct evidence distribution\n\n")
        f.write(direct_evidence.to_markdown(index=False))
        f.write("\n\n")

        f.write("## Inference score distribution\n\n")
        f.write(inference_score.to_markdown(index=False))
        f.write("\n\n")

        f.write("## Quality and interpretation flags\n\n")
        f.write(flags.to_markdown(index=False))
        f.write("\n\n")

    print(f"[REPORT] {report_path}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Exploratory analysis of CTD chemical-disease output.")

    parser.add_argument(
        "--input",
        default="./data/processed/ctd_aggregate_chemical_disease.csv",
        help="Path to final CTD processed CSV file.",
    )
    parser.add_argument(
        "--outdir",
        default="./analysis/ctd_analysis",
        help="Output directory for CTD analysis tables and figures.",
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

    print_section("CTD ANALYSIS")
    print(f"[STARTED] {now_utc()}")

    df = pd.read_csv(input_path, dtype=str, low_memory=False)

    metadata = {
        "analysis_name": "CTD Analysis",
        "started_at_utc": now_utc(),
        "interpretation": (
            "CTD records represent chemical-disease associations, not necessarily "
            "therapeutic drug indications."
        ),
    }

    metadata["input_validation"] = validate_input(df, input_path)

    df = add_normalised_columns(df)

    overview = overview_stats(df)
    top_chemicals = top_chemicals_by_unique_diseases(df, top_n=args.top_n)
    top_diseases = top_diseases_by_unique_chemicals(df, top_n=args.top_n)
    common_pairs = most_common_chemical_disease_pairs(df, top_n=args.top_n)
    direct_evidence = direct_evidence_distribution(df)
    inference_score = inference_score_distribution(df)
    top_genes = top_inferred_gene_symbols(df, top_n=args.top_n)
    pubmed_support = pubmed_support_analysis(df)
    keyword_landscape = disease_keyword_landscape(df)
    flags = quality_flags(df)

    save_table(overview, tables_dir / "overview_statistics.csv")
    save_table(top_chemicals, tables_dir / "top_chemicals_by_unique_diseases.csv")
    save_table(top_diseases, tables_dir / "top_diseases_by_unique_chemicals.csv")
    save_table(common_pairs, tables_dir / "most_common_chemical_disease_pairs.csv")
    save_table(direct_evidence, tables_dir / "direct_evidence_distribution.csv")
    save_table(inference_score, tables_dir / "inference_score_distribution.csv")
    save_table(top_genes, tables_dir / "top_inference_gene_symbols.csv")
    save_table(pubmed_support, tables_dir / "pubmed_support_summary.csv")
    save_table(keyword_landscape, tables_dir / "broad_disease_keyword_landscape.csv")
    save_table(flags, tables_dir / "quality_interpretation_flags.csv")

    make_figures(
        outdir=outdir,
        top_chemicals=top_chemicals,
        top_diseases=top_diseases,
        direct_evidence=direct_evidence,
        keyword_landscape=keyword_landscape,
    )

    write_markdown_report(
        outdir=outdir,
        overview=overview,
        top_chemicals=top_chemicals,
        top_diseases=top_diseases,
        direct_evidence=direct_evidence,
        inference_score=inference_score,
        flags=flags,
    )

    metadata["finished_at_utc"] = now_utc()
    metadata["status"] = "success"

    metadata_path = outdir / "ctd_analysis_metadata.json"
    with open(metadata_path, "w", encoding="utf-8") as f:
        json.dump(metadata, f, indent=2)

    print_section("DONE")
    print(f"[OUTPUT DIRECTORY] {outdir}")
    print(f"[METADATA] {metadata_path}")


if __name__ == "__main__":
    main()