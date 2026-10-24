#!/usr/bin/env python3
"""
SIDER Analysis: Exploratory analysis of SIDER drug-disease / indication output.

Input:
  data/processed/sider_drug_disease.csv

Expected columns:
  drug_id
  drug_name
  disease_id
  disease_name
  relationship
  meddra_concept_type
  source_concept_id
  source_concept_name
  extraction_method
  source
  internal_source
  evidence_type

Purpose:
  1. Validate processed SIDER output.
  2. Summarise drug-indication coverage.
  3. Identify drugs linked to the most indications.
  4. Identify indications linked to the most drugs.
  5. Analyse MedDRA concept types, e.g. PT and LLT.
  6. Analyse extraction method distribution.
  7. Export tables, figures, report, and metadata.

Important interpretation:
  SIDER records are label-derived drug-indication / phenotype mappings.
  They should not be interpreted as direct therapeutic efficacy evidence.
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
    "relationship",
    "meddra_concept_type",
    "source_concept_id",
    "source_concept_name",
    "extraction_method",
    "source",
    "internal_source",
    "evidence_type",
]


def now_utc() -> str:
    return datetime.now(timezone.utc).isoformat()


def print_section(title: str) -> None:
    print("\n" + "=" * 100)
    print(title)
    print("=" * 100)


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

    plot_df = df.head(top_n).copy().iloc[::-1]

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
    preview_cols = [
        "drug_id",
        "drug_name",
        "disease_id",
        "disease_name",
        "meddra_concept_type",
        "extraction_method",
        "source",
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
        if col in df.columns:
            df[col] = df[col].fillna("").astype(str).str.strip()

    df["drug_id_norm"] = df["drug_id"].map(normalise_text)
    df["drug_name_norm"] = df["drug_name"].map(normalise_text)
    df["disease_id_norm"] = df["disease_id"].map(normalise_text)
    df["disease_name_norm"] = df["disease_name"].map(normalise_text)
    df["meddra_concept_type_clean"] = df["meddra_concept_type"].replace("", "Missing")
    df["extraction_method_clean"] = df["extraction_method"].replace("", "Missing")

    df["has_drug_name"] = df["drug_name"].ne("")
    df["has_disease_name"] = df["disease_name"].ne("")
    df["is_pt"] = df["meddra_concept_type_clean"].str.upper().eq("PT")
    df["is_llt"] = df["meddra_concept_type_clean"].str.upper().eq("LLT")

    print("[HEAD AFTER NORMALISATION]")
    show_cols = [
        "drug_id",
        "drug_name",
        "disease_id",
        "disease_name",
        "meddra_concept_type_clean",
        "extraction_method_clean",
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
        "unique_drug_disease_pairs": df[["drug_id_norm", "disease_id_norm"]].drop_duplicates().shape[0],
        "records_with_drug_names": int(df["has_drug_name"].sum()),
        "records_with_disease_names": int(df["has_disease_name"].sum()),
        "preferred_term_records_PT": int(df["is_pt"].sum()),
        "lowest_level_term_records_LLT": int(df["is_llt"].sum()),
        "exact_duplicate_rows": int(df.duplicated().sum()),
    }

    out = pd.DataFrame(list(stats.items()), columns=["metric", "value"])
    print(out.to_string(index=False))
    return out


def top_drugs_by_unique_diseases(df: pd.DataFrame, top_n: int = 100) -> pd.DataFrame:
    print_section("4. DRUGS LINKED TO THE LARGEST NUMBER OF DISTINCT INDICATIONS")

    out = (
        df.groupby(["drug_id", "drug_name"])
        .agg(
            unique_disease_ids=("disease_id_norm", "nunique"),
            unique_disease_names=("disease_name_norm", "nunique"),
            total_records=("drug_id", "size"),
            pt_records=("is_pt", "sum"),
            llt_records=("is_llt", "sum"),
            extraction_methods=("extraction_method_clean", lambda x: "; ".join(sorted(set(x.astype(str)))[:20])),
        )
        .reset_index()
        .sort_values(["unique_disease_ids", "total_records"], ascending=False)
        .head(top_n)
    )

    print(out.head(25).to_string(index=False))
    return out


def top_diseases_by_unique_drugs(df: pd.DataFrame, top_n: int = 100) -> pd.DataFrame:
    print_section("5. INDICATIONS LINKED TO THE LARGEST NUMBER OF DISTINCT DRUGS")

    out = (
        df.groupby(["disease_id", "disease_name"])
        .agg(
            unique_drug_ids=("drug_id_norm", "nunique"),
            unique_drug_names=("drug_name_norm", "nunique"),
            total_records=("disease_id", "size"),
            pt_records=("is_pt", "sum"),
            llt_records=("is_llt", "sum"),
            extraction_methods=("extraction_method_clean", lambda x: "; ".join(sorted(set(x.astype(str)))[:20])),
        )
        .reset_index()
        .sort_values(["unique_drug_ids", "total_records"], ascending=False)
        .head(top_n)
    )

    print(out.head(25).to_string(index=False))
    return out


def most_common_pairs(df: pd.DataFrame, top_n: int = 100) -> pd.DataFrame:
    print_section("6. MOST COMMON DRUG-INDICATION PAIRS")

    out = (
        df.groupby(["drug_id", "drug_name", "disease_id", "disease_name"])
        .agg(
            total_records=("drug_id", "size"),
            meddra_concept_types=("meddra_concept_type_clean", lambda x: "; ".join(sorted(set(x.astype(str)))[:20])),
            extraction_methods=("extraction_method_clean", lambda x: "; ".join(sorted(set(x.astype(str)))[:20])),
            source_concept_ids=("source_concept_id", lambda x: "; ".join(sorted(set(x.astype(str)))[:20])),
            source_concept_names=("source_concept_name", lambda x: "; ".join(sorted(set(x.astype(str)))[:20])),
        )
        .reset_index()
        .sort_values(["total_records"], ascending=False)
        .head(top_n)
    )

    print(out.head(30).to_string(index=False))
    return out


def meddra_type_distribution(df: pd.DataFrame) -> pd.DataFrame:
    print_section("7. MEDDRA CONCEPT TYPE DISTRIBUTION")

    out = (
        df.groupby("meddra_concept_type_clean")
        .agg(
            total_records=("drug_id", "size"),
            unique_drugs=("drug_id_norm", "nunique"),
            unique_diseases=("disease_id_norm", "nunique"),
            unique_pairs=("drug_id_norm", lambda x: None),
        )
        .reset_index()
        .drop(columns=["unique_pairs"])
        .sort_values("total_records", ascending=False)
    )

    pair_counts = (
        df.groupby("meddra_concept_type_clean")
        .apply(
            lambda x: x[["drug_id_norm", "disease_id_norm"]].drop_duplicates().shape[0],
            include_groups=False,
        )
        .reset_index(name="unique_pairs")
    )

    out = out.merge(pair_counts, on="meddra_concept_type_clean", how="left")
    print(out.to_string(index=False))
    return out


def extraction_method_distribution(df: pd.DataFrame) -> pd.DataFrame:
    print_section("8. EXTRACTION METHOD DISTRIBUTION")

    out = (
        df.groupby("extraction_method_clean")
        .agg(
            total_records=("drug_id", "size"),
            unique_drugs=("drug_id_norm", "nunique"),
            unique_diseases=("disease_id_norm", "nunique"),
            unique_pairs=("drug_id_norm", lambda x: None),
        )
        .reset_index()
        .drop(columns=["unique_pairs"])
        .sort_values("total_records", ascending=False)
    )

    pair_counts = (
        df.groupby("extraction_method_clean")
        .apply(
            lambda x: x[["drug_id_norm", "disease_id_norm"]].drop_duplicates().shape[0],
            include_groups=False,
        )
        .reset_index(name="unique_pairs")
    )

    out = out.merge(pair_counts, on="extraction_method_clean", how="left")
    print(out.to_string(index=False))
    return out


def disease_keyword_landscape(df: pd.DataFrame) -> pd.DataFrame:
    print_section("9. BROAD DISEASE KEYWORD LANDSCAPE")

    keyword_patterns = {
        "cancer_or_tumour": r"cancer|tumou?r|neoplasm|carcinoma|lymphoma|leukemia|melanoma|sarcoma|myeloma",
        "infection": r"infection|infectious|hiv|hepatitis|covid|tuberculosis|pneumonia|influenza|malaria|viral|bacterial|sepsis",
        "cardiovascular": r"heart|cardiac|cardio|hypertension|coronary|stroke|vascular|myocardial|infarction",
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
                "pt_records": int(sub["is_pt"].sum()),
                "llt_records": int(sub["is_llt"].sum()),
            }
        )

    out = pd.DataFrame(rows).sort_values("records", ascending=False)
    print(out.to_string(index=False))
    return out


def quality_flags(df: pd.DataFrame) -> pd.DataFrame:
    print_section("10. QUALITY AND INTERPRETATION FLAGS")

    flags = {
        "missing_drug_id": df["drug_id"].eq(""),
        "missing_drug_name": df["drug_name"].eq(""),
        "missing_disease_id": df["disease_id"].eq(""),
        "missing_disease_name": df["disease_name"].eq(""),
        "disease_id_is_LLT_or_PT_BAD": df["disease_id"].isin(["LLT", "PT"]),
        "drug_id_looks_like_UMLS_BAD": df["drug_id"].str.match(r"^C\d+", na=False),
        "drug_id_not_CID": ~df["drug_id"].str.startswith("CID", na=False),
        "missing_extraction_method": df["extraction_method_clean"].eq("Missing"),
        "missing_meddra_concept_type": df["meddra_concept_type_clean"].eq("Missing"),
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


def make_figures(
    outdir: Path,
    top_drugs: pd.DataFrame,
    top_diseases: pd.DataFrame,
    meddra_dist: pd.DataFrame,
    method_dist: pd.DataFrame,
    keyword_landscape: pd.DataFrame,
) -> None:
    print_section("11. GENERATING FIGURES")

    figdir = outdir / "figures"
    figdir.mkdir(parents=True, exist_ok=True)

    barplot(
        top_drugs,
        x_col="drug_name",
        y_col="unique_disease_ids",
        title="Top SIDER drugs by number of distinct indications",
        out_path=figdir / "top_drugs_by_unique_indications.png",
    )

    barplot(
        top_diseases,
        x_col="disease_name",
        y_col="unique_drug_ids",
        title="Top SIDER indications by number of distinct drugs",
        out_path=figdir / "top_indications_by_unique_drugs.png",
    )

    barplot(
        meddra_dist,
        x_col="meddra_concept_type_clean",
        y_col="total_records",
        title="SIDER MedDRA concept type distribution",
        out_path=figdir / "meddra_concept_type_distribution.png",
    )

    barplot(
        method_dist,
        x_col="extraction_method_clean",
        y_col="total_records",
        title="SIDER extraction method distribution",
        out_path=figdir / "extraction_method_distribution.png",
    )

    barplot(
        keyword_landscape,
        x_col="category",
        y_col="records",
        title="Broad SIDER disease keyword landscape",
        out_path=figdir / "disease_keyword_landscape.png",
    )


def write_markdown_report(
    outdir: Path,
    overview: pd.DataFrame,
    top_drugs: pd.DataFrame,
    top_diseases: pd.DataFrame,
    meddra_dist: pd.DataFrame,
    method_dist: pd.DataFrame,
    flags: pd.DataFrame,
) -> None:
    report_path = outdir / "SIDER_Analysis_Report.md"

    with open(report_path, "w", encoding="utf-8") as f:
        f.write("# SIDER Drug-Indication Analysis Report\n\n")
        f.write(
            "This report summarises exploratory analysis of SIDER-derived "
            "drug-indication / phenotype mappings.\n\n"
        )
        f.write(
            "Important interpretation: SIDER records are derived from marketed "
            "drug labels and package inserts. They should be interpreted as "
            "label-derived indication/phenotype evidence, not direct therapeutic "
            "efficacy evidence.\n\n"
        )

        f.write("## Overview statistics\n\n")
        f.write(overview.to_markdown(index=False))
        f.write("\n\n")

        f.write("## Top drugs by number of distinct indications\n\n")
        f.write(top_drugs.head(20).to_markdown(index=False))
        f.write("\n\n")

        f.write("## Top indications by number of distinct drugs\n\n")
        f.write(top_diseases.head(20).to_markdown(index=False))
        f.write("\n\n")

        f.write("## MedDRA concept type distribution\n\n")
        f.write(meddra_dist.to_markdown(index=False))
        f.write("\n\n")

        f.write("## Extraction method distribution\n\n")
        f.write(method_dist.to_markdown(index=False))
        f.write("\n\n")

        f.write("## Quality and interpretation flags\n\n")
        f.write(flags.to_markdown(index=False))
        f.write("\n\n")

    print(f"[REPORT] {report_path}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Exploratory analysis of SIDER drug-indication output."
    )

    parser.add_argument(
        "--input",
        default="./data/processed/sider_drug_disease.csv",
        help="Path to final SIDER processed CSV file.",
    )

    parser.add_argument(
        "--outdir",
        default="./analysis/sider_analysis",
        help="Output directory for SIDER analysis tables and figures.",
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

    print_section("SIDER ANALYSIS")
    print(f"[STARTED] {now_utc()}")

    df = pd.read_csv(input_path, dtype=str, low_memory=False)

    metadata = {
        "analysis_name": "SIDER Drug-Indication Analysis",
        "started_at_utc": now_utc(),
        "interpretation": (
            "SIDER records represent label-derived drug-indication/phenotype "
            "mappings and should not be interpreted as direct therapeutic efficacy evidence."
        ),
    }

    metadata["input_validation"] = validate_input(df, input_path)

    df = add_normalised_columns(df)

    overview = overview_stats(df)
    top_drugs = top_drugs_by_unique_diseases(df, top_n=args.top_n)
    top_diseases = top_diseases_by_unique_drugs(df, top_n=args.top_n)
    common_pairs = most_common_pairs(df, top_n=args.top_n)
    meddra_dist = meddra_type_distribution(df)
    method_dist = extraction_method_distribution(df)
    keyword_landscape = disease_keyword_landscape(df)
    flags = quality_flags(df)

    save_table(overview, tables_dir / "overview_statistics.csv")
    save_table(top_drugs, tables_dir / "top_drugs_by_unique_indications.csv")
    save_table(top_diseases, tables_dir / "top_indications_by_unique_drugs.csv")
    save_table(common_pairs, tables_dir / "most_common_drug_indication_pairs.csv")
    save_table(meddra_dist, tables_dir / "meddra_concept_type_distribution.csv")
    save_table(method_dist, tables_dir / "extraction_method_distribution.csv")
    save_table(keyword_landscape, tables_dir / "broad_disease_keyword_landscape.csv")
    save_table(flags, tables_dir / "quality_interpretation_flags.csv")

    make_figures(
        outdir=outdir,
        top_drugs=top_drugs,
        top_diseases=top_diseases,
        meddra_dist=meddra_dist,
        method_dist=method_dist,
        keyword_landscape=keyword_landscape,
    )

    write_markdown_report(
        outdir=outdir,
        overview=overview,
        top_drugs=top_drugs,
        top_diseases=top_diseases,
        meddra_dist=meddra_dist,
        method_dist=method_dist,
        flags=flags,
    )

    metadata["finished_at_utc"] = now_utc()
    metadata["status"] = "success"

    metadata["outputs"] = {
        "outdir": str(outdir),
        "tables_dir": str(tables_dir),
        "report": str(outdir / "SIDER_Analysis_Report.md"),
        "figures_dir": str(outdir / "figures"),
    }

    metadata_path = outdir / "sider_analysis_metadata.json"

    with open(metadata_path, "w", encoding="utf-8") as f:
        json.dump(metadata, f, indent=2, default=str)

    print_section("DONE")
    print(f"[OUTPUT DIRECTORY] {outdir}")
    print(f"[TABLES] {tables_dir}")
    print(f"[REPORT] {outdir / 'SIDER_Analysis_Report.md'}")
    print(f"[METADATA] {metadata_path}")


if __name__ == "__main__":
    main()