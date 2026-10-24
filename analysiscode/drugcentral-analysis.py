#!/usr/bin/env python3
"""
DrugCentral Analysis: Exploratory analysis of DrugCentral drug-disease/concept relationships.

Input:
  data/processed/drugcentral_drug_disease.csv

Expected columns:
  drug_id
  drug_name
  disease_id
  disease_name
  relationship
  concept_kind
  source
  internal_source
  evidence_type

Purpose:
  1. Validate final processed DrugCentral output.
  2. Summarise DrugCentral relationship coverage.
  3. Identify relationship classes: indication, contraindication, off-label use, etc.
  4. Identify drugs linked to the most diseases/concepts.
  5. Identify diseases/concepts linked to the most drugs.
  6. Generate indication-only and contraindication-only summaries.
  7. Export tables, figures, report, and metadata.

Important interpretation:
  DrugCentral records represent heterogeneous clinical relationships.
  They include contraindications, indications, off-label uses and other labels.
  Therefore, they must not be interpreted as a homogeneous set of therapeutic indications.

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
    "relationship",
    "concept_kind",
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
    if df.empty:
        print(f"[SKIP FIGURE] Empty table: {title}")
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

    empty_df = (
        df.fillna("")
        .eq("")
        .sum()
        .reset_index()
        .rename(columns={"index": "column", 0: "empty_count"})
    )
    empty_df["empty_percent"] = (empty_df["empty_count"] / len(df) * 100).round(3)

    print("\n[MISSING VALUES]")
    print(missing_df.to_string(index=False))

    print("\n[EMPTY STRING VALUES]")
    print(empty_df.to_string(index=False))

    return {
        "input_file": str(input_path),
        "input_size_bytes": input_path.stat().st_size,
        "input_sha256": sha256_file(input_path),
        "n_rows": int(df.shape[0]),
        "n_columns": int(df.shape[1]),
        "columns": list(df.columns),
        "missing_values": missing_df.to_dict(orient="records"),
        "empty_values": empty_df.to_dict(orient="records"),
    }


def add_normalised_columns(df: pd.DataFrame) -> pd.DataFrame:
    print_section("2. CLEANING / NORMALISED ANALYSIS COLUMNS")

    df = df.copy()

    for col in REQUIRED_COLUMNS:
        df[col] = df[col].fillna("").astype(str).str.strip()

    df["drug_id_norm"] = df["drug_id"].map(normalise_text)
    df["drug_name_norm"] = df["drug_name"].map(normalise_text)
    df["disease_id_norm"] = df["disease_id"].map(normalise_text)
    df["disease_name_norm"] = df["disease_name"].map(normalise_text)
    df["relationship_norm"] = df["relationship"].map(normalise_text).replace("", "missing")
    df["concept_kind_norm"] = df["concept_kind"].map(normalise_text).replace("", "missing")

    df["is_indication"] = df["relationship_norm"].eq("indication")
    df["is_contraindication"] = df["relationship_norm"].eq("contraindication")
    df["is_off_label"] = df["relationship_norm"].eq("off-label use")
    df["is_other_relationship"] = ~(
        df["is_indication"] | df["is_contraindication"] | df["is_off_label"]
    )

    print("[HEAD AFTER NORMALISATION]")
    show_cols = [
        "drug_id",
        "drug_name",
        "drug_name_norm",
        "disease_id",
        "disease_name",
        "disease_name_norm",
        "relationship",
        "relationship_norm",
        "concept_kind_norm",
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
        "indication_records": int(df["is_indication"].sum()),
        "contraindication_records": int(df["is_contraindication"].sum()),
        "off_label_use_records": int(df["is_off_label"].sum()),
        "other_relationship_records": int(df["is_other_relationship"].sum()),
        "exact_duplicate_rows": int(df.duplicated().sum()),
    }

    out = pd.DataFrame(list(stats.items()), columns=["metric", "value"])
    print(out.to_string(index=False))
    return out


def relationship_distribution(df: pd.DataFrame) -> pd.DataFrame:
    print_section("4. RELATIONSHIP DISTRIBUTION")

    out = (
        df.groupby("relationship_norm")
        .agg(
            total_records=("relationship_norm", "size"),
            unique_drug_ids=("drug_id_norm", "nunique"),
            unique_drug_names=("drug_name_norm", "nunique"),
            unique_disease_ids=("disease_id_norm", "nunique"),
            unique_disease_names=("disease_name_norm", "nunique"),
        )
        .reset_index()
        .sort_values("total_records", ascending=False)
    )

    print(out.to_string(index=False))
    return out


def concept_kind_distribution(df: pd.DataFrame) -> pd.DataFrame:
    print_section("5. CONCEPT KIND DISTRIBUTION")

    out = (
        df.groupby("concept_kind_norm")
        .agg(
            total_records=("concept_kind_norm", "size"),
            unique_drug_ids=("drug_id_norm", "nunique"),
            unique_disease_names=("disease_name_norm", "nunique"),
        )
        .reset_index()
        .sort_values("total_records", ascending=False)
    )

    print(out.to_string(index=False))
    return out


def top_drugs_by_unique_diseases(df: pd.DataFrame, top_n: int = 50) -> pd.DataFrame:
    print_section("6. DRUGS LINKED TO THE LARGEST NUMBER OF DISTINCT DISEASES/CONCEPTS")

    out = (
        df.groupby(["drug_id_norm", "drug_name_norm"])
        .agg(
            unique_disease_names=("disease_name_norm", "nunique"),
            unique_disease_ids=("disease_id_norm", "nunique"),
            total_records=("drug_id_norm", "size"),
            indication_records=("is_indication", "sum"),
            contraindication_records=("is_contraindication", "sum"),
            off_label_use_records=("is_off_label", "sum"),
        )
        .reset_index()
        .sort_values(["unique_disease_names", "total_records"], ascending=False)
        .head(top_n)
    )

    print(out.head(25).to_string(index=False))
    return out


def top_diseases_by_unique_drugs(df: pd.DataFrame, top_n: int = 50) -> pd.DataFrame:
    print_section("7. DISEASES/CONCEPTS LINKED TO THE LARGEST NUMBER OF DISTINCT DRUGS")

    out = (
        df.groupby(["disease_id_norm", "disease_name_norm"])
        .agg(
            unique_drug_names=("drug_name_norm", "nunique"),
            unique_drug_ids=("drug_id_norm", "nunique"),
            total_records=("disease_name_norm", "size"),
            indication_records=("is_indication", "sum"),
            contraindication_records=("is_contraindication", "sum"),
            off_label_use_records=("is_off_label", "sum"),
        )
        .reset_index()
        .sort_values(["unique_drug_names", "total_records"], ascending=False)
        .head(top_n)
    )

    print(out.head(25).to_string(index=False))
    return out


def most_common_drug_disease_pairs(df: pd.DataFrame, top_n: int = 100) -> pd.DataFrame:
    print_section("8. MOST COMMON DRUG–DISEASE/CONCEPT PAIRS")

    out = (
        df.groupby(["drug_id_norm", "drug_name_norm", "disease_id_norm", "disease_name_norm"])
        .agg(
            total_records=("drug_id_norm", "size"),
            relationship_labels=("relationship_norm", lambda x: "; ".join(sorted(set(x.astype(str))))),
            indication_records=("is_indication", "sum"),
            contraindication_records=("is_contraindication", "sum"),
            off_label_use_records=("is_off_label", "sum"),
        )
        .reset_index()
        .sort_values(["total_records", "indication_records", "contraindication_records"], ascending=False)
        .head(top_n)
    )

    print(out.head(30).to_string(index=False))
    return out


def indication_only_summary(df: pd.DataFrame, top_n: int = 50) -> tuple[pd.DataFrame, pd.DataFrame]:
    print_section("9. INDICATION-ONLY SUMMARY")

    sub = df[df["is_indication"]].copy()

    top_drugs = (
        sub.groupby(["drug_id_norm", "drug_name_norm"])
        .agg(
            unique_indication_names=("disease_name_norm", "nunique"),
            unique_indication_ids=("disease_id_norm", "nunique"),
            total_indication_records=("drug_id_norm", "size"),
        )
        .reset_index()
        .sort_values(["unique_indication_names", "total_indication_records"], ascending=False)
        .head(top_n)
    )

    top_diseases = (
        sub.groupby(["disease_id_norm", "disease_name_norm"])
        .agg(
            unique_drug_names=("drug_name_norm", "nunique"),
            unique_drug_ids=("drug_id_norm", "nunique"),
            total_indication_records=("disease_name_norm", "size"),
        )
        .reset_index()
        .sort_values(["unique_drug_names", "total_indication_records"], ascending=False)
        .head(top_n)
    )

    print("\n[TOP DRUGS BY INDICATION COUNT]")
    print(top_drugs.head(25).to_string(index=False))

    print("\n[TOP INDICATIONS BY DRUG COUNT]")
    print(top_diseases.head(25).to_string(index=False))

    return top_drugs, top_diseases


def contraindication_only_summary(df: pd.DataFrame, top_n: int = 50) -> tuple[pd.DataFrame, pd.DataFrame]:
    print_section("10. CONTRAINDICATION-ONLY SUMMARY")

    sub = df[df["is_contraindication"]].copy()

    top_drugs = (
        sub.groupby(["drug_id_norm", "drug_name_norm"])
        .agg(
            unique_contraindication_names=("disease_name_norm", "nunique"),
            unique_contraindication_ids=("disease_id_norm", "nunique"),
            total_contraindication_records=("drug_id_norm", "size"),
        )
        .reset_index()
        .sort_values(["unique_contraindication_names", "total_contraindication_records"], ascending=False)
        .head(top_n)
    )

    top_diseases = (
        sub.groupby(["disease_id_norm", "disease_name_norm"])
        .agg(
            unique_drug_names=("drug_name_norm", "nunique"),
            unique_drug_ids=("drug_id_norm", "nunique"),
            total_contraindication_records=("disease_name_norm", "size"),
        )
        .reset_index()
        .sort_values(["unique_drug_names", "total_contraindication_records"], ascending=False)
        .head(top_n)
    )

    print("\n[TOP DRUGS BY CONTRAINDICATION COUNT]")
    print(top_drugs.head(25).to_string(index=False))

    print("\n[TOP CONTRAINDICATIONS BY DRUG COUNT]")
    print(top_diseases.head(25).to_string(index=False))

    return top_drugs, top_diseases


def disease_keyword_landscape(df: pd.DataFrame) -> pd.DataFrame:
    print_section("11. BROAD DISEASE/CONCEPT KEYWORD LANDSCAPE")

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
        "diagnostic_or_imaging": r"tomography|imaging|diagnosis|diagnostic|scan|radiograph|spect|pet|mri",
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
                "indication_records": int(sub["is_indication"].sum()),
                "contraindication_records": int(sub["is_contraindication"].sum()),
                "off_label_use_records": int(sub["is_off_label"].sum()),
            }
        )

    out = pd.DataFrame(rows).sort_values("records", ascending=False)
    print(out.to_string(index=False))
    return out


def quality_flags(df: pd.DataFrame) -> pd.DataFrame:
    print_section("12. QUALITY AND INTERPRETATION FLAGS")

    flags = {
        "missing_disease_id": df["disease_id_norm"].eq(""),
        "missing_concept_kind": df["concept_kind_norm"].eq("missing"),
        "exact_duplicate_rows": df.duplicated(),
        "contraindication_records": df["is_contraindication"],
        "indication_records": df["is_indication"],
        "off_label_use_records": df["is_off_label"],
        "non_indication_records": ~df["is_indication"],
        "diagnostic_or_imaging_concepts": df["disease_name_norm"].str.contains(
            r"tomography|imaging|diagnosis|diagnostic|scan|radiograph|spect|pet|mri",
            regex=True,
            na=False,
        ),
    }

    rows = []

    for flag_name, mask in flags.items():
        sub = df[mask]
        rows.append(
            {
                "flag": flag_name,
                "records": int(len(sub)),
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
    relationship: pd.DataFrame,
    top_drugs: pd.DataFrame,
    top_diseases: pd.DataFrame,
    indication_top_diseases: pd.DataFrame,
    contraindication_top_diseases: pd.DataFrame,
    keyword_landscape: pd.DataFrame,
) -> None:
    print_section("13. GENERATING FIGURES")

    figdir = outdir / "figures"
    figdir.mkdir(parents=True, exist_ok=True)

    barplot(
        relationship,
        x_col="relationship_norm",
        y_col="total_records",
        title="DrugCentral relationship distribution",
        out_path=figdir / "relationship_distribution.png",
    )

    barplot(
        top_drugs,
        x_col="drug_name_norm",
        y_col="unique_disease_names",
        title="Top DrugCentral drugs by number of linked disease/concept names",
        out_path=figdir / "top_drugs_by_unique_diseases.png",
    )

    barplot(
        top_diseases,
        x_col="disease_name_norm",
        y_col="unique_drug_names",
        title="Top DrugCentral disease/concept names by number of linked drugs",
        out_path=figdir / "top_diseases_by_unique_drugs.png",
    )

    barplot(
        indication_top_diseases,
        x_col="disease_name_norm",
        y_col="unique_drug_names",
        title="Top indication concepts by number of drugs",
        out_path=figdir / "top_indications_by_unique_drugs.png",
    )

    barplot(
        contraindication_top_diseases,
        x_col="disease_name_norm",
        y_col="unique_drug_names",
        title="Top contraindication concepts by number of drugs",
        out_path=figdir / "top_contraindications_by_unique_drugs.png",
    )

    barplot(
        keyword_landscape,
        x_col="category",
        y_col="records",
        title="Broad DrugCentral disease/concept keyword landscape",
        out_path=figdir / "disease_keyword_landscape.png",
    )


def write_markdown_report(
    outdir: Path,
    overview: pd.DataFrame,
    relationship: pd.DataFrame,
    top_drugs: pd.DataFrame,
    top_diseases: pd.DataFrame,
    indication_top_diseases: pd.DataFrame,
    contraindication_top_diseases: pd.DataFrame,
    flags: pd.DataFrame,
) -> None:
    report_path = outdir / "DrugCentral_Analysis_Report.md"

    with open(report_path, "w", encoding="utf-8") as f:
        f.write("# DrugCentral Analysis Report\n\n")
        f.write(
            "This report summarises exploratory analysis of DrugCentral-derived "
            "drug-disease/concept relationship records.\n\n"
        )
        f.write(
            "Important interpretation: DrugCentral includes heterogeneous relationship "
            "types such as contraindications, indications and off-label uses. These "
            "records should not be treated as a homogeneous set of therapeutic indications.\n\n"
        )

        f.write("## Overview statistics\n\n")
        f.write(overview.to_markdown(index=False))
        f.write("\n\n")

        f.write("## Relationship distribution\n\n")
        f.write(relationship.to_markdown(index=False))
        f.write("\n\n")

        f.write("## Top drugs by number of linked disease/concept names\n\n")
        f.write(top_drugs.head(20).to_markdown(index=False))
        f.write("\n\n")

        f.write("## Top disease/concept names by number of linked drugs\n\n")
        f.write(top_diseases.head(20).to_markdown(index=False))
        f.write("\n\n")

        f.write("## Top indication concepts by number of drugs\n\n")
        f.write(indication_top_diseases.head(20).to_markdown(index=False))
        f.write("\n\n")

        f.write("## Top contraindication concepts by number of drugs\n\n")
        f.write(contraindication_top_diseases.head(20).to_markdown(index=False))
        f.write("\n\n")

        f.write("## Quality and interpretation flags\n\n")
        f.write(flags.to_markdown(index=False))
        f.write("\n\n")

    print(f"[REPORT] {report_path}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Exploratory analysis of DrugCentral output.")

    parser.add_argument(
        "--input",
        default="./data/processed/drugcentral_drug_disease.csv",
        help="Path to final DrugCentral processed CSV file.",
    )
    parser.add_argument(
        "--outdir",
        default="./analysis/drugcentral_analysis",
        help="Output directory for DrugCentral analysis tables and figures.",
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

    print_section("DRUGCENTRAL ANALYSIS")
    print(f"[STARTED] {now_utc()}")

    df = pd.read_csv(input_path, dtype=str, low_memory=False)

    metadata = {
        "analysis_name": "DrugCentral Analysis",
        "started_at_utc": now_utc(),
        "interpretation": (
            "DrugCentral records include heterogeneous clinical relationship types "
            "and should not be treated as homogeneous drug indications."
        ),
    }

    metadata["input_validation"] = validate_input(df, input_path)

    df = add_normalised_columns(df)

    overview = overview_stats(df)
    relationship = relationship_distribution(df)
    concept_kind = concept_kind_distribution(df)
    top_drugs = top_drugs_by_unique_diseases(df, top_n=args.top_n)
    top_diseases = top_diseases_by_unique_drugs(df, top_n=args.top_n)
    common_pairs = most_common_drug_disease_pairs(df, top_n=args.top_n)
    indication_top_drugs, indication_top_diseases = indication_only_summary(df, top_n=args.top_n)
    contraindication_top_drugs, contraindication_top_diseases = contraindication_only_summary(df, top_n=args.top_n)
    keyword_landscape = disease_keyword_landscape(df)
    flags = quality_flags(df)

    save_table(overview, tables_dir / "overview_statistics.csv")
    save_table(relationship, tables_dir / "relationship_distribution.csv")
    save_table(concept_kind, tables_dir / "concept_kind_distribution.csv")
    save_table(top_drugs, tables_dir / "top_drugs_by_unique_diseases.csv")
    save_table(top_diseases, tables_dir / "top_diseases_by_unique_drugs.csv")
    save_table(common_pairs, tables_dir / "most_common_drug_disease_pairs.csv")
    save_table(indication_top_drugs, tables_dir / "indication_top_drugs.csv")
    save_table(indication_top_diseases, tables_dir / "indication_top_diseases.csv")
    save_table(contraindication_top_drugs, tables_dir / "contraindication_top_drugs.csv")
    save_table(contraindication_top_diseases, tables_dir / "contraindication_top_diseases.csv")
    save_table(keyword_landscape, tables_dir / "broad_disease_keyword_landscape.csv")
    save_table(flags, tables_dir / "quality_interpretation_flags.csv")

    make_figures(
        outdir=outdir,
        relationship=relationship,
        top_drugs=top_drugs,
        top_diseases=top_diseases,
        indication_top_diseases=indication_top_diseases,
        contraindication_top_diseases=contraindication_top_diseases,
        keyword_landscape=keyword_landscape,
    )

    write_markdown_report(
        outdir=outdir,
        overview=overview,
        relationship=relationship,
        top_drugs=top_drugs,
        top_diseases=top_diseases,
        indication_top_diseases=indication_top_diseases,
        contraindication_top_diseases=contraindication_top_diseases,
        flags=flags,
    )

    metadata["finished_at_utc"] = now_utc()
    metadata["status"] = "success"

    metadata_path = outdir / "drugcentral_analysis_metadata.json"
    with open(metadata_path, "w", encoding="utf-8") as f:
        json.dump(metadata, f, indent=2)

    print_section("DONE")
    print(f"[OUTPUT DIRECTORY] {outdir}")
    print(f"[METADATA] {metadata_path}")


if __name__ == "__main__":
    main()