#!/usr/bin/env python3

import polars as pl
from pathlib import Path

# =========================
# Helper: normalize schema
# =========================
def normalize_df(
    df: pl.DataFrame,
    *,
    drug_name_col: str,
    disease_col: str,
    source: str,
    internal_source: str,
    drug_id_col: str | None = None,
) -> pl.DataFrame:

    return df.select([
        pl.col(drug_name_col).cast(pl.Utf8).alias("drug_name"),

        (
            pl.col(drug_id_col).cast(pl.Utf8)
            if drug_id_col and drug_id_col in df.columns
            else pl.lit(None).cast(pl.Utf8)
        ).alias("drug_id"),

        pl.col(disease_col).cast(pl.Utf8).alias("disease"),
        pl.lit(source).alias("source"),
        pl.lit(internal_source).alias("internal_source"),
    ])


# =========================
# Load files
# =========================
base = Path(".")

aact        = pl.read_csv(base / "aact_drug_disease.csv")
ctgov       = pl.read_csv(base / "ctgov_drug_disease.csv")
chembl      = pl.read_csv(base / "chembl.csv")
ctd         = pl.read_csv(base / "ctd_drug_disease.csv")
drugcentral = pl.read_csv(base / "drugcentral_drug_disease.csv")
opentargets = pl.read_csv(base / "opentargets_drug_disease.csv")
sider       = pl.read_csv(base / "sider_drug_indication.csv")

# =========================
# Normalize all sources
# =========================
dfs = [

    normalize_df(
        aact,
        drug_name_col="drug_name",
        disease_col="disease_name",
        source="AACT",
        internal_source="ClinicalTrials.gov",
    ),

    normalize_df(
        ctgov,
        drug_name_col="drug_name",
        disease_col="disease_name",
        source="CTGOV",
        internal_source="ClinicalTrials.gov",
    ),

    normalize_df(
        chembl,
        drug_name_col="drug_name",
        drug_id_col="drug_id",
        disease_col="disease_name",
        source="ChEMBL",
        internal_source="ChEMBL",
    ),

    normalize_df(
        ctd,
        drug_name_col="drug_name",
        drug_id_col="drug_id",
        disease_col="disease_name",
        source="CTD",
        internal_source="CTD_aggregate",
    ),

    normalize_df(
        drugcentral,
        drug_name_col="drug_name",
        drug_id_col="drug_id",
        disease_col="disease_name",
        source="DrugCentral",
        internal_source="DrugCentral",
    ),

    normalize_df(
        opentargets,
        drug_name_col="drug_name",
        drug_id_col="drug_id",
        disease_col="disease_name",
        source="OpenTargets",
        internal_source="OpenTargets",
    ),

    normalize_df(
        sider,
        drug_name_col="drug_name",
        drug_id_col="drug_id",
        disease_col="indication_name",
        source="SIDER",
        internal_source="SIDER4.1",
    ),
]

# =========================
# Merge + normalize text
# =========================
merged = pl.concat(dfs, how="vertical")

merged = merged.with_columns([
    pl.col("drug_name").str.to_lowercase().str.strip_chars(),
    pl.col("disease").str.to_lowercase().str.strip_chars(),
    pl.col("drug_id").str.strip_chars(),
])

# =========================
# DEDUPLICATION (AUTHORITATIVE)
# =========================
dedup = merged.unique(
    subset=["drug_name", "drug_id", "disease"],
    keep="first",
)

# =========================
# Save
# =========================
out_file = "ALL_SOURCES_drug_disease_deduplicated.csv"
dedup.write_csv(out_file)

print(f"[DONE] Rows after deduplication: {dedup.height:,}")
print(f"[SAVED] {out_file}")
