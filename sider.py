#!/usr/bin/env python3
"""
SIDER 4.1 -> Drug ↔ Indication (Disease/Phenotype) mapping

Downloads:
  - drug_names.tsv
  - meddra_all_indications.tsv.gz

Outputs:
  - sider_drug_indication.csv with columns:
      drug_id, drug_name, indication_id, indication_name, source

Notes:
  - SIDER uses STITCH compound identifiers (derived from PubChem).
  - The MedDRA-derived files are licensed CC BY-NC-SA (non-commercial). Check SIDER download page.

Requirements:
  pip install pandas requests
"""

from __future__ import annotations

import argparse
from pathlib import Path
import gzip
import requests
import pandas as pd


BASE_URL = "http://sideeffects.embl.de/media/download"
FILES = {
    "drug_names": "drug_names.tsv",
    "indications": "meddra_all_indications.tsv.gz",
}


def download(url: str, out_path: Path, timeout: int = 300) -> Path:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    if out_path.exists() and out_path.stat().st_size > 0:
        return out_path
    with requests.get(url, stream=True, timeout=timeout) as r:
        r.raise_for_status()
        tmp = out_path.with_suffix(out_path.suffix + ".part")
        with open(tmp, "wb") as f:
            for chunk in r.iter_content(chunk_size=1024 * 1024):
                if chunk:
                    f.write(chunk)
        tmp.replace(out_path)
    return out_path


def read_tsv_auto(path: Path, gz: bool = False) -> pd.DataFrame:
    if gz:
        return pd.read_csv(path, sep="\t", header=None, compression="gzip", dtype=str)
    return pd.read_csv(path, sep="\t", header=None, dtype=str)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workdir", default="./sider_download", help="Where to download SIDER files")
    ap.add_argument("--out", default="sider_drug_indication.csv", help="Output CSV path")
    args = ap.parse_args()

    workdir = Path(args.workdir).resolve()
    out_csv = Path(args.out).resolve()
    out_csv.parent.mkdir(parents=True, exist_ok=True)

    # Download
    drug_names_path = download(f"{BASE_URL}/{FILES['drug_names']}", workdir / FILES["drug_names"])
    ind_path = download(f"{BASE_URL}/{FILES['indications']}", workdir / FILES["indications"])

    # Read drug_names.tsv
    # Usually: stitch_id \t name
    dn = read_tsv_auto(drug_names_path, gz=False)
    if dn.shape[1] < 2:
        raise SystemExit(f"drug_names.tsv unexpected columns: {dn.shape[1]}")
    dn = dn.rename(columns={0: "drug_id", 1: "drug_name"})
    dn["drug_id"] = dn["drug_id"].astype(str).str.strip()
    dn["drug_name"] = dn["drug_name"].astype(str).str.strip()

    # If multiple names per drug_id, keep the first (or you can aggregate)
    dn = dn.dropna(subset=["drug_id"]).drop_duplicates(subset=["drug_id"], keep="first")

    # Read indications file
    ind = read_tsv_auto(ind_path, gz=True)

    # SIDER formats can vary slightly; we handle by position.
    # Most commonly, meddra_all_indications.tsv.gz has 4 columns like:
    #   drug_id (STITCH) | indication_id (UMLS/MedDRA) | ... | indication_name
    #
    # We’ll map:
    #   col0 -> drug_id
    #   last column -> indication_name
    #   col1 (if exists) -> indication_id
    #
    # If your file has more columns, we keep them but still extract these key ones.
    ncols = ind.shape[1]
    if ncols < 2:
        raise SystemExit(f"meddra_all_indications.tsv.gz unexpected columns: {ncols}")

    ind = ind.copy()
    ind.columns = [f"col{i}" for i in range(ncols)]
    ind["drug_id"] = ind["col0"].astype(str).str.strip()
    ind["indication_id"] = ind["col1"].astype(str).str.strip() if ncols >= 2 else ""
    ind["indication_name"] = ind[f"col{ncols-1}"].astype(str).str.strip()

    # Join names
    out = ind.merge(dn, on="drug_id", how="left")

    # Clean + select
    out["source"] = "SIDER4.1"
    out = out[["drug_id", "drug_name", "indication_id", "indication_name", "source"]]

    # Drop empties
    out = out.dropna(subset=["drug_id", "indication_name"])
    out = out[out["drug_id"].astype(str).str.len() > 0]
    out = out[out["indication_name"].astype(str).str.len() > 0]

    # Deduplicate
    out = out.drop_duplicates(subset=["drug_id", "indication_id", "indication_name"])

    out.to_csv(out_csv, index=False)
    print("✅ Done")
    print(f"Downloaded: {drug_names_path}")
    print(f"Downloaded: {ind_path}")
    print(f"CSV:        {out_csv}")
    print(f"Rows:       {len(out):,}")
    print("\nSanity checks:")
    print(out.head(10).to_string(index=False))


if __name__ == "__main__":
    main()
