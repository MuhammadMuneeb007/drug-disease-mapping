#!/usr/bin/env python3
"""
CTD (Comparative Toxicogenomics Database) -> Chemical/Drug ? Disease mapping

Downloads CTD chemical-disease associations and exports a CSV.

Default uses the *aggregate* file (bigger, includes inference score).
You can switch to curated with --subset curated.

Output columns:
  drug_name, drug_id, disease_name, disease_id, direct_evidence, inference_score, source

Notes:
- CTD is chemical?disease, not strictly drug?disease; many chemicals are drugs.
- CTD TSVs have a comment block and the header line may start with '# '.
"""

from __future__ import annotations

import argparse
import gzip
from pathlib import Path
import requests
import pandas as pd


CTD_URLS = {
    "aggregate": "https://ctdbase.org/reports/CTD_chemicals_diseases.tsv.gz",
    "curated":   "https://ctdbase.org/reports/CTD_curated_chemicals_diseases.tsv.gz",
}


def download(url: str, out_path: Path, timeout: int = 600) -> Path:
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


def find_header_line_gz(path: Path) -> int:
    """
    Find the header line index in a CTD gzipped TSV.

    CTD files:
      - start with many '#' comment lines
      - the TRUE header line may ALSO start with '# ' (e.g. '# ChemicalName\tChemicalID...')

    We locate the header by searching for required column names.
    Return the 0-based line index of that header line.
    """
    required = {"ChemicalName", "ChemicalID", "DiseaseName", "DiseaseID"}

    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as f:
        for i, line in enumerate(f):
            candidate = line.lstrip("#").strip()
            if not candidate:
                continue
            cols = candidate.split("\t")
            if required.issubset(set(cols)):
                return i

    raise RuntimeError("Could not find CTD header line containing expected columns.")


def normalize_columns(df: pd.DataFrame) -> pd.DataFrame:
    """
    Normalize CTD column names:
      - strip whitespace
      - remove leading '# ' if present
    """
    df.columns = [str(c).strip().lstrip("#").strip() for c in df.columns]
    return df


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--subset", choices=["aggregate", "curated"], default="aggregate")
    ap.add_argument("--workdir", default="./ctd_download")
    ap.add_argument("--out", default="ctd_drug_disease.csv")
    ap.add_argument("--chunksize", type=int, default=250_000)
    args = ap.parse_args()

    workdir = Path(args.workdir).resolve()
    out_csv = Path(args.out).resolve()
    out_csv.parent.mkdir(parents=True, exist_ok=True)

    url = CTD_URLS[args.subset]
    gz_path = workdir / Path(url).name

    print(f"Downloading CTD ({args.subset})...")
    download(url, gz_path)
    print(f"Downloaded: {gz_path}")

    header_line = find_header_line_gz(gz_path)
    print(f"[INFO] header_line={header_line} (0-based)")

    wrote_header = False
    total = 0

    # Stream read: skip everything before header_line, then treat next line as header
    for chunk in pd.read_csv(
        gz_path,
        sep="\t",
        compression="gzip",
        skiprows=header_line,  # skips comment lines BEFORE the header
        header=0,              # the header line becomes the header
        dtype=str,
        chunksize=args.chunksize,
        low_memory=False,
    ):
        chunk = normalize_columns(chunk)

        need = ["ChemicalName", "ChemicalID", "DiseaseName", "DiseaseID"]
        missing = [c for c in need if c not in chunk.columns]
        if missing:
            raise RuntimeError(
                f"Missing expected CTD columns: {missing}\nFound: {list(chunk.columns)}"
            )

        # Optional fields
        direct = chunk["DirectEvidence"] if "DirectEvidence" in chunk.columns else ""
        infer  = chunk["InferenceScore"] if "InferenceScore" in chunk.columns else ""

        out = pd.DataFrame({
            "drug_name": chunk["ChemicalName"].fillna("").astype(str),
            "drug_id": chunk["ChemicalID"].fillna("").astype(str),      # MeSH ID (usually)
            "disease_name": chunk["DiseaseName"].fillna("").astype(str),
            "disease_id": chunk["DiseaseID"].fillna("").astype(str),    # MeSH or OMIM
            "direct_evidence": pd.Series(direct).fillna("").astype(str),
            "inference_score": pd.Series(infer).fillna("").astype(str),
            "source": f"CTD_{args.subset}",
        })

        # Drop empties
        out = out[(out["drug_id"].str.len() > 0) & (out["disease_id"].str.len() > 0)]

        # De-dup within chunk (avoid massive memory use across all chunks)
        out = out.drop_duplicates(
            subset=["drug_id", "disease_id", "direct_evidence", "inference_score"]
        )

        out.to_csv(out_csv, mode="a", index=False, header=(not wrote_header))
        wrote_header = True
        total += len(out)

        print(f"  wrote {len(out):,} rows (total {total:,})")

    print("\n? Done")
    print(f"CTD file: {gz_path}")
    print(f"CSV:      {out_csv}")
    print(f"Rows:     {total:,}")


if __name__ == "__main__":
    main()
