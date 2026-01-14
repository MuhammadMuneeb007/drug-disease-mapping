#!/usr/bin/env python3
"""
Open Targets Platform -> Complete DRUG-DISEASE mapping with drug names

Downloads all drug-disease pairs from OpenTargets and adds drug names from local ChEMBL file.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path
from typing import List, Dict, Any

import requests
import pandas as pd
import pyarrow.parquet as pq


UA = "Mozilla/5.0 (OpenTargetsIndicationExporter/1.0)"
BASE = "https://ftp.ebi.ac.uk/pub/databases/opentargets/platform"


def fetch_text(url: str, timeout: int = 60) -> str:
    r = requests.get(url, headers={"User-Agent": UA}, timeout=timeout)
    r.raise_for_status()
    return r.text


def list_parquet_files(dir_url: str) -> List[str]:
    html = fetch_text(dir_url)
    files = re.findall(r'href="([^"]+\.parquet)"', html)
    files = [f for f in files if f.endswith(".parquet")]
    return [dir_url.rstrip("/") + "/" + f for f in files]


def download(url: str, out_path: Path, timeout: int = 300) -> Path:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    if out_path.exists() and out_path.stat().st_size > 0:
        print(f"  [CACHED] {out_path.name}")
        return out_path

    print(f"  [DOWNLOAD] {url.split('/')[-1]}")
    with requests.get(url, headers={"User-Agent": UA}, stream=True, timeout=timeout) as r:
        r.raise_for_status()
        tmp = out_path.with_suffix(out_path.suffix + ".part")
        with open(tmp, "wb") as f:
            for chunk in r.iter_content(chunk_size=1024 * 1024):
                if chunk:
                    f.write(chunk)
        tmp.replace(out_path)
    return out_path


def load_chembl_drug_names(chembl_csv_path: Path) -> Dict[str, str]:
    """Load drug ID -> drug name mapping from local chembl.csv"""
    print(f"[INFO] Loading drug names from: {chembl_csv_path}")
    
    if not chembl_csv_path.exists():
        raise FileNotFoundError(f"ChEMBL file not found: {chembl_csv_path}")
    
    df = pd.read_csv(chembl_csv_path)
    
    if 'drug_id' not in df.columns or 'drug_name' not in df.columns:
        raise ValueError(f"chembl.csv must have 'drug_id' and 'drug_name' columns")
    
    drug_names = {}
    for _, row in df.iterrows():
        drug_id_raw = str(row['drug_id']).strip()
        drug_name = str(row['drug_name']).strip()
        
        if drug_id_raw.startswith('CHEMBL'):
            drug_names[drug_id_raw] = drug_name
        else:
            drug_names[f'CHEMBL{drug_id_raw}'] = drug_name
            drug_names[drug_id_raw] = drug_name
    
    print(f"[INFO] Loaded {len(drug_names):,} drug names")
    return drug_names


def flatten_drug_indications(table: pq.Table, drug_names: Dict[str, str]) -> List[Dict[str, Any]]:
    """Extract all drug-disease pairs from OpenTargets parquet"""
    df = table.to_pandas()
    rows = []
    
    for _, row in df.iterrows():
        drug_id = str(row['id'])
        drug_name = drug_names.get(drug_id, drug_id)
        
        indications = row.get('indications', [])
        if indications is None:
            continue
        
        if not isinstance(indications, list):
            try:
                indications = list(indications)
            except:
                continue
        
        if len(indications) == 0:
            continue
            
        for indication in indications:
            if indication is None or not hasattr(indication, 'get'):
                continue
            
            disease_id = indication.get('disease', '')
            disease_name = indication.get('efoName', '')
            max_phase = indication.get('maxPhaseForIndication', None)
            references = indication.get('references', None)
            
            if max_phase is not None:
                try:
                    max_phase = int(float(max_phase))
                except:
                    max_phase = 0
            else:
                max_phase = 0
            
            sources = []
            if references is not None:
                try:
                    if not isinstance(references, list):
                        references = list(references)
                    for ref in references:
                        if ref is not None and hasattr(ref, 'get'):
                            source = ref.get('source', None)
                            if source:
                                sources.append(str(source))
                except:
                    pass
            
            sources_str = ';'.join(set(sources)) if sources else ''
            
            rows.append({
                'drug_id': drug_id,
                'drug_name': drug_name,
                'disease_id': str(disease_id) if disease_id else '',
                'disease_name': str(disease_name) if disease_name else '',
                'max_phase': max_phase,
                'references': sources_str
            })
    
    return rows


def main():
    ap = argparse.ArgumentParser(description="Download OpenTargets drug-disease mappings with drug names")
    ap.add_argument("--chembl", default="./chembl.csv", help="Path to chembl.csv file")
    ap.add_argument("--out", default="opentargets_drug_disease_complete.csv", help="Output CSV file")
    ap.add_argument("--workdir", default="./opentargets_download", help="Download directory")
    args = ap.parse_args()

    chembl_csv = Path(args.chembl).resolve()
    out_csv = Path(args.out).resolve()
    workdir = Path(args.workdir).resolve()

    # Load drug names
    drug_names = load_chembl_drug_names(chembl_csv)
    
    # Download OpenTargets data
    dir_url = f"{BASE}/latest/output/drug_indication/"
    print(f"\n[INFO] Downloading from OpenTargets: {dir_url}")
    
    parquet_urls = list_parquet_files(dir_url)
    if not parquet_urls:
        raise SystemExit(f"No parquet files found at {dir_url}")
    
    print(f"[INFO] Found {len(parquet_urls)} file(s)")
    
    local_paths = []
    for u in parquet_urls:
        fname = u.split("/")[-1]
        local_paths.append(download(u, workdir / "drug_indication" / fname))
    
    # Process all files
    print(f"\n[INFO] Processing data...")
    all_rows = []
    
    for p in local_paths:
        print(f"[INFO] Reading {p.name}")
        table = pq.read_table(p)
        rows = flatten_drug_indications(table, drug_names)
        all_rows.extend(rows)
        print(f"  -> Extracted {len(rows):,} drug-disease pairs")
    
    # Create DataFrame and save
    print(f"\n[INFO] Creating final dataset...")
    df = pd.DataFrame(all_rows)
    df = df.drop_duplicates(subset=["drug_id", "disease_id", "max_phase"])
    
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_csv, index=False, quoting=1)
    
    print(f"\n? DONE!")
    print(f"Output: {out_csv}")
    print(f"Total rows: {len(df):,}")
    print(f"Unique drugs: {df['drug_id'].nunique():,}")
    print(f"Unique diseases: {df['disease_id'].nunique():,}")
    
    has_names = df[df['drug_name'] != df['drug_id']]
    print(f"Drugs with names: {len(has_names):,} / {len(df):,} ({len(has_names)/len(df)*100:.1f}%)")


if __name__ == "__main__":
    main()