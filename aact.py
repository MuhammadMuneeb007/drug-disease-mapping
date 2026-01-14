#!/usr/bin/env python3
"""
AACT flat files -> Drug ? Disease mapping (Drug interventions ? Conditions)

Downloads the AACT flat export zip, extracts only interventions + conditions (+ optional studies),
then uses DuckDB to join and write a CSV.

Output columns:
  nct_id, drug_name, disease_name, intervention_type, study_phase, overall_status

Requirements:
  pip install duckdb requests
"""

from __future__ import annotations

import argparse
import zipfile
from pathlib import Path
import requests
import duckdb


AACT_FLAT_ZIP_URL_20260114 = "https://ctti-aact.nyc3.digitaloceanspaces.com/pw7s52pt0ighmd1qcb5hbl8hikcs"


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


def extract_selected(zip_path: Path, out_dir: Path, want: set[str]) -> dict[str, Path]:
    """
    Extract files whose *basename* matches entries in want.
    Returns mapping basename -> extracted path.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    found: dict[str, Path] = {}

    with zipfile.ZipFile(zip_path, "r") as z:
        for member in z.namelist():
            base = Path(member).name
            if base in want:
                dest = out_dir / base
                if not dest.exists() or dest.stat().st_size == 0:
                    z.extract(member, out_dir)
                    # zip may contain subfolders; find the extracted file
                    extracted = next(out_dir.rglob(base))
                    extracted.replace(dest)
                found[base] = dest

    missing = want - set(found.keys())
    if missing:
        raise RuntimeError(
            "Missing expected files in zip: "
            + ", ".join(sorted(missing))
            + "\nTip: list zip contents: python -c \"import zipfile; z=zipfile.ZipFile('file.zip'); print('\\n'.join(z.namelist()[:200]))\""
        )
    return found


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--url", default=AACT_FLAT_ZIP_URL_20260114, help="AACT flat-files zip URL")
    ap.add_argument("--workdir", default="./aact_20260114", help="Working directory")
    ap.add_argument("--out", default="aact_drug_disease.csv", help="Output CSV path")
    ap.add_argument("--keep_all_intervention_types", action="store_true",
                    help="If set, do NOT filter to intervention_type='Drug'")
    args = ap.parse_args()

    workdir = Path(args.workdir).resolve()
    zip_path = workdir / "aact_export.zip"
    extract_dir = workdir / "extract"
    out_csv = Path(args.out).resolve()
    out_csv.parent.mkdir(parents=True, exist_ok=True)

    print(f"Downloading AACT flat files zip...\n  {args.url}")
    download(args.url, zip_path)

    # Common filenames in AACT flat exports
    # If your zip uses different names, list it and adjust here.
    needed = {"interventions.txt", "conditions.txt", "studies.txt"}

    print("Extracting interventions/conditions/studies only...")
    files = extract_selected(zip_path, extract_dir, needed)
    interventions = files["interventions.txt"]
    conditions = files["conditions.txt"]
    studies = files["studies.txt"]

    print("Building mapping with DuckDB (streaming, low RAM)...")
    con = duckdb.connect(database=":memory:")

    # AACT flat files are pipe-delimited.
    # We set delim='|' and header=True.
    con.execute(f"""
        CREATE VIEW interventions AS
        SELECT * FROM read_csv_auto('{interventions.as_posix()}', delim='|', header=True, quote='"', escape='"');
    """)
    con.execute(f"""
        CREATE VIEW conditions AS
        SELECT * FROM read_csv_auto('{conditions.as_posix()}', delim='|', header=True, quote='"', escape='"');
    """)
    con.execute(f"""
        CREATE VIEW studies AS
        SELECT * FROM read_csv_auto('{studies.as_posix()}', delim='|', header=True, quote='"', escape='"');
    """)

    # Inspect expected columns (helps if AACT changes slightly)
    # con.execute("DESCRIBE interventions").fetchall()

    where_clause = ""
    if not args.keep_all_intervention_types:
        where_clause = "WHERE lower(i.intervention_type) = 'drug'"

    # Many AACT tables use: interventions.name, conditions.name, studies.phase, studies.overall_status
    # If a column differs, run: duckdb -c "DESCRIBE interventions" etc.
    sql = f"""
        COPY (
            SELECT
                i.nct_id,
                i.name AS drug_name,
                c.name AS disease_name,
                i.intervention_type,
                COALESCE(s.phase, '') AS study_phase,
                COALESCE(s.overall_status, '') AS overall_status
            FROM interventions i
            JOIN conditions c USING (nct_id)
            LEFT JOIN studies s USING (nct_id)
            {where_clause}
        ) TO '{out_csv.as_posix()}' (HEADER, DELIMITER ',');
    """
    con.execute(sql)
    con.close()

    print("\n? Done")
    print(f"Zip: {zip_path}")
    print(f"CSV: {out_csv}")
    print("\nQuick checks:")
    print(f"  head -n 5 {out_csv}")


if __name__ == "__main__":
    main()
