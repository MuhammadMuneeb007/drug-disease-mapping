#!/usr/bin/env python3
"""
ChEMBL SQLite tar.gz  ->  Drug <-> Disease mapping CSV (schema-aware, chunked)

What this does
--------------
1) Finds the SQLite DB inside chembl_XX_sqlite.tar.gz (auto)
2) Extracts ONLY the DB to a temp folder (or keep it if you want)
3) Introspects ChEMBL schema to figure out the correct join keys
4) Exports drug_indication -> (drug, disease) mapping to a CSV in chunks

Works on ChEMBL 36 schema differences (e.g., drug_indication may use molecule_id/molregno not molecule_chembl_id).

Usage
-----
python chembl_sqlite_to_csv.py /path/to/chembl_36_sqlite.tar.gz --out chembl_drug_disease.csv
python chembl_sqlite_to_csv.py /path/to/chembl_36_sqlite.tar.gz --out chembl_drug_disease.csv --keep_db

Notes
-----
- Output includes: drug_id, drug_name, disease_id, disease_name, mesh_id, mesh_heading, max_phase
- Keeps duplicates across sources? Here it's only ChEMBL, but we keep all rows from drug_indication.
"""

from __future__ import annotations

import argparse
import sys
import tarfile
import tempfile
from pathlib import Path
import sqlite3
import pandas as pd


# -----------------------------
# TAR / extraction helpers
# -----------------------------

def find_sqlite_member(tar: tarfile.TarFile) -> tarfile.TarInfo:
    """
    Locate the SQLite database file inside the tar.gz.
    Choose the largest .db/.sqlite/.sqlite3 file as the DB.
    """
    candidates = []
    for m in tar.getmembers():
        if not m.isfile():
            continue
        lname = m.name.lower()
        if lname.endswith(".db") or lname.endswith(".sqlite") or lname.endswith(".sqlite3"):
            candidates.append(m)

    if not candidates:
        raise RuntimeError("No .db/.sqlite file found inside the tar.gz. Is this the ChEMBL sqlite tarball?")

    candidates.sort(key=lambda x: x.size, reverse=True)
    return candidates[0]


def extract_member(tar_path: Path, member: tarfile.TarInfo, out_dir: Path) -> Path:
    """
    Extract a single member from tar into out_dir and return the extracted file path.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    with tarfile.open(tar_path, "r:gz") as tar:
        tar.extract(member, path=out_dir)

    extracted = out_dir / member.name
    if extracted.exists() and extracted.is_file():
        return extracted

    # In case it extracted under nested dirs and we need to find it
    db_files = list(out_dir.rglob("*.db")) + list(out_dir.rglob("*.sqlite")) + list(out_dir.rglob("*.sqlite3"))
    if not db_files:
        raise RuntimeError(f"Extracted member but could not locate DB file under: {out_dir}")
    return max(db_files, key=lambda p: p.stat().st_size)


# -----------------------------
# SQLite schema helpers
# -----------------------------

def table_exists(conn: sqlite3.Connection, table: str) -> bool:
    q = "SELECT 1 FROM sqlite_master WHERE type='table' AND name=? LIMIT 1;"
    return conn.execute(q, (table,)).fetchone() is not None


def get_columns(conn: sqlite3.Connection, table: str) -> set[str]:
    df = pd.read_sql_query(f"PRAGMA table_info({table});", conn)
    return set(df["name"].astype(str))


def choose_first_present(cols: set[str], *candidates: str) -> str | None:
    for c in candidates:
        if c in cols:
            return c
    return None


# -----------------------------
# Core export
# -----------------------------

def build_chembl_indication_sql(conn: sqlite3.Connection) -> str:
    """
    Build a schema-aware SQL query for ChEMBL drug_indication -> drug/disease mapping.
    Handles schema differences across releases.
    """
    if not table_exists(conn, "drug_indication"):
        raise RuntimeError("Missing table: drug_indication")

    di_cols = get_columns(conn, "drug_indication")

    md_exists = table_exists(conn, "molecule_dictionary")
    md_cols = get_columns(conn, "molecule_dictionary") if md_exists else set()

    dd_exists = table_exists(conn, "disease_dictionary")
    dd_cols = get_columns(conn, "disease_dictionary") if dd_exists else set()

    print("[INFO] drug_indication columns contain:", ", ".join(sorted(list(di_cols))[:40]) + (" ..." if len(di_cols) > 40 else ""))
    if md_exists:
        print("[INFO] molecule_dictionary columns contain:", ", ".join(sorted(list(md_cols))[:40]) + (" ..." if len(md_cols) > 40 else ""))
    if dd_exists:
        print("[INFO] disease_dictionary columns contain:", ", ".join(sorted(list(dd_cols))[:40]) + (" ..." if len(dd_cols) > 40 else ""))

    # -------------------------
    # Drug linking: figure join + output drug_id/drug_name
    # -------------------------
    drug_join = None
    drug_id_expr = None

    # Prefer returning a real ChEMBL drug ID if possible
    # Case A: di has molecule_chembl_id directly
    if "molecule_chembl_id" in di_cols:
        drug_id_expr = "di.molecule_chembl_id"
        if md_exists and "molecule_chembl_id" in md_cols:
            drug_join = "di.molecule_chembl_id = md.molecule_chembl_id"

    # Case B: di has molecule_id -> md.molecule_id
    if drug_id_expr is None and ("molecule_id" in di_cols) and md_exists and ("molecule_id" in md_cols):
        drug_join = "di.molecule_id = md.molecule_id"
        if "molecule_chembl_id" in md_cols:
            drug_id_expr = "md.molecule_chembl_id"
        else:
            drug_id_expr = "CAST(md.molecule_id AS TEXT)"

    # Case C: di has molregno -> md.molregno
    if drug_id_expr is None and ("molregno" in di_cols) and md_exists and ("molregno" in md_cols):
        drug_join = "di.molregno = md.molregno"
        if "molecule_chembl_id" in md_cols:
            drug_id_expr = "md.molecule_chembl_id"
        else:
            drug_id_expr = "CAST(md.molregno AS TEXT)"

    # Fallback: no md join possible
    if drug_id_expr is None:
        if "molecule_chembl_id" in di_cols:
            drug_id_expr = "di.molecule_chembl_id"
        elif "molecule_id" in di_cols:
            drug_id_expr = "CAST(di.molecule_id AS TEXT)"
        elif "molregno" in di_cols:
            drug_id_expr = "CAST(di.molregno AS TEXT)"
        else:
            drug_id_expr = "''"

    # Drug name expression
    if md_exists and "pref_name" in md_cols:
        drug_name_expr = "COALESCE(md.pref_name, '')"
    elif "molecule_pref_name" in di_cols:
        drug_name_expr = "COALESCE(di.molecule_pref_name, '')"
    else:
        drug_name_expr = "''"

    # -------------------------
    # Disease fields
    # -------------------------
    # efo_id/efo_term typically exist in drug_indication
    disease_id_expr = "di.efo_id" if "efo_id" in di_cols else "''"

    # disease name from either disease_dictionary.pref_name or di.efo_term
    di_efo_term = "di.efo_term" if "efo_term" in di_cols else "''"

    disease_join = None
    if dd_exists and ("efo_id" in di_cols) and ("efo_id" in dd_cols):
        disease_join = "di.efo_id = dd.efo_id"
        if "pref_name" in dd_cols:
            disease_name_expr = f"COALESCE(dd.pref_name, {di_efo_term}, '')"
        else:
            disease_name_expr = f"COALESCE({di_efo_term}, '')"
    else:
        disease_name_expr = f"COALESCE({di_efo_term}, '')"

    mesh_id_expr = "di.mesh_id" if "mesh_id" in di_cols else "''"
    mesh_heading_expr = "di.mesh_heading" if "mesh_heading" in di_cols else "''"

    phase_col = choose_first_present(di_cols, "max_phase_for_ind", "max_phase", "max_phase_for_indication")
    phase_expr = f"di.{phase_col}" if phase_col else "0"

    # -------------------------
    # Assemble SQL
    # -------------------------
    sql = f"""
    SELECT
        {drug_id_expr}        AS drug_id,
        {drug_name_expr}      AS drug_name,
        {disease_id_expr}     AS disease_id,
        {disease_name_expr}   AS disease_name,
        {mesh_id_expr}        AS mesh_id,
        {mesh_heading_expr}   AS mesh_heading,
        {phase_expr}          AS max_phase
    FROM drug_indication di
    """

    if md_exists and drug_join:
        sql += f"\nLEFT JOIN molecule_dictionary md ON {drug_join}\n"
    if dd_exists and disease_join:
        sql += f"LEFT JOIN disease_dictionary dd ON {disease_join}\n"

    return sql


def export_drug_disease(db_path: Path, out_csv: Path, chunksize: int = 200_000) -> None:
    conn = sqlite3.connect(str(db_path))
    try:
        sql = build_chembl_indication_sql(conn)
        out_csv.parent.mkdir(parents=True, exist_ok=True)

        wrote_header = False
        total = 0

        for chunk in pd.read_sql_query(sql, conn, chunksize=chunksize):
            # Clean
            for c in ["drug_id", "drug_name", "disease_id", "disease_name", "mesh_id", "mesh_heading"]:
                chunk[c] = chunk[c].fillna("").astype(str).str.strip()

            chunk["max_phase"] = pd.to_numeric(chunk["max_phase"], errors="coerce").fillna(0).astype(int)

            # Drop empty core fields
            chunk = chunk[(chunk["drug_id"] != "") & (chunk["disease_name"] != "")]

            chunk.to_csv(out_csv, mode="a", header=(not wrote_header), index=False)
            wrote_header = True

            total += len(chunk)
            print(f"  wrote {total:,} rows -> {out_csv}")

        print("\n? Done")
        print(f"CSV: {out_csv}")
        print(f"Total exported rows: {total:,}")

    finally:
        conn.close()


# -----------------------------
# Main
# -----------------------------

def main():
    ap = argparse.ArgumentParser(description="Export ChEMBL drug_indication (drug<->disease) to CSV from sqlite tar.gz")
    ap.add_argument("chembl_sqlite_tar_gz", type=str, help="Path to chembl_XX_sqlite.tar.gz")
    ap.add_argument("--out", type=str, required=True, help="Output CSV path")
    ap.add_argument("--chunksize", type=int, default=200_000, help="Rows per chunk")
    ap.add_argument("--keep_db", action="store_true", help="Keep extracted DB beside output (default: temp dir)")
    args = ap.parse_args()

    tar_path = Path(args.chembl_sqlite_tar_gz).expanduser().resolve()
    out_csv = Path(args.out).expanduser().resolve()

    if not tar_path.exists():
        print(f"ERROR: tarball not found: {tar_path}", file=sys.stderr)
        sys.exit(1)

    print("=" * 80)
    print("ChEMBL SQLite tar.gz -> Drug-Disease CSV (schema-aware)")
    print("=" * 80)
    print(f"Tarball: {tar_path}")
    print(f"Output:  {out_csv}")

    if args.keep_db:
        extract_dir = out_csv.parent / "_chembl_sqlite_extracted"
        extract_dir.mkdir(parents=True, exist_ok=True)
        temp_ctx = None
    else:
        temp_ctx = tempfile.TemporaryDirectory(prefix="chembl_sqlite_")
        extract_dir = Path(temp_ctx.name)

    # Identify DB inside tar
    with tarfile.open(tar_path, "r:gz") as tar:
        member = find_sqlite_member(tar)
        print(f"Found DB inside tar: {member.name} ({member.size/1e9:.2f} GB)")

    db_path = extract_member(tar_path, member, extract_dir)
    print(f"Using DB: {db_path}")

    export_drug_disease(db_path=db_path, out_csv=out_csv, chunksize=args.chunksize)

    if temp_ctx is not None:
        try:
            temp_ctx.cleanup()
        except Exception:
            pass


if __name__ == "__main__":
    main()
