#!/usr/bin/env python3
"""
DrugCentral Postgres -> drug_disease CSV (FIXED: fills drug_name via join)

Exports DrugCentral relationships such as:
  - indication
  - contraindication
  - treatment
  - etc.

Key fix:
  joins drug IDs to a drug table (structures/drug/drugs) to populate drug_name.

Run (public instance):
  python drugcentral_to_csv_fixed.py --out drugcentral_drug_disease.csv

Requires:
  pip install psycopg2-binary
"""

from __future__ import annotations

import argparse
import csv
import sys
from typing import Optional, List, Dict, Tuple

import psycopg2
import psycopg2.extras


# -----------------------------
# Introspection helpers
# -----------------------------

def pg_table_exists(cur, table: str, schema: str = "public") -> bool:
    cur.execute(
        """
        SELECT 1
        FROM information_schema.tables
        WHERE table_schema=%s AND table_name=%s
        LIMIT 1
        """,
        (schema, table),
    )
    return cur.fetchone() is not None


def pg_columns(cur, table: str, schema: str = "public") -> List[str]:
    cur.execute(
        """
        SELECT column_name
        FROM information_schema.columns
        WHERE table_schema=%s AND table_name=%s
        ORDER BY ordinal_position
        """,
        (schema, table),
    )
    return [r[0] for r in cur.fetchall()]


def pick_first_present(cols: List[str], *cands: str) -> Optional[str]:
    s = set(cols)
    for c in cands:
        if c in s:
            return c
    return None


def discover_drug_table(cur) -> Tuple[Optional[str], Optional[str], Optional[str]]:
    """
    Find a DrugCentral drug table + join keys.
    Returns (table, id_col, name_col).
    """
    for t in ["structures", "drug", "drugs"]:
        if not pg_table_exists(cur, t):
            continue
        cols = pg_columns(cur, t)
        id_col = pick_first_present(cols, "id", "struct_id", "drug_id")
        name_col = pick_first_present(cols, "name", "drug_name", "inn", "preferred_name")
        if id_col and name_col:
            return (t, id_col, name_col)
    return (None, None, None)


# -----------------------------
# Query builders
# -----------------------------

def build_query_from_omop_relationship(cur) -> Optional[str]:
    """
    Best case: use omop_relationship.
    We join to structures/drug/drugs to fill drug_name.
    """
    if not pg_table_exists(cur, "omop_relationship"):
        return None

    r_cols = pg_columns(cur, "omop_relationship")
    struct_id = pick_first_present(r_cols, "struct_id", "drug_id", "id")
    rel = pick_first_present(r_cols, "relationship_name", "relationship")
    dis_id = pick_first_present(r_cols, "snomed_conceptid", "concept_id", "snomed_id", "disease_id")
    dis_name = pick_first_present(r_cols, "snomed_full_name", "concept_name", "disease_name")

    if not (struct_id and dis_name):
        return None

    drug_table, drug_id_col, drug_name_col = discover_drug_table(cur)

    # Relationship filter: keep only clinically meaningful edges
    # (you can remove this WHERE if you literally want everything)
    where = ""
    if rel:
        where = f"WHERE {rel} IS NOT NULL"

    join = ""
    drug_name_expr = "''"
    if drug_table and drug_id_col and drug_name_col:
        join = f"LEFT JOIN {drug_table} d ON r.{struct_id} = d.{drug_id_col}"
        drug_name_expr = f"COALESCE(d.{drug_name_col}, '')"

    disease_id_expr = f"r.{dis_id}" if dis_id else "''"
    relationship_expr = f"r.{rel}" if rel else "''"

    # Some schemas have a column describing the concept type (disease/procedure/etc.)
    concept_kind = pick_first_present(r_cols, "concept_class", "concept_type", "domain_id", "vocabulary_id")
    concept_kind_expr = f"r.{concept_kind}" if concept_kind else "''"

    return f"""
        SELECT
            CAST(r.{struct_id} AS TEXT)      AS drug_id,
            {drug_name_expr}                AS drug_name,
            {disease_id_expr}               AS disease_id,
            r.{dis_name}                    AS disease_name,
            {relationship_expr}             AS relationship,
            {concept_kind_expr}             AS concept_kind,
            'DrugCentral'                   AS source
        FROM omop_relationship r
        {join}
        {where}
    """


def build_query_from_indication_table(cur) -> Optional[str]:
    """
    Fallback if omop_relationship doesn't exist.
    Try indication-like tables and join to drug table.
    """
    candidates = ["indication", "indications", "drug_indication", "drug_indications", "use", "uses"]
    drug_table, drug_id_col, drug_name_col = discover_drug_table(cur)

    for t in candidates:
        if not pg_table_exists(cur, t):
            continue

        cols = pg_columns(cur, t)

        struct_id = pick_first_present(cols, "struct_id", "drug_id", "structure_id", "id")
        disease_name = pick_first_present(cols, "concept_name", "disease_name", "indication", "name")
        disease_id = pick_first_present(cols, "omop_concept_id", "concept_id", "snomed_conceptid", "mesh_id", "efo_id")
        rel = pick_first_present(cols, "relationship_name", "relationship", "type")

        if not (struct_id and disease_name):
            continue

        join = ""
        drug_name_expr = "''"
        if drug_table and drug_id_col and drug_name_col:
            join = f"LEFT JOIN {drug_table} d ON i.{struct_id} = d.{drug_id_col}"
            drug_name_expr = f"COALESCE(d.{drug_name_col}, '')"

        disease_id_expr = f"i.{disease_id}" if disease_id else "''"
        relationship_expr = f"i.{rel}" if rel else "'indication'"

        return f"""
            SELECT
                CAST(i.{struct_id} AS TEXT)  AS drug_id,
                {drug_name_expr}            AS drug_name,
                {disease_id_expr}           AS disease_id,
                i.{disease_name}            AS disease_name,
                {relationship_expr}         AS relationship,
                ''                          AS concept_kind,
                'DrugCentral'               AS source
            FROM {t} i
            {join}
        """

    return None


def choose_sql(cur) -> str:
    sql = build_query_from_omop_relationship(cur)
    if sql:
        return sql
    sql = build_query_from_indication_table(cur)
    if sql:
        return sql

    cur.execute(
        """
        SELECT table_name
        FROM information_schema.tables
        WHERE table_schema='public'
        ORDER BY table_name
        """
    )
    tables = [r[0] for r in cur.fetchall()]
    raise RuntimeError(
        "Could not find usable DrugCentral mapping tables.\n"
        "First 80 tables:\n" + "\n".join(tables[:80])
    )


# -----------------------------
# Export
# -----------------------------

def export_to_csv(conn, sql: str, out_csv: str, fetch_size: int = 100_000) -> None:
    with conn.cursor(name="dc_stream", cursor_factory=psycopg2.extras.DictCursor) as cur:
        cur.itersize = fetch_size
        cur.execute(sql)

        with open(out_csv, "w", newline="", encoding="utf-8") as f:
            w = csv.writer(f)
            w.writerow(["drug_id", "drug_name", "disease_id", "disease_name", "relationship", "concept_kind", "source"])

            n = 0
            for row in cur:
                drug_id = (row["drug_id"] or "").strip()
                drug_name = (row["drug_name"] or "").strip()
                disease_id = (str(row["disease_id"]) if row["disease_id"] is not None else "").strip()
                disease_name = (row["disease_name"] or "").strip()
                relationship = (row["relationship"] or "").strip()
                concept_kind = (row["concept_kind"] or "").strip()
                source = (row["source"] or "DrugCentral").strip()

                if not disease_name:
                    continue

                w.writerow([drug_id, drug_name, disease_id, disease_name, relationship, concept_kind, source])
                n += 1
                if n % 500_000 == 0:
                    print(f"  wrote {n:,} rows -> {out_csv}")

    print(f"\n? Export complete: {out_csv}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--host", default="unmtid-dbs.net", help="Postgres host")
    ap.add_argument("--port", type=int, default=5433, help="Postgres port")
    ap.add_argument("--user", default="drugman", help="Postgres user")
    ap.add_argument("--password", default="dosage", help="Postgres password")
    ap.add_argument("--dbname", default="drugcentral", help="Database name")
    ap.add_argument("--out", required=True, help="Output CSV path")
    ap.add_argument("--fetch_size", type=int, default=100_000, help="Fetch size for server-side cursor")
    args = ap.parse_args()

    print("=" * 80)
    print("DrugCentral -> Drug-Disease CSV (FIXED drug_name join)")
    print("=" * 80)
    print(f"Connecting: {args.host}:{args.port} db={args.dbname} user={args.user}")
    print(f"Output: {args.out}")

    conn = psycopg2.connect(
        host=args.host,
        port=args.port,
        user=args.user,
        password=args.password,
        dbname=args.dbname,
        connect_timeout=30,
    )

    try:
        with conn.cursor() as cur:
            sql = choose_sql(cur)
            # Print the chosen route for transparency
            print("\n[INFO] Using SQL based on:")
            if "FROM omop_relationship" in sql:
                print("  - omop_relationship (best)")
            else:
                print("  - indication-like table (fallback)")

        export_to_csv(conn, sql, args.out, fetch_size=args.fetch_size)

    finally:
        conn.close()


if __name__ == "__main__":
    main()
