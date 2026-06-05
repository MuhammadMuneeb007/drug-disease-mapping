#!/usr/bin/env python3
"""
DrugCentral PostgreSQL -> unified drug-disease association CSV

Purpose
-------
Extract DrugCentral drug-condition relationships into a publication-ready,
provenance-aware CSV for a multi-source drug-disease association dataset.

Source
------
DrugCentral public PostgreSQL instance:
    host:     unmtid-dbs.net
    port:     5433
    database: drugcentral
    user:     drugman
    password: dosage

DrugCentral download page:
    https://drugcentral.org/download

Interpretation
--------------
DrugCentral contains structured drug clinical relationships such as
indications, contraindications, off-label uses, treatment relationships,
and related OMOP/SNOMED concept mappings depending on the database version.

These records should be interpreted according to the original relationship
field. Do not force all DrugCentral rows to approved indications.

Unified output columns
----------------------
drug_name
drug_identifier
drug_identifier_type
disease_or_condition_name
disease_or_condition_identifier
disease_or_condition_identifier_type
relationship_type
evidence_type
source
internal_source

DrugCentral-specific metadata columns
-------------------------------------
drugcentral_struct_id
drugcentral_relationship
concept_kind
source_table
source_relationship_raw

Requirements
------------
pip install psycopg2-binary pandas
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

import pandas as pd
import psycopg2
import psycopg2.extras


# ============================================================
# CONFIG
# ============================================================

DEFAULT_HOST = "unmtid-dbs.net"
DEFAULT_PORT = 5433
DEFAULT_DBNAME = "drugcentral"
DEFAULT_USER = "drugman"
DEFAULT_PASSWORD = "dosage"

DRUGCENTRAL_DOWNLOAD_PAGE = "https://drugcentral.org/download"
DRUGCENTRAL_DUMP_URL = "https://unmtid-dbs.net/download/drugcentral.dump.11012023.sql.gz"

OUTPUT_COLUMNS = [
    # Unified schema
    "drug_name",
    "drug_identifier",
    "drug_identifier_type",
    "disease_or_condition_name",
    "disease_or_condition_identifier",
    "disease_or_condition_identifier_type",
    "relationship_type",
    "evidence_type",
    "source",
    "internal_source",

    # DrugCentral-specific metadata
    "drugcentral_struct_id",
    "drugcentral_relationship",
    "concept_kind",
    "source_table",
    "source_relationship_raw",
]


# ============================================================
# HELPERS
# ============================================================

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


def clean_text(x) -> str:
    if x is None:
        return ""
    return str(x).strip()


def pg_table_exists(cur, table: str, schema: str = "public") -> bool:
    cur.execute(
        """
        SELECT 1
        FROM information_schema.tables
        WHERE table_schema = %s
        AND table_name = %s
        LIMIT 1
        """,
        (schema, table),
    )
    return cur.fetchone() is not None


def pg_columns(cur, table: str, schema: str = "public") -> list[str]:
    cur.execute(
        """
        SELECT column_name
        FROM information_schema.columns
        WHERE table_schema = %s
        AND table_name = %s
        ORDER BY ordinal_position
        """,
        (schema, table),
    )
    return [r[0] for r in cur.fetchall()]


def pick_first_present(columns: list[str], *candidates: str) -> Optional[str]:
    available = set(columns)

    for c in candidates:
        if c in available:
            return c

    return None


def infer_disease_identifier_type(identifier: str) -> str:
    x = clean_text(identifier)

    if not x:
        return ""

    upper = x.upper()

    if upper.startswith("SNOMED"):
        return "SNOMED"

    if upper.startswith("MESH"):
        return "MeSH"

    if upper.startswith("OMOP"):
        return "OMOP"

    if upper.startswith("UMLS"):
        return "UMLS"

    # Most DrugCentral omop_relationship concept IDs are numeric SNOMED/OMOP-like IDs.
    if x.isdigit():
        return "SNOMED/OMOP"

    return "DrugCentral_concept"


def normalise_relationship_type(raw_relationship: str) -> str:
    """
    Preserve the source meaning but make the relationship_type machine-friendly.

    Examples:
        indication -> indication
        contraindication -> contraindication
        off-label use -> off_label_use
    """
    r = clean_text(raw_relationship).lower()

    if not r:
        return "drug_condition_relationship"

    r = r.replace("-", " ")
    r = r.replace("/", " ")
    r = "_".join(r.split())

    return r


def choose_evidence_type(raw_relationship: str) -> str:
    """
    Keep DrugCentral relationship evidence as DrugCentral-derived clinical annotation.
    """
    r = normalise_relationship_type(raw_relationship)

    if "contraindication" in r:
        return "drugcentral_contraindication"

    if "off" in r and "label" in r:
        return "drugcentral_off_label_use"

    if "indication" in r:
        return "drugcentral_indication"

    if "treatment" in r or "treat" in r:
        return "drugcentral_treatment"

    return "drugcentral_clinical_relationship"


# ============================================================
# DATABASE DISCOVERY
# ============================================================

def discover_drug_table(cur) -> tuple[Optional[str], Optional[str], Optional[str]]:
    """
    Find a DrugCentral drug table and join fields.

    Returns:
        table_name, id_column, name_column
    """
    candidates = ["structures", "drug", "drugs"]

    for table in candidates:
        if not pg_table_exists(cur, table):
            continue

        cols = pg_columns(cur, table)

        id_col = pick_first_present(
            cols,
            "id",
            "struct_id",
            "drug_id",
        )

        name_col = pick_first_present(
            cols,
            "name",
            "drug_name",
            "inn",
            "preferred_name",
        )

        if id_col and name_col:
            return table, id_col, name_col

    return None, None, None


def inspect_database(cur, head_n: int = 5) -> dict:
    print_section("Inspecting DrugCentral database")

    cur.execute(
        """
        SELECT table_name
        FROM information_schema.tables
        WHERE table_schema = 'public'
        ORDER BY table_name
        """
    )

    tables = [r[0] for r in cur.fetchall()]

    print(f"[TABLE COUNT] {len(tables):,}")
    print("[FIRST 80 TABLES]")
    for t in tables[:80]:
        print(f"  - {t}")

    inspected = {}

    for table in ["omop_relationship", "structures", "drug", "drugs", "indication", "indications"]:
        if not pg_table_exists(cur, table):
            continue

        cols = pg_columns(cur, table)

        cur.execute(f'SELECT COUNT(*) FROM "{table}"')
        n = cur.fetchone()[0]

        print()
        print("-" * 100)
        print(f"[TABLE] {table}")
        print(f"[ROWS] {n:,}")
        print(f"[COLUMNS] {cols}")

        try:
            cur.execute(f'SELECT * FROM "{table}" LIMIT %s', (head_n,))
            rows = cur.fetchall()
            colnames = [desc[0] for desc in cur.description]
            head = [dict(zip(colnames, row)) for row in rows]
        except Exception:
            head = []

        inspected[table] = {
            "rows": int(n),
            "columns": cols,
            "head_preview": head,
        }

    drug_table, drug_id_col, drug_name_col = discover_drug_table(cur)

    print()
    print("[DRUG TABLE DISCOVERY]")
    print(f"  table:     {drug_table}")
    print(f"  id_col:    {drug_id_col}")
    print(f"  name_col:  {drug_name_col}")

    return {
        "table_count": len(tables),
        "tables_preview": tables[:100],
        "inspected_tables": inspected,
        "drug_table_discovery": {
            "table": drug_table,
            "id_col": drug_id_col,
            "name_col": drug_name_col,
        },
    }


# ============================================================
# QUERY BUILDERS
# ============================================================

def build_query_from_omop_relationship(cur) -> Optional[tuple[str, str]]:
    """
    Best route: use omop_relationship.

    We join to structures/drug/drugs to fill drug_name.
    """
    if not pg_table_exists(cur, "omop_relationship"):
        return None

    r_cols = pg_columns(cur, "omop_relationship")

    struct_id = pick_first_present(
        r_cols,
        "struct_id",
        "drug_id",
        "id",
    )

    rel = pick_first_present(
        r_cols,
        "relationship_name",
        "relationship",
    )

    disease_id = pick_first_present(
        r_cols,
        "snomed_conceptid",
        "concept_id",
        "snomed_id",
        "disease_id",
        "umls_cui",
    )

    disease_name = pick_first_present(
        r_cols,
        "snomed_full_name",
        "concept_name",
        "disease_name",
    )

    concept_kind = pick_first_present(
        r_cols,
        "concept_class",
        "concept_type",
        "domain_id",
        "vocabulary_id",
    )

    if not (struct_id and disease_name):
        return None

    drug_table, drug_id_col, drug_name_col = discover_drug_table(cur)

    join = ""
    drug_name_expr = "''"

    if drug_table and drug_id_col and drug_name_col:
        join = f'LEFT JOIN "{drug_table}" d ON r."{struct_id}" = d."{drug_id_col}"'
        drug_name_expr = f'COALESCE(d."{drug_name_col}", \'\')'

    disease_id_expr = f'CAST(r."{disease_id}" AS TEXT)' if disease_id else "''"
    relationship_expr = f'COALESCE(CAST(r."{rel}" AS TEXT), \'\')' if rel else "''"
    concept_kind_expr = f'COALESCE(CAST(r."{concept_kind}" AS TEXT), \'\')' if concept_kind else "''"

    sql = f"""
        SELECT
            CAST(r."{struct_id}" AS TEXT) AS drugcentral_struct_id,
            {drug_name_expr} AS drug_name,
            {disease_id_expr} AS disease_or_condition_identifier,
            COALESCE(CAST(r."{disease_name}" AS TEXT), '') AS disease_or_condition_name,
            {relationship_expr} AS source_relationship_raw,
            {concept_kind_expr} AS concept_kind,
            'omop_relationship' AS source_table
        FROM omop_relationship r
        {join}
        WHERE r."{disease_name}" IS NOT NULL
        AND TRIM(CAST(r."{disease_name}" AS TEXT)) <> ''
    """

    return sql, "omop_relationship"


def build_query_from_indication_table(cur) -> Optional[tuple[str, str]]:
    """
    Fallback route if omop_relationship is not available.
    """
    candidates = [
        "indication",
        "indications",
        "drug_indication",
        "drug_indications",
        "use",
        "uses",
    ]

    drug_table, drug_id_col, drug_name_col = discover_drug_table(cur)

    for table in candidates:
        if not pg_table_exists(cur, table):
            continue

        cols = pg_columns(cur, table)

        struct_id = pick_first_present(
            cols,
            "struct_id",
            "drug_id",
            "structure_id",
            "id",
        )

        disease_name = pick_first_present(
            cols,
            "concept_name",
            "disease_name",
            "indication",
            "name",
        )

        disease_id = pick_first_present(
            cols,
            "omop_concept_id",
            "concept_id",
            "snomed_conceptid",
            "mesh_id",
            "efo_id",
            "umls_cui",
        )

        rel = pick_first_present(
            cols,
            "relationship_name",
            "relationship",
            "type",
        )

        if not (struct_id and disease_name):
            continue

        join = ""
        drug_name_expr = "''"

        if drug_table and drug_id_col and drug_name_col:
            join = f'LEFT JOIN "{drug_table}" d ON i."{struct_id}" = d."{drug_id_col}"'
            drug_name_expr = f'COALESCE(d."{drug_name_col}", \'\')'

        disease_id_expr = f'CAST(i."{disease_id}" AS TEXT)' if disease_id else "''"
        relationship_expr = f'COALESCE(CAST(i."{rel}" AS TEXT), \'indication\')' if rel else "'indication'"

        sql = f"""
            SELECT
                CAST(i."{struct_id}" AS TEXT) AS drugcentral_struct_id,
                {drug_name_expr} AS drug_name,
                {disease_id_expr} AS disease_or_condition_identifier,
                COALESCE(CAST(i."{disease_name}" AS TEXT), '') AS disease_or_condition_name,
                {relationship_expr} AS source_relationship_raw,
                '' AS concept_kind,
                '{table}' AS source_table
            FROM "{table}" i
            {join}
            WHERE i."{disease_name}" IS NOT NULL
            AND TRIM(CAST(i."{disease_name}" AS TEXT)) <> ''
        """

        return sql, table

    return None


def choose_sql(cur) -> tuple[str, str]:
    route = build_query_from_omop_relationship(cur)

    if route is not None:
        return route

    route = build_query_from_indication_table(cur)

    if route is not None:
        return route

    cur.execute(
        """
        SELECT table_name
        FROM information_schema.tables
        WHERE table_schema = 'public'
        ORDER BY table_name
        """
    )

    tables = [r[0] for r in cur.fetchall()]

    raise RuntimeError(
        "Could not find usable DrugCentral mapping tables.\n"
        "First 100 tables:\n"
        + "\n".join(tables[:100])
    )


# ============================================================
# EXPORT
# ============================================================

def export_to_csv(
    conn,
    sql: str,
    source_table: str,
    out_csv: Path,
    fetch_size: int = 100_000,
) -> dict:
    print_section("Exporting DrugCentral records")

    out_csv.parent.mkdir(parents=True, exist_ok=True)

    if out_csv.exists():
        out_csv.unlink()

    n_written = 0
    n_skipped_no_drug = 0
    n_skipped_no_disease = 0

    relationship_counts = {}
    evidence_counts = {}
    source_table_counts = {}

    with conn.cursor(
        name="drugcentral_stream",
        cursor_factory=psycopg2.extras.DictCursor,
    ) as cur:
        cur.itersize = fetch_size
        cur.execute(sql)

        with open(out_csv, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=OUTPUT_COLUMNS)
            writer.writeheader()

            for row in cur:
                drugcentral_struct_id = clean_text(row["drugcentral_struct_id"])
                drug_name = clean_text(row["drug_name"])

                disease_id = clean_text(row["disease_or_condition_identifier"])
                disease_name = clean_text(row["disease_or_condition_name"])

                source_relationship_raw = clean_text(row["source_relationship_raw"])
                concept_kind = clean_text(row["concept_kind"])
                source_table_value = clean_text(row["source_table"]) or source_table

                if not drugcentral_struct_id and not drug_name:
                    n_skipped_no_drug += 1
                    continue

                if not disease_name:
                    n_skipped_no_disease += 1
                    continue

                relationship_type = normalise_relationship_type(source_relationship_raw)
                evidence_type = choose_evidence_type(source_relationship_raw)

                disease_identifier_type = infer_disease_identifier_type(disease_id)

                output_row = {
                    "drug_name": drug_name,
                    "drug_identifier": drugcentral_struct_id,
                    "drug_identifier_type": "DrugCentral_struct_id",
                    "disease_or_condition_name": disease_name,
                    "disease_or_condition_identifier": disease_id,
                    "disease_or_condition_identifier_type": disease_identifier_type,
                    "relationship_type": relationship_type,
                    "evidence_type": evidence_type,
                    "source": "DrugCentral",
                    "internal_source": "DrugCentral PostgreSQL public instance",
                    "drugcentral_struct_id": drugcentral_struct_id,
                    "drugcentral_relationship": relationship_type,
                    "concept_kind": concept_kind,
                    "source_table": source_table_value,
                    "source_relationship_raw": source_relationship_raw,
                }

                writer.writerow(output_row)
                n_written += 1

                relationship_counts[relationship_type] = relationship_counts.get(relationship_type, 0) + 1
                evidence_counts[evidence_type] = evidence_counts.get(evidence_type, 0) + 1
                source_table_counts[source_table_value] = source_table_counts.get(source_table_value, 0) + 1

                if n_written % 100_000 == 0:
                    print(f"[WRITE] {n_written:,} rows -> {out_csv}")

    if not out_csv.exists() or out_csv.stat().st_size == 0:
        raise RuntimeError("Output CSV was not created or is empty.")

    print_section("Output validation")

    head = pd.read_csv(out_csv, dtype=str, nrows=10)

    print("[OUTPUT HEAD]")
    print(head.to_string(index=False))

    full = pd.read_csv(out_csv, dtype=str, low_memory=False).fillna("")

    exact_duplicate_rows = int(full.duplicated().sum())

    unique_drugs = int(full["drug_identifier"].nunique(dropna=True))
    unique_drug_names = int(full["drug_name"].nunique(dropna=True))
    unique_diseases = int(full["disease_or_condition_identifier"].nunique(dropna=True))
    unique_disease_names = int(full["disease_or_condition_name"].nunique(dropna=True))

    unique_pairs = int(
        full[
            [
                "drug_identifier",
                "disease_or_condition_name",
                "relationship_type",
            ]
        ]
        .drop_duplicates()
        .shape[0]
    )

    print()
    print("[OUTPUT SUMMARY]")
    print(f"Rows written:              {len(full):,}")
    print(f"Skipped missing drug:      {n_skipped_no_drug:,}")
    print(f"Skipped missing disease:   {n_skipped_no_disease:,}")
    print(f"Exact duplicate rows:      {exact_duplicate_rows:,}")
    print(f"Unique drug IDs:           {unique_drugs:,}")
    print(f"Unique drug names:         {unique_drug_names:,}")
    print(f"Unique disease IDs:        {unique_diseases:,}")
    print(f"Unique disease names:      {unique_disease_names:,}")
    print(f"Unique drug-disease-rel pairs: {unique_pairs:,}")
    print(f"Output SHA256:             {sha256_file(out_csv)}")

    print()
    print("[RELATIONSHIP TYPE COUNTS]")
    for k, v in sorted(relationship_counts.items(), key=lambda x: x[1], reverse=True):
        print(f"  - {k}: {v:,}")

    return {
        "output_file": str(out_csv),
        "output_size_bytes": out_csv.stat().st_size,
        "output_sha256": sha256_file(out_csv),
        "rows_written": int(len(full)),
        "rows_stream_written": int(n_written),
        "rows_skipped_no_drug": int(n_skipped_no_drug),
        "rows_skipped_no_disease": int(n_skipped_no_disease),
        "exact_duplicate_rows": exact_duplicate_rows,
        "unique_drug_identifiers": unique_drugs,
        "unique_drug_names": unique_drug_names,
        "unique_disease_identifiers": unique_diseases,
        "unique_disease_names": unique_disease_names,
        "unique_drug_disease_relationship_pairs": unique_pairs,
        "relationship_type_counts": relationship_counts,
        "evidence_type_counts": evidence_counts,
        "source_table_counts": source_table_counts,
        "output_columns": list(full.columns),
        "output_head": head.to_dict(orient="records"),
    }


# ============================================================
# MAIN
# ============================================================

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Export DrugCentral clinical drug-condition relationships into unified schema."
    )

    parser.add_argument("--host", default=DEFAULT_HOST)
    parser.add_argument("--port", type=int, default=DEFAULT_PORT)
    parser.add_argument("--user", default=DEFAULT_USER)
    parser.add_argument("--password", default=DEFAULT_PASSWORD)
    parser.add_argument("--dbname", default=DEFAULT_DBNAME)

    parser.add_argument(
        "--out",
        default="./data/processed/drugcentral_drug_disease.csv",
        help="Output CSV path.",
    )

    parser.add_argument(
        "--metadata",
        default="./data/metadata/drugcentral_processing_metadata.json",
        help="Metadata JSON path.",
    )

    parser.add_argument(
        "--fetch_size",
        type=int,
        default=100_000,
        help="Server-side cursor fetch size.",
    )

    parser.add_argument(
        "--head_n",
        type=int,
        default=5,
        help="Number of rows to preview during inspection.",
    )

    args = parser.parse_args()

    out_csv = Path(args.out).resolve()
    metadata_path = Path(args.metadata).resolve()

    metadata = {
        "script": "drugcentral.py",
        "source": "DrugCentral",
        "source_page": DRUGCENTRAL_DOWNLOAD_PAGE,
        "database_dump_url": DRUGCENTRAL_DUMP_URL,
        "public_postgres": {
            "host": args.host,
            "port": args.port,
            "dbname": args.dbname,
            "user": args.user,
        },
        "started_at_utc": now_utc(),
        "interpretation": (
            "DrugCentral records are structured drug-condition clinical annotations. "
            "The original relationship label is preserved and normalised into "
            "relationship_type. Rows may include indications, contraindications, "
            "off-label uses, treatment relationships, and related clinical mappings "
            "depending on the source table."
        ),
        "unified_output_columns": OUTPUT_COLUMNS,
    }

    conn = None

    try:
        print_section("DrugCentral -> unified drug-disease association processing")
        print(f"[HOST] {args.host}:{args.port}")
        print(f"[DB]   {args.dbname}")
        print(f"[USER] {args.user}")
        print(f"[OUT]  {out_csv}")

        conn = psycopg2.connect(
            host=args.host,
            port=args.port,
            user=args.user,
            password=args.password,
            dbname=args.dbname,
            connect_timeout=60,
        )

        with conn.cursor() as cur:
            metadata["database_inspection"] = inspect_database(cur, head_n=args.head_n)

            sql, source_table = choose_sql(cur)

            metadata["chosen_source_table"] = source_table
            metadata["sql_used"] = sql

            print_section("Chosen DrugCentral extraction route")
            print(f"[SOURCE TABLE] {source_table}")
            print("[SQL USED]")
            print(sql)

        metadata["processing"] = export_to_csv(
            conn=conn,
            sql=sql,
            source_table=source_table,
            out_csv=out_csv,
            fetch_size=args.fetch_size,
        )

        metadata["finished_at_utc"] = now_utc()
        metadata["status"] = "success"

    except Exception as e:
        metadata["finished_at_utc"] = now_utc()
        metadata["status"] = "failed"
        metadata["error"] = repr(e)
        raise

    finally:
        if conn is not None:
            conn.close()

        metadata_path.parent.mkdir(parents=True, exist_ok=True)

        with open(metadata_path, "w", encoding="utf-8") as f:
            json.dump(metadata, f, indent=2, default=str)

        print_section("METADATA WRITTEN")
        print(metadata_path)


if __name__ == "__main__":
    main()