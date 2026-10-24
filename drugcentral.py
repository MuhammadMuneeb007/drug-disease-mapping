
#!/usr/bin/env python3
"""DrugCentral PostgreSQL -> frozen raw snapshot -> merge-compatible evidence.

First run downloads a query-selected snapshot from the public DrugCentral database.
Later runs reuse that snapshot. Pass --refresh-source to fetch again, or --offline
for an explicitly network-free run. No existing processed files are deleted early.

Indications, contraindications and off-label uses are NOT interchangeable.
This is a source-evidence table, not an ontology-resolved/independent-efficacy table.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(os.environ.get(
    "DRUG_DISEASE_ROOT", Path(__file__).resolve().parent
)).resolve()

RAW_DIR = ROOT / "data/raw/drugcentral"
OUT_DIR = ROOT / "data/processed"
META_DIR = ROOT / "data/metadata"

RAW_NAME = "drugcentral_omop_relationship_snapshot.csv"
MANIFEST_NAME = "drugcentral_omop_relationship_snapshot.json"

CORE = [
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
]

EXTRA = [
    "drugcentral_struct_id",
    "drugcentral_relationship",
    "concept_kind",
    "source_table",
    "source_relationship_raw",
    "relationship_category",
    "drugcentral_record_id",
    "drugcentral_source_row_number",
    "disease_identifier_source_column",
    "snomed_conceptid",
    "umls_cui",
    "omop_concept_id",
    "concept_name",
    "snomed_full_name",
    "drug_name_lookup_status",
    "disease_name_lookup_status",
    "evidence_granularity",
]

for entity in ("drug", "disease"):
    EXTRA += [
        f"{entity}_identifier_{s}"
        for s in (
            "original",
            "namespace",
            "canonical",
            "validation",
        )
    ]

COLUMNS = CORE + EXTRA

RAW_COLUMNS = [
    "drugcentral_source_row_number",
    "drugcentral_record_id",
    "drugcentral_struct_id",
    "drug_name",
    "snomed_conceptid",
    "umls_cui",
    "omop_concept_id",
    "concept_name",
    "snomed_full_name",
    "disease_or_condition_name",
    "source_relationship_raw",
    "concept_kind",
    "source_table",
]

REJECT_COLS = [
    "source_row_number",
    "source_record_id",
    "reason",
    "raw_record_json",
]

DEDUP_EXTRA = [
    "drugcentral_supporting_record_count",
    "drugcentral_supporting_source_rows",
]

RELATIONS = {
    "indication": (
        "indication",
        "drugcentral_indication",
    ),
    "contraindication": (
        "contraindication",
        "drugcentral_contraindication",
    ),
    "off_label_use": (
        "off_label_use",
        "drugcentral_off_label_use",
    ),
    "symptomatic_treatment": (
        "symptomatic_treatment",
        "drugcentral_symptomatic_treatment",
    ),
    "reduce_risk": (
        "risk_reduction",
        "drugcentral_risk_reduction",
    ),
    "diagnosis": (
        "diagnostic_use",
        "drugcentral_diagnostic_use",
    ),
}


def now():
    return datetime.now(timezone.utc).isoformat()


def heading(label):
    print(
        "\n" + "=" * 84
        + "\n" + label
        + "\n" + "=" * 84,
        flush=True,
    )


def clean(x):
    return "" if x is None else str(x).strip()


def digest(path):
    h = hashlib.sha256()

    with Path(path).open("rb") as f:
        for b in iter(
            lambda: f.read(4 * 1024 * 1024),
            b"",
        ):
            h.update(b)

    return h.hexdigest()


def stage_json(path, content):
    path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    tmp = path.with_name(
        path.name + ".processing.tmp"
    )

    tmp.write_text(
        json.dumps(
            content,
            indent=2,
            ensure_ascii=False,
            default=str,
        ) + "\n",
        encoding="utf-8",
    )

    return tmp


def identifier(raw, namespace):
    original = clean(raw)

    if not original:
        return original, namespace or "", "", "missing"

    v = original.upper()

    if namespace == "DrugCentral":
        good = bool(
            re.fullmatch(r"\d+", v)
        )

        return (
            original,
            namespace,
            "DRUGCENTRAL:" + v if good else "",
            "syntax_valid" if good else "invalid_syntax",
        )

    if namespace == "SNOMEDCT":
        good = bool(
            re.fullmatch(r"\d{6,18}", v)
        )

        return (
            original,
            namespace,
            "SNOMEDCT:" + v if good else "",
            "syntax_valid" if good else "invalid_syntax",
        )

    if namespace == "UMLS":
        good = bool(
            re.fullmatch(
                r"(?:UMLS:)?C\d{7}",
                v,
            )
        )

        return (
            original,
            namespace,
            "UMLS:" + v.removeprefix("UMLS:") if good else "",
            "syntax_valid" if good else "invalid_syntax",
        )

    if namespace == "OMOP":
        good = bool(
            re.fullmatch(r"\d+", v)
        )

        return (
            original,
            namespace,
            "OMOP:" + v if good else "",
            "syntax_valid" if good else "invalid_syntax",
        )

    return (
        original,
        namespace or "unresolved",
        "",
        "unrecognised_namespace",
    )


def relationship(raw):
    token = re.sub(
        r"[^a-z0-9]+",
        "_",
        clean(raw).casefold(),
    ).strip("_")

    rel, ev = RELATIONS.get(
        token,
        (
            "unclassified_drug_condition_relationship",
            "drugcentral_unclassified_relationship",
        ),
    )

    return rel, ev, token or "missing"


def from_raw(row):
    drug_id = clean(
        row["drugcentral_struct_id"]
    )

    drug_name = clean(
        row["drug_name"]
    )

    if drug_name:
        drug_lookup = "structures_name"

    elif drug_id:
        drug_name = drug_id
        drug_lookup = "struct_id_placeholder"

    else:
        drug_lookup = "missing_drug_identity"

    ids = [
        ("snomed_conceptid", "SNOMEDCT"),
        ("umls_cui", "UMLS"),
        ("omop_concept_id", "OMOP"),
    ]

    disease_id = ""
    disease_ns = ""
    id_col = ""

    for col, namespace in ids:
        if clean(row[col]):
            disease_id = clean(row[col])
            disease_ns = namespace
            id_col = col
            break

    disease_name = clean(
        row["disease_or_condition_name"]
    )

    if disease_name:
        disease_lookup = "source_concept_name"

    elif disease_id:
        disease_name = disease_id
        disease_lookup = "identifier_placeholder"

    else:
        disease_lookup = "missing_disease_identity"

    reasons = []

    if not drug_id:
        reasons.append("missing_struct_id")

    if not disease_name:
        reasons.append(
            "missing_disease_name_and_identifier"
        )

    if reasons:
        return None, "|".join(reasons)

    rel, evidence, normalized = relationship(
        row["source_relationship_raw"]
    )

    drug_fields = identifier(
        drug_id,
        "DrugCentral",
    )

    disease_fields = identifier(
        disease_id,
        disease_ns,
    )

    out = {
        "drug_name": drug_name,
        "drug_identifier": drug_id,
        "drug_identifier_type": "DrugCentral",

        "disease_or_condition_name":
            disease_name,

        "disease_or_condition_identifier":
            disease_id,

        "disease_or_condition_identifier_type":
            disease_ns,

        "relationship_type": rel,
        "evidence_type": evidence,

        "source": "DrugCentral",

        "internal_source":
            "DrugCentral omop_relationship",

        "drugcentral_struct_id": drug_id,

        "drugcentral_relationship": normalized,

        "concept_kind": clean(
            row["concept_kind"]
        ),

        "source_table": "omop_relationship",

        "source_relationship_raw": clean(
            row["source_relationship_raw"]
        ),

        "relationship_category": rel,

        "drugcentral_record_id": clean(
            row["drugcentral_record_id"]
        ),

        "drugcentral_source_row_number": clean(
            row["drugcentral_source_row_number"]
        ),

        "disease_identifier_source_column":
            id_col,

        "snomed_conceptid": clean(
            row["snomed_conceptid"]
        ),

        "umls_cui": clean(
            row["umls_cui"]
        ),

        "omop_concept_id": clean(
            row["omop_concept_id"]
        ),

        "concept_name": clean(
            row["concept_name"]
        ),

        "snomed_full_name": clean(
            row["snomed_full_name"]
        ),

        "drug_name_lookup_status":
            drug_lookup,

        "disease_name_lookup_status":
            disease_lookup,

        "evidence_granularity":
            "source_omop_relationship_row",
    }

    for entity, data in (
        ("drug", drug_fields),
        ("disease", disease_fields),
    ):
        for suffix, val in zip(
            (
                "original",
                "namespace",
                "canonical",
                "validation",
            ),
            data,
        ):
            out[
                f"{entity}_identifier_{suffix}"
            ] = val

    return out, ""


def safe_ident(s):
    if not re.fullmatch(
        r"[a-z_][a-z0-9_]*",
        s,
    ):
        raise RuntimeError(
            f"Unexpected PostgreSQL identifier {s!r}"
        )

    return '"' + s + '"'


def pg_columns(cur, table):
    cur.execute(
        "SELECT column_name "
        "FROM information_schema.columns "
        "WHERE table_schema='public' "
        "AND table_name=%s "
        "ORDER BY ordinal_position",
        (table,),
    )

    return [
        r[0]
        for r in cur.fetchall()
    ]


def text_sql(alias, name):
    if not name:
        return "''"

    return (
        f"COALESCE(TRIM(CAST("
        f"{alias}.{safe_ident(name)} "
        f"AS TEXT)), '')"
    )


def fetch_snapshot(args, temp_raw):
    try:
        import psycopg2
    except ImportError as exc:
        raise RuntimeError(
            "Install psycopg2-binary to fetch "
            "the DrugCentral PostgreSQL snapshot"
        ) from exc

    heading(
        "DATABASE SOURCE INSPECTION "
        "AND FROZEN EXTRACTION"
    )

    conn = psycopg2.connect(
        host=args.host,
        port=args.port,
        dbname=args.dbname,
        user=args.user,
        password=args.password,
        connect_timeout=45,
        application_name="drugcentral_source_event_export",
    )

    try:
        conn.set_session(
            readonly=True,
            isolation_level="REPEATABLE READ",
            autocommit=False,
        )

        with conn.cursor() as c:
            cols = pg_columns(
                c,
                "omop_relationship",
            )

            if not cols:
                raise RuntimeError(
                    "public.omop_relationship absent; "
                    "refusing unverified alternative table"
                )

            required = {
                "struct_id",
                "relationship_name",
            }

            if (
                not required.issubset(cols)
                or not (
                    {
                        "snomed_full_name",
                        "concept_name",
                    } & set(cols)
                )
            ):
                raise RuntimeError(
                    "Unrecognised omop_relationship "
                    f"schema: {cols}"
                )

            print(
                "[SOURCE HEADER] omop_relationship:",
                " | ".join(cols),
            )

            struct_cols = pg_columns(
                c,
                "structures",
            )

            if (
                struct_cols
                and {
                    "id",
                    "name",
                }.issubset(struct_cols)
            ):
                c.execute(
                    'SELECT COUNT(*) FROM '
                    '(SELECT id FROM public.structures '
                    'GROUP BY id HAVING COUNT(*) > 1) x'
                )

                duplicate_ids = c.fetchone()[0]

                if duplicate_ids:
                    raise RuntimeError(
                        f"structures.id has "
                        f"{duplicate_ids} duplicate keys"
                    )

                join = (
                    'LEFT JOIN public.structures d '
                    'ON r."struct_id" = d."id"'
                )

                drug_name = text_sql(
                    "d",
                    "name",
                )

                name_strategy = (
                    "structures.id unique LEFT JOIN "
                    "omop_relationship.struct_id"
                )

            else:
                join = ""
                drug_name = "''"

                name_strategy = (
                    "no structures.id/name; "
                    "use flagged struct_id placeholder"
                )

            print(
                "[STRUCTURES HEADER]",
                struct_cols,
            )

            print(
                "[DRUG NAME JOIN]",
                name_strategy,
            )

            c.execute(
                "SELECT COUNT(*) "
                "FROM public.omop_relationship"
            )

            original_count = int(
                c.fetchone()[0]
            )

            dbversion = []

            if pg_columns(
                c,
                "dbversion",
            ):
                c.execute(
                    "SELECT * FROM public.dbversion"
                )

                dbcols = [
                    desc[0]
                    for desc in c.description
                ]

                dbversion = [
                    dict(zip(dbcols, r))
                    for r in c.fetchall()
                ]

            print(
                "[ORIGINAL SOURCE ROWS]",
                f"{original_count:,}",
            )

            print(
                "[DBVERSION]",
                dbversion,
            )

            def optional(name):
                return text_sql(
                    "r",
                    name if name in cols else None,
                )

            name_expr = (
                f"COALESCE("
                f"NULLIF({optional('snomed_full_name')}, ''), "
                f"{optional('concept_name')})"
            )

            query = f"""
                SELECT
                    {optional('id')} AS drugcentral_record_id,
                    {text_sql('r', 'struct_id')} AS drugcentral_struct_id,
                    {drug_name} AS drug_name,
                    {optional('snomed_conceptid')} AS snomed_conceptid,
                    {optional('umls_cui')} AS umls_cui,
                    {optional('concept_id')} AS omop_concept_id,
                    {optional('concept_name')} AS concept_name,
                    {optional('snomed_full_name')} AS snomed_full_name,
                    {name_expr} AS disease_or_condition_name,
                    {text_sql('r', 'relationship_name')} AS source_relationship_raw,
                    {optional('cui_semantic_type')} AS concept_kind,
                    'omop_relationship' AS source_table
                FROM public.omop_relationship r
                {join}
                ORDER BY r."struct_id", r."relationship_name"
                {', r."concept_id" NULLS FIRST' if 'concept_id' in cols else ''}
            """

            c.execute(
                f"SELECT COUNT(*) FROM ({query}) q"
            )

            joined_count = int(
                c.fetchone()[0]
            )

            if joined_count != original_count:
                raise RuntimeError(
                    "Join changes row count "
                    f"{original_count} -> {joined_count}"
                )

        count = 0

        temp_raw.parent.mkdir(
            parents=True,
            exist_ok=True,
        )

        with temp_raw.open(
            "w",
            encoding="utf-8",
            newline="",
        ) as output:

            writer = csv.DictWriter(
                output,
                fieldnames=RAW_COLUMNS,
            )

            writer.writeheader()

            with conn.cursor(
                name="drugcentral_snapshot_stream"
            ) as stream:

                stream.itersize = args.fetch_size
                stream.execute(query)

                for fields in iter(
                    lambda: stream.fetchmany(
                        args.fetch_size
                    ),
                    [],
                ):
                    for source_row in fields:
                        count += 1

                        row = dict(
                            zip(
                                RAW_COLUMNS[1:],
                                source_row,
                            )
                        )

                        writer.writerow({
                            "drugcentral_source_row_number":
                                count,
                            **row,
                        })

                    if count % 100000 < len(fields):
                        print(
                            "[SNAPSHOT EXPORTED]",
                            f"{count:,}",
                            flush=True,
                        )

        if (
            count != joined_count
            or count == 0
        ):
            raise RuntimeError(
                "Snapshot row conservation failure: "
                f"{count} != {joined_count}"
            )

        manifest = {
            "source": "DrugCentral",

            "source_table":
                "public.omop_relationship",

            "captured_at_utc": now(),

            "dbversion": dbversion,

            "source_table_header": cols,

            "structures_header": struct_cols,

            "name_join": name_strategy,

            "source_rows": original_count,

            "rows_after_name_join": joined_count,

            "raw_snapshot_rows": count,

            "raw_snapshot_columns": RAW_COLUMNS,

            "query": query,

            "raw_snapshot_sha256":
                digest(temp_raw),

            "access_host": args.host,

            "access_port": args.port,

            "access_database": args.dbname,
        }

        print(
            "[FROZEN RAW SNAPSHOT]",
            f"{count:,}",
            "rows",
            manifest["raw_snapshot_sha256"],
        )

        return manifest

    finally:
        conn.close()


def process_source(
    path,
    expected_rows,
    out_tmp,
    dedup_tmp,
    rejects_tmp,
):
    heading(
        "PROCESS FROZEN DRUGCENTRAL SOURCE ROWS"
    )

    counts = Counter()
    relations = Counter()
    raw_rels = Counter()
    id_sources = Counter()
    drug_names = Counter()
    diseases = Counter()
    validations = Counter()

    unique_pairs = set()
    seen = {}

    example = None

    with (
        path.open(
            "r",
            encoding="utf-8-sig",
            newline="",
        ) as stream,

        out_tmp.open(
            "w",
            encoding="utf-8",
            newline="",
        ) as out,

        rejects_tmp.open(
            "w",
            encoding="utf-8",
            newline="",
        ) as reject,
    ):
        reader = csv.reader(
            stream,
            strict=True,
        )

        header = next(reader)

        if header != RAW_COLUMNS:
            raise RuntimeError(
                "Frozen source header differs: "
                f"{header}"
            )

        writer = csv.DictWriter(
            out,
            fieldnames=COLUMNS,
        )

        bad_writer = csv.DictWriter(
            reject,
            fieldnames=REJECT_COLS,
        )

        writer.writeheader()
        bad_writer.writeheader()

        for ordinal, fields in enumerate(
            reader,
            1,
        ):
            if len(fields) != len(RAW_COLUMNS):
                raise RuntimeError(
                    "Malformed frozen row "
                    f"{ordinal}: {len(fields)} fields"
                )

            item = dict(
                zip(RAW_COLUMNS, fields)
            )

            if (
                item["drugcentral_source_row_number"]
                != str(ordinal)
            ):
                raise RuntimeError(
                    "Unexpected source record "
                    f"ordinal at row {ordinal}"
                )

            counts["input_rows"] += 1

            rec, reason = from_raw(item)

            if reason:
                counts["rejected_rows"] += 1

                bad_writer.writerow({
                    "source_row_number": ordinal,

                    "source_record_id":
                        item["drugcentral_record_id"],

                    "reason": reason,

                    "raw_record_json": json.dumps(
                        item,
                        ensure_ascii=False,
                    ),
                })

                continue

            writer.writerow(rec)

            counts["retained_rows"] += 1

            if example is None:
                example = rec

            relations[
                rec["relationship_type"]
            ] += 1

            raw_rels[
                rec["source_relationship_raw"]
                or "[missing]"
            ] += 1

            id_sources[
                rec["disease_identifier_source_column"]
                or "[missing]"
            ] += 1

            drug_names[
                rec["drug_name_lookup_status"]
            ] += 1

            diseases[
                rec["disease_name_lookup_status"]
            ] += 1

            validations[
                rec["disease_identifier_validation"]
            ] += 1

            unique_pairs.add((
                rec["drug_identifier"],
                rec["disease_or_condition_identifier"],
                rec["disease_or_condition_identifier_type"],
                rec["relationship_type"],
            ))

            # Dedup is a separate reporting view.
            # The main output preserves every event.

            key = tuple(
                rec[x]
                for x in (
                    "drugcentral_struct_id",
                    "snomed_conceptid",
                    "umls_cui",
                    "omop_concept_id",
                    "source_relationship_raw",
                    "concept_kind",
                    "concept_name",
                    "snomed_full_name",
                )
            )

            if key not in seen:
                seen[key] = [rec, []]

            seen[key][1].append(
                rec["drugcentral_source_row_number"]
            )

    if (
        counts["input_rows"] != expected_rows
        or (
            counts["retained_rows"]
            + counts["rejected_rows"]
        ) != expected_rows
    ):
        raise RuntimeError(
            "Source row accounting failed: "
            f"{dict(counts)}, expected {expected_rows}"
        )

    if not counts["retained_rows"]:
        raise RuntimeError(
            "No eligible DrugCentral rows; "
            "old processed output retained"
        )

    with dedup_tmp.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as f:
        writer = csv.DictWriter(
            f,
            fieldnames=COLUMNS + DEDUP_EXTRA,
        )

        writer.writeheader()

        for rec, refs in seen.values():
            writer.writerow({
                **rec,

                "drugcentral_supporting_record_count":
                    len(refs),

                "drugcentral_supporting_source_rows":
                    "|".join(refs),
            })

    if sum(
        len(v[1])
        for v in seen.values()
    ) != counts["retained_rows"]:
        raise RuntimeError(
            "Deduplication lost source support"
        )

    print(
        "[INPUT ROWS]",
        f"{counts['input_rows']:,}",
    )

    print(
        "[RETAINED / REJECTED]",
        f"{counts['retained_rows']:,}",
        f"/ {counts['rejected_rows']:,}",
    )

    print(
        "[DEDUP VIEW ASSERTIONS]",
        f"{len(seen):,}",
    )

    print(
        "[RELATIONSHIP COUNTS]",
        dict(relations),
    )

    print(
        "[RAW RELATIONSHIPS]",
        dict(raw_rels),
    )

    print(
        "[DISEASE IDENTIFIER SOURCE]",
        dict(id_sources),
    )

    print(
        "[DRUG NAME LOOKUP]",
        dict(drug_names),
    )

    print(
        "[DISEASE NAME LOOKUP]",
        dict(diseases),
    )

    print(
        "[IDENTIFIER SYNTAX]",
        dict(validations),
    )

    print(
        "[OUTPUT COLUMNS]",
        " | ".join(COLUMNS),
    )

    print(
        "[REAL OUTPUT EXAMPLE]",
        json.dumps(
            example,
            ensure_ascii=False,
            indent=2,
        )[:2200],
    )

    return {
        "counts": dict(counts),

        "deduplicated_assertions":
            len(seen),

        "distinct_drug_disease_relationship_keys":
            len(unique_pairs),

        "relationship_type_counts":
            dict(relations),

        "source_relationship_raw_counts":
            dict(raw_rels),

        "disease_identifier_source_column_counts":
            dict(id_sources),

        "drug_name_lookup_counts":
            dict(drug_names),

        "disease_name_lookup_counts":
            dict(diseases),

        "disease_identifier_validation_counts":
            dict(validations),

        "output_header": COLUMNS,

        "first_output_record": example,
    }


def main():
    ap = argparse.ArgumentParser(
        description=__doc__
    )

    ap.add_argument(
        "--raw-dir",
        type=Path,
        default=RAW_DIR,
    )

    ap.add_argument(
        "--out",
        type=Path,
        default=(
            OUT_DIR /
            "drugcentral_drug_disease.csv"
        ),
    )

    ap.add_argument(
        "--metadata",
        type=Path,
        default=(
            META_DIR /
            "drugcentral_processing_metadata.json"
        ),
    )

    ap.add_argument(
        "--offline",
        action="store_true",
        help="Only reuse verified cached raw snapshot",
    )

    ap.add_argument(
        "--refresh-source",
        action="store_true",
        help=(
            "Explicitly replace frozen source "
            "from live PostgreSQL"
        ),
    )

    ap.add_argument(
        "--host",
        default=os.getenv(
            "DRUGCENTRAL_HOST",
            "unmtid-dbs.net",
        ),
    )

    ap.add_argument(
        "--port",
        type=int,
        default=int(
            os.getenv(
                "DRUGCENTRAL_PORT",
                "5433",
            )
        ),
    )

    ap.add_argument(
        "--dbname",
        default=os.getenv(
            "DRUGCENTRAL_DBNAME",
            "drugcentral",
        ),
    )

    ap.add_argument(
        "--user",
        default=os.getenv(
            "DRUGCENTRAL_USER",
            "drugman",
        ),
    )

    ap.add_argument(
        "--password",
        default=os.getenv(
            "DRUGCENTRAL_PASSWORD",
            "dosage",
        ),
    )

    ap.add_argument(
        "--fetch-size",
        type=int,
        default=20000,
    )

    args = ap.parse_args()

    if (
        args.offline
        and args.refresh_source
    ):
        ap.error(
            "--offline and --refresh-source "
            "are mutually exclusive"
        )

    if args.fetch_size < 1:
        ap.error(
            "--fetch-size must be positive"
        )

    raw_dir = args.raw_dir.resolve()
    out = args.out.resolve()
    meta = args.metadata.resolve()

    snapshot = raw_dir / RAW_NAME
    manifest_path = raw_dir / MANIFEST_NAME

    dedup = out.with_name(
        "drugcentral_drug_disease_deduplicated.csv"
    )

    rejects = out.with_name(
        "drugcentral_rejected_rows.csv"
    )

    for p in (
        snapshot,
        manifest_path,
        out,
        dedup,
        rejects,
        meta,
    ):
        p.parent.mkdir(
            parents=True,
            exist_ok=True,
        )

    heading(
        "1. RAW FILE INVENTORY AND REUSE RULE"
    )

    for path in sorted(
        raw_dir.rglob("*")
    ):
        if path.is_file():
            print(
                "  ",
                path.relative_to(raw_dir),
                f"({path.stat().st_size:,} bytes)",
            )

    temp_paths = []
    staged = {}

    try:
        if (
            snapshot.is_file()
            and manifest_path.is_file()
            and not args.refresh_source
        ):
            manifest = json.loads(
                manifest_path.read_text(
                    encoding="utf-8"
                )
            )

            if (
                manifest.get("raw_snapshot_sha256")
                != digest(snapshot)
                or manifest.get("raw_snapshot_columns")
                != RAW_COLUMNS
            ):
                raise RuntimeError(
                    "Frozen snapshot/manifest hash "
                    "or schema mismatch; do not merge"
                )

            source = snapshot

            print(
                "[REUSING VERIFIED RAW SNAPSHOT]",
                source,
            )

        else:
            if args.offline:
                raise FileNotFoundError(
                    "Offline mode needs both verified "
                    "snapshot CSV and manifest JSON"
                )

            if (
                snapshot.exists()
                and not args.refresh_source
            ):
                raise RuntimeError(
                    "Raw CSV exists without a verified "
                    "manifest; use --refresh-source explicitly"
                )

            source = snapshot.with_name(
                snapshot.name + ".processing.tmp"
            )

            temp_paths.append(source)

            manifest = fetch_snapshot(
                args,
                source,
            )

            manifest_tmp = stage_json(
                manifest_path,
                manifest,
            )

            staged[snapshot] = source
            staged[manifest_path] = manifest_tmp

            temp_paths.append(manifest_tmp)

        expected = int(
            manifest["raw_snapshot_rows"]
        )

        output_tmp = out.with_name(
            out.name + ".processing.tmp"
        )

        dedup_tmp = dedup.with_name(
            dedup.name + ".processing.tmp"
        )

        rejects_tmp = rejects.with_name(
            rejects.name + ".processing.tmp"
        )

        for p in (
            output_tmp,
            dedup_tmp,
            rejects_tmp,
        ):
            temp_paths.append(p)

        result = process_source(
            source,
            expected,
            output_tmp,
            dedup_tmp,
            rejects_tmp,
        )

        staged[out] = output_tmp
        staged[dedup] = dedup_tmp
        staged[rejects] = rejects_tmp

        report = {
            "status": "success",

            "finished_at_utc": now(),

            "source": "DrugCentral",

            "source_table":
                "public.omop_relationship",

            "source_page":
                "https://drugcentral.org/download",

            "raw_snapshot": str(snapshot),

            "raw_snapshot_metadata":
                str(manifest_path),

            "raw_snapshot_sha256":
                digest(source),

            "database_release":
                manifest.get("dbversion", []),

            "source_columns":
                manifest.get("source_table_header", []),

            "raw_columns": RAW_COLUMNS,

            "source_name_join":
                manifest.get("name_join", ""),

            "source_sql":
                manifest.get("query", ""),

            "source_row_accounting": result,

            "deduplication": (
                "Full source evidence retained; "
                "separate exact-assertion view "
                "with supporting raw row numbers"
            ),

            "identifiers": (
                "DrugCentral structure ID; condition "
                "primary ID column-specific "
                "SNOMEDCT > UMLS > OMOP; "
                "all IDs retained"
            ),

            "caveat": (
                "Clinical annotations contain "
                "contraindications, indications and "
                "off-label uses; overlapping upstream "
                "evidence is not independent validation"
            ),

            "output_files": {
                "full": str(out),

                "deduplicated": str(dedup),

                "rejections": str(rejects),
            },

            "output_sha256": {
                "full": digest(output_tmp),

                "deduplicated": digest(dedup_tmp),

                "rejections": digest(rejects_tmp),
            },
        }

        metadata_tmp = stage_json(
            meta,
            report,
        )

        temp_paths.append(metadata_tmp)
        staged[meta] = metadata_tmp

        heading(
            "OVERWRITE ORIGINAL PATHS "
            "ONLY AFTER VALIDATION"
        )

        for destination, tmp in staged.items():
            tmp.replace(destination)

            print(
                "[WRITTEN]",
                destination,
            )

        print(
            "SUCCESS: Merge the FULL "
            "drugcentral_drug_disease.csv, "
            "not the deduplicated view"
        )

    finally:
        for p in temp_paths:
            p.unlink(missing_ok=True)


if __name__ == "__main__":
    main()
