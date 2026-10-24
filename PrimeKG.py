
#!/usr/bin/env python3
"""PrimeKG v2.1 kg.csv -> source-preserving drug--disease CSVs.

Read a pinned local kg.csv by default. Stream the raw CSV into a temporary,
indexed SQLite store so bidirectional graph edges can be reconciled without
loading the ~1 GB graph into memory. A reverse edge with a forward counterpart
is recorded as provenance, not counted as independent evidence.

No download unless --download-missing. Existing outputs are replaced only
following source parsing, count reconciliation, and staged CSV validation.
These knowledge-graph relations do not independently establish approval.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
import shutil
import sqlite3
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(os.environ.get(
    "DRUG_DISEASE_ROOT",
    Path(__file__).resolve().parent
)).resolve()

RAW_DIR = ROOT / "data/raw/primekg"
OUTPUT_DIR = ROOT / "data/processed/primekg"

DATA_URL = (
    "https://dataverse.harvard.edu/api/access/datafile/6180620"
)

REQUIRED = [
    "relation",
    "display_relation",
    "x_index",
    "x_id",
    "x_type",
    "x_name",
    "x_source",
    "y_index",
    "y_id",
    "y_type",
    "y_name",
    "y_source",
]

NODE_FIELDS = ("index", "id", "type", "name", "source")

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
    "primekg_relation",
    "primekg_display_relation",
    "source_node_type",
    "target_node_type",
    "evidence_text",
    "primekg_x_index",
    "primekg_x_id",
    "primekg_x_type",
    "primekg_x_source",
    "primekg_y_index",
    "primekg_y_id",
    "primekg_y_type",
    "primekg_y_source",
    "orientation_in_source",
    "mondo_grouped_member_ids_unverified",
    "primekg_source_row_number",
    "primekg_source_line_number",
    "primekg_source_record_json",
    "primekg_original_x_id",
    "primekg_original_y_id",
    "primekg_original_x_source",
    "primekg_original_y_source",
    "primekg_reverse_support_json",
    "primekg_reverse_support_count",
    "primekg_name_resolution_drug",
    "primekg_name_resolution_disease",
    "primekg_evidence_granularity",
]

for entity in ("drug", "disease"):
    EXTRA.extend(
        f"{entity}_identifier_{s}"
        for s in (
            "original",
            "namespace",
            "canonical",
            "validation",
        )
    )

COLUMNS = CORE + EXTRA

DEDUP_EXTRA = [
    "primekg_supporting_record_count",
    "primekg_supporting_source_rows",
]

REJECT_COLS = [
    "source_row_number",
    "physical_line_number",
    "reason",
    "source_record_json",
]

RELATIONS = {
    "indication": (
        "indication",
        "primekg_indication",
    ),
    "contraindication": (
        "contraindication",
        "primekg_contraindication",
    ),
    "off label use": (
        "off_label_use",
        "primekg_off_label_use",
    ),
}


def stamp():
    return datetime.now(timezone.utc).isoformat()


def clean(value):
    return "" if value is None else str(value).strip()


def title(text):
    print(
        "\n" + "=" * 88 +
        "\n" + text +
        "\n" + "=" * 88,
        flush=True,
    )


def checksum(path, algorithm="sha256"):
    h = hashlib.new(algorithm)

    with Path(path).open("rb") as handle:
        for block in iter(
            lambda: handle.read(4 * 1024 * 1024),
            b"",
        ):
            h.update(block)

    return h.hexdigest()


def stable_json(x):
    return json.dumps(
        x,
        ensure_ascii=False,
        separators=(",", ":"),
        sort_keys=True,
    )


def canonical_rel(raw):
    return re.sub(
        r"[\s_\-/]+",
        " ",
        clean(raw).casefold(),
    ).strip()


def id_annotation(raw, source, entity):
    """Validate identifier syntax without inferring ontology equivalence."""

    original = clean(raw)
    namespace = clean(source)

    key = re.sub(
        r"[\s_\-]+",
        "",
        namespace,
    ).casefold()

    token = original.upper()

    if not original:
        return (
            original,
            namespace or "unresolved",
            "",
            "missing",
        )

    if key == "drugbank":
        good = bool(
            re.fullmatch(r"DB\d{5,}", token)
        )

        return (
            original,
            "DrugBank",
            "DRUGBANK:" + token if good else "",
            "syntax_valid" if good else "invalid_syntax",
        )

    if key == "chembl":
        good = bool(
            re.fullmatch(
                r"(?:CHEMBL:)?CHEMBL\d+",
                token,
            )
        )

        return (
            original,
            "ChEMBL",
            "CHEMBL:" + token.removeprefix("CHEMBL:")
            if good else "",
            "syntax_valid" if good else "invalid_syntax",
        )

    if key == "mondogrouped":
        return (
            original,
            "MONDO_grouped",
            "",
            "grouped_concept_unresolved",
        )

    if key == "mondo":
        match = re.fullmatch(
            r"(?:MONDO[:_])?(\d{1,7})",
            token,
        )

        return (
            original,
            "MONDO",
            "MONDO:" + match[1].zfill(7)
            if match else "",
            "syntax_valid_zero_padded"
            if match else "invalid_syntax",
        )

    if key in {"doid", "diseaseontology"}:
        match = re.fullmatch(
            r"(?:DOID[:_])?(\d+)",
            token,
        )

        return (
            original,
            "DOID",
            "DOID:" + match[1] if match else "",
            "syntax_valid" if match else "invalid_syntax",
        )

    if key in {"mesh", "medicalsubjectheadings"}:
        match = re.fullmatch(
            r"(?:MESH[:_])?([CD]\d{6,9})",
            token,
        )

        return (
            original,
            "MeSH",
            "MESH:" + match[1] if match else "",
            "syntax_valid" if match else "invalid_syntax",
        )

    if key in {"umls", "umlsconcept"}:
        match = re.fullmatch(
            r"(?:UMLS:)?(C\d{7})",
            token,
        )

        return (
            original,
            "UMLS",
            "UMLS:" + match[1] if match else "",
            "syntax_valid" if match else "invalid_syntax",
        )

    if key in {"omim", "mim"}:
        match = re.fullmatch(
            r"(?:OMIM:)?(\d{6})",
            token,
        )

        return (
            original,
            "OMIM",
            "OMIM:" + match[1] if match else "",
            "syntax_valid" if match else "invalid_syntax",
        )

    if key == "efo":
        match = re.fullmatch(
            r"(?:EFO[:_])?(\d{1,7})",
            token,
        )

        return (
            original,
            "EFO",
            "EFO:" + match[1].zfill(7)
            if match else "",
            "syntax_valid" if match else "invalid_syntax",
        )

    if key in {"hpo", "hp"}:
        match = re.fullmatch(
            r"(?:HP[:_])?(\d{1,7})",
            token,
        )

        return (
            original,
            "HPO",
            "HP:" + match[1].zfill(7)
            if match else "",
            "syntax_valid" if match else "invalid_syntax",
        )

    return (
        original,
        namespace or "unresolved",
        "",
        "unrecognised_namespace",
    )


def grouped_mondo_members(raw):
    pieces = clean(raw).split("_")

    if (
        not pieces
        or not all(
            re.fullmatch(r"\d{1,7}", x)
            for x in pieces
        )
    ):
        return ""

    return "|".join(
        "MONDO:" + x.zfill(7)
        for x in pieces
    )


def source_csv_header_ok(path):
    if (
        not path.is_file()
        or path.stat().st_size == 0
    ):
        return False

    try:
        with path.open(
            "r",
            encoding="utf-8-sig",
            newline="",
        ) as f:
            header = next(csv.reader(f))

        return (
            len(header) == len(set(header))
            and set(REQUIRED).issubset(header)
        )

    except (
        OSError,
        UnicodeError,
        csv.Error,
        StopIteration,
    ):
        return False


def obtain_source(raw, download_missing):
    raw.mkdir(
        parents=True,
        exist_ok=True,
    )

    target = raw / "kg.csv"

    if target.exists():
        if not source_csv_header_ok(target):
            raise RuntimeError(
                "Local kg.csv is not the required "
                f"PrimeKG CSV: {target}"
            )

        print("[REUSE]", target)
        return target

    manual = [
        raw / "6180620",
        ROOT / "6180620",
        Path.cwd() / "6180620",
    ]

    for source in manual:
        if (
            source.is_file()
            and source_csv_header_ok(source)
        ):
            if source.resolve() != target.resolve():
                stage = target.with_name(
                    target.name + ".copy.tmp"
                )

                try:
                    shutil.copyfile(
                        source,
                        stage,
                    )

                    if not source_csv_header_ok(stage):
                        raise RuntimeError(
                            "Manual source did not "
                            "pass header validation"
                        )

                    stage.replace(target)

                finally:
                    stage.unlink(missing_ok=True)

            print(
                "[COPIED EXISTING MANUAL DATAVERSE FILE]",
                source,
            )

            return target

    if not download_missing:
        raise FileNotFoundError(
            f"Missing {target} and valid manually "
            "downloaded 6180620; use "
            "--download-missing if you want to fetch it"
        )

    try:
        import requests

    except ImportError as exc:
        raise RuntimeError(
            "Install requests only if downloading: "
            "pip install requests"
        ) from exc

    temporary = target.with_name(
        "kg.csv.download.tmp"
    )

    try:
        print("[DOWNLOAD]", DATA_URL)

        with requests.get(
            DATA_URL,
            stream=True,
            allow_redirects=True,
            timeout=(30, 600),
        ) as resp:
            resp.raise_for_status()

            with temporary.open("wb") as f:
                for block in resp.iter_content(
                    4 * 1024 * 1024
                ):
                    if block:
                        f.write(block)

        if not source_csv_header_ok(temporary):
            raise RuntimeError(
                "Dataverse response was not "
                "a PrimeKG kg.csv"
            )

        temporary.replace(target)

    finally:
        temporary.unlink(missing_ok=True)

    return target


def oriented_key(row, direction, relation):
    if direction == "forward":
        drug, disease = "x", "y"
    else:
        drug, disease = "y", "x"

    return stable_json([
        relation,
        clean(row[drug + "_index"]),
        clean(row[drug + "_id"]),
        clean(row[drug + "_source"]),
        clean(row[disease + "_index"]),
        clean(row[disease + "_id"]),
        clean(row[disease + "_source"]),
    ])


def create_database(path):
    con = sqlite3.connect(path)

    con.execute("PRAGMA journal_mode=DELETE")
    con.execute("PRAGMA synchronous=NORMAL")
    con.execute("PRAGMA temp_store=FILE")

    con.execute("""
        CREATE TABLE candidate (
            source_row INTEGER PRIMARY KEY,
            physical_line INTEGER NOT NULL,
            pair_key TEXT NOT NULL,
            direction TEXT NOT NULL,
            relation_key TEXT NOT NULL,
            raw_json TEXT NOT NULL
        )
    """)

    con.execute("""
        CREATE TABLE assertion_first (
            assertion_key TEXT PRIMARY KEY,
            first_json TEXT NOT NULL
        )
    """)

    con.execute("""
        CREATE TABLE assertion_support (
            assertion_key TEXT NOT NULL,
            source_row INTEGER NOT NULL
        )
    """)

    return con


def index_raw(
    path,
    conn,
    reject_writer,
    batch_size,
):
    title(
        "2. STREAM RAW GRAPH, VALIDATE COLUMNS, "
        "INDEX SELECTED EDGES"
    )

    all_relations = Counter()
    counters = Counter()
    buffer = []
    first = None

    with path.open(
        "r",
        encoding="utf-8-sig",
        newline="",
    ) as f:

        reader = csv.reader(
            f,
            strict=True,
        )

        header = next(reader)

        if (
            len(header) != len(set(header))
            or not set(REQUIRED).issubset(header)
        ):
            raise RuntimeError(
                "kg.csv missing expected columns "
                "or contains repeated headers: "
                + str(header)
            )

        print(
            "Original kg.csv header:",
            " | ".join(header),
        )

        for source_row, cells in enumerate(
            reader,
            1,
        ):

            if len(cells) != len(header):
                raise RuntimeError(
                    f"Malformed CSV at data row {source_row}, "
                    f"file line {reader.line_num}: "
                    f"{len(cells)} cells != "
                    f"{len(header)} columns"
                )

            row = dict(
                zip(header, cells)
            )

            if first is None:
                first = row

            counters["graph_rows"] += 1

            relation_raw = clean(
                row["relation"]
            )

            all_relations[
                relation_raw or "[blank]"
            ] += 1

            rel = canonical_rel(
                relation_raw
            )

            if rel not in RELATIONS:
                continue

            counters[
                "selected_relation_rows"
            ] += 1

            xtype = clean(
                row["x_type"]
            ).casefold()

            ytype = clean(
                row["y_type"]
            ).casefold()

            if xtype == "drug" and ytype == "disease":
                direction = "forward"

            elif (
                xtype == "disease"
                and ytype == "drug"
            ):
                direction = "reverse"

            else:
                counters[
                    "excluded_wrong_node_types"
                ] += 1

                reject_writer.writerow({
                    "source_row_number": source_row,
                    "physical_line_number":
                        reader.line_num,
                    "reason": (
                        "selected_relationship_but_"
                        "not_drug_disease_nodes"
                    ),
                    "source_record_json":
                        stable_json(row),
                })

                continue

            counters[
                "indexed_" + direction + "_rows"
            ] += 1

            buffer.append((
                source_row,
                reader.line_num,
                oriented_key(
                    row,
                    direction,
                    rel,
                ),
                direction,
                rel,
                stable_json(row),
            ))

            if len(buffer) >= batch_size:
                conn.executemany(
                    "INSERT INTO candidate "
                    "VALUES (?,?,?,?,?,?)",
                    buffer,
                )

                conn.commit()
                buffer.clear()

            if source_row % 1000000 == 0:
                print(
                    "[SCANNED]",
                    f"{source_row:,}",
                    flush=True,
                )

    if buffer:
        conn.executemany(
            "INSERT INTO candidate "
            "VALUES (?,?,?,?,?,?)",
            buffer,
        )

        conn.commit()

    conn.execute(
        "CREATE INDEX idx_candidate_pair_direction "
        "ON candidate(pair_key, direction)"
    )

    conn.commit()

    print(
        "Graph rows:",
        f"{counters['graph_rows']:,}",
    )

    print(
        "Selected relation rows:",
        f"{counters['selected_relation_rows']:,}",
    )

    print(
        "Indexed drug->disease:",
        f"{counters['indexed_forward_rows']:,}",
    )

    print(
        "Indexed disease->drug:",
        f"{counters['indexed_reverse_rows']:,}",
    )

    print(
        "Wrong node types rejected:",
        f"{counters['excluded_wrong_node_types']:,}",
    )

    print(
        "First source row:",
        stable_json(first)[:1100],
    )

    return (
        header,
        dict(all_relations),
        counters,
        first,
    )


def make_record(
    raw,
    source_row,
    physical_line,
    direction,
    reverse_rows,
):
    d, c = (
        ("x", "y")
        if direction == "forward"
        else ("y", "x")
    )

    drug_id = clean(
        raw[d + "_id"]
    )

    disease_id = clean(
        raw[c + "_id"]
    )

    missing = []

    if not drug_id:
        missing.append("missing_drug_id")

    if not disease_id:
        missing.append("missing_disease_id")

    if missing:
        return None, "|".join(missing)

    drug_source = clean(
        raw[d + "_source"]
    )

    disease_source = clean(
        raw[c + "_source"]
    )

    drug_name = clean(
        raw[d + "_name"]
    )

    disease_name = clean(
        raw[c + "_name"]
    )

    drug_name_status = (
        "provided"
        if drug_name
        else "identifier_placeholder"
    )

    disease_name_status = (
        "provided"
        if disease_name
        else "identifier_placeholder"
    )

    if not drug_name:
        drug_name = drug_id

    if not disease_name:
        disease_name = disease_id

    rel_raw = clean(
        raw["relation"]
    )

    normalized_relation = canonical_rel(
        rel_raw
    )

    relationship_type, evidence_type = (
        RELATIONS[normalized_relation]
    )

    dmeta = id_annotation(
        drug_id,
        drug_source,
        "drug",
    )

    cmeta = id_annotation(
        disease_id,
        disease_source,
        "disease",
    )

    if direction == "reverse":
        orientation = (
            "disease_to_drug_only_flipped"
        )

    elif reverse_rows:
        orientation = (
            "drug_to_disease_with_reverse"
        )

    else:
        orientation = (
            "drug_to_disease_only"
        )

    row = {
        "drug_name": drug_name,
        "drug_identifier": drug_id,
        "drug_identifier_type": dmeta[1],

        "disease_or_condition_name":
            disease_name,

        "disease_or_condition_identifier":
            disease_id,

        "disease_or_condition_identifier_type":
            cmeta[1],

        "relationship_type": relationship_type,
        "evidence_type": evidence_type,
        "source": "PrimeKG",
        "internal_source": "PrimeKG v2.1 kg.csv",

        "primekg_relation": rel_raw,

        "primekg_display_relation":
            clean(raw["display_relation"]),

        "source_node_type": "drug",
        "target_node_type": "disease",

        "evidence_text": (
            f"PrimeKG relation={rel_raw}; "
            f"display_relation="
            f"{clean(raw['display_relation'])}; "
            f"orientation={orientation}; "
            "not independent approval evidence"
        ),

        "primekg_x_index":
            clean(raw[d + "_index"]),

        "primekg_x_id": drug_id,
        "primekg_x_type": "drug",
        "primekg_x_source": drug_source,

        "primekg_y_index":
            clean(raw[c + "_index"]),

        "primekg_y_id": disease_id,
        "primekg_y_type": "disease",
        "primekg_y_source": disease_source,

        "orientation_in_source": orientation,

        "mondo_grouped_member_ids_unverified": (
            grouped_mondo_members(disease_id)
            if disease_source.casefold()
            == "mondo_grouped"
            else ""
        ),

        "primekg_source_row_number":
            str(source_row),

        "primekg_source_line_number":
            str(physical_line),

        "primekg_source_record_json":
            stable_json(raw),

        "primekg_original_x_id":
            clean(raw["x_id"]),

        "primekg_original_y_id":
            clean(raw["y_id"]),

        "primekg_original_x_source":
            clean(raw["x_source"]),

        "primekg_original_y_source":
            clean(raw["y_source"]),

        "primekg_reverse_support_json":
            stable_json(reverse_rows),

        "primekg_reverse_support_count":
            str(len(reverse_rows)),

        "primekg_name_resolution_drug":
            drug_name_status,

        "primekg_name_resolution_disease":
            disease_name_status,

        "primekg_evidence_granularity":
            "oriented_primekg_source_graph_edge",
    }

    for entity, metadata in (
        ("drug", dmeta),
        ("disease", cmeta),
    ):
        for key, item in zip(
            (
                "original",
                "namespace",
                "canonical",
                "validation",
            ),
            metadata,
        ):
            row[
                f"{entity}_identifier_{key}"
            ] = item

    return row, ""


def write_processed(
    con,
    full_tmp,
    dedup_tmp,
    reject_writer,
):
    title(
        "3. ORIENT EDGES, PRESERVE "
        "REVERSE SUPPORT, WRITE EVIDENCE"
    )

    counts = Counter()
    relations = Counter()
    orientations = Counter()
    namespace_counts = Counter()
    validation_counts = Counter()

    examples = []

    canonical_mondo_raw = {}
    mondo_collisions = set()

    pair_cur = con.cursor()

    output_query = """
        SELECT
            c.source_row,
            c.physical_line,
            c.pair_key,
            c.direction,
            c.raw_json

        FROM candidate c

        WHERE c.direction = 'forward'
           OR (
               c.direction = 'reverse'
               AND NOT EXISTS (
                   SELECT 1
                   FROM candidate p
                   WHERE p.pair_key = c.pair_key
                   AND p.direction = 'forward'
               )
           )

        ORDER BY c.source_row
    """

    with full_tmp.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as f:

        writer = csv.DictWriter(
            f,
            fieldnames=COLUMNS,
        )

        writer.writeheader()

        for (
            source_row,
            physical_line,
            pair_key,
            direction,
            original,
        ) in con.execute(output_query):

            raw = json.loads(original)
            reverse_rows = []

            if direction == "forward":
                for rev_row, rev_json in pair_cur.execute(
                    """
                    SELECT source_row, raw_json
                    FROM candidate
                    WHERE pair_key = ?
                    AND direction = 'reverse'
                    ORDER BY source_row
                    """,
                    (pair_key,),
                ):
                    reverse_rows.append({
                        "source_row_number": rev_row,
                        "source_record":
                            json.loads(rev_json),
                    })

            rec, reject_reason = make_record(
                raw,
                source_row,
                physical_line,
                direction,
                reverse_rows,
            )

            if reject_reason:
                counts[
                    "identity_rejected_rows"
                ] += 1

                reject_writer.writerow({
                    "source_row_number": source_row,
                    "physical_line_number":
                        physical_line,
                    "reason": reject_reason,
                    "source_record_json": original,
                })

                continue

            writer.writerow(rec)
            counts["output_rows"] += 1

            relations[
                rec["relationship_type"]
            ] += 1

            orientations[
                rec["orientation_in_source"]
            ] += 1

            namespace_counts[
                rec["disease_identifier_namespace"]
            ] += 1

            validation_counts[
                rec["disease_identifier_validation"]
            ] += 1

            if len(examples) < 2:
                examples.append(rec)

            if (
                rec["disease_identifier_namespace"]
                == "MONDO"
            ):
                canon = rec[
                    "disease_identifier_canonical"
                ]

                original_id = rec[
                    "disease_identifier_original"
                ]

                if canon:
                    if (
                        canon in canonical_mondo_raw
                        and canonical_mondo_raw[canon]
                        != original_id
                    ):
                        mondo_collisions.add(canon)

                    else:
                        canonical_mondo_raw[
                            canon
                        ] = original_id

            # Keep the main file complete.
            # Deduplication only builds a
            # separate assertion-level view.

            fields = (
                "drug_identifier",
                "drug_identifier_type",
                "disease_or_condition_identifier",
                "disease_or_condition_identifier_type",
                "relationship_type",
                "evidence_type",
                "primekg_display_relation",
            )

            key = stable_json([
                rec[name]
                for name in fields
            ])

            con.execute(
                """
                INSERT OR IGNORE
                INTO assertion_first
                VALUES (?, ?)
                """,
                (
                    key,
                    stable_json(rec),
                ),
            )

            con.execute(
                """
                INSERT INTO assertion_support
                VALUES (?, ?)
                """,
                (
                    key,
                    source_row,
                ),
            )

            if counts["output_rows"] % 10000 == 0:
                con.commit()

                print(
                    "[EVIDENCE WRITTEN]",
                    f"{counts['output_rows']:,}",
                    flush=True,
                )

    con.commit()

    if not counts["output_rows"]:
        raise RuntimeError(
            "No valid drug-disease associations: "
            "outputs were not replaced"
        )

    con.execute(
        "CREATE INDEX idx_assertion_support "
        "ON assertion_support("
        "assertion_key, source_row)"
    )

    con.commit()

    assertions = 0
    support_sum = 0

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

        for key, first_json in con.execute(
            """
            SELECT assertion_key, first_json
            FROM assertion_first
            ORDER BY rowid
            """
        ):

            supports = [
                str(r[0])
                for r in pair_cur.execute(
                    """
                    SELECT source_row
                    FROM assertion_support
                    WHERE assertion_key = ?
                    ORDER BY source_row
                    """,
                    (key,),
                )
            ]

            row = json.loads(first_json)

            row[
                "primekg_supporting_record_count"
            ] = len(supports)

            row[
                "primekg_supporting_source_rows"
            ] = "|".join(supports)

            writer.writerow(row)

            assertions += 1
            support_sum += len(supports)

    if support_sum != counts["output_rows"]:
        raise RuntimeError(
            "Dedup support mismatch: "
            f"{support_sum} != "
            f"{counts['output_rows']}"
        )

    print(
        "Output rows:",
        f"{counts['output_rows']:,}",
    )

    print(
        "Dedup assertions:",
        f"{assertions:,}",
    )

    print(
        "Relationship counts:",
        dict(relations),
    )

    print(
        "Source orientations:",
        dict(orientations),
    )

    print(
        "Disease namespaces:",
        dict(namespace_counts),
    )

    print(
        "Identifier checks:",
        dict(validation_counts),
    )

    print(
        "MONDO normalisation collisions:",
        len(mondo_collisions),
    )

    print(
        "Output schema:",
        " | ".join(COLUMNS),
    )

    print(
        "First real output record:",
        json.dumps(
            examples[0],
            ensure_ascii=False,
            indent=2,
        )[:3200],
    )

    return {
        "counts": dict(counts),
        "deduplicated_assertions": assertions,

        "relationship_counts":
            dict(relations),

        "orientation_counts":
            dict(orientations),

        "disease_namespace_counts":
            dict(namespace_counts),

        "disease_identifier_validation_counts":
            dict(validation_counts),

        "mondo_canonical_collision_count":
            len(mondo_collisions),

        "mondo_canonical_collision_examples":
            sorted(mondo_collisions)[:20],

        "output_columns": COLUMNS,

        "examples": examples,
    }


def validate_csv_header(path, expected):
    with path.open(
        "r",
        encoding="utf-8",
        newline="",
    ) as f:
        header = next(csv.reader(f))

    if header != expected:
        raise RuntimeError(
            f"Staged CSV header mismatch for "
            f"{path}: {header}"
        )


def main():
    parser = argparse.ArgumentParser(
        description=__doc__
    )

    parser.add_argument(
        "--raw-dir",
        type=Path,
        default=RAW_DIR,
    )

    parser.add_argument(
        "--output-dir",
        type=Path,
        default=OUTPUT_DIR,
    )

    parser.add_argument(
        "--download-missing",
        action="store_true",
    )

    parser.add_argument(
        "--batch-size",
        type=int,
        default=10000,
    )

    parser.add_argument(
        "--expected-md5",
        default="",
        help=(
            "Optional checksum from pinned "
            "Dataverse file metadata"
        ),
    )

    args = parser.parse_args()

    if args.batch_size < 1:
        parser.error(
            "--batch-size must be >= 1"
        )

    raw_dir = (
        args.raw_dir.expanduser().resolve()
    )

    out_dir = (
        args.output_dir.expanduser().resolve()
    )

    kg = obtain_source(
        raw_dir,
        args.download_missing,
    )

    out_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    title(
        "1. PINNED SOURCE AND "
        "LOCAL RAW INVENTORY"
    )

    inventory = []

    for path in sorted(raw_dir.rglob("*")):
        if path.is_file():
            record = {
                "relative_path": str(
                    path.relative_to(raw_dir)
                ),
                "bytes": path.stat().st_size,
                "input": (
                    path.resolve()
                    == kg.resolve()
                ),
            }

            inventory.append(record)

            print(
                "[INPUT]"
                if record["input"]
                else "[OTHER]",
                record,
            )

    raw_sha256 = checksum(kg)
    raw_md5 = checksum(kg, "md5")

    print(
        "[SOURCE]",
        kg,
        "bytes",
        f"{kg.stat().st_size:,}",
    )

    print(
        "[SHA256]",
        raw_sha256,
    )

    print(
        "[MD5]",
        raw_md5,
    )

    if (
        args.expected_md5
        and raw_md5.casefold()
        != args.expected_md5.strip().casefold()
    ):
        raise RuntimeError(
            "Raw dataset MD5 differs from "
            "requested pinned snapshot"
        )

    full = (
        out_dir / "primekg_drug_disease.csv"
    )

    dedup = (
        out_dir /
        "primekg_drug_disease_deduplicated.csv"
    )

    rejects = (
        out_dir / "primekg_rejected_rows.csv"
    )

    summary = (
        out_dir / "primekg_summary.json"
    )

    dbpath = (
        out_dir /
        ".primekg_source_index.processing.tmp.sqlite3"
    )

    destinations = (
        full,
        dedup,
        rejects,
        summary,
    )

    stages = {
        p: p.with_name(
            p.name + ".processing.tmp"
        )
        for p in destinations
    }

    con = None

    try:
        if dbpath.exists():
            dbpath.unlink()

        con = create_database(dbpath)

        with stages[rejects].open(
            "w",
            encoding="utf-8",
            newline="",
        ) as bad_f:

            bad_writer = csv.DictWriter(
                bad_f,
                fieldnames=REJECT_COLS,
            )

            bad_writer.writeheader()

            (
                source_header,
                relation_counts,
                input_counts,
                original_sample,
            ) = index_raw(
                kg,
                con,
                bad_writer,
                args.batch_size,
            )

            if not input_counts[
                "selected_relation_rows"
            ]:
                raise RuntimeError(
                    "No PrimeKG drug--disease "
                    "relationship labels found"
                )

            result = write_processed(
                con,
                stages[full],
                stages[dedup],
                bad_writer,
            )

        title(
            "4. ORIENTATION AND "
            "RECORD CONSERVATION"
        )

        counts = input_counts

        reverse_only = (
            result["orientation_counts"].get(
                "disease_to_drug_only_flipped",
                0,
            )
        )

        mirrored = (
            counts["indexed_reverse_rows"]
            - reverse_only
        )

        excluded = (
            counts["excluded_wrong_node_types"]
            + result["counts"].get(
                "identity_rejected_rows",
                0,
            )
        )

        observed = (
            result["counts"]["output_rows"]
            + mirrored
            + excluded
        )

        print(
            "Selected raw relation rows:",
            f"{counts['selected_relation_rows']:,}",
        )

        print(
            "Processed source events:",
            f"{result['counts']['output_rows']:,}",
        )

        print(
            "Mirrored reverse rows "
            "not double-counted:",
            f"{mirrored:,}",
        )

        print(
            "Rejected rows:",
            f"{excluded:,}",
        )

        print(
            "Reconciled:",
            observed == counts[
                "selected_relation_rows"
            ],
        )

        if (
            observed != counts[
                "selected_relation_rows"
            ]
            or mirrored < 0
        ):
            raise RuntimeError(
                "PrimeKG orientation / "
                "rejection accounting mismatch"
            )

        for path in (
            full,
            dedup,
            rejects,
        ):
            validate_csv_header(
                stages[path],
                COLUMNS + (
                    DEDUP_EXTRA
                    if path == dedup
                    else []
                )
                if path != rejects
                else REJECT_COLS,
            )

        report = {
            "status": "success",
            "timestamp_utc": stamp(),
            "source": "PrimeKG",

            "source_artifact": (
                "PrimeKG v2.1 kg.csv "
                "(Dataverse datafile 6180620)"
            ),

            "dataset_doi":
                "10.7910/DVN/IXA7BM",

            "publication_doi":
                "10.1038/s41597-023-01960-3",

            "source_download_url": DATA_URL,

            "raw_file_inventory": inventory,

            "raw_file_sha256": raw_sha256,

            "raw_file_md5": raw_md5,

            "raw_file_bytes":
                kg.stat().st_size,

            "original_source_header":
                source_header,

            "original_source_first_record":
                original_sample,

            "all_raw_relation_counts":
                relation_counts,

            "raw_scan_counts":
                dict(input_counts),

            "final_result": result,

            "selected_source_rows_accounted_for":
                observed,

            "reverse_rows_considered_mirrors":
                mirrored,

            "deduplication": (
                "Full output preserves each kept "
                "original oriented source edge; "
                "separate assertion view includes "
                "supporting source row numbers"
            ),

            "node_identity": (
                "Source node index, original ID "
                "and declared namespace; "
                "no cross-source equivalence inferred"
            ),

            "licensing_note": (
                "PrimeKG dataset license and "
                "upstream source terms must be "
                "reviewed separately for "
                "redistribution; code license "
                "is not evidence license"
            ),

            "interpretation": (
                "PrimeKG indication != verified "
                "regulatory approval; "
                "contraindication and off-label "
                "use are separate evidence categories"
            ),

            "output_files": {
                p.name: str(p)
                for p in destinations
            },

            "output_sha256": {
                p.name: checksum(stages[p])
                for p in (
                    full,
                    dedup,
                    rejects,
                )
            },
        }

        stages[summary].write_text(
            json.dumps(
                report,
                indent=2,
                ensure_ascii=False,
            ) + "\n",
            encoding="utf-8",
        )

        title(
            "5. REPLACING EXISTING "
            "PROCESSED FILES "
            "AFTER VALIDATION"
        )

        if con is not None:
            con.close()
            con = None

        for destination in destinations:
            stages[destination].replace(
                destination
            )

            print(
                "[WRITTEN]",
                destination,
            )

        print(
            "SUCCESS. Merge "
            "primekg_drug_disease.csv; "
            "keep relationship_type "
            "and evidence_type."
        )

    finally:
        if con is not None:
            con.close()

        for temp in stages.values():
            temp.unlink(
                missing_ok=True
            )

        dbpath.unlink(
            missing_ok=True
        )


if __name__ == "__main__":
    main()
