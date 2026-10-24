
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Analysis 07: locked source audit, schema validation, and name QC.

Read-only for all locked inputs. The merge manifest, not a directory glob or a
filename guess, defines the eleven included sources. Optional output tables go
to tables/. All scientific counts print as tab-separated tables on stdout.

Run:
    POLARS_MAX_THREADS=16 python Analysis_07_Source_Schema_Name_QC.py --save-tables
    POLARS_MAX_THREADS=16 python Analysis_07_Source_Schema_Name_QC.py --verify-hashes
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import polars as pl
from tqdm import tqdm

ROOT = Path(__file__).resolve().parent
CORE = [
    "drug_name", "drug_identifier", "drug_identifier_type",
    "disease_or_condition_name", "disease_or_condition_identifier",
    "disease_or_condition_identifier_type", "relationship_type",
    "evidence_type", "source", "internal_source",
]
SOURCES = {
    "AACT", "Broad Drug Repurposing Hub", "ChEMBL", "CTD", "DrugCentral",
    "Hetionet", "MEDI", "OpenTargets", "PrimeKG", "repoDB", "SIDER",
}
FINAL = {
    "events": "ALL_INCLUDED_SOURCES_drug_disease_merged.csv",
    "assertions": "ALL_INCLUDED_SOURCES_drug_disease_deduplicated.csv",
    "pairs": "ALL_INCLUDED_SOURCES_drug_disease_source_collapsed.csv",
}
EXPECTED_ROWS = {
    "events": 10_895_395,
    "assertions": 4_534_385,
    "pairs": 4_470_601,
}
EXPECTED_COLUMNS = {
    "events": 200,
    "assertions": 14,
    "pairs": 18,
}
EXPECTED_CONFLICTS = 87

PLACEHOLDERS = {
    "na", "n/a", "nan", "none", "null", "missing", "unknown",
    "not available", "not_applicable", "missing_drug_name",
    "missing_disease_or_condition_name",
}

# Identifier-shaped names are review flags, not automatic errors.
ID_SHAPED = (
    r"(?i)^(?:DB[0-9]{5}|CHEMBL[0-9]+|C[0-9]{7}|D[0-9]{6}|"
    r"(?:EFO|MONDO|HP|DOID|ORPHA)[_:][0-9]+)$"
)


def must(label, actual, expected):
    if actual != expected:
        raise RuntimeError(
            f"[FAIL] {label}: {actual!r} != {expected!r}"
        )

    desc = (
        f"{actual:,}"
        if isinstance(actual, int)
        else "OK"
    )
    print(f"[PASS] {label}: {desc}", flush=True)


def warn(msg, warnings):
    warnings.append(msg)
    print(f"[REVIEW] {msg}", flush=True)


def print_table(title, frame):
    print("\n" + "=" * 108)
    print(f"{title}  [{frame.height:,} rows]")
    print("=" * 108)

    if not frame.height:
        print("[EMPTY]")
        return

    print("\t".join(frame.columns))

    for row in frame.iter_rows():
        values = []
        for item in row:
            if item is None:
                val = ""
            elif isinstance(item, float):
                val = f"{item:.3f}"
            else:
                val = str(item)

            values.append(
                val.replace("\t", " ")
                   .replace("\n", " ")
                   .replace("\r", " ")
            )

        print("\t".join(values))


def csv_header(path):
    if not path.is_file() or path.stat().st_size == 0:
        raise RuntimeError(
            f"Missing or empty CSV: {path}"
        )

    with path.open(
        "r",
        encoding="utf-8-sig",
        newline="",
    ) as handle:
        names = next(
            csv.reader(handle, strict=True),
            [],
        )

    if not names or len(names) != len(set(names)):
        raise RuntimeError(
            f"Empty or duplicate CSV header: {path}"
        )

    return names


def scan(path, selected=None):
    q = pl.scan_csv(
        str(path),
        infer_schema=False,
        ignore_errors=False,
        truncate_ragged_lines=False,
        encoding="utf8",
    )
    return q.select(selected) if selected is not None else q


def string_key(c):
    return (
        pl.col(c)
        .fill_null("")
        .str.strip_chars()
        .str.replace_all(r"\s+", " ")
        .str.to_lowercase()
    )


def hash_file(path):
    h = hashlib.sha256()

    with path.open("rb") as f:
        for block in iter(
            lambda: f.read(8 * 1024 * 1024),
            b"",
        ):
            h.update(block)

    return h.hexdigest()


def source_file_path(project_root, processed, spec):
    """Use only exact paths recorded by the merger."""

    stated = spec.get("path")

    if not isinstance(stated, str) or not stated.strip():
        raise RuntimeError(
            f"No input path recorded for {spec.get('source')}"
        )

    raw = Path(stated).expanduser()

    if raw.is_absolute():
        candidates = [raw]
    else:
        candidates = [
            processed / raw,
            project_root / raw,
            project_root / "data" / "processed" / raw,
        ]

    for p in candidates:
        if p.is_file():
            return p.resolve()

    raise RuntimeError(
        f"Included source {spec.get('source')} not found. "
        f"Manifest path={stated!r}; "
        f"checked={[str(x) for x in candidates]}"
    )


def file_schema_audit(
    label,
    path,
    fields,
    n_columns,
    schema_from_manifest=None,
):
    header = csv_header(path)

    if header[:len(fields)] != fields:
        raise RuntimeError(
            f"{label}: required column order differs: "
            f"{header[:len(fields)]}"
        )

    if len(header) != n_columns:
        raise RuntimeError(
            f"{label}: found {len(header)} columns; "
            f"expected {n_columns}"
        )

    if (
        schema_from_manifest
        and header != schema_from_manifest
    ):
        raise RuntimeError(
            f"{label}: CSV header differs from merge manifest"
        )

    return {
        "file_role": label,
        "file": str(path),
        "columns": len(header),
        "required_prefix_columns": len(fields),
        "schema_result": "PASS",
    }


def exclusion_audit(
    project_root,
    processed,
    manifest_sources,
):
    """TTD/openFDA must not enter the frozen 11-source dataset."""

    audit = []

    for label, folder in [
        ("TTD", "ttd"),
        ("openFDA", "openfda"),
    ]:
        if any(
            label.lower() in s.lower()
            for s in manifest_sources
        ):
            raise RuntimeError(
                f"Excluded source {label} in merge manifest"
            )

        places = [
            project_root / folder,
            project_root / "data" / "raw" / folder,
            project_root / "data" / "processed" / folder,
        ]

        matching = [
            p for p in processed.glob("*.csv")
            if folder in p.name.lower()
        ]

        present = [
            str(p)
            for p in places
            if p.exists()
        ]

        audit.append({
            "excluded_source": label,
            "included_in_manifest": "NO",
            "path_hints_present": (
                " | ".join(present[:8])
                or "none in inspected locations"
            ),
            "processed_root_matching_csv": len(matching),
            "status": "EXCLUDED_FROM_LOCKED_11",
        })

    return pl.DataFrame(audit)


def name_flags(field, prefix):
    value = string_key(field)

    return [
        (value == "").alias(
            prefix + "_blank"
        ),

        value.is_in(
            sorted(PLACEHOLDERS)
        ).alias(
            prefix + "_placeholder_token"
        ),

        (
            (value != "")
            & (value.str.len_chars() <= 2)
        ).alias(
            prefix + "_short_1_2"
        ),

        value.str.contains(
            ID_SHAPED
        ).alias(
            prefix + "_identifier_shaped"
        ),
    ]


def main():
    ap = argparse.ArgumentParser(
        description=__doc__
    )

    ap.add_argument(
        "--processed-dir",
        type=Path,
        default=ROOT / "data" / "processed",
    )

    ap.add_argument(
        "--save-tables",
        action="store_true",
    )

    ap.add_argument(
        "--verify-hashes",
        action="store_true",
        help="Verify all 11 source and 3 output SHA256 hashes",
    )

    ap.add_argument(
        "--examples",
        type=int,
        default=15,
    )

    args = ap.parse_args()

    if not 1 <= args.examples <= 100:
        ap.error("--examples must be 1..100")

    base = args.processed_dir.expanduser().resolve()

    root = (
        base.parent.parent
        if base.name == "processed"
        and base.parent.name == "data"
        else ROOT
    )

    manifest_file = (
        base / "ALL_INCLUDED_SOURCES_merge_summary.json"
    )

    if not manifest_file.is_file():
        raise RuntimeError(
            f"Missing merge manifest: {manifest_file}"
        )

    manifest = json.loads(
        manifest_file.read_text(encoding="utf8")
    )

    must(
        "locked merger status",
        manifest.get("status"),
        "success",
    )

    specs = manifest.get("sources", [])

    manifest_sources = [
        x.get("source")
        for x in specs
    ]

    must(
        "11 source manifest entries",
        len(specs),
        11,
    )

    must(
        "exact included sources",
        set(manifest_sources),
        SOURCES,
    )

    if len(set(manifest_sources)) != len(manifest_sources):
        raise RuntimeError(
            "Duplicate source entries in manifest"
        )

    warnings = []

    print("\nANALYSIS 07 | SOURCE / SCHEMA / NAME QC")
    print("Source membership is determined by the merge manifest.")
    print("Name flags are QC findings, not clinical conclusions.")

    # ========================================================
    # STEP 1 - SOURCE FILE AUDIT
    # ========================================================

    print(
        "\n[STEP 1] Auditing 11 original processed inputs",
        flush=True,
    )

    source_rows = []

    for spec in tqdm(
        specs,
        desc="Included source file checks",
        unit="source",
    ):
        source = spec["source"]

        path = source_file_path(
            root,
            base,
            spec,
        )

        cols = csv_header(path)

        missing = [
            c for c in CORE
            if c not in cols
        ]

        if missing or cols[:10] != CORE:
            raise RuntimeError(
                f"{source}: invalid harmonised schema. "
                f"Missing={missing}; first columns={cols[:10]}"
            )

        observed = int(
            scan(path)
            .select(pl.len())
            .collect(engine="streaming")
            .item()
        )

        expected = int(spec["observed"])

        must(
            f"source rows: {source}",
            observed,
            expected,
        )

        if args.verify_hashes:
            expected_hash = spec.get("sha256")

            if not expected_hash:
                raise RuntimeError(
                    f"No source SHA256 recorded: {source}"
                )

            must(
                f"source SHA256: {source}",
                hash_file(path),
                expected_hash,
            )

        source_rows.append({
            "source": source,
            "source_csv": str(path),
            "rows": observed,
            "columns": len(cols),
            "required_fields": "10/10",
            "schema": "PASS",
            "sha256": (
                "VERIFIED"
                if args.verify_hashes
                else "NOT_CHECKED"
            ),
        })

    source_audit = (
        pl.DataFrame(source_rows)
        .sort("rows", descending=True)
    )

    must(
        "total source evidence events",
        int(source_audit["rows"].sum()),
        EXPECTED_ROWS["events"],
    )

    # ========================================================
    # STEP 2 - FINAL INTEGRATED SCHEMAS
    # ========================================================

    print(
        "\n[STEP 2] Validating all 3 frozen CSV schemas",
        flush=True,
    )

    outputs = {
        k: base / v
        for k, v in FINAL.items()
    }

    schema = manifest.get("schema", {})

    schema_keys = {
        "events": "full",
        "assertions": "assertions",
        "pairs": "pairs",
    }

    schemas = []

    for label in [
        "events",
        "assertions",
        "pairs",
    ]:
        required_prefix = (
            CORE
            if label != "pairs"
            else CORE[:6]
        )

        stated_schema = (
            schema.get(schema_keys[label])
            if isinstance(schema, dict)
            else None
        )

        schemas.append(
            file_schema_audit(
                label,
                outputs[label],
                required_prefix,
                EXPECTED_COLUMNS[label],
                stated_schema,
            )
        )

    schemas_df = pl.DataFrame(schemas)

    event_cols = set(csv_header(outputs["events"]))
    assert_cols = set(csv_header(outputs["assertions"]))
    pair_cols = set(csv_header(outputs["pairs"]))

    more_required = {
        "events": {
            "merge_event_id",
            "merge_assertion_id",
            "merge_pair_id",
            "merge_source_row_number",
            "merge_evidence_family",
            "merge_original_source",
            "merge_drug_name_status",
            "merge_disease_name_status",
        },

        "assertions": {
            "merge_assertion_id",
            "merge_pair_id",
            "supporting_event_count",
            "representative_event_id",
        },

        "pairs": {
            "merge_pair_id",
            "supporting_event_count",
            "supporting_assertion_count",
            "supporting_source_count",
            "supporting_sources",
            "relationship_types",
            "evidence_types",
            "potential_conflict",
        },
    }

    observed_columns = {
        "events": event_cols,
        "assertions": assert_cols,
        "pairs": pair_cols,
    }

    for label, expected in more_required.items():
        observed = observed_columns[label]

        if not expected.issubset(observed):
            raise RuntimeError(
                f"{label}: missing fields "
                f"{sorted(expected - observed)}"
            )

    for k, manifest_key in [
        ("events", "source_events"),
        ("assertions", "assertions"),
        ("pairs", "exact_identity_pairs"),
    ]:
        must(
            f"manifest {k}",
            int(manifest[manifest_key]),
            EXPECTED_ROWS[k],
        )

    if args.verify_hashes:
        map_hash = {
            "events": "merged",
            "assertions": "dedup",
            "pairs": "pairs",
        }

        hashes = manifest.get(
            "output_sha256",
            {},
        )

        for k, hash_key in map_hash.items():
            expected_hash = hashes.get(hash_key)

            if not expected_hash:
                raise RuntimeError(
                    f"Missing output checksum: {hash_key}"
                )

            must(
                f"output SHA256 {k}",
                hash_file(outputs[k]),
                expected_hash,
            )

    # ========================================================
    # STEP 3 - ASSERTION-LEVEL NAME QUALITY
    # ========================================================

    print(
        "\n[STEP 3] Validating drug and disease names",
        flush=True,
    )

    assertion_cols = [
        "source",
        "drug_name",
        "disease_or_condition_name",
        "supporting_event_count",
    ]

    assertion_lf = (
        scan(
            outputs["assertions"],
            assertion_cols,
        )
        .with_columns(
            pl.col("supporting_event_count")
            .cast(pl.Int64, strict=True)
            .alias("_events")
        )
        .with_columns(
            name_flags("drug_name", "drug")
            + name_flags(
                "disease_or_condition_name",
                "disease",
            )
        )
    )

    name_source = (
        assertion_lf
        .group_by("source")
        .agg([
            pl.len().alias("assertions"),

            pl.col("_events")
            .sum()
            .alias("events"),

            pl.col("drug_blank")
            .sum()
            .alias("drug_blank"),

            pl.col("disease_blank")
            .sum()
            .alias("disease_blank"),

            pl.col("drug_placeholder_token")
            .sum()
            .alias("drug_placeholder_tokens"),

            pl.col("disease_placeholder_token")
            .sum()
            .alias("disease_placeholder_tokens"),

            pl.col("drug_short_1_2")
            .sum()
            .alias("drug_short_1_2"),

            pl.col("disease_short_1_2")
            .sum()
            .alias("disease_short_1_2"),

            pl.col("drug_identifier_shaped")
            .sum()
            .alias("drug_identifier_shaped"),

            pl.col("disease_identifier_shaped")
            .sum()
            .alias("disease_identifier_shaped"),

            pl.when(
                pl.col("drug_blank")
                | pl.col("disease_blank")
            )
            .then(pl.col("_events"))
            .otherwise(0)
            .sum()
            .alias("events_with_blank_name"),
        ])
        .collect(engine="streaming")
        .sort("events", descending=True)
    )

    must(
        "assertion source set",
        set(name_source["source"].to_list()),
        SOURCES,
    )

    must(
        "source assertions",
        int(name_source["assertions"].sum()),
        EXPECTED_ROWS["assertions"],
    )

    must(
        "source-linked events",
        int(name_source["events"].sum()),
        EXPECTED_ROWS["events"],
    )

    actual_source_counts = {
        x["source"]: x["events"]
        for x in name_source.to_dicts()
    }

    expected_source_counts = {
        x["source"]: int(x["observed"])
        for x in specs
    }

    must(
        "assertion-to-input source conservation",
        actual_source_counts,
        expected_source_counts,
    )

    must(
        "empty drug-name fields",
        int(name_source["drug_blank"].sum()),
        0,
    )

    must(
        "empty disease-name fields",
        int(name_source["disease_blank"].sum()),
        0,
    )

    must(
        "events linked to empty names",
        int(name_source["events_with_blank_name"].sum()),
        0,
    )

    for attr in [
        "drug_placeholder_tokens",
        "disease_placeholder_tokens",
        "drug_short_1_2",
        "disease_short_1_2",
        "drug_identifier_shaped",
        "disease_identifier_shaped",
    ]:
        n = int(name_source[attr].sum())

        if n:
            warn(
                f"{attr}: {n:,} assertions require review",
                warnings,
            )

    name_source = name_source.with_columns([
        (
            100
            * pl.col("drug_placeholder_tokens")
            / pl.col("assertions")
        )
        .round(4)
        .alias("pct_drug_placeholder"),

        (
            100
            * pl.col("disease_placeholder_tokens")
            / pl.col("assertions")
        )
        .round(4)
        .alias("pct_disease_placeholder"),
    ])

    # ========================================================
    # STEP 4 - ORIGINAL EVENT NAME PROVENANCE
    # ========================================================

    print(
        "\n[STEP 4] Validating event provenance "
        "and final pair support",
        flush=True,
    )

    event_status = (
        scan(
            outputs["events"],
            [
                "source",
                "merge_drug_name_status",
                "merge_disease_name_status",
            ],
        )
        .group_by([
            "source",
            "merge_drug_name_status",
            "merge_disease_name_status",
        ])
        .agg(
            pl.len().alias("events")
        )
        .collect(engine="streaming")
        .sort("events", descending=True)
    )

    must(
        "event ledger name-status conservation",
        int(event_status["events"].sum()),
        EXPECTED_ROWS["events"],
    )

    observed_event_counts = {
        row["source"]: sum(
            int(x["events"])
            for x in event_status.to_dicts()
            if x["source"] == row["source"]
        )
        for row in source_rows
    }

    must(
        "event source counts vs manifest",
        observed_event_counts,
        expected_source_counts,
    )

    # The source_name status was assigned during integration.
    # It cannot prove that every source processor obtained
    # a genuine human-readable label.

    pair_summary = (
        scan(
            outputs["pairs"],
            [
                "supporting_event_count",
                "supporting_assertion_count",
                "potential_conflict",
            ],
        )
        .with_columns([
            pl.col("supporting_event_count")
            .cast(pl.Int64, strict=True),

            pl.col("supporting_assertion_count")
            .cast(pl.Int64, strict=True),
        ])
        .select([
            pl.len().alias("pairs"),

            pl.col("supporting_event_count")
            .sum()
            .alias("events"),

            pl.col("supporting_assertion_count")
            .sum()
            .alias("assertions"),

            (
                pl.col("potential_conflict")
                .str.to_lowercase()
                == "true"
            )
            .sum()
            .alias("potential_conflicts"),
        ])
        .collect(engine="streaming")
        .row(0, named=True)
    )

    for attr, expected in [
        ("pairs", EXPECTED_ROWS["pairs"]),
        ("events", EXPECTED_ROWS["events"]),
        ("assertions", EXPECTED_ROWS["assertions"]),
        ("potential_conflicts", EXPECTED_CONFLICTS),
    ]:
        must(
            f"pair summary {attr}",
            int(pair_summary[attr]),
            expected,
        )

    # ========================================================
    # STEP 5 - EXCLUDED AND AUXILIARY FILE AUDIT
    # ========================================================

    print(
        "\n[STEP 5] Auditing exclusions and auxiliary CSVs",
        flush=True,
    )

    excluded = exclusion_audit(
        root,
        base,
        SOURCES,
    )

    included_paths = {
        str(source_file_path(root, base, s))
        for s in specs
    }

    final_paths = {
        str(p.resolve())
        for p in outputs.values()
    }

    auxiliary_rows = []

    for file in sorted(base.glob("*.csv")):
        absolute = str(file.resolve())

        if (
            absolute in included_paths
            or absolute in final_paths
        ):
            continue

        try:
            columns = len(csv_header(file))
            status = "HEADER_READ_OK"

        except Exception as exc:
            columns = 0
            status = (
                f"HEADER_ERROR_{type(exc).__name__}"
            )

            warn(
                f"Auxiliary CSV needs review: "
                f"{file.name}: {exc}",
                warnings,
            )

        auxiliary_rows.append({
            "file": file.name,
            "columns": columns,
            "classification":
                "NOT_PART_OF_LOCKED_SOURCE_MANIFEST",
            "header_read": status,
        })

    if auxiliary_rows:
        auxiliaries = pl.DataFrame(
            auxiliary_rows
        )
    else:
        auxiliaries = pl.DataFrame({
            "file": [],
            "columns": [],
            "classification": [],
            "header_read": [],
        })

    # ========================================================
    # SUSPICIOUS NAME EXAMPLES
    # ========================================================

    suspect_count = sum(
        int(name_source[name].sum())
        for name in [
            "drug_placeholder_tokens",
            "disease_placeholder_tokens",
            "drug_identifier_shaped",
            "disease_identifier_shaped",
            "drug_short_1_2",
            "disease_short_1_2",
        ]
    )

    examples = pl.DataFrame({
        "source": [],
        "drug_name": [],
        "disease_or_condition_name": [],
        "supporting_event_count": [],
    })

    if suspect_count:
        examples = (
            assertion_lf
            .filter(
                pl.col("drug_placeholder_token")
                | pl.col("disease_placeholder_token")
                | pl.col("drug_identifier_shaped")
                | pl.col("disease_identifier_shaped")
                | pl.col("drug_short_1_2")
                | pl.col("disease_short_1_2")
            )
            .select([
                "source",
                "drug_name",
                "disease_or_condition_name",
                "supporting_event_count",
            ])
            .limit(args.examples)
            .collect(engine="streaming")
        )

    # ========================================================
    # OVERALL SUMMARY
    # ========================================================

    overall = pl.DataFrame([
        {
            "metric": "included_source_files",
            "value": 11,
        },
        {
            "metric": "full_evidence_events",
            "value": EXPECTED_ROWS["events"],
        },
        {
            "metric": "source_specific_assertions",
            "value": EXPECTED_ROWS["assertions"],
        },
        {
            "metric": "exact_identity_pairs",
            "value": EXPECTED_ROWS["pairs"],
        },
        {
            "metric": "potential_conflict_pairs",
            "value": EXPECTED_CONFLICTS,
        },
        {
            "metric": "full_merged_columns",
            "value": EXPECTED_COLUMNS["events"],
        },
        {
            "metric": "assertion_columns",
            "value": EXPECTED_COLUMNS["assertions"],
        },
        {
            "metric": "pair_columns",
            "value": EXPECTED_COLUMNS["pairs"],
        },
        {
            "metric": "assertions_with_blank_drug_name",
            "value": int(
                name_source["drug_blank"].sum()
            ),
        },
        {
            "metric": "assertions_with_blank_disease_name",
            "value": int(
                name_source["disease_blank"].sum()
            ),
        },
        {
            "metric": "drug_name_placeholder_token_assertions",
            "value": int(
                name_source["drug_placeholder_tokens"].sum()
            ),
        },
        {
            "metric": "disease_name_placeholder_token_assertions",
            "value": int(
                name_source["disease_placeholder_tokens"].sum()
            ),
        },
        {
            "metric": "drug_id_shaped_name_assertions",
            "value": int(
                name_source["drug_identifier_shaped"].sum()
            ),
        },
        {
            "metric": "disease_id_shaped_name_assertions",
            "value": int(
                name_source["disease_identifier_shaped"].sum()
            ),
        },
        {
            "metric": "excluded_sources_in_manifest",
            "value": 0,
        },
        {
            "metric": "non_final_root_csv_files",
            "value": auxiliaries.height,
        },
        {
            "metric": "review_warnings",
            "value": len(warnings),
        },
    ])

    # ========================================================
    # PRINT ALL TABLES
    # ========================================================

    print_table(
        "1. FINAL LOCKED DATASET SUMMARY",
        overall,
    )

    print_table(
        "2. SOURCE FILES IN FROZEN MANIFEST",
        source_audit,
    )

    print_table(
        "3. FINAL OUTPUT SCHEMA VALIDATION",
        schemas_df,
    )

    print_table(
        "4. SOURCE-LEVEL DRUG/DISEASE NAME QC",
        name_source,
    )

    print_table(
        "5. EVENT NAME PROVENANCE STATUS DISTRIBUTION",
        event_status,
    )

    print_table(
        "6. EXCLUDED DATABASE AUDIT",
        excluded,
    )

    print_table(
        "7. OTHER CSV FILES (NOT SOURCE FAILURES)",
        auxiliaries,
    )

    print_table(
        "8. FLAGGED NAME EXAMPLES FOR REVIEW",
        examples,
    )

    # ========================================================
    # INTERPRETATION
    # ========================================================

    print("\n[INTERPRETATION]")

    print(
        "All 11 source snapshots and all 3 integrated "
        "output structures passed the tested checks."
    )

    print(
        "Blank names are hard QC failures. Short, "
        "placeholder-like and identifier-shaped names "
        "are review flags."
    )

    print(
        "Source-name provenance status is merger metadata, "
        "not independent verification of original names."
    )

    print(
        "Extra CSV files are not additional included sources "
        "and are not automatically schema failures."
    )

    print(
        "TTD and openFDA remain excluded by the manifest, "
        "regardless of whether raw files exist on disk."
    )

    print(
        "Passing schema validation does not establish "
        "clinical validity or correct ontology mapping."
    )

    print(
        f"[PASS] ANALYSIS 07 STRUCTURAL AUDIT; "
        f"{len(warnings)} review warning(s)"
    )

    # ========================================================
    # OPTIONAL SAVE
    # ========================================================

    if args.save_tables:
        output = root / "tables"

        output.mkdir(
            parents=True,
            exist_ok=True,
        )

        tables = {
            "overall": overall,
            "included_source_audit": source_audit,
            "final_schema": schemas_df,
            "name_quality_by_source": name_source,
            "name_provenance": event_status,
            "excluded_sources": excluded,
            "auxiliary_csv_files": auxiliaries,
            "name_flag_examples": examples,
        }

        for label, frame in tables.items():
            dst = (
                output
                / f"Analysis_07_{label}.csv"
            )

            frame.write_csv(dst)

            print(f"[SAVED] {dst}")

        summary = {
            "snapshot": "paper11-locked-v1",
            "polars_version": pl.__version__,
            "results": overall.to_dicts(),
            "warnings": warnings,
            "sha256_verified": bool(
                args.verify_hashes
            ),
        }

        (
            output / "Analysis_07_summary.json"
        ).write_text(
            json.dumps(
                summary,
                indent=2,
            ) + "\n",
            encoding="utf8",
        )

        print(
            "[SAVED]",
            output / "Analysis_07_summary.json",
        )


if __name__ == "__main__":
    main()
