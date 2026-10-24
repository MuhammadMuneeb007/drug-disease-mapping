#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Analysis 23:
WHO Essential Medicines / Approved-Drug Coverage Validation

Purpose
-------
Evaluate whether drugs listed in external reference datasets are represented
in the locked, harmonised 11-source drug-disease database.

Reference files
---------------
data/reference/who_essential_medicines.csv
data/reference/approved_drugs.csv

At least one reference file must exist.

Accepted drug-name columns
--------------------------
medicine
drug_name
name
active_ingredient
struct_name
ingredient

Matching hierarchy
------------------
1. EXACT_NAME
2. CONTROLLED_SYNONYM
3. COMBINATION_LABEL_MATCH
4. FORM_VARIANT_REVIEW
5. COMPONENTS_PRESENT_ONLY
6. PARTIAL_COMPONENTS
7. NOT_FOUND

Only categories 1-3 are counted as covered.

IMPORTANT
---------
This is DRUG-NAME COVERAGE, not indication-specific validation.

A reference drug appearing in at least one dataset record does not establish:
- Regulatory approval
- Correct indication
- Clinical efficacy
- Clinical safety

The original locked files are never modified.

Outputs
-------
Print all principal results directly to stdout.

Optional:
--save-tables
--print-all
--top N

Run
---
POLARS_MAX_THREADS=16 \
python Analysis_23_WHO_EML_or_Approved_Drug_Coverage.py --save-tables

Requirements
------------
pip install polars
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re

from collections import Counter, defaultdict
from pathlib import Path

import polars as pl


# ============================================================
# CONFIGURATION
# ============================================================

ROOT = Path(__file__).resolve().parent

ASSERT_FILE = (
    "ALL_INCLUDED_SOURCES_drug_disease_deduplicated.csv"
)

MANIFEST_FILE = (
    "ALL_INCLUDED_SOURCES_merge_summary.json"
)

REF_FILES = {
    "WHO EML": "who_essential_medicines.csv",
    "Approved drugs": "approved_drugs.csv",
}

NAME_COLUMNS = [
    "medicine",
    "drug_name",
    "name",
    "active_ingredient",
    "struct_name",
    "ingredient",
]

SOURCES = {
    "AACT",
    "Broad Drug Repurposing Hub",
    "ChEMBL",
    "CTD",
    "DrugCentral",
    "Hetionet",
    "MEDI",
    "OpenTargets",
    "PrimeKG",
    "repoDB",
    "SIDER",
}

EXPECTED_ASSERTIONS = 4_534_385
EXPECTED_EVENTS = 10_895_395
EXPECTED_DRUG_KEYS = 158_364


# ============================================================
# CONSERVATIVE DRUG SYNONYMS
# ============================================================

# Explicit, limited equivalents.
#
# Avoid broad fuzzy matching, because different medicines
# can have similar names.

ALIAS_TO_CANONICAL = {
    "acetaminophen": "paracetamol",
    "acyclovir": "aciclovir",
    "salbutamol": "albuterol",
    "adrenaline": "epinephrine",
    "noradrenaline": "norepinephrine",
    "glyceryl trinitrate": "nitroglycerin",
    "hydroxycarbamide": "hydroxyurea",
}


# ============================================================
# FORMULATION REVIEW TERMS
# ============================================================

# These are used for REVIEW ONLY.
#
# Removing these terms can change the meaning of a drug name.
# Consequently, such matches do not count as strict coverage.

FORM_SUFFIXES = {
    "hydrochloride",
    "hcl",
    "sodium",
    "potassium",
    "calcium",
    "magnesium",
    "phosphate",
    "sulfate",
    "sulphate",
    "citrate",
    "maleate",
    "mesylate",
    "besylate",
    "besilate",
    "acetate",
    "tartrate",
    "succinate",
    "fumarate",
    "nitrate",
    "bromide",
    "chloride",
    "monohydrate",
    "dihydrate",
    "trihydrate",
}


# ============================================================
# REGEX
# ============================================================

COMBO_DELIMITER = re.compile(
    r"\s*(?:\+|;|/|\s+and\s+)\s*",
    re.I,
)

MULTISPACE = re.compile(r"\s+")


# ============================================================
# NAME NORMALISATION
# ============================================================

def norm(value: str) -> str:
    """
    Same conservative key as Analyses 01, 04 and 08.

    Lowercase, trim, and collapse whitespace.

    Do not remove punctuation or merge ontology concepts.
    """
    return MULTISPACE.sub(
        " ",
        str(value or "").strip(),
    ).lower()


def canon(value: str) -> str:
    """
    Apply only explicitly curated aliases.
    """
    key = norm(value)

    return ALIAS_TO_CANONICAL.get(
        key,
        key,
    )


def components(value: str) -> tuple[str, ...]:
    """
    Identify explicit combination-medicine labels.

    Examples:
        amoxicillin + clavulanic acid
        amoxicillin/clavulanic acid

    Returns a sorted tuple of components.

    This does NOT imply that separately available ingredients
    constitute a covered combination medicine.
    """
    parts = [
        canon(x)
        for x in COMBO_DELIMITER.split(norm(value))
    ]

    if len(parts) <= 1:
        return ()

    return tuple(
        sorted(
            x for x in parts if x
        )
    )


def form_variant_key(value: str) -> str:
    """
    Remove trailing salt/form terms for review-only matching.
    """
    parts = canon(value).split()

    while (
        len(parts) > 1
        and parts[-1] in FORM_SUFFIXES
    ):
        parts.pop()

    return " ".join(parts)


# ============================================================
# VALIDATION HELPERS
# ============================================================

def check(
    label: str,
    actual,
    expected,
) -> None:

    if actual != expected:
        raise RuntimeError(
            f"[FAIL] {label}: "
            f"{actual!r}, expected {expected!r}"
        )

    suffix = (
        f"{actual:,}"
        if isinstance(actual, int)
        else "OK"
    )

    print(
        f"[PASS] {label}: {suffix}",
        flush=True,
    )


def table(
    title: str,
    rows: list[dict],
    limit: int | None = None,
) -> None:
    """
    Print tab-separated data without Polars truncation.
    """
    print("\n" + "=" * 110)
    print(f"{title} [{len(rows):,} rows]")
    print("=" * 110)

    if not rows:
        print("[NO ROWS]")
        return

    columns = list(rows[0])

    print("\t".join(columns))

    selected = (
        rows
        if limit is None
        else rows[:limit]
    )

    for row in selected:

        values = []

        for key in columns:
            value = row.get(key, "")

            if value is None:
                value = ""

            value = (
                str(value)
                .replace("\t", " ")
                .replace("\r", " ")
                .replace("\n", " ")
            )

            values.append(value)

        print("\t".join(values))

    if (
        limit is not None
        and len(rows) > limit
    ):
        print(
            f"[SHOWING FIRST {limit:,}; "
            "use --print-all or --save-tables]"
        )


def write_csv(
    path: Path,
    rows: list[dict],
    columns=None,
) -> None:
    """
    Save supplementary tables when requested.
    """
    if not rows and not columns:
        return

    cols = (
        columns
        if columns
        else list(rows[0])
    )

    with path.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:

        writer = csv.DictWriter(
            handle,
            fieldnames=cols,
            extrasaction="ignore",
        )

        writer.writeheader()
        writer.writerows(rows)

    print(f"[SAVED] {path}")


def file_sha256(path: Path) -> str:
    """
    Calculate SHA256 for reference-file provenance.
    """
    digest = hashlib.sha256()

    with path.open("rb") as handle:

        for chunk in iter(
            lambda: handle.read(1024 * 1024),
            b"",
        ):
            digest.update(chunk)

    return digest.hexdigest()


# ============================================================
# REFERENCE DATA LOADING
# ============================================================

def read_reference(
    path: Path,
    label: str,
) -> tuple[list[dict], dict]:
    """
    Read a reference CSV.

    Count distinct listed medicine names rather than
    blindly treating every source row as a unique drug.
    """

    if not path.exists():

        print(
            f"[NOT PROVIDED] {label}: {path}"
        )

        return [], {}

    if not path.stat().st_size:
        raise RuntimeError(
            f"Empty reference file: {path}"
        )

    entries = {}

    original_rows = 0
    blanks = 0

    with path.open(
        "r",
        encoding="utf-8-sig",
        newline="",
    ) as handle:

        reader = csv.DictReader(handle)

        headers = reader.fieldnames or []

        lookup = {
            h.strip().lower(): h
            for h in headers
            if h
        }

        field = next(
            (
                lookup[k]
                for k in NAME_COLUMNS
                if k in lookup
            ),
            None,
        )

        if not field:
            raise RuntimeError(
                f"{path.name}: no drug-name column.\n"
                f"Accepted: {NAME_COLUMNS}\n"
                f"Observed: {headers}"
            )

        for record in reader:

            original_rows += 1

            raw_name = MULTISPACE.sub(
                " ",
                (record.get(field) or "").strip(),
            )

            key = norm(raw_name)

            if not key:
                blanks += 1
                continue

            if key not in entries:

                entries[key] = {
                    "reference_list": label,
                    "reference_name": raw_name,
                    "reference_key": key,
                }

    result = sorted(
        entries.values(),
        key=lambda x: x["reference_key"],
    )

    metadata = {
        "reference_list": label,
        "file": str(path),
        "file_sha256": file_sha256(path),
        "input_rows": original_rows,
        "blank_rows": blanks,
        "distinct_reference_names": len(result),
        "name_column": field,
    }

    if not result:
        raise RuntimeError(
            f"{path.name}: no usable reference names"
        )

    return result, metadata


# ============================================================
# BUILD DATASET DRUG NAME INDEX
# ============================================================

def make_index(
    summary: pl.DataFrame,
):
    """
    Build lookup indexes from grouped drug names.

    The complete 4.5-million-assertion dataset is NOT
    converted into a Python list or dictionary.

    Only the much smaller aggregated drug-name vocabulary
    is indexed in memory.
    """

    names: dict[str, dict] = {}

    for row in summary.iter_rows(named=True):

        key = row["_name"]

        if not key:
            continue

        state = names.setdefault(
            key,
            {
                "sources": set(),
                "relationships": set(),
                "evidence_types": set(),
            },
        )

        state["sources"].add(
            row["source"]
        )

        state["relationships"].add(
            row["relationship_type"]
        )

        state["evidence_types"].add(
            row["evidence_type"]
        )

    canonical = defaultdict(set)
    variants = defaultdict(set)
    combos = defaultdict(set)

    for key in names:

        canonical[
            canon(key)
        ].add(key)

        if not components(key):

            variants[
                form_variant_key(key)
            ].add(key)

        else:

            combos[
                components(key)
            ].add(key)

    return (
        names,
        canonical,
        variants,
        combos,
    )


# ============================================================
# MATCH ONE REFERENCE MEDICINE
# ============================================================

def resolve(
    entry: dict,
    names,
    canonical,
    variants,
    combos,
) -> dict:
    """
    Match a reference medicine against the dataset.

    Covered:
        EXACT_NAME
        CONTROLLED_SYNONYM
        COMBINATION_LABEL_MATCH

    Not covered without manual review:
        FORM_VARIANT_REVIEW
        COMPONENTS_PRESENT_ONLY
        PARTIAL_COMPONENTS
        NOT_FOUND
    """

    key = entry["reference_key"]

    parts = components(key)

    direct = (
        {key}
        if key in names
        else set()
    )

    # --------------------------------------------------------
    # CASE 1: EXACT NAME
    # --------------------------------------------------------

    if direct:

        matched = direct

        status = "EXACT_NAME"

        covered = True

    # --------------------------------------------------------
    # CASE 2: COMBINATION MEDICINE
    # --------------------------------------------------------

    elif parts:

        matched = combos.get(
            parts,
            set(),
        )

        if matched:

            status = "COMBINATION_LABEL_MATCH"

            covered = True

        else:

            component_hits = [
                set(
                    canonical.get(
                        part,
                        set(),
                    )
                )
                for part in parts
            ]

            matched = set().union(
                *component_hits
            )

            if all(component_hits):

                status = "COMPONENTS_PRESENT_ONLY"

            elif any(component_hits):

                status = "PARTIAL_COMPONENTS"

            else:

                status = "NOT_FOUND"

            covered = False

    # --------------------------------------------------------
    # CASE 3: CONTROLLED SYNONYM
    # --------------------------------------------------------

    else:

        matched = set(
            canonical.get(
                canon(key),
                set(),
            )
        )

        if matched:

            status = "CONTROLLED_SYNONYM"

            covered = True

        # ----------------------------------------------------
        # CASE 4: FORMULATION VARIANT
        # ----------------------------------------------------

        else:

            matched = set(
                variants.get(
                    form_variant_key(key),
                    set(),
                )
            )

            if matched:

                status = "FORM_VARIANT_REVIEW"

            else:

                status = "NOT_FOUND"

            covered = False

    # --------------------------------------------------------
    # COLLECT SOURCE AND EVIDENCE METADATA
    # --------------------------------------------------------

    src = set()
    rel = set()
    evid = set()

    for label in matched:

        record = names[label]

        src.update(
            record["sources"]
        )

        rel.update(
            record["relationships"]
        )

        evid.update(
            record["evidence_types"]
        )

    return {
        **entry,

        "combination": (
            "yes"
            if parts
            else "no"
        ),

        "component_count": (
            len(parts)
            if parts
            else 1
        ),

        "status": status,

        "strict_covered": (
            "yes"
            if covered
            else "no"
        ),

        "matched_dataset_name_keys": (
            " | ".join(
                sorted(matched)[:15]
            )
        ),

        "matching_key_count": len(matched),

        "n_sources": (
            len(src)
            if covered
            else 0
        ),

        "sources": (
            " | ".join(sorted(src))
            if covered
            else ""
        ),

        "relationship_types": (
            " | ".join(sorted(rel))
            if covered
            else ""
        ),

        "evidence_types": (
            " | ".join(sorted(evid))
            if covered
            else ""
        ),

        "review_sources_or_components": (
            " | ".join(sorted(src))
            if not covered
            else ""
        ),
    }


# ============================================================
# SUMMARY STATISTICS
# ============================================================

def summarise(
    results: list[dict],
):
    """
    Generate source-aware supplementary tables.

    Each reference name is counted at most once per
    source, relationship type or evidence type.
    """

    by_list = defaultdict(list)

    for row in results:

        by_list[
            row["reference_list"]
        ].append(row)

    overall = []
    by_status = []
    by_source = []
    by_rel = []
    by_evid = []

    source_names = sorted(SOURCES)

    for label, rows in sorted(
        by_list.items()
    ):

        total = len(rows)

        covered = [
            x for x in rows
            if x["strict_covered"] == "yes"
        ]

        by_status_count = Counter(
            x["status"]
            for x in rows
        )

        # ----------------------------------------------------
        # OVERALL COVERAGE
        # ----------------------------------------------------

        overall.append({
            "reference_list": label,

            "reference_names": total,

            "strict_covered": len(covered),

            "not_strictly_covered": (
                total - len(covered)
            ),

            "strict_coverage_percent": round(
                100 * len(covered) / total,
                3,
            ),

            "exact_name": (
                by_status_count["EXACT_NAME"]
            ),

            "controlled_synonym": (
                by_status_count["CONTROLLED_SYNONYM"]
            ),

            "combination_label": (
                by_status_count["COMBINATION_LABEL_MATCH"]
            ),

            "review_only": sum(
                n
                for status, n in by_status_count.items()
                if status in {
                    "FORM_VARIANT_REVIEW",
                    "COMPONENTS_PRESENT_ONLY",
                    "PARTIAL_COMPONENTS",
                }
            ),

            "not_found": (
                by_status_count["NOT_FOUND"]
            ),
        })

        # ----------------------------------------------------
        # MATCHING STATUS
        # ----------------------------------------------------

        for status, n in sorted(
            by_status_count.items()
        ):

            by_status.append({
                "reference_list": label,
                "status": status,
                "reference_names": n,

                "pct_of_reference_list": round(
                    100 * n / total,
                    3,
                ),
            })

        # ----------------------------------------------------
        # SOURCE-LEVEL COVERAGE
        # ----------------------------------------------------

        for source in source_names:

            n = sum(
                source in x["sources"].split(" | ")
                for x in covered
            )

            by_source.append({
                "reference_list": label,
                "source": source,

                "covered_reference_names": n,

                "pct_of_reference_list": round(
                    100 * n / total,
                    3,
                ),
            })

        # ----------------------------------------------------
        # RELATIONSHIP AND EVIDENCE COVERAGE
        # ----------------------------------------------------

        relations = Counter()
        evidences = Counter()

        for row in covered:

            relations.update(
                row["relationship_types"].split(" | ")
            )

            evidences.update(
                row["evidence_types"].split(" | ")
            )

        for relation, n in sorted(
            relations.items(),
            key=lambda pair: (
                -pair[1],
                pair[0],
            ),
        ):

            by_rel.append({
                "reference_list": label,
                "relationship_type": relation,
                "covered_reference_names_with_type": n,
            })

        for evidence, n in sorted(
            evidences.items(),
            key=lambda pair: (
                -pair[1],
                pair[0],
            ),
        ):

            by_evid.append({
                "reference_list": label,
                "evidence_type": evidence,
                "covered_reference_names_with_type": n,
            })

    return (
        overall,
        by_status,
        by_source,
        by_rel,
        by_evid,
    )


# ============================================================
# MAIN
# ============================================================

def main() -> None:

    parser = argparse.ArgumentParser(
        description=__doc__
    )

    parser.add_argument(
        "--processed-dir",
        type=Path,
        default=ROOT / "data" / "processed",
    )

    parser.add_argument(
        "--reference-dir",
        type=Path,
        default=ROOT / "data" / "reference",
    )

    parser.add_argument(
        "--output-dir",
        type=Path,
        default=ROOT / "tables",
    )

    parser.add_argument(
        "--save-tables",
        action="store_true",
    )

    parser.add_argument(
        "--print-all",
        action="store_true",
        help="Print every reference entry",
    )

    parser.add_argument(
        "--top",
        type=int,
        default=60,
    )

    args = parser.parse_args()

    if not 1 <= args.top <= 100000:
        parser.error(
            "--top must be positive"
        )

    base = (
        args.processed_dir
        .expanduser()
        .resolve()
    )

    reference = (
        args.reference_dir
        .expanduser()
        .resolve()
    )

    assertion_file = base / ASSERT_FILE

    manifest_path = base / MANIFEST_FILE

    # ========================================================
    # STEP 0: VALIDATE LOCKED INPUT FILES
    # ========================================================

    if not manifest_path.is_file():

        raise RuntimeError(
            f"Missing locked manifest: {manifest_path}"
        )

    if (
        not assertion_file.is_file()
        or not assertion_file.stat().st_size
    ):

        raise RuntimeError(
            f"Missing assertion file: {assertion_file}"
        )

    with assertion_file.open(
        "r",
        encoding="utf-8-sig",
        newline="",
    ) as handle:

        header = next(
            csv.reader(
                handle,
                strict=True,
            ),
            [],
        )

    needed = {
        "source",
        "drug_name",
        "relationship_type",
        "evidence_type",
        "supporting_event_count",
    }

    if (
        not needed.issubset(header)
        or len(header) != len(set(header))
    ):

        raise RuntimeError(
            "Invalid assertion header. Missing: "
            f"{sorted(needed - set(header))}"
        )

    manifest = json.loads(
        manifest_path.read_text(
            encoding="utf-8"
        )
    )

    check(
        "locked manifest status",
        manifest.get("status"),
        "success",
    )

    check(
        "manifest assertions",
        int(manifest["assertions"]),
        EXPECTED_ASSERTIONS,
    )

    check(
        "manifest source events",
        int(manifest["source_events"]),
        EXPECTED_EVENTS,
    )

    # ========================================================
    # STEP 1: READ REFERENCE LISTS
    # ========================================================

    references = []
    audits = []

    for label, filename in REF_FILES.items():

        rows, audit = read_reference(
            reference / filename,
            label,
        )

        references.extend(rows)

        if audit:
            audits.append(audit)

    if not references:

        raise SystemExit(
            "[NO VALIDATION PERFORMED] "
            "Provide at least one reference CSV:\n"
            f"{reference / REF_FILES['WHO EML']}\n"
            "or\n"
            f"{reference / REF_FILES['Approved drugs']}"
        )

    print(
        "\nANALYSIS 23 | "
        "WHO EML / APPROVED-DRUG NAME COVERAGE"
    )

    print(
        "[SCOPE] Drug-name retrieval only; "
        "not indication-specific verification"
    )

    print(
        "[RULE] Component matches and salt/form variants "
        "are review-only"
    )

    # ========================================================
    # STEP 2: STREAM AND AGGREGATE ASSERTIONS
    # ========================================================

    print(
        "[STEP 1] Stream assertion fields "
        "and aggregate source/name/type counts",
        flush=True,
    )

    name = (
        pl.col("drug_name")
        .fill_null("")
        .str.strip_chars()
        .str.replace_all(
            r"\s+",
            " ",
        )
        .str.to_lowercase()
    )

    aggregated = (
        pl.scan_csv(
            str(assertion_file),
            infer_schema=False,
            ignore_errors=False,
            truncate_ragged_lines=False,
            encoding="utf8",
        )

        .select(
            sorted(needed)
        )

        .with_columns([
            name.alias("_name"),

            pl.col("supporting_event_count")
            .cast(
                pl.Int64,
                strict=True,
            ),
        ])

        .group_by([
            "source",
            "_name",
            "relationship_type",
            "evidence_type",
        ])

        .agg([
            pl.len().alias("assertions"),

            pl.col("supporting_event_count")
            .sum()
            .alias("events"),
        ])

        .collect(
            engine="streaming"
        )
    )

    # ========================================================
    # VALIDATE LOCKED DATASET COUNTS
    # ========================================================

    check(
        "locked assertions",
        int(
            aggregated["assertions"].sum()
        ),
        EXPECTED_ASSERTIONS,
    )

    check(
        "linked evidence events",
        int(
            aggregated["events"].sum()
        ),
        EXPECTED_EVENTS,
    )

    check(
        "all 11 sources",
        set(
            aggregated["source"].to_list()
        ),
        SOURCES,
    )

    check(
        "drug search-name keys",
        aggregated["_name"].n_unique(),
        EXPECTED_DRUG_KEYS,
    )

    if aggregated.filter(
        pl.col("_name") == ""
    ).height:

        raise RuntimeError(
            "Blank drug name in locked assertions"
        )

    actual_source_counts = {
        row["source"]: row["events"]

        for row in (
            aggregated
            .group_by("source")
            .agg(
                pl.col("events").sum()
            )
            .to_dicts()
        )
    }

    manifest_counts = {
        row["source"]: int(row["observed"])

        for row in manifest.get(
            "sources",
            [],
        )
    }

    check(
        "per-source event conservation",
        actual_source_counts,
        manifest_counts,
    )

    # ========================================================
    # STEP 3: BUILD LOOKUP INDEXES
    # ========================================================

    print(
        "[STEP 2] Build exact, synonym, "
        "combination and review indexes",
        flush=True,
    )

    (
        names,
        canonical,
        variants,
        combos,
    ) = make_index(aggregated)

    # ========================================================
    # STEP 4: EVALUATE REFERENCE MEDICINES
    # ========================================================

    print(
        "[STEP 3] Evaluate all reference names",
        flush=True,
    )

    outcomes = [
        resolve(
            entry,
            names,
            canonical,
            variants,
            combos,
        )

        for entry in references
    ]

    (
        overall,
        status,
        by_source,
        by_rel,
        by_evid,
    ) = summarise(outcomes)

    # ========================================================
    # VALIDATE COVERAGE COUNTS
    # ========================================================

    check(
        "reference-name classification conservation",
        len(outcomes),

        sum(
            row["reference_names"]
            for row in overall
        ),
    )

    for row in overall:

        check(
            f"{row['reference_list']} "
            "covered/not covered",

            (
                row["strict_covered"]
                + row["not_strictly_covered"]
            ),

            row["reference_names"],
        )

    for row in outcomes:

        if (
            row["strict_covered"] == "yes"
            and not row["sources"]
        ):

            raise RuntimeError(
                "Covered name without source: "
                f"{row}"
            )

    # ========================================================
    # PREPARE SUPPLEMENTARY RESULTS
    # ========================================================

    unverified = [
        row for row in outcomes
        if row["strict_covered"] != "yes"
    ]

    combos_all = [
        row for row in outcomes
        if row["combination"] == "yes"
    ]

    aliases = [
        row for row in outcomes
        if row["status"] == "CONTROLLED_SYNONYM"
    ]

    source_status = sorted(
        by_source,

        key=lambda row: (
            row["reference_list"],
            -row["covered_reference_names"],
            row["source"],
        ),
    )

    # ========================================================
    # PRINT COMPLETE SUMMARY TABLES
    # ========================================================

    table(
        "1. REFERENCE INPUT AUDIT AND PROVENANCE",
        audits,
    )

    table(
        "2. STRICT DRUG-NAME COVERAGE BY REFERENCE LIST",
        overall,
    )

    table(
        "3. COVERAGE CLASSIFICATION DISTRIBUTION",
        status,
    )

    table(
        "4. COVERAGE BY INCLUDED SOURCE",
        source_status,
    )

    table(
        "5. RELATIONSHIP TYPES OF MATCHED NAMES "
        "(NON-EXCLUSIVE)",
        by_rel,
    )

    table(
        "6. EVIDENCE TYPES OF MATCHED NAMES "
        "(NON-EXCLUSIVE)",
        by_evid,
    )

    table(
        "7. UNVERIFIED REFERENCE NAMES / REVIEW CASES",
        unverified,

        None
        if args.print_all
        else args.top,
    )

    table(
        "8. REFERENCE COMBINATION MEDICINES",
        combos_all,

        None
        if args.print_all
        else args.top,
    )

    table(
        "9. CONTROLLED SYNONYM MATCHES",
        aliases,

        None
        if args.print_all
        else args.top,
    )

    table(
        "10. ALL REFERENCE MEDICINE MATCH RESULTS",
        outcomes,

        None
        if args.print_all
        else args.top,
    )

    # ========================================================
    # MANUSCRIPT-READY INTERPRETATION
    # ========================================================

    print("\n" + "=" * 110)
    print("MANUSCRIPT-READY COVERAGE SUMMARY")
    print("=" * 110)

    for row in overall:

        print(
            f"{row['reference_list']}: "
            f"{row['strict_covered']:,} / "
            f"{row['reference_names']:,} "
            "distinct reference name entries "
            f"({row['strict_coverage_percent']:.3f}%) "
            "had exact, controlled-synonym or "
            "single-label combination matches "
            "in the integrated database."
        )

    print()

    print(
        "Missing and component-only reference entries "
        "require manual review."
    )

    print()

    print(
        "A drug found anywhere in the integrated database "
        "does not validate its approved indication, "
        "essential-medicine status, clinical efficacy "
        "or safety."
    )

    print()

    print(
        "The reference CSV edition, provenance and "
        "SHA256 checksum should be reported alongside "
        "the manuscript results."
    )

    # ========================================================
    # SAVE SUPPLEMENTARY TABLES
    # ========================================================

    if args.save_tables:

        output = (
            args.output_dir
            .expanduser()
            .resolve()
        )

        output.mkdir(
            parents=True,
            exist_ok=True,
        )

        results = {
            "reference_audit": audits,

            "overall_coverage": overall,

            "match_status": status,

            "by_source": source_status,

            "by_relationship_type": by_rel,

            "by_evidence_type": by_evid,

            "not_strictly_covered": unverified,

            "combination_entries": combos_all,

            "controlled_synonyms": aliases,

            "all_reference_drugs": outcomes,
        }

        for name, rows in results.items():

            destination = (
                output
                / f"Analysis_23_{name}.csv"
            )

            write_csv(
                destination,
                rows,
            )

    print(
        "\n[PASS] Analysis 23 completed; "
        "locked merged files unchanged"
    )


if __name__ == "__main__":
    main()