
#!/usr/bin/env python3
"""Process local repoDB full export into merge-compatible evidence and audited views.

Run: python repoDB.py

No source download. Uses existing:
    data/raw/repodb/repodb_full_download.tsv

Writes to the existing data/processed/repodb/ files after source validation.

repoDB status is benchmark/development evidence, not an independently verified
clinical efficacy judgment. repoDB also derives records from AACT and DrugCentral.
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

import pandas as pd


# ============================================================
# CONFIGURATION
# ============================================================

ROOT = Path(
    os.environ.get(
        "DRUG_DISEASE_ROOT",
        Path(__file__).resolve().parent,
    )
).resolve()

RAW = ROOT / "data/raw/repodb/repodb_full_download.tsv"
OUT_DIR = ROOT / "data/processed/repodb"

SOURCE_URL = "https://unmtid-shinyapps.net/shiny/repodb/"
DOI = "10.1038/sdata.2017.29"

NA = {"", "na", "n/a", "none", "null", "nan"}

REQUIRED = {
    "drug_name",
    "drugbank_id",
    "ind_name",
    "ind_id",
    "status",
    "phase",
    "NCT",
}

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
    "status",
    "development_status",
    "outcome_evidence_class",
    "phase",
    "nct_id",
    "nct_id_normalised",
    "detailed_status",
    "detailed_status_normalised",
    "evidence_text",
    "repodb_source_row_number",
    "repodb_source_record_json",
    "repodb_source_line_number",
    "repodb_name_resolution_drug",
    "repodb_name_resolution_disease",
    "repodb_status_normalised",
    "repodb_nct_ids_extracted_json",
    "repodb_evidence_granularity",
]

for entity in ("drug", "disease"):
    EXTRA.extend(
        f"{entity}_identifier_{suffix}"
        for suffix in (
            "original",
            "namespace",
            "canonical",
            "validation",
        )
    )

COLUMNS = CORE + EXTRA

REJECT_COLUMNS = [
    "repodb_source_row_number",
    "source_line_number",
    "reason",
    "source_record_json",
]

NAME_AUDIT_COLUMNS = [
    "repodb_source_row_number",
    "drug_name_state",
    "disease_name_state",
    "drug_identifier_state",
    "disease_identifier_state",
    "action",
]

MALFORMED_COLUMNS = [
    "data_record_number",
    "source_line_number",
    "fields_found",
    "expected_fields",
    "record_preview",
]


# ============================================================
# GENERAL UTILITIES
# ============================================================

def stamp():
    return datetime.now(timezone.utc).isoformat()


def display(title):
    print(
        "\n" + "=" * 88
        + "\n" + title
        + "\n" + "=" * 88,
        flush=True,
    )


def sha256(path):
    h = hashlib.sha256()

    with Path(path).open("rb") as f:
        for piece in iter(
            lambda: f.read(4 * 1024 * 1024),
            b"",
        ):
            h.update(piece)

    return h.hexdigest()


def raw_value(x):
    return "" if x is None else str(x).strip()


def valid_value(x):
    s = raw_value(x)

    return "" if s.casefold() in NA else s


def value_state(x):
    s = raw_value(x)

    if not s:
        return "empty"

    if s.casefold() in NA:
        return "na_literal"

    return "present"


# ============================================================
# REPO DB STATUS CLASSIFICATION
# ============================================================

def normalise_status(status):
    s = valid_value(status).casefold()

    s = re.sub(
        r"[\s/_-]+",
        "_",
        s,
    ).strip("_")

    aliases = {
        "approve": "approved",
        "approval": "approved",
        "approved": "approved",

        "program_terminated": "program_terminated",
        "trial_halted": "trial_halted",

        "not_approved": "not_approved",
        "notapprove": "not_approved",

        "terminated": "terminated",
        "withdrawn": "withdrawn",
        "suspended": "suspended",
        "failed": "failed",
    }

    return aliases.get(
        s,
        s or "unspecified",
    )


def classify_status(status):
    """
    Preserve the source classification.

    Development discontinuation and repoDB's
    nonapproved benchmark class must not
    automatically imply lack of efficacy.
    """
    s = normalise_status(status)

    if s == "approved":
        return (
            s,
            "indication",
            "repodb_approved_indication",
            "repoDB_approved",
        )

    if s in {
        "program_terminated",
        "trial_halted",
        "terminated",
        "withdrawn",
        "suspended",
    }:
        return (
            s,
            "discontinued_development",
            f"repodb_{s}",
            "development_discontinued_or_halted",
        )

    if s in {
        "not_approved",
        "failed",
    }:
        return (
            s,
            "nonapproved_development_record",
            f"repodb_{s}",
            "repoDB_nonapproved_benchmark",
        )

    if s == "unspecified":
        return (
            s,
            "drug_indication_with_status",
            "repodb_status_unspecified",
            "status_not_reported",
        )

    return (
        s,
        "drug_indication_with_status",
        "repodb_status_other",
        "other_reported_status",
    )


# ============================================================
# DETAILED STATUS INTERPRETATION
# ============================================================

REASONS = (
    (
        "possible_efficacy_reason",
        r"\black of (?:clinical )?efficacy\b"
        r"|\binsufficient efficacy\b"
        r"|\bfutility\b"
        r"|\bfailed to meet (?:the )?"
        r"(?:primary|secondary|efficacy) endpoint\b"
        r"|\bno clinical benefit\b",
    ),
    (
        "possible_safety_reason",
        r"\bsafety\b"
        r"|\badverse event\b"
        r"|\btoxicit"
        r"|\btolerab",
    ),
    (
        "possible_operational_reason",
        r"\brecruit"
        r"|\benrol"
        r"|\baccrual\b"
        r"|\bfund"
        r"|\bfinancial\b"
        r"|\bbusiness\b"
        r"|\bstrateg"
        r"|\bcommercial\b"
        r"|\bsponsor\b"
        r"|\blogistic",
    ),
)


def outcome_reason(status, detail):
    """
    Text-based reason flags are candidates for review,
    not verified negative clinical outcomes.
    """
    if normalise_status(status) == "approved":
        return "not_applicable_approved"

    text = valid_value(detail).casefold()

    if not text:
        return "reason_not_reported"

    matches = [
        label
        for label, pattern in REASONS
        if re.search(pattern, text)
    ]

    return (
        "|".join(matches)
        if matches
        else "other_or_unclassified_reason"
    )


# ============================================================
# IDENTIFIER VALIDATION
# ============================================================

def identifier(x, entity):
    """
    Return:
       original, namespace, canonical, validation

    Validation is syntactic only. No external
    cross-database equivalence is inferred.
    """
    orig = raw_value(x)
    token = valid_value(orig).upper()

    if not token:
        return (
            orig,
            "",
            "",
            "missing",
        )

    if entity == "drug":

        if re.fullmatch(r"DB\d{5}", token):
            return (
                orig,
                "DrugBank",
                token,
                "syntax_valid",
            )

        return (
            orig,
            "unresolved",
            "",
            "unrecognised_or_invalid",
        )

    if re.fullmatch(
        r"(?:UMLS:)?C\d{7}",
        token,
    ):
        return (
            orig,
            "UMLS",
            "UMLS:" + token.removeprefix("UMLS:"),
            "syntax_valid",
        )

    for pat, label, prefix in (
        (
            r"MESH:[CD]\d{6,9}",
            "MeSH",
            "MESH",
        ),
        (
            r"OMIM:\d{6}",
            "OMIM",
            "OMIM",
        ),
        (
            r"DOID:\d+",
            "DOID",
            "DOID",
        ),
        (
            r"MONDO:\d{1,7}",
            "MONDO",
            "MONDO",
        ),
    ):

        if re.fullmatch(pat, token):

            if prefix == "MONDO":
                token = (
                    "MONDO:"
                    + token.split(":")[1].zfill(7)
                )

            return (
                orig,
                label,
                token,
                "syntax_valid",
            )

    return (
        orig,
        "unresolved",
        "",
        "unrecognised_or_invalid",
    )


# ============================================================
# INPUT FORMAT VALIDATION
# ============================================================

def detect_delimiter(path):
    with path.open(
        "r",
        encoding="utf-8-sig",
        newline="",
    ) as f:
        head = f.readline()

    return (
        "\t"
        if head.count("\t") >= head.count(",")
        else ","
    )


def parse_original(path, delimiter):
    """
    Strict parser with source record and physical
    line tracking.

    Malformed records are not silently skipped.
    """
    rows = []
    bad = []

    with path.open(
        "r",
        encoding="utf-8-sig",
        newline="",
    ) as f:

        reader = csv.reader(
            f,
            delimiter=delimiter,
            quotechar='"',
            strict=True,
        )

        header = [
            c.strip()
            for c in next(reader)
        ]

        if (
            not REQUIRED.issubset(header)
            or len(header) != len(set(header))
        ):
            raise RuntimeError(
                f"Invalid repoDB header: {header}; "
                f"missing={sorted(REQUIRED - set(header))}"
            )

        for number, fields in enumerate(
            reader,
            1,
        ):

            if len(fields) != len(header):

                bad.append({
                    "data_record_number": number,
                    "source_line_number": reader.line_num,
                    "fields_found": len(fields),
                    "expected_fields": len(header),
                    "record_preview": str(fields)[:300],
                })

                continue

            rows.append(
                (
                    number,
                    reader.line_num,
                    dict(zip(header, fields)),
                )
            )

    return header, rows, bad


# ============================================================
# TRANSFORM ONE REPO DB RECORD
# ============================================================

def make_record(source_number, physical_line, source):

    dn = raw_value(source["drug_name"])
    di = raw_value(source["drugbank_id"])

    cn = raw_value(source["ind_name"])
    ci = raw_value(source["ind_id"])

    ds = value_state(dn)
    cs = value_state(cn)

    did = identifier(di, "drug")
    cid = identifier(ci, "disease")

    dname = valid_value(dn)
    cname = valid_value(cn)

    # Preserve missing-name records where an
    # original identifier is available.
    #
    # A placeholder is NOT a verified drug name.

    if not dname and valid_value(di):

        dname = valid_value(di)
        dname_resolution = "identifier_placeholder"

    else:

        dname_resolution = (
            "provided_name"
            if dname
            else "unresolved"
        )

    if not cname and valid_value(ci):

        cname = valid_value(ci)
        cname_resolution = "identifier_placeholder"

    else:

        cname_resolution = (
            "provided_name"
            if cname
            else "unresolved"
        )

    audit = {
        "repodb_source_row_number": source_number,
        "drug_name_state": ds,
        "disease_name_state": cs,
        "drug_identifier_state": value_state(di),
        "disease_identifier_state": value_state(ci),
        "action": (
            "retained_with_placeholder"
            if (
                not valid_value(dn)
                or not valid_value(cn)
            )
            else "retained"
        ),
    }

    # Reject only when there is no usable
    # name AND no usable identifier.

    if not dname or not cname:

        audit["action"] = (
            "rejected_no_name_or_identifier"
        )

        return (
            None,
            audit,
            {
                "repodb_source_row_number": source_number,
                "source_line_number": physical_line,
                "reason": (
                    "missing_drug_identity"
                    if not dname
                    else "missing_disease_identity"
                ),
                "source_record_json": json.dumps(
                    source,
                    ensure_ascii=False,
                ),
            },
        )

    status = raw_value(source["status"])
    phase = raw_value(source["phase"])
    nct = raw_value(source["NCT"])

    detailed = raw_value(
        source.get("DetailedStatus", "")
    )

    (
        norm_status,
        relationship,
        evidence,
        dev_status,
    ) = classify_status(status)

    # Extract NCT identifiers without turning one
    # source row into several independent trial
    # observations.

    nct_ids = list(dict.fromkeys(
        x.upper()
        for x in re.findall(
            r"\bNCT\d{8}\b",
            nct,
            flags=re.I,
        )
    ))

    normalized_nct = (
        nct_ids[0]
        if len(nct_ids) == 1
        else ""
    )

    parts = [
        f"status={status}",
        f"development_category={dev_status}",
    ]

    for k, v in (
        ("phase", phase),
        ("NCT", nct),
        ("DetailedStatus", detailed),
    ):
        if valid_value(v):
            parts.append(f"{k}={v}")

    output = dict(
        drug_name=dname,

        drug_identifier=valid_value(di),

        drug_identifier_type=did[1],

        disease_or_condition_name=cname,

        disease_or_condition_identifier=valid_value(ci),

        disease_or_condition_identifier_type=cid[1],

        relationship_type=relationship,

        evidence_type=evidence,

        source="repoDB",

        internal_source="repoDB full database",

        status=status,

        development_status=dev_status,

        outcome_evidence_class=outcome_reason(
            status,
            detailed,
        ),

        phase=phase,

        nct_id=nct,

        nct_id_normalised=normalized_nct,

        detailed_status=detailed,

        detailed_status_normalised=valid_value(detailed),

        evidence_text="; ".join(parts),

        repodb_source_row_number=str(source_number),

        repodb_source_record_json=json.dumps(
            source,
            ensure_ascii=False,
        ),

        repodb_source_line_number=str(physical_line),

        repodb_name_resolution_drug=dname_resolution,

        repodb_name_resolution_disease=cname_resolution,

        repodb_status_normalised=norm_status,

        repodb_nct_ids_extracted_json=json.dumps(
            nct_ids
        ),

        repodb_evidence_granularity=(
            "repoDB_source_status_record"
        ),
    )

    for entity, vals in (
        ("drug", did),
        ("disease", cid),
    ):

        for suffix, val in zip(
            (
                "original",
                "namespace",
                "canonical",
                "validation",
            ),
            vals,
        ):

            output[
                f"{entity}_identifier_{suffix}"
            ] = val

    return output, audit, None


# ============================================================
# DEDUPLICATION WITH SOURCE PROVENANCE
# ============================================================

def dedup_with_provenance(full):
    """
    Preserve the full source-event CSV unchanged.

    The additional dedup view groups matching
    assertions while retaining original row links.
    """
    keys = [
        "drug_name",
        "drug_identifier",
        "drug_identifier_type",
        "disease_or_condition_name",
        "disease_or_condition_identifier",
        "disease_or_condition_identifier_type",
        "relationship_type",
        "evidence_type",
        "status",
        "phase",
        "nct_id",
        "detailed_status",
        "outcome_evidence_class",
    ]

    groups = full.groupby(
        keys,
        sort=False,
        dropna=False,
    )

    base = (
        groups.first()
        .reset_index()[COLUMNS]
    )

    counts = (
        groups.size()
        .rename("repodb_supporting_record_count")
        .reset_index()
    )

    links = (
        groups["repodb_source_row_number"]
        .agg("|".join)
        .rename("repodb_supporting_source_rows")
        .reset_index()
    )

    dedup = (
        base.merge(
            counts,
            on=keys,
            validate="one_to_one",
        )
        .merge(
            links,
            on=keys,
            validate="one_to_one",
        )
    )

    if int(
        dedup[
            "repodb_supporting_record_count"
        ].sum()
    ) != len(full):
        raise RuntimeError(
            "Dedup provenance row accounting failed"
        )

    return dedup, keys


def legacy_duplicate_audit(full):
    """
    Show where a broad association-level collapse
    would have hidden different trials or reasons.
    """
    legacy = [
        "drug_name",
        "drug_identifier",
        "disease_or_condition_name",
        "disease_or_condition_identifier",
        "relationship_type",
        "evidence_type",
        "status",
        "phase",
    ]

    items = []

    for key, g in full.groupby(
        legacy,
        sort=False,
        dropna=False,
    ):

        if len(g) < 2:
            continue

        obj = dict(zip(legacy, key))

        obj.update(
            record_count=len(g),

            source_rows="|".join(
                g["repodb_source_row_number"]
            ),

            distinct_nct_count=(
                g["nct_id"].nunique()
            ),

            distinct_detailed_status_count=(
                g["detailed_status"].nunique()
            ),

            warning=(
                "distinct trial/reason information"
                if (
                    g["nct_id"].nunique() > 1
                    or g["detailed_status"].nunique() > 1
                )
                else "repeated same trial/reason assertion"
            ),
        )

        items.append(obj)

    return pd.DataFrame(
        items,
        columns=legacy + [
            "record_count",
            "source_rows",
            "distinct_nct_count",
            "distinct_detailed_status_count",
            "warning",
        ],
    )


# ============================================================
# SAFE OUTPUT STAGING
# ============================================================

def stage_csv(data, path):
    path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    tmp = path.with_name(
        path.name + ".processing.tmp"
    )

    data.to_csv(
        tmp,
        index=False,
    )

    check = (
        pd.read_csv(
            tmp,
            nrows=0,
        )
        .columns.tolist()
    )

    if check != list(data.columns):
        raise RuntimeError(
            "Staged CSV schema failed: " + str(path)
        )

    return tmp


# ============================================================
# MAIN PROCESSING
# ============================================================

def run(args):

    inp = args.input.resolve()
    dest = args.output_dir.resolve()

    if (
        not inp.is_file()
        or inp.stat().st_size == 0
    ):
        raise FileNotFoundError(
            f"Existing repoDB raw input missing/empty: {inp}"
        )

    display(
        "1. RAW FILE INVENTORY / SOURCE DETAILS"
    )

    for path in sorted(
        x for x in inp.parent.rglob("*")
        if x.is_file()
    ):

        print(
            "  ",
            "USED" if path == inp else "OTHER",
            path.name,
            f"{path.stat().st_size:,} bytes",
        )

    print(
        "Source:",
        inp,
        "SHA256:",
        sha256(inp),
    )

    delimiter = detect_delimiter(inp)

    header, rows, malformed = (
        parse_original(
            inp,
            delimiter,
        )
    )

    print(
        "Delimiter:",
        repr(delimiter),
        "Header:",
        " | ".join(header),
    )

    print(
        "Parsed data records:",
        len(rows),
        "Malformed rows:",
        len(malformed),
    )

    if rows:
        print(
            "First raw record:",
            json.dumps(
                rows[0][2],
                ensure_ascii=False,
            )[:1600],
        )

    if malformed:

        bad_path = (
            dest / "repodb_malformed_lines.csv"
        )

        dest.mkdir(
            parents=True,
            exist_ok=True,
        )

        pd.DataFrame(
            malformed,
            columns=MALFORMED_COLUMNS,
        ).to_csv(
            bad_path,
            index=False,
        )

        raise RuntimeError(
            f"Malformed repoDB rows: {len(malformed)}; "
            f"audit: {bad_path}. "
            "Existing main outputs unchanged"
        )

    display(
        "2. CLASSIFICATION AND SOURCE-TO-MERGE RULES"
    )

    print(
        "ONE FILE, NO TABLE JOIN. "
        "Preserve every repoDB source status event."
    )

    print(
        "DrugBank ID and UMLS/other disease ID "
        "remain original; optional syntax-only "
        "canonical fields."
    )

    print(
        "Approved means repoDB's approved classification, "
        "not new independent approval verification."
    )

    print(
        "Not Approved, Trial Halted, Program Terminated "
        "are NOT automatically lack of efficacy."
    )

    print(
        "Missing names with present IDs use flagged "
        "ID placeholders; missing both identity "
        "fields are rejected."
    )

    print(
        "Underlying repoDB integrates AACT and DrugCentral: "
        "overlap must be handled in the final merger."
    )

    full_rows = []
    rejects = []
    name_audit = []

    statuses = Counter()

    for number, line, raw in rows:

        converted, audit, rejected = (
            make_record(
                number,
                line,
                raw,
            )
        )

        statuses[
            normalise_status(raw["status"])
        ] += 1

        if audit["action"] != "retained":
            name_audit.append(audit)

        if rejected is not None:
            rejects.append(rejected)

        else:
            full_rows.append(converted)

    full = pd.DataFrame(
        full_rows,
        columns=COLUMNS,
    )

    if (
        full.empty
        or len(full) + len(rejects) != len(rows)
    ):
        raise RuntimeError(
            "No valid repoDB output or "
            "source row conservation failed"
        )

    dedup, dedup_keys = (
        dedup_with_provenance(full)
    )

    audit = legacy_duplicate_audit(full)

    rejected_df = pd.DataFrame(
        rejects,
        columns=REJECT_COLUMNS,
    )

    name_df = pd.DataFrame(
        name_audit,
        columns=NAME_AUDIT_COLUMNS,
    )

    malformed_df = pd.DataFrame(
        [],
        columns=MALFORMED_COLUMNS,
    )

    display(
        "3. COUNTS, STATUSES, SCHEMA, AND OUTPUT PREVIEW"
    )

    print(
        "Source records:",
        len(rows),
        "Retained:",
        len(full),
        "Rejected:",
        len(rejects),
    )

    print(
        "Deduplicated assertions:",
        len(dedup),
        "Repeated evidence:",
        len(full) - len(dedup),
    )

    print(
        "Source status counts:",
        dict(statuses),
    )

    print(
        "Relationship counts:",
        full[
            "relationship_type"
        ].value_counts().to_dict(),
    )

    print(
        "Evidence counts:",
        full[
            "evidence_type"
        ].value_counts().to_dict(),
    )

    print(
        "Drug name resolution:",
        full[
            "repodb_name_resolution_drug"
        ].value_counts().to_dict(),
    )

    print(
        "Disease name resolution:",
        full[
            "repodb_name_resolution_disease"
        ].value_counts().to_dict(),
    )

    print(
        "Drug identifier validation:",
        full[
            "drug_identifier_validation"
        ].value_counts().to_dict(),
    )

    print(
        "Disease identifier validation:",
        full[
            "disease_identifier_validation"
        ].value_counts().to_dict(),
    )

    print(
        "Fields:",
        " | ".join(COLUMNS),
    )

    print(
        "First output row:",
        json.dumps(
            full_rows[0],
            ensure_ascii=False,
            indent=2,
        )[:2500],
    )

    # --------------------------------------------------------
    # WRITE TO THE SAME ORIGINAL PROCESSED DIRECTORY
    # --------------------------------------------------------

    output_files = {
        "repodb_drug_disease.csv": full,

        "repodb_drug_disease_deduplicated.csv":
            dedup,

        "repodb_duplicate_audit.csv": audit,

        "repodb_name_audit.csv": name_df,

        "repodb_malformed_lines.csv":
            malformed_df,

        "repodb_rejected_rows.csv":
            rejected_df,
    }

    report = {
        "status": "success",
        "timestamp_utc": stamp(),
        "source": "repoDB",
        "source_url": SOURCE_URL,
        "citation_doi": DOI,

        "source_sha256": sha256(inp),
        "source_file": str(inp),
        "source_header": header,

        "source_data_rows": len(rows),

        "valid_source_rows": len(full),

        "rejected_identity_rows":
            len(rejects),

        "deduplicated_assertions":
            len(dedup),

        "repeated_evidence_rows":
            len(full) - len(dedup),

        "status_counts":
            dict(statuses),

        "relationship_type_counts":
            full["relationship_type"]
            .value_counts()
            .to_dict(),

        "evidence_type_counts":
            full["evidence_type"]
            .value_counts()
            .to_dict(),

        "drug_name_lookup_counts":
            full["repodb_name_resolution_drug"]
            .value_counts()
            .to_dict(),

        "disease_name_lookup_counts":
            full["repodb_name_resolution_disease"]
            .value_counts()
            .to_dict(),

        "duplicate_audit_groups":
            len(audit),

        "missing_name_audit_rows":
            len(name_df),

        "unique_drug_disease_id_pairs":
            len(
                full[
                    [
                        "drug_identifier",
                        "disease_or_condition_identifier",
                    ]
                ].drop_duplicates()
            ),

        "output_columns":
            COLUMNS,

        "deduplication_keys":
            dedup_keys,

        "scientific_interpretation": (
            "Source reported status from repoDB; "
            "terminated/not-approved/halted does not "
            "by itself establish lack of efficacy. "
            "repoDB reuses AACT and DrugCentral evidence."
        ),

        "record_identity": (
            "repodb_source_row_number; "
            "raw record JSON preserved"
        ),

        "output_files": {
            name: str(dest / name)
            for name in output_files
        },
    }

    staged = {}

    try:

        for name, frame in output_files.items():

            staged[name] = stage_csv(
                frame,
                dest / name,
            )

        report["output_sha256"] = {
            name: sha256(temp)
            for name, temp in staged.items()
        }

        summary_path = (
            dest / "repodb_summary.json"
        )

        summary_tmp = summary_path.with_name(
            summary_path.name + ".processing.tmp"
        )

        summary_tmp.write_text(
            json.dumps(
                report,
                indent=2,
                ensure_ascii=False,
            ) + "\n",
            encoding="utf-8",
        )

        staged[
            "repodb_summary.json"
        ] = summary_tmp

        for name, tmp in staged.items():

            tmp.replace(
                dest / name
            )

            print(
                "[OVERWRITTEN]",
                dest / name,
            )

    finally:

        for tmp in staged.values():
            tmp.unlink(
                missing_ok=True
            )

    print(
        "SUCCESS. Merge the full repodb_drug_disease.csv; "
        "preserve status and NCT evidence."
    )


# ============================================================
# ENTRY POINT
# ============================================================

def main():

    p = argparse.ArgumentParser(
        description=__doc__,
    )

    p.add_argument(
        "--input",
        type=Path,
        default=RAW,
    )

    p.add_argument(
        "--output-dir",
        type=Path,
        default=OUT_DIR,
    )

    args = p.parse_args()

    run(args)


if __name__ == "__main__":
    main()
