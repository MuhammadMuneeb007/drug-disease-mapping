
#!/usr/bin/env python3
"""TTD drug-disease status records -> merge-compatible CSV.

Default: reuse existing P1-05 and P1-03 in data/raw/ttd.
Overwrite data/processed/ttd outputs after validation.
TTD-reported drug development status is not independent
proof of approval or therapeutic efficacy.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
import shutil
import urllib.request
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(
    os.environ.get(
        "DRUG_DISEASE_ROOT",
        Path(__file__).resolve().parent,
    )
).resolve()

RAW = ROOT / "data/raw/ttd"
OUT = ROOT / "data/processed/ttd"

BASE = (
    "https://db.idrblab.net/ttd/"
    "sites/default/files/ttd_database"
)

FILES = (
    "P1-05-Drug_disease.txt",
    "P1-03-TTD_crossmatching.txt",
)

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
    "ttd_drug_id",
    "ttd_drug_id_original",
    "ttd_drug_id_validation",
    "ttd_drug_name_lookup_status",
    "ttd_indication_original",
    "ttd_indication_cells_json",
    "ttd_source_line_number",
    "ttd_source_record_layout",
    "ttd_status_original",
    "clinical_status",
    "ttd_status_class",
    "ttd_status_parse_method",
    "ttd_indication_parse_status",
    "ttd_icd_matches_json",
    "ttd_icd_mapping_status",
    "pubchem_cid",
    "pubchem_cids_json",
    "pubchem_mapping_status",
    "cas",
    "cas_all_json",
    "chebi",
    "chebi_all_json",
    "atc",
    "atc_all_json",
    "crossmatch_record_found",
    "crossmatch_drug_name",
    "evidence_text",
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
    "ttd_drug_id",
    "ttd_source_line_number",
    "reason",
    "indication_cells_json",
    "indication_original",
]

PARSER_AUDIT_COLUMNS = [
    "file",
    "physical_line",
    "reason",
    "line_preview",
]

STATUS_TERMS = [
    "Discontinued in Preclinical",
    "Discontinued in Phase 1/2",
    "Discontinued in Phase 2/3",
    "Discontinued in Phase 3/4",
    "Discontinued in Phase 1",
    "Discontinued in Phase 2",
    "Discontinued in Phase 3",
    "Discontinued in Phase 4",
    "Withdrawn from market",
    "Approved (orphan drug)",
    "Approved in China",
    "Approved in EU",
    "Approved in USA",
    "Application submitted",
    "IND submitted",
    "Registration",
    "Preregistration",
    "Preclinical",
    "Investigative",
    "Patented",
    "Clinical trial",
    "Discontinued",
    "Terminated",
    "Withdrawn",
    "Suspended",
    "Approved",
    "Phase 1/2",
    "Phase 2/3",
    "Phase 3/4",
    "Phase 0",
    "Phase 1",
    "Phase 2",
    "Phase 3",
    "Phase 4",
    "Phase I/II",
    "Phase II/III",
    "Phase I",
    "Phase II",
    "Phase III",
    "Phase IV",
]

STATUS_TERMS.sort(
    key=len,
    reverse=True,
)

STATUS_PATTERN = "|".join(
    re.escape(x)
    for x in STATUS_TERMS
)

STATUS_EXACT = re.compile(
    rf"^(?:{STATUS_PATTERN})$",
    re.IGNORECASE,
)

STATUS_PREFIX = re.compile(
    rf"^\s*(?P<status>{STATUS_PATTERN})"
    rf"(?=$|\s*[|,:;\-–])",
    re.IGNORECASE,
)

STATUS_SUFFIX = re.compile(
    rf"(?<!\w)(?P<status>{STATUS_PATTERN})\s*$",
    re.IGNORECASE,
)

ICD_PATTERN = re.compile(
    r"\[\s*ICD\s*[-_]?\s*"
    r"(9|10|11)\s*:?\s*([^\]]+)\]",
    re.IGNORECASE,
)

TAG_PATTERN = re.compile(
    r"^[A-Z][A-Z0-9_]{2,24}$"
)

IDENTITY_TAGS = {
    "TTDDRUID",
    "TTD_DRUG_ID",
    "DRUGID",
    "DRUG_ID",
}

DATA_TAGS = IDENTITY_TAGS | {
    "DRUGNAME",
    "INDICATI",
    "PUBCHCID",
    "CASNUMBE",
    "CHEBI_ID",
    "SUPDRATC",
}

NA_VALUES = {
    "",
    "#N/A",
    "N/A",
    "NA",
    "NONE",
    "NULL",
    "NAN",
    "-",
}


def timestamp():
    return datetime.now(timezone.utc).isoformat()


def heading(title):
    print(
        "\n" + "=" * 94
        + "\n" + title
        + "\n" + "=" * 94,
        flush=True,
    )


def clean(value):
    if value is None:
        return ""

    return str(value).strip()


def valid(value):
    result = clean(value)

    if result.upper() in NA_VALUES:
        return ""

    return result


def sha256(path):
    h = hashlib.sha256()

    with Path(path).open("rb") as stream:
        for block in iter(
            lambda: stream.read(4 * 1024 * 1024),
            b"",
        ):
            h.update(block)

    return h.hexdigest()


def obtain_file(path, allow_download):
    if path.is_file() and path.stat().st_size:
        print("[REUSE]", path)
        return

    if not allow_download:
        raise FileNotFoundError(
            f"Missing {path}. "
            "Pass --download-missing explicitly."
        )

    path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    tmp = path.with_name(
        path.name + ".download"
    )

    try:
        url = BASE + "/" + path.name

        print("[DOWNLOAD]", url)

        with (
            urllib.request.urlopen(
                url,
                timeout=180,
            ) as source,
            tmp.open("wb") as target,
        ):
            shutil.copyfileobj(
                source,
                target,
            )

        if tmp.stat().st_size == 0:
            raise RuntimeError(
                "Downloaded an empty source file"
            )

        tmp.replace(path)

    finally:
        tmp.unlink(missing_ok=True)


def looks_like_tag(value):
    return bool(
        TAG_PATTERN.fullmatch(value)
    )


def parse_file(path):
    """
    Supports both source structures:

        DRUG_ID  TAG  VALUE
        TAG      VALUE

    Multiple value columns are retained.

    In tag-first format, TTDDRUID introduces
    the current drug record.

    Unrecognised records are audited.
    """
    records = defaultdict(
        lambda: defaultdict(list)
    )

    entries = []
    problems = []

    counts = Counter()
    examples = []

    current_id = ""
    layout = ""

    with path.open(
        "r",
        encoding="utf-8-sig",
        errors="replace",
    ) as source:

        for lineno, line in enumerate(
            source,
            start=1,
        ):
            raw = line.rstrip("\r\n")

            if not raw.strip():
                counts["blank_lines"] += 1
                continue

            pieces = [
                clean(x)
                for x in raw.split("\t")
            ]

            if len(pieces) < 2:
                counts["non_tab_lines"] += 1

                if len(examples) < 10:
                    examples.append(raw[:160])

                continue

            if pieces[0].upper() in DATA_TAGS:
                this_layout = "tag_first"

                tag = pieces[0].upper()
                cells = pieces[1:]

                if tag in IDENTITY_TAGS:
                    current_id = (
                        valid(cells[0])
                        if cells
                        else ""
                    )

                record_id = current_id

            elif (
                len(pieces) >= 3
                and looks_like_tag(
                    pieces[1].upper()
                )
                and pieces[0]
            ):
                this_layout = "id_first"

                record_id = pieces[0]
                tag = pieces[1].upper()
                cells = pieces[2:]

            elif (
                looks_like_tag(
                    pieces[0].upper()
                )
                and pieces[0].upper()
                not in DATA_TAGS
            ):
                this_layout = "tag_first"

                record_id = current_id
                tag = pieces[0].upper()
                cells = pieces[1:]

            else:
                counts["unparsed_tab_lines"] += 1

                problems.append({
                    "file": path.name,
                    "physical_line": lineno,
                    "reason": "unrecognised_layout",
                    "line_preview": raw[:250],
                })

                continue

            if (
                layout
                and this_layout != layout
            ):
                raise RuntimeError(
                    f"{path.name}:{lineno}: "
                    "inconsistent source layouts "
                    f"{layout}/{this_layout}"
                )

            layout = this_layout

            if not record_id:
                if tag in DATA_TAGS:
                    problems.append({
                        "file": path.name,
                        "physical_line": lineno,
                        "reason":
                            "tag_without_ttd_drug_id",
                        "line_preview": raw[:250],
                    })

                    counts["unlinked_tag_lines"] += 1

                continue

            counts["tag_lines"] += 1
            counts["tag_" + tag] += 1

            entry = {
                "drug_id": record_id,
                "tag": tag,
                "cells": cells,
                "physical_line": lineno,
                "layout": this_layout,
                "original": raw,
            }

            entries.append(entry)
            records[record_id][tag].append(
                entry
            )

    if not layout or not entries:
        raise RuntimeError(
            "No structured TTD records could "
            f"be read from {path}"
        )

    print(
        f"[{path.name}] "
        f"layout={layout} "
        f"unique_drugs={len(records):,} "
        f"tag_lines={counts['tag_lines']:,} "
        f"indications={counts['tag_INDICATI']:,}"
    )

    tag_counts = Counter({
        k[4:]: value
        for k, value in counts.items()
        if k.startswith("tag_")
    })

    print(
        "Tag distribution:",
        dict(tag_counts.most_common(30)),
    )

    if examples:
        print(
            "Non-tab preamble examples:",
            examples[:3],
        )

    return (
        records,
        counts,
        problems,
        layout,
    )


def record_values(fields, tag):
    return [
        cell
        for entry in fields.get(tag, [])
        for cell in entry["cells"]
        if valid(cell)
    ]


def resolve_name(fields):
    names = {
        valid(x)
        for x in record_values(
            fields,
            "DRUGNAME",
        )
        if valid(x)
    }

    if len(names) == 1:
        return (
            next(iter(names)),
            "unambiguous_source_name",
        )

    if names:
        return (
            "",
            "conflicting_source_names",
        )

    return (
        "",
        "missing_source_name",
    )


def crossmatch_fields(record):
    pub = record_values(
        record,
        "PUBCHCID",
    )

    cas = record_values(
        record,
        "CASNUMBE",
    )

    chebi = record_values(
        record,
        "CHEBI_ID",
    )

    atc = record_values(
        record,
        "SUPDRATC",
    )

    pubchem_valid = set()

    for item in pub:
        candidate = re.sub(
            r"^(?:PUBCHEM(?:\s*CID|\s*_CID)?|CID)"
            r"\s*[:# ]\s*",
            "",
            item,
            flags=re.IGNORECASE,
        ).strip()

        if re.fullmatch(
            r"\d+",
            candidate,
        ):
            pubchem_valid.add(
                str(int(candidate))
            )

    pubchem_id = (
        next(iter(pubchem_valid))
        if len(pubchem_valid) == 1
        else ""
    )

    if len(pubchem_valid) == 1:
        pubchem_status = (
            "unique_numeric_pubchem_cid"
        )

    elif len(pubchem_valid) > 1:
        pubchem_status = (
            "ambiguous_multiple_cids"
        )

    elif pub:
        pubchem_status = "present_unparsed"

    else:
        pubchem_status = "missing"

    cas_clean = sorted({
        re.sub(
            r"^CAS\s*:?[ ]*",
            "",
            x,
            flags=re.IGNORECASE,
        ).strip()
        for x in cas
    })

    chebi_clean = sorted({
        x.upper()
        for x in chebi
    })

    atc_clean = sorted(set(atc))

    return {
        "pubchem_cid": pubchem_id,
        "pubchem_cids_json":
            json.dumps(pub, ensure_ascii=False),
        "pubchem_mapping_status":
            pubchem_status,

        "cas": (
            cas_clean[0]
            if len(cas_clean) == 1
            else ""
        ),

        "cas_all_json":
            json.dumps(cas, ensure_ascii=False),

        "chebi": (
            chebi_clean[0]
            if len(chebi_clean) == 1
            else ""
        ),

        "chebi_all_json":
            json.dumps(chebi, ensure_ascii=False),

        "atc": (
            atc_clean[0]
            if len(atc_clean) == 1
            else ""
        ),

        "atc_all_json":
            json.dumps(atc, ensure_ascii=False),
    }


def classify_status(raw):
    """
    Clinical statuses are TTD-reported,
    not independently verified approvals
    or efficacy outcomes.
    """
    status = clean(raw).casefold()

    if not status:
        return (
            "status_not_reported",
            "ttd_status_unreported",
            "ttd_status_unreported",
        )

    if status.startswith("approved"):
        return (
            "ttd_reported_approved",
            "ttd_reported_approved_indication",
            "ttd_curated_approved_status",
        )

    if any(
        term in status
        for term in (
            "withdrawn",
            "discontinued",
            "terminated",
            "suspended",
        )
    ):
        return (
            "development_halted_or_withdrawn",
            "ttd_discontinued_development",
            "ttd_" + re.sub(
                r"\W+",
                "_",
                status,
            ).strip("_"),
        )

    if (
        status.startswith("phase")
        or status.startswith("clinical trial")
        or status in {
            "registration",
            "preregistration",
            "application submitted",
        }
    ):
        return (
            "clinical_development",
            "ttd_clinical_development",
            "ttd_" + re.sub(
                r"\W+",
                "_",
                status,
            ).strip("_"),
        )

    if (
        status.startswith("preclinical")
        or status.startswith("ind submitted")
    ):
        return (
            "preclinical_or_ind",
            "ttd_preclinical_development",
            "ttd_" + re.sub(
                r"\W+",
                "_",
                status,
            ).strip("_"),
        )

    if status in {
        "patented",
        "investigative",
    }:
        return (
            "early_research",
            "ttd_early_research",
            "ttd_" + re.sub(
                r"\W+",
                "_",
                status,
            ).strip("_"),
        )

    return (
        "unclassified_status",
        "ttd_unclassified_status",
        "ttd_unclassified_status",
    )


def parse_indication(entry):
    """
    Parse disease, status and ICD annotation.

    Only unambiguous ICD annotations enter
    the standard identifier field.

    All original source cells remain available
    in the source evidence columns.
    """
    values = [
        valid(x)
        for x in entry["cells"]
    ]

    values = [
        x for x in values if x
    ]

    if not values:
        return (
            None,
            "empty_indication",
        )

    disease_text = values[0]

    status = ""
    method = "none"

    if (
        len(values) > 1
        and STATUS_EXACT.fullmatch(disease_text)
    ):
        status = disease_text
        disease_text = values[1]
        method = "leading_explicit_column"

    if len(values) > 1:
        status_cells = [
            x
            for x in values[1:]
            if STATUS_EXACT.fullmatch(x)
        ]

        if (
            len(set(status_cells)) == 1
            and status_cells
            and not status
        ):
            status = status_cells[0]
            method = "separate_explicit_column"

        elif (
            status_cells
            and any(
                s.casefold() != status.casefold()
                for s in status_cells
            )
        ):
            return (
                None,
                "conflicting_explicit_statuses",
            )

        elif len(set(status_cells)) > 1:
            return (
                None,
                "conflicting_explicit_statuses",
            )

    if not status:
        suffix = STATUS_SUFFIX.search(
            disease_text
        )

        prefix = STATUS_PREFIX.match(
            disease_text
        )

        if suffix and suffix.start() > 0:
            status = suffix.group("status")

            disease_text = (
                disease_text[:suffix.start()]
                .strip(" \t,;|:-–")
            )

            method = (
                "anchored_suffix_in_indication"
            )

        elif (
            prefix
            and prefix.end() < len(disease_text)
        ):
            status = prefix.group("status")

            disease_text = (
                disease_text[prefix.end():]
                .strip(" \t,;|:-–")
            )

            method = (
                "anchored_prefix_in_indication"
            )

    icd_matches = []

    for cell in values:
        for match in ICD_PATTERN.finditer(
            cell
        ):
            icd_matches.append({
                "version":
                    "ICD-" + match.group(1),
                "code_raw":
                    match.group(2).strip(),
            })

    disease_name = " ".join(
        ICD_PATTERN.sub(
            "",
            disease_text,
        ).split()
    ).strip(" ,;|:-–")

    if (
        not disease_name
        or disease_name.upper() in NA_VALUES
    ):
        return (
            None,
            "missing_or_invalid_disease_name",
        )

    unique_icd = {
        (
            x["version"],
            x["code_raw"],
        )
        for x in icd_matches
        if x["code_raw"]
    }

    if len(unique_icd) == 1:
        icd_version, icd_code = next(
            iter(unique_icd)
        )

        mapping_status = (
            "single_raw_source_icd_annotation_unverified"
        )

    else:
        icd_version = ""
        icd_code = ""

        mapping_status = (
            "multiple_icd_annotations_ambiguous"
            if unique_icd
            else "no_icd_annotation"
        )

    extra_values = [
        x
        for x in values[1:]
        if (
            x != disease_text
            and not STATUS_EXACT.fullmatch(x)
            and not ICD_PATTERN.fullmatch(x)
        )
    ]

    parse_status = (
        "review_extra_source_cells"
        if extra_values
        else "parsed"
    )

    if not status:
        parse_status = (
            "review_extra_source_cells_and_status"
            if extra_values
            else "unrecognised_or_missing_status"
        )

    return ({
        "disease_name": disease_name,
        "icd_code": icd_code,
        "icd_version": icd_version,
        "icd_mapping_status":
            mapping_status,
        "icd_matches": icd_matches,
        "status": status,
        "status_parse_method": method,
        "indication_parse_status":
            parse_status,
    }, None)


def process_records(
    drug_records,
    crossmatch_records,
):
    full = []
    rejected = []

    counts = Counter()

    for drug_id, fields in drug_records.items():
        source_name, source_name_status = (
            resolve_name(fields)
        )

        crossmatch = crossmatch_records.get(
            drug_id
        )

        if crossmatch is not None:
            cross_name, _ = resolve_name(
                crossmatch
            )
        else:
            cross_name = ""

        if source_name:
            drug_name = source_name
            name_status = source_name_status

        elif cross_name:
            drug_name = cross_name
            name_status = (
                "unambiguous_crossmatch_name"
            )

        else:
            drug_name = drug_id
            name_status = (
                "identifier_placeholder"
            )

        external = crossmatch_fields(
            crossmatch or {}
        )

        indications = fields.get(
            "INDICATI",
            [],
        )

        if not indications:
            counts[
                "drugs_without_indication"
            ] += 1

        for entry in indications:
            counts[
                "source_indication_rows"
            ] += 1

            parsed, reason = (
                parse_indication(entry)
            )

            if reason:
                rejected.append({
                    "ttd_drug_id": drug_id,
                    "ttd_source_line_number":
                        entry["physical_line"],
                    "reason": reason,
                    "indication_cells_json":
                        json.dumps(
                            entry["cells"],
                            ensure_ascii=False,
                        ),
                    "indication_original":
                        entry["original"],
                })

                counts[
                    "rejected_indications"
                ] += 1

                continue

            (
                status_class,
                relationship,
                evidence,
            ) = classify_status(
                parsed["status"]
            )

            if status_class in {
                "unclassified_status",
                "status_not_reported",
            }:
                counts[
                    "unclassified_or_missing_status"
                ] += 1

            if (
                parsed["indication_parse_status"]
                != "parsed"
            ):
                counts[
                    "indication_review_flags"
                ] += 1

            icd_namespace = parsed[
                "icd_version"
            ]

            icd_code = parsed[
                "icd_code"
            ]

            row = {
                "drug_name": drug_name,
                "drug_identifier": drug_id,
                "drug_identifier_type": "TTD",

                "disease_or_condition_name":
                    parsed["disease_name"],

                "disease_or_condition_identifier":
                    icd_code,

                "disease_or_condition_identifier_type":
                    icd_namespace,

                "relationship_type":
                    relationship,

                "evidence_type":
                    evidence,

                "source": "TTD",

                "internal_source":
                    "TTD_P1-05_Drug_disease",

                "ttd_drug_id": drug_id,
                "ttd_drug_id_original": drug_id,

                "ttd_drug_id_validation":
                    "source_key_present_not_externally_verified",

                "ttd_drug_name_lookup_status":
                    name_status,

                "ttd_indication_original":
                    entry["original"],

                "ttd_indication_cells_json":
                    json.dumps(
                        entry["cells"],
                        ensure_ascii=False,
                    ),

                "ttd_source_line_number":
                    str(entry["physical_line"]),

                "ttd_source_record_layout":
                    entry["layout"],

                "ttd_status_original":
                    parsed["status"],

                "clinical_status":
                    parsed["status"],

                "ttd_status_class":
                    status_class,

                "ttd_status_parse_method":
                    parsed["status_parse_method"],

                "ttd_indication_parse_status":
                    parsed["indication_parse_status"],

                "ttd_icd_matches_json":
                    json.dumps(
                        parsed["icd_matches"],
                        ensure_ascii=False,
                    ),

                "ttd_icd_mapping_status":
                    parsed["icd_mapping_status"],

                **external,

                "crossmatch_record_found":
                    str(
                        crossmatch is not None
                    ),

                "crossmatch_drug_name":
                    cross_name,

                "evidence_text": (
                    "TTD-reported clinical status: "
                    f"{parsed['status'] or 'unreported'}; "
                    "source INDICATI line: "
                    f"{entry['physical_line']}"
                ),

                "drug_identifier_original":
                    drug_id,

                "drug_identifier_namespace":
                    "TTD",

                "drug_identifier_canonical": "",

                "drug_identifier_validation":
                    "source_specific_only",

                "disease_identifier_original":
                    icd_code,

                "disease_identifier_namespace":
                    icd_namespace,

                "disease_identifier_canonical": "",

                "disease_identifier_validation": (
                    "raw_ICD_annotation_not_ontology_validated"
                    if icd_code
                    else "missing_or_ambiguous"
                ),
            }

            full.append(row)

            counts[
                "retained_indications"
            ] += 1

    if (
        counts["source_indication_rows"]
        != len(full) + len(rejected)
    ):
        raise RuntimeError(
            "TTD source-row accounting failed"
        )

    if not full:
        raise RuntimeError(
            "No eligible TTD records; "
            "existing output is unchanged"
        )

    return full, rejected, counts


def deduplicate_with_provenance(full):
    """
    The main CSV preserves all evidence events.

    This separate view collapses identical
    source assertions while retaining all
    supporting original line references.
    """
    keys = [
        "ttd_drug_id",
        "disease_or_condition_name",
        "disease_or_condition_identifier",
        "disease_or_condition_identifier_type",
        "relationship_type",
        "evidence_type",
        "ttd_status_original",
        "ttd_indication_cells_json",
    ]

    groups = {}

    for record in full:
        key = tuple(
            record[field]
            for field in keys
        )

        if key not in groups:
            groups[key] = {
                **record,

                "ttd_supporting_source_lines": [],

                "ttd_supporting_record_count": 0,
            }

        groups[key][
            "ttd_supporting_record_count"
        ] += 1

        groups[key][
            "ttd_supporting_source_lines"
        ].append(
            record["ttd_source_line_number"]
        )

    deduplicated = []

    for row in groups.values():
        row["ttd_supporting_source_lines"] = (
            "|".join(
                row["ttd_supporting_source_lines"]
            )
        )

        deduplicated.append(row)

    support_total = sum(
        row["ttd_supporting_record_count"]
        for row in deduplicated
    )

    if support_total != len(full):
        raise RuntimeError(
            "Deduplicated support accounting failed"
        )

    return deduplicated, keys


def stage_csv(destination, columns, rows):
    destination.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    temp = destination.with_name(
        destination.name + ".processing.tmp"
    )

    with temp.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as stream:
        writer = csv.DictWriter(
            stream,
            fieldnames=columns,
            extrasaction="ignore",
        )

        writer.writeheader()
        writer.writerows(rows)

    with temp.open(
        "r",
        encoding="utf-8",
        newline="",
    ) as stream:
        actual = next(
            csv.reader(stream)
        )

    if actual != columns:
        raise RuntimeError(
            "Staged output schema mismatch"
        )

    return temp


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
    )

    parser.add_argument(
        "--workdir",
        type=Path,
        default=RAW,
    )

    parser.add_argument(
        "--output-dir",
        type=Path,
        default=OUT,
    )

    parser.add_argument(
        "--download-missing",
        action="store_true",
    )

    args = parser.parse_args()

    raw_directory = (
        args.workdir.expanduser().resolve()
    )

    output_directory = (
        args.output_dir.expanduser().resolve()
    )

    paths = [
        raw_directory / filename
        for filename in FILES
    ]

    for path in paths:
        obtain_file(
            path,
            args.download_missing,
        )

    heading(
        "1. ORIGINAL RAW FILE INVENTORY"
    )

    inventory = []

    for path in sorted(
        raw_directory.rglob("*")
    ):
        if not path.is_file():
            continue

        item = {
            "path": str(
                path.relative_to(raw_directory)
            ),
            "size_bytes":
                path.stat().st_size,
            "selected":
                path in paths,
        }

        inventory.append(item)

        print(
            "INPUT" if item["selected"]
            else "OTHER",
            item["path"],
            f"{item['size_bytes']:,} bytes",
        )

    hashes = {
        path.name: sha256(path)
        for path in paths
    }

    heading(
        "2. PARSE TTD MULTILINE RECORDS"
    )

    (
        drug_records,
        drug_stats,
        drug_audit,
        drug_layout,
    ) = parse_file(paths[0])

    (
        crossmatch_records,
        crossmatch_stats,
        crossmatch_audit,
        crossmatch_layout,
    ) = parse_file(paths[1])

    if drug_stats["tag_INDICATI"] == 0:
        raise RuntimeError(
            "P1-05 contains no INDICATI records. "
            "Inspect the actual source format."
        )

    parse_audit = (
        drug_audit + crossmatch_audit
    )

    if parse_audit:
        audit_path = (
            output_directory
            / "ttd_parse_audit.csv"
        )

        temp = stage_csv(
            audit_path,
            PARSER_AUDIT_COLUMNS,
            parse_audit,
        )

        temp.replace(audit_path)

        raise RuntimeError(
            f"{len(parse_audit)} ambiguous "
            "tagged source lines. Review "
            f"{audit_path}. "
            "Main output not overwritten."
        )

    heading(
        "3. TTD DRUG ID JOIN AND "
        "SCIENTIFIC EVIDENCE RULES"
    )

    print(
        "JOIN: P1-05 TTD drug ID "
        "-> P1-03 TTD drug ID"
    )

    print(
        "Crossmatches enrich records; "
        "they do not multiply associations."
    )

    print(
        "Primary source ID: TTD drug ID. "
        "PubChem/CAS/ChEBI/ATC retained separately."
    )

    print(
        "ONE output row per original "
        "INDICATI record."
    )

    print(
        "Approval, clinical development, "
        "preclinical and discontinued statuses "
        "are distinct."
    )

    print(
        "Unknown statuses are retained "
        "as unclassified rather than "
        "called therapeutic indications."
    )

    full, rejected, counts = (
        process_records(
            drug_records,
            crossmatch_records,
        )
    )

    deduplicated, dedup_keys = (
        deduplicate_with_provenance(full)
    )

    heading(
        "4. RECORD COUNTS AND "
        "MERGE-COMPATIBLE OUTPUT"
    )

    print(
        "Original indication records:",
        counts["source_indication_rows"],
    )

    print(
        "Retained records:",
        len(full),
    )

    print(
        "Rejected records:",
        len(rejected),
    )

    print(
        "Deduplicated assertions:",
        len(deduplicated),
    )

    print(
        "Source records supported by dedup:",
        sum(
            row["ttd_supporting_record_count"]
            for row in deduplicated
        ),
    )

    print(
        "Clinical status classes:",
        dict(Counter(
            row["ttd_status_class"]
            for row in full
        )),
    )

    print(
        "Indication parsing statuses:",
        dict(Counter(
            row["ttd_indication_parse_status"]
            for row in full
        )),
    )

    print(
        "ICD mapping statuses:",
        dict(Counter(
            row["ttd_icd_mapping_status"]
            for row in full
        )),
    )

    print(
        "PubChem mapping statuses:",
        dict(Counter(
            row["pubchem_mapping_status"]
            for row in full
        )),
    )

    print(
        "Crossmatch-linked evidence rows:",
        sum(
            row["crossmatch_record_found"]
            == "True"
            for row in full
        ),
    )

    print(
        "FULL OUTPUT SCHEMA:",
        " | ".join(COLUMNS),
    )

    print(
        "FIRST ACTUAL OUTPUT ROW:"
    )

    print(
        json.dumps(
            full[0],
            indent=2,
            ensure_ascii=False,
        )[:3500]
    )

    outputs = {
        "ttd_drug_disease.csv": (
            COLUMNS,
            full,
        ),

        "ttd_drug_disease_deduplicated.csv": (
            COLUMNS + [
                "ttd_supporting_record_count",
                "ttd_supporting_source_lines",
            ],
            deduplicated,
        ),

        "ttd_rejected_indications.csv": (
            REJECT_COLUMNS,
            rejected,
        ),

        "ttd_parse_audit.csv": (
            PARSER_AUDIT_COLUMNS,
            parse_audit,
        ),
    }

    summary = {
        "status": "success",
        "timestamp_utc": timestamp(),

        "source": "TTD",
        "citation_doi": "10.1093/nar/gkaf1154",

        "license_note": (
            "Verify TTD redistribution terms. "
            "Free access does not itself "
            "establish open redistribution rights."
        ),

        "source_inventory": inventory,
        "input_sha256": hashes,

        "input_tagged_lines": {
            "P1-05": dict(drug_stats),
            "P1-03": dict(crossmatch_stats),
        },

        "source_layouts": {
            "P1-05": drug_layout,
            "P1-03": crossmatch_layout,
        },

        "unique_drug_records":
            len(drug_records),

        "unique_crossmatch_records":
            len(crossmatch_records),

        "source_indication_events":
            counts["source_indication_rows"],

        "retained_rows":
            len(full),

        "rejected_indications":
            len(rejected),

        "deduplicated_assertions":
            len(deduplicated),

        "dedup_keys":
            dedup_keys,

        "clinical_status_counts":
            dict(Counter(
                row["clinical_status"]
                or "[missing]"
                for row in full
            )),

        "status_class_counts":
            dict(Counter(
                row["ttd_status_class"]
                for row in full
            )),

        "relation_counts":
            dict(Counter(
                row["relationship_type"]
                for row in full
            )),

        "parse_status_counts":
            dict(Counter(
                row["ttd_indication_parse_status"]
                for row in full
            )),

        "icd_mapping_counts":
            dict(Counter(
                row["ttd_icd_mapping_status"]
                for row in full
            )),

        "pubchem_mapping_counts":
            dict(Counter(
                row["pubchem_mapping_status"]
                for row in full
            )),

        "drug_name_lookup_status":
            dict(Counter(
                row["ttd_drug_name_lookup_status"]
                for row in full
            )),

        "output_columns":
            COLUMNS,

        "source_id_join": (
            "P1-05 TTD drug ID == "
            "P1-03 TTD drug ID; "
            "zero-to-one enrichment"
        ),

        "scientific_interpretation": (
            "TTD-reported development statuses, "
            "not independently verified "
            "regulatory approval or efficacy. "
            "ICD annotations are source values, "
            "not externally ontology-validated."
        ),
    }

    staged = {}

    summary_path = (
        output_directory / "ttd_summary.json"
    )

    try:
        for filename, (
            columns,
            records,
        ) in outputs.items():
            destination = (
                output_directory / filename
            )

            staged[destination] = stage_csv(
                destination,
                columns,
                records,
            )

        summary["output_sha256"] = {
            path.name: sha256(temp)
            for path, temp in staged.items()
        }

        summary["output_files"] = [
            str(path)
            for path in staged
        ]

        summary_path.parent.mkdir(
            parents=True,
            exist_ok=True,
        )

        summary_temp = (
            summary_path.with_name(
                summary_path.name
                + ".processing.tmp"
            )
        )

        summary_temp.write_text(
            json.dumps(
                summary,
                indent=2,
                ensure_ascii=False,
            ) + "\n",
            encoding="utf-8",
        )

        staged[
            summary_path
        ] = summary_temp

        heading(
            "5. OVERWRITING ORIGINAL "
            "PROCESSED TTD FILES"
        )

        for destination, temp in (
            staged.items()
        ):
            temp.replace(destination)

            print(
                "[WRITTEN]",
                destination,
            )

    finally:
        for temp in staged.values():
            temp.unlink(missing_ok=True)

    print(
        "SUCCESS: use the FULL "
        "ttd_drug_disease.csv in "
        "the evidence-preserving merger."
    )


if __name__ == "__main__":
    main()
