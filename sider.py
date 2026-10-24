
#!/usr/bin/env python3
"""SIDER 4.1 -> merge-compatible label phenotype evidence with strict indication view.

Default: reuse pinned raw files in data/raw/sider; overwrite existing processed
files in data/processed/sider only after source validations succeed.

NLP_indication, NLP_precondition, and text_mention are distinct. Only the first
is included in the strict SIDER indication view; no row proves FDA approval.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
import re
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

ROOT = Path(
    os.environ.get(
        "DRUG_DISEASE_ROOT",
        Path(__file__).resolve().parent,
    )
).resolve()

RAW = ROOT / "data/raw/sider"
OUT = ROOT / "data/processed/sider"

URL_PREFIX = "https://sideeffects.embl.de/media/download"

FILE_NAMES = (
    "drug_names.tsv",
    "meddra_all_indications.tsv.gz",
)

SOURCE = "SIDER4.1"

RAW_KEYS = [
    f"raw_col_{i}"
    for i in range(7)
]

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
    "sider_drug_id",
    "source_concept_id",
    "source_concept_name",
    "meddra_concept_type",
    "extraction_method",
    "relationship",
    "evidence_text",
    *RAW_KEYS,
    "extraction_evidence_quality",
    "meddra_llt_same_cui_as_pt",
    "drug_name_match_status",
    "sider_source_row_index",
    "sider_source_line_number",
    "sider_disease_name_match_status",
    "sider_record_granularity",
]

for entity in ("drug", "disease"):
    EXTRA += [
        f"{entity}_identifier_{field}"
        for field in (
            "original",
            "namespace",
            "canonical",
            "validation",
        )
    ]

COLS = CORE + EXTRA

REJECT_COLS = [
    "sider_source_row_index",
    "sider_source_line_number",
    "reason",
    *RAW_KEYS,
]

METHODS = {
    "nlp_indication": (
        "label_nlp_indication",
        "sider_nlp_classified_indication",
        "nlp_classified_indication",
    ),
    "nlp_precondition": (
        "label_precondition",
        "sider_nlp_classified_precondition",
        "nlp_classified_precondition",
    ),
    "text_mention": (
        "label_text_mention",
        "sider_unclassified_text_mention",
        "unclassified_text_mention",
    ),
}


def ts():
    return datetime.now(timezone.utc).isoformat()


def section(label):
    print(
        "\n" + "=" * 88
        + "\n" + label
        + "\n" + "=" * 88,
        flush=True,
    )


def hash_file(path):
    h = hashlib.sha256()

    with Path(path).open("rb") as stream:
        for part in iter(
            lambda: stream.read(4 * 1024 * 1024),
            b"",
        ):
            h.update(part)

    return h.hexdigest()


def normal(value):
    return (
        str(value).strip()
        if value is not None
        else ""
    )


def id_info(identifier, entity):
    """
    Syntax validation only.

    No external chemical or disease mapping
    is inferred here.
    """
    original = normal(identifier)
    upper = original.upper()

    if not original:
        return original, "", "", "missing"

    if entity == "drug":
        valid = bool(
            re.fullmatch(
                r"CID[01]\d{8}",
                upper,
            )
        )

        return (
            original,
            "STITCH" if valid else "unresolved",
            "STITCH:" + upper if valid else "",
            "syntax_valid"
            if valid
            else "unrecognised_or_invalid",
        )

    valid = bool(
        re.fullmatch(
            r"C\d{7}",
            upper,
        )
    )

    return (
        original,
        "UMLS" if valid else "unresolved",
        "UMLS:" + upper if valid else "",
        "syntax_valid"
        if valid
        else "unrecognised_or_invalid",
    )


def obtain(
    path,
    download_missing,
    force_download,
):
    if (
        path.is_file()
        and path.stat().st_size
        and not force_download
    ):
        print("[REUSE]", path)
        return

    if not (
        download_missing
        or force_download
    ):
        raise FileNotFoundError(
            f"Required pinned input missing: {path}; "
            "use --download-missing"
        )

    import requests

    path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    url = URL_PREFIX + "/" + path.name

    tmp = path.with_name(
        path.name + ".download"
    )

    try:
        print("[DOWNLOAD]", url)

        with requests.get(
            url,
            stream=True,
            timeout=(30, 600),
        ) as response:
            response.raise_for_status()

            with tmp.open("wb") as file:
                for block in response.iter_content(
                    4 * 1024 * 1024
                ):
                    if block:
                        file.write(block)

        if not tmp.stat().st_size:
            raise RuntimeError(
                "Downloaded source was empty: "
                + path.name
            )

        if path.suffix == ".gz":
            with gzip.open(tmp, "rb") as g:
                if not g.read(1):
                    raise RuntimeError(
                        "Downloaded gzip file was empty"
                    )
        else:
            with tmp.open("rb") as f:
                if not f.read(1):
                    raise RuntimeError(
                        "Downloaded name table was empty"
                    )

        tmp.replace(path)

    finally:
        tmp.unlink(missing_ok=True)


def names_lookup(path):
    """
    Build an ID-to-name mapping only when
    one identifier has an unambiguous name.

    Conflicting names are audited.
    """
    names = defaultdict(set)
    rows = 0
    example = None

    with path.open(
        "r",
        encoding="utf-8-sig",
        newline="",
    ) as file:

        reader = csv.reader(
            file,
            delimiter="\t",
            strict=True,
        )

        for row in reader:
            rows += 1

            if len(row) != 2:
                raise RuntimeError(
                    f"drug_names.tsv record {rows}: "
                    "expected 2 fields, "
                    f"got {len(row)}"
                )

            if example is None:
                example = row

            drug, label = map(normal, row)

            if drug and label:
                names[drug].add(label)

    if not rows:
        raise RuntimeError(
            "drug_names.tsv contains zero records"
        )

    unique = {
        drug: next(iter(values))
        for drug, values in names.items()
        if len(values) == 1
    }

    conflicts = {
        drug: sorted(values)
        for drug, values in names.items()
        if len(values) > 1
    }

    print(
        f"[NAME TABLE] source={rows:,}; "
        f"IDs={len(names):,}; "
        f"unambiguous={len(unique):,}; "
        f"conflicts={len(conflicts):,}"
    )

    print("[NAME RAW EXAMPLE]", example)

    stats = {
        "source_rows": rows,
        "distinct_ids": len(names),
        "conflicting_ids": len(conflicts),
        "conflict_examples": dict(
            list(conflicts.items())[:20]
        ),
    }

    return unique, conflicts, stats


def read_indications(path):
    """
    Strict seven-column headerless TSV reader.

    Preserve each original record and its location.
    """
    records = []
    first = None

    with gzip.open(
        path,
        "rt",
        encoding="utf-8-sig",
        newline="",
    ) as handle:

        reader = csv.reader(
            handle,
            delimiter="\t",
            strict=True,
        )

        for number, cells in enumerate(
            reader,
            1,
        ):
            if len(cells) != 7:
                raise RuntimeError(
                    f"indications source row {number}: "
                    "expected exactly 7 columns, "
                    f"got {len(cells)}"
                )

            if first is None:
                first = cells.copy()

            records.append(
                (
                    number,
                    reader.line_num,
                    cells,
                )
            )

    if not records:
        raise RuntimeError(
            "The indication file has zero rows"
        )

    print(
        "[INDICATION TABLE] source rows:",
        f"{len(records):,}",
    )

    print(
        "[ORIGINAL SCHEMA]",
        RAW_KEYS,
    )

    print(
        "[FIRST ORIGINAL RECORD]",
        dict(zip(RAW_KEYS, first)),
    )

    return records


def make_row(
    index,
    line,
    cells,
    names,
    conflicts,
):
    (
        drug,
        source_id,
        method,
        source_name,
        level,
        disease,
        disease_name,
    ) = map(normal, cells)

    missing = [
        key
        for key, field in (
            ("drug_id", drug),
            ("disease_id", disease),
        )
        if not field
    ]

    if missing:
        return None, (
            dict(zip(RAW_KEYS, cells))
            | {
                "sider_source_row_index": str(index),
                "sider_source_line_number": str(line),
                "reason": "|".join(
                    "missing_" + key
                    for key in missing
                ),
            }
        )

    # Do not choose arbitrary names
    # where an identifier has conflicts.

    if drug in conflicts:
        drug_label = drug
        match = (
            "ambiguous_names_identifier_placeholder"
        )

    elif drug in names:
        drug_label = names[drug]
        match = "matched"

    else:
        drug_label = drug
        match = (
            "unmatched_identifier_placeholder"
        )

    disease_label = disease_name or disease

    disease_status = (
        "provided"
        if disease_name
        else "identifier_placeholder"
    )

    (
        relation,
        evidence,
        quality,
    ) = METHODS.get(
        method.casefold(),
        (
            "label_unclassified",
            "sider_unrecognised_extraction_method",
            "unrecognised_method"
            if method
            else "missing_method",
        ),
    )

    drug_id = id_info(
        drug,
        "drug",
    )

    disease_id = id_info(
        disease,
        "disease",
    )

    row = {
        "drug_name": drug_label,
        "drug_identifier": drug,
        "drug_identifier_type":
            drug_id[1] or "unresolved",

        "disease_or_condition_name":
            disease_label,

        "disease_or_condition_identifier":
            disease,

        "disease_or_condition_identifier_type":
            disease_id[1] or "unresolved",

        "relationship_type": relation,
        "evidence_type": evidence,

        "source": SOURCE,
        "internal_source":
            "SIDER_meddra_all_indications",

        "sider_drug_id": drug,
        "source_concept_id": source_id,
        "source_concept_name": source_name,
        "meddra_concept_type": level,
        "extraction_method": method,
        "relationship": relation,

        "evidence_text": (
            f"source_concept_id={source_id}; "
            f"extraction_method={method}; "
            f"MedDRA_level={level}"
        ),

        **{
            f"raw_col_{i}": cells[i]
            for i in range(7)
        },

        "extraction_evidence_quality": quality,
        "meddra_llt_same_cui_as_pt": "False",
        "drug_name_match_status": match,
        "sider_source_row_index": str(index),
        "sider_source_line_number": str(line),

        "sider_disease_name_match_status":
            disease_status,

        "sider_record_granularity":
            "original_SIDER_indication_file_row",
    }

    for entity, vals in (
        ("drug", drug_id),
        ("disease", disease_id),
    ):
        for field, val in zip(
            (
                "original",
                "namespace",
                "canonical",
                "validation",
            ),
            vals,
        ):
            row[
                f"{entity}_identifier_{field}"
            ] = val

    return row, None


def deduplicated(df):
    """
    Assertion-level reporting view.

    All original source rows remain
    in the full CSV.

    Different source concepts, methods
    and MedDRA levels remain distinct.
    """
    keys = [
        "drug_identifier",
        "disease_or_condition_identifier",
        "disease_or_condition_name",
        "source_concept_id",
        "source_concept_name",
        "meddra_concept_type",
        "extraction_method",
        "relationship_type",
        "evidence_type",
    ]

    groups = df.groupby(
        keys,
        sort=False,
        dropna=False,
    )

    output = (
        groups.first()
        .reset_index()
        .reindex(columns=COLS)
    )

    support = (
        groups["sider_source_row_index"]
        .agg("|".join)
        .rename("sider_supporting_source_rows")
        .reset_index()
    )

    counts = (
        groups.size()
        .rename("sider_supporting_record_count")
        .reset_index()
    )

    output = output.merge(
        support,
        on=keys,
        validate="one_to_one",
    )

    output = output.merge(
        counts,
        on=keys,
        validate="one_to_one",
    )

    if int(
        output[
            "sider_supporting_record_count"
        ].sum()
    ) != len(df):
        raise RuntimeError(
            "Deduplicated support count does not "
            "match full records"
        )

    return output, keys


def write_staged(frame, dest):
    """
    Stage output before overwriting
    the original destination.
    """
    dest.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    temp = dest.with_name(
        dest.name + ".processing.tmp"
    )

    frame.to_csv(
        temp,
        index=False,
    )

    with temp.open(
        "r",
        encoding="utf-8",
        newline="",
    ) as f:
        if next(csv.reader(f)) != list(frame.columns):
            raise RuntimeError(
                f"Output header mismatch: {dest}"
            )

    return temp


def main():
    p = argparse.ArgumentParser(
        description=__doc__
    )

    p.add_argument(
        "--workdir",
        type=Path,
        default=RAW,
    )

    p.add_argument(
        "--out",
        type=Path,
        default=OUT / "sider_drug_disease.csv",
    )

    p.add_argument(
        "--dedup_out",
        type=Path,
        default=(
            OUT /
            "sider_drug_disease_deduplicated.csv"
        ),
    )

    p.add_argument(
        "--pt_out",
        type=Path,
        default=OUT / "sider_drug_disease_pt_only.csv",
    )

    p.add_argument(
        "--unmatched_out",
        type=Path,
        default=OUT / "sider_unmatched_drug_ids.csv",
    )

    p.add_argument(
        "--metadata",
        type=Path,
        default=OUT / "sider_summary.json",
    )

    p.add_argument(
        "--download-missing",
        action="store_true",
    )

    p.add_argument(
        "--force_download",
        "--force-download",
        action="store_true",
    )

    args = p.parse_args()

    raw = args.workdir.resolve()

    paths = {
        name: raw / name
        for name in FILE_NAMES
    }

    # ========================================================
    # 1. FILE INVENTORY
    # ========================================================

    section("1. RAW INPUT INVENTORY")

    for path in paths.values():
        obtain(
            path,
            args.download_missing,
            args.force_download,
        )

    inventory = []

    for path in sorted(raw.rglob("*")):
        if path.is_file():
            info = {
                "path": str(
                    path.relative_to(raw)
                ),
                "bytes": path.stat().st_size,
                "input": (
                    path.resolve()
                    in {
                        x.resolve()
                        for x in paths.values()
                    }
                ),
            }

            inventory.append(info)

            print(
                "[INPUT]"
                if info["input"]
                else "[OTHER]",
                info,
            )

    hashes = {
        name: hash_file(path)
        for name, path in paths.items()
    }

    # ========================================================
    # 2. ORIGINAL SCHEMA AND ROW COUNTS
    # ========================================================

    section(
        "2. ORIGINAL SCHEMAS "
        "AND RECORD COUNTS"
    )

    names, conflicts, names_stats = (
        names_lookup(
            paths["drug_names.tsv"]
        )
    )

    sources = read_indications(
        paths["meddra_all_indications.tsv.gz"]
    )

    # ========================================================
    # 3. SCIENTIFIC MAPPING RULES
    # ========================================================

    section(
        "3. MERGE RULE / "
        "EVIDENCE INTERPRETATION"
    )

    print(
        "Join drug_names TSV column 0 to "
        "indications TSV column 0 "
        "by exact STITCH compound ID."
    )

    print(
        "Seven original indication columns "
        "retained; every source event "
        "has a row index."
    )

    print(
        "NLP_indication = label NLP classification."
    )

    print(
        "NLP_precondition = label precondition."
    )

    print(
        "text_mention = unclassified label mention."
    )

    print(
        "Unknown methods = unclassified."
    )

    print(
        "Source records are NOT FDA approvals "
        "or validated therapeutic efficacy."
    )

    # ========================================================
    # 4. TRANSFORM
    # ========================================================

    out_rows = []
    rejects = []

    for row_index, line, cells in sources:
        record, bad = make_row(
            row_index,
            line,
            cells,
            names,
            conflicts,
        )

        if bad is None:
            out_rows.append(record)
        else:
            rejects.append(bad)

    full = pd.DataFrame(
        out_rows,
        columns=COLS,
    )

    rejected = pd.DataFrame(
        rejects,
        columns=REJECT_COLS,
    )

    if (
        full.empty
        or len(full) + len(rejected) != len(sources)
    ):
        raise RuntimeError(
            "No valid evidence or source "
            "row accounting failed"
        )

    # Identify LLT records with a matching
    # PT-level CUI, without removing them.

    pt_rows = full.loc[
        full["meddra_concept_type"].eq("PT")
    ]

    pt_keys = set(zip(
        pt_rows["drug_identifier"],
        pt_rows["source_concept_id"],
        pt_rows["extraction_method"],
        pt_rows["disease_or_condition_identifier"],
    ))

    full["meddra_llt_same_cui_as_pt"] = [
        str(
            level == "LLT"
            and (
                drug,
                source,
                method,
                disease,
            ) in pt_keys
        )
        for (
            drug,
            source,
            method,
            disease,
            level,
        ) in full[
            [
                "drug_identifier",
                "source_concept_id",
                "extraction_method",
                "disease_or_condition_identifier",
                "meddra_concept_type",
            ]
        ].itertuples(
            index=False,
            name=None,
        )
    ]

    dedup, dedup_keys = deduplicated(full)

    # Existing PT-only view retains
    # every extraction-method category.

    pt_only = full.loc[
        full["meddra_concept_type"].eq("PT")
    ].copy()

    # This additional view is more conservative:
    # only NLP-classified indication records
    # with PT-level disease concepts.

    nlp_pt = pt_only.loc[
        pt_only[
            "extraction_method"
        ].str.casefold().eq("nlp_indication")
    ].copy()

    unmatched = (
        full.loc[
            full["drug_name_match_status"].ne("matched")
        ]
        .groupby(
            [
                "drug_identifier",
                "drug_name_match_status",
            ],
            sort=False,
        )
        .size()
        .rename("n_indication_rows")
        .reset_index()
    )

    # ========================================================
    # 5. FINAL COUNTS AND QUALITY CHECKS
    # ========================================================

    section(
        "4. ROW CONSERVATION, "
        "METHODS, AND OUTPUT SCHEMA"
    )

    print(
        "Source indication rows:",
        len(sources),
    )

    print(
        "Full output rows:",
        len(full),
        "Rejected rows:",
        len(rejected),
    )

    print(
        "Deduplicated assertion rows:",
        len(dedup),
    )

    print(
        "PT rows:",
        len(pt_only),
    )

    print(
        "Strict NLP indication PT rows:",
        len(nlp_pt),
    )

    print(
        "Method counts:",
        full["extraction_method"]
        .value_counts()
        .to_dict(),
    )

    print(
        "Relationship counts:",
        full["relationship_type"]
        .value_counts()
        .to_dict(),
    )

    print(
        "MedDRA levels:",
        full["meddra_concept_type"]
        .value_counts()
        .to_dict(),
    )

    print(
        "Unmatched drug-name rows:",
        int(
            full[
                "drug_name_match_status"
            ].ne("matched").sum()
        ),
    )

    print(
        "Drug ID validation:",
        full["drug_identifier_validation"]
        .value_counts()
        .to_dict(),
    )

    print(
        "Disease ID validation:",
        full["disease_identifier_validation"]
        .value_counts()
        .to_dict(),
    )

    print(
        "Core fields:",
        " | ".join(CORE),
    )

    print(
        "Full output fields:",
        " | ".join(COLS),
    )

    print(
        "First mapped record:\n",
        json.dumps(
            full.iloc[0].to_dict(),
            ensure_ascii=False,
            indent=2,
        )[:2800],
    )

    # ========================================================
    # 6. OUTPUT FILES
    # ========================================================

    outputs = {
        args.out.resolve(): full,

        args.dedup_out.resolve(): dedup,

        args.pt_out.resolve(): pt_only,

        args.unmatched_out.resolve(): unmatched,

        args.out.resolve().with_name(
            "sider_drug_disease_nlp_indication_pt.csv"
        ): nlp_pt,

        args.out.resolve().with_name(
            "sider_rejected_rows.csv"
        ): rejected,
    }

    if len(outputs) != 6:
        raise RuntimeError(
            "Duplicate configured output file paths"
        )

    summary = {
        "status": "success",
        "timestamp_utc": ts(),
        "source": SOURCE,

        "source_repository":
            "https://sideeffects.embl.de/download/",

        "license_note": (
            "Check primary SIDER distribution "
            "and MedDRA terms before redistribution; "
            "do not assume unrestricted "
            "downstream licensing."
        ),

        "source_files": inventory,
        "source_sha256": hashes,

        "name_table_stats": names_stats,
        "indications_raw_columns": RAW_KEYS,

        "source_rows": len(sources),
        "output_rows": len(full),
        "rejected_rows": len(rejected),
        "deduplicated_rows": len(dedup),

        "pt_only_rows": len(pt_only),
        "nlp_indication_pt_rows": len(nlp_pt),

        "drug_name_match_status":
            full["drug_name_match_status"]
            .value_counts()
            .to_dict(),

        "extraction_method_counts":
            full["extraction_method"]
            .value_counts()
            .to_dict(),

        "relationship_type_counts":
            full["relationship_type"]
            .value_counts()
            .to_dict(),

        "meddra_level_counts":
            full["meddra_concept_type"]
            .value_counts()
            .to_dict(),

        "llt_rows_same_cui_as_pt":
            int(
                full[
                    "meddra_llt_same_cui_as_pt"
                ].eq("True").sum()
            ),

        "distinct_compound_disease_pairs":
            len(
                full[
                    [
                        "drug_identifier",
                        "disease_or_condition_identifier",
                    ]
                ].drop_duplicates()
            ),

        "drug_identifier_validation":
            full["drug_identifier_validation"]
            .value_counts()
            .to_dict(),

        "disease_identifier_validation":
            full["disease_identifier_validation"]
            .value_counts()
            .to_dict(),

        "output_columns": COLS,
        "dedup_keys": dedup_keys,

        "merge_instructions": (
            "Use full CSV with relationship_type "
            "and evidence_type; never treat "
            "preconditions or text mentions "
            "as therapeutic indications. "
            "Use the NLP_indication PT view "
            "for conservative exploratory analysis."
        ),

        "name_join": (
            "Exact STITCH CID key: "
            "indications.raw_col_0 "
            "-> drug_names.column_0"
        ),

        "pt_only_view_semantics": (
            "ALL extraction methods, PT level; "
            "NOT a strict indication-only dataset"
        ),
    }

    staged = {}

    try:
        # Stage every output before replacement.

        for dest, frame in outputs.items():
            staged[dest] = write_staged(
                frame,
                dest,
            )

        summary["output_files_sha256"] = {
            str(dest): hash_file(tmp)
            for dest, tmp in staged.items()
        }

        meta_path = args.metadata.resolve()

        if meta_path in staged:
            raise RuntimeError(
                "Metadata path overlaps a CSV"
            )

        meta_path.parent.mkdir(
            parents=True,
            exist_ok=True,
        )

        meta_tmp = meta_path.with_name(
            meta_path.name + ".processing.tmp"
        )

        meta_tmp.write_text(
            json.dumps(
                summary,
                indent=2,
                ensure_ascii=False,
            ) + "\n",
            encoding="utf-8",
        )

        staged[meta_path] = meta_tmp

        section(
            "5. REPLACING FILES "
            "IN ORIGINAL PROCESSED DIRECTORY"
        )

        for dest, tmp in staged.items():
            tmp.replace(dest)

            print(
                "[OVERWRITTEN]",
                dest,
            )

        print(
            "SUCCESS: use sider_drug_disease.csv "
            "for the evidence-preserving merge."
        )

        print(
            "For conservative candidate indications, "
            "use sider_drug_disease_nlp_indication_pt.csv."
        )

    finally:
        for tmp in staged.values():
            tmp.unlink(missing_ok=True)


if __name__ == "__main__":
    main()
