
#!/usr/bin/env python3
"""MEDI-2 original ICD/UMLS files -> merge-compatible medication-indication evidence.

Run: python MEDI.py
Reuses data/raw/medi/MEDI-2.csv and MEDI-2_UMLS.csv; no network by default.
Overwrites original data/processed/medi/ outputs, without processed_v2.
The main CSV contains all eligible source records; dedup view retains support counts.
MEDI annotations are not evidence of regulatory approval or independent validation.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import urllib.request
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

ROOT = Path(
    os.environ.get(
        "DRUG_DISEASE_ROOT",
        Path(__file__).resolve().parent,
    )
).resolve()

RAW = ROOT / "data/raw/medi"
OUT = ROOT / "data/processed/medi"

URLS = {
    "MEDI-2.csv":
        "https://www.vumc.org/wei-lab/sites/default/files/public_files/MEDI-2.csv",
    "MEDI-2_UMLS.csv":
        "https://www.vumc.org/wei-lab/sites/default/files/public_files/MEDI-2_UMLS.csv",
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
    "high_precision_subset",
    "high_precision_subset_raw",
    "number_of_resources",
    "resources_mentioning",
    "medi_source_flags",
    "medi_confidence_tier",
    "medi_low_confidence_flag",
    "medi_indication_text_key",
    "source_record_type",
    "medi_source_file",
    "medi_source_row_number",
    "medi_source_record_json",
    "evidence_text",
    "medi_original_vocabulary",
    "medi_support_count_used",
    "medi_support_count_basis",
    "medi_source_semantics",
]

for entity in ("drug", "disease"):
    EXTRA += [
        f"{entity}_identifier_{x}"
        for x in (
            "original",
            "namespace",
            "canonical",
            "validation",
        )
    ]

COLUMNS = CORE + EXTRA

ICD_REQUIRED = [
    "RXCUI",
    "DRUG_DESCRIPTION",
    "VOCABULARY",
    "CODE",
    "INDICATION_DESCRIPTION",
]

UMLS_REQUIRED = [
    "RXCUI",
    "DRUG_DESC",
    "CUI",
    "STR",
]

HPS = "HIGH_PRECISION_SUBSET"
N_RES = "NUMBER_OF_RESOURCES_MENTIONED"

TRUE = {"TRUE", "T", "1", "YES", "Y"}
FALSE = {"FALSE", "F", "0", "NO", "N"}
NULL_TOKENS = {"", "NA", "N/A", "NULL", "NONE", "NAN"}


def section(title):
    print(
        "\n" + "=" * 90
        + "\n" + title
        + "\n" + "=" * 90,
        flush=True,
    )


def utc():
    return datetime.now(timezone.utc).isoformat()


def checksum(path):
    h = hashlib.sha256()

    with open(path, "rb") as f:
        for b in iter(
            lambda: f.read(4 * 1024 * 1024),
            b"",
        ):
            h.update(b)

    return h.hexdigest()


def norm_text(value):
    return re.sub(
        r"[^a-z0-9]+",
        " ",
        str(value).casefold(),
    ).strip()


def blank(value):
    return str(value).strip().upper() in NULL_TOKENS


def truth(value):
    s = str(value).strip().upper()

    if s in TRUE:
        return True

    if s in FALSE:
        return False

    return None


def confidence(value, present):
    if not present:
        return (
            "",
            "not_reported",
            "medi_indication_precision_not_reported",
        )

    result = truth(value)

    if result is True:
        return (
            "True",
            "high_precision",
            "medi_high_precision_indication",
        )

    if result is False:
        return (
            "False",
            "ensemble_only",
            "medi_ensemble_indication",
        )

    return (
        "",
        "not_reported",
        "medi_indication_precision_not_reported",
    )


def namespace_of(vocab):
    """
    Conservative namespace assignment.
    ICD9 is not assumed to mean ICD9CM.
    """
    v = re.sub(
        r"[ _-]",
        "",
        str(vocab).upper(),
    )

    return {
        "ICD9CM": "ICD9CM",
        "ICD10CM": "ICD10CM",
        "ICD9": "ICD9",
        "ICD10": "ICD10",
        "UMLS": "UMLS",
        "SNOMEDCT": "SNOMEDCT",
    }.get(v, "unresolved")


def id_fields(identifier, ns):
    """
    Return original ID, namespace,
    optional canonical CURIE and validation status.

    Syntax validation does not prove a code is active
    or establish equivalent ontology concepts.
    """
    raw = str(identifier).strip()

    if blank(raw):
        return "", ns, "", "missing"

    up = raw.upper()

    if ns == "RxNorm":
        good = bool(re.fullmatch(r"\d+", raw))

        return (
            raw,
            ns,
            f"RXNORM:{raw}" if good else "",
            "syntax_valid" if good else "invalid_syntax",
        )

    if ns == "UMLS":
        good = bool(
            re.fullmatch(r"C\d{7}", up)
        )

        return (
            raw,
            ns,
            f"UMLS:{up}" if good else "",
            "syntax_valid" if good else "invalid_syntax",
        )

    if ns in {"ICD9CM", "ICD10CM"}:
        if ns == "ICD9CM":
            pat = (
                r"(?:"
                r"\d{3}(?:\.\d{1,2})?"
                r"|V\d{2}(?:\.\d{1,2})?"
                r"|E\d{3}(?:\.\d)?"
                r")"
            )
        else:
            pat = (
                r"[A-Z]\d[A-Z0-9]"
                r"(?:\.[A-Z0-9]{1,4})?"
            )

        good = bool(re.fullmatch(pat, up))

        return (
            raw,
            ns,
            f"{ns}:{up}" if good else "",
            "syntax_valid" if good else "invalid_syntax",
        )

    return raw, ns, "", "unverified_namespace"


def write_temp_csv(df, destination):
    destination.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    temp = destination.with_name(
        destination.name + ".processing.tmp"
    )

    df.to_csv(temp, index=False)

    actual = list(
        pd.read_csv(temp, nrows=0).columns
    )

    if actual != list(df.columns):
        raise RuntimeError(
            f"Unexpected written columns: {temp}"
        )

    return temp


def write_temp_json(content, destination):
    destination.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    temp = destination.with_name(
        destination.name + ".processing.tmp"
    )

    temp.write_text(
        json.dumps(
            content,
            indent=2,
            ensure_ascii=False,
            default=str,
        ) + "\n",
        encoding="utf-8",
    )

    return temp


def obtain_file(path, url, download_missing):
    if path.is_file() and path.stat().st_size:
        print(f"[REUSE] {path}")
        return

    if not download_missing:
        raise FileNotFoundError(
            f"Missing source file {path}; "
            "use --download_missing intentionally"
        )

    path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    temp = path.with_name(
        path.name + ".download"
    )

    try:
        urllib.request.urlretrieve(
            url,
            temp,
        )

        if not temp.stat().st_size:
            raise RuntimeError(
                f"Empty downloaded source: {url}"
            )

        temp.replace(path)

    finally:
        temp.unlink(missing_ok=True)


def read_source(path, required):
    df = pd.read_csv(
        path,
        dtype=str,
        keep_default_na=False,
        encoding="utf-8-sig",
        on_bad_lines="error",
    )

    missing = sorted(
        set(required) - set(df.columns)
    )

    if missing:
        raise RuntimeError(
            f"{path.name} missing required columns: "
            f"{missing}; observed {list(df.columns)}"
        )

    if df.columns.duplicated().any():
        raise RuntimeError(
            f"Duplicate headers in {path.name}"
        )

    print(
        f"{path.name}: {len(df):,} source records "
        f"| {len(df.columns)} fields"
    )

    print(
        "FULL HEADER:",
        " | ".join(df.columns),
    )

    if len(df):
        print("FIRST RAW RECORD:")

        for k, v in df.iloc[0].items():
            print(
                f"  {k:34} "
                f"{str(v)[:100]}"
            )

    return df


def process_source(df, file_name, mode):
    section(f"3. PROCESS {file_name}")

    configs = {
        "ICD-coded": (
            "DRUG_DESCRIPTION",
            "CODE",
            "INDICATION_DESCRIPTION",
            "VOCABULARY",
        ),
        "UMLS-coded": (
            "DRUG_DESC",
            "CUI",
            "STR",
            None,
        ),
    }

    drug_col, code_col, disease_col, vocab_col = (
        configs[mode]
    )

    support_cols = [
        c
        for c in df.columns
        if c.startswith("IN_")
        and c not in {HPS, N_RES}
    ]

    print(
        "SUPPORT FLAGS (kept unaltered):",
        support_cols,
    )

    unknown_flags = {
        c: sorted(
            set(
                df[c].str.upper().str.strip()
            ) - TRUE - FALSE - NULL_TOKENS
        )[:20]
        for c in support_cols
    }

    unknown_flags = {
        k: v
        for k, v in unknown_flags.items()
        if v
    }

    if unknown_flags:
        raise RuntimeError(
            "Unexpected nonboolean source-support "
            f"flag values in {file_name}: "
            f"{unknown_flags}"
        )

    hps_exists = HPS in df.columns

    rows = []
    rejected = []

    confidence_count = Counter()
    id_validation = Counter()

    for source_row_number, record in enumerate(
        df.to_dict("records"),
        1,
    ):
        drug = str(
            record[drug_col]
        ).strip()

        rxcui = str(
            record["RXCUI"]
        ).strip()

        disease = str(
            record[disease_col]
        ).strip()

        code = str(
            record[code_col]
        ).strip()

        vocab = (
            str(record[vocab_col]).strip()
            if vocab_col
            else "UMLS"
        )

        missing = [
            field
            for field, value in (
                ("drug_name", drug),
                ("rxcui", rxcui),
                ("disease_name", disease),
                ("disease_code", code),
            )
            if blank(value)
        ]

        if missing:
            rejected.append({
                "source_file": file_name,
                "source_row_number": source_row_number,
                "rejection_reason":
                    "missing:" + "|".join(missing),
                "original_record_json":
                    json.dumps(
                        record,
                        ensure_ascii=False,
                    ),
            })

            continue

        raw_hps = (
            str(
                record.get(HPS, "")
            ).strip()
            if hps_exists
            else ""
        )

        hp, tier, evidence = confidence(
            raw_hps,
            hps_exists,
        )

        ns = namespace_of(vocab)

        drug_id = id_fields(
            rxcui,
            "RxNorm",
        )

        disease_id = id_fields(
            code,
            ns,
        )

        support_true = [
            flag
            for flag in support_cols
            if truth(record[flag]) is True
        ]

        reported_support = str(
            record.get(N_RES, "")
        ).strip()

        if (
            reported_support
            and not blank(reported_support)
        ):
            if not re.fullmatch(
                r"\d+",
                reported_support,
            ):
                raise RuntimeError(
                    f"Noninteger {N_RES} at "
                    f"{file_name}:{source_row_number}: "
                    f"{reported_support!r}"
                )

            support_count = int(
                reported_support
            )

            support_basis = N_RES

        elif support_cols:
            support_count = len(
                support_true
            )

            support_basis = "observed_IN_flags"

        else:
            support_count = None

            support_basis = "not_available"

        low = (
            tier == "ensemble_only"
            and support_count is not None
            and support_count <= 1
        )

        low_label = (
            "Unknown"
            if support_count is None
            else str(low)
        )

        flags = "|".join(
            f"{c}={str(record[c]).strip()}"
            for c in support_cols
        )

        out = dict(
            drug_name=drug,
            drug_identifier=rxcui,
            drug_identifier_type="RxCUI",
            disease_or_condition_name=disease,
            disease_or_condition_identifier=code,
            disease_or_condition_identifier_type=ns,
            relationship_type="indication",
            evidence_type=evidence,
            source="MEDI",
            internal_source=f"MEDI-2 {mode}",
            high_precision_subset=hp,
            high_precision_subset_raw=raw_hps,
            number_of_resources=reported_support,
            resources_mentioning="|".join(
                support_true
            ),
            medi_source_flags=flags,
            medi_confidence_tier=tier,
            medi_low_confidence_flag=low_label,
            medi_indication_text_key=(
                f"{rxcui}|{norm_text(disease)}"
            ),
            source_record_type=mode,
            medi_source_file=file_name,
            medi_source_row_number=str(
                source_row_number
            ),
            medi_source_record_json=json.dumps(
                record,
                ensure_ascii=False,
            ),
            evidence_text=(
                f"{mode}; high_precision={tier}; "
                f"source_support="
                f"{support_count if support_count is not None else 'unknown'}"
            ),
            medi_original_vocabulary=vocab,
            medi_support_count_used=(
                ""
                if support_count is None
                else str(support_count)
            ),
            medi_support_count_basis=support_basis,
            medi_source_semantics=(
                "MEDI ensemble medication-indication "
                "annotation; not regulatory approval"
            ),
        )

        for entity, id_data in (
            ("drug", drug_id),
            ("disease", disease_id),
        ):
            for suffix, value in zip(
                (
                    "original",
                    "namespace",
                    "canonical",
                    "validation",
                ),
                id_data,
            ):
                out[
                    f"{entity}_identifier_{suffix}"
                ] = value

        rows.append(out)

        confidence_count[tier] += 1

        id_validation[
            (ns, disease_id[3])
        ] += 1

    result = pd.DataFrame(
        rows,
        columns=COLUMNS,
    )

    rejects = pd.DataFrame(
        rejected,
        columns=[
            "source_file",
            "source_row_number",
            "rejection_reason",
            "original_record_json",
        ],
    )

    if len(result) + len(rejects) != len(df):
        raise RuntimeError(
            f"Conservation failure in {file_name}"
        )

    print(
        f"Source rows: {len(df):,} | "
        f"retained: {len(result):,} | "
        f"audited rejects: {len(rejects):,}"
    )

    print(
        "Precision tiers:",
        dict(confidence_count),
    )

    print(
        "Vocabulary/validation:",
        {
            f"{k[0]}:{k[1]}": v
            for k, v in id_validation.items()
        },
    )

    stats = {
        "input_rows": len(df),
        "accepted": len(result),
        "rejected": len(rejects),
        "original_header": list(df.columns),
        "resource_flag_columns": support_cols,
        "precision_tiers":
            dict(confidence_count),
    }

    return result, rejects, stats


def audited_dedup(out):
    """
    Retain one assertion per identical metadata key,
    together with every supporting source-row locator.

    The full source evidence CSV is never deduplicated.
    """
    keys = [
        "drug_identifier",
        "drug_name",
        "disease_or_condition_identifier",
        "disease_or_condition_identifier_type",
        "disease_or_condition_name",
        "relationship_type",
        "evidence_type",
        "source_record_type",
        "high_precision_subset_raw",
        "number_of_resources",
        "medi_source_flags",
    ]

    annotated = out.copy()

    annotated["_source_locator"] = (
        annotated["medi_source_file"]
        + ":"
        + annotated["medi_source_row_number"]
    )

    grouping = annotated.groupby(
        keys,
        dropna=False,
        sort=False,
    )

    dedup = grouping.first().reset_index()

    links = (
        grouping["_source_locator"]
        .agg("|".join)
        .rename("medi_supporting_source_rows")
        .reset_index()
    )

    counts = (
        grouping.size()
        .rename("medi_supporting_record_count")
        .reset_index()
    )

    dedup = dedup.merge(
        links,
        on=keys,
        validate="one_to_one",
    )

    dedup = dedup.merge(
        counts,
        on=keys,
        validate="one_to_one",
    )

    dedup = dedup.reindex(
        columns=COLUMNS + [
            "medi_supporting_record_count",
            "medi_supporting_source_rows",
        ]
    )

    if int(
        dedup[
            "medi_supporting_record_count"
        ].sum()
    ) != len(out):
        raise RuntimeError(
            "Deduplication provenance "
            "reconciliation failed"
        )

    return dedup


def redundancy_audit(out, crosswalk=None):
    """
    Measure lexical cross-code overlap.

    A matching description alone does NOT establish
    equivalence between ICD and UMLS concepts.
    """
    group = out.groupby(
        "medi_indication_text_key"
    )["disease_identifier_namespace"].agg(
        lambda x: sorted(set(x))
    )

    result = {
        "same_drug_and_text_groups":
            len(group),

        "text_groups_with_multiple_coding_namespaces_unverified":
            int(
                (group.map(len) > 1).sum()
            ),

        "interpretation": (
            "Text-based overlap is not verified "
            "concept equivalence and not "
            "independent evidence"
        ),

        "crosswalk": "not supplied",
    }

    if not crosswalk:
        return result

    cross = pd.read_csv(
        crosswalk,
        dtype=str,
        keep_default_na=False,
    )

    required = {
        "source_curie",
        "target_curie",
        "relation",
        "verified",
    }

    if not required.issubset(cross.columns):
        raise RuntimeError(
            f"Crosswalk needs columns "
            f"{sorted(required)}"
        )

    if "entity" in cross.columns:
        cross = cross.loc[
            cross["entity"].str.casefold()
            == "disease"
        ]

    cross = cross.loc[
        cross["verified"].str.upper().isin(TRUE)
        &
        cross["relation"].str.casefold().isin(
            {
                "exact",
                "equivalent",
                "exact_match",
                "skos:exactmatch",
            }
        )
    ]

    cgroup = cross.groupby(
        "source_curie"
    )["target_curie"].agg(
        lambda x: set(x)
    )

    ambiguous = {
        key
        for key, v in cgroup.items()
        if len(v) != 1
    }

    trusted = {
        key: next(iter(v))
        for key, v in cgroup.items()
        if len(v) == 1
    }

    umls = set(
        zip(
            out.loc[
                out["disease_identifier_namespace"]
                == "UMLS",
                "drug_identifier",
            ],
            out.loc[
                out["disease_identifier_namespace"]
                == "UMLS",
                "disease_identifier_canonical",
            ],
        )
    )

    icd = out.loc[
        out["disease_identifier_namespace"].isin(
            {
                "ICD9CM",
                "ICD10CM",
                "ICD9",
                "ICD10",
            }
        )
    ]

    matching = sum(
        (drug, trusted.get(source, "")) in umls
        for drug, source in zip(
            icd["drug_identifier"],
            icd["disease_identifier_canonical"],
        )
    )

    result["crosswalk"] = {
        "path": str(crosswalk),

        "sha256": checksum(crosswalk),

        "verified_one_to_one_source_mappings":
            len(trusted),

        "ambiguous_source_mappings_omitted":
            len(ambiguous),

        "ICD_rows_with_explicit_exact_UMLS_equivalent_for_same_RxCUI":
            int(matching),

        "note": (
            "Audit only: no source records are "
            "merged or overwritten by crosswalk"
        ),
    }

    return result


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
    )

    ap.add_argument(
        "--raw_dir",
        default=str(RAW),
    )

    ap.add_argument(
        "--output_dir",
        default=str(OUT),
    )

    ap.add_argument(
        "--download_missing",
        action="store_true",
    )

    ap.add_argument(
        "--crosswalk",
        default="",
        help=(
            "Optional reviewed exact-equivalence "
            "crosswalk (audit only)"
        ),
    )

    args = ap.parse_args()

    raw = Path(args.raw_dir).resolve()
    destination = Path(
        args.output_dir
    ).resolve()

    paths = {
        name: raw / name
        for name in URLS
    }

    # ========================================================
    # 1. RAW FILE INVENTORY
    # ========================================================

    section("1. SOURCE FILE INVENTORY")

    for name, url in URLS.items():
        obtain_file(
            paths[name],
            url,
            args.download_missing,
        )

    inventory = []

    for p in sorted(raw.rglob("*")):
        if p.is_file():
            info = {
                "name": str(
                    p.relative_to(raw)
                ),
                "bytes": p.stat().st_size,
                "used": (
                    p.resolve()
                    in {
                        x.resolve()
                        for x in paths.values()
                    }
                ),
            }

            inventory.append(info)

            print(
                f"  [{'USED' if info['used'] else 'OTHER'}] "
                f"{info['name']} | "
                f"{info['bytes']:,} bytes"
            )

    print(
        f"Physical raw files: {len(inventory)} "
        "| required processing files: 2"
    )

    # ========================================================
    # 2. ACTUAL INPUT SCHEMAS
    # ========================================================

    section(
        "2. INPUT SCHEMAS, ORIGINAL "
        "HEADERS AND EXAMPLES"
    )

    icd = read_source(
        paths["MEDI-2.csv"],
        ICD_REQUIRED,
    )

    umls = read_source(
        paths["MEDI-2_UMLS.csv"],
        UMLS_REQUIRED,
    )

    section(
        "DATA MAPPING / JOIN RULES"
    )

    print(
        "No SQL join: ICD-coded and UMLS-coded "
        "files are transformed separately "
        "and concatenated vertically."
    )

    print(
        "ICD: RXCUI, DRUG_DESCRIPTION, "
        "VOCABULARY, CODE, "
        "INDICATION_DESCRIPTION "
        "-> common fields."
    )

    print(
        "UMLS: RXCUI, DRUG_DESC, CUI, STR "
        "-> common fields."
    )

    print(
        "Original vocabulary and source_record_type "
        "are retained. ICD and UMLS records are "
        "NOT independent confirmations."
    )

    print(
        "HIGH_PRECISION_SUBSET describes confidence, "
        "not a different biological relationship "
        "or proof of approval."
    )

    # ========================================================
    # 3. PROCESS EACH ORIGINAL FILE
    # ========================================================

    x, xbad, xstats = process_source(
        icd,
        "MEDI-2.csv",
        "ICD-coded",
    )

    y, ybad, ystats = process_source(
        umls,
        "MEDI-2_UMLS.csv",
        "UMLS-coded",
    )

    # ========================================================
    # 4. COMBINE WITHOUT DISCARDING SOURCE EVIDENCE
    # ========================================================

    out = pd.concat(
        [x, y],
        axis=0,
        ignore_index=True,
    ).reindex(columns=COLUMNS)

    rejects = pd.concat(
        [xbad, ybad],
        axis=0,
        ignore_index=True,
    )

    if (
        out.empty
        or len(out) + len(rejects)
        != len(icd) + len(umls)
    ):
        raise RuntimeError(
            "Combined source conservation failure; "
            "refusing overwrite"
        )

    dedup = audited_dedup(out)

    low = out.loc[
        out["medi_low_confidence_flag"]
        == "True"
    ].copy()

    redundancy = redundancy_audit(
        out,
        args.crosswalk or None,
    )

    # ========================================================
    # 5. HUMAN-READABLE SUMMARY
    # ========================================================

    section(
        "4. COUNTS / SCHEMA / RECONCILIATION"
    )

    print(
        f"Original ICD + UMLS: "
        f"{len(icd):,} + {len(umls):,} "
        f"= {len(icd) + len(umls):,}"
    )

    print(
        f"Accepted evidence events: {len(out):,} "
        f"| audited rejected rows: {len(rejects):,}"
    )

    print(
        f"Assertion-level records: {len(dedup):,} "
        f"| repeated evidence rows: "
        f"{len(out) - len(dedup):,}"
    )

    print(
        f"Low-confidence rows: {len(low):,}"
    )

    print(
        "Precision tiers:",
        out["medi_confidence_tier"]
        .value_counts()
        .to_dict(),
    )

    print(
        "Disease vocabularies:",
        out[
            "disease_or_condition_identifier_type"
        ]
        .value_counts()
        .to_dict(),
    )

    print(
        "Cross-system coding redundancy:",
        redundancy,
    )

    print(
        f"Output field count: {len(COLUMNS)}"
    )

    print(
        "FULL MERGE SCHEMA:",
        " | ".join(COLUMNS),
    )

    print(
        "FIRST ACCEPTED RECORD:"
    )

    for key, value in out.iloc[0].items():
        print(
            f"  {key:41} "
            f"{str(value)[:130]}"
        )

    # ========================================================
    # 6. WRITE VALIDATED DATA TO SAME DIRECTORY
    # ========================================================

    section(
        "5. OVERWRITE OF ORIGINAL "
        "PROCESSED DIRECTORY"
    )

    outputs = {
        "medi_drug_disease.csv": out,

        "medi_drug_disease_deduplicated.csv":
            dedup,

        "medi_low_confidence_audit.csv":
            low,

        "medi_rejected_rows.csv":
            rejects,
    }

    temps = {}

    try:
        # Stage and validate every CSV before
        # replacing the existing processed files.
        for filename, frame in outputs.items():
            temps[filename] = write_temp_csv(
                frame,
                destination / filename,
            )

        summary_path = (
            destination / "medi_summary.json"
        )

        summary = {
            "status": "success",
            "timestamp_utc": utc(),
            "source": "MEDI-2",

            "interpretation": (
                "MEDI ensemble medication-indication "
                "associations; no automatic "
                "regulatory approval inference"
            ),

            "licence": (
                "CC BY-NC-SA 3.0; review original "
                "terms before redistributing "
                "derived data"
            ),

            "source_file_inventory":
                inventory,

            "input_headers": {
                "MEDI-2.csv": list(icd.columns),
                "MEDI-2_UMLS.csv":
                    list(umls.columns),
            },

            "original_file_sha256": {
                name: checksum(path)
                for name, path in paths.items()
            },

            "input_rows": {
                "ICD": len(icd),
                "UMLS": len(umls),
            },

            "source_processing": {
                "ICD": xstats,
                "UMLS": ystats,
            },

            "join": (
                "none; separate source processing, "
                "then vertical concatenation"
            ),

            "source_row_identity": (
                "medi_source_file + "
                "medi_source_row_number"
            ),

            "rows_full": len(out),
            "rows_rejected": len(rejects),
            "rows_dedup": len(dedup),
            "rows_low_confidence": len(low),

            "precision_tiers":
                out["medi_confidence_tier"]
                .value_counts()
                .to_dict(),

            "coding_redundancy":
                redundancy,

            "output_columns":
                COLUMNS,

            "output_files": {
                name: str(destination / name)
                for name in outputs
            },

            "output_file_sha256": {
                name: checksum(tmp)
                for name, tmp in temps.items()
            },

            "licence_note": (
                "Data are subject to MEDI "
                "redistribution limitations; "
                "distributing only scripts "
                "may be preferable"
            ),
        }

        summary_tmp = write_temp_json(
            summary,
            summary_path,
        )

        temps["medi_summary.json"] = (
            summary_tmp
        )

        # Replace the current processed files.
        for filename, tmp in temps.items():
            tmp.replace(
                destination / filename
            )

            print(
                "[OVERWRITTEN]",
                destination / filename,
            )

        print(
            "SUCCESS. Use medi_drug_disease.csv "
            "(full source evidence) in the "
            "updated central merger."
        )

    finally:
        for tmp in temps.values():
            tmp.unlink(missing_ok=True)


if __name__ == "__main__":
    main()
