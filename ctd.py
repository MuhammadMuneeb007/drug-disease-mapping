
#!/usr/bin/env python3
"""CTD chemical-disease source-event processor (merge-compatible).

Default: reuse pinned data/raw/ctd/CTD_chemicals_diseases.tsv.gz and overwrite
         data/processed/ctd_aggregate_chemical_disease.csv only after success.
No automatic download, no full-data Pandas load, no source-event deduplication.

The ten common merge columns come first. CTD chemical associations are NOT
synonymous with clinically approved drug indications. Each source event is
retained with its line index, evidence class, identifiers, genes and citations.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import urlparse

import pandas as pd
import requests

ROOT = Path(__file__).resolve().parent

URLS = {
    "aggregate": "https://ctdbase.org/reports/CTD_chemicals_diseases.tsv.gz",
    "curated": "https://ctdbase.org/reports/CTD_curated_chemicals_diseases.tsv.gz",
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
    "chemical_name",
    "chemical_id",
    "cas_rn",
    "disease_name",
    "disease_id",
    "direct_evidence",
    "inference_gene_symbol",
    "inference_score",
    "omim_ids",
    "pubmed_ids",
    "ctd_evidence_class",
    "ctd_has_direct_evidence",
    "ctd_has_inference_evidence",
    "ctd_relationship_family",
    "ctd_source_row_index",
    "ctd_source_dataset",
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

SOURCE_FIELDS = [
    "ChemicalName",
    "ChemicalID",
    "CasRN",
    "DiseaseName",
    "DiseaseID",
    "DirectEvidence",
    "InferenceGeneSymbol",
    "InferenceScore",
    "OmimIDs",
    "PubMedIDs",
]

REQUIRED = {
    "ChemicalName",
    "ChemicalID",
    "DiseaseName",
    "DiseaseID",
}


def now():
    return datetime.now(timezone.utc).isoformat()


def heading(title):
    print(
        "\n" + "=" * 94
        + "\n" + title
        + "\n" + "=" * 94,
        flush=True,
    )


def sha256(path):
    h = hashlib.sha256()

    with open(path, "rb") as f:
        for chunk in iter(
            lambda: f.read(4 * 1024 * 1024),
            b"",
        ):
            h.update(chunk)

    return h.hexdigest()


def atomic_json(data, destination):
    destination.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    tmp = destination.with_name(
        destination.name + ".tmp"
    )

    try:
        tmp.write_text(
            json.dumps(
                data,
                indent=2,
                ensure_ascii=False,
                default=str,
            ) + "\n",
            encoding="utf-8",
        )
        tmp.replace(destination)
    finally:
        tmp.unlink(missing_ok=True)


def download_if_needed(
    path,
    url,
    allow=False,
    force=False,
):
    if (
        path.is_file()
        and path.stat().st_size
        and not force
    ):
        print(f"[EXISTING RAW DATA] {path}")
        return

    if not allow:
        raise FileNotFoundError(
            f"Missing raw CTD file: {path}; "
            "use --download_missing"
        )

    path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    tmp = path.with_name(
        path.name + ".download"
    )

    try:
        with requests.get(
            url,
            stream=True,
            timeout=1800,
        ) as r:
            r.raise_for_status()

            with tmp.open("wb") as output:
                for data in r.iter_content(
                    4 * 1024 * 1024
                ):
                    if data:
                        output.write(data)

        if not tmp.stat().st_size:
            raise RuntimeError(
                "Downloaded CTD file is empty"
            )

        with gzip.open(tmp, "rb") as f:
            if not f.read(1):
                raise RuntimeError(
                    "Downloaded CTD gzip contains no data"
                )

        tmp.replace(path)

    finally:
        tmp.unlink(missing_ok=True)


def scan_raw_directory(
    directory,
    input_path,
):
    heading("1. ORIGINAL DOWNLOADED FILES")

    files = []

    for item in sorted(
        directory.rglob("*")
    ):
        if not item.is_file():
            continue

        status = (
            "USED"
            if item.resolve() == input_path.resolve()
            else "OTHER"
        )

        files.append({
            "path": str(
                item.relative_to(directory)
            ),
            "bytes": item.stat().st_size,
            "processing_input": status == "USED",
        })

        print(
            f"  [{status:5}] "
            f"{item.relative_to(directory)} | "
            f"{item.stat().st_size:,} bytes"
        )

    print(
        f"Total local files: {len(files)}; "
        "selected source files: 1"
    )

    return files


def check_gzip(path):
    heading("2. GZIP INTEGRITY")

    byte_count = 0

    with gzip.open(path, "rb") as f:
        for b in iter(
            lambda: f.read(4 * 1024 * 1024),
            b"",
        ):
            byte_count += len(b)

    print(
        "GZIP CRC/decompression: PASS; "
        f"uncompressed bytes: {byte_count:,}"
    )

    return byte_count


def inspect_header(path):
    heading("3. ORIGINAL CTD TABULAR SCHEMA")

    comments = 0

    with gzip.open(
        path,
        "rt",
        encoding="utf-8-sig",
        errors="strict",
        newline="",
    ) as f:

        for index, line in enumerate(f):
            candidate = (
                line.rstrip("\r\n")
                .lstrip("#")
                .strip()
            )

            cols = [
                x.strip()
                for x in candidate.split("\t")
            ]

            if REQUIRED.issubset(cols):

                if len(cols) != len(set(cols)):
                    raise RuntimeError(
                        "Duplicate source column names"
                    )

                print(
                    f"Header is file line {index + 1}, "
                    f"following {comments} preamble lines"
                )

                print(
                    f"Columns ({len(cols)}): "
                    + " | ".join(cols)
                )

                return index, cols

            comments += 1

    raise RuntimeError(
        "Could not find CTD header "
        "in compressed source"
    )


def class_of_evidence(
    direct,
    gene,
    score,
):
    """Preserve direct and inferred evidence semantics."""

    tokens = {
        x.strip().casefold()
        for x in direct.split("|")
        if x.strip()
    }

    if tokens:
        known = {
            "therapeutic",
            "marker/mechanism",
        }

        if tokens - known:
            return "direct_unrecognised"

        if tokens == known:
            return (
                "direct_therapeutic_and_"
                "marker_mechanism"
            )

        if "therapeutic" in tokens:
            return "direct_therapeutic"

        return "direct_marker_mechanism"

    if gene:
        return "inferred_via_gene"

    if score:
        return "inferred_score_without_gene"

    return "unclassified_association"


def disease_identifier_metadata(ids):
    """
    Check identifier namespace and syntax.

    This does not establish ontology equivalence.
    """
    upper = (
        ids.str.upper()
        .str.replace(
            "_",
            ":",
            n=1,
            regex=False,
        )
    )

    prefix = (
        upper.str.extract(
            r"^([A-Z][A-Z0-9]*):",
            expand=False,
        )
        .fillna("")
    )

    namespace_by_prefix = {
        "MESH": "MeSH",
        "OMIM": "OMIM",
        "DOID": "DOID",
        "MONDO": "MONDO",
        "EFO": "EFO",
        "HP": "HPO",
        "ORPHA": "Orphanet",
        "ORPHANET": "Orphanet",
    }

    namespace = (
        prefix.map(namespace_by_prefix)
        .fillna("unresolved")
    )

    canonical = upper.copy()

    canonical = canonical.str.replace(
        r"^ORPHANET:",
        "ORPHA:",
        regex=True,
    )

    for p in ("MONDO", "EFO", "HP"):

        mask = canonical.str.fullmatch(
            p + r":\d{1,7}"
        )

        if mask.any():
            values = (
                canonical.loc[mask]
                .str.split(":", n=1)
                .str[-1]
            )

            canonical.loc[mask] = (
                p + ":" + values.str.zfill(7)
            )

    valid = (
        upper.str.fullmatch(
            r"MESH:[CD]\d{6,9}"
        )
        | upper.str.fullmatch(
            r"OMIM:\d{6}"
        )
        | upper.str.fullmatch(
            r"DOID:\d+"
        )
        | upper.str.fullmatch(
            r"(?:EFO|MONDO|HP):\d{1,7}"
        )
        | upper.str.fullmatch(
            r"(?:ORPHA|ORPHANET):\d+"
        )
    )

    canonical = canonical.where(
        valid, ""
    )

    namespace = namespace.where(
        valid, "unresolved"
    )

    validation = valid.map({
        True: "syntax_valid",
        False: "unresolved_or_invalid",
    })

    return namespace, canonical, validation


def transform(chunk, subset):
    """
    Preserve one output record per eligible
    original CTD source record.
    """
    missing = REQUIRED - set(chunk.columns)

    if missing:
        raise RuntimeError(
            "CTD source columns missing: "
            f"{sorted(missing)}"
        )

    for col in SOURCE_FIELDS:

        if col not in chunk:
            chunk[col] = ""

        chunk[col] = (
            chunk[col]
            .fillna("")
            .astype(str)
            .str.strip()
        )

    empty_core = {
        src: int(chunk[src].eq("").sum())
        for src in (
            "ChemicalName",
            "ChemicalID",
            "DiseaseName",
            "DiseaseID",
        )
    }

    mask = (
        chunk["ChemicalName"].ne("")
        & chunk["ChemicalID"].ne("")
        & chunk["DiseaseName"].ne("")
        & chunk["DiseaseID"].ne("")
    )

    out = chunk.loc[mask].copy()

    if out.empty:
        return (
            pd.DataFrame(columns=COLUMNS),
            empty_core,
            int((~mask).sum()),
        )

    out["ctd_source_row_index"] = (
        (out.index + 1).astype(str)
    )

    out["ctd_source_dataset"] = subset

    assignments = {
        "chemical_name": "ChemicalName",
        "chemical_id": "ChemicalID",
        "cas_rn": "CasRN",
        "disease_name": "DiseaseName",
        "disease_id": "DiseaseID",
        "direct_evidence": "DirectEvidence",
        "inference_gene_symbol":
            "InferenceGeneSymbol",
        "inference_score": "InferenceScore",
        "omim_ids": "OmimIDs",
        "pubmed_ids": "PubMedIDs",
    }

    for target, source in assignments.items():
        out[target] = out[source]

    out["drug_name"] = out["chemical_name"]

    out["drug_identifier"] = out["chemical_id"]

    out["disease_or_condition_name"] = (
        out["disease_name"]
    )

    out["disease_or_condition_identifier"] = (
        out["disease_id"]
    )

    direct_flag = out["direct_evidence"].ne("")

    inferred_flag = (
        out["inference_gene_symbol"].ne("")
        | out["inference_score"].ne("")
    )

    out["ctd_has_direct_evidence"] = direct_flag

    out["ctd_has_inference_evidence"] = (
        inferred_flag
    )

    # Classify each distinct combination once.
    keys = list(zip(
        out["direct_evidence"],
        out["inference_gene_symbol"].ne(""),
        out["inference_score"].ne(""),
    ))

    lookup = {
        key: class_of_evidence(*key)
        for key in set(keys)
    }

    out["ctd_evidence_class"] = [
        lookup[key]
        for key in keys
    ]

    # Important:
    # The current merger deduplicates by nine core
    # fields, including relationship_type/evidence_type.
    # These categories prevent direct therapeutic
    # evidence from being treated as equivalent to
    # inferred or marker/mechanism evidence.
    #
    # The merger will still require modification to
    # preserve multiple events within one class.

    out["relationship_type"] = (
        "ctd_" + out["ctd_evidence_class"]
    )

    out["evidence_type"] = (
        "ctd_" + out["ctd_evidence_class"]
    )

    out["ctd_relationship_family"] = (
        "chemical_disease_association"
    )

    out["source"] = "CTD"

    if subset == "aggregate":
        out["internal_source"] = (
            "CTD_aggregate_chemicals_diseases"
        )
    else:
        out["internal_source"] = (
            "CTD_curated_chemicals_diseases"
        )

    # Original CTD chemical identifiers are retained.
    # A separate canonical MeSH identifier is added
    # only when syntax can be recognised.

    ids = out["chemical_id"].str.upper()

    valid = ids.str.fullmatch(
        r"[CD]\d{6,9}"
    )

    out["drug_identifier_original"] = (
        out["chemical_id"]
    )

    out["drug_identifier_namespace"] = (
        valid.map({
            True: "MeSH",
            False: "unresolved",
        })
    )

    out["drug_identifier_canonical"] = (
        "MESH:" + ids
    ).where(valid, "")

    out["drug_identifier_validation"] = (
        valid.map({
            True: "syntax_valid",
            False: "unresolved_or_invalid",
        })
    )

    out["drug_identifier_type"] = (
        out["drug_identifier_namespace"]
    )

    out["disease_identifier_original"] = (
        out["disease_id"]
    )

    namespace, canonical, status = (
        disease_identifier_metadata(
            out["disease_id"]
        )
    )

    out["disease_identifier_namespace"] = namespace

    out["disease_identifier_canonical"] = canonical

    out["disease_identifier_validation"] = status

    out["disease_or_condition_identifier_type"] = (
        namespace
    )

    return (
        out.loc[:, COLUMNS],
        empty_core,
        int((~mask).sum()),
    )


def main():

    ap = argparse.ArgumentParser(
        description=__doc__
    )

    ap.add_argument(
        "--subset",
        choices=list(URLS),
        default="aggregate",
    )

    ap.add_argument(
        "--workdir",
        default=str(
            ROOT / "data/raw/ctd"
        ),
    )

    ap.add_argument(
        "--out",
        default=None,
    )

    ap.add_argument(
        "--metadata",
        default=None,
    )

    ap.add_argument(
        "--url",
        default=None,
    )

    ap.add_argument(
        "--chunksize",
        type=int,
        default=100_000,
    )

    ap.add_argument(
        "--head_n",
        type=int,
        default=3,
    )

    ap.add_argument(
        "--download_missing",
        action="store_true",
    )

    ap.add_argument(
        "--force_download",
        action="store_true",
    )

    ap.add_argument(
        "--skip_gzip_test",
        action="store_true",
    )

    args = ap.parse_args()

    if args.chunksize < 1:
        ap.error(
            "--chunksize must be positive"
        )

    subset = args.subset

    workdir = (
        Path(args.workdir)
        .expanduser()
        .resolve()
    )

    output = (
        Path(args.out)
        .expanduser()
        .resolve()
        if args.out
        else (
            ROOT
            / "data/processed"
            / f"ctd_{subset}_chemical_disease.csv"
        )
    )

    metadata = (
        Path(args.metadata)
        .expanduser()
        .resolve()
        if args.metadata
        else (
            ROOT
            / "data/metadata"
            / f"ctd_{subset}_processing_metadata.json"
        )
    )

    source_url = args.url or URLS[subset]

    rawname = Path(
        urlparse(source_url).path
    ).name

    if not rawname.endswith(".tsv.gz"):
        ap.error(
            "CTD source URL must end with .tsv.gz"
        )

    workdir.mkdir(
        parents=True,
        exist_ok=True,
    )

    infile = workdir / rawname

    temp = output.with_name(
        output.name + ".processing.tmp"
    )

    report = {
        "status": "started",
        "started_at_utc": now(),
        "source": "CTD",
        "subset": subset,
        "download_url": source_url,
        "source_semantics": (
            "chemical-disease associations, "
            "not proof of drug approval"
        ),
        "columns": COLUMNS,
    }

    try:

        download_if_needed(
            infile,
            source_url,
            args.download_missing
            or args.force_download,
            args.force_download,
        )

        report["raw_file_inventory"] = (
            scan_raw_directory(
                workdir,
                infile,
            )
        )

        report["raw_file"] = {
            "path": str(infile),
            "size_bytes": infile.stat().st_size,
            "sha256": sha256(infile),
        }

        report["gzip_uncompressed_bytes"] = (
            None
            if args.skip_gzip_test
            else check_gzip(infile)
        )

        header_line, original_cols = (
            inspect_header(infile)
        )

        report["header_file_line"] = (
            header_line + 1
        )

        report["original_columns"] = original_cols

        raw_head = pd.read_csv(
            infile,
            compression="gzip",
            sep="\t",
            skiprows=header_line + 1,
            header=None,
            names=original_cols,
            nrows=max(1, args.head_n),
            dtype=str,
            keep_default_na=False,
            on_bad_lines="error",
            encoding="utf-8-sig",
        )

        print(
            "Original source-record preview "
            "(first data row):"
        )

        if not raw_head.empty:
            for field, value in (
                raw_head.iloc[0].items()
            ):
                print(
                    f"  {field:35} "
                    f"{str(value)[:140]}"
                )

        report["raw_record_preview"] = (
            raw_head.head(
                args.head_n
            ).to_dict("records")
        )

        heading(
            "4. SOURCE-TO-OUTPUT MAPPING, "
            "EVIDENCE AND JOIN STRATEGY"
        )

        print(
            "ONE FILE; NO SQL JOIN. Each eligible "
            "input row becomes one evidence record."
        )

        print(
            "ChemicalName/ChemicalID -> "
            "drug_name/identifier "
            "(chemical != approved drug)"
        )

        print(
            "DiseaseName/DiseaseID -> "
            "disease_or_condition_name/identifier"
        )

        print(
            "DirectEvidence + inference fields -> "
            "ctd_evidence_class"
        )

        print(
            "Core relationship/evidence types "
            "preserve direct vs inferred distinctions."
        )

        print(
            "Ten common merge fields first; "
            "original evidence and identifiers follow."
        )

        print(
            "Source row numbers begin at 1 "
            "after the header, across chunks."
        )

        print(
            "The source-event output is NOT "
            "prematurely deduplicated."
        )

        counts = Counter()
        classes = Counter()
        direct_values = Counter()
        drug_validation = Counter()
        disease_validation = Counter()
        disease_ns = Counter()
        missing = Counter()

        sample_rows = []

        output.parent.mkdir(
            parents=True,
            exist_ok=True,
        )

        # Temporary output, committed only on success.
        with temp.open(
            "w",
            encoding="utf-8",
            newline="",
        ) as f:
            import csv
            csv.writer(f).writerow(
                COLUMNS
            )

        reader = pd.read_csv(
            infile,
            compression="gzip",
            sep="\t",
            skiprows=header_line + 1,
            header=None,
            names=original_cols,
            dtype=str,
            keep_default_na=False,
            chunksize=args.chunksize,
            on_bad_lines="error",
            encoding="utf-8-sig",
            low_memory=False,
        )

        for n, chunk in enumerate(
            reader,
            start=1,
        ):

            counts["rows_raw"] += len(chunk)

            out, core_missing, rejected = (
                transform(
                    chunk,
                    subset,
                )
            )

            missing.update(
                core_missing
            )

            counts[
                "rows_rejected_missing_core"
            ] += rejected

            counts["rows_output"] += len(out)

            if not out.empty:

                if len(sample_rows) < args.head_n:
                    sample_rows.extend(
                        out.head(
                            args.head_n
                            - len(sample_rows)
                        ).to_dict("records")
                    )

                classes.update(
                    out["ctd_evidence_class"]
                    .value_counts()
                    .to_dict()
                )

                direct_values.update(
                    out["direct_evidence"]
                    .replace("", "[empty]")
                    .value_counts()
                    .to_dict()
                )

                drug_validation.update(
                    out["drug_identifier_validation"]
                    .value_counts()
                    .to_dict()
                )

                disease_validation.update(
                    out["disease_identifier_validation"]
                    .value_counts()
                    .to_dict()
                )

                disease_ns.update(
                    out["disease_identifier_namespace"]
                    .value_counts()
                    .to_dict()
                )

                counts[
                    "rows_direct_and_inferred"
                ] += int((
                    out["ctd_has_direct_evidence"]
                    & out["ctd_has_inference_evidence"]
                ).sum())

                counts[
                    "rows_with_pubmed_ids"
                ] += int(
                    out["pubmed_ids"].ne("").sum()
                )

                out.to_csv(
                    temp,
                    mode="a",
                    header=False,
                    index=False,
                )

            print(
                f"[CHUNK {n:04d}] "
                f"raw={len(chunk):,} "
                f"accepted={len(out):,} "
                f"cumulative={counts['rows_output']:,}",
                flush=True,
            )

        if (
            counts["rows_raw"]
            != counts["rows_output"]
            + counts["rows_rejected_missing_core"]
        ):
            raise RuntimeError(
                "Source-to-output row counts "
                "do not reconcile"
            )

        if counts["rows_output"] == 0:
            raise RuntimeError(
                "No CTD rows passed validation; "
                "existing output remains untouched"
            )

        heading(
            "5. FINAL DATASET COUNTS "
            "AND EVIDENCE DISTRIBUTION"
        )

        for label, number in sorted(
            counts.items()
        ):
            print(
                f"  {label:38} {number:,}"
            )

        print(
            "  Nonexclusive missing-core counts:"
        )

        for label, number in sorted(
            missing.items()
        ):
            print(
                f"    {label:36} {number:,}"
            )

        print("  CTD evidence classes:")

        for label, number in (
            classes.most_common()
        ):
            print(
                f"    {label:48} {number:,}"
            )

        print("  Original DirectEvidence values:")

        for label, number in (
            direct_values.most_common()
        ):
            print(
                f"    {label:48} {number:,}"
            )

        print(
            "  Disease identifier namespace counts:"
        )

        for label, number in (
            disease_ns.most_common()
        ):
            print(
                f"    {label:36} {number:,}"
            )

        heading(
            "6. COMPLETE MERGE SCHEMA "
            "AND REAL OUTPUT SAMPLE"
        )

        print(
            f"Output fields ({len(COLUMNS)}): "
            + " | ".join(COLUMNS)
        )

        if sample_rows:
            print(
                "First accepted source record:"
            )

            for k, v in (
                sample_rows[0].items()
            ):
                print(
                    f"  {k:39} "
                    f"{str(v)[:130]}"
                )

        report.update({
            "input_rows":
                counts["rows_raw"],

            "output_rows":
                counts["rows_output"],

            "filtered_rows":
                counts["rows_rejected_missing_core"],

            "counts":
                dict(counts),

            "missing_input_fields_nonexclusive":
                dict(missing),

            "direct_evidence_counts":
                dict(direct_values),

            "ctd_evidence_class_counts":
                dict(classes),

            "drug_identifier_validation":
                dict(drug_validation),

            "disease_identifier_validation":
                dict(disease_validation),

            "disease_namespace_counts":
                dict(disease_ns),

            "sample_rows":
                sample_rows,

            "processing": {
                "join": "none",
                "combine": (
                    "one input row -> "
                    "one source evidence event"
                ),
                "filter": (
                    "four nonempty identity fields"
                ),
                "deduplication": (
                    "none at evidence-event stage"
                ),
                "merger_requirement": (
                    "preserve evidence events; "
                    "do not deduplicate on only "
                    "nine common fields"
                ),
            },
        })

        # Commit the complete output only after
        # all processing and row checks pass.
        temp.replace(output)

        report["output_file"] = str(output)

        report["output_file_sha256"] = (
            sha256(output)
        )

        report["finished_at_utc"] = now()

        report["status"] = "success"

        atomic_json(
            report,
            metadata,
        )

        heading(
            "SUCCESS: CTD OUTPUTS OVERWRITTEN "
            "IN ORIGINAL DIRECTORIES"
        )

        print("Main CSV:", output)
        print("Metadata:", metadata)
        print(
            "Rows:",
            f"{counts['rows_output']:,}",
        )

        print(
            "The main merger must still be "
            "updated before the final "
            "integrated rebuild."
        )

    except Exception as exc:

        report["finished_at_utc"] = now()
        report["status"] = "failed"
        report["error"] = repr(exc)

        failure = metadata.with_name(
            f"ctd_{subset}_failed_run.json"
        )

        atomic_json(
            report,
            failure,
        )

        print(
            "FAILED. Check the diagnostic JSON. "
            "Error:",
            repr(exc),
        )

        raise

    finally:
        temp.unlink(missing_ok=True)


if __name__ == "__main__":
    main()
