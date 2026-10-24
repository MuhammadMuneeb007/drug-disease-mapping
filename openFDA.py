
#!/usr/bin/env python3
"""Inspect (do not merge) pinned openFDA Human Drug Label partitions.

Examples:
  python inspect_openfda_fields.py --download-missing
  python inspect_openfda_fields.py --all-partitions --max-records 25 --download-missing
  python inspect_openfda_fields.py --all-partitions --full-scan --download-missing

A sampled inspection cannot demonstrate that a field is absent
from all labels.

FDA indication sections contain narrative text, not automatically
verified drug-disease ontology relationships.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json

from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import urlparse
from zipfile import ZipFile

import requests

try:
    import ijson
except ImportError as exc:
    raise SystemExit(
        "Missing dependency ijson. Install with: "
        "python -m pip install ijson requests"
    ) from exc


# ============================================================
# CONFIGURATION
# ============================================================

ROOT = Path(__file__).resolve().parent

RAW = ROOT / "data/raw/openfda_drug_label"

INSPECT = (
    ROOT / "data/processed/openfda/inspection"
)

MANIFEST_URL = "https://api.fda.gov/download.json"

SECTIONS = (
    "indications_and_usage",
    "purpose",
    "uses",
    "contraindications",
    "warnings",
    "warnings_and_precautions",
    "boxed_warning",
    "adverse_reactions",
    "clinical_studies",
    "description",
    "dosage_and_administration",
    "use_in_specific_populations",
)

IDENTIFIERS = (
    "brand_name",
    "generic_name",
    "substance_name",
    "rxcui",
    "product_ndc",
    "application_number",
    "spl_id",
    "spl_set_id",
    "product_type",
)


# ============================================================
# UTILITIES
# ============================================================

def now():
    return datetime.now(timezone.utc).isoformat()


def section(name):
    print(
        "\n" + "=" * 90
        + "\n" + name
        + "\n" + "=" * 90,
        flush=True,
    )


def sha256(path):
    h = hashlib.sha256()

    with path.open("rb") as f:
        for block in iter(
            lambda: f.read(4 * 1024 * 1024),
            b"",
        ):
            h.update(block)

    return h.hexdigest()


def write_json(obj, path):
    path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    temp = path.with_name(
        path.name + ".tmp"
    )

    try:
        temp.write_text(
            json.dumps(
                obj,
                indent=2,
                ensure_ascii=False,
                default=str,
            ) + "\n",
            encoding="utf-8",
        )

        temp.replace(path)

    finally:
        temp.unlink(missing_ok=True)


# ============================================================
# PINNED MANIFEST
# ============================================================

def pinned_manifest(path, refresh=False):

    if path.is_file() and not refresh:
        print(
            "[REUSE PINNED MANIFEST]",
            path,
        )

        return json.loads(
            path.read_text(
                encoding="utf-8"
            )
        )

    print(
        "[FETCH MANIFEST]",
        MANIFEST_URL,
    )

    r = requests.get(
        MANIFEST_URL,
        timeout=(20, 120),
    )

    r.raise_for_status()

    manifest = r.json()

    label = (
        manifest.get("results", {})
        .get("drug", {})
        .get("label", {})
    )

    if (
        not isinstance(
            label.get("partitions"),
            list,
        )
        or not label["partitions"]
    ):
        raise RuntimeError(
            "Manifest lacks "
            "results.drug.label.partitions"
        )

    # Avoid mixing older downloaded partitions
    # with a newly refreshed manifest.
    if refresh and path.is_file():

        old = json.loads(
            path.read_text(
                encoding="utf-8"
            )
        )

        old_label = (
            old.get("results", {})
            .get("drug", {})
            .get("label", {})
        )

        changed = (
            old_label.get("export_date")
            != label.get("export_date")
            or old_label.get("partitions")
            != label.get("partitions")
        )

        if changed:
            existing_zips = list(
                path.parent.glob("*.zip")
            )

            if existing_zips:
                raise RuntimeError(
                    "Manifest release changed, "
                    "but old ZIPs exist. "
                    "Do not mix snapshots. "
                    "Keep the pinned manifest "
                    "or archive/move existing ZIPs "
                    "before refreshing."
                )

    write_json(
        manifest,
        path,
    )

    return manifest


# ============================================================
# PARTITION DOWNLOAD
# ============================================================

def part_filename(part):

    url = str(
        part.get("file", "")
    )

    parsed = urlparse(url)

    if parsed.hostname != "download.open.fda.gov":
        raise RuntimeError(
            "Unexpected openFDA partition host: "
            + url
        )

    name = Path(
        parsed.path
    ).name

    if (
        not name.endswith(".json.zip")
        or name in {".", ".."}
    ):
        raise RuntimeError(
            "Unexpected partition filename: "
            + name
        )

    return name


def download_partition(
    part,
    raw,
    download_missing,
):

    name = part_filename(part)

    dest = raw / name

    if (
        dest.is_file()
        and dest.stat().st_size
    ):
        print(
            "[REUSE PARTITION]",
            name,
        )

        return dest

    if not download_missing:
        raise FileNotFoundError(
            f"Missing {dest}; "
            "use --download-missing "
            "to obtain this pinned partition"
        )

    url = str(
        part["file"]
    )

    if url.startswith(
        "http://download.open.fda.gov/"
    ):
        url = (
            "https://"
            + url[len("http://"):]
        )

    tmp = dest.with_name(
        dest.name + ".download"
    )

    try:
        print(
            "[DOWNLOAD]",
            url,
        )

        with requests.get(
            url,
            stream=True,
            timeout=(30, 1200),
        ) as r:

            r.raise_for_status()

            with tmp.open("wb") as output:
                for chunk in r.iter_content(
                    chunk_size=4 * 1024 * 1024
                ):
                    if chunk:
                        output.write(chunk)

        with ZipFile(tmp) as z:
            if not any(
                x.endswith(".json")
                for x in z.namelist()
            ):
                raise RuntimeError(
                    "Downloaded ZIP contains "
                    "no JSON member: " + name
                )

        tmp.replace(dest)

    finally:
        tmp.unlink(missing_ok=True)

    return dest


# ============================================================
# RECORD INSPECTION
# ============================================================

def preview(value, limit=650):

    if isinstance(
        value,
        (list, dict),
    ):
        value = json.dumps(
            value,
            ensure_ascii=False,
        )

    text = " ".join(
        str(value).split()
    )

    return (
        text[:limit]
        + (
            "…"
            if len(text) > limit
            else ""
        )
    )


def looks_like_structured_disease_key(path):
    """
    Detect possible structured disease-code
    field names for manual review.

    This does NOT verify disease ontology mapping.
    """

    tokens = path.lower()

    related = any(
        word in tokens
        for word in (
            "disease",
            "diagnos",
            "condition",
            "indication",
        )
    )

    explicit = any(
        word in tokens
        for word in (
            "_id",
            "_code",
            "_cui",
            "ontology",
            "_mesh",
            "_icd",
            "_snomed",
        )
    )

    return related and explicit


# ============================================================
# STREAMING ZIP INSPECTION
# ============================================================

def scan_zip(
    path,
    max_records,
    full_scan,
    verify_zip,
):

    totals = Counter()

    fields = Counter()
    openfda_fields = Counter()

    field_types = defaultdict(Counter)

    section_examples = defaultdict(list)
    identifier_examples = defaultdict(list)

    candidate_disease_fields = Counter()

    first_meta = None
    sample_record = None

    if verify_zip:
        with ZipFile(path) as z:
            bad = z.testzip()

            if bad:
                raise RuntimeError(
                    f"ZIP CRC failed for "
                    f"{path.name}: {bad}"
                )

        print(
            "[ZIP CRC] PASS",
            path.name,
        )

    with ZipFile(path) as z:

        zip_members = [
            info.filename
            for info in z.infolist()
            if info.filename.endswith(".json")
        ]

        if not zip_members:
            raise RuntimeError(
                f"No JSON entries in {path}"
            )

        for member in zip_members:

            # Read only the small metadata object.
            with z.open(member) as stream:

                meta_iter = ijson.items(
                    stream,
                    "meta",
                )

                if first_meta is None:
                    first_meta = next(
                        meta_iter,
                        {},
                    )

            # Stream result objects without
            # loading the entire JSON document.
            with z.open(member) as stream:

                for record in ijson.items(
                    stream,
                    "results.item",
                ):

                    if not isinstance(
                        record,
                        dict,
                    ):
                        raise RuntimeError(
                            "Non-object label in "
                            f"{path.name}:{member}"
                        )

                    totals["records"] += 1

                    if sample_record is None:
                        sample_record = {
                            k: preview(v, 200)
                            for k, v in record.items()
                            if k != "openfda"
                        }

                    locator = {
                        "partition": path.name,
                        "record_index":
                            totals["records"],
                        "id": str(
                            record.get("id", "")
                        ),
                        "set_id": str(
                            record.get("set_id", "")
                        ),
                        "effective_time": str(
                            record.get(
                                "effective_time",
                                "",
                            )
                        ),
                    }

                    # --------------------------
                    # All original record fields
                    # --------------------------

                    for key, value in record.items():

                        fields[key] += 1

                        field_types[key][
                            type(value).__name__
                        ] += 1

                        if looks_like_structured_disease_key(
                            key
                        ):
                            candidate_disease_fields[
                                key
                            ] += 1

                    # --------------------------
                    # Drug label text sections
                    # --------------------------

                    for key in SECTIONS:

                        val = record.get(key)

                        if val not in (
                            None,
                            "",
                            [],
                        ):
                            totals[
                                "section_" + key
                            ] += 1

                            if (
                                len(
                                    section_examples[key]
                                ) < 3
                            ):
                                section_examples[
                                    key
                                ].append({
                                    **locator,
                                    "text_preview":
                                        preview(val),
                                })

                    # --------------------------
                    # openfda harmonised fields
                    # --------------------------

                    annotation = record.get(
                        "openfda",
                        {},
                    )

                    if not isinstance(
                        annotation,
                        dict,
                    ):
                        totals[
                            "malformed_openfda"
                        ] += 1

                        annotation = {}

                    for key in annotation:

                        openfda_fields[
                            key
                        ] += 1

                        path_key = (
                            "openfda." + key
                        )

                        if looks_like_structured_disease_key(
                            path_key
                        ):
                            candidate_disease_fields[
                                path_key
                            ] += 1

                    for key in IDENTIFIERS:

                        val = annotation.get(key)

                        if val not in (
                            None,
                            "",
                            [],
                        ):
                            totals[
                                "identifier_" + key
                            ] += 1

                            if (
                                len(
                                    identifier_examples[key]
                                ) < 3
                            ):
                                identifier_examples[
                                    key
                                ].append({
                                    **locator,
                                    "value_preview":
                                        preview(
                                            val,
                                            200,
                                        ),
                                })

                    if (
                        not full_scan
                        and totals["records"]
                        >= max_records
                    ):
                        break

            if (
                not full_scan
                and totals["records"]
                >= max_records
            ):
                break

    return {
        "zip_members": zip_members,

        "zip_crc_verified":
            bool(verify_zip),

        "meta": first_meta,

        "sample_count":
            int(totals["records"]),

        "scan_complete":
            full_scan,

        "counters":
            dict(totals),

        "record_field_counts":
            dict(fields),

        "openfda_field_counts":
            dict(openfda_fields),

        "record_field_types": {
            k: dict(v)
            for k, v in field_types.items()
        },

        "section_examples":
            dict(section_examples),

        "identifier_examples":
            dict(identifier_examples),

        "candidate_disease_field_names":
            dict(candidate_disease_fields),

        "first_record_preview":
            sample_record,
    }


# ============================================================
# MAIN
# ============================================================

def main():

    p = argparse.ArgumentParser(
        description=__doc__,
    )

    p.add_argument(
        "--raw-dir",
        type=Path,
        default=RAW,
    )

    p.add_argument(
        "--inspect-dir",
        type=Path,
        default=INSPECT,
    )

    group = p.add_mutually_exclusive_group()

    group.add_argument(
        "--all-partitions",
        action="store_true",
    )

    group.add_argument(
        "--partition",
        type=int,
        default=1,
        help=(
            "1-based partition number; "
            "default: first partition"
        ),
    )

    p.add_argument(
        "--max-records",
        type=int,
        default=25,
        help=(
            "Sample size per partition; "
            "ignored by --full-scan"
        ),
    )

    p.add_argument(
        "--full-scan",
        action="store_true",
        help=(
            "Stream ALL records in selected ZIPs"
        ),
    )

    p.add_argument(
        "--download-missing",
        action="store_true",
    )

    p.add_argument(
        "--refresh-manifest",
        action="store_true",
    )

    p.add_argument(
        "--verify-zip",
        action="store_true",
        help=(
            "Verify ZIP CRC by fully "
            "decompressing selected archives"
        ),
    )

    args = p.parse_args()

    if args.max_records < 1:
        p.error(
            "--max-records must be >= 1"
        )

    raw = args.raw_dir.resolve()
    inspect = args.inspect_dir.resolve()

    raw.mkdir(
        parents=True,
        exist_ok=True,
    )

    inspect.mkdir(
        parents=True,
        exist_ok=True,
    )

    # --------------------------------------------------------
    # 1. Manifest
    # --------------------------------------------------------

    section(
        "1. PINNED OPENFDA "
        "DRUG-LABEL DOWNLOAD MANIFEST"
    )

    manifest_path = (
        raw / "openfda_download_manifest.json"
    )

    manifest = pinned_manifest(
        manifest_path,
        args.refresh_manifest,
    )

    endpoint = (
        manifest["results"]["drug"]["label"]
    )

    parts = endpoint["partitions"]

    names = [
        part_filename(x)
        for x in parts
    ]

    if len(names) != len(set(names)):
        raise RuntimeError(
            "Two manifest partitions resolve "
            "to the same ZIP filename"
        )

    print(
        "Export date:",
        endpoint.get(
            "export_date",
            "not reported",
        ),
    )

    print(
        "Partitions:",
        len(parts),
        "| reported records:",
        endpoint.get(
            "total_records",
            "not reported",
        ),
    )

    manifest_rows = []

    for index, (
        part,
        name,
    ) in enumerate(
        zip(parts, names),
        1,
    ):

        item = {
            "partition_index": index,
            "filename": name,
            "records":
                part.get("records"),
            "size_mb":
                part.get("size_mb"),
            "url":
                part["file"],
        }

        manifest_rows.append(item)

        print(
            f"  {index:2d}. {name}"
            f" | records={part.get('records')}"
            f" | size_mb={part.get('size_mb')}"
        )

    manifest_csv = (
        inspect
        / "openfda_label_partition_manifest.csv"
    )

    with manifest_csv.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as stream:

        writer = csv.DictWriter(
            stream,
            fieldnames=list(
                manifest_rows[0]
            ),
        )

        writer.writeheader()
        writer.writerows(
            manifest_rows
        )

    # --------------------------------------------------------
    # 2. Partition selection
    # --------------------------------------------------------

    chosen = (
        list(range(len(parts)))
        if args.all_partitions
        else [args.partition - 1]
    )

    if (
        min(chosen) < 0
        or max(chosen) >= len(parts)
    ):
        p.error(
            "--partition must be between "
            f"1 and {len(parts)}"
        )

    # --------------------------------------------------------
    # 3. Local raw files
    # --------------------------------------------------------

    section(
        "2. LOCAL RAW FILE INVENTORY"
    )

    local = []

    for file in sorted(
        raw.rglob("*")
    ):

        if file.is_file():

            item = {
                "path": str(
                    file.relative_to(raw)
                ),
                "bytes":
                    file.stat().st_size,
            }

            local.append(item)

            print(
                f"  {item['path']}"
                f" | {item['bytes']:,} bytes"
            )

    # --------------------------------------------------------
    # 4. Inspect selected partitions
    # --------------------------------------------------------

    aggregate_fields = Counter()
    aggregate_openfda = Counter()

    aggregate_sections = Counter()
    aggregate_identifiers = Counter()

    aggregate_candidates = Counter()

    partition_reports = []

    section(
        "3. RECORD FIELDS, LABEL SECTIONS "
        "AND PRODUCT IDENTIFIERS"
    )

    for idx in chosen:

        part = parts[idx]

        zip_path = download_partition(
            part,
            raw,
            args.download_missing,
        )

        report = scan_zip(
            zip_path,
            args.max_records,
            args.full_scan,
            args.verify_zip,
        )

        manifest_count = int(
            part["records"]
        )

        if (
            args.full_scan
            and report["sample_count"]
            != manifest_count
        ):
            raise RuntimeError(
                f"Record-count mismatch in "
                f"{zip_path.name}: "
                f"observed {report['sample_count']}, "
                f"manifest {manifest_count}. "
                "Check for corrupted or "
                "mixed-snapshot files."
            )

        report["partition_index"] = (
            idx + 1
        )

        report["filename"] = (
            zip_path.name
        )

        report["manifest_records"] = (
            manifest_count
        )

        report["manifest_size_mb"] = (
            part.get("size_mb")
        )

        report["zip_size_bytes"] = (
            zip_path.stat().st_size
        )

        report["zip_sha256"] = (
            sha256(zip_path)
        )

        partition_reports.append(
            report
        )

        aggregate_fields.update(
            report["record_field_counts"]
        )

        aggregate_openfda.update(
            report["openfda_field_counts"]
        )

        aggregate_sections.update({
            k.removeprefix("section_"): v
            for k, v
            in report["counters"].items()
            if k.startswith("section_")
        })

        aggregate_identifiers.update({
            k.removeprefix("identifier_"): v
            for k, v
            in report["counters"].items()
            if k.startswith("identifier_")
        })

        aggregate_candidates.update(
            report[
                "candidate_disease_field_names"
            ]
        )

        print(
            f"[{idx + 1}/{len(parts)}] "
            f"{zip_path.name}: "
            f"inspected="
            f"{report['sample_count']:,}"
            f" of {manifest_count:,}; "
            f"full_scan={args.full_scan}"
        )

        print(
            "  Indications-and-usage present:",
            report["counters"].get(
                "section_indications_and_usage",
                0,
            ),
        )

        print(
            "  Record keys:",
            ", ".join(
                sorted(
                    report[
                        "record_field_counts"
                    ]
                )
            ),
        )

        print(
            "  openfda keys:",
            ", ".join(
                sorted(
                    report[
                        "openfda_field_counts"
                    ]
                )
            ),
        )

        examples = report[
            "section_examples"
        ].get(
            "indications_and_usage",
            [],
        )

        for item in examples[:1]:
            print(
                "  Example indication text:",
                item["text_preview"][:300],
            )

    # --------------------------------------------------------
    # 5. Dataset summary
    # --------------------------------------------------------

    observed = sum(
        x["sample_count"]
        for x in partition_reports
    )

    full_census = bool(
        args.full_scan
        and len(chosen) == len(parts)
    )

    print(
        "\nTotal records inspected:",
        f"{observed:,}",
    )

    print(
        "All partitions inspected:",
        len(chosen) == len(parts),
    )

    print(
        "Complete dataset field census:",
        full_census,
    )

    print(
        "\nLabel section frequencies:",
        dict(aggregate_sections),
    )

    print(
        "Drug/product identifier coverage:",
        dict(aggregate_identifiers),
    )

    print(
        "Candidate structured disease-ID "
        "field names:",
        dict(aggregate_candidates),
    )

    print(
        "\nIMPORTANT: No indication-to-ontology "
        "mappings are verified by this inspector."
    )

    print(
        "Do not treat contraindications, "
        "adverse reactions, pharmacologic "
        "classes or incidental disease mentions "
        "as therapeutic indications."
    )

    # --------------------------------------------------------
    # 6. Complete inspection report
    # --------------------------------------------------------

    report = {
        "status": "success",

        "timestamp_utc": now(),

        "report_type":
            "inspection_only_no_association_extraction",

        "manifest_url":
            MANIFEST_URL,

        "manifest_snapshot_file":
            str(manifest_path),

        "manifest_sha256":
            sha256(manifest_path),

        "export_date":
            endpoint.get("export_date"),

        "manifest_partitions":
            len(parts),

        "manifest_total_records":
            endpoint.get("total_records"),

        "manifest_partition_rows":
            manifest_rows,

        "physical_raw_file_inventory_before_optional_downloads":
            local,

        "partitions_inspected":
            len(partition_reports),

        "records_inspected":
            observed,

        "full_dataset_field_census":
            full_census,

        "sampling_rule": (
            "all records per selected partition"
            if args.full_scan
            else
            "first N records per selected partition"
        ),

        "record_field_counts":
            dict(aggregate_fields),

        "openfda_field_counts":
            dict(aggregate_openfda),

        "section_counts":
            dict(aggregate_sections),

        "product_identifier_counts":
            dict(aggregate_identifiers),

        "potential_structured_disease_field_names_not_verified":
            dict(aggregate_candidates),

        "has_verified_structured_indication_disease_mapping":
            False,

        "scientific_interpretation": (
            "Human drug label sections contain "
            "narrative indication text. "
            "Automatic drug-disease assertions "
            "require additional text extraction, "
            "negation/section handling, and "
            "validated ontology mapping. "
            "Label presence alone does not "
            "establish FDA-approved product status."
        ),

        "partition_inspections":
            partition_reports,
    }

    report_file = (
        inspect
        / "openfda_field_inspection_report.json"
    )

    write_json(
        report,
        report_file,
    )

    section(
        "4. INSPECTION COMPLETE "
        "(NO DRUG-DISEASE CSV CREATED)"
    )

    print(
        "[MANIFEST CSV]",
        manifest_csv,
    )

    print(
        "[FIELD REPORT]",
        report_file,
    )

    print(
        "[RAW ZIP FILES]",
        raw,
    )


if __name__ == "__main__":
    main()
