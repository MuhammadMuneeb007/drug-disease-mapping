
#!/usr/bin/env python3
"""Open Targets Platform 26.03 clinical indications -> merge-compatible CSV.

Read existing Parquet files, check all partitions against the official release
listing (unless --offline), and overwrite original processed output only after
successful validation. No downloads unless --download-missing is passed.

Clinical indications aggregate clinical reports. Their evidence is not independent
of underlying AACT, ChEMBL, EMA, PMDA, etc. APPROVAL != PHASE_4.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import urljoin, urlparse

import pyarrow.parquet as pq
import requests

ROOT = Path(__file__).resolve().parent
BASE = "https://ftp.ebi.ac.uk/pub/databases/opentargets/platform"
RELEASE = "26.03"
DATASETS = ("clinical_indication", "drug_molecule", "disease")

CORE = [
    "drug_name", "drug_identifier", "drug_identifier_type",
    "disease_or_condition_name", "disease_or_condition_identifier",
    "disease_or_condition_identifier_type", "relationship_type", "evidence_type",
    "source", "internal_source",
]
META = [
    "drug_id", "disease_id", "disease_name",
    "max_clinical_stage_for_indication", "clinical_phase_raw",
    "clinical_stage_validation", "clinical_phase_numeric",
    "clinical_stage_is_explicit_approval", "clinical_report_ids",
    "clinical_report_ids_original", "clinical_report_count",
    "opentargets_indication_id", "drug_name_lookup_status",
    "disease_name_lookup_status", "dataset", "release",
    "opentargets_source_parquet", "opentargets_source_row_number",
    "evidence_granularity",
]
for entity in ("drug", "disease"):
    META += [f"{entity}_identifier_{suffix}" for suffix in
             ("original", "namespace", "canonical", "validation")]
COLUMNS = CORE + META
REJECT_COLUMNS = ["source_parquet", "source_row_number", "indication_id",
                  "drug_id", "disease_id", "rejection_reason"]

STAGES = {
    "UNKNOWN", "PRECLINICAL", "IND", "EARLY_PHASE_1", "PHASE_1",
    "PHASE_1_2", "PHASE_2", "PHASE_2_3", "PHASE_3", "PREAPPROVAL",
    "APPROVAL", "PHASE_4", "WITHDRAWAL",
}
PHASES = {f"PHASE_{n}": str(n) for n in range(1, 5)}


def now():
    return datetime.now(timezone.utc).isoformat()


def section(s):
    print("\n" + "=" * 90 + "\n" + s + "\n" + "=" * 90, flush=True)


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for b in iter(lambda: stream.read(4 * 1024 * 1024), b""):
            h.update(b)
    return h.hexdigest()


def value(x):
    return "" if x is None else str(x).strip()


def list_items(x):
    if x is None:
        return []
    if hasattr(x, "tolist"):
        x = x.tolist()
    if isinstance(x, (list, tuple)):
        return [value(v) for v in x if value(v)]
    if value(x):
        return [value(x)]
    return []


def normalize_stage(raw):
    """Do not conflate marketing authorisation with Phase IV trials."""
    original = value(raw)
    key = re.sub(r"[\s/-]+", "_", original.upper()).strip("_")
    variants = {"APPROVED": "APPROVAL", "LAUNCHED": "APPROVAL",
                "PHASE_IV": "PHASE_4", "PHASE_III": "PHASE_3",
                "PHASE_II": "PHASE_2", "PHASE_I": "PHASE_1",
                "EARLY_PHASE_I": "EARLY_PHASE_1", "PHASE_I_II": "PHASE_1_2",
                "PHASE_II_III": "PHASE_2_3", "PRE_APPROVAL": "PREAPPROVAL"}
    key = variants.get(key, key)
    if not original:
        return "", "missing", "", "False"
    if key not in STAGES:
        return original, "unrecognised", "", "False"
    return key, "recognised", PHASES.get(key, ""), str(key == "APPROVAL")


def identifier(x, entity):
    """Syntax-only validation; canonical syntax does not imply ontology mapping."""
    original = value(x)
    token = original.upper().replace("_", ":")
    if not original:
        return ("", "", "", "missing")
    if entity == "drug":
        if re.fullmatch(r"CHEMBL\d+", token):
            return (original, "ChEMBL", f"CHEMBL:{token}", "syntax_valid")
        return (original, "unresolved", "", "unrecognised_namespace")
    for pattern, prefix, label, padding in (
        (r"EFO:(\d+)", "EFO", "EFO", 7),
        (r"MONDO:(\d+)", "MONDO", "MONDO", 7),
        (r"HP:(\d+)", "HP", "HPO", 7),
        (r"DOID:(\d+)", "DOID", "DOID", 0),
        (r"(?:ORPHA|ORPHANET):(\d+)", "ORPHA", "Orphanet", 0),
        (r"OMIM:(\d+)", "OMIM", "OMIM", 0),
        (r"MESH:([CD]\d{6,9})", "MESH", "MeSH", 0),
    ):
        m = re.fullmatch(pattern, token)
        if m:
            suffix = m.group(1).zfill(padding) if padding else m.group(1)
            return (original, label, f"{prefix}:{suffix}", "syntax_valid")
    return (original, "unresolved", "", "unrecognised_namespace")


def source_listing(release, dataset):
    """Discover official Parquet partition names, not source data values."""
    base = f"{BASE}/{release}/output/{dataset}/"
    response = requests.get(
        base,
        timeout=(25, 120),
        headers={"User-Agent": "OT-Research-Exporter/3.0"},
    )
    response.raise_for_status()
    hrefs = re.findall(r'href=["\']([^"\']+\.parquet)["\']', response.text, re.I)
    result = {}
    for href in hrefs:
        full = urljoin(base, href)
        parsed = urlparse(full)
        if (
            parsed.netloc != "ftp.ebi.ac.uk"
            or parsed.path.rsplit("/", 1)[0]
            != urlparse(base).path.rstrip("/")
        ):
            continue
        filename = Path(parsed.path).name
        if filename in result and result[filename] != full:
            raise ValueError("Duplicate partition filename with different URLs")
        result[filename] = full
    if not result:
        raise RuntimeError(
            f"No parquet partitions found in official directory: {base}"
        )
    return result


def obtain(url, dest):
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_name(dest.name + ".download")
    try:
        with requests.get(url, stream=True, timeout=(30, 1800)) as r:
            r.raise_for_status()
            with tmp.open("wb") as f:
                for piece in r.iter_content(4 * 1024 * 1024):
                    if piece:
                        f.write(piece)
        if not tmp.is_file() or not tmp.stat().st_size:
            raise RuntimeError(f"Empty downloaded parquet: {url}")
        pq.ParquetFile(tmp)
        tmp.replace(dest)
    finally:
        tmp.unlink(missing_ok=True)


def choose_files(workdir, release, dataset, offline, download_missing):
    directory = workdir / release / dataset
    local = sorted(
        p for p in directory.glob("*.parquet") if p.is_file()
    )

    if offline:
        if not local:
            raise FileNotFoundError(
                f"No offline parquet for {dataset} in {directory}"
            )
        print(
            f"[OFFLINE/NOT INVENTORY VERIFIED] "
            f"{dataset}: {len(local)} local files"
        )
        return local, False, []

    listing = source_listing(release, dataset)
    expected = set(listing)
    existing = {p.name for p in local}
    extras = sorted(existing - expected)

    if extras:
        raise RuntimeError(
            f"Unexpected {dataset} parquet files: {extras[:15]}"
        )

    missing = sorted(expected - existing)

    if missing and not download_missing:
        raise FileNotFoundError(
            f"{dataset}: {len(missing)} missing partitions: {missing[:15]}. "
            "Pass --download-missing to fetch them, or --offline to inspect "
            "existing files (not lockable)."
        )

    for name in missing:
        print(f"[DOWNLOAD MISSING] {dataset}/{name}")
        obtain(listing[name], directory / name)

    found = sorted(directory / name for name in expected)

    for p in found:
        if p.stat().st_size == 0:
            raise RuntimeError(f"Empty parquet {p}")
        pq.ParquetFile(p)

    print(
        f"[VERIFIED INVENTORY] {dataset}: "
        f"{len(found)} official parquet files, missing=0"
    )
    return found, True, sorted(expected)


def inspect_source(dataset, paths):
    info = []
    required = (
        {"id", "name"}
        if dataset != "clinical_indication"
        else {
            "id", "drugId", "diseaseId",
            "maxClinicalStage", "clinicalReportIds",
        }
    )

    total = 0

    for path in paths:
        p = pq.ParquetFile(path)
        columns = p.schema_arrow.names
        missing = required - set(columns)

        if missing:
            raise RuntimeError(
                f"{path.name}: {dataset} missing source columns "
                f"{sorted(missing)}"
            )

        n = p.metadata.num_rows
        total += n

        info.append({
            "file": str(path),
            "rows": n,
            "columns": columns,
            "bytes": path.stat().st_size,
            "sha256": sha256(path),
        })

        print(
            f"  {path.name}: rows={n:,}; "
            f"fields={len(columns)}; bytes={path.stat().st_size:,}"
        )
        print("    HEADER: " + " | ".join(columns))

    print(
        f"[{dataset}] {len(paths)} partitions; "
        f"{total:,} records"
    )

    return info, total


def name_lookup(paths, kind):
    """Map unique names only; use identifier placeholder if conflicting."""
    names = defaultdict(set)

    for path in paths:
        reader = pq.ParquetFile(path)

        for batch in reader.iter_batches(
            batch_size=40000,
            columns=["id", "name"],
        ):
            for row in batch.to_pylist():
                id_ = value(row["id"])
                name = value(row["name"])
                if id_ and name:
                    names[id_].add(name)

    mapping = {
        k: next(iter(v))
        for k, v in names.items()
        if len(v) == 1
    }

    conflicts = {
        k: sorted(v)
        for k, v in names.items()
        if len(v) > 1
    }

    print(
        f"[{kind} LOOKUP] "
        f"unambiguous={len(mapping):,}, "
        f"conflicting IDs={len(conflicts):,}"
    )

    return mapping, conflicts


def optional_chembl_lookup(path):
    if not path:
        return {}, {}

    names = defaultdict(set)

    with Path(path).open(
        encoding="utf-8-sig",
        newline="",
    ) as f:
        reader = csv.DictReader(f)

        idcol = (
            "drug_id"
            if "drug_id" in reader.fieldnames
            else "drug_identifier"
        )

        if (
            idcol not in reader.fieldnames
            or "drug_name" not in reader.fieldnames
        ):
            raise RuntimeError(
                "ChEMBL fallback CSV needs "
                "drug_id/drug_identifier and drug_name"
            )

        for row in reader:
            drug_id = value(row[idcol])
            name = value(row["drug_name"])

            if drug_id and name and name != drug_id:
                names[drug_id].add(name)

    clean = {
        k: next(iter(v))
        for k, v in names.items()
        if len(v) == 1
    }

    conflicts = {
        k: sorted(v)
        for k, v in names.items()
        if len(v) > 1
    }

    print(
        f"[OPTIONAL ChEMBL MAP] "
        f"unambiguous={len(clean):,} "
        f"conflicts={len(conflicts):,}"
    )

    return clean, conflicts


def form_record(
    raw,
    drug_names,
    disease_names,
    drug_conflicts,
    disease_conflicts,
    chembl_names,
    chembl_conflicts,
    release,
    partition,
    row_index,
):
    drug = value(raw.get("drugId"))
    disease = value(raw.get("diseaseId"))
    indication = value(raw.get("id"))

    errors = []

    if not drug:
        errors.append("missing_drug_id")

    if not disease:
        errors.append("missing_disease_id")

    if errors:
        return None, {
            "source_parquet": partition,
            "source_row_number": row_index,
            "indication_id": indication,
            "drug_id": drug,
            "disease_id": disease,
            "rejection_reason": "|".join(errors),
        }

    dorig, dns, dcanonical, dvalid = identifier(drug, "drug")
    sorig, sns, scanonical, svalid = identifier(disease, "disease")

    if drug in drug_conflicts:
        drug_name, drug_lookup = (
            drug,
            "conflicting_names_identifier_placeholder",
        )
    elif drug in drug_names:
        drug_name, drug_lookup = drug_names[drug], "drug_molecule"
    elif drug in chembl_conflicts:
        drug_name, drug_lookup = (
            drug,
            "conflicting_chembl_names_identifier_placeholder",
        )
    elif drug in chembl_names:
        drug_name, drug_lookup = (
            chembl_names[drug],
            "optional_chembl_csv",
        )
    else:
        drug_name, drug_lookup = (
            drug,
            "unmatched_identifier_placeholder",
        )

    if disease in disease_conflicts:
        disease_name, disease_lookup = (
            disease,
            "conflicting_names_identifier_placeholder",
        )
    elif disease in disease_names:
        disease_name, disease_lookup = (
            disease_names[disease],
            "disease",
        )
    else:
        disease_name, disease_lookup = (
            disease,
            "unmatched_identifier_placeholder",
        )

    original_stage = value(raw.get("maxClinicalStage"))

    stage, stage_validation, phase_number, approval = normalize_stage(
        original_stage
    )

    refs = list_items(raw.get("clinicalReportIds"))
    normalized_refs = sorted(set(refs))

    out = {
        "drug_name": drug_name,
        "drug_identifier": drug,
        "drug_identifier_type": dns if dns else "unresolved",
        "disease_or_condition_name": disease_name,
        "disease_or_condition_identifier": disease,
        "disease_or_condition_identifier_type": sns if sns else "unresolved",
        "relationship_type": "clinical_indication",
        "evidence_type": "opentargets_clinical_indication",
        "source": "OpenTargets",
        "internal_source": (
            f"Open Targets Platform {release} clinical_indication"
        ),
        "drug_id": drug,
        "disease_id": disease,
        "disease_name": disease_name,
        "max_clinical_stage_for_indication": stage,
        "clinical_phase_raw": original_stage,
        "clinical_stage_validation": stage_validation,
        "clinical_phase_numeric": phase_number,
        "clinical_stage_is_explicit_approval": approval,
        "clinical_report_ids": ";".join(normalized_refs),
        "clinical_report_ids_original": json.dumps(
            refs,
            ensure_ascii=False,
        ),
        "clinical_report_count": len(normalized_refs),
        "opentargets_indication_id": indication,
        "drug_name_lookup_status": drug_lookup,
        "disease_name_lookup_status": disease_lookup,
        "dataset": "clinical_indication",
        "release": release,
        "opentargets_source_parquet": partition,
        "opentargets_source_row_number": row_index,
        "evidence_granularity": "aggregated_clinical_indication",
        "drug_identifier_original": dorig,
        "drug_identifier_namespace": dns,
        "drug_identifier_canonical": dcanonical,
        "drug_identifier_validation": dvalid,
        "disease_identifier_original": sorig,
        "disease_identifier_namespace": sns,
        "disease_identifier_canonical": scanonical,
        "disease_identifier_validation": svalid,
    }

    return out, None


def run(args):
    release = args.release

    if release != RELEASE:
        raise ValueError(
            f"This script pins release {RELEASE}; got {release!r}"
        )

    workdir = args.workdir.resolve()
    outfile = args.out.resolve()
    metapath = args.metadata.resolve()

    outfile.parent.mkdir(parents=True, exist_ok=True)
    metapath.parent.mkdir(parents=True, exist_ok=True)

    tmp = outfile.with_name(
        outfile.name + ".processing.tmp"
    )

    rejected_path = outfile.with_name(
        "opentargets_rejected_rows.csv"
    )

    reject_tmp = rejected_path.with_name(
        rejected_path.name + ".processing.tmp"
    )

    report = {
        "source": "OpenTargets",
        "release": release,
        "started_at_utc": now(),
        "status": "started",
        "complete_partition_inventory_verified": False,
        "interpretation": (
            "Aggregate drug-disease clinical indication; reports "
            "not independent of constituent sources"
        ),
        "core_columns": CORE,
        "output_columns": COLUMNS,
    }

    try:
        section("1. ORIGINAL RELEASE FILE INVENTORY")

        all_paths = {}
        inventories = {}
        input_counts = {}
        verified = True

        for dataset in DATASETS:
            paths, flag, names = choose_files(
                workdir,
                release,
                dataset,
                args.offline,
                args.download_missing,
            )

            all_paths[dataset] = paths

            info, count = inspect_source(dataset, paths)

            inventories[dataset] = {
                "partition_listing_verified": flag,
                "official_parquet_filenames": names,
                "files": info,
            }

            input_counts[dataset] = count
            verified &= flag

        report["datasets"] = inventories
        report["input_rows"] = input_counts
        report["complete_partition_inventory_verified"] = verified

        section("2. EXACT SOURCE LOOKUPS / NAMING")

        drug_names, drug_conflicts = name_lookup(
            all_paths["drug_molecule"], "drug"
        )

        disease_names, disease_conflicts = name_lookup(
            all_paths["disease"], "disease"
        )

        chembl_names, chembl_conflicts = optional_chembl_lookup(
            args.chembl
        )

        report["name_lookup"] = {
            "drug_unambiguous": len(drug_names),
            "disease_unambiguous": len(disease_names),
            "drug_conflicts": len(drug_conflicts),
            "disease_conflicts": len(disease_conflicts),
            "drug_conflict_examples": dict(
                list(drug_conflicts.items())[:15]
            ),
            "disease_conflict_examples": dict(
                list(disease_conflicts.items())[:15]
            ),
            "optional_chembl_unambiguous": len(chembl_names),
            "optional_chembl_conflicts": len(chembl_conflicts),
        }

        print(
            "JOIN: clinical_indication.drugId -> "
            "drug_molecule.id (unique lookup)"
        )

        print(
            "JOIN: clinical_indication.diseaseId -> "
            "disease.id (unique lookup)"
        )

        print(
            "AMBIGUOUS/MISSING NAMES: use identifier placeholder, "
            "never arbitrary last name"
        )

        print(
            "REPORT IDs: exact source list JSON plus "
            "sorted unique list; not independent evidence"
        )

        print(
            "STAGES: APPROVAL separate from PHASE_4 and WITHDRAWAL"
        )

        section(
            "3. STREAM INDICATION ROWS AND KEEP REPORT PROVENANCE"
        )

        counters = Counter()
        stages = Counter()
        name_status = Counter()
        disease_status = Counter()
        ns_counts = Counter()
        validation_counts = Counter()

        pairs = set()
        indication_ids = Counter()
        sample = []

        with (
            tmp.open(
                "w",
                encoding="utf-8",
                newline="",
            ) as full_f,
            reject_tmp.open(
                "w",
                encoding="utf-8",
                newline="",
            ) as rej_f,
        ):
            writer = csv.DictWriter(
                full_f,
                fieldnames=COLUMNS,
            )
            writer.writeheader()

            rejects = csv.DictWriter(
                rej_f,
                fieldnames=REJECT_COLUMNS,
            )
            rejects.writeheader()

            for path in all_paths["clinical_indication"]:
                parquet = pq.ParquetFile(path)
                source_row = 0

                columns = [
                    "id",
                    "drugId",
                    "diseaseId",
                    "maxClinicalStage",
                    "clinicalReportIds",
                ]

                for batch in parquet.iter_batches(
                    batch_size=args.batch_size,
                    columns=columns,
                ):
                    for raw in batch.to_pylist():
                        source_row += 1
                        counters["input_rows"] += 1

                        record, reject = form_record(
                            raw,
                            drug_names,
                            disease_names,
                            drug_conflicts,
                            disease_conflicts,
                            chembl_names,
                            chembl_conflicts,
                            release,
                            path.name,
                            source_row,
                        )

                        if reject:
                            rejects.writerow(reject)
                            counters["rejected_rows"] += 1
                            continue

                        writer.writerow(record)
                        counters["output_rows"] += 1

                        counters[
                            "clinical_report_links_unique_per_aggregate_sum"
                        ] += int(record["clinical_report_count"])

                        if record["clinical_report_count"] == 0:
                            counters[
                                "indications_without_report_ids"
                            ] += 1

                        stages[
                            record["max_clinical_stage_for_indication"]
                            or "[missing]"
                        ] += 1

                        name_status[
                            record["drug_name_lookup_status"]
                        ] += 1

                        disease_status[
                            record["disease_name_lookup_status"]
                        ] += 1

                        ns_counts[
                            record["disease_identifier_namespace"]
                        ] += 1

                        validation_counts[
                            record["disease_identifier_validation"]
                        ] += 1

                        pairs.add((
                            record["drug_identifier"],
                            record["disease_or_condition_identifier"],
                        ))

                        if record["opentargets_indication_id"]:
                            indication_ids[
                                record["opentargets_indication_id"]
                            ] += 1
                        else:
                            counters[
                                "missing_indication_ids"
                            ] += 1

                        if len(sample) < max(1, args.head_n):
                            sample.append(record)

                print(
                    f"  {path.name}: source rows={source_row:,}"
                )

                if source_row != parquet.metadata.num_rows:
                    raise RuntimeError(
                        f"Parquet row count mismatch in {path}"
                    )

        if (
            counters["input_rows"]
            != input_counts["clinical_indication"]
        ):
            raise RuntimeError(
                "Clinical-indication partition rows do not reconcile"
            )

        if (
            counters["input_rows"]
            != counters["output_rows"] + counters["rejected_rows"]
        ):
            raise RuntimeError(
                "Rejected + retained != input count"
            )

        if counters["output_rows"] == 0:
            raise RuntimeError(
                "No eligible Open Targets records; refusing overwrite"
            )

        duplicates = sum(
            n - 1
            for n in indication_ids.values()
            if n > 1
        )

        report.update({
            "row_counts": dict(counters),
            "distinct_drug_disease_pairs": len(pairs),
            "duplicate_nonempty_indication_ids": duplicates,
            "clinical_stage_counts": dict(stages),
            "drug_name_lookup_status_counts": dict(name_status),
            "disease_name_lookup_status_counts": dict(disease_status),
            "disease_namespace_counts": dict(ns_counts),
            "disease_validation_counts": dict(validation_counts),
            "sample_records": sample,
            "merge_caveat": (
                "The clinical_indication evidence is derived from "
                "clinical reports and overlaps with AACT/ChEMBL; "
                "retain clinical_report_ids and never assume independence."
            ),
        })

        section("4. OUTPUT VALIDATION AND COUNTS")

        for k, v in counters.items():
            print(f"  {k}: {v:,}")

        print(
            "  distinct drug-disease pairs:",
            len(pairs),
        )

        print(
            "  repeated indication IDs:",
            duplicates,
        )

        print("  stages:", dict(stages))
        print("  drug names:", dict(name_status))
        print("  disease names:", dict(disease_status))
        print("  disease namespace:", dict(ns_counts))

        print(
            "  output columns:",
            " | ".join(COLUMNS),
        )

        print(
            "  first record:",
            json.dumps(sample[0], indent=2)[:4500],
        )

        with tmp.open(
            encoding="utf-8",
            newline="",
        ) as f:
            if next(csv.reader(f)) != COLUMNS:
                raise RuntimeError("Output schema mismatch")

        # Commit only after all source reads, count reconciliation
        # and schema checks.
        tmp.replace(outfile)
        reject_tmp.replace(rejected_path)

        report["output_file"] = str(outfile)
        report["rejected_file"] = str(rejected_path)
        report["output_sha256"] = sha256(outfile)
        report["rejected_sha256"] = sha256(rejected_path)
        report["finished_at_utc"] = now()

        report["status"] = (
            "success"
            if verified
            else "offline_unverified_inventory"
        )

        report_tmp = metapath.with_name(
            metapath.name + ".processing.tmp"
        )

        report_tmp.write_text(
            json.dumps(
                report,
                indent=2,
                ensure_ascii=False,
                default=str,
            ) + "\n",
            encoding="utf-8",
        )

        report_tmp.replace(metapath)

        print("[OUTPUT OVERWRITTEN]", outfile)
        print("[AUDIT OVERWRITTEN]", rejected_path)
        print("[METADATA OVERWRITTEN]", metapath)

        if not verified:
            print(
                "[WARNING] --offline used: release partition "
                "completeness NOT independently verified. "
                "Do not lock yet."
            )

    except Exception as exc:
        failure = metapath.with_name(
            "opentargets_failed_run.json"
        )

        report.update({
            "status": "failed",
            "error": repr(exc),
            "finished_at_utc": now(),
        })

        failure.parent.mkdir(
            parents=True,
            exist_ok=True,
        )

        failure.write_text(
            json.dumps(
                report,
                indent=2,
                default=str,
            ),
            encoding="utf-8",
        )

        raise

    finally:
        tmp.unlink(missing_ok=True)
        reject_tmp.unlink(missing_ok=True)


def main():
    p = argparse.ArgumentParser(
        description=__doc__,
    )

    p.add_argument(
        "--release",
        default=RELEASE,
    )

    p.add_argument(
        "--workdir",
        type=Path,
        default=ROOT / "data/raw/opentargets",
    )

    p.add_argument(
        "--out",
        type=Path,
        default=ROOT / "data/processed/opentargets_drug_disease.csv",
    )

    p.add_argument(
        "--metadata",
        type=Path,
        default=ROOT / "data/metadata/opentargets_processing_metadata.json",
    )

    p.add_argument(
        "--chembl",
        type=Path,
        default=None,
        help="Optional unambiguous drug name fallback from processed ChEMBL",
    )

    p.add_argument(
        "--download-missing",
        action="store_true",
        help="Download only missing official partitions",
    )

    p.add_argument(
        "--offline",
        action="store_true",
        help=(
            "Use cached parquet without online inventory check; "
            "result NOT lockable until verified"
        ),
    )

    p.add_argument(
        "--batch-size",
        type=int,
        default=30000,
    )

    p.add_argument(
        "--head-n",
        type=int,
        default=3,
    )

    args = p.parse_args()

    if args.batch_size < 1:
        p.error("--batch-size must be >=1")

    if args.offline and args.download_missing:
        p.error(
            "--offline and --download-missing cannot both be used"
        )

    run(args)


if __name__ == "__main__":
    main()
