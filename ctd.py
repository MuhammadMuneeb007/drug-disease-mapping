#!/usr/bin/env python3
"""
CTD -> unified chemical/drug–disease association CSV

Purpose
-------
Download and process CTD chemical–disease association data into a
publication-ready, provenance-aware CSV for a multi-source drug–disease
association dataset.

Recommended main input
----------------------
aggregate:
  CTD_chemicals_diseases.tsv.gz

Optional conservative input
---------------------------
curated:
  CTD_curated_chemicals_diseases.tsv.gz

Important interpretation
------------------------
CTD records represent chemical–disease associations, not necessarily
therapeutic drug indications. Some chemicals are drugs, but CTD should be
kept distinct from curated indication sources such as ChEMBL and from
clinical-trial co-occurrence sources such as AACT/ClinicalTrials.gov.

Unified output columns
----------------------
drug_name
drug_identifier
drug_identifier_type
disease_or_condition_name
disease_or_condition_identifier
disease_or_condition_identifier_type
relationship_type
evidence_type
source
internal_source

CTD-specific metadata columns
-----------------------------
chemical_name
chemical_id
cas_rn
disease_name
disease_id
direct_evidence
inference_gene_symbol
inference_score
omim_ids
pubmed_ids

Notes
-----
- CTD chemical names and identifiers are retained as CTD chemical-level fields.
- They are also copied into drug-compatible unified fields for harmonisation.
- This does NOT imply every CTD chemical is an approved drug.
- Deduplication is limited to exact duplicate rows after metadata retention.

Requirements
------------
pip install pandas requests
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import requests


# ============================================================
# CONFIG
# ============================================================

CTD_URLS = {
    "aggregate": "https://ctdbase.org/reports/CTD_chemicals_diseases.tsv.gz",
    "curated": "https://ctdbase.org/reports/CTD_curated_chemicals_diseases.tsv.gz",
}


CTD_OUTPUT_COLUMNS = [
    # Unified schema columns
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

    # CTD-specific retained metadata
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
]


# ============================================================
# HELPERS
# ============================================================

def now_utc() -> str:
    return datetime.now(timezone.utc).isoformat()


def print_section(title: str) -> None:
    print("\n" + "=" * 100)
    print(title)
    print("=" * 100)


def human_size(num_bytes: int) -> str:
    size = float(num_bytes)

    for unit in ["B", "KB", "MB", "GB", "TB"]:
        if size < 1024:
            return f"{size:.2f} {unit}"
        size /= 1024

    return f"{size:.2f} PB"


def sha256_file(path: Path, chunk_size: int = 1024 * 1024) -> str:
    h = hashlib.sha256()

    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(chunk_size), b""):
            h.update(chunk)

    return h.hexdigest()


def download_file(
    url: str,
    out_path: Path,
    force: bool = False,
    timeout: int = 1200,
) -> Path:
    out_path.parent.mkdir(parents=True, exist_ok=True)

    if out_path.exists() and out_path.stat().st_size > 0 and not force:
        print(f"[CACHE] Using existing file: {out_path}")
        print(f"        Size:   {human_size(out_path.stat().st_size)}")
        print(f"        SHA256: {sha256_file(out_path)}")
        return out_path

    print(f"[DOWNLOAD] {url}")

    tmp = out_path.with_suffix(out_path.suffix + ".part")

    with requests.get(url, stream=True, timeout=timeout) as r:
        r.raise_for_status()

        with open(tmp, "wb") as f:
            for chunk in r.iter_content(chunk_size=1024 * 1024):
                if chunk:
                    f.write(chunk)

    tmp.replace(out_path)

    print(f"[DOWNLOADED] {out_path}")
    print(f"             Size:   {human_size(out_path.stat().st_size)}")
    print(f"             SHA256: {sha256_file(out_path)}")

    return out_path


def test_gzip_integrity(path: Path) -> None:
    """
    Read through gzip stream to ensure the file is not corrupted.
    """
    print_section("GZIP integrity check")

    total_uncompressed = 0

    with gzip.open(path, "rb") as f:
        while True:
            chunk = f.read(1024 * 1024)

            if not chunk:
                break

            total_uncompressed += len(chunk)

    print("[GZIP CHECK] PASS")
    print(f"[UNCOMPRESSED BYTES READ] {total_uncompressed:,}")


def find_header_line_gz(path: Path) -> tuple[int, list[str], list[str]]:
    """
    CTD files begin with comment lines. The true header line may also begin
    with '#'. This function identifies the header by finding required columns.

    Returns:
      header_line_index, raw_comment_preview, header_columns
    """
    required = {"ChemicalName", "ChemicalID", "DiseaseName", "DiseaseID"}
    comment_preview: list[str] = []

    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as f:
        for i, line in enumerate(f):
            stripped = line.rstrip("\n")

            if i < 30:
                comment_preview.append(stripped)

            candidate = stripped.lstrip("#").strip()

            if not candidate:
                continue

            cols = [
                c.strip().lstrip("#").strip()
                for c in candidate.split("\t")
            ]

            if required.issubset(set(cols)):
                return i, comment_preview, cols

    raise RuntimeError(
        "Could not find CTD header line containing required columns."
    )


def normalize_columns(df: pd.DataFrame) -> pd.DataFrame:
    df.columns = [
        str(c).strip().lstrip("#").strip()
        for c in df.columns
    ]
    return df


def clean_series(s: pd.Series) -> pd.Series:
    return s.fillna("").astype(str).str.strip()


def get_column_or_empty(
    df: pd.DataFrame,
    candidates: list[str],
) -> pd.Series:
    """
    Return the first matching column as a clean string Series.
    Otherwise return empty strings.
    """
    for c in candidates:
        if c in df.columns:
            return clean_series(df[c])

    return pd.Series([""] * len(df), index=df.index, dtype="string")


def infer_ctd_disease_identifier_type(disease_id: str) -> str:
    """
    Infer disease identifier type from CTD DiseaseID string.

    Common examples:
      MESH:D003920
      OMIM:123456
    """
    x = str(disease_id).strip()

    if not x:
        return ""

    upper = x.upper()

    if upper.startswith("MESH:"):
        return "MeSH"

    if upper.startswith("OMIM:"):
        return "OMIM"

    if upper.startswith("DOID:"):
        return "DOID"

    if upper.startswith("MONDO:"):
        return "MONDO"

    return "CTD_DiseaseID"


def inspect_input_file(
    path: Path,
    header_line: int,
    header_columns: list[str],
    head_n: int,
) -> dict:
    print_section("Inspecting CTD input file")

    print(f"[INPUT] {path}")
    print(f"[SIZE] {human_size(path.stat().st_size)}")
    print(f"[SHA256] {sha256_file(path)}")
    print(f"[HEADER LINE] {header_line}")
    print(f"[HEADER COLUMNS] {header_columns}")

    head = pd.read_csv(
        path,
        sep="\t",
        compression="gzip",
        skiprows=header_line,
        header=0,
        dtype=str,
        nrows=head_n,
        low_memory=False,
    )

    head = normalize_columns(head)

    print("\n[INPUT HEAD]")
    print(head.to_string(index=False))

    return {
        "input_file": str(path),
        "input_size_bytes": path.stat().st_size,
        "input_sha256": sha256_file(path),
        "header_line": header_line,
        "header_columns": header_columns,
        "head_preview": head.to_dict(orient="records"),
    }


def transform_ctd_chunk(chunk: pd.DataFrame, subset: str) -> pd.DataFrame:
    """
    Transform one CTD input chunk into the unified output schema.

    Handles both aggregate and curated CTD files. Some fields exist only in
    the aggregate file, so they are filled with empty values when absent.
    """
    chunk = normalize_columns(chunk)

    required = [
        "ChemicalName",
        "ChemicalID",
        "DiseaseName",
        "DiseaseID",
    ]

    missing = [c for c in required if c not in chunk.columns]

    if missing:
        raise RuntimeError(
            f"Missing required CTD columns: {missing}\n"
            f"Available columns: {list(chunk.columns)}"
        )

    chemical_name = get_column_or_empty(chunk, ["ChemicalName"])
    chemical_id = get_column_or_empty(chunk, ["ChemicalID"])
    cas_rn = get_column_or_empty(chunk, ["CasRN", "CASRN", "CAS_RN"])

    disease_name = get_column_or_empty(chunk, ["DiseaseName"])
    disease_id = get_column_or_empty(chunk, ["DiseaseID"])

    direct_evidence = get_column_or_empty(chunk, ["DirectEvidence"])
    inference_gene_symbol = get_column_or_empty(
        chunk,
        ["InferenceGeneSymbol", "InferenceGeneSymbols"],
    )
    inference_score = get_column_or_empty(chunk, ["InferenceScore"])
    omim_ids = get_column_or_empty(chunk, ["OmimIDs", "OMIMIDs", "OMIMID"])
    pubmed_ids = get_column_or_empty(chunk, ["PubMedIDs", "PubMedID"])

    out = pd.DataFrame(
        {
            "chemical_name": chemical_name,
            "chemical_id": chemical_id,
            "cas_rn": cas_rn,
            "disease_name": disease_name,
            "disease_id": disease_id,
            "direct_evidence": direct_evidence,
            "inference_gene_symbol": inference_gene_symbol,
            "inference_score": inference_score,
            "omim_ids": omim_ids,
            "pubmed_ids": pubmed_ids,
        }
    )

    # Drop empty essential fields.
    out = out[
        (out["chemical_name"].astype(str).str.strip() != "")
        & (out["chemical_id"].astype(str).str.strip() != "")
        & (out["disease_name"].astype(str).str.strip() != "")
        & (out["disease_id"].astype(str).str.strip() != "")
    ].copy()

    if len(out) == 0:
        return pd.DataFrame(columns=CTD_OUTPUT_COLUMNS)

    # ------------------------------------------------------------
    # Unified schema aliases
    # ------------------------------------------------------------
    # Important: CTD is chemical-level. These unified columns are for
    # harmonisation only and do not imply every chemical is an approved drug.
    out["drug_name"] = out["chemical_name"]
    out["drug_identifier"] = out["chemical_id"]
    out["drug_identifier_type"] = "CTD_ChemicalID"

    out["disease_or_condition_name"] = out["disease_name"]
    out["disease_or_condition_identifier"] = out["disease_id"]
    out["disease_or_condition_identifier_type"] = out["disease_id"].apply(
        infer_ctd_disease_identifier_type
    )

    out["relationship_type"] = "chemical_disease_association"

    out["source"] = "CTD"

    if subset == "aggregate":
        out["internal_source"] = "CTD_aggregate_chemicals_diseases"
        out["evidence_type"] = "chemical_disease_association"

    elif subset == "curated":
        out["internal_source"] = "CTD_curated_chemicals_diseases"
        out["evidence_type"] = "curated_chemical_disease_association"

    else:
        out["internal_source"] = f"CTD_{subset}"
        out["evidence_type"] = "chemical_disease_association"

    # Exact duplicate removal only after preserving all useful metadata.
    out = out.drop_duplicates()

    return out[CTD_OUTPUT_COLUMNS]


def process_ctd_file(
    gz_path: Path,
    out_csv: Path,
    subset: str,
    header_line: int,
    chunksize: int,
) -> dict:
    print_section("Processing CTD file")

    out_csv.parent.mkdir(parents=True, exist_ok=True)

    if out_csv.exists():
        out_csv.unlink()

    total_input_rows = 0
    total_output_rows = 0
    wrote_header = False
    chunk_index = 0

    missing_counts = {
        "chemical_name": 0,
        "chemical_id": 0,
        "disease_name": 0,
        "disease_id": 0,
        "direct_evidence": 0,
        "inference_gene_symbol": 0,
        "inference_score": 0,
        "omim_ids": 0,
        "pubmed_ids": 0,
    }

    for chunk in pd.read_csv(
        gz_path,
        sep="\t",
        compression="gzip",
        skiprows=header_line,
        header=0,
        dtype=str,
        chunksize=chunksize,
        low_memory=False,
    ):
        chunk_index += 1
        chunk = normalize_columns(chunk)
        total_input_rows += len(chunk)

        out = transform_ctd_chunk(chunk, subset=subset)

        for col in missing_counts:
            if col in out.columns:
                missing_counts[col] += int(
                    (
                        out[col]
                        .fillna("")
                        .astype(str)
                        .str.strip()
                        == ""
                    ).sum()
                )

        if len(out) > 0:
            out.to_csv(
                out_csv,
                mode="a",
                index=False,
                header=(not wrote_header),
            )
            wrote_header = True

        total_output_rows += len(out)

        print(
            f"[CHUNK {chunk_index}] input={len(chunk):,} "
            f"output={len(out):,} total_output={total_output_rows:,}"
        )

    if not out_csv.exists() or out_csv.stat().st_size == 0:
        raise RuntimeError("Output CSV was not created or is empty.")

    print_section("Output validation")

    head = pd.read_csv(out_csv, dtype=str, nrows=10)

    print("[OUTPUT HEAD]")
    print(head.to_string(index=False))

    # Read full output for summary. If memory becomes an issue, replace later
    # with streaming summaries.
    full = pd.read_csv(out_csv, dtype=str, low_memory=False).fillna("")

    exact_duplicate_rows = int(full.duplicated().sum())

    unique_chemicals = int(full["chemical_id"].nunique(dropna=True))
    unique_chemical_names = int(full["chemical_name"].nunique(dropna=True))
    unique_diseases = int(full["disease_id"].nunique(dropna=True))
    unique_disease_names = int(full["disease_name"].nunique(dropna=True))

    unique_pairs = int(
        full[
            [
                "drug_identifier",
                "disease_or_condition_identifier",
            ]
        ]
        .drop_duplicates()
        .shape[0]
    )

    evidence_counts = (
        full["direct_evidence"]
        .fillna("")
        .replace("", "Missing")
        .value_counts()
        .head(30)
        .reset_index()
    )
    evidence_counts.columns = ["direct_evidence", "n"]

    relationship_type_counts = (
        full["relationship_type"]
        .value_counts()
        .to_dict()
    )

    evidence_type_counts = (
        full["evidence_type"]
        .value_counts()
        .to_dict()
    )

    disease_identifier_type_counts = (
        full["disease_or_condition_identifier_type"]
        .value_counts()
        .to_dict()
    )

    print("\n[OUTPUT SUMMARY]")
    print(f"Input rows processed:          {total_input_rows:,}")
    print(f"Output rows:                   {len(full):,}")
    print(f"Exact duplicate rows:          {exact_duplicate_rows:,}")
    print(f"Unique chemical IDs:           {unique_chemicals:,}")
    print(f"Unique chemical names:         {unique_chemical_names:,}")
    print(f"Unique disease IDs:            {unique_diseases:,}")
    print(f"Unique disease names:          {unique_disease_names:,}")
    print(f"Unique chemical-disease pairs: {unique_pairs:,}")
    print(f"Output SHA256:                 {sha256_file(out_csv)}")

    print("\n[DIRECT EVIDENCE COUNTS TOP 30]")
    print(evidence_counts.to_string(index=False))

    return {
        "output_file": str(out_csv),
        "output_size_bytes": out_csv.stat().st_size,
        "output_sha256": sha256_file(out_csv),
        "input_rows_processed": int(total_input_rows),
        "output_rows": int(len(full)),
        "exact_duplicate_rows": exact_duplicate_rows,
        "unique_chemical_ids": unique_chemicals,
        "unique_chemical_names": unique_chemical_names,
        "unique_disease_ids": unique_diseases,
        "unique_disease_names": unique_disease_names,
        "unique_chemical_disease_pairs": unique_pairs,
        "missing_counts": missing_counts,
        "direct_evidence_counts_top30": evidence_counts.to_dict(
            orient="records"
        ),
        "relationship_type_counts": relationship_type_counts,
        "evidence_type_counts": evidence_type_counts,
        "disease_identifier_type_counts": disease_identifier_type_counts,
        "output_columns": list(full.columns),
        "output_head": head.to_dict(orient="records"),
    }


# ============================================================
# MAIN
# ============================================================

def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Download and process CTD chemical-disease associations into "
            "a unified drug-disease-compatible schema."
        )
    )

    parser.add_argument(
        "--subset",
        choices=["aggregate", "curated"],
        default="aggregate",
        help=(
            "CTD subset to process. Use aggregate for broad coverage; "
            "curated for conservative curated-only associations."
        ),
    )

    parser.add_argument(
        "--url",
        default=None,
        help=(
            "Optional custom CTD URL. If omitted, uses the official URL "
            "for the selected subset."
        ),
    )

    parser.add_argument(
        "--workdir",
        default="./data/raw/ctd",
        help="Directory for downloaded CTD files.",
    )

    parser.add_argument(
        "--out",
        default=None,
        help="Output CSV path. If omitted, subset-specific default is used.",
    )

    parser.add_argument(
        "--metadata",
        default=None,
        help="Metadata JSON path. If omitted, subset-specific default is used.",
    )

    parser.add_argument(
        "--chunksize",
        type=int,
        default=250_000,
        help="Rows per chunk for streaming processing.",
    )

    parser.add_argument(
        "--head_n",
        type=int,
        default=10,
        help="Number of rows to print in previews.",
    )

    parser.add_argument(
        "--force_download",
        action="store_true",
        help="Force re-download even if file already exists.",
    )

    parser.add_argument(
        "--skip_gzip_test",
        action="store_true",
        help=(
            "Skip full gzip integrity test. Useful for speed if the file "
            "was already checked."
        ),
    )

    args = parser.parse_args()

    subset = args.subset
    url = args.url or CTD_URLS[subset]

    workdir = Path(args.workdir).resolve()
    workdir.mkdir(parents=True, exist_ok=True)

    out_csv = (
        Path(args.out).resolve()
        if args.out
        else Path(
            f"./data/processed/ctd_{subset}_chemical_disease.csv"
        ).resolve()
    )

    metadata_path = (
        Path(args.metadata).resolve()
        if args.metadata
        else Path(
            f"./data/metadata/ctd_{subset}_processing_metadata.json"
        ).resolve()
    )

    gz_path = workdir / Path(url).name

    metadata = {
        "script": "ctd_to_drug_disease.py",
        "source": "CTD",
        "subset": subset,
        "source_url": url,
        "relationship_type": "chemical_disease_association",
        "started_at_utc": now_utc(),
        "interpretation": (
            "CTD records represent chemical-disease associations and are not "
            "necessarily therapeutic drug indications. CTD chemical identifiers "
            "are copied into drug-compatible unified fields only for harmonisation."
        ),
        "unified_output_columns": CTD_OUTPUT_COLUMNS,
    }

    try:
        print_section("CTD -> unified chemical/drug-disease association processing")
        print(f"[SUBSET] {subset}")
        print(f"[URL] {url}")
        print(f"[WORKDIR] {workdir}")
        print(f"[OUTPUT] {out_csv}")

        gz_path = download_file(
            url=url,
            out_path=gz_path,
            force=args.force_download,
            timeout=1800,
        )

        metadata["downloaded_file"] = {
            "path": str(gz_path),
            "size_bytes": gz_path.stat().st_size,
            "sha256": sha256_file(gz_path),
        }

        if not args.skip_gzip_test:
            test_gzip_integrity(gz_path)
            metadata["gzip_integrity_check"] = "pass"
        else:
            print("[SKIP] GZIP integrity check skipped by user.")
            metadata["gzip_integrity_check"] = "skipped"

        header_line, comment_preview, header_columns = find_header_line_gz(
            gz_path
        )

        metadata["comment_preview"] = comment_preview
        metadata["header_line"] = header_line
        metadata["header_columns"] = header_columns

        metadata["input_inspection"] = inspect_input_file(
            path=gz_path,
            header_line=header_line,
            header_columns=header_columns,
            head_n=args.head_n,
        )

        metadata["processing"] = process_ctd_file(
            gz_path=gz_path,
            out_csv=out_csv,
            subset=subset,
            header_line=header_line,
            chunksize=args.chunksize,
        )

        metadata["finished_at_utc"] = now_utc()
        metadata["status"] = "success"

    except Exception as e:
        metadata["finished_at_utc"] = now_utc()
        metadata["status"] = "failed"
        metadata["error"] = repr(e)
        raise

    finally:
        metadata_path.parent.mkdir(parents=True, exist_ok=True)

        with open(metadata_path, "w", encoding="utf-8") as f:
            json.dump(metadata, f, indent=2)

        print_section("METADATA WRITTEN")
        print(metadata_path)


if __name__ == "__main__":
    main()