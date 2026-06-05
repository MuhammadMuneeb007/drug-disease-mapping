#!/usr/bin/env python3
"""
Hetionet v1.0 -> unified compound/drug-disease association CSV

Purpose
-------
Process Hetionet v1.0 Compound-Disease edges into a publication-ready,
provenance-aware CSV for a multi-source drug-disease association dataset.

Source repository:
    https://github.com/hetio/hetionet

Files used:
    hetionet-v1.0-nodes.tsv
    hetionet-v1.0-edges.sif.gz

Input URLs:
    https://github.com/hetio/hetionet/raw/main/hetnet/tsv/hetionet-v1.0-nodes.tsv
    https://github.com/hetio/hetionet/raw/main/hetnet/tsv/hetionet-v1.0-edges.sif.gz

Citation:
    Himmelstein DS et al. Systematic integration of biomedical knowledge
    prioritizes drugs for repurposing. eLife. 2017;6:e26726.
    doi:10.7554/eLife.26726

License:
    CC0 1.0 Universal / Public Domain

Important interpretation
------------------------
Hetionet records are knowledge-graph Compound-Disease relationships.
They should not be treated as the same evidence type as approved drug
indications from curated clinical drug resources.

Drug-disease metaedges retained:
    CtD = Compound-treats-Disease
    CpD = Compound-palliates-Disease

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

Hetionet-specific metadata columns
----------------------------------
metaedge
source_node_id
target_node_id
evidence_text

Requirements
------------
pip install pandas
"""

from __future__ import annotations

import gzip
import hashlib
import json
import os
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd


# ============================================================
# CONFIG
# ============================================================

HETIONET_BASE_URL = "https://github.com/hetio/hetionet/raw/main/hetnet/tsv"

NODES_URL = f"{HETIONET_BASE_URL}/hetionet-v1.0-nodes.tsv"
EDGES_URL = f"{HETIONET_BASE_URL}/hetionet-v1.0-edges.sif.gz"

RAW_DIR = (
    "/data/ascher02/uqmmune1/ANNOVAR/drug_disease_data/"
    "drug-disease-mapping/data/raw/hetionet"
)

OUTPUT_DIR = (
    "/data/ascher02/uqmmune1/ANNOVAR/drug_disease_data/"
    "drug-disease-mapping/data/processed/hetionet"
)

NODES_FILE = os.path.join(RAW_DIR, "hetionet-v1.0-nodes.tsv")
EDGES_FILE = os.path.join(RAW_DIR, "hetionet-v1.0-edges.sif.gz")

OUTPUT_FILE = os.path.join(OUTPUT_DIR, "hetionet_drug_disease.csv")
DEDUP_FILE = os.path.join(OUTPUT_DIR, "hetionet_drug_disease_deduplicated.csv")
SUMMARY_FILE = os.path.join(OUTPUT_DIR, "hetionet_summary.json")


# IMPORTANT:
# Do not call CtD "approved_indication".
# Hetionet CtD means Compound-treats-Disease, not necessarily approved indication.
DRUG_DISEASE_METAEDGES = {
    "CtD": (
        "compound_treats_disease",
        "hetionet_compound_treats_disease",
    ),
    "CpD": (
        "compound_palliates_disease",
        "hetionet_compound_palliates_disease",
    ),
}


OUTPUT_COLUMNS = [
    # Unified schema
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

    # Hetionet-specific metadata
    "metaedge",
    "source_node_id",
    "target_node_id",
    "evidence_text",
]


# ============================================================
# HELPERS
# ============================================================

def now_utc() -> str:
    return datetime.now(timezone.utc).isoformat()


def clean_text(x) -> str:
    if pd.isna(x):
        return ""
    return str(x).strip()


def sha256_file(path: str | Path, chunk_size: int = 1024 * 1024) -> str:
    path = Path(path)
    h = hashlib.sha256()

    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(chunk_size), b""):
            h.update(chunk)

    return h.hexdigest()


def human_size(num_bytes: int) -> str:
    size = float(num_bytes)

    for unit in ["B", "KB", "MB", "GB", "TB"]:
        if size < 1024:
            return f"{size:.2f} {unit}"
        size /= 1024

    return f"{size:.2f} PB"


def is_probable_git_lfs_pointer(path: str | Path) -> bool:
    """
    Detect accidental Git LFS pointer downloads.

    Git LFS pointer files are small text files that often contain:
        version https://git-lfs.github.com/spec/v1
        oid sha256:...
        size ...
    """
    path = Path(path)

    if not path.exists() or path.stat().st_size == 0:
        return True

    # Most real Hetionet files are much larger than 1 KB.
    # If small, inspect the beginning.
    if path.stat().st_size < 1024:
        try:
            text = path.read_text(encoding="utf-8", errors="ignore")
            if "git-lfs.github.com/spec" in text or "oid sha256:" in text:
                return True
        except Exception:
            return True

    return False


def test_gzip_file(path: str | Path) -> None:
    """
    Ensure gzip file can be read.
    """
    path = Path(path)

    total = 0
    with gzip.open(path, "rb") as f:
        while True:
            chunk = f.read(1024 * 1024)
            if not chunk:
                break
            total += len(chunk)

    print(f"[GZIP CHECK] PASS: {path} | uncompressed bytes={total:,}")


def download_if_missing(url: str, path: str, force: bool = False) -> None:
    """
    Download file if missing.

    If a Git LFS pointer is downloaded instead of real content, raise an error.
    """
    path_obj = Path(path)

    if path_obj.exists() and path_obj.stat().st_size > 0 and not force:
        print(f"[CACHE] Already present: {path}")
        print(f"        Size:   {human_size(path_obj.stat().st_size)}")
        print(f"        SHA256: {sha256_file(path_obj)}")

        if is_probable_git_lfs_pointer(path_obj):
            raise RuntimeError(
                f"Existing file appears to be a Git LFS pointer or empty: {path}\n"
                "Delete it and re-download manually from:\n"
                "https://github.com/hetio/hetionet/tree/main/hetnet/tsv"
            )

        return

    os.makedirs(path_obj.parent, exist_ok=True)

    print(f"[DOWNLOAD] {url}")
    urllib.request.urlretrieve(url, path)

    print(f"[DOWNLOADED] {path}")
    print(f"             Size:   {human_size(path_obj.stat().st_size)}")
    print(f"             SHA256: {sha256_file(path_obj)}")

    if is_probable_git_lfs_pointer(path_obj):
        raise RuntimeError(
            f"Downloaded file appears to be a Git LFS pointer or empty: {path}\n"
            "Download manually from:\n"
            "https://github.com/hetio/hetionet/tree/main/hetnet/tsv"
        )


def extract_identifier(node_id: str) -> str:
    """
    Hetionet node IDs are formatted as 'Kind::identifier'.

    Examples:
        Compound::DB00945  -> DB00945
        Disease::DOID:1612 -> DOID:1612
    """
    node_id = clean_text(node_id)

    if "::" in node_id:
        return node_id.split("::", 1)[1]

    return node_id


def infer_disease_identifier_type(identifier: str) -> str:
    identifier = clean_text(identifier)

    if not identifier:
        return ""

    upper = identifier.upper()

    if upper.startswith("DOID:"):
        return "DOID"

    if upper.startswith("MESH:"):
        return "MeSH"

    if upper.startswith("OMIM:"):
        return "OMIM"

    if upper.startswith("MONDO:"):
        return "MONDO"

    return "Hetionet_DiseaseID"


# ============================================================
# MAIN
# ============================================================

def main() -> None:
    print("=" * 100)
    print("HETIONET v1.0 DRUG-DISEASE PROCESSING")
    print("=" * 100)
    print("[SOURCE] Hetionet v1.0")
    print(f"[NODES URL] {NODES_URL}")
    print(f"[EDGES URL] {EDGES_URL}")
    print("[DOI] 10.7554/eLife.26726")
    print("[LICENSE] CC0 1.0 Universal")
    print("=" * 100)

    os.makedirs(RAW_DIR, exist_ok=True)
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    download_if_missing(NODES_URL, NODES_FILE)
    download_if_missing(EDGES_URL, EDGES_FILE)

    # Validate gzip edge file before reading.
    test_gzip_file(EDGES_FILE)

    # ------------------------------------------------------------
    # Load nodes
    # ------------------------------------------------------------
    print(f"[INFO] Reading nodes: {NODES_FILE}")

    nodes = pd.read_csv(
        NODES_FILE,
        sep="\t",
        dtype=str,
        keep_default_na=False,
    )

    print(f"[INFO] Nodes loaded: {len(nodes):,}")
    print(f"[INFO] Node columns: {list(nodes.columns)}")

    expected_node_cols = {"id", "name", "kind"}
    missing_node_cols = expected_node_cols - set(nodes.columns)

    if missing_node_cols:
        raise ValueError(
            f"Missing required node columns: {sorted(missing_node_cols)}"
        )

    id_to_name = dict(zip(nodes["id"], nodes["name"]))
    id_to_kind = dict(zip(nodes["id"], nodes["kind"]))

    n_compound = sum(1 for k in id_to_kind.values() if k == "Compound")
    n_disease = sum(1 for k in id_to_kind.values() if k == "Disease")

    print(f"[INFO] Compound nodes: {n_compound:,}")
    print(f"[INFO] Disease nodes:  {n_disease:,}")

    # ------------------------------------------------------------
    # Load edges
    # ------------------------------------------------------------
    print(f"[INFO] Reading edges: {EDGES_FILE}")

    edges = pd.read_csv(
        EDGES_FILE,
        sep="\t",
        dtype=str,
        keep_default_na=False,
        compression="gzip",
    )

    print(f"[INFO] Total edges loaded: {len(edges):,}")
    print(f"[INFO] Edge columns: {list(edges.columns)}")

    expected_edge_cols = {"source", "metaedge", "target"}
    missing_edge_cols = expected_edge_cols - set(edges.columns)

    if missing_edge_cols:
        raise ValueError(
            f"Missing required edge columns: {sorted(missing_edge_cols)}"
        )

    # ------------------------------------------------------------
    # Filter drug-disease edges
    # ------------------------------------------------------------
    edges_dd = edges[
        edges["metaedge"].isin(DRUG_DISEASE_METAEDGES.keys())
    ].copy()

    print(f"[INFO] Drug-disease edges retained: {len(edges_dd):,}")

    metaedge_counts = edges_dd["metaedge"].value_counts().to_dict()
    print(f"[INFO] Metaedge counts: {metaedge_counts}")

    # ------------------------------------------------------------
    # Build unified records
    # ------------------------------------------------------------
    rows = []

    skipped_wrong_kind = 0
    skipped_missing_name = 0

    for _, edge in edges_dd.iterrows():
        source_node_id = clean_text(edge["source"])
        target_node_id = clean_text(edge["target"])
        metaedge = clean_text(edge["metaedge"])

        # CtD and CpD are Compound -> Disease in Hetionet.
        # Keep this defensive check to avoid accidental direction errors.
        if id_to_kind.get(source_node_id) != "Compound":
            skipped_wrong_kind += 1
            continue

        if id_to_kind.get(target_node_id) != "Disease":
            skipped_wrong_kind += 1
            continue

        drug_name = clean_text(id_to_name.get(source_node_id, ""))
        disease_name = clean_text(id_to_name.get(target_node_id, ""))

        if not drug_name or not disease_name:
            skipped_missing_name += 1
            continue

        drug_identifier = extract_identifier(source_node_id)
        disease_identifier = extract_identifier(target_node_id)

        relationship_type, evidence_type = DRUG_DISEASE_METAEDGES[metaedge]

        rows.append(
            {
                "drug_name": drug_name,
                "drug_identifier": drug_identifier,
                "drug_identifier_type": "DrugBank",
                "disease_or_condition_name": disease_name,
                "disease_or_condition_identifier": disease_identifier,
                "disease_or_condition_identifier_type": infer_disease_identifier_type(
                    disease_identifier
                ),
                "relationship_type": relationship_type,
                "evidence_type": evidence_type,
                "source": "Hetionet",
                "internal_source": "Hetionet v1.0",
                "metaedge": metaedge,
                "source_node_id": source_node_id,
                "target_node_id": target_node_id,
                "evidence_text": f"metaedge={metaedge}",
            }
        )

    out = pd.DataFrame(rows, columns=OUTPUT_COLUMNS)

    print(f"[INFO] Standardised rows: {len(out):,}")
    print(f"[INFO] Skipped wrong kind/direction: {skipped_wrong_kind:,}")
    print(f"[INFO] Skipped missing name: {skipped_missing_name:,}")

    if len(out) == 0:
        raise RuntimeError("No Hetionet drug-disease rows were produced.")

    print()
    print("[OUTPUT PREVIEW]")
    print(out.head(10).to_string(index=False))

    out.to_csv(OUTPUT_FILE, index=False)

    # ------------------------------------------------------------
    # Deduplicate
    # ------------------------------------------------------------
    dedup_cols = [
        "drug_name",
        "drug_identifier",
        "disease_or_condition_name",
        "disease_or_condition_identifier",
        "relationship_type",
    ]

    dedup = out.drop_duplicates(subset=dedup_cols)
    dedup.to_csv(DEDUP_FILE, index=False)

    # ------------------------------------------------------------
    # Summary
    # ------------------------------------------------------------
    relationship_type_counts = out["relationship_type"].value_counts().to_dict()
    evidence_type_counts = out["evidence_type"].value_counts().to_dict()

    summary = {
        "timestamp_utc": now_utc(),
        "source": "Hetionet v1.0",
        "source_repository": "https://github.com/hetio/hetionet",
        "nodes_url": NODES_URL,
        "edges_url": EDGES_URL,
        "citation_doi": "10.7554/eLife.26726",
        "license": "CC0 1.0 Universal",
        "interpretation": (
            "Hetionet records are knowledge-graph Compound-Disease edges. "
            "CtD is interpreted as compound_treats_disease and CpD as "
            "compound_palliates_disease. These records are not labelled as "
            "approved indications."
        ),
        "nodes_file": NODES_FILE,
        "edges_file": EDGES_FILE,
        "nodes_file_sha256": sha256_file(NODES_FILE),
        "edges_file_sha256": sha256_file(EDGES_FILE),
        "nodes_loaded": int(len(nodes)),
        "compound_nodes": int(n_compound),
        "disease_nodes": int(n_disease),
        "total_edges_loaded": int(len(edges)),
        "drug_disease_edges_filtered": int(len(edges_dd)),
        "metaedge_counts": metaedge_counts,
        "rows_written": int(len(out)),
        "deduplicated_rows_written": int(len(dedup)),
        "skipped_wrong_kind_or_direction": int(skipped_wrong_kind),
        "skipped_missing_name": int(skipped_missing_name),
        "relationship_type_counts": relationship_type_counts,
        "evidence_type_counts": evidence_type_counts,
        "unique_drugs": int(out["drug_identifier"].nunique()),
        "unique_drug_names": int(out["drug_name"].nunique()),
        "unique_diseases": int(out["disease_or_condition_identifier"].nunique()),
        "unique_disease_names": int(out["disease_or_condition_name"].nunique()),
        "output_file": OUTPUT_FILE,
        "deduplicated_file": DEDUP_FILE,
        "output_file_sha256": sha256_file(OUTPUT_FILE),
        "deduplicated_file_sha256": sha256_file(DEDUP_FILE),
        "output_columns": OUTPUT_COLUMNS,
        "deduplication_columns": dedup_cols,
    }

    with open(SUMMARY_FILE, "w") as f:
        json.dump(summary, f, indent=2)

    print("=" * 100)
    print("DONE")
    print("=" * 100)
    print(f"[OUTPUT]  {OUTPUT_FILE}")
    print(f"[DEDUP]   {DEDUP_FILE}")
    print(f"[SUMMARY] {SUMMARY_FILE}")
    print("=" * 100)


if __name__ == "__main__":
    main()