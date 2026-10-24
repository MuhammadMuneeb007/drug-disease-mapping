
#!/usr/bin/env python3
"""Hetionet v1.0 Compound-Disease edges -> source-event CSV.

Reuses existing data/raw/hetionet/ files. Run: python hetionet.py
Overwrites the SAME data/processed/hetionet output files only after checks.
CtD = treats; CpD = palliates; neither implies regulatory approval.
TSV lacks upstream license annotations: inspect Hetionet upstream terms before redistributing.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
import shutil
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

ROOT = Path(os.environ.get("DRUG_DISEASE_ROOT", Path(__file__).resolve().parent)).resolve()
RAW = ROOT / "data/raw/hetionet"
OUT = ROOT / "data/processed/hetionet"
BASE = "https://github.com/hetio/hetionet/raw/main/hetnet/tsv"
NODES_NAME = "hetionet-v1.0-nodes.tsv"
EDGES_NAME = "hetionet-v1.0-edges.sif.gz"

CORE = [
    "drug_name", "drug_identifier", "drug_identifier_type",
    "disease_or_condition_name", "disease_or_condition_identifier",
    "disease_or_condition_identifier_type", "relationship_type",
    "evidence_type", "source", "internal_source",
]
EXTRA = [
    "metaedge", "source_node_id", "target_node_id", "evidence_text",
    "hetionet_source_edge_row_number",
]
COLUMNS = CORE + EXTRA
RELATIONS = {
    "CtD": ("compound_treats_disease", "hetionet_compound_treats_disease"),
    "CpD": ("compound_palliates_disease", "hetionet_compound_palliates_disease"),
}


def heading(text):
    print("\n" + "=" * 94 + "\n" + text + "\n" + "=" * 94, flush=True)


def now():
    return datetime.now(timezone.utc).isoformat()


def checksum(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(4 * 1024 * 1024), b""):
            h.update(b)
    return h.hexdigest()


def assert_not_lfs_pointer(path):
    if not path.is_file() or path.stat().st_size == 0:
        raise ValueError(f"Empty or absent input: {path}")
    with open(path, "rb") as f:
        prefix = f.read(256)
    if b"git-lfs.github.com/spec" in prefix or b"oid sha256:" in prefix:
        raise ValueError(f"Git LFS pointer, not actual Hetionet data: {path}")


def ensure_file(path, download_missing):
    if path.is_file() and path.stat().st_size:
        assert_not_lfs_pointer(path)
        print(f"[REUSE] {path}")
        return
    if not download_missing:
        raise FileNotFoundError(f"Required file not present: {path}; use --download_missing")
    path.parent.mkdir(parents=True, exist_ok=True)
    url = f"{BASE}/{path.name}"
    tmp = path.with_name(path.name + ".download")
    try:
        print("[DOWNLOAD]", url)
        urllib.request.urlretrieve(url, tmp)
        assert_not_lfs_pointer(tmp)
        tmp.replace(path)
    finally:
        tmp.unlink(missing_ok=True)


def gzip_check(path):
    total = 0
    with gzip.open(path, "rb") as f:
        for chunk in iter(lambda: f.read(4 * 1024 * 1024), b""):
            total += len(chunk)
    print(f"[GZIP CRC CHECK] PASS; decompressed={total:,} bytes")
    return total


def clean_series(s):
    return s.fillna("").astype(str).str.strip()


def extract_id(s):
    return s.str.split("::", n=1).str[-1]


def disease_type(s):
    v = str(s).strip().upper()
    for prefix, namespace in (
        ("DOID:", "DOID"), ("MESH:", "MeSH"),
        ("OMIM:", "OMIM"), ("MONDO:", "MONDO"),
        ("EFO:", "EFO"), ("HP:", "HPO"),
        ("ORPHA:", "Orphanet"),
    ):
        if v.startswith(prefix):
            return namespace
    return "Hetionet_DiseaseID"


def write_csv_atomic(df, dest):
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_name(dest.name + ".tmp")
    try:
        df.to_csv(tmp, index=False)
        if list(pd.read_csv(tmp, nrows=0).columns) != list(df.columns):
            raise RuntimeError("Written CSV schema does not match expected columns")
        tmp.replace(dest)
    finally:
        tmp.unlink(missing_ok=True)


def write_json_atomic(obj, dest):
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_name(dest.name + ".tmp")
    try:
        tmp.write_text(json.dumps(obj, indent=2, ensure_ascii=False, default=str) + "\n", encoding="utf-8")
        tmp.replace(dest)
    finally:
        tmp.unlink(missing_ok=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--raw_dir", default=str(RAW))
    ap.add_argument("--output_dir", default=str(OUT))
    ap.add_argument("--download_missing", action="store_true")
    args = ap.parse_args()

    raw = Path(args.raw_dir).resolve()
    dest = Path(args.output_dir).resolve()
    nodes_path, edges_path = raw / NODES_NAME, raw / EDGES_NAME
    output = dest / "hetionet_drug_disease.csv"
    dedup_output = dest / "hetionet_drug_disease_deduplicated.csv"
    summary_output = dest / "hetionet_summary.json"

    heading("1. DOWNLOADED SOURCE FILE INVENTORY")
    for file in (nodes_path, edges_path):
        ensure_file(file, args.download_missing)
    inventory = []
    for file in sorted(raw.rglob("*")):
        if file.is_file():
            record = {
                "file": str(file.relative_to(raw)),
                "bytes": file.stat().st_size,
                "used": file.resolve() in {nodes_path.resolve(), edges_path.resolve()},
            }
            inventory.append(record)
            print("  ", "INPUT" if record["used"] else "OTHER", record["file"], f"{record['bytes']:,} bytes")
    print(f"Physical raw files: {len(inventory)}; processing inputs: 2")

    compressed_bytes = gzip_check(edges_path)
    heading("2. ACTUAL RAW HEADERS AND SOURCE RECORD COUNTS")
    nodes = pd.read_csv(nodes_path, sep="\t", dtype=str, keep_default_na=False, on_bad_lines="error")
    edges = pd.read_csv(edges_path, sep="\t", dtype=str, keep_default_na=False, compression="gzip", on_bad_lines="error")
    required_nodes = {"id", "name", "kind"}
    required_edges = {"source", "metaedge", "target"}
    if not required_nodes.issubset(nodes.columns) or not required_edges.issubset(edges.columns):
        raise RuntimeError(f"Wrong schemas: nodes={list(nodes.columns)}, edges={list(edges.columns)}")
    if nodes["id"].duplicated().any() or nodes["id"].eq("").any():
        raise RuntimeError("Node IDs are duplicated or empty: join would be ambiguous")
    if edges[list(required_edges)].isna().any().any():
        raise RuntimeError("Unexpected null graph edge fields")

    for name, frame, path in (("NODES", nodes, nodes_path), ("EDGES", edges, edges_path)):
        print(f"\n[{name}] {len(frame):,} rows; {len(frame.columns)} columns")
        print("Original header:", " | ".join(frame.columns))
        if len(frame):
            for col, value in frame.iloc[0].items():
                print(f"  {col:18} {str(value)[:160]}")
    print("\nOriginal metaedge counts: ")
    print(edges["metaedge"].value_counts().to_string())

    heading("3. GRAPH NODE LOOKUP AND ASSOCIATION MAPPING")
    print("No cross-database join. Within-source mapping: edges.source -> nodes.id for the compound; edges.target -> nodes.id for the disease.")
    print("Retain CtD (compound treats disease) and CpD (compound palliates disease); other metaedges excluded.")
    print("Raw node identifier Kind::ID -> remove Kind:: in unified identifier, keep full original source/target node IDs.")
    print("Keep source-edge row numbers, source evidence type, and original relationship labels.")
    print("These graph edges are not proof of regulatory approval or independent clinical validation.")

    selected = edges.loc[edges["metaedge"].isin(RELATIONS)].copy()
    selected["hetionet_source_edge_row_number"] = (selected.index + 1).astype(str)
    kind_map = nodes.set_index("id")["kind"]
    name_map = nodes.set_index("id")["name"]
    selected["source_kind"] = selected["source"].map(kind_map)
    selected["target_kind"] = selected["target"].map(kind_map)
    wrong = selected.loc[(selected["source_kind"] != "Compound") | (selected["target_kind"] != "Disease")]
    if len(wrong):
        raise RuntimeError(f"{len(wrong):,} CtD/CpD edges have missing/wrong-direction node types; first rows: {wrong.head(5).to_dict('records')}")

    selected["drug_name"] = clean_series(selected["source"].map(name_map))
    selected["disease_or_condition_name"] = clean_series(selected["target"].map(name_map))
    if (selected["drug_name"].eq("") | selected["disease_or_condition_name"].eq("")).any():
        raise RuntimeError("CtD/CpD source nodes have missing drug or disease names")

    out = pd.DataFrame(index=selected.index)
    out["drug_name"] = selected["drug_name"]
    out["drug_identifier"] = extract_id(selected["source"])
    out["drug_identifier_type"] = "DrugBank"
    out["disease_or_condition_name"] = selected["disease_or_condition_name"]
    out["disease_or_condition_identifier"] = extract_id(selected["target"])
    out["disease_or_condition_identifier_type"] = out["disease_or_condition_identifier"].map(disease_type)
    out["relationship_type"] = selected["metaedge"].map(lambda m: RELATIONS[m][0])
    out["evidence_type"] = selected["metaedge"].map(lambda m: RELATIONS[m][1])
    out["source"] = "Hetionet"
    out["internal_source"] = "Hetionet v1.0"
    out["metaedge"] = selected["metaedge"]
    out["source_node_id"] = selected["source"]
    out["target_node_id"] = selected["target"]
    out["evidence_text"] = "metaedge=" + selected["metaedge"]
    out["hetionet_source_edge_row_number"] = selected["hetionet_source_edge_row_number"]
    out = out.loc[:, COLUMNS].reset_index(drop=True)
    if out.empty:
        raise RuntimeError("No CtD/CpD edges found; refusing to overwrite outputs")
    if len(out) != len(selected):
        raise RuntimeError("Row conservation failed; evidence events unexpectedly lost")
    if out["drug_identifier"].eq("").any() or out["disease_or_condition_identifier"].eq("").any():
        raise RuntimeError("An output identifier is empty")
    if not out["drug_identifier"].str.fullmatch(r"DB\d{5}").all():
        raise RuntimeError("Hetionet Compound IDs not in expected DrugBank DB##### format; investigate source nodes")

    # Full file = all source edges. Dedup file = identical assertions only;
    # retain count and source-edge indices as extra metadata.
    assertion_key = ["drug_identifier", "disease_or_condition_identifier", "relationship_type", "evidence_type"]
    duplicate_rows = int(out.duplicated(subset=assertion_key).sum())
    dedup = out.drop_duplicates(subset=assertion_key).copy()
    support = out.groupby(assertion_key, sort=False)["hetionet_source_edge_row_number"].agg("|".join).rename("supporting_edge_rows")
    supports = out.groupby(assertion_key, sort=False).size().rename("supporting_edge_count")
    dedup = dedup.merge(support.reset_index(), on=assertion_key, how="left", validate="one_to_one")
    dedup = dedup.merge(supports.reset_index(), on=assertion_key, how="left", validate="one_to_one")
    if int(dedup["supporting_edge_count"].sum()) != len(out):
        raise RuntimeError("Deduplicated support links fail row accounting")

    heading("4. CONSERVATION, DEDUP AND OUTPUT SCHEMA")
    print(f"Graph nodes:               {len(nodes):,}")
    print(f"Graph edges:               {len(edges):,}")
    print(f"CtD/CpD source edges:      {len(selected):,}")
    print(f"Full evidence output:      {len(out):,}")
    print(f"Unique assertions:         {len(dedup):,}")
    print(f"Repeated assertions:       {duplicate_rows:,}")
    print("Relation distribution:", out["metaedge"].value_counts().to_dict())
    print("Output fields:", " | ".join(COLUMNS))
    print("First actual output row:")
    for col, value in out.iloc[0].items():
        print(f"  {col:38} {str(value)[:150]}")

    heading("5. OVERWRITING SAME PROCESSED DIRECTORY")
    write_csv_atomic(out, output)
    write_csv_atomic(dedup, dedup_output)
    summary = {
        "status": "success", "timestamp_utc": now(), "source": "Hetionet", "graph_version": "v1.0",
        "interpretation": "Knowledge-graph CtD treats / CpD palliates, not approved indications",
        "license_caveat": "Repository's original content CC0; heterogeneous upstream sources have their own licensing; TSV lacks per-edge license properties. Review source-specific licensing before redistribution.",
        "source_repository": "https://github.com/hetio/hetionet",
        "input_files_inventory": inventory,
        "input_headers": {"nodes": list(nodes.columns), "edges": list(edges.columns)},
        "nodes_file": str(nodes_path), "edges_file": str(edges_path),
        "nodes_sha256": checksum(nodes_path), "edges_sha256": checksum(edges_path),
        "edges_uncompressed_bytes": compressed_bytes,
        "nodes_count": len(nodes), "edges_count": len(edges),
        "compound_node_count": int(nodes["kind"].eq("Compound").sum()),
        "disease_node_count": int(nodes["kind"].eq("Disease").sum()),
        "ctd_cpd_edges": len(selected), "full_output_rows": len(out),
        "deduplicated_rows": len(dedup), "duplicate_assertion_rows": duplicate_rows,
        "metaedge_counts": out["metaedge"].value_counts().to_dict(),
        "join": "edges.source -> nodes.id, edges.target -> nodes.id; validated Compound -> Disease",
        "identity": "Raw node IDs retained; DrugBank and DOID identifiers extracted after Kind::",
        "core_schema": CORE, "output_schema": COLUMNS,
        "source_event_deduplication": "None in full output; assertion-level dedup retains source-edge support links",
        "output_file": str(output), "deduplicated_file": str(dedup_output),
        "output_sha256": checksum(output), "dedup_sha256": checksum(dedup_output),
    }
    write_json_atomic(summary, summary_output)
    for path in (output, dedup_output, summary_output):
        print("[WRITTEN]", path)
    print("SUCCESS. Feed hetionet_drug_disease.csv (FULL evidence) to the eventual merger.")


if __name__ == "__main__":
    main()
