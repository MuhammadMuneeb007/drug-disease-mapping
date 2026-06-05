#!/usr/bin/env python3
"""
TTD (Therapeutic Target Database) drug-disease processing.

Source:     https://idrblab.org/ttd/
Citation:   Zhang YT et al. Therapeutic target database 2026: facilitating
            targeted therapies and precision medicine. Nucleic Acids Research.
            54(D1): D1692-D1701 (2026). PMID: 41243978.
Licence:    Default copyright with citation-required free-access norm
            (no explicit open licence published with the data).

Two source files are processed:
  P1-05-Drug_disease.txt   - drug-to-disease mapping with ICD-11 codes
  P1-03-TTD_crossmatching.txt - PubChem/CAS/ChEBI/ATC IDs for TTD drugs

File format (both files):
  Multi-line per-drug records. Each line has three tab-separated parts:
      <TTD_DRUG_ID>\\t<FIELD_TAG>\\t<VALUE>
  Each drug has multiple lines (one per field). Common fields:
      TTDDRUID, DRUGNAME, INDICATI, PUBCHCID, CASNUMBE, CHEBI_ID, SUPDRATC

INDICATI value contains: disease name + [ICD-11: code] + clinical status,
in varying orders across versions. Parsed defensively with regex.
"""

import os
import re
import json
import urllib.request
from collections import defaultdict
from datetime import datetime

import pandas as pd


# ============================================================
# CONFIG
# ============================================================

TTD_BASE_URL = "https://db.idrblab.net/ttd/sites/default/files/ttd_database"
DRUG_DISEASE_URL = f"{TTD_BASE_URL}/P1-05-Drug_disease.txt"
CROSSMATCH_URL   = f"{TTD_BASE_URL}/P1-03-TTD_crossmatching.txt"

RAW_DIR = "/data/ascher02/uqmmune1/ANNOVAR/drug_disease_data/drug-disease-mapping/data/raw/ttd"
OUTPUT_DIR = "/data/ascher02/uqmmune1/ANNOVAR/drug_disease_data/drug-disease-mapping/data/processed/ttd"

DRUG_DISEASE_FILE = os.path.join(RAW_DIR, "P1-05-Drug_disease.txt")
CROSSMATCH_FILE   = os.path.join(RAW_DIR, "P1-03-TTD_crossmatching.txt")

OUTPUT_FILE  = os.path.join(OUTPUT_DIR, "ttd_drug_disease.csv")
DEDUP_FILE   = os.path.join(OUTPUT_DIR, "ttd_drug_disease_deduplicated.csv")
SUMMARY_FILE = os.path.join(OUTPUT_DIR, "ttd_summary.json")


# Order matters: longer/more specific keywords first, otherwise "Phase 1"
# would match before "Phase 1/2".
CLINICAL_STATUS_KEYWORDS = [
    "Phase 1/2", "Phase 2/3", "Phase 3/4",
    "Phase 0", "Phase 1", "Phase 2", "Phase 3", "Phase 4",
    "Approved", "Preclinical", "Investigative",
    "Discontinued in Phase 1", "Discontinued in Phase 2",
    "Discontinued in Phase 3", "Discontinued in Preclinical",
    "Discontinued", "Terminated", "Withdrawn from market", "Withdrawn",
    "Patented", "Clinical trial",
]

ICD_RE = re.compile(r"\[ICD-?(9|10|11)\s*:?\s*([^\]]+)\]", re.IGNORECASE)


# ============================================================
# HELPERS
# ============================================================

def clean_text(x):
    if x is None:
        return ""
    return str(x).strip()


def download_if_missing(url, path):
    if os.path.exists(path):
        size_kb = os.path.getsize(path) / 1024
        print(f"[INFO] Already present: {path} ({size_kb:.1f} KB)")
        return
    os.makedirs(os.path.dirname(path), exist_ok=True)
    print(f"[INFO] Downloading: {url}")
    urllib.request.urlretrieve(url, path)
    size_kb = os.path.getsize(path) / 1024
    print(f"[INFO] Saved to: {path} ({size_kb:.1f} KB)")


def parse_ttd_multiline(file_path):
    """
    Parse TTD multi-line record format.
    Returns: dict mapping ID -> dict of field -> list of values.
    """
    records = defaultdict(lambda: defaultdict(list))
    with open(file_path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.rstrip("\n").rstrip("\r")
            if not line:
                continue
            # Skip metadata / decoration lines
            if line.startswith("---") or line.startswith("Abbreviations"):
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            rec_id, field = parts[0].strip(), parts[1].strip()
            value = "\t".join(parts[2:]).strip()
            if rec_id and field and value:
                records[rec_id][field].append(value)
    return {k: dict(v) for k, v in records.items()}


def parse_indicati(text):
    """
    Extract (clinical_status, disease_name, icd_code, icd_version) from
    an INDICATI value. Robust to varying field order across TTD versions.
    """
    icd_code = ""
    icd_version = ""

    icd_match = ICD_RE.search(text)
    if icd_match:
        icd_version = f"ICD-{icd_match.group(1)}"
        icd_code = icd_match.group(2).strip()
        text = ICD_RE.sub("", text).strip()

    clinical_status = ""
    for kw in CLINICAL_STATUS_KEYWORDS:
        m = re.search(rf"(?<!\w){re.escape(kw)}(?!\w)", text, re.IGNORECASE)
        if m:
            clinical_status = kw
            text = (text[:m.start()] + " " + text[m.end():]).strip()
            break

    disease_name = re.sub(r"\s+", " ", text).strip(" ,;.")
    return clinical_status, disease_name, icd_code, icd_version


def classify_status(status):
    """Map TTD clinical status to (relationship_type, evidence_type)."""
    s = status.lower().strip()
    if not s:
        return "indication", "ttd_indication"
    if "approved" in s:
        return "approved_indication", "ttd_approved"
    if "withdrawn" in s:
        return "withdrawn_indication", "ttd_withdrawn"
    if "discontin" in s or "terminat" in s:
        return "discontinued_indication", "ttd_discontinued"
    if "preclinic" in s:
        return "preclinical_indication", "ttd_preclinical"
    if "patented" in s:
        return "patented_indication", "ttd_patented"
    if "phase" in s or "clinical trial" in s or "investig" in s:
        tag = re.sub(r"\W+", "_", s).strip("_")
        return "clinical_indication", f"ttd_{tag}"
    return "indication", "ttd_indication"


def build_crossmatch_map(crossmatch_records):
    """Build TTD drug ID -> external IDs lookup."""
    cm = {}
    for drug_id, fields in crossmatch_records.items():
        pubchem_raw = fields.get("PUBCHCID", [""])[0]
        cas_raw     = fields.get("CASNUMBE", [""])[0]
        chebi       = fields.get("CHEBI_ID", [""])[0]
        atc         = fields.get("SUPDRATC", [""])[0]
        drugname    = fields.get("DRUGNAME", [""])[0]

        # CASNUMBE values look like "CAS 65807-02-5"; strip the prefix
        cas = re.sub(r"^CAS\s*", "", cas_raw).strip() if cas_raw else ""
        pubchem = pubchem_raw.strip()

        cm[drug_id] = {
            "drugname": drugname,
            "pubchem_cid": pubchem,
            "cas": cas,
            "chebi": chebi,
            "atc": atc,
        }
    return cm


# ============================================================
# MAIN
# ============================================================

def main():
    print("=" * 100)
    print("TTD PROCESSING")
    print("=" * 100)
    print("[CITATION] Zhang YT et al. Therapeutic target database 2026. NAR 54(D1):D1692-D1701.")
    print("=" * 100)

    os.makedirs(OUTPUT_DIR, exist_ok=True)

    download_if_missing(DRUG_DISEASE_URL, DRUG_DISEASE_FILE)
    download_if_missing(CROSSMATCH_URL,   CROSSMATCH_FILE)

    # --- Parse cross-matching file -----------------------------------------
    print(f"[INFO] Parsing crossmatch: {CROSSMATCH_FILE}")
    crossmatch_records = parse_ttd_multiline(CROSSMATCH_FILE)
    print(f"[INFO] Crossmatch drug records: {len(crossmatch_records)}")
    crossmap = build_crossmatch_map(crossmatch_records)

    cm_with_pubchem = sum(1 for v in crossmap.values() if v["pubchem_cid"])
    cm_with_cas     = sum(1 for v in crossmap.values() if v["cas"])
    print(f"[INFO]   Drugs with PubChem CID: {cm_with_pubchem}")
    print(f"[INFO]   Drugs with CAS:         {cm_with_cas}")

    # --- Parse drug-disease file -------------------------------------------
    print(f"[INFO] Parsing drug-disease: {DRUG_DISEASE_FILE}")
    dd_records = parse_ttd_multiline(DRUG_DISEASE_FILE)
    print(f"[INFO] Drug records in P1-05: {len(dd_records)}")

    rows = []
    skipped_no_indication = 0
    skipped_no_disease = 0

    for drug_id, fields in dd_records.items():
        drug_name = clean_text(fields.get("DRUGNAME", [""])[0])
        if not drug_name and drug_id in crossmap:
            drug_name = clean_text(crossmap[drug_id]["drugname"])
        if not drug_name:
            continue

        indications = fields.get("INDICATI", [])
        if not indications:
            skipped_no_indication += 1
            continue

        cm = crossmap.get(drug_id, {})
        pubchem_cid = cm.get("pubchem_cid", "")
        cas         = cm.get("cas", "")
        chebi       = cm.get("chebi", "")
        atc         = cm.get("atc", "")

        # Prefer PubChem CID as primary identifier; fall back to TTD ID
        if pubchem_cid:
            drug_identifier = pubchem_cid
            drug_identifier_type = "PubChem"
        else:
            drug_identifier = drug_id
            drug_identifier_type = "TTD"

        for indication in indications:
            clinical_status, disease_name, icd_code, icd_version = parse_indicati(indication)

            if not disease_name:
                skipped_no_disease += 1
                continue

            relationship_type, evidence_type = classify_status(clinical_status)

            evidence_parts = []
            if clinical_status:
                evidence_parts.append(f"clinical_status={clinical_status}")
            if cas:
                evidence_parts.append(f"CAS={cas}")
            if chebi:
                evidence_parts.append(f"ChEBI={chebi}")
            if atc:
                evidence_parts.append(f"ATC={atc}")
            evidence_parts.append(f"TTD_drug_id={drug_id}")

            rows.append({
                "drug_name": drug_name,
                "drug_identifier": drug_identifier,
                "drug_identifier_type": drug_identifier_type,
                "disease_or_condition_name": disease_name,
                "disease_or_condition_identifier": icd_code,
                "disease_or_condition_identifier_type": icd_version,
                "relationship_type": relationship_type,
                "evidence_type": evidence_type,
                "source": "TTD",
                "internal_source": "TTD (drug-disease mapping with ICD)",
                "ttd_drug_id": drug_id,
                "clinical_status": clinical_status,
                "pubchem_cid": pubchem_cid,
                "cas": cas,
                "chebi": chebi,
                "atc": atc,
                "evidence_text": "; ".join(evidence_parts),
            })

    out = pd.DataFrame(rows)
    print(f"[INFO] Rows skipped (no indication):    {skipped_no_indication}")
    print(f"[INFO] Rows skipped (no disease parsed): {skipped_no_disease}")
    print(f"[INFO] Standardised rows: {len(out)}")

    out.to_csv(OUTPUT_FILE, index=False)

    dedup_cols = [
        "drug_name",
        "drug_identifier",
        "disease_or_condition_name",
        "disease_or_condition_identifier",
        "relationship_type",
    ]
    dedup = out.drop_duplicates(subset=dedup_cols)
    dedup.to_csv(DEDUP_FILE, index=False)

    # --- Summary ------------------------------------------------------------
    icd_version_counts = (
        out["disease_or_condition_identifier_type"].value_counts().to_dict()
        if len(out) else {}
    )

    summary = {
        "timestamp": datetime.utcnow().isoformat(),
        "source": "Therapeutic Target Database (TTD)",
        "source_url": "https://idrblab.org/ttd/",
        "data_url_drug_disease": DRUG_DISEASE_URL,
        "data_url_crossmatch":   CROSSMATCH_URL,
        "citation": (
            "Zhang YT, Zhou Y, Xu HW, et al. Therapeutic target database 2026. "
            "Nucleic Acids Research. 54(D1):D1692-D1701 (2026). PMID: 41243978."
        ),
        "licence": (
            "Default copyright with citation-required free-access norm; "
            "no explicit open licence published with the data."
        ),
        "drug_disease_file": DRUG_DISEASE_FILE,
        "crossmatch_file":   CROSSMATCH_FILE,
        "ttd_drugs_in_crossmatch": int(len(crossmatch_records)),
        "ttd_drugs_in_drug_disease": int(len(dd_records)),
        "rows_skipped_no_indication": int(skipped_no_indication),
        "rows_skipped_no_disease_parsed": int(skipped_no_disease),
        "rows_written": int(len(out)),
        "deduplicated_rows_written": int(len(dedup)),
        "unique_drugs_ttd_id": int(out["ttd_drug_id"].nunique()) if len(out) else 0,
        "unique_drugs_pubchem": (
            int(out[out["pubchem_cid"] != ""]["pubchem_cid"].nunique())
            if len(out) else 0
        ),
        "unique_diseases_text": int(out["disease_or_condition_name"].nunique()) if len(out) else 0,
        "unique_icd_codes": (
            int(out[out["disease_or_condition_identifier"] != ""]
                ["disease_or_condition_identifier"].nunique())
            if len(out) else 0
        ),
        "icd_version_counts": icd_version_counts,
        "clinical_status_counts": (
            out["clinical_status"].value_counts().to_dict() if len(out) else {}
        ),
        "relationship_type_counts": (
            out["relationship_type"].value_counts().to_dict() if len(out) else {}
        ),
        "output_file": OUTPUT_FILE,
        "deduplicated_file": DEDUP_FILE,
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