#!/usr/bin/env python3
"""Human-readable, READ-ONLY inspector for the eleven drug-disease sources.

Runs after download (even without data/processed); never imports source processors.
Only reads directory listings, CSV/TSV first rows, parquet schemas/small previews,
ZIP central directories, JSON metadata, and optional SQLite table schemas.
Does not scan entire raw or processed datasets.

Examples:
  python Analysis_0.1_InspectFiles.py
  python Analysis_0.1_InspectFiles.py --source CTD
  python Analysis_0.1_InspectFiles.py --markdown Supplementary_S0_Inspection.md
  python Analysis_0.1_InspectFiles.py --source AACT --archive-members
  python Analysis_0.1_InspectFiles.py --source ChEMBL --sqlite /path/to/chembl_37.db
"""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import sqlite3
import textwrap
import zipfile
from pathlib import Path
from urllib.parse import quote

WIDTH = 104
COMMON = [
    "drug_name", "drug_identifier", "drug_identifier_type",
    "disease_or_condition_name", "disease_or_condition_identifier",
    "disease_or_condition_identifier_type", "relationship_type",
    "evidence_type", "source", "internal_source",
]

# These are descriptions of the current repository processors, NOT validation of
# biological equivalence. Do not infer additional cross-namespace matches.
SOURCES = [
    dict(name="AACT", folder="aact_20260114", parser="aact.py",
         inputs=["extract/interventions.txt", "extract/conditions.txt", "extract/studies.txt"],
         processed="aact/aact_drug_disease.csv", method="Relational joins",
         steps=[
             "Read the three pipe-separated tables with their first row as the header.",
             "INNER JOIN interventions with conditions USING nct_id; this pairs each eligible intervention with conditions listed in the same trial.",
             "LEFT JOIN studies USING nct_id to attach phase, status and study titles.",
             "Default filter retains DRUG interventions with nonempty intervention name, condition name and nct_id.",
             "Create trial drug-condition co-occurrence records; this is NOT evidence of efficacy or approval."],
         maps=[("interventions.name", "drug_name"), ("conditions.name", "disease_or_condition_name"),
               ("interventions.nct_id", "nct_id; join key"),
               ("interventions.intervention_type", "intervention_type; filter"),
               ("studies.phase", "study_phase"), ("studies.overall_status", "overall_status"),
               ("studies.brief_title", "brief_title"), ("studies.official_title", "official_title")],
         caveat="Only three of the 49 ZIP members are used by this parser. Raw samples come from individual tables, not one joined trial."),
    dict(name="Broad Drug Repurposing Hub", folder="broad_repurposing", parser="Broad.py",
         inputs=["repurposing_drugs_20200324.txt"],
         processed="broad_repurposing/broad_repurposing_drug_disease.csv", method="Single-file expansion",
         steps=["Read the tab-separated file after skipping initial metadata lines beginning with !.",
                "Use pert_iname as compound name. Split indication at | into separate drug-indication rows.",
                "Preserve clinical_phase, disease_area, moa, target and original indication text.",
                "Exclude source rows with missing drug name or indication. No table JOIN occurs."],
         maps=[("pert_iname", "drug_name"), ("indication", "disease_or_condition_name; split |"),
               ("clinical_phase", "clinical_phase"), ("disease_area", "disease_area"),
               ("moa", "moa"), ("target", "target")],
         caveat="Clinical phase remains a source annotation; Launched, Preclinical and Withdrawn are not equivalent."),
    dict(name="ChEMBL", folder="chembl", parser="chembl.py",
         inputs=["chembl_37_sqlite.tar.gz"], processed="chembl_drug_disease.csv",
         method="SQL extraction and conditional joins",
         steps=["Use the ChEMBL 37 SQLite database inside the downloaded archive. ChEMBL 36 is an additional local snapshot, not the processed release.",
                "Read drug_indication (alias di) as the association table.",
                "Conditionally LEFT JOIN molecule_dictionary for drug ID/name using the molecule key present in this release.",
                "Conditionally LEFT JOIN disease_dictionary on efo_id only if it exists; otherwise use the indication row efo_term.",
                "Extract indication-level phase, MeSH annotation, drug ID/name and disease ID/name."],
         maps=[("drug_indication.molregno / molecule ID", "drug_id (resolved ChEMBL ID)"),
               ("molecule_dictionary.pref_name", "drug_name, if available"),
               ("drug_indication.efo_id", "disease_id"),
               ("drug_indication.efo_term / dictionary", "disease_name"),
               ("drug_indication.mesh_id", "mesh_id"),
               ("drug_indication.mesh_heading", "mesh_heading"),
               ("drug_indication.max_phase_for_ind", "max_phase (indication-specific)")],
         caveat="The raw TAR.GZ is not a flat table. Its actual SQL tables/headers cannot be confirmed without an extracted SQLite file; use --sqlite for a read-only schema check."),
    dict(name="CTD", folder="ctd", parser="ctd.py",
         inputs=["CTD_chemicals_diseases.tsv.gz"], processed="ctd_aggregate_chemical_disease.csv",
         method="Single-file field mapping",
         steps=["Detect the real tabular header after the CTD comment preamble; read the compressed TSV in chunks.",
                "Require ChemicalName, ChemicalID, DiseaseName and DiseaseID; remove records lacking these required fields.",
                "Copy chemical and disease identity fields to the common columns and retain direct evidence, inferred gene, score, OMIM and PubMed fields.",
                "Remove exact duplicate full records within processing chunks; later merge-level deduplication uses a narrower key."],
         maps=[("ChemicalName", "chemical_name; drug_name"), ("ChemicalID", "chemical_id; drug_identifier"),
               ("CasRN", "cas_rn"), ("DiseaseName", "disease_name; disease_or_condition_name"),
               ("DiseaseID", "disease_id; disease_or_condition_identifier"),
               ("DirectEvidence", "direct_evidence"), ("InferenceGeneSymbol", "inference_gene_symbol"),
               ("InferenceScore", "inference_score"), ("OmimIDs", "omim_ids"), ("PubMedIDs", "pubmed_ids")],
         caveat="CTD contains chemical–disease associations, including inferred/nontherapeutic ones. The generic evidence_type does not encode DirectEvidence."),
    dict(name="DrugCentral", folder="drugcentral", parser="drugcentral.py",
         inputs=[], processed="drugcentral_drug_disease.csv",
         method="Live PostgreSQL query",
         steps=["No raw DrugCentral file has been identified under data/raw/drugcentral. The processor queries PostgreSQL.",
                "Prefer omop_relationship; discover a suitable drug/structures table dynamically.",
                "LEFT JOIN on the discovered structure/drug identifier to obtain the drug name.",
                "If unavailable, use a supported indication table as the fallback source; map native relationship to indication, contraindication or other categories."],
         maps=[("omop_relationship.struct_id", "drugcentral_struct_id; drug_identifier"),
               ("discovered drug table name", "drug_name"),
               ("relationship_name / relationship", "source_relationship_raw; relationship_type"),
               ("snomed_conceptid / concept_id", "disease_or_condition_identifier"),
               ("snomed_full_name / concept_name", "disease_or_condition_name"),
               ("concept_class / concept_type", "concept_kind")],
         caveat="PostgreSQL original headers and exact query path need a snapshot/query log or read-only connection. Root-level FDA archives should NOT be assumed to be DrugCentral input."),
    dict(name="Hetionet", folder="hetionet", parser="hetionet.py",
         inputs=["hetionet-v1.0-nodes.tsv", "hetionet-v1.0-edges.sif.gz"],
         processed="hetionet/hetionet_drug_disease.csv", method="Edge-to-node lookup",
         steps=["Read nodes(id, name, kind) and edges(source, metaedge, target) as tab-separated tables.",
                "Use edge.source and edge.target to look up node names/types using nodes.id.",
                "Keep only compound-to-disease CtD and CpD edges, checking the endpoint node kinds.",
                "Map CtD to treats disease and CpD to palliates disease; retain the original metaedge."],
         maps=[("edges.source", "source_node_id; drug_identifier via node ID"),
               ("edges.target", "target_node_id; disease_identifier via node ID"),
               ("edges.metaedge", "metaedge; relationship_type; evidence_type"),
               ("nodes.id", "node lookup key"), ("nodes.name", "drug_name / disease name"),
               ("nodes.kind", "direction and entity-type filter")],
         caveat="Original graph contains many unrelated edge types; only CtD and CpD are processed."),
    dict(name="MEDI", folder="medi", parser="MEDI.py",
         inputs=["MEDI-2.csv", "MEDI-2_UMLS.csv"],
         processed="medi/medi_drug_disease.csv", method="Two streams concatenated",
         steps=["Read MEDI-2.csv as RxCUI–ICD-9/ICD-10 indication rows, retaining precision and resource flags.",
                "Read MEDI-2_UMLS.csv separately as RxCUI–UMLS indication rows.",
                "Convert each stream to the same output columns and append its rows; no row-wise JOIN between the two input files.",
                "Retain source_record_type so ICD-coded and UMLS-coded records remain distinguishable."],
         maps=[("MEDI-2.csv: RXCUI", "drug_identifier (RxCUI)"),
               ("MEDI-2.csv: DRUG_DESCRIPTION", "drug_name"),
               ("MEDI-2.csv: VOCABULARY / CODE", "disease ID namespace / disease ID"),
               ("MEDI-2.csv: INDICATION_DESCRIPTION", "disease_or_condition_name"),
               ("MEDI-2.csv: HIGH_PRECISION_SUBSET", "high_precision_subset; evidence_type"),
               ("MEDI-2.csv: IN_*", "resources_mentioning"),
               ("MEDI-2_UMLS.csv: RXCUI / DRUG_DESC", "drug ID / name"),
               ("MEDI-2_UMLS.csv: CUI / STR", "UMLS disease ID / disease name")],
         caveat="Concatenating ICD and UMLS streams does not assert that two differently coded records are independent indications."),
    dict(name="OpenTargets", folder="opentargets", parser="opentargets.py",
         inputs=["26.03/clinical_indication/clinical_indication.parquet", "26.03/disease/disease.parquet", "26.03/drug_molecule/*.parquet"],
         processed="opentargets_drug_disease.csv", method="Partition concatenation and ID lookups",
         steps=["Read all five drug_molecule Parquet partitions and concatenate them into one molecule lookup table.",
                "Read disease.parquet and create an ID-to-name lookup from id and name.",
                "Read clinical_indication.parquet as the primary association table.",
                "Resolve clinical drugId through drug_molecule.id; resolve diseaseId through disease.id. Optionally use a ChEMBL-name fallback.",
                "Retain clinicalReportIds, indication ID, maximum clinical stage and source release."],
         maps=[("clinical_indication.drugId", "drug_id; drug_identifier"),
               ("clinical_indication.diseaseId", "disease_id; disease_or_condition_identifier"),
               ("clinical_indication.id", "opentargets_indication_id"),
               ("clinical_indication.maxClinicalStage", "clinical_phase_raw; max_phase"),
               ("clinical_indication.clinicalReportIds", "clinical_report_ids"),
               ("drug_molecule.id / name", "drug name lookup"),
               ("disease.id / name", "disease name lookup")],
         caveat="Drug-molecule partitions are pieces of one table, not five independent evidence sources. Clinical stage is indication-specific."),
    dict(name="PrimeKG", folder="primekg", parser="PrimeKG.py",
         inputs=["kg.csv"], processed="primekg/primekg_drug_disease.csv",
         method="Single-graph edge filtering",
         steps=["Read the kg.csv knowledge-graph edge table.",
                "Keep only indication, contraindication and off-label-use relationships.",
                "Filter to the drug -> disease direction (x_type=drug, y_type=disease), dropping reverse orientation.",
                "Map x-side IDs/names/namespaces to drug and y-side to disease; preserve the original edge relation."],
         maps=[("relation", "primekg_relation; relationship_type; evidence_type"),
               ("display_relation", "primekg_display_relation"),
               ("x_id / x_name / x_source", "drug_identifier / drug_name / drug namespace"),
               ("y_id / y_name / y_source", "disease identifier / name / namespace"),
               ("x_type / y_type", "direction and node-type filters")],
         caveat="MONDO_grouped labels may not denote a single resolvable MONDO identifier; keep source ID and type."),
    dict(name="repoDB", folder="repodb", parser="repoDB.py",
         inputs=["repodb_full_download.tsv"], processed="repodb/repodb_drug_disease.csv",
         method="Single-file record mapping",
         steps=["Read repodb_full_download.tsv; despite its extension, the downloaded file is comma-separated.",
                "Map drugbank_id to drug identifier; ind_id to disease identifier; retain NCT, status, phase and detailed status.",
                "Convert Approved and terminated/withdrawn/suspended statuses into distinct relationship and evidence categories.",
                "The local R scripts/RData are supplementary repository materials; the parser does not JOIN them to the primary file."],
         maps=[("drug_name", "drug_name"), ("drugbank_id", "drug_identifier (DrugBank)"),
               ("ind_name", "disease_or_condition_name"), ("ind_id", "disease_or_condition_identifier (UMLS)"),
               ("NCT", "nct_id"), ("status", "status; relationship_type; evidence_type"),
               ("phase", "phase"), ("DetailedStatus", "detailed_status")],
         caveat="Terminated/withdrawn does not necessarily indicate lack of efficacy. Final merged filtering also removes missing-name records."),
    dict(name="SIDER", folder="sider", parser="sider.py",
         inputs=["drug_names.tsv", "meddra_all_indications.tsv.gz"],
         processed="sider/sider_drug_disease.csv", method="Drug-name lookup by ID",
         steps=["Read drug_names.tsv as a headerless 2-column file: drug_id, drug_name.",
                "Read meddra_all_indications.tsv.gz as a headerless seven-column file; assign meanings by position.",
                "Map each indication drug ID to its drug name using the first file.",
                "Use column 5 for the disease ID, not the MedDRA concept type in column 4; preserve extraction method and LLT/PT level."],
         maps=[("drug_names.tsv col 0 / 1", "drug ID -> drug name lookup"),
               ("indications col 0", "sider_drug_id; drug_identifier"),
               ("indications col 1 / 3", "source_concept_id / source_concept_name"),
               ("indications col 2", "extraction_method"),
               ("indications col 4", "meddra_concept_type"),
               ("indications col 5 / 6", "disease identifier / name")],
         caveat="Headerless fields are positional names ASSIGNED by sider.py, not original file headers. Different extraction methods and MedDRA levels need later scientific validation."),
]

# Short, human-readable verification against the local parser file. These checks
# compare recognizable operations, never execute the parser or print code dumps.
CODE_CHECKS = {
    "AACT": [("Trial join uses nct_id", "USING (nct_id)"),
             ("Default filter selects drug interventions", "lower(i.intervention_type) = 'drug'")],
    "Broad Drug Repurposing Hub": [("Indication strings are expanded", "split_indications(indication)"),
                                   ("Raw file read as TSV", 'sep="\\t"')],
    "ChEMBL": [("drug_indication is primary SQL source", "FROM drug_indication di"),
               ("Molecule names added via LEFT JOIN", "LEFT JOIN molecule_dictionary"),
               ("Disease dictionary conditional", 'table_exists(conn, "disease_dictionary")')],
    "CTD": [("Preserves direct evidence", 'out["evidence_type"]'),
            ("Retains inferred gene field", "InferenceGeneSymbol"),
            ("Reads input by chunks", "chunksize=chunksize")],
    "DrugCentral": [("Uses omop_relationship route", "FROM omop_relationship r"),
                    ("Joins a discovered drug table", "LEFT JOIN")],
    "Hetionet": [("Filters CtD/CpD edges", 'DRUG_DISEASE_METAEDGES'),
                 ("Looks up node names", "id_to_name")],
    "MEDI": [("ICD and UMLS records concatenated", "icd_rows + umls_rows"),
             ("Retains supporting resources", "resources_mentioning")],
    "OpenTargets": [("Combines drug molecule partitions", "pd.concat(dfs"),
                    ("Resolves drug ID through name map", "drug_name_map.get"),
                    ("Resolves disease ID through name map", "disease_name_map.get")],
    "PrimeKG": [("Filters selected relations", "DRUG_DISEASE_RELATIONS"),
                ("Keeps drug-to-disease orientation", '== "drug"')],
    "repoDB": [("Detects delimiter", "detect_separator"),
               ("Maps development status", "classify_relationship(status)")],
    "SIDER": [("Drug name lookup keyed by ID", "drug_name_map.get"),
              ("Reads disease ID from column 5", 'df["raw_col_5"]'),
              ("Preserves extraction method", 'df["raw_col_2"]')],
}

HEADERLESS = {
    "drug_names.tsv": ["drug_id", "drug_name"],
    "meddra_all_indications.tsv.gz": ["drug_id", "source_concept_id",
                                     "extraction_method", "source_concept_name",
                                     "meddra_concept_type", "disease_id", "disease_name"],
}


def snip(value, n=125):
    t = str(value).replace("\r", " ").replace("\n", " ")
    return t if len(t) <= n else t[:n - 3] + "..."


def readable_size(n):
    a = float(n)
    for unit in ("B", "KiB", "MiB", "GiB", "TiB"):
        if a < 1024 or unit == "TiB":
            return ("{:,.0f} {}" if unit == "B" else "{:,.1f} {}").format(a, unit)
        a /= 1024


def grouped(items, width=WIDTH, prefix="      "):
    """Wrap fields at token boundaries; never omit/truncate a header."""
    lines, current = [], prefix
    for token in items:
        token = str(token)
        candidate = token if current == prefix else "  |  " + token
        if current != prefix and len(current + candidate) > width:
            lines.append(current)
            current = prefix + token
        else:
            current += candidate
    if current != prefix:
        lines.append(current)
    return lines or [prefix + "(none)"]


def get_root(explicit=None):
    if explicit:
        p = Path(explicit).expanduser().resolve()
        if (p / "data/raw").is_dir() or (p / "Analysis_00_MergeDatabases.py").exists():
            return p
        raise SystemExit("Root does not contain data/raw or Analysis_00_MergeDatabases.py: " + str(p))
    for p in (Path.cwd(), *Path.cwd().parents):
        if (p / "data/raw").is_dir() or (p / "Analysis_00_MergeDatabases.py").exists():
            return p
    raise SystemExit("Run from inside drug-disease-mapping, or pass --root /path/to/repository")


def inspect_csv(path):
    name = path.name.lower()
    opener = gzip.open if name.endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8-sig", errors="replace", newline="") as f:
        if name == "repurposing_drugs_20200324.txt":
            first = next((s for s in f if s.strip() and not s.startswith("!")), "")
        elif name == "ctd_chemicals_diseases.tsv.gz":
            first = next((s.lstrip("# ").strip() for s in f
                          if "ChemicalName" in s and "DiseaseID" in s), "")
        else:
            first = f.readline()
        if not first:
            return {"kind": "table", "status": "empty file / no header", "fields": [], "sample": {}}
        if name in ("interventions.txt", "conditions.txt", "studies.txt"):
            delim = "|"
        elif name == "repurposing_drugs_20200324.txt":
            delim = "\t"
        else:
            # repoDB file carries .tsv extension but has comma-separated CSV.
            delim = "\t" if "\t" in first else ","
        fields = next(csv.reader([first], delimiter=delim), [])
        if name in HEADERLESS:
            sample_values = fields
            fields = HEADERLESS[name]
            mode = "Headerless: positional schema assigned by parser"
        else:
            mode = "Headered"
            reader = csv.reader(f, delimiter=delim)
            sample_values = next(reader, [])
            if name == "ctd_chemicals_diseases.tsv.gz":
                while sample_values and (len(sample_values) == 1 or sample_values[0].startswith("#")):
                    sample_values = next(reader, [])
        return {"kind": "table", "status": mode,
                "delimiter": {"|": "pipe (|)", "\t": "tab", ",": "comma"}[delim],
                "fields": fields, "sample": dict(zip(fields, sample_values))}


def inspect_parquet(path):
    try:
        import pyarrow.parquet as pq
    except ImportError:
        return {"kind": "parquet", "status": "pyarrow not installed; schema unavailable", "fields": [], "sample": {}}
    p = pq.ParquetFile(path)
    flds = [(f.name, str(f.type)) for f in p.schema_arrow]
    sample = {}
    if p.metadata.num_row_groups and p.metadata.num_rows:
        try:
            columns = [k for k, _ in flds if k in
                       {"id", "drugId", "diseaseId", "name", "maxClinicalStage", "clinicalReportIds", "code", "inchiKey"}]
            if columns:
                values = p.read_row_group(0, columns=columns).slice(0, 1).to_pylist()
                if values:
                    sample = {k: snip(v, 80) for k, v in values[0].items()}
        except Exception as e:
            sample = {"note": "preview unavailable: " + snip(e, 120)}
    return {"kind": "parquet", "status": "Parquet", "rows": p.metadata.num_rows,
            "fields": [k for k, _ in flds], "types": flds, "sample": sample}


def inspect_input(path):
    name = path.name.lower()
    try:
        if name.endswith(".parquet"):
            return inspect_parquet(path)
        if name.endswith(".zip"):
            with zipfile.ZipFile(path) as z:
                members = [n for n in z.namelist() if not n.endswith("/")]
            return {"kind": "archive", "status": "ZIP archive", "members": members,
                    "fields": [], "sample": {}}
        if name.endswith((".tar.gz", ".tar", ".sqlite", ".db", ".rdata", ".rds")):
            return {"kind": "binary", "status": "Archive/database: no single flat-file header",
                    "fields": [], "sample": {}}
        if name.endswith((".csv", ".tsv", ".sif.gz")) or name in (
            "interventions.txt", "conditions.txt", "studies.txt", "repurposing_drugs_20200324.txt",
            "ctd_chemicals_diseases.tsv.gz"):
            return inspect_csv(path)
        return {"kind": "aux", "status": "Auxiliary documentation/code/binary; no table header",
                "fields": [], "sample": {}}
    except Exception as e:
        return {"kind": "error", "status": type(e).__name__ + ": " + snip(e, 160),
                "fields": [], "sample": {}}


def get_processed(path):
    if not path or not path.is_file():
        return None
    try:
        with path.open("r", encoding="utf-8-sig", errors="replace", newline="") as f:
            reader = csv.DictReader(f)
            fields = reader.fieldnames or []
            sample = next(reader, None) or {}
        return {"path": path, "fields": fields, "sample": sample}
    except Exception as e:
        return {"path": path, "error": type(e).__name__ + ": " + snip(e)}


def read_json(path):
    try:
        return json.loads(path.read_text(encoding="utf-8")) if path.is_file() else {}
    except (OSError, ValueError):
        return {}


def find_primary(d, patterns):
    names = set()
    paths = []
    for pat in patterns:
        items = sorted(d.glob(pat))
        for p in items:
            if p.is_file() and p not in names:
                paths.append(p)
                names.add(p)
    return paths


def get_record(root, source, merge_report):
    d = root / "data/raw" / source["folder"]
    allfiles = sorted((p for p in d.rglob("*") if p.is_file()), key=lambda p: str(p)) if d.is_dir() else []
    primary = find_primary(d, source["inputs"])
    merge_entry = next((r for r in merge_report.get("source_reports", [])
                        if r.get("source") == source["name"]), {})
    selected = root / "data/processed" / source["processed"]
    if not selected.exists() and merge_entry.get("selected_file"):
        candidate = Path(merge_entry["selected_file"])
        if candidate.exists():
            selected = candidate
    parser_path = root / source["parser"]
    return {"config": source, "dir": d, "files": allfiles, "primary": primary,
            "inputs": {p: inspect_input(p) for p in primary},
            "processed": get_processed(selected), "merge": merge_entry,
            "parser": parser_path, "parser_exists": parser_path.is_file()}


def check_local_parser(rec):
    path = rec["parser"]
    if not path.is_file():
        return [("Local parser is present", False)]
    try:
        code = path.read_text(encoding="utf-8", errors="replace").casefold()
    except Exception:
        return [("Local parser can be read", False)]
    return [(label, token.casefold() in code)
            for label, token in CODE_CHECKS.get(rec["config"]["name"], [])]


def print_archive_members(record, root):
    # Optional: inspect ZIP members without extracting and without scanning their rows.
    for archive in record["files"]:
        if archive.suffix.lower() != ".zip":
            continue
        print("\n  ZIP ARCHIVE MEMBER INSPECTION:", archive.relative_to(root))
        try:
            with zipfile.ZipFile(archive) as z:
                for item in z.infolist():
                    if item.is_dir():
                        continue
                    print("    •", item.filename)
                    if not item.filename.lower().endswith((".txt", ".tsv", ".csv")):
                        print("       Auxiliary/non-tabular member; no header inspected.")
                        continue
                    try:
                        with z.open(item) as f:
                            first = f.readline(131072).decode("utf-8-sig", "replace").strip()
                        sep = ("|" if "|" in first else
                               "\t" if "\t" in first else
                               "," if "," in first else None)
                        fields = next(csv.reader([first], delimiter=sep), []) if sep else []
                        if len(fields) > 1:
                            print("       Header:")
                            for line in grouped(fields, prefix="         "):
                                print(line)
                        else:
                            print("       No confidently detected tabular header.")
                    except Exception as ex:
                        print("       Header unavailable:", snip(ex, 110))
        except Exception as ex:
            print("    ZIP directory could not be read:", snip(ex, 110))


def print_terminal(records, root, merge_report):
    print("\n" + "═" * WIDTH)
    print("  DRUG–DISEASE SOURCE INSPECTOR | RAW FILES → FIELD MAPPING → PROCESSED SCHEMA")
    print("═" * WIDTH)
    print("  Repository :", root)
    print("  Scope      : inspection of actual local files; no downloads, parsers or merges executed")
    print("  Evidence   : file headers/examples read from disk; workflow described from existing parser design")
    print("  Sample rows from different input files are INDEPENDENT, not joined examples.")
    print("\n  QUICK INVENTORY")
    print("  " + "─" * (WIDTH - 4))
    for i, rec in enumerate(records, 1):
        c = rec["config"]
        print("  {:>2}. {:<29} {:>2} local files   {:>2} processing inputs   {}".format(
            i, c["name"], len(rec["files"]), len(rec["primary"]),
            "processed ✓" if rec["processed"] else "not processed"))
    print("\n  NOTE: archived files count once. AACT ZIP contains additional members; ChEMBL SQL schema")
    print("        needs an extracted SQLite DB; DrugCentral live PostgreSQL fields are not on disk.\n")
    for index, rec in enumerate(records, 1):
        c = rec["config"]
        print("\n" + "═" * WIDTH)
        print("  {:02d} / {:02d}  {}".format(index, len(records), c["name"].upper()))
        print("═" * WIDTH)
        print("  RAW FILES      : {} local files total; {} used by this processor".format(len(rec["files"]), len(rec["primary"])))
        print("  RAW FOLDER     : {}".format(rec["dir"].relative_to(root)))
        print("  PARSER         : {} ({})".format(c["parser"], "present" if rec["parser_exists"] else "not present"))
        print("  COMBINATION    : {}".format(c["method"]))
        print("\n  A. FILE INVENTORY")
        if not rec["files"]:
            print("     No local files in raw directory.")
        for p in rec["files"]:
            role = "INPUT" if p in rec["primary"] else "OTHER"
            print("     [{:<5}] {:<74} {}".format(role, snip(str(p.relative_to(rec["dir"])), 72), readable_size(p.stat().st_size)))
        for j, p in enumerate(rec["primary"], 1):
            x = rec["inputs"][p]
            print("\n  B{}. DOWNLOADED INPUT: {}".format(j, p.relative_to(rec["dir"])))
            print("      FORMAT   : {}".format(x["status"]))
            if x.get("delimiter"):
                print("      DELIMITER: {}".format(x["delimiter"]))
            if "rows" in x:
                print("      ROWS     : {:,} (Parquet metadata)".format(x["rows"]))
            if x["kind"] == "archive":
                print("      CONTENTS : {} member files (not separately counted)".format(len(x["members"])))
                print("      FIRST MEMBERS:")
                for name in x["members"][:14]:
                    print("       - " + name)
                if len(x["members"]) > 14:
                    print("       ... {} more; see --archive-members".format(len(x["members"]) - 14))
            if x.get("fields"):
                print("      HEADER   : {} {}".format(len(x["fields"]), "assigned positional fields" if
                                                   "Headerless" in x["status"] else "original columns"))
                for line in grouped(x["fields"], prefix="         "):
                    print(line)
            if x.get("types"):
                print("      DATA TYPES (Parquet):")
                for field, typ in x["types"]:
                    print("       - {}: {}".format(field, textwrap.shorten(typ, width=93, placeholder="…")))
            if x.get("sample"):
                print("      ONE RAW RECORD (selected fields):")
                for k, v in list(x["sample"].items())[:8]:
                    print("       {:<30} {}".format(snip(k, 29), snip(v or "[empty]", 65)))
                if len(x["sample"]) > 8:
                    print("       ... {} other fields in header".format(len(x["sample"]) - 8))
        if not rec["primary"]:
            print("\n  B. NO LOCAL PROCESSING INPUT AVAILABLE")
            print("     PostgreSQL source must be documented from the original database/query snapshot.")
        missing = [x for x in c["inputs"] if not list(rec["dir"].glob(x))]
        if missing:
            print("\n  EXPECTED BUT MISSING:")
            for s in missing:
                print("     - " + s)
        print("\n  C. HOW THE FILES ARE COMBINED")
        for num, step in enumerate(c["steps"], 1):
            lines = textwrap.wrap(step, width=WIDTH - 10, break_long_words=False)
            for ix, ln in enumerate(lines):
                print(("     {}. ".format(num) if ix == 0 else "        ") + ln)
        print("\n  D. SOURCE → OUTPUT FIELD MAPPING")
        for origin, dest in c["maps"]:
            text = "{}  →  {}".format(origin, dest)
            for ln in textwrap.wrap(text, width=WIDTH - 10, break_long_words=False):
                print("     " + ln)
        print("\n  E. OUTPUT SCHEMA AND ACTUAL OUTPUT RECORD")
        proc = rec["processed"]
        if not proc:
            print("     Not processed yet — rerun this inspector after the source processor finishes.")
        elif "error" in proc:
            print("     Could not read the output header: " + proc["error"])
        else:
            print("     PROCESSED FILE: {}".format(proc["path"].relative_to(root) if
                  proc["path"].is_relative_to(root) else proc["path"]))
            mr = rec["merge"]
            if mr:
                print("     MERGE RECORDS : {:,} source rows → {:,} admitted".format(
                      int(mr.get("rows_original") or 0), int(mr.get("rows_after_core_filter") or 0)))
                if merge_report.get("per_source_rows_deduplicated", {}).get(c["name"]) is not None:
                    print("     AFTER DEDUP  : {:,}".format(merge_report["per_source_rows_deduplicated"][c["name"]]))
            print("     FULL HEADER ({} columns):".format(len(proc["fields"])))
            for line in grouped(proc["fields"], prefix="         "):
                print(line)
            print("     CORE FIELDS:")
            for line in grouped([f for f in proc["fields"] if f in COMMON], prefix="         "):
                print(line)
            print("     SOURCE-SPECIFIC FIELDS:")
            for line in grouped([f for f in proc["fields"] if f not in COMMON], prefix="         "):
                print(line)
            print("     ONE ACTUAL PROCESSED ROW (all columns):")
            for k in proc["fields"]:
                print("       {:<39} {}".format(snip(k, 38), snip(proc["sample"].get(k) or "[empty]", 58)))
        print("\n  F. CHECKS AGAINST LOCAL PROCESSOR (NO CODE DUMPS)")
        for label, passed in check_local_parser(rec):
            print("     {}  {}".format("[FOUND]" if passed else "[REVIEW]", label))
        print("     [NOTE] These are static code-text checks, not proof that the processor ran correctly.")
        print("\n  G. IMPORTANT INTERPRETATION")
        for line in textwrap.wrap(c["caveat"], width=WIDTH - 10, break_long_words=False):
            print("     " + line)
    if len(records) > 1 and merge_report:
        print("\n" + "═" * WIDTH)
        print("  INTEGRATION SUMMARY (FROM EXISTING MERGE JSON, IF AVAILABLE)")
        print("═" * WIDTH)
        r = [x for x in merge_report.get("source_reports", []) if x.get("status") == "loaded"]
        before = sum(int(x.get("rows_original", 0)) for x in r)
        admitted = sum(int(x.get("rows_after_core_filter", 0)) for x in r)
        print("  Original processed source rows : {:,}".format(before))
        print("  After source-level core filters: {:,}".format(admitted))
        print("  Rows removed by these filters : {:,}".format(before-admitted))
        for x in r:
            delta = int(x.get("rows_original", 0)) - int(x.get("rows_after_core_filter", 0))
            if delta:
                print("    - {}: {:,} removed".format(x.get("source"), delta))
        for key in ("rows_merged", "rows_deduplicated", "rows_collapsed_drug_disease_pairs"):
            val = merge_report.get(key)
            if val is not None:
                print("  {:<34} {:,}".format(key + ":", int(val)))
        print("\n  HOW ELEVEN OUTPUT FILES ARE INTEGRATED")
        print("   1. Add missing harmonised fields, standardise source labels and require nonempty names.")
        print("   2. Vertically CONCATENATE the eleven CSV tables; this is not cross-database entity resolution.")
        print("   3. Deduplicate using nine identity/source/relationship/evidence columns.")
        print("   4. The separate collapsed output groups by six name/ID fields; relationship labels are aggregated.")
        print("   5. Distinct trials, papers or source records can be lost from the nine-key deduplicated view.")
    print("\n" + "═" * WIDTH)
    print("  END OF REPORT — NO DATASETS MODIFIED, DOWNLOADED OR REPROCESSED")
    print("═" * WIDTH + "\n")


def md_cell(s):
    return str(s).replace("|", "\\|").replace("\n", " ").replace("\r", " ")


def markdown_report(records, root, merge_report):
    out = ["# Supplementary S0 — Source file and transformation inspection", "",
           "This is an **inspection report**, not a completed validation. Headers and sample values are read from local files. Workflow descriptions reflect the current source parser design. No data was changed.", "",
           "## Source inventory", "",
           "| Database | Local raw files | Processing inputs found | Processed file |",
           "|---|---:|---:|---|"]
    for rec in records:
        name = rec["config"]["name"]
        out.append("| {} | {} | {} | {} |".format(
            md_cell(name), len(rec["files"]), len(rec["primary"]),
            "Present" if rec["processed"] else "Not yet available"))
    for idx, rec in enumerate(records, 1):
        c = rec["config"]
        out += ["", "## {}. {}".format(idx, c["name"]), "",
                "**Local files:** {}. **Processing inputs:** {}. **Processor:** `{}`.".format(
                    len(rec["files"]), len(rec["primary"]), c["parser"]), "",
                "### Files and raw schemas", "",
                "| File | Role | Format | Size |", "|---|---|---|---:|"]
        for p in rec["files"]:
            role = "Processing input" if p in rec["primary"] else "Auxiliary / archive"
            x = rec["inputs"].get(p, {})
            out.append("| `{}` | {} | {} | {} |".format(
                md_cell(str(p.relative_to(rec["dir"]))), role,
                md_cell(x.get("status", "Not inspected")), readable_size(p.stat().st_size)))
        if not rec["files"]:
            out.append("| No local raw files | PostgreSQL input | Schema not locally accessible | — |")
        for p in rec["primary"]:
            x = rec["inputs"][p]
            out += ["", "**Input:** `{}`".format(p.relative_to(rec["dir"])), ""]
            if x.get("fields"):
                out += ["{} columns ({}):".format(len(x["fields"]), x.get("delimiter", x["status"])), "",
                        "```text", ", ".join(x["fields"]), "```", ""]
            elif x["kind"] == "binary":
                out += ["Binary/database archive. The underlying SQL table schemas have not been inspected.", ""]
            if x.get("types"):
                out += ["**Parquet field types:**", ""]
                for key, typ in x["types"]:
                    out.append("- `{}`: `{}`".format(key, typ))
                out.append("")
            if x.get("sample"):
                out += ["**Actual raw-record example (selected fields):**", "",
                        "| Field | Value |", "|---|---|"]
                for k, v in list(x["sample"].items())[:8]:
                    out.append("| `{}` | {} |".format(md_cell(k), md_cell(snip(v or "[empty]", 95))))
        out += ["", "### How the source information is combined", ""]
        for i, step in enumerate(c["steps"], 1):
            out.append("{}. {}".format(i, step))
        out += ["", "### Source-to-output field mapping", "",
                "| Original field | Harmonised/output field |", "|---|---|"]
        for k, v in c["maps"]:
            out.append("| `{}` | {} |".format(md_cell(k), md_cell(v)))
        out += ["", "### Final processed header and example", ""]
        p = rec["processed"]
        if not p:
            out.append("Not yet available; re-inspect after source processing.")
        elif "error" in p:
            out.append("Could not read processed file: " + p["error"])
        else:
            out += ["**Processed file:** `{}`".format(p["path"].name), "",
                    "**Complete header ({} columns):**".format(len(p["fields"])), "",
                    "```text", ", ".join(p["fields"]), "```", "",
                    "**Actual processed record:**", "",
                    "| Field | Value |", "|---|---|"]
            for k in p["fields"]:
                out.append("| `{}` | {} |".format(md_cell(k), md_cell(snip(p["sample"].get(k) or "[empty]", 105))))
        out += ["", "**Checks against the local processor (static only):**", ""]
        for label, passed in check_local_parser(rec):
            out.append("- {} — {}".format("Located" if passed else "Needs verification", label))
        out += ["", "**Interpretation/limitation:** " + c["caveat"]]
    if len(records) > 1 and merge_report:
        rows = [r for r in merge_report.get("source_reports", []) if r.get("status") == "loaded"]
        before = sum(int(r.get("rows_original", 0)) for r in rows)
        admitted = sum(int(r.get("rows_after_core_filter", 0)) for r in rows)
        out += ["", "## Source integration (existing release)", "",
                "The processed source CSVs are vertically concatenated after common-field harmonisation and required-name filtering. This is not identity resolution across chemical/disease ontologies.", "",
                "| Stage | Rows |", "|---|---:|",
                "| Sum of processed source files | {:,} |".format(before),
                "| After required-name filter | {:,} |".format(admitted),
                "| Final deduplicated assertions | {:,} |".format(merge_report.get("rows_deduplicated", 0)),
                "| Separate exact-schema collapsed view | {:,} |".format(merge_report.get("rows_collapsed_drug_disease_pairs", 0)),
                "", "The nine-field deduplication key omits many original evidence fields. The six-field collapse aggregates relationship labels and should not be described as ontology-resolved matching."]
    return "\n".join(out) + "\n"


def inspect_sqlite(path):
    if not path.is_file():
        print("\nSQLite file not found:", path)
        return
    print("\n  CHEmBL SQLITE DATABASE (READ-ONLY SCHEMA)")
    connection = None
    try:
        uri = "file:" + quote(str(path.resolve()), safe="/:") + "?mode=ro&immutable=1"
        connection = sqlite3.connect(uri, uri=True)
        for name in ("drug_indication", "molecule_dictionary", "disease_dictionary"):
            found = connection.execute("SELECT name FROM sqlite_master WHERE type=\"table\" AND name=?", (name,)).fetchone()
            if not found:
                print("    {}: NOT PRESENT".format(name))
                continue
            print("    TABLE: {}".format(name))
            rows = connection.execute('PRAGMA table_info("{}")'.format(name)).fetchall()
            for r in rows:
                print("     - {} [{}]".format(r[1], r[2]))
    except Exception as e:
        print("    SQLite inspection failed: " + snip(e))
    finally:
        if connection is not None:
            connection.close()


def main():
    ap = argparse.ArgumentParser(description="Readably inspect raw source files, their headers and mapping, and processed outputs")
    ap.add_argument("--root", help="Path to drug-disease-mapping repository (otherwise auto-detected)")
    ap.add_argument("--source", help="Only one source, e.g. AACT, CTD, ChEMBL, OpenTargets, repoDB")
    ap.add_argument("--markdown", metavar="FILE", help="Additionally create a readable Markdown table report (the only optional write)")
    ap.add_argument("--archive-members", action="store_true", help="Print all file names inside ZIP archives (no extraction)")
    ap.add_argument("--sqlite", help="Optional path to existing extracted ChEMBL SQLite database; schema only")
    args = ap.parse_args()
    root = get_root(args.root)
    chosen = SOURCES
    if args.source:
        chosen = [s for s in SOURCES if s["name"].lower().replace(" ", "") == args.source.lower().replace(" ", "")]
        if not chosen:
            ap.error("Unknown source. Choices: " + ", ".join(s["name"] for s in SOURCES))
    merge_report = read_json(root / "data/processed/ALL_INCLUDED_SOURCES_merge_summary.json")
    records = [get_record(root, s, merge_report) for s in chosen]
    print_terminal(records, root, merge_report)
    if args.archive_members:
        for r in records:
            print_archive_members(r, root)
    if args.sqlite:
        path = Path(args.sqlite).expanduser()
        inspect_sqlite(path if path.is_absolute() else root / path)
    if args.markdown:
        target = Path(args.markdown).expanduser()
        target.write_text(markdown_report(records, root, merge_report), encoding="utf-8")
        print("\nSupplementary inspection Markdown saved to:", target.resolve())


if __name__ == "__main__":
    main()
