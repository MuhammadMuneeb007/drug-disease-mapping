```markdown
# Drug-Disease Association Data Pipeline

A comprehensive pipeline for downloading, processing, and merging drug-disease associations from multiple public biomedical databases. This repository integrates data from 7 major sources to create a unified drug-disease association dataset for pharmaceutical research and computational biology applications.

## Data (Zenodo)

- **DOI:** https://doi.org/10.5281/zenodo.18308460  
- **Zenodo record:** https://zenodo.org/records/18308460  

---

## Table of Contents

- [Overview](#overview)
- [Data Sources](#data-sources)
- [Pipeline Architecture](#pipeline-architecture)
- [File Structure](#file-structure)
- [Data Processing Scripts](#data-processing-scripts)
- [Output Files](#output-files)
- [Merging Methodology](#merging-methodology)
- [Usage](#usage)
- [Requirements](#requirements)
- [License & Citations](#license--citations)
- [Contact & Support](#contact--support)

---

## Overview

This pipeline systematically:
1. **Downloads** drug-disease association data from 7 major public databases
2. **Processes** each source into a standardized schema
3. **Merges** all sources into unified datasets
4. **Deduplicates** entries while preserving source provenance
5. **Generates** analysis-ready CSV files for downstream research

**Final Dataset Statistics:**
- **5,107,064 rows** in merged dataset
- **4,148,850 rows** after deduplication
- **7 data sources** integrated
- **688,747 associations** from clinical trials (AACT)
- **3,564,540 associations** from chemical-disease databases (CTD)

---

## Data Sources

### 1. **AACT (Aggregate Analysis of ClinicalTrials.gov)**
- **File:** `aact_drug_disease.csv` (688,747 rows)
- **Source:** Clinical trials registry data
- **Content:** Drug interventions and their associated conditions from clinical trials
- **Columns:** `nct_id`, `drug_name`, `disease_name`, `intervention_type`, `study_phase`, `overall_status`

### 2. **ChEMBL**
- **File:** `chembl.csv` (59,787 rows)
- **Source:** European Bioinformatics Institute (EBI)
- **Content:** Manually curated drug-disease indications with clinical phase information
- **Columns:** `drug_id`, `drug_name`, `disease_id`, `disease_name`, `mesh_id`, `mesh_heading`, `max_phase`

### 3. **CTD (Comparative Toxicogenomics Database)**
- **File:** `ctd_drug_disease.csv` (3,564,540 rows)
- **Source:** Curated and inferred chemical-disease associations
- **Content:** Chemical-disease relationships with evidence types and inference scores
- **Columns:** `drug_name`, `drug_id`, `disease_name`, `disease_id`, `direct_evidence`, `inference_score`, `source`

### 4. **DrugCentral**
- **File:** `drugcentral_drug_disease.csv` (38,431 rows)
- **Source:** University of New Mexico DrugCentral database
- **Content:** Curated drug-disease indications and contraindications
- **Columns:** `drug_id`, `drug_name`, `disease_id`, `disease_name`, `relationship`, `concept_kind`, `source`

### 5. **Open Targets Platform**
- **File:** `opentargets_drug_disease.csv` (65,198 rows)
- **Source:** Open Targets Platform integration of multiple sources
- **Content:** Drug-disease associations with clinical phase and reference information
- **Columns:** `drug_id`, `drug_name`, `disease_id`, `disease_name`, `max_phase`, `references`, `source`

### 6. **SIDER (Side Effect Resource)**
- **File:** `sider_drug_indication.csv` (19,002 rows)
- **Source:** Drug indications from package inserts
- **Content:** Approved drug indications extracted from FDA labels
- **Columns:** `drug_id`, `drug_name`, `indication_id`, `indication_name`, `source`

### 7. **ClinicalTrials.gov (Direct)**
- **File:** `ctgov_drug_disease.csv` (671,365 rows)
- **Source:** Direct processing of ClinicalTrials.gov data
- **Content:** Drug interventions mapped to clinical trial conditions
- **Columns:** `nct_id`, `drug_name`, `disease_name`, `intervention_type`, `study_phase`, `overall_status`

---

## Pipeline Architecture

```

┌─────────────────────────────────────────────────────────────────┐
│                    DATA ACQUISITION LAYER                        │
├─────────────────────────────────────────────────────────────────┤
│  aact.py │ chembl.py │ ctd.py │ drugcentral.py │ opentargets.py │
│                    sider.py                                      │
└────────────────────┬────────────────────────────────────────────┘
│
▼
┌─────────────────────────────────────────────────────────────────┐
│                  STANDARDIZATION LAYER                           │
├─────────────────────────────────────────────────────────────────┤
│  • Field normalization (drug_name, disease_name, IDs)           │
│  • Source tracking (source, internal_source)                    │
│  • Schema harmonization                                         │
└────────────────────┬────────────────────────────────────────────┘
│
▼
┌─────────────────────────────────────────────────────────────────┐
│                    MERGING LAYER                                 │
├─────────────────────────────────────────────────────────────────┤
│  mergeDrugs.py - Vertical concatenation with Polars             │
│  • Text normalization (lowercase, strip whitespace)             │
│  • Source annotation preservation                               │
└────────────────────┬────────────────────────────────────────────┘
│
├──────────────────────────────────┐
▼                                  ▼
┌────────────────────────┐      ┌─────────────────────────┐
│   MERGED DATASET       │      │  DEDUPLICATED DATASET   │
├────────────────────────┤      ├─────────────────────────┤
│ 5,107,064 rows         │      │ 4,148,850 rows          │
│ All associations       │      │ Unique combinations     │
└────────────────────────┘      └─────────────────────────┘

```

---

## File Structure

```

AllDiseasesToDrugs/
├─ DownloadAllDrugs/
│  ├─ aact.py
│  ├─ chembl.py
│  ├─ ctd.py
│  ├─ drugcentral.py
│  ├─ opentargets.py
│  ├─ sider.py
│  └─ (downloaded raw source files / caches as applicable)
├─ mergeDrugs.py
├─ AllDiseasesToDrugs.py
├─ README.md
└─ output/
├─ aact_drug_disease.csv
├─ chembl.csv
├─ ctd_drug_disease.csv
├─ ctgov_drug_disease.csv
├─ drugcentral_drug_disease.csv
├─ opentargets_drug_disease.csv
├─ sider_drug_indication.csv
├─ ALL_SOURCES_drug_disease_merged.csv
├─ ALL_SOURCES_drug_disease_deduplicated.csv
└─ migraine_drugs.csv

```

---

## Data Processing Scripts

### Downloaders (per-source)
- `DownloadAllDrugs/aact.py` — downloads/derives trial intervention–condition mappings from AACT / ClinicalTrials.gov
- `DownloadAllDrugs/chembl.py` — extracts drug–indication mappings and phase information from ChEMBL SQLite
- `DownloadAllDrugs/ctd.py` — processes CTD chemical–disease associations (curated + inferred)
- `DownloadAllDrugs/drugcentral.py` — queries DrugCentral (PostgreSQL-backed) for indications/contraindications
- `DownloadAllDrugs/opentargets.py` — pulls Open Targets drug–disease associations (may use `chembl.csv` for name harmonization)
- `DownloadAllDrugs/sider.py` — parses SIDER indications extracted from labels

### Merging
- `mergeDrugs.py` — standardizes, concatenates, annotates sources, and exports merged + deduplicated CSVs

### Bulk pipeline runner
- `AllDiseasesToDrugs.py` — end-to-end orchestration:
  1. download all sources
  2. normalize/standardize
  3. generate unified SQLite (if configured)
  4. export CSV/Parquet outputs

---

## Output Files

### Primary Outputs

#### **ALL_SOURCES_drug_disease_merged.csv** (5,107,064 rows)

```

Purpose:     Complete merged dataset with all associations
Schema:      drug_name, disease_name, source, internal_source, drug_id,
disease_id, extra_metadata
Features:

* Preserves all source-specific metadata in JSON format
* Includes duplicate associations from multiple sources
* Annotates each row with source provenance

```

#### **ALL_SOURCES_drug_disease_deduplicated.csv** (4,148,850 rows)

```

Purpose:     Deduplicated dataset for analysis
Schema:      drug_name, drug_id, disease, source, internal_source
Features:

* Unique (drug_name, drug_id, disease) combinations
* Removes ~20% redundancy from merged dataset
* Maintains first-occurrence source attribution

```

### Source-Specific Files

| File | Rows | Source Database |
|------|------|-----------------|
| `aact_drug_disease.csv` | 688,747 | AACT/ClinicalTrials.gov |
| `chembl.csv` | 59,787 | ChEMBL |
| `ctd_drug_disease.csv` | 3,564,540 | CTD |
| `ctgov_drug_disease.csv` | 671,365 | ClinicalTrials.gov Direct |
| `drugcentral_drug_disease.csv` | 38,431 | DrugCentral |
| `opentargets_drug_disease.csv` | 65,198 | Open Targets |
| `sider_drug_indication.csv` | 19,002 | SIDER |

### Domain-Specific Subsets

#### **migraine_drugs.csv** (5,049 rows)

```

Purpose:     Migraine-specific drug associations
Schema:      drug_name, drug_id, n_migraine_rows, unique_disease_terms,
sources, internal_sources
Features:

* Aggregates all migraine-related associations
* Counts occurrences across sources
* Lists all associated sources per drug

```

---

## Data Quality Metrics

```

Total Sources:              7
Original Row Count:         14,433,231 (all individual CSVs)
After Source Processing:    5,107,064
After Deduplication:        4,148,850
Deduplication Rate:         18.8%
Unique Drugs:               ~200,000
Unique Diseases:            ~50,000
Sources Per Drug (avg):     2.3

````

---

## Merging Methodology

1. **Standardize schemas**
   - Normalize column names and types across sources
   - Ensure presence of core fields (drug, disease, source provenance)

2. **Normalize text**
   - Lowercase drug/disease strings
   - Strip whitespace / clean common formatting issues

3. **Preserve provenance**
   - Add `source` and `internal_source` annotations
   - Store source-specific extra fields in `extra_metadata` (JSON) where applicable

4. **Merge**
   - Vertical concatenation of standardized tables (Polars)

5. **Deduplicate**
   - Collapse duplicates to unique `(drug_name, drug_id, disease)` (or equivalent harmonized key)
   - Keep first-occurrence provenance fields for traceability

---

## Usage

### Running Individual Downloaders

```bash
# Download AACT data
cd DownloadAllDrugs
python aact.py --out aact_drug_disease.csv

# Download ChEMBL data (requires chembl_XX_sqlite.tar.gz)
python chembl.py /path/to/chembl_36_sqlite.tar.gz --out chembl.csv

# Download CTD data
python ctd.py --out ctd_drug_disease.csv

# Download DrugCentral data (requires PostgreSQL connection)
python drugcentral.py --out drugcentral_drug_disease.csv

# Download Open Targets data (requires chembl.csv for drug names)
python opentargets.py --chembl chembl.csv --out opentargets_drug_disease.csv

# Download SIDER data
python sider.py --out sider_drug_indication.csv
````

### Merging All Sources

```bash
# Ensure all source CSV files are in current directory
python mergeDrugs.py

# Output: ALL_SOURCES_drug_disease_deduplicated.csv
```

### Complete Pipeline (Bulk Download + Merge)

```bash
# Run full pipeline
python AllDiseasesToDrugs.py ./output_directory

# This will:
# 1. Download all sources
# 2. Process and normalize
# 3. Create unified SQLite database
# 4. Export CSV and Parquet files
```

---

## Requirements

### Python Dependencies

```bash
pip install requests pandas polars duckdb psycopg2-binary pyarrow
```

### Individual Package Requirements

* **aact.py:** `duckdb`, `requests`
* **chembl.py:** `sqlite3`, `pandas`, `tarfile`
* **ctd.py:** `pandas`, `requests`, `gzip`
* **drugcentral.py:** `psycopg2-binary`
* **opentargets.py:** `requests`, `pandas`, `pyarrow`
* **sider.py:** `pandas`, `requests`
* **mergeDrugs.py:** `polars`

---

## License & Citations

### Data Licenses

Each data source has its own license terms:

* **AACT/ClinicalTrials.gov:** Public domain (U.S. Government)
* **ChEMBL:** Creative Commons Attribution-ShareAlike 3.0 Unported License
* **CTD:** Free for academic and commercial use with attribution
* **DrugCentral:** Creative Commons Attribution-ShareAlike 4.0
* **Open Targets:** Creative Commons Attribution 4.0
* **SIDER:** Creative Commons Attribution-NonCommercial-ShareAlike 4.0
* **ClinicalTrials.gov:** Public domain

### Dataset citation (Zenodo)

```bibtex
@dataset{muneeb_drug_disease_2026,
  author       = {Muneeb, Muhammad},
  title        = {Drug-Disease Association Data Pipeline and Unified Dataset},
  year         = {2026},
  publisher    = {Zenodo},
  doi          = {10.5281/zenodo.18308460},
  url          = {https://doi.org/10.5281/zenodo.18308460}
}
```

### Citations (data sources)

Please cite the original data sources when using this pipeline:

```bibtex
@article{chembl2024,
  title   = {ChEMBL: towards direct deposition of bioassay data},
  journal = {Nucleic Acids Research},
  year    = {2024}
}

@article{ctd2023,
  title   = {The Comparative Toxicogenomics Database: update 2023},
  journal = {Nucleic Acids Research},
  year    = {2023}
}

@article{drugcentral2023,
  title   = {DrugCentral 2023 extends human clinical data and integrates veterinary drugs},
  journal = {Nucleic Acids Research},
  year    = {2023}
}

@article{opentargets2023,
  title   = {Open Targets Platform: supporting systematic drug-target identification and prioritisation},
  journal = {Nucleic Acids Research},
  year    = {2023}
}

@article{sider2016,
  title   = {The SIDER database of drugs and side effects},
  journal = {Nucleic Acids Research},
  year    = {2016}
}
```

---

## Contact & Support

For questions, issues, or contributions, please refer to the project repository or contact the research team.

**Muhammad Muneeb**

* Email: [muneebsiddique007@gmail.com](mailto:muneebsiddique007@gmail.com)
* Email: [m.muneeb@uq.edu.au](mailto:m.muneeb@uq.edu.au)
* GitHub: [https://github.com/MuhammadMuneeb007/drug-disease-mapping](https://github.com/MuhammadMuneeb007/drug-disease-mapping)

```
::contentReference[oaicite:0]{index=0}
```
