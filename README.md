# Drug-Disease Association Data Pipeline

A comprehensive pipeline for downloading, processing, and merging drug-disease associations from multiple public biomedical databases. This repository integrates data from 7 major sources to create a unified drug-disease association dataset for pharmaceutical research and computational biology applications.

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
.
├── README.md                           # This file
├── AllDiseasesToDrugs.py              # Master orchestration script (bulk download)
├── DownloadAllDrugs/                  # Individual source downloaders
│   ├── aact.py                        # AACT/ClinicalTrials.gov downloader
│   ├── chembl.py                      # ChEMBL SQLite processor
│   ├── ctd.py                         # CTD chemical-disease downloader
│   ├── drugcentral.py                 # DrugCentral PostgreSQL connector
│   ├── opentargets.py                 # Open Targets Platform downloader
│   ├── sider.py                       # SIDER indication processor
│   └── mergeDrugs.py                  # Multi-source merger with Polars
├── Data Files (Generated)
│   ├── aact_drug_disease.csv          # 688,747 rows
│   ├── chembl.csv                     # 59,787 rows
│   ├── ctd_drug_disease.csv           # 3,564,540 rows
│   ├── ctgov_drug_disease.csv         # 671,365 rows
│   ├── drugcentral_drug_disease.csv   # 38,431 rows
│   ├── opentargets_drug_disease.csv   # 65,198 rows
│   ├── sider_drug_indication.csv      # 19,002 rows
│   ├── ALL_SOURCES_drug_disease_merged.csv           # 5,107,064 rows
│   ├── ALL_SOURCES_drug_disease_deduplicated.csv     # 4,148,850 rows
│   └── migraine_drugs.csv             # 5,049 rows (disease-specific subset)
└── Analysis Scripts
    ├── predict4-GeneDifferentialAnalysisSixMethods.py
    ├── predict4.1.2.11-DrugFinder.py
    └── ... (additional analysis pipelines)
```

---

## Data Processing Scripts

### Download Scripts (`DownloadAllDrugs/`)

#### **1. aact.py**
```
Purpose:     Download and process AACT flat files from ClinicalTrials.gov
Input:       AACT flat-files ZIP (interventions.txt, conditions.txt, studies.txt)
Processing:  
  - Downloads latest AACT snapshot
  - Extracts interventions and conditions tables
  - Joins on nct_id using DuckDB
  - Filters to intervention_type='Drug' (optional)
Output:      aact_drug_disease.csv
Columns:     nct_id, drug_name, disease_name, intervention_type, 
             study_phase, overall_status
```

#### **2. chembl.py**
```
Purpose:     Extract drug-disease indications from ChEMBL SQLite database
Input:       chembl_XX_sqlite.tar.gz
Processing:  
  - Auto-detects SQLite DB in tar.gz
  - Introspects schema for column names
  - Joins drug_indication → molecule_dictionary → disease tables
  - Handles schema variations (ChEMBL 36+ compatibility)
  - Exports in chunks for memory efficiency
Output:      chembl.csv
Columns:     drug_id, drug_name, disease_id, disease_name, mesh_id, 
             mesh_heading, max_phase
```

#### **3. ctd.py**
```
Purpose:     Download and parse CTD chemical-disease associations
Input:       CTD_chemicals_diseases.tsv.gz (aggregate or curated)
Processing:  
  - Downloads from CTD FTP
  - Auto-detects header line (handles comment blocks)
  - Normalizes column names (removes '# ' prefix)
  - Processes in chunks (default 250k rows)
  - Preserves inference scores and direct evidence flags
Output:      ctd_drug_disease.csv
Columns:     drug_name, drug_id, disease_name, disease_id, 
             direct_evidence, inference_score, source
```

#### **4. drugcentral.py**
```
Purpose:     Connect to DrugCentral PostgreSQL and export relationships
Input:       DrugCentral public PostgreSQL instance
Processing:  
  - Connects to unmtid-shinyapps.net PostgreSQL
  - Introspects schema to find relationship tables
  - Discovers drug table (structures/drug/drugs) and join keys
  - Builds query for omop_relationship or indication tables
  - Joins to populate drug_name
  - Filters relationship types (indication, contraindication, etc.)
Output:      drugcentral_drug_disease.csv
Columns:     drug_id, drug_name, disease_id, disease_name, relationship, 
             concept_kind, source
```

#### **5. opentargets.py**
```
Purpose:     Download Open Targets drug-disease associations
Input:       Open Targets Platform Parquet files
Processing:  
  - Lists available Parquet files from FTP
  - Downloads drug indication datasets
  - Loads ChEMBL drug names from local chembl.csv
  - Flattens nested indications structure
  - Maps CHEMBL IDs to drug names
  - Extracts max_phase and references
Output:      opentargets_drug_disease.csv
Columns:     drug_id, drug_name, disease_id, disease_name, max_phase, 
             references, source
```

#### **6. sider.py**
```
Purpose:     Download and process SIDER drug indication data
Input:       drug_names.tsv, meddra_all_indications.tsv.gz
Processing:  
  - Downloads from SIDER EMBL server
  - Reads drug names mapping (STITCH → name)
  - Parses MedDRA indication file
  - Handles variable column positions
  - Joins drug IDs to names
  - Deduplicates by drug-indication pair
Output:      sider_drug_indication.csv
Columns:     drug_id, drug_name, indication_id, indication_name, source
```

#### **7. mergeDrugs.py**
```
Purpose:     Merge all source files into unified datasets
Input:       All 7 source CSV files
Processing:  
  1. Schema Normalization
     - Maps source-specific columns to standard schema
     - Creates drug_name, drug_id, disease, source, internal_source
  
  2. Text Normalization
     - Converts to lowercase
     - Strips leading/trailing whitespace
     - Standardizes drug and disease names
  
  3. Vertical Concatenation
     - Uses Polars for high-performance merging
     - Concatenates all sources vertically
  
  4. Deduplication
     - Removes duplicate (drug_name, drug_id, disease) tuples
     - Keeps first occurrence (preserves source priority)
  
Output:      
  - ALL_SOURCES_drug_disease_deduplicated.csv (4,148,850 rows)
  
Technology:  Polars DataFrame library for memory-efficient processing
```

### Master Orchestration Script

#### **AllDiseasesToDrugs.py**
```
Purpose:     Complete end-to-end bulk download and unification
Input:       None (downloads from public APIs)
Processing:  
  - Orchestrates downloads from 13+ sources (7 core + 6 extended)
  - Calls individual downloaders with error handling
  - Parses heterogeneous file formats (TSV, CSV, Parquet, SQLite, ZIP)
  - Creates unified schema across all sources
  - Exports to CSV, Parquet, and SQLite
  - Creates database indices for efficient querying
  
Extended Sources (beyond core 7):
  - TTD (Therapeutic Target Database)
  - BioSNAP Drug-Disease Network
  - KEGG Drug Database
  - Drugs@FDA
  - DisGeNET
  - DGIdb (Drug-Gene Interaction Database)
  
Output:      unified_drug_disease.csv, unified_drug_disease.sqlite
```

---

## Output Files

### Primary Outputs

#### **ALL_SOURCES_drug_disease_merged.csv** (5,107,064 rows)
```
Purpose:     Complete merged dataset with all associations
Schema:      drug_name, disease_name, source, internal_source, drug_id, 
             disease_id, extra_metadata
Features:    
  - Preserves all source-specific metadata in JSON format
  - Includes duplicate associations from multiple sources
  - Annotates each row with source provenance
```

#### **ALL_SOURCES_drug_disease_deduplicated.csv** (4,148,850 rows)
```
Purpose:     Deduplicated dataset for analysis
Schema:      drug_name, drug_id, disease, source, internal_source
Features:    
  - Unique (drug_name, drug_id, disease) combinations
  - Removes ~20% redundancy from merged dataset
  - Maintains first-occurrence source attribution
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
  - Aggregates all migraine-related associations
  - Counts occurrences across sources
  - Lists all associated sources per drug
```

---

## Merging Methodology

### Phase 1: Schema Normalization

Each source is normalized to a common schema:

```python
Standard Schema:
  - drug_name       (str)  : Standardized drug name
  - drug_id         (str)  : Source-specific drug identifier
  - disease         (str)  : Standardized disease/condition name
  - source          (str)  : Primary source (AACT, ChEMBL, CTD, etc.)
  - internal_source (str)  : Detailed source provenance
```

**Normalization Rules:**
1. **Drug Names:** Extracted from source-specific columns
   - `aact` → `drug_name`
   - `chembl` → `drug_name` (from molecule_dictionary)
   - `ctd` → `ChemicalName`
   - `drugcentral` → `name` (from structures table)
   - `opentargets` → mapped from ChEMBL IDs
   - `sider` → `drug_name` (from drug_names.tsv)

2. **Disease Names:** Mapped from various columns
   - `aact` → `disease_name` (condition)
   - `chembl` → `mesh_heading` / `disease_name`
   - `ctd` → `DiseaseName`
   - `drugcentral` → disease concept name
   - `opentargets` → `disease_name`
   - `sider` → `indication_name` (MedDRA)

3. **Identifiers:** Preserved when available
   - Drug IDs: ChEMBL IDs, PubChem CIDs, STITCH IDs
   - Disease IDs: MeSH, EFO, UMLS, OMIM

### Phase 2: Text Normalization

```python
Transformations Applied:
  1. Convert to lowercase
  2. Strip leading/trailing whitespace
  3. Normalize Unicode characters
  4. Remove extra internal spaces
```

### Phase 3: Vertical Concatenation

```python
Method: Polars.concat(dfs, how="vertical")

Process:
  1. Load all 7 source DataFrames
  2. Ensure schema compatibility
  3. Concatenate vertically (row-wise union)
  4. Result: Single DataFrame with 5,107,064 rows
```

### Phase 4: Deduplication

```python
Method: Polars.unique(subset=["drug_name", "drug_id", "disease"], keep="first")

Logic:
  - Identifies duplicate (drug, disease) pairs
  - Keeps first occurrence (source priority implicit)
  - Removes 958,214 duplicate rows (18.8% reduction)
  - Preserves source attribution for kept rows

Result: 4,148,850 unique associations
```

### Source Priority (Implicit in "keep first")

Order of concatenation determines priority:
1. AACT (Clinical trial evidence)
2. ClinicalTrials.gov Direct
3. ChEMBL (Curated bioactivity)
4. CTD (Toxicogenomics)
5. DrugCentral (Regulatory indications)
6. Open Targets (Multi-source integration)
7. SIDER (Package insert data)

**Rationale:** Clinical trial and curated sources prioritized over inferred associations.

### Data Quality Metrics

```
Total Sources:              7
Original Row Count:         14,433,231 (all individual CSVs)
After Source Processing:    5,107,064
After Deduplication:        4,148,850
Deduplication Rate:         18.8%
Unique Drugs:               ~200,000
Unique Diseases:            ~50,000
Sources Per Drug (avg):     2.3
```

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
```

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

### Querying the Data

```python
# Python example
import polars as pl

df = pl.read_csv("ALL_SOURCES_drug_disease_deduplicated.csv")

# Find all diseases for aspirin
aspirin = df.filter(pl.col("drug_name").str.contains("aspirin"))
print(aspirin)

# Count associations per source
source_counts = df.group_by("source").count()
print(source_counts)
```

```sql
-- SQLite example (if using AllDiseasesToDrugs.py output)
sqlite3 unified_drug_disease.sqlite

-- Find all diseases for a drug
SELECT disease_name, source, max_phase 
FROM drug_disease 
WHERE drug_name LIKE '%aspirin%' 
ORDER BY max_phase DESC;

-- Count associations by source
SELECT source, COUNT(*) as n 
FROM drug_disease 
GROUP BY source 
ORDER BY n DESC;
```

---

## Requirements

### Python Dependencies

```bash
pip install requests pandas polars duckdb psycopg2-binary pyarrow
```

### Individual Package Requirements

- **aact.py:** `duckdb`, `requests`
- **chembl.py:** `sqlite3`, `pandas`, `tarfile`
- **ctd.py:** `pandas`, `requests`, `gzip`
- **drugcentral.py:** `psycopg2-binary`
- **opentargets.py:** `requests`, `pandas`, `pyarrow`
- **sider.py:** `pandas`, `requests`
- **mergeDrugs.py:** `polars`

### System Requirements

- **Memory:** 16GB+ RAM recommended (for processing large datasets)
- **Disk Space:** 50GB+ for raw downloads and processed files
- **Network:** Stable internet connection for bulk downloads

---

## License & Citations

### Data Licenses

Each data source has its own license terms:

- **AACT/ClinicalTrials.gov:** Public domain (U.S. Government)
- **ChEMBL:** Creative Commons Attribution-ShareAlike 3.0 Unported License
- **CTD:** Free for academic and commercial use with attribution
- **DrugCentral:** Creative Commons Attribution-ShareAlike 4.0
- **Open Targets:** Creative Commons Attribution 4.0
- **SIDER:** Creative Commons Attribution-NonCommercial-ShareAlike 4.0
- **ClinicalTrials.gov:** Public domain

### Citations

Please cite the original data sources when using this pipeline:

```bibtex
@article{chembl2024,
  title={ChEMBL: towards direct deposition of bioassay data},
  journal={Nucleic Acids Research},
  year={2024}
}

@article{ctd2023,
  title={The Comparative Toxicogenomics Database: update 2023},
  journal={Nucleic Acids Research},
  year={2023}
}

@article{drugcentral2023,
  title={DrugCentral 2023 extends human clinical data and integrates veterinary drugs},
  journal={Nucleic Acids Research},
  year={2023}
}

@article{opentargets2023,
  title={Open Targets Platform: supporting systematic drug-target identification and prioritisation},
  journal={Nucleic Acids Research},
  year={2023}
}

@article{sider2016,
  title={The SIDER database of drugs and side effects},
  journal={Nucleic Acids Research},
  year={2016}
}
```

---

## Contact & Support

For questions, issues, or contributions, please refer to the project repository or contact the research team.

**Repository Maintainer:** PhD Research Project - Migraine Drug Discovery

**Last Updated:** January 2026

---

## Appendix: Field Mappings

### Merged File Schema

```
ALL_SOURCES_drug_disease_merged.csv:
  - drug_name        : Standardized drug name (lowercase, trimmed)
  - disease_name     : Standardized disease name (lowercase, trimmed)
  - source           : Primary source identifier (AACT, ChEMBL, etc.)
  - internal_source  : Detailed source provenance
  - drug_id          : Source-specific drug identifier (may be empty)
  - disease_id       : Source-specific disease identifier (may be empty)
  - extra_metadata   : JSON-encoded source-specific fields
```

### Deduplicated File Schema

```
ALL_SOURCES_drug_disease_deduplicated.csv:
  - drug_name        : Standardized drug name (lowercase, trimmed)
  - drug_id          : Source-specific drug identifier (may be empty)
  - disease          : Standardized disease name (lowercase, trimmed)
  - source           : Primary source identifier
  - internal_source  : Detailed source provenance
```

### Extra Metadata Examples

```json
// AACT/ClinicalTrials.gov
{
  "nct_id": "NCT00391339",
  "intervention_type": "DRUG",
  "study_phase": "PHASE2/PHASE3",
  "overall_status": "UNKNOWN"
}

// ChEMBL
{
  "mesh_id": "D001172",
  "max_phase": 4,
  "efo_id": "EFO:0000685"
}

// CTD
{
  "direct_evidence": "therapeutic",
  "inference_score": 4.08
}
```

---

*This README provides comprehensive documentation for the drug-disease association data pipeline. For additional details on specific scripts or methodologies, please refer to inline code documentation.*
