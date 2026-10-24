
#!/usr/bin/env python3
"""AACT clinical-trial drug/condition co-occurrence extraction.

Uses a pinned local AACT ZIP snapshot and writes to data/processed/aact/
(overwriting the existing processed files). The resulting records are NOT
validated therapeutic indications. Run on an HPC compute node, not a login node.

Default: DRUG interventions. Review the intervention-type audit before
including BIOLOGICAL or other intervention types.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import shutil
import zipfile
from datetime import datetime, timezone
from pathlib import Path

import duckdb
import pandas as pd
import requests

ROOT = Path(
    "/data/ascher02/uqmmune1/ANNOVAR/drug_disease_data/drug-disease-mapping"
)
RAW = ROOT / "data/raw/aact_20260114"
OUT = ROOT / "data/processed/aact"

URL = (
    "https://ctti-aact.nyc3.digitaloceanspaces.com/"
    "pw7s52pt0ighmd1qcb5hbl8hikcs"
)

REQUIRED = {
    "interventions.txt",
    "conditions.txt",
    "studies.txt",
}

OPTIONAL = {
    "design_groups.txt",
    "design_group_interventions.txt",
}

SCHEMA = {
    "interventions.txt": {
        "id", "nct_id", "intervention_type", "name"
    },
    "conditions.txt": {
        "id", "nct_id", "name"
    },
    "studies.txt": {
        "nct_id", "study_type", "phase",
        "overall_status", "brief_title", "official_title"
    },
    "design_groups.txt": {
        "id", "nct_id", "group_type"
    },
    "design_group_interventions.txt": {
        "design_group_id", "intervention_id"
    },
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

EXTRAS = [
    "nct_id",
    "aact_intervention_id",
    "aact_condition_id",
    "intervention_name_original",
    "condition_name_original",
    "intervention_type",
    "study_type",
    "study_phase",
    "overall_status",
    "study_source_original",
    "arm_group_types",
    "intervention_role",
    "intervention_role_basis",
    "brief_title",
    "official_title",
    "clinicaltrials_url",
    "source_url",
]

COLUMNS = CORE + EXTRAS

ROLE_PATTERNS = [
    (
        "placebo_or_sham",
        r"\b(placebo|sham|vehicle|dummy)\b",
    ),
    (
        "no_intervention_or_standard_of_care",
        r"standard of care|standard care|usual care|"
        r"best supportive care|no intervention|no treatment|"
        r"observation only|watchful waiting",
    ),
    (
        "possibly_nontherapeutic",
        r"contrast agent|contrast medium|radiotracer|"
        r"\btracer\b|imaging agent|questionnaire|"
        r"\bsurvey\b|blood sampl|blood draw|\bbiopsy\b",
    ),
]

# Only unambiguous arm assignments enter this candidate view.
CANDIDATE_ROLES = {
    "experimental",
    "active_comparator",
}


def now():
    return datetime.now(timezone.utc).isoformat()


def section(title):
    print(
        "\n" + "=" * 92 +
        "\n" + title +
        "\n" + "=" * 92,
        flush=True,
    )


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(
            lambda: f.read(4 * 1024 * 1024), b""
        ):
            h.update(b)
    return h.hexdigest()


def save_csv(df, path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    tmp = path.with_name(path.name + ".tmp")

    try:
        df.to_csv(tmp, index=False)
        tmp.replace(path)
    finally:
        if tmp.exists():
            tmp.unlink()


def save_json(data, path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    tmp = path.with_name(path.name + ".tmp")

    try:
        tmp.write_text(
            json.dumps(
                data,
                indent=2,
                ensure_ascii=False,
                default=str,
            ) + "\n"
        )
        tmp.replace(path)
    finally:
        if tmp.exists():
            tmp.unlink()


def ensure_zip(
    zip_path,
    url,
    force_download=False,
    download_missing=False,
):
    if zip_path.is_file() and not force_download:
        with zipfile.ZipFile(zip_path) as z:
            bad = REQUIRED - {
                Path(x).name for x in z.namelist()
            }
            if bad:
                raise RuntimeError(
                    f"Cached ZIP is missing required files: "
                    f"{sorted(bad)}"
                )

        print(
            f"[REUSE] {zip_path} "
            f"({zip_path.stat().st_size:,} bytes)"
        )
        return

    if not (force_download or download_missing):
        raise FileNotFoundError(
            "AACT ZIP not found. Supply --download_missing "
            "to fetch it; the default NEVER silently "
            "downloads a different snapshot."
        )

    zip_path.parent.mkdir(parents=True, exist_ok=True)
    tmp = zip_path.with_name(
        zip_path.name + ".download"
    )

    try:
        with requests.get(
            url,
            stream=True,
            timeout=1200,
        ) as r:
            r.raise_for_status()

            with open(tmp, "wb") as f:
                for chunk in r.iter_content(1024 * 1024):
                    if chunk:
                        f.write(chunk)

        with zipfile.ZipFile(tmp) as z:
            bad = REQUIRED - {
                Path(x).name for x in z.namelist()
            }
            if bad:
                raise RuntimeError(
                    f"Downloaded ZIP missing: {sorted(bad)}"
                )

        tmp.replace(zip_path)

    finally:
        if tmp.exists():
            tmp.unlink()

    print(f"[DOWNLOADED] {zip_path}")


def inventory_and_extract(
    zip_path,
    extract_dir,
    force_extract=False,
):
    section(
        "1. DOWNLOADED AACT INVENTORY "
        "(INCLUDING EVERY ZIP MEMBER)"
    )

    with zipfile.ZipFile(zip_path) as z:
        members = [
            x for x in z.infolist()
            if not x.is_dir()
        ]

        print(f"ZIP file: {zip_path}")
        print(f"ZIP members: {len(members):,}")

        for member in members:
            filename = Path(member.filename).name

            label = (
                "USED"
                if filename in (REQUIRED | OPTIONAL)
                else "NOT USED"
            )

            print(
                f"  [{label:8}] {member.filename} "
                f"({member.file_size:,} bytes uncompressed)"
            )

        found = {}
        extract_dir.mkdir(
            parents=True,
            exist_ok=True,
        )

        for filename in sorted(REQUIRED | OPTIONAL):
            matches = [
                m for m in members
                if Path(m.filename).name == filename
            ]

            if not matches:
                if filename in REQUIRED:
                    raise RuntimeError(
                        f"ZIP missing required table: {filename}"
                    )

                print(f"  [OPTIONAL MISSING] {filename}")
                continue

            if len(matches) != 1:
                raise RuntimeError(
                    f"Ambiguous ZIP member basename: {filename}"
                )

            dest = extract_dir / filename

            if (
                force_extract
                or not dest.is_file()
                or not dest.stat().st_size
            ):
                tmp = dest.with_name(
                    dest.name + ".tmp"
                )

                try:
                    with (
                        z.open(matches[0]) as src,
                        open(tmp, "wb") as dst,
                    ):
                        shutil.copyfileobj(
                            src,
                            dst,
                            1024 * 1024,
                        )

                    tmp.replace(dest)

                    print(f"  [EXTRACTED] {filename}")

                finally:
                    if tmp.exists():
                        tmp.unlink()

            else:
                print(f"  [REUSE] {dest}")

            found[filename] = dest

    section(
        "2. PHYSICAL RAW FILES AND ACTUAL HEADERS"
    )

    physical = sorted(
        p for p in zip_path.parent.rglob("*")
        if p.is_file()
    )

    for p in physical:
        print(
            f"  {p.relative_to(zip_path.parent)} "
            f"({p.stat().st_size:,} bytes)"
        )

    print(
        f"Physical files in raw directory: "
        f"{len(physical):,}"
    )

    details = {}

    for filename, path in sorted(found.items()):
        with path.open(
            "r",
            encoding="utf-8-sig",
            errors="replace",
            newline="",
        ) as f:
            reader = csv.reader(
                f,
                delimiter="|",
                quotechar='"',
            )

            headers = next(reader)
            sample = next(reader, [])

        missing = SCHEMA[filename] - set(headers)

        if missing:
            raise RuntimeError(
                f"{filename}: missing columns {sorted(missing)}"
            )

        section(
            f"INPUT: {filename} | "
            f"{len(headers)} columns"
        )

        print(
            "Full header: " +
            " | ".join(headers)
        )

        print(
            "First raw row "
            "(independent example, not a joined record):"
        )

        for field, value in list(
            zip(headers, sample)
        )[:12]:
            val = str(value).replace("\n", " ")

            print(
                f"  {field:34} "
                f"{val[:115]}"
            )

        if len(headers) > 12:
            print(
                f"  ... ({len(headers)-12} remaining "
                f"columns listed in full header above)"
            )

        details[filename] = {
            "file": str(path),
            "size_bytes": path.stat().st_size,
            "sha256": sha256(path),
            "columns": headers,
        }

    return (
        found,
        details,
        len(physical),
        len(members),
    )


def classify_role(name, arm_string):
    """
    Conservative annotation.

    Missing arm metadata does not establish
    therapeutic intent.
    """
    arms = {
        x.strip().upper()
        for x in str(arm_string).split("|")
        if x.strip()
    }

    lname = str(name).lower()

    name_role = next(
        (
            label
            for label, pattern in ROLE_PATTERNS
            if re.search(pattern, lname)
        ),
        None,
    )

    if len(arms) > 1:
        return (
            "ambiguous_or_mixed_role",
            "multiple_arm_types",
        )

    if name_role:
        if (
            arms
            and name_role == "placebo_or_sham"
            and arms == {"EXPERIMENTAL"}
        ):
            return (
                "ambiguous_or_mixed_role",
                "name_arm_conflict",
            )

        return (
            name_role,
            "intervention_name_heuristic",
        )

    if arms == {"PLACEBO_COMPARATOR"}:
        return (
            "placebo_or_sham",
            "arm_group_type",
        )

    if arms == {"SHAM_COMPARATOR"}:
        return (
            "placebo_or_sham",
            "arm_group_type",
        )

    if arms == {"NO_INTERVENTION"}:
        return (
            "no_intervention_or_standard_of_care",
            "arm_group_type",
        )

    if arms == {"ACTIVE_COMPARATOR"}:
        return (
            "active_comparator",
            "arm_group_type",
        )

    if arms == {"EXPERIMENTAL"}:
        return (
            "experimental",
            "arm_group_type",
        )

    if arms:
        return (
            "uncertain_other_arm_type",
            "arm_group_type",
        )

    return (
        "unknown_role",
        "no_arm_evidence",
    )


def sql_csv_view(con, view, path):
    p = str(path.resolve()).replace("'", "''")

    con.execute(
        f"""
        CREATE VIEW {view} AS
        SELECT *
        FROM read_csv_auto(
            '{p}',
            delim='|',
            header=true,
            quote='"',
            escape='"',
            all_varchar=true,
            ignore_errors=false,
            max_line_size=10000000
        )
        """
    )


def check_unique(con, table, col):
    n = con.execute(
        f"""
        SELECT
            COUNT(*) - COUNT(DISTINCT {col})
        FROM {table}
        """
    ).fetchone()[0]

    if n:
        raise RuntimeError(
            f"{table}.{col} contains {n:,} "
            f"duplicated or missing keys"
        )


def run_processing(
    found,
    args,
    input_details,
    physical_count,
    zip_members,
):
    section(
        "3. JOIN DESIGN / NAMING / "
        "SCIENTIFIC INTERPRETATION"
    )

    print(
        "Association: interventions INNER JOIN "
        "conditions USING (nct_id)"
    )

    print(
        "Metadata:    LEFT JOIN studies "
        "USING (nct_id)"
    )

    print(
        "Arm role:    LEFT JOIN aggregated "
        "design_groups/design_group_interventions "
        "by intervention ID"
    )

    print(
        "Naming:      drug_name = "
        "trim(interventions.name); "
        "condition = trim(conditions.name)"
    )

    print(
        "Identifiers: none inferred from free-text names; "
        "AACT IDs and NCT are retained"
    )

    print(
        "Study source field is saved as "
        "study_source_original; core source remains AACT"
    )

    print(
        "MEANING: trial intervention–condition "
        "CO-OCCURRENCE, NOT efficacy or approval"
    )

    con = duckdb.connect(database=":memory:")

    try:
        for name, table in [
            ("interventions.txt", "interventions"),
            ("conditions.txt", "conditions"),
            ("studies.txt", "studies"),
        ]:
            sql_csv_view(
                con,
                table,
                found[name],
            )

        check_unique(
            con, "interventions", "id"
        )
        check_unique(
            con, "conditions", "id"
        )
        check_unique(
            con, "studies", "nct_id"
        )

        design = all(
            k in found for k in OPTIONAL
        )

        if design:
            sql_csv_view(
                con,
                "design_groups",
                found["design_groups.txt"],
            )

            sql_csv_view(
                con,
                "design_group_interventions",
                found["design_group_interventions.txt"],
            )

            con.execute(
                """
                CREATE TEMP VIEW arm_types AS
                SELECT
                    dgi.intervention_id,
                    string_agg(
                        DISTINCT upper(trim(dg.group_type)),
                        '|'
                        ORDER BY upper(trim(dg.group_type))
                    ) AS arm_group_types
                FROM design_group_interventions dgi
                JOIN design_groups dg
                    ON dgi.design_group_id = dg.id
                WHERE
                    dg.group_type IS NOT NULL
                    AND trim(dg.group_type) <> ''
                GROUP BY dgi.intervention_id
                """
            )

        else:
            con.execute(
                """
                CREATE TEMP VIEW arm_types AS
                SELECT
                    CAST(NULL AS VARCHAR) intervention_id,
                    CAST(NULL AS VARCHAR) arm_group_types
                WHERE false
                """
            )

        counts = {
            table: int(
                con.execute(
                    f"SELECT COUNT(*) FROM {table}"
                ).fetchone()[0]
            )
            for table in [
                "interventions",
                "conditions",
                "studies",
            ]
        }

        section(
            "4. SOURCE TABLE COUNTS "
            "AND INTERVENTION TYPES"
        )

        for k, v in counts.items():
            print(
                f"{k:16} {v:>12,} rows"
            )

        type_audit = con.execute(
            """
            SELECT
                upper(coalesce(intervention_type, ''))
                    AS intervention_type,
                count(*) AS intervention_rows,
                count(distinct nct_id) AS trials,
                count(distinct lower(trim(name)))
                    AS distinct_names
            FROM interventions
            GROUP BY 1
            ORDER BY 2 DESC
            """
        ).df()

        save_csv(
            type_audit,
            args.intervention_audit_out,
        )

        print(
            type_audit.to_string(index=False)
        )

        if args.keep_all_intervention_types:
            type_where = "TRUE"
            chosen_types = ["ALL"]

        else:
            chosen_types = sorted(
                set(
                    x.strip().upper()
                    for x in
                    args.include_intervention_types.split(",")
                    if x.strip()
                )
            )

            if (
                not chosen_types
                or any(
                    not re.fullmatch(r"[A-Z_ ]+", x)
                    for x in chosen_types
                )
            ):
                raise ValueError(
                    "Invalid --include_intervention_types value"
                )

            type_where = (
                "upper(i.intervention_type) IN (" +
                ", ".join(
                    "'" + x + "'"
                    for x in chosen_types
                ) +
                ")"
            )

        print(
            "Included intervention types:",
            ", ".join(chosen_types),
        )
        print(
            "Design tables used:", design
        )

        safe_url = args.url.replace("'", "''")

        query = f"""
            SELECT
                trim(i.name) AS drug_name,
                '' AS drug_identifier,
                '' AS drug_identifier_type,

                trim(c.name)
                    AS disease_or_condition_name,
                '' AS disease_or_condition_identifier,
                '' AS disease_or_condition_identifier_type,

                CASE
                    WHEN upper(i.intervention_type) = 'DRUG'
                    THEN 'clinical_trial_drug_condition'
                    ELSE 'clinical_trial_intervention_condition'
                END AS relationship_type,

                CASE
                    WHEN upper(i.intervention_type) = 'DRUG'
                    THEN
                    'aact_clinical_trial_drug_condition_cooccurrence'
                    ELSE
                    'aact_clinical_trial_intervention_condition_cooccurrence'
                END AS evidence_type,

                'AACT' AS source,
                'ClinicalTrials.gov' AS internal_source,

                i.nct_id,
                i.id AS aact_intervention_id,
                c.id AS aact_condition_id,

                i.name AS intervention_name_original,
                c.name AS condition_name_original,

                coalesce(i.intervention_type,'')
                    AS intervention_type,

                coalesce(s.study_type,'')
                    AS study_type,

                coalesce(s.phase,'')
                    AS study_phase,

                coalesce(s.overall_status,'')
                    AS overall_status,

                coalesce(s.source,'')
                    AS study_source_original,

                coalesce(a.arm_group_types,'')
                    AS arm_group_types,

                coalesce(s.brief_title,'')
                    AS brief_title,

                coalesce(s.official_title,'')
                    AS official_title,

                'https://clinicaltrials.gov/study/'
                    || i.nct_id
                    AS clinicaltrials_url,

                '{safe_url}'
                    AS source_url

            FROM interventions i

            JOIN conditions c
                USING (nct_id)

            LEFT JOIN studies s
                USING (nct_id)

            LEFT JOIN arm_types a
                ON a.intervention_id = i.id

            WHERE
                {type_where}

                AND i.name IS NOT NULL
                AND trim(i.name) <> ''

                AND c.name IS NOT NULL
                AND trim(c.name) <> ''

                AND i.nct_id IS NOT NULL
                AND trim(i.nct_id) <> ''
        """

        expected = con.execute(
            f"""
            WITH selected AS (
                SELECT
                    nct_id,
                    count(*) AS n_interventions
                FROM interventions i
                WHERE
                    {type_where}
                    AND i.name IS NOT NULL
                    AND trim(i.name) <> ''
                    AND i.nct_id IS NOT NULL
                    AND trim(i.nct_id) <> ''
                GROUP BY nct_id
            ),

            conditions_ok AS (
                SELECT
                    nct_id,
                    count(*) AS n_conditions
                FROM conditions
                WHERE
                    name IS NOT NULL
                    AND trim(name) <> ''
                GROUP BY nct_id
            )

            SELECT
                s.nct_id,
                s.n_interventions,
                c.n_conditions,
                s.n_interventions * c.n_conditions
                    AS expected_cartesian_combinations

            FROM selected s
            JOIN conditions_ok c
                USING (nct_id)
            """
        ).df()

        section(
            "5. JOIN EXECUTION AND "
            "PER-TRIAL QUALITY CHECK"
        )

        df = con.execute(query).df()

        for col in COLUMNS:
            if col not in df:
                df[col] = ""

            df[col] = (
                df[col]
                .fillna("")
                .astype(str)
            )

        df["drug_name"] = (
            df["drug_name"].str.strip()
        )

        df["disease_or_condition_name"] = (
            df["disease_or_condition_name"].str.strip()
        )

        roles = [
            classify_role(n, a)
            for n, a in zip(
                df["intervention_name_original"],
                df["arm_group_types"],
            )
        ]

        df["intervention_role"] = [
            x[0] for x in roles
        ]

        df["intervention_role_basis"] = [
            x[1] for x in roles
        ]

        df = df[COLUMNS]

        observed = (
            df.groupby("nct_id")
            .size()
            .rename("observed_combinations")
            .reset_index()
        )

        combos = expected.merge(
            observed,
            on="nct_id",
            how="outer",
        ).fillna(0)

        combos["matches"] = (
            combos["expected_cartesian_combinations"]
            == combos["observed_combinations"]
        )

        mismatches = combos.loc[
            ~combos["matches"]
        ]

        print(
            "Expected combinations:",
            f"{int(combos['expected_cartesian_combinations'].sum()):,}",
        )

        print(
            "Observed combinations:",
            f"{len(df):,}",
        )

        print(
            "Mismatching trials:",
            f"{len(mismatches):,}",
        )

        if len(mismatches):
            print(
                mismatches.head(10).to_string(index=False)
            )

            raise RuntimeError(
                "Join multiplicity did not equal filtered "
                "Cartesian count. Outputs not overwritten."
            )

        print(
            "Role distribution:\n" +
            df["intervention_role"]
            .value_counts()
            .to_string()
        )

        section(
            "6. WRITING EXISTING PROCESSED "
            "DIRECTORY (OVERWRITES)"
        )

        save_csv(df, args.out)

        dedup_cols = [
            "nct_id",
            "aact_intervention_id",
            "aact_condition_id",
        ]

        duplicates = df.duplicated(
            subset=dedup_cols
        )

        dedup = df.loc[
            ~duplicates
        ].copy()

        save_csv(
            dedup,
            args.dedup_out,
        )

        # Conservative candidate therapeutic co-occurrence view.
        # Not evidence of efficacy or regulatory approval.
        candidate = df[
            (
                df["study_type"].str.upper()
                == "INTERVENTIONAL"
            )
            &
            (
                df["intervention_role"].isin(
                    CANDIDATE_ROLES
                )
            )
        ].copy()

        save_csv(
            candidate,
            args.therapeutic_view_out,
        )

        combo_path = (
            Path(args.out)
            .with_name(
                "aact_trial_combination_counts.csv"
            )
        )

        save_csv(
            combos,
            combo_path,
        )

        section(
            "7. OUTPUT SCHEMA, COUNTS, "
            "AND ACTUAL SAMPLE"
        )

        print(
            "Full output columns:", len(COLUMNS)
        )

        print(
            " | ".join(COLUMNS)
        )

        print(
            f"Full co-occurrences:        {len(df):,}"
        )

        print(
            f"Deduped intervention events:{len(dedup):,}"
        )

        print(
            f"Candidate co-occurrences:  {len(candidate):,} "
            "(NOT proven indications)"
        )

        print(
            f"Unique trials:             "
            f"{df['nct_id'].nunique():,}"
        )

        print("Example processed row:")

        if len(df):
            for key, value in df.iloc[0].items():
                print(
                    f"  {key:36} "
                    f"{str(value)[:100]}"
                )

        print("Files overwritten:")

        for f in [
            args.out,
            args.dedup_out,
            args.therapeutic_view_out,
            args.intervention_audit_out,
            combo_path,
        ]:
            print(" ", f)

        return {
            "source": "AACT",
            "interpretation": (
                "trial intervention-condition "
                "co-occurrence only"
            ),
            "started_at_utc": args.started,
            "finished_at_utc": now(),
            "source_url": args.url,
            "zip_file": str(
                Path(args.workdir) /
                "aact_export.zip"
            ),
            "zip_members": zip_members,
            "physical_raw_files": physical_count,
            "source_inputs": input_details,
            "input_row_counts": counts,
            "design_tables_used": design,
            "included_intervention_types": chosen_types,
            "join": (
                "interventions INNER JOIN conditions "
                "ON nct_id; LEFT JOIN studies ON nct_id; "
                "LEFT JOIN arm_types ON intervention id"
            ),
            "name_mapping": {
                "drug_name": (
                    "trim(interventions.name)"
                ),
                "disease_or_condition_name": (
                    "trim(conditions.name)"
                ),
            },
            "identity": (
                "nct_id + interventions.id + conditions.id"
            ),
            "drug_identifier_strategy": (
                "not inferred (empty)"
            ),
            "output_columns": COLUMNS,
            "output_file": str(args.out),
            "output_rows": len(df),
            "dedup_file": str(args.dedup_out),
            "dedup_keys": dedup_cols,
            "dedup_rows": len(dedup),
            "candidate_view_file": str(
                args.therapeutic_view_out
            ),
            "candidate_definition": (
                "INTERVENTIONAL study and unambiguous "
                "EXPERIMENTAL or ACTIVE_COMPARATOR arm; "
                "never efficacy proof"
            ),
            "candidate_rows": len(candidate),
            "trial_combinations_file": str(combo_path),
            "trial_combinations_verified": True,
            "intervention_audit_file": str(
                args.intervention_audit_out
            ),
            "intervention_role_counts": (
                df["intervention_role"]
                .value_counts()
                .to_dict()
            ),
            "intervention_type_counts": (
                df["intervention_type"]
                .value_counts()
                .to_dict()
            ),
            "study_type_counts": (
                df["study_type"]
                .value_counts()
                .to_dict()
            ),
            "unique_trials": (
                df["nct_id"].nunique()
            ),
            "output_sha256": sha256(args.out),
            "dedup_sha256": sha256(args.dedup_out),
        }

    finally:
        con.close()


def main():
    parser = argparse.ArgumentParser(
        description=__doc__
    )

    parser.add_argument(
        "--url",
        default=URL,
    )

    parser.add_argument(
        "--workdir",
        default=str(RAW),
    )

    parser.add_argument(
        "--out",
        default=str(
            OUT / "aact_drug_disease.csv"
        ),
    )

    parser.add_argument(
        "--dedup_out",
        default=str(
            OUT / "aact_drug_disease_deduplicated.csv"
        ),
    )

    parser.add_argument(
        "--therapeutic_view_out",
        default=str(
            OUT / "aact_drug_disease_therapeutic_view.csv"
        ),
    )

    parser.add_argument(
        "--intervention_audit_out",
        default=str(
            OUT / "aact_intervention_type_audit.csv"
        ),
    )

    parser.add_argument(
        "--summary",
        default=str(
            OUT / "aact_summary.json"
        ),
    )

    parser.add_argument(
        "--include_intervention_types",
        default="DRUG",
    )

    parser.add_argument(
        "--keep_all_intervention_types",
        action="store_true",
    )

    parser.add_argument(
        "--force_download",
        action="store_true",
    )

    parser.add_argument(
        "--download_missing",
        action="store_true",
    )

    parser.add_argument(
        "--force_extract",
        action="store_true",
    )

    parser.add_argument(
        "--head_n",
        type=int,
        default=5,
        help="Accepted for compatibility",
    )

    args = parser.parse_args()
    args.started = now()

    zip_path = (
        Path(args.workdir) /
        "aact_export.zip"
    )

    try:
        ensure_zip(
            zip_path,
            args.url,
            args.force_download,
            args.download_missing,
        )

        files, details, nphysical, nmembers = (
            inventory_and_extract(
                zip_path,
                Path(args.workdir) / "extract",
                args.force_extract,
            )
        )

        result = run_processing(
            files,
            args,
            details,
            nphysical,
            nmembers,
        )

        result["status"] = "success"

        save_json(
            result,
            args.summary,
        )

        print(
            f"\nSUCCESS — summary overwritten: "
            f"{args.summary}"
        )

    except Exception as e:
        failure = (
            Path(args.summary)
            .with_name("aact_failed_run.json")
        )

        save_json(
            {
                "status": "failed",
                "error": repr(e),
                "time": now(),
            },
            failure,
        )

        print(
            "\nFAILED — original summary untouched; "
            f"diagnostic: {failure}"
        )

        raise


if __name__ == "__main__":
    main()
