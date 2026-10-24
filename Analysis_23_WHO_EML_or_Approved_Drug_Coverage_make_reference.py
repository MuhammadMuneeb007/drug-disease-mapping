
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Analysis 23 reference preparation: FDA Orange Book ingredient audit.

Produces data/reference/approved_drugs.csv for Analysis 23.
Default reference includes ONLY active ingredients appearing in at least one
single-ingredient RX/OTC Orange Book product (Type=RX/OTC).

A second complete provenance table records ingredients found in combination
products. Use --reference-scope all-components only if you deliberately want
components found solely in combinations in the reference denominator.

Source names are retained verbatim apart from whitespace/case normalisation.
Salt/prodrug/isotope changes are NEVER automatically applied to the reference.
All outputs are derived; the frozen 11-source merger remains unchanged.

Runs with Python standard library only.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import os
import re
import shutil
import tempfile
import unicodedata
import zipfile

from collections import Counter, defaultdict
from contextlib import contextmanager
from datetime import datetime, timezone
from pathlib import Path


# ============================================================
# CONFIGURATION
# ============================================================

ROOT = Path(__file__).resolve().parent

EXPECTED_HEADER = (
    "Ingredient",
    "Type",
    "Appl_Type",
    "Appl_No",
    "Product_No",
)

ACCEPTED_TYPES = {
    "RX",
    "OTC",
    "DISCN",
}

ACCEPTED_APPLICATIONS = {
    "A",
    "N",
}

# These terms are only used to generate potential
# name equivalences for manual review.
#
# They are NEVER automatically removed from approved_drugs.csv.

REVIEW_SUFFIXES = {
    "HYDROCHLORIDE",
    "HCL",
    "SODIUM",
    "POTASSIUM",
    "CALCIUM",
    "MAGNESIUM",
    "PHOSPHATE",
    "SULFATE",
    "SULPHATE",
    "NITRATE",
    "CITRATE",
    "MESYLATE",
    "MALEATE",
    "TARTRATE",
    "ACETATE",
    "SUCCINATE",
    "FUMARATE",
    "BESYLATE",
    "BESILATE",
    "BROMIDE",
    "CHLORIDE",
    "DIHYDRATE",
    "MONOHYDRATE",
    "TRIHYDRATE",
    "HYDRATE",
    "DISODIUM",
    "MONOSODIUM",
    "HEMIFUMARATE",
    "DIMALEATE",
    "LACTATE",
    "PROPIONATE",
    "TOSYLATE",
    "BITARTRATE",
}

# These standalone fragments are suspicious in an
# ingredient field and require inspection.

SUSPICIOUS_STANDALONE = {
    "BITARTRATE",
    "CARBONATE",
    "DIMETHYL",
    "DIBASIC",
    "TRIBASIC",
    "MONOBASIC",
    "ACETATE",
    "SULFATE",
    "PHOSPHATE",
    "SODIUM",
    "POTASSIUM",
    "FERRIC",
    "FERROUS",
    "CUPRIC",
    "CHROMIC",
}

WS = re.compile(r"\s+")

AMBIGUOUS_SEPARATORS = re.compile(
    r"\s+AND\s+|\s*\+\s*|/",
    re.I,
)


# ============================================================
# NORMALISATION
# ============================================================

def norm(value: str | None) -> str:
    """
    Preserve the original drug-substance identity.

    Do not remove salts, isotopes, prodrug names
    or other chemically meaningful name components.
    """
    return WS.sub(
        " ",
        unicodedata.normalize(
            "NFKC",
            value or "",
        ).strip(),
    ).upper()


# ============================================================
# CHECKSUM
# ============================================================

def digest(path: Path) -> str:
    sha = hashlib.sha256()

    with path.open("rb") as handle:
        for part in iter(
            lambda: handle.read(8 * 1024 * 1024),
            b"",
        ):
            sha.update(part)

    return sha.hexdigest()


# ============================================================
# LOCATE FDA PRODUCTS FILE
# ============================================================

def find_products_file(
    custom: Path | None,
) -> Path:

    if custom is not None:
        path = custom.expanduser().resolve()

        if not path.is_file():
            raise FileNotFoundError(
                f"--products not found: {path}"
            )

        return path

    choices = [
        ROOT / "products.txt",
        ROOT / "Products.txt",

        ROOT / "data/raw/orange_book/products.txt",
        ROOT / "data/raw/orange_book/Products.txt",

        ROOT / "data/reference/products.txt",

        ROOT / "eobzip.zip",
        ROOT / "data/raw/orange_book/eobzip.zip",
    ]

    for path in choices:
        if path.is_file() and path.stat().st_size:
            return path.resolve()

    raise FileNotFoundError(
        "FDA Orange Book products.txt (or zip) not found. "
        "Use --products /path/to/products.txt\nChecked:\n"
        + "\n".join(str(x) for x in choices)
    )


# ============================================================
# OPEN FDA FILE OR ZIP ARCHIVE
# ============================================================

@contextmanager
def open_products(path: Path):
    """
    Open a tilde-delimited FDA Products.txt file.

    Also supports the original ZIP archive directly.
    """

    if path.suffix.lower() == ".zip":

        with zipfile.ZipFile(path) as archive:

            candidates = [
                x
                for x in archive.namelist()
                if Path(x).name.lower() == "products.txt"
            ]

            if len(candidates) != 1:
                raise RuntimeError(
                    "ZIP must contain exactly one products.txt; "
                    f"found {candidates}"
                )

            with archive.open(candidates[0]) as binary:
                with io.TextIOWrapper(
                    binary,
                    encoding="utf-8-sig",
                    errors="strict",
                    newline="",
                ) as handle:
                    yield handle

    else:

        with path.open(
            "r",
            encoding="utf-8-sig",
            errors="strict",
            newline="",
        ) as handle:
            yield handle


# ============================================================
# TERMINAL OUTPUT
# ============================================================

def print_tsv(
    title: str,
    rows: list[dict],
    limit: int = 25,
) -> None:

    print("\n" + "=" * 100)
    print(f"{title}: {len(rows):,} rows")
    print("=" * 100)

    if not rows:
        print("[NO ROWS]")
        return

    columns = list(rows[0])

    print("\t".join(columns))

    for row in rows[:limit]:
        values = []

        for column in columns:

            value = str(
                row.get(column, "")
            )

            value = (
                value
                .replace("\t", " ")
                .replace("\n", " ")
            )

            values.append(value)

        print("\t".join(values))

    if len(rows) > limit:
        print(
            f"[SHOWING {limit:,} OF {len(rows):,}; "
            "saved CSV has all rows]"
        )


# ============================================================
# SAFE CSV WRITING
# ============================================================

def write_csv_atomic(
    path: Path,
    rows: list[dict],
    fieldnames: list[str],
) -> None:
    """
    Avoid partially written supplementary CSV files.
    """

    path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    with tempfile.NamedTemporaryFile(
        "w",
        newline="",
        encoding="utf-8",
        dir=path.parent,
        prefix=".ob_",
        suffix=".tmp",
        delete=False,
    ) as handle:

        tmp = Path(handle.name)

        try:

            writer = csv.DictWriter(
                handle,
                fieldnames=fieldnames,
                extrasaction="ignore",
            )

            writer.writeheader()
            writer.writerows(rows)

        except Exception:
            tmp.unlink(missing_ok=True)
            raise

    os.replace(tmp, path)

    print(
        f"[SAVED] {path} ({len(rows):,} rows)"
    )


# ============================================================
# MAIN
# ============================================================

def main() -> None:

    parser = argparse.ArgumentParser(
        description=__doc__,
    )

    parser.add_argument(
        "--products",
        type=Path,
        default=None,
        help=(
            "FDA products.txt, Products.txt "
            "or eobzip.zip"
        ),
    )

    parser.add_argument(
        "--reference-dir",
        type=Path,
        default=ROOT / "data/reference",
    )

    parser.add_argument(
        "--reference-scope",
        choices=[
            "mono-only",
            "all-components",
        ],
        default="mono-only",
        help=(
            "Default is scientifically stricter; "
            "combination-only does not establish "
            "single-agent approval"
        ),
    )

    parser.add_argument(
        "--include-discontinued",
        action="store_true",
        help=(
            "Include Type=DISCN historical/"
            "discontinued entries"
        ),
    )

    parser.add_argument(
        "--nda-only",
        action="store_true",
        help=(
            "Use NDA (Appl_Type=N), "
            "rather than NDA+ANDA"
        ),
    )

    parser.add_argument(
        "--top",
        type=int,
        default=30,
    )

    parser.add_argument(
        "--dry-run",
        action="store_true",
    )

    args = parser.parse_args()

    if not 1 <= args.top <= 1000:
        parser.error(
            "--top must be 1..1000"
        )

    input_path = find_products_file(
        args.products
    )

    outdir = (
        args.reference_dir
        .expanduser()
        .resolve()
    )

    outfile = outdir / "approved_drugs.csv"

    print("\nFDA ORANGE BOOK REFERENCE BUILDER (AUDITED)")
    print("[SOURCE]", input_path)

    print(
        "[REFERENCE SCOPE]",
        args.reference_scope,
    )

    print(
        "[TYPES]",
        (
            "RX / OTC / DISCN"
            if args.include_discontinued
            else "RX / OTC"
        ),
    )

    print(
        "[APPLICATIONS]",
        (
            "NDA only"
            if args.nda_only
            else "NDA + ANDA"
        ),
    )

    print(
        "[POLICY] FDA ingredient substance names remain "
        "unchanged; no automatic salt stripping"
    )

    print(
        "[POLICY] Split combinations on the FDA's "
        "semicolon delimiter only"
    )

    # ========================================================
    # INGREDIENT AND COMBINATION INVENTORIES
    # ========================================================

    ingredients = defaultdict(
        lambda: {
            "monotherapy_product_rows": 0,
            "combination_product_rows": 0,
            "nda_product_rows": 0,
            "anda_product_rows": 0,
            "products": set(),
            "trade_names": set(),
        }
    )

    combos = defaultdict(
        lambda: {
            "product_rows": 0,
            "products": set(),
            "trade_names": set(),
        }
    )

    issues = Counter()

    review = []

    rows_read = 0
    rows_eligible = 0
    rows_mono = 0
    rows_combo = 0

    excluded_types = Counter()
    other_apps = Counter()

    eligible_types = {
        "RX",
        "OTC",
    }

    if args.include_discontinued:
        eligible_types.add("DISCN")

    # ========================================================
    # READ FDA ORANGE BOOK PRODUCTS.TXT
    # ========================================================

    with open_products(input_path) as handle:

        reader = csv.DictReader(
            handle,
            delimiter="~",
            restkey="__extra__",
        )

        headers = reader.fieldnames or []

        lookup = {
            norm(h).replace(" ", "_"): h
            for h in headers
            if h
        }

        needed = {
            h.upper(): h
            for h in EXPECTED_HEADER
        }

        missing = [
            h
            for h in needed
            if h not in lookup
        ]

        if missing:
            raise RuntimeError(
                f"Missing FDA Orange Book fields: {missing}; "
                f"observed: {headers}\n"
                "Verify this is products.txt, not patent.txt."
            )

        cols = {
            h: lookup[h]
            for h in needed
        }

        trade_col = lookup.get("TRADE_NAME")

        for lineno, row in enumerate(
            reader,
            start=2,
        ):

            rows_read += 1

            extra = row.get("__extra__")

            if extra and any(
                str(x).strip()
                for x in extra
            ):
                raise RuntimeError(
                    "Extra unexpected tilde-delimited "
                    f"cells on line {lineno}"
                )

            # ------------------------------------------------
            # FDA PRODUCT CATEGORY
            # ------------------------------------------------

            status = norm(
                row.get(cols["TYPE"])
            )

            if status not in ACCEPTED_TYPES:
                raise RuntimeError(
                    "Unknown Orange Book Type "
                    f"on line {lineno}: {status!r}"
                )

            if status not in eligible_types:
                excluded_types[status] += 1
                continue

            # ------------------------------------------------
            # APPLICATION TYPE
            # ------------------------------------------------

            app = norm(
                row.get(cols["APPL_TYPE"])
            )

            if app not in ACCEPTED_APPLICATIONS:
                other_apps[
                    app or "BLANK"
                ] += 1
                continue

            if args.nda_only and app != "N":
                continue

            # ------------------------------------------------
            # INGREDIENT AND PRODUCT INFORMATION
            # ------------------------------------------------

            raw = norm(
                row.get(cols["INGREDIENT"])
            )

            appl_no = norm(
                row.get(cols["APPL_NO"])
            )

            prod_no = norm(
                row.get(cols["PRODUCT_NO"])
            )

            trade = (
                norm(row.get(trade_col))
                if trade_col
                else ""
            )

            if not raw:

                issues["EMPTY_INGREDIENT"] += 1

                review.append({
                    "line": lineno,
                    "issue": "EMPTY_INGREDIENT",
                    "ingredient": raw,
                    "application": app + appl_no,
                    "product_no": prod_no,
                    "type": status,
                })

                continue

            # FDA documents semicolons as the separator
            # between multiple active ingredients.

            parts = [
                norm(part)
                for part in raw.split(";")
            ]

            # ------------------------------------------------
            # AMBIGUOUS INGREDIENT DELIMITERS
            # ------------------------------------------------

            if (
                any(not p for p in parts)
                or any(
                    AMBIGUOUS_SEPARATORS.search(p)
                    for p in parts
                )
            ):

                issues["AMBIGUOUS_DELIMITER"] += 1

                review.append({
                    "line": lineno,
                    "issue": "AMBIGUOUS_DELIMITER",
                    "ingredient": raw,
                    "application": app + appl_no,
                    "product_no": prod_no,
                    "type": status,
                })

                continue

            if len(set(parts)) != len(parts):

                issues["REPEATED_COMPONENT"] += 1

                review.append({
                    "line": lineno,
                    "issue": "REPEATED_COMPONENT",
                    "ingredient": raw,
                    "application": app + appl_no,
                    "product_no": prod_no,
                    "type": status,
                })

                continue

            # ------------------------------------------------
            # SUSPICIOUS STANDALONE FRAGMENTS
            # ------------------------------------------------

            if any(
                p in SUSPICIOUS_STANDALONE
                or len(p) < 3
                for p in parts
            ):

                issues[
                    "SUSPICIOUS_STANDALONE_INGREDIENT"
                ] += 1

                review.append({
                    "line": lineno,
                    "issue": (
                        "SUSPICIOUS_STANDALONE_INGREDIENT"
                    ),
                    "ingredient": raw,
                    "application": app + appl_no,
                    "product_no": prod_no,
                    "type": status,
                })

                continue

            product_id = (
                f"{app}:{appl_no}:{prod_no}"
            )

            rows_eligible += 1

            # ------------------------------------------------
            # CLASSIFY SINGLE VS COMBINATION PRODUCT
            # ------------------------------------------------

            if len(parts) == 1:

                rows_mono += 1

            else:

                rows_combo += 1

                signature = "; ".join(
                    sorted(parts)
                )

                group = combos[signature]

                group["product_rows"] += 1
                group["products"].add(product_id)

                if trade:
                    group["trade_names"].add(trade)

            # ------------------------------------------------
            # AGGREGATE EACH COMPONENT'S PROVENANCE
            # ------------------------------------------------

            for part in parts:

                item = ingredients[part]

                if len(parts) == 1:
                    item[
                        "monotherapy_product_rows"
                    ] += 1
                else:
                    item[
                        "combination_product_rows"
                    ] += 1

                if app == "N":
                    item["nda_product_rows"] += 1
                else:
                    item["anda_product_rows"] += 1

                item["products"].add(product_id)

                if trade:
                    item["trade_names"].add(trade)

    # ========================================================
    # CONSERVATION CHECKS
    # ========================================================

    if not rows_eligible or not ingredients:
        raise RuntimeError(
            "No usable product rows; check Type/Appl_Type "
            "and the input data"
        )

    if rows_mono + rows_combo != rows_eligible:
        raise RuntimeError(
            "Failed to reconcile eligible product counts"
        )

    # ========================================================
    # BUILD REFERENCE DRUG SETS
    # ========================================================

    mono = {
        name
        for name, obj in ingredients.items()
        if obj["monotherapy_product_rows"] > 0
    }

    all_components = set(
        ingredients
    )

    combo_only = (
        all_components - mono
    )

    if args.reference_scope == "mono-only":
        selected = mono
    else:
        selected = all_components

    approved = [
        {
            "active_ingredient": name
        }
        for name in sorted(selected)
    ]

    # No semicolon-separated combination can survive
    # into approved_drugs.csv.

    if any(
        ";" in r["active_ingredient"]
        for r in approved
    ):
        raise RuntimeError(
            "Unexpected combination delimiter "
            "in final reference"
        )

    if len(approved) != len({
        r["active_ingredient"]
        for r in approved
    }):
        raise RuntimeError(
            "Duplicate ingredient name in reference"
        )

    # ========================================================
    # INGREDIENT-LEVEL PROVENANCE
    # ========================================================

    provenance = []

    base_review = []

    for name, obj in sorted(
        ingredients.items()
    ):

        provenance.append({
            "active_ingredient": name,

            "listed_in_single_ingredient_product":
                obj["monotherapy_product_rows"] > 0,

            "single_ingredient_product_rows":
                obj["monotherapy_product_rows"],

            "combination_product_rows":
                obj["combination_product_rows"],

            "nda_product_rows":
                obj["nda_product_rows"],

            "anda_product_rows":
                obj["anda_product_rows"],

            "distinct_application_product_ids":
                len(obj["products"]),

            "included_in_approved_drugs_csv":
                name in selected,

            "trade_name_examples":
                " | ".join(
                    sorted(obj["trade_names"])[:5]
                ),
        })

        # Salt/formulation stripping is only suggested
        # for human review, never automatically applied.

        words = name.split()

        if (
            len(words) > 1
            and words[-1] in REVIEW_SUFFIXES
        ):

            candidate = " ".join(
                words[:-1]
            )

            base_review.append({
                "fda_ingredient": name,

                "proposed_base_for_manual_review":
                    candidate,

                "warning": (
                    "NOT_AUTOMATIC_EQUIVALENCE; "
                    "salt/prodrug identity may matter"
                ),

                "candidate_is_another_FDA_ingredient":
                    candidate in ingredients,

                "candidate_is_in_main_reference":
                    candidate in selected,
            })

    # ========================================================
    # COMBINATION PRODUCT INVENTORY
    # ========================================================

    combo_rows = []

    for name, obj in combos.items():

        combo_rows.append({
            "full_combination_ingredient_field": name,

            "n_components":
                len(name.split("; ")),

            "product_rows":
                obj["product_rows"],

            "distinct_application_product_ids":
                len(obj["products"]),

            "trade_name_examples":
                " | ".join(
                    sorted(obj["trade_names"])[:5]
                ),
        })

    combo_rows.sort(
        key=lambda r: (
            -r["product_rows"],
            r["full_combination_ingredient_field"],
        )
    )

    combo_only_rows = [
        row
        for row in provenance
        if row["active_ingredient"] in combo_only
    ]

    # ========================================================
    # PRINT SOURCE COUNTS
    # ========================================================

    print(
        f"[INPUT PRODUCT ROWS] {rows_read:,}"
    )

    print(
        f"[ELIGIBLE PRODUCT ROWS] {rows_eligible:,}"
    )

    print(
        f"[SINGLE-INGREDIENT PRODUCT ROWS] {rows_mono:,}"
    )

    print(
        f"[COMBINATION PRODUCT ROWS] {rows_combo:,}"
    )

    print(
        f"[MONOTHERAPY INGREDIENT NAMES] {len(mono):,}"
    )

    print(
        "[ALL INGREDIENT COMPONENT NAMES] "
        f"{len(all_components):,}"
    )

    print(
        "[COMBINATION-ONLY INGREDIENT NAMES] "
        f"{len(combo_only):,}"
    )

    print(
        f"[MAIN REFERENCE NAMES] {len(approved):,}"
    )

    print(
        f"[EXCLUDED TYPES] {dict(excluded_types)}"
    )

    print(
        f"[SKIPPED APPLICATION TYPES] {dict(other_apps)}"
    )

    print(
        f"[REVIEW ISSUES] {dict(issues)}"
    )

    # ========================================================
    # PRINT PREVIEW TABLES
    # ========================================================

    print_tsv(
        "FIRST REFERENCE ENTRIES",
        approved,
        args.top,
    )

    print_tsv(
        "COMBINATION-ONLY INGREDIENTS "
        "(NOT MONOTHERAPY)",
        combo_only_rows,
        args.top,
    )

    print_tsv(
        "SALT/FORM BASE-NAME SUGGESTIONS: "
        "REVIEW ONLY",
        base_review,
        args.top,
    )

    print_tsv(
        "ROWS REQUIRING PARSING REVIEW",
        review,
        args.top,
    )

    # ========================================================
    # COMPARE WITH EXISTING APPROVED DRUG REFERENCE
    # ========================================================

    before_names = set()

    if outfile.exists():

        with outfile.open(
            encoding="utf-8-sig",
            newline="",
        ) as handle:

            reader = csv.DictReader(handle)

            if (
                "active_ingredient"
                in (reader.fieldnames or [])
            ):

                before_names = {
                    norm(x.get("active_ingredient"))
                    for x in reader
                    if norm(x.get("active_ingredient"))
                }

    removed = sorted(
        before_names - selected
    )

    added = sorted(
        selected - before_names
    )

    print(
        f"[EXISTING REFERENCE] {len(before_names):,} names"
    )

    print(
        f"[WILL ADD] {len(added):,} names"
    )

    print(
        f"[WILL REMOVE] {len(removed):,} names"
    )

    print_tsv(
        "PREVIOUS-ONLY NAME EXAMPLES",
        [
            {"active_ingredient": x}
            for x in removed
        ],
        args.top,
    )

    # ========================================================
    # OPTIONAL DRY RUN
    # ========================================================

    if args.dry_run:

        print(
            "[DRY RUN] No outputs or backup written"
        )

        return

    outdir.mkdir(
        parents=True,
        exist_ok=True,
    )

    # ========================================================
    # BACK UP PREVIOUS REFERENCE
    # ========================================================

    backup = (
        outdir / "approved_drugs.before_rebuild.csv"
    )

    if (
        outfile.exists()
        and not backup.exists()
    ):

        shutil.copy2(
            outfile,
            backup,
        )

        print(
            f"[BACKUP] {backup}"
        )

    # ========================================================
    # SAVE REFERENCE AND AUDIT TABLES
    # ========================================================

    write_csv_atomic(
        outfile,
        approved,
        ["active_ingredient"],
    )

    write_csv_atomic(
        outdir / "approved_drugs_provenance.csv",
        provenance,
        list(provenance[0]),
    )

    write_csv_atomic(
        outdir / "approved_drug_combinations.csv",
        combo_rows,
        [
            "full_combination_ingredient_field",
            "n_components",
            "product_rows",
            "distinct_application_product_ids",
            "trade_name_examples",
        ],
    )

    write_csv_atomic(
        outdir / "approved_drug_combination_only_components.csv",
        combo_only_rows,
        list(provenance[0]),
    )

    write_csv_atomic(
        outdir / "approved_drug_salt_form_review.csv",
        base_review,
        [
            "fda_ingredient",
            "proposed_base_for_manual_review",
            "warning",
            "candidate_is_another_FDA_ingredient",
            "candidate_is_in_main_reference",
        ],
    )

    write_csv_atomic(
        outdir / "approved_drug_parse_review.csv",
        review,
        [
            "line",
            "issue",
            "ingredient",
            "application",
            "product_no",
            "type",
        ],
    )

    diff_rows = (
        [
            {
                "active_ingredient": x,
                "change": "added",
            }
            for x in added
        ]
        +
        [
            {
                "active_ingredient": x,
                "change": "removed",
            }
            for x in removed
        ]
    )

    write_csv_atomic(
        outdir / "approved_drug_reference_diff.csv",
        diff_rows,
        [
            "active_ingredient",
            "change",
        ],
    )

    # ========================================================
    # SAVE REFERENCE MANIFEST
    # ========================================================

    metadata = {
        "source": "FDA Orange Book Products.txt",

        "source_path": str(input_path),

        "source_sha256": digest(input_path),

        "generation_utc":
            datetime.now(timezone.utc).isoformat(),

        "reference_scope": args.reference_scope,

        "include_discontinued":
            args.include_discontinued,

        "nda_only": args.nda_only,

        "counts": {
            "input_product_rows":
                rows_read,

            "eligible_product_rows":
                rows_eligible,

            "single_product_rows":
                rows_mono,

            "combination_product_rows":
                rows_combo,

            "monotherapy_ingredient_names":
                len(mono),

            "all_component_names":
                len(all_components),

            "combination_only_component_names":
                len(combo_only),

            "approved_drugs_csv_names":
                len(approved),

            "review_rows":
                len(review),
        },

        "excluded_product_types":
            dict(excluded_types),

        "skipped_application_types":
            dict(other_apps),

        "parse_review_issues":
            dict(issues),

        "approved_drugs_csv_sha256":
            digest(outfile),

        "interpretation": (
            "Ingredient-level representation in Orange Book "
            "non-discontinued RX/OTC products; not verification "
            "of each drug-disease indication, clinical benefit "
            "or currently marketed status."
        ),
    }

    manifest_file = (
        outdir / "approved_drugs_reference_manifest.json"
    )

    with tempfile.NamedTemporaryFile(
        "w",
        encoding="utf-8",
        dir=outdir,
        prefix=".ob_manifest_",
        suffix=".tmp",
        delete=False,
    ) as handle:

        tmp = Path(handle.name)

        json.dump(
            metadata,
            handle,
            indent=2,
        )

        handle.write("\n")

    os.replace(
        tmp,
        manifest_file,
    )

    print(
        f"[SAVED] {manifest_file}"
    )

    # ========================================================
    # FINAL REPORT
    # ========================================================

    print(
        "[PASS] Orange Book reference preparation complete"
    )

    print(
        "[NOTE] Changing the reference list changes "
        "Analysis 23 denominator and coverage."
    )

    print(
        "[NOTE] FDA Orange Book is not a complete list "
        "of all FDA-approved biologics."
    )


if __name__ == "__main__":
    main()
