import pandas as pd
import glob
import os
import re

# TCRCompare, Version 1.0
#
# This script compares TCR clonotype datasets between two groups belonging to a specimen (or multiple specimens).
# For example, you may have clonotype data for a specimen named "wildtype" in two groups: "neg" and "pos". This script
# will identify clonotypes common to both groups and compare their prevalence in each group by absolute count and by
# prevalence (rank) within each group. You can also use this to identify clonotypes unique to one group.
#
# Note 1: A clonotype ID and a count of 0 indicates that the clonotype was not found in that group. Please ensure that
# the input files do not use 0 as a valid clonotype ID and do not include any rows with a count of 0.
# Note 2: Ensure your input files are placed in your working folder (the same folder where this script is located).
# Note 3: At present you must have exactly two groups per specimen. If the working folder contains more than two groups
# for a specimen, the script will skip that specimen and print a warning message.
#
# Written by Mel Symeonides at the University of Vermont, July 2025.


# Edit the variable definition below to specify your filename template, using the {specimen} and {group} placeholders.
# Example: Assume you have the datasets "TCR-clonotypes-wildtype-pos.xlsx" and "TCR-clonotypes-wildtype-neg.xlsx".
# Define your template as "TCR-clonotypes-{specimen}-{group}.xlsx"; this will extract the specimen "wildtype" and
# groups "pos" and "neg" from those filenames and will compare "pos" to "neg" within "wildtype".
# If you also have the same pair of datasets but for the specimen "mutant", this same template will serve to compare
# those two groups for that specimen as well, but the comparison will be run independently for "wildtype" and "mutant".
FILENAME_TEMPLATE = "PBMC-{specimen}-*-{group}-*_new.xlsx"


# The variable below lists the columns of interest that will be compared between groups. Ensure that these columns are
# present in your input files and named exactly as shown. The columns "Clonotype_ID" and "Count" (again named exactly
# as shown) are also expected to be present. The column "d_gene_b" is ignored as the D segment of the beta chain is
# entirely contained within the CDR3 region.
cols_of_interest = ["v_gene_a", "j_gene_a", "cdr3_a", "cdr3_b", "v_gene_b", "j_gene_b"]


def infer_separator(template):
    m = re.search(r"\{specimen}(\W)", template)
    if m:
        return m.group(1)
    m = re.search(r"\{group}(\W)", template)
    if m:
        return m.group(1)
    if '-' in template:
        return '-'
    elif '_' in template:
        return '_'
    else:
        raise ValueError("Could not infer separator from template.")

def get_template_positions(template, sep):
    template_no_ext = os.path.splitext(template)[0]
    parts = template_no_ext.split(sep)
    specimen_idx = None
    group_idx = None
    for i, part in enumerate(parts):
        if "{specimen}" in part:
            specimen_idx = i
        if "{group}" in part:
            group_idx = i
    return specimen_idx, group_idx

def extract_specimen_and_group(template, filename):
    sep = infer_separator(template)
    filename_no_ext = os.path.splitext(filename)[0]
    fname_parts = filename_no_ext.split(sep)
    specimen_idx, group_idx = get_template_positions(template, sep)
    try:
        specimen_part = fname_parts[specimen_idx]
        group = fname_parts[group_idx]
        template_no_ext = os.path.splitext(template)[0]
        template_parts = template_no_ext.split(sep)
        specimen_template = template_parts[specimen_idx]
        if "{specimen}" in specimen_template:
            prefix = specimen_template.split("{specimen}")[0]
            if prefix and specimen_part.startswith(prefix):
                specimen_part = specimen_part[len(prefix):]
        return specimen_part, group
    except (IndexError, TypeError):
        return None, None

def build_glob_pattern(template):
    return template.replace("{specimen}", "*").replace("{group}", "*")

def load_files():
    pattern = build_glob_pattern(FILENAME_TEMPLATE)
    files = glob.glob(pattern)
    print(f"Files found with pattern '{pattern}': {files}")  # Debug print
    specimens = {}
    for f in files:
        fname = os.path.basename(f)
        specimen, group = extract_specimen_and_group(FILENAME_TEMPLATE, fname)
        print(f"Extracted specimen: '{specimen}', group: '{group}' from file: '{fname}'")  # Debug print
        if specimen and group:
            specimens.setdefault(specimen, {})[group] = f
    if not specimens:
        print("No input files detected. Please check your FILENAME_TEMPLATE and input files.")
    return specimens

def compare_specimen(file_a, file_b, group_a, group_b, specimen):
    df_a = pd.read_excel(file_a, dtype=str).fillna("")
    df_b = pd.read_excel(file_b, dtype=str).fillna("")

    missing_cols_a = [col for col in cols_of_interest if col not in df_a.columns]
    missing_cols_b = [col for col in cols_of_interest if col not in df_b.columns]
    if missing_cols_a or missing_cols_b:
        print(f"Skipping specimen {specimen}: missing columns in {group_a}: {missing_cols_a}, {group_b}: {missing_cols_b}")
        return

    def valid_row(row):
        for col in cols_of_interest:
            val = str(row[col]).strip()
            if val == "" or val.lower() == "na":
                return False
        return True

    df_a = df_a[df_a.apply(valid_row, axis=1)].copy()
    df_b = df_b[df_b.apply(valid_row, axis=1)].copy()

    df_a["Count_int"] = df_a["Count"].apply(lambda x: int(x) if str(x).isdigit() else 0)
    df_a = df_a.sort_values("Count_int", ascending=False)
    df_a["rank"] = range(1, len(df_a) + 1)
    rank_map_a = {tuple(row[col] for col in cols_of_interest): row["rank"] for _, row in df_a.iterrows()}

    df_b["Count_int"] = df_b["Count"].apply(lambda x: int(x) if str(x).isdigit() else 0)
    df_b = df_b.sort_values("Count_int", ascending=False)
    df_b["rank"] = range(1, len(df_b) + 1)
    rank_map_b = {tuple(row[col] for col in cols_of_interest): row["rank"] for _, row in df_b.iterrows()}

    def make_key(row):
        return tuple(row[col] for col in cols_of_interest)

    rows_a = {make_key(row): row for _, row in df_a.iterrows()}
    rows_b = {make_key(row): row for _, row in df_b.iterrows()}

    all_keys = set(rows_a.keys()).union(rows_b.keys())
    output_rows = []

    for key in all_keys:
        in_a = key in rows_a
        in_b = key in rows_b
        common = in_a and in_b

        row_a = rows_a.get(key, {})
        row_b = rows_b.get(key, {})

        def to_int(val):
            try:
                return int(val)
            except (ValueError, TypeError):
                return 0

        rank_a = int(rank_map_a.get(key, 0))
        rank_b = int(rank_map_b.get(key, 0))
        rank_diff = rank_b - rank_a

        count_a = to_int(row_a.get("Count", 0)) if in_a else 0
        count_b = to_int(row_b.get("Count", 0)) if in_b else 0
        count_eq = (count_a == count_b == 1)

        # Set presence value
        presence = "common" if common else "unique"

        output_rows.append({
            f"{group_a} ID": to_int(row_a.get("Clonotype_ID", 0)) if in_a else 0,
            f"{group_b} ID": to_int(row_b.get("Clonotype_ID", 0)) if in_b else 0,
            "presence": presence,
            **{col: key[i] for i, col in enumerate(cols_of_interest)},
            f"{group_a} count": count_a,
            f"{group_b} count": count_b,
            "count=1 in both": count_eq,
            f"{group_a} rank": rank_a,
            f"{group_b} rank": rank_b,
            "rank diff": rank_diff
        })

    col_order = [
        f"{group_a} ID", f"{group_b} ID", "presence",
        *cols_of_interest,
        f"{group_a} count", f"{group_b} count", "count=1 in both",
        f"{group_a} rank", f"{group_b} rank", "rank diff"
    ]
    out_name = f"TCRCompare_{specimen}-{group_a}-{group_b}.xlsx"
    pd.DataFrame(output_rows)[col_order].to_excel(out_name, index=False)
    print(f"Comparison for specimen {specimen} ({group_a} vs {group_b}) written to {out_name}")

def main():
    specimens = load_files()
    for specimen, groups in specimens.items():
        group_names = sorted(groups.keys())
        if len(group_names) != 2:
            print(f"Skipping specimen '{specimen}': found {len(group_names)} groups ({group_names}), but exactly two are required.")
            continue
        group_a, group_b = group_names
        compare_specimen(groups[group_a], groups[group_b], group_a, group_b, specimen)

if __name__ == "__main__":
    main()
