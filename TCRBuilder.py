import pandas as pd
import subprocess
import sys

# TCRBuilder, version 1.1
#
# This script converts TCR clonotype data from an Excel file into a format suitable for the stitchr/thimble pipeline,
# runs thimble to generate TCR constructs, and outputs final TCR chain assemblies with homology arms for HiFi cloning.
# Written by Mel Symeonides at the University of Vermont, October 2025.
# Stitchr documentation can be found at https://jamieheather.github.io/stitchr/

# --- Configuration ---
basename = "PBMC-Rota-D33-a4b7pos-0-10_new"   # Your input Excel filename (without the .xlsx extension)
input_file = f"{basename}.xlsx"
output_file = f"{basename}_converted.tsv"  # prepared for thimble input
thimble_output = f"{basename}_thimble.tsv"  # output from thimble
final_excel = f"{basename}_final-assembly.xlsx"  # final output Excel with HiFi constructs
cluster_analysis_output = f"{basename}_cluster_analysis.xlsx" # final output for cluster analysis
species = "human"  # Your dataset species for thimble, e.g., "human" or "mouse"
trbc = "TRBC2" # Default TRBC to use, set to either TRBC1 or TRBC2

def check_stitchr_installed():
    try:
        subprocess.run(["thimble", "--help"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
    except (subprocess.CalledProcessError, FileNotFoundError):
        print("Aborting - please install stitchr and IMGTgeneDL (pip install stitchr IMGTgeneDL) to proceed.")
        sys.exit(1)

def convert_clonotypes_to_thimble_base():
    df = pd.read_excel(input_file, dtype=str)

    # Columns to check for NA/NaN
    cols_to_check = ["v_gene_a", "j_gene_a", "cdr3_a", "cdr3_b", "v_gene_b", "j_gene_b"]

    # Normalize strings for reliable filtering
    df_cleaned = df[cols_to_check].map(lambda x: str(x).strip().lower() if pd.notnull(x) else x)

    # Drop rows with any 'na' or actual NaN values in the specified columns
    mask = ~((df_cleaned.isin(['na']) | df_cleaned.isna()).any(axis=1))
    df = df[mask].copy()

    # Add constant TRAC and TRBC only to remaining rows
    df["TRAC"] = "TRAC"
    df["TRBC"] = trbc

    thimble_columns = [
        "TCR_name", "TRAV", "TRAJ", "TRA_CDR3",
        "TRBV", "TRBJ", "TRB_CDR3",
        "TRAC", "TRBC",
        "TRA_leader", "TRB_leader",
        "Linker", "Link_order",
        "TRA_5_prime_seq", "TRA_3_prime_seq",
        "TRB_5_prime_seq", "TRB_3_prime_seq"
    ]

    thimble_rows = []

    for _, row in df.iterrows():
        thimble_row = {
            "TCR_name": row["Clonotype_ID"],
            "TRAV": row.get("v_gene_a", ""),
            "TRAJ": row.get("j_gene_a", ""),
            "TRA_CDR3": row.get("cdr3_a", ""),
            "TRBV": row.get("v_gene_b", ""),
            "TRBJ": row.get("j_gene_b", ""),
            "TRB_CDR3": row.get("cdr3_b", ""),
            "TRAC": row["TRAC"],
            "TRBC": row["TRBC"],
            "TRA_leader": "",
            "TRB_leader": "",
            "Linker": "",
            "Link_order": "",
            "TRA_5_prime_seq": "",
            "TRA_3_prime_seq": "",
            "TRB_5_prime_seq": "",
            "TRB_3_prime_seq": ""
        }
        thimble_rows.append(thimble_row)

    df_thimble = pd.DataFrame(thimble_rows)[thimble_columns]
    df_thimble.to_csv(output_file, sep="\t", index=False)

    print(f"Output written to '{output_file}' with {len(df_thimble)} rows.")

    # Run thimble as a subprocess
    try:
        subprocess.run([
            "thimble",
            "-i", output_file,
            "-o", thimble_output,
            "-s", species
        ], check=True)
        print(f"Thimble run complete. Output TSV: '{thimble_output}'")
    except Exception as e:
        print(f"Thimble run failed: {e}")
        return

    # Read thimble output TSV
    thimble_df = pd.read_csv(thimble_output, sep="\t")

    # Prepare final assembly output
    prefix = "GAGTCGCCCGGGGGGGATCGCTCGAGACC"    # 5` homology arm for HiFi cloning into pHR vector
    suffix = "GGATCCGGAGCTACTAACTTCAGCCT"       # 3` homology arm for HiFi cloning into pHR vector
    final_rows = []

    for _, row in thimble_df.iterrows():
        tcr_name = row.get("TCR_name", "")
        tra_nt = row.get("TRA_nt", "")
        trb_nt = row.get("TRB_nt", "")
        final_rows.append({
            "Name": f"{basename}_TCR{tcr_name}",   # Change TCR naming format as needed
            "TRA_construct": f"{prefix}{tra_nt}{suffix}" if pd.notna(tra_nt) and tra_nt else "",
            "TRB_construct": f"{prefix}{trb_nt}{suffix}" if pd.notna(trb_nt) and trb_nt else ""
        })

    final_df = pd.DataFrame(final_rows, columns=["Name", "TRA_construct", "TRB_construct"])
    final_df.to_excel(final_excel, index=False)
    print(f"Final assembly Excel written: '{final_excel}'")

def calculate_cluster_frequencies():
    """
    Calculates clonotype frequencies within clusters and saves the analysis to an Excel file.
    """
    try:
        # Read both sheets from the input Excel file
        xls = pd.read_excel(input_file, sheet_name=None)
        df_main = xls[list(xls.keys())[0]]  # Main data is the first sheet
        df_clusters = xls.get("Clusters")

        if df_clusters is None:
            print("'Clusters' worksheet not found in the input file. Skipping cluster frequency analysis.")
            return

        # Prepare cluster total cells for lookup
        # Assumes cluster names in 'Clusters' sheet are like 'C0', 'C1', etc.
        # The user mentioned the column is "Cluster", so we will use that.
        # And that the format is just the number as text, e.g. "0", "1", etc.
        df_clusters['Cluster'] = df_clusters['Cluster'].astype(str)
        df_clusters['Cells'] = pd.to_numeric(df_clusters['Cells'])
        cluster_totals = df_clusters.set_index("Cluster")["Cells"].to_dict()

        # New columns to be added
        new_cols = {
            "Primary cluster": [], "Primary cluster %": [],
            "Secondary cluster": [], "Secondary cluster %": [],
            "Tertiary cluster": [], "Tertiary cluster %": []
        }

        for _, row in df_main.iterrows():
            cluster_id_str = row.get("cluster_ID_w_count", "")
            clusters = []
            if isinstance(cluster_id_str, str) and cluster_id_str:
                # E.g., "C0(113),C3(27),C5(6)"
                parts = cluster_id_str.split(',')
                for part in parts:
                    try:
                        cluster_name = part.split('(')[0]
                        cell_count = int(part.split('(')[1][:-1])
                        clusters.append((cluster_name, cell_count))
                    except (ValueError, IndexError):
                        continue # Skip malformed entries

            # Fill cluster data for up to 3 clusters
            for i in range(3):
                if i < len(clusters):
                    cluster_name, cell_count = clusters[i]
                    cluster_num = ''.join(filter(str.isdigit, cluster_name))
                    total_cells = cluster_totals.get(cluster_num)

                    if total_cells and total_cells > 0:
                        percentage = (cell_count / total_cells)
                    else:
                        percentage = 0.0

                    if i == 0:
                        new_cols["Primary cluster"].append(cluster_num)
                        new_cols["Primary cluster %"].append(percentage)
                    elif i == 1:
                        new_cols["Secondary cluster"].append(cluster_num)
                        new_cols["Secondary cluster %"].append(percentage)
                    elif i == 2:
                        new_cols["Tertiary cluster"].append(cluster_num)
                        new_cols["Tertiary cluster %"].append(percentage)
                else:
                    # No data for this cluster level
                    if i == 0:
                        new_cols["Primary cluster"].append(None)
                        new_cols["Primary cluster %"].append(None)
                    elif i == 1:
                        new_cols["Secondary cluster"].append(None)
                        new_cols["Secondary cluster %"].append(None)
                    elif i == 2:
                        new_cols["Tertiary cluster"].append(None)
                        new_cols["Tertiary cluster %"].append(None)

        # Add new columns to the main dataframe
        df_main["Primary cluster"] = new_cols["Primary cluster"]
        df_main["Primary cluster %"] = new_cols["Primary cluster %"]
        df_main["Secondary cluster"] = new_cols["Secondary cluster"]
        df_main["Secondary cluster %"] = new_cols["Secondary cluster %"]
        df_main["Tertiary cluster"] = new_cols["Tertiary cluster"]
        df_main["Tertiary cluster %"] = new_cols["Tertiary cluster %"]

        # Save to a new Excel file with both sheets
        with pd.ExcelWriter(cluster_analysis_output) as writer:
            df_main.to_excel(writer, sheet_name=list(xls.keys())[0], index=False)
            df_clusters.to_excel(writer, sheet_name="Clusters", index=False)

        print(f"Cluster frequency analysis written to '{cluster_analysis_output}'")

    except Exception as e:
        print(f"An error occurred during cluster frequency analysis: {e}")


if __name__ == "__main__":
    check_stitchr_installed()
    convert_clonotypes_to_thimble_base()
    calculate_cluster_frequencies()
