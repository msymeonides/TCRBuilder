# TCR plasmid sequencing analysis

## Setup

Install Python and BioPython: `py -m pip install biopython`

Put `analyze_tcr_plasmids.py` in a working folder with:

- `layouts/IDT-eBlocks_v4.xlsx`
- `layouts/Miniprep sequencing.xlsx`
- `10282026_1198368/` (or another `DATE_PLATEID` folder), containing `1198368_summary.csv`, `1198368_FASTA_FILES/` and `1198368_PER_BASE_DATA/`
- `10282026_1198373/` (or another `DATE_PLATEID` folder), containing the equivalent TCRβ files

Run `py analyze_tcr_plasmids.py` to process all detected datasets. To process just one, use `py analyze_tcr_plasmids.py --data-folder 10282026_1198373`.

## Outputs (two per plate)

- **`results/TCRa_REPORT.html` / `results/TCRb_REPORT.html`**: the primary, browser-readable action plan, grouped by decision, including preferred colonies, important warnings, and an interactive gallery of three 96-well layouts. Open with any web browser; it also prints to PDF.
- **`results/TCRa_ACTIONS.csv` / `results/TCRb_ACTIONS.csv`**: sortable data with one row per intended eBlock; includes alternative good and purification candidates.

The old text summary is removed when you rerun the script. Live sample-by-sample terminal progress remains enabled. Optional detailed colony-level output is available with `--diagnostics`.

## Plate gallery

The HTML report has three switchable plate views:

1. **Sequencing**: locations on the contiguous sequencing submission plate, color-coded by the result of each colony. Failed reads, incorrect/highly mixed, correct-but-needing-purification, correct pure and correct unexpected-insert results are shown separately.
2. **Miniprep**: the same colony results at their actual miniprep positions, read from the **third worksheet for TCRa** and **fourth worksheet for TCRb** of `layouts/Miniprep sequencing.xlsx`. The report refuses to display an unverified mapping if identifiers are missing or inconsistent.
3. **Gibson assembly**: the downstream decision for each intended construct at its Gibson assembly well, honoring the IDT relocation instructions for rows D/H. Each well shows the original IDT position and preferred colony or needed action.

Hover over a well to see more information. Use the tabs above the plate to switch views. The HTML works offline. **The ZIP includes the corrected workbook**, which should be saved as `layouts/Miniprep sequencing.xlsx`, replacing the earlier version. Keep `layouts/IDT-eBlocks_v4.xlsx` as before.

## Decisions

- **Ready to use**: accepted full-length insert, median insert consensus strictly above 98%; proceed from the indicated patch colony *only if backbone QC has no unresolved warning*.
- **Correct insert found in an unexpected colony**: accepted insert from another mapped original well; confirm sample/patch provenance before use.
- **Purify and screen**: accepted insert with median consensus strictly above 90%, but not above 98%; pick and sequence a single colony from the source sample.
- **Screen more colonies**: no accepted candidate among the sequenced colonies; screen others from the transformation if available, otherwise retransform the existing Gibson assembly.
- **Retransform Gibson**: no acceptable clone was recovered; retransform the existing Gibson assembly and screen fresh colonies.
- **Check sequencing inputs**: the sequencing files or identifiers were missing/inconsistent.

Each eBlock is assigned one prioritized action: a ready colony takes precedence over a purification candidate. The HTML flags selected colonies if any 500-bp backbone window (excluding insert bases) has median consensus below 99%. Such colonies must be followed up before progression. This consensus test does **not** by itself verify full-plasmid structure.

Approved sequence exceptions include the exact TCRβ junction G→A and synonymous substitutions from the specified start codon, but no nonsynonymous changes or indels. Provider-reported `Fail` is treated as plasmid absent. Thresholds and plate-to-chain mapping are editable near the top of the script. New plates can also be assigned by `--plate-chain 1234567=TCRa`.

### Unexpected-colony result wording

The report heading **Correct insert found in an unexpected colony** means that the colony contains a match to an eBlock other than the one assigned to its source sample. Recover it using the original colony identifier on the patch plate and confirm its identity by sequencing an independently recovered isolate before proceeding. If purity is insufficient, isolate and sequence a pure clone first.
