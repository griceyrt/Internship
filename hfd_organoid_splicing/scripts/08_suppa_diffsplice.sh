#!/bin/bash
#SBATCH --job-name=hfd_organoid_suppa_diffsplice
#SBATCH --partition=Lake
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --output=%j_hfd_organoid_suppa_diffsplice.log

source ~/miniconda3/etc/profile.d/conda.sh
conda activate biotools

# =============================================================================
# Phase 8 — SUPPA2 diffSplice (HFD organoid splicing project)
# Author: Gricey
# Description: The actual statistical comparison -- for each of the 4 HFD
#              timepoints vs the WT_SD baseline, and for each of the 7 event
#              types, tests whether PSI differs significantly. This is the
#              step that actually answers the project's question.
#
# SUPPA2's diffSplice compares exactly two groups at a time, so the TPM
# matrix and each of the 7 .psi files first get split into one file per
# condition (5 conditions total: WT_SD, W8/18/24/42_WHFD), matching the
# proven approach from orthogonal_validation/scripts/run_suppa.sh (there:
# WT vs KO, 2 groups; here: 5 groups, but only compared pairwise against the
# WT_SD baseline -- not all 10 possible pairs).
#
# Splitting done in Python (not awk like run_suppa.sh) because the header
# row and data rows in our TPM matrix / .psi files have different column
# counts (header = N sample labels only, data rows = N labels + 1 leading
# id column) -- awk's NR==1 vs NR>1 special-casing gets fragile with that
# offset, plain Python string splitting is more robust and readable.
#
# diffSplice flags match the PROVEN working command from run_suppa.sh:
#   --method empirical -gc (gene-corrected p-values)
#
# STORAGE: everything on $SCRATCH, only this script lives in $HOME.
# =============================================================================

set -euo pipefail

CODE_DIR="$HOME/hfd_organoid_splicing"
SCRATCH_DIR="/scratch/Lake/$USER/hfd_organoid_splicing"
SAMPLE_SHEET="$CODE_DIR/meta/sample_sheet.tsv"
TPM_MATRIX="$SCRATCH_DIR/results/suppa/tpm_matrix.tab"
PSI_DIR="$SCRATCH_DIR/results/suppa/psi"
EVENTS_DIR="$SCRATCH_DIR/results/suppa/events"
SPLIT_DIR="$SCRATCH_DIR/results/suppa/split_by_condition"
DIFF_DIR="$SCRATCH_DIR/results/suppa/diffsplice"
mkdir -p "$SPLIT_DIR" "$DIFF_DIR"

SUPPA="python3 $CODE_DIR/SUPPA/suppa.py"

echo "=== Splitting TPM matrix + all 7 PSI files by condition ==="
python3 - "$SAMPLE_SHEET" "$TPM_MATRIX" "$PSI_DIR" "$SPLIT_DIR" << 'PYEOF'
import sys, csv, glob, os

sample_sheet_path, tpm_path, psi_dir, split_dir = sys.argv[1:5]

# Map: condition -> ordered list of sample labels (e.g. WT_SD -> [WT_SD_1, WT_SD_2, WT_SD_3])
conditions = {}
with open(sample_sheet_path) as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        label = f'{row["condition"]}_{row["replicate"]}'
        conditions.setdefault(row["condition"], []).append(label)

print("Conditions found:", {c: len(v) for c, v in conditions.items()})

def split_matrix(in_path, out_prefix):
    """Split a TPM/PSI matrix (header = labels only, data rows = id + values)
    into one file per condition, preserving row order."""
    with open(in_path) as f:
        header_labels = f.readline().rstrip("\n").split("\t")
        label_to_col = {label: i for i, label in enumerate(header_labels)}

        # Open one output file handle per condition, write its header
        out_files = {}
        for cond, labels in conditions.items():
            out_path = f"{out_prefix}_{cond}.tab"
            out_files[cond] = open(out_path, "w")
            out_files[cond].write("\t".join(labels) + "\n")

        for line in f:
            fields = line.rstrip("\n").split("\t")
            row_id = fields[0]
            values = fields[1:]
            for cond, labels in conditions.items():
                cols = [values[label_to_col[label]] for label in labels]
                out_files[cond].write(row_id + "\t" + "\t".join(cols) + "\n")

        for fh in out_files.values():
            fh.close()

# Split the TPM matrix
split_matrix(tpm_path, f"{split_dir}/tpm")
print(f"Split TPM matrix -> {split_dir}/tpm_<condition>.tab")

# Split every .psi file
for psi_path in sorted(glob.glob(f"{psi_dir}/*.psi")):
    event_name = os.path.basename(psi_path)[:-4]  # strip .psi
    split_matrix(psi_path, f"{split_dir}/{event_name}")
    print(f"Split {event_name}.psi -> {split_dir}/{event_name}_<condition>.tab")
PYEOF

echo ""
echo "=== Running diffSplice: each HFD timepoint vs WT_SD, all 7 event types ==="

for TIMEPOINT in W8_WHFD W18_WHFD W24_WHFD W42_WHFD; do
    for ioefile in "$EVENTS_DIR"/*.ioe; do
        event_name="$(basename "$ioefile" .ioe)"
        comparison="${TIMEPOINT}_vs_WT_SD"

        WT_PSI="$SPLIT_DIR/${event_name}_WT_SD.tab"
        TP_PSI="$SPLIT_DIR/${event_name}_${TIMEPOINT}.tab"
        WT_TPM="$SPLIT_DIR/tpm_WT_SD.tab"
        TP_TPM="$SPLIT_DIR/tpm_${TIMEPOINT}.tab"

        echo "=== diffSplice: ${comparison} / ${event_name} ==="
        $SUPPA diffSplice \
            --method empirical \
            --input "$ioefile" \
            --psi "$WT_PSI" "$TP_PSI" \
            --tpm "$WT_TPM" "$TP_TPM" \
            -gc \
            -o "$DIFF_DIR/${comparison}_${event_name}"
    done
done

echo ""
echo "=== Done. diffSplice results: $DIFF_DIR ==="
ls -lh "$DIFF_DIR" | head -20
echo "..."
echo "Total files: $(ls "$DIFF_DIR" | wc -l)"
