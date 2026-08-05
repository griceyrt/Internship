#!/bin/bash
#SBATCH --job-name=hfd_organoid_tpm_matrix
#SBATCH --partition=Lake-short
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=%j_hfd_organoid_tpm_matrix.log

source ~/miniconda3/etc/profile.d/conda.sh
conda activate biotools

# =============================================================================
# Combine per-sample TPM into one matrix (HFD organoid splicing project)
# Author: Gricey
# Description: SUPPA2's psiPerEvent needs ONE table with all 15 samples' TPM
#              values side by side, not 15 separate quant.sf files. This is
#              a step Desmond's pipeline never needed (he worked per-sample
#              via tximeta/DESeq2 in R) -- built fresh for the splicing
#              analysis, sits between Phase 6 (generateEvents) and Phase 7
#              (psiPerEvent).
#
# Column naming matches Desmond's own convention exactly (from
# rna_seq_pipeline.r: "WT_SD_1", "W8_WHFD_1", etc. -- condition_replicate),
# read from meta/sample_sheet.tsv rather than hardcoded, so it can't drift
# out of sync with the confirmed Lib#->condition mapping.
#
# Pure Python standard library only (no pandas) -- deliberately avoiding a
# new dependency given how much friction this env's packages have already
# caused (see reference_biotools_env_gotchas memory).
#
# VERSION-STRIPPING (added 2026-08-05, after first run failed): salmon's
# quant.sf "Name" column carries the transcript ID WITH version suffix
# (e.g. ENSMUST00000200568.2), because that's what's in the cDNA FASTA
# headers. But the GTF's transcript_id attribute -- and therefore every
# .ioe event file SUPPA2 generated -- does NOT carry the version suffix
# (just ENSMUST00000200568). psiPerEvent failed to find almost every
# transcript as a result. Fix: strip everything from the first "." onward
# before writing each transcript ID, so both files agree.
#
# PREREQUISITE: meta/sample_sheet.tsv must be synced to
# $HOME/hfd_organoid_splicing/meta/ (same rsync pattern as scripts/) --
# not yet done as of writing this script, do it before running.
#
# STORAGE: everything on $SCRATCH except the small sample sheet, which
# lives in $HOME alongside the scripts (tracked in git, treated as
# code/config, not data).
# =============================================================================

set -euo pipefail

CODE_DIR="$HOME/hfd_organoid_splicing"
SCRATCH_DIR="/scratch/Lake/$USER/hfd_organoid_splicing"
SAMPLE_SHEET="$CODE_DIR/meta/sample_sheet.tsv"
QUANT_DIR="$SCRATCH_DIR/results/salmon"
OUT_DIR="$SCRATCH_DIR/results/suppa"
mkdir -p "$OUT_DIR"

if [ ! -f "$SAMPLE_SHEET" ]; then
    echo "ERROR: sample sheet not found at $SAMPLE_SHEET"
    echo "Sync it up first, from your Mac:"
    echo "  rsync -avP \"/Users/gricey/Desktop/Internship/hfd_organoid_splicing/meta/\" cl5218comp1:hfd_organoid_splicing/meta/"
    exit 1
fi

python3 - "$SAMPLE_SHEET" "$QUANT_DIR" "$OUT_DIR/tpm_matrix.tab" << 'PYEOF'
import sys, csv

sample_sheet_path, quant_dir, out_path = sys.argv[1], sys.argv[2], sys.argv[3]

# Read sample sheet: sample, condition, replicate, fastq_file
samples = []
with open(sample_sheet_path) as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        samples.append({
            "sample": row["sample"],
            "label": f'{row["condition"]}_{row["replicate"]}',
        })

print(f"Found {len(samples)} samples in sample sheet.")

# Read each sample's quant.sf, extract Name (col 0) + TPM (col 3)
transcript_ids = None
tpm_columns = {}  # label -> list of TPM values, same order as transcript_ids

for s in samples:
    quant_path = f'{quant_dir}/{s["sample"]}/quant.sf'
    ids = []
    tpms = []
    with open(quant_path) as f:
        f.readline()  # skip header
        for line in f:
            fields = line.rstrip("\n").split("\t")
            ids.append(fields[0].split(".")[0])  # strip version suffix
            tpms.append(fields[3])

    if transcript_ids is None:
        transcript_ids = ids
    elif ids != transcript_ids:
        sys.exit(f"ERROR: transcript ID order mismatch in {quant_path} vs "
                  f"first sample -- quant.sf files aren't from the same "
                  f"salmon index, or something else is wrong.")

    tpm_columns[s["label"]] = tpms
    print(f'  {s["sample"]} -> column "{s["label"]}" ({len(tpms)} transcripts)')

# Write combined matrix: header = sample labels, rows = transcript_id + TPMs
labels = [s["label"] for s in samples]
with open(out_path, "w") as out:
    out.write("\t".join(labels) + "\n")
    for i, tid in enumerate(transcript_ids):
        row = [tid] + [tpm_columns[label][i] for label in labels]
        out.write("\t".join(row) + "\n")

print(f"Wrote combined TPM matrix: {out_path} "
      f"({len(transcript_ids)} transcripts x {len(labels)} samples)")
PYEOF

echo ""
echo "=== Done. TPM matrix: $OUT_DIR/tpm_matrix.tab ==="
head -2 "$OUT_DIR/tpm_matrix.tab" | cut -c1-200
wc -l "$OUT_DIR/tpm_matrix.tab"
