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

# Author: Gricey
#
# Combines the 15 per-sample quant.sf files into one TPM matrix, as needed
# by SUPPA2's psiPerEvent. Column labels (condition_replicate) come from
# meta/sample_sheet.tsv, not hardcoded. Transcript IDs are stripped of their
# Ensembl version suffix (e.g. ENSMUST00000200568.2 -> ENSMUST00000200568)
# since the GTF/.ioe files don't carry it -- without stripping, psiPerEvent
# fails to match almost every transcript.
# Pure Python standard library, no pandas.

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
    echo "  rsync -avP meta/ cl5218comp1:hfd_organoid_splicing/meta/"
    exit 1
fi

python3 - "$SAMPLE_SHEET" "$QUANT_DIR" "$OUT_DIR/tpm_matrix.tab" << 'PYEOF'
import sys, csv

sample_sheet_path, quant_dir, out_path = sys.argv[1], sys.argv[2], sys.argv[3]

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
wc -l "$OUT_DIR/tpm_matrix.tab"
