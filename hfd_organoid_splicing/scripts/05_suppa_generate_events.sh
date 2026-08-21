#!/bin/bash
#SBATCH --job-name=hfd_organoid_suppa_events
#SBATCH --partition=Lake-short
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=%j_hfd_organoid_suppa_events.log

source ~/miniconda3/etc/profile.d/conda.sh
conda activate biotools

# Author: Gricey
#
# Generates local AS events from the GTF. SUPPA2 isn't a conda/pip package
# in biotools -- git-cloned once into $HOME, reused by later SUPPA2 scripts.
# -e SE SS MX RI FL expands to all 7 event types (SS -> A5+A3, FL -> AF+AL).
# -b S = strict boundary mode.
#
# suppa.py unconditionally imports its clustering module at the top of the
# file, which needs scikit-learn/scipy -- broken in this env and unrelated
# to anything I use (generateEvents/psiPerEvent/diffSplice, not
# clusterEvents). The block below patches that one import to be optional
# rather than fighting the dependency; idempotent, safe to re-run.

set -euo pipefail

CODE_DIR="$HOME/hfd_organoid_splicing"
SCRATCH_DIR="/scratch/Lake/$USER/hfd_organoid_splicing"
GTF="$SCRATCH_DIR/data/reference/Mus_musculus.GRCm39.116.gtf"
EVENTS_DIR="$SCRATCH_DIR/results/suppa/events"
mkdir -p "$EVENTS_DIR"

echo "=== Setting up SUPPA2 (clone once, reused by future scripts) ==="
if [ ! -f "$CODE_DIR/SUPPA/suppa.py" ]; then
    git clone https://github.com/comprna/SUPPA.git "$CODE_DIR/SUPPA"
else
    echo "SUPPA already present at $CODE_DIR/SUPPA, skipping clone."
fi
SUPPA="python3 $CODE_DIR/SUPPA/suppa.py"

echo "=== Patching suppa.py: make sklearn/clusterEvents import optional ==="
if ! grep -q "clusterAnalysis = None" "$CODE_DIR/SUPPA/suppa.py"; then
    python3 - "$CODE_DIR/SUPPA/suppa.py" << 'PYEOF'
import sys
path = sys.argv[1]
with open(path) as f:
    content = f.read()

content = content.replace(
    "import eventClusterer as clusterAnalysis\n",
    "try:\n    import eventClusterer as clusterAnalysis\n"
    "except Exception:\n    clusterAnalysis = None\n"
)

content = content.replace(
    '''# eventClusterer parser
eventClustererSubparser = subparsers.add_parser(
    "clusterEvents", parents=[clusterAnalysis.parser],
    help="Calculates clusters of events across conditions.")
eventClustererSubparser.set_defaults(which="clusterEvents")''',
    '''# eventClusterer parser (skipped if sklearn/clusterAnalysis unavailable)
if clusterAnalysis is not None:
    eventClustererSubparser = subparsers.add_parser(
        "clusterEvents", parents=[clusterAnalysis.parser],
        help="Calculates clusters of events across conditions.")
    eventClustererSubparser.set_defaults(which="clusterEvents")'''
)

with open(path, "w") as f:
    f.write(content)

print("Patched suppa.py successfully.")
PYEOF
else
    echo "Already patched, skipping."
fi

echo "=== Running generateEvents ==="
$SUPPA generateEvents \
    -i "$GTF" \
    -o "$EVENTS_DIR/events" \
    -f ioe \
    -e SE SS MX RI FL \
    -b S

echo ""
echo "=== Done. Event files: $EVENTS_DIR ==="
ls -lh "$EVENTS_DIR"
