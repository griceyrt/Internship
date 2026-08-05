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

# =============================================================================
# Phase 6 — SUPPA2 generateEvents (HFD organoid splicing project)
# Author: Gricey
# Description: Generates local alternative splicing events from the
#              release-116 GTF. SUPPA2 isn't pip-installed anywhere in the
#              biotools env (checked, nothing found) -- it's used as a
#              git-cloned script in every prior project
#              (https://github.com/comprna/SUPPA.git), so this script clones
#              it fresh into this project too, matching that convention.
#
# Flags copied exactly from the PROVEN working command in
# orthogonal_validation/scripts/run_suppa.sh:
#   -e SE SS MX RI FL  -- 5 input codes, expands to 7 output event types
#                          (SS -> A5 + A3, FL -> AF + AL)
#   -b S                -- strict boundary mode (matches every prior
#                           project's default; project 1's "variable"
#                           boundary mode was specific to that project's own
#                           boundary-analysis question, not relevant here)
#   -f ioe              -- local-event .ioe output format (what psiPerEvent/
#                           diffSplice expect downstream)
#
# STORAGE: SUPPA2 itself is code, so it lives in $HOME (cloned once, reused
# by every future SUPPA2 script). Its output (the actual event files) goes
# on $SCRATCH like all other data. Only this script's own contents are in
# $HOME/hfd_organoid_splicing/scripts/ as usual.
#
# PATCH (added 2026-08-05, fixed same day): suppa.py unconditionally imports
# its clustering module (eventClusterer) at the top of the file, which
# requires sklearn -- and sklearn's scipy dependency has a real bug/
# incompatibility in this biotools env (ValueError from
# scipy.special._multiufuncs, survived two attempts to fix numpy/scipy/
# scikit-learn versions via conda). We never use SUPPA2's clusterEvents
# command (only generateEvents/psiPerEvent/diffSplice), so rather than keep
# fighting that dependency, this script patches the clone to make the
# import optional. NOTE: must catch `except Exception`, not
# `except ImportError` -- the actual failure is a ValueError raised deep in
# scipy's own module-level code during the import, not a clean "module not
# found", so a narrower except clause silently fails to catch it (this bit
# us on the first attempt). Idempotent (checks if already patched before
# re-patching), so safe to re-run.
# =============================================================================

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
echo "=== (we never use clusterEvents; unblocks generateEvents/psiPerEvent/diffSplice"
echo "===  from an unrelated sklearn/scipy incompatibility in this env)"
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
    '''# eventClusterer parser (skipped if sklearn/clusterAnalysis unavailable --
# not needed for generateEvents/psiPerEvent/diffSplice, only for
# clusterEvents, which this project doesn't use)
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
