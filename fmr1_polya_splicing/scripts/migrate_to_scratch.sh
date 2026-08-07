#!/bin/bash
#SBATCH --job-name=migrate_to_scratch
#SBATCH --partition=Lake
#SBATCH --time=02:00:00
#SBATCH --mem=4G
#SBATCH --output=migrate_to_scratch_%j.log

# Moves large data OUT of $HOME (shared 41T igfl group quota, 76% full)
# INTO $SCRATCH/Lake (473T free, no quota) and leaves a symlink behind
# at the original path, so no existing script needs its paths edited.
#
# SCOPE: only touches these two subfolders inside ~/fmr1_polya_splicing.
# Nothing outside this project (e.g. ~/hfd_organoid_splicing) is touched.
#
# Safety: rsync-copies first, verifies byte-identical size, and only
# THEN deletes the original + creates the symlink. If sizes don't match,
# the original is left untouched and the mismatch is printed so it can
# be checked manually -- nothing is deleted on an unverified copy.

set -e

PROJECT_HOME=~/fmr1_polya_splicing
SCRATCH_BASE=/scratch/Lake/grangel/fmr1_polya_splicing

mkdir -p "$SCRATCH_BASE"

migrate() {
    local name="$1"
    local src="$PROJECT_HOME/$2"
    local dst="$SCRATCH_BASE/$2"

    if [ -L "$src" ]; then
        echo "[$name] already a symlink, skipping."
        return
    fi
    if [ ! -e "$src" ]; then
        echo "[$name] not found at $src, skipping."
        return
    fi

    echo "[$name] copying $src -> $dst ..."
    mkdir -p "$(dirname "$dst")"
    rsync -a --info=progress2 "$src"/ "$dst"/

    src_size=$(du -sb "$src" | cut -f1)
    dst_size=$(du -sb "$dst" | cut -f1)

    if [ "$src_size" -eq "$dst_size" ]; then
        echo "[$name] verified ($src_size bytes match). Replacing with symlink."
        rm -rf "$src"
        ln -s "$dst" "$src"
        echo "[$name] done: $src -> $dst"
    else
        echo "[$name] SIZE MISMATCH (src=$src_size dst=$dst_size)."
        echo "[$name] NOT deleting original. Check $dst manually before retrying."
    fi
}

migrate "minimap2 BAMs (103G)"   "results/gse188840_minimap2_2026-08-03"
migrate "fastq temp (42G)"       "temp"

# --- Straight delete (not migrated): oarfish_test, superseded by the real
# 6-sample production run already saved in results/gse188840_oarfish_2026-08-04.
# Confirmed by Gricey 2026-08-07 as safe to delete outright, not just move.
OARFISH_TEST="$PROJECT_HOME/results/oarfish_test"
if [ -d "$OARFISH_TEST" ]; then
    echo "Deleting superseded test run: $OARFISH_TEST"
    rm -rf "$OARFISH_TEST"
    echo "Deleted."
else
    echo "$OARFISH_TEST not found, skipping."
fi

echo ""
echo "=== Migration finished. Disk state: ==="
df -h /home/grangel /scratch/Lake
