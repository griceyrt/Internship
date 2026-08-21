#!/bin/bash
#SBATCH --job-name=migrate_fast5
#SBATCH --partition=Lake
#SBATCH --time=12:00:00
#SBATCH --mem=4G
#SBATCH --output=%j_migrate_fast5.log

# Run this AFTER download_fast5_japan.sh has fully finished and its final
# `du -sh` output looks right (~3.4TB). Moves the fast5 data off the
# periodically-wiped /Xnfs/igfldb staging area into its real home,
# /scratch/Lake/grangel/fmr1_polya_splicing/fast5/, verifying byte size
# before deleting the igfldb copy.
#
# NOTE on compression: fast5 (HDF5) files are often already internally
# compressed (VBZ or gzip, depending on the basecaller/kit version) --
# gzip/tar.gz on top typically saves very little (sometimes nothing) for
# real cost in CPU time on 3.4TB. This script does a plain rsync copy, no
# extra compression layer. If you want to check whether compression is
# worth it here, test on ONE file first:
#   du -h somefile.fast5; gzip -k somefile.fast5; du -h somefile.fast5.gz
# and only bother compressing everything if that shows a real saving.

SRC_DIR=/Xnfs/igfldb/grangel/fmr1_fast5_raw
DST_DIR=/scratch/Lake/grangel/fmr1_polya_splicing/fast5

mkdir -p "$DST_DIR"

echo "Copying $SRC_DIR -> $DST_DIR ..."
rsync -a --info=progress2 "$SRC_DIR"/ "$DST_DIR"/

src_size=$(du -sb "$SRC_DIR" | cut -f1)
dst_size=$(du -sb "$DST_DIR" | cut -f1)

if [ "$src_size" -eq "$dst_size" ]; then
    echo "Verified: sizes match ($src_size bytes). Safe to delete the igfldb copy."
    rm -rf "$SRC_DIR"
    echo "Deleted $SRC_DIR. fast5 data now lives only at $DST_DIR."
    echo "Reminder: \$SCRATCH/Lake has NO backup/snapshots -- this is working"
    echo "space, not archival. Keep Kiran's/Japan's copy as the backup of record."
else
    echo "SIZE MISMATCH (src=$src_size dst=$dst_size bytes)."
    echo "NOT deleting the igfldb copy. Check $DST_DIR manually before retrying."
fi

df -h /Xnfs/igfldb /scratch/Lake
