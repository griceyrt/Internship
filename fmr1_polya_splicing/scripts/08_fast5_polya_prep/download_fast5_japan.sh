#!/bin/bash
#SBATCH --job-name=fast5_download
#SBATCH --partition=Lake
#SBATCH --time=7-00:00:00
#SBATCH --mem=4G
#SBATCH --output=%j_fast5_download.log

# Downloads the ~3.4TB fast5 raw-signal dataset from Kiran's Japan
# collaborator (kero.hgc.jp).
#
# LANDING ZONE: /Xnfs/igfldb, a shared lab-group NFS mount that gets
# periodically cleared -- treat as temporary staging, not a final
# destination. Once the download is verified complete, run
# migrate_fast5_from_igfldb.sh to move it to its permanent home on
# /scratch/Lake.
#
# SERVER QUIRKS (why this script looks the way it does): kero.hgc.jp
# doesn't redirect a no-trailing-slash directory URL to its slash-
# terminated form the way well-behaved servers do, so wget must be given
# each directory explicitly with a trailing slash -- otherwise it tries to
# save the listing as a FILE at a path where a real directory already
# exists, fails with "Is a directory", and silently never recurses into
# it. This applies at every directory level (sample folder, then
# fast5_pass/fast5_fail within it), hence the nested loop below rather
# than one recursive crawl from the top. The server also doesn't handle
# HTTP Range requests cleanly (partial-file resume gets stuck retrying
# indefinitely) and never sends Last-Modified headers, so -N/-c (wget's
# normal resume flags) don't work here -- using -nc (no-clobber) instead,
# which does a clean full download for anything not already present.
#
# Because -nc trusts "file already exists" without checking completeness,
# the pre-cleanup step below deletes any locally-truncated .fast5 file
# (below 50% of its own folder's median size) before each run, so
# partial files from an earlier interrupted job get a fresh download
# instead of being silently skipped forever. It also clears cached
# index.html listing pages, since -nc would otherwise skip re-fetching
# those too and the crawl might never re-discover new files.
#
# VERIFY EACH RUN ACTUALLY PROGRESSED: check `du -sh` on $DEST_DIR before
# and after -- it should grow. If a run finishes in well under a minute
# moving no new bytes, something upstream of this logic has changed and
# needs a fresh look at the log for "Is a directory" or "Retrying" patterns.

DEST_DIR=/Xnfs/igfldb/grangel/fmr1_fast5_raw
mkdir -p "$DEST_DIR"

SAMPLES=(210126_DrLies_WT1_CTX 210126_DrLies_WT2_CTX 210126_DrLies_WT3_CTX
         210126_DrLies_KO1_CTX 210126_DrLies_KO2_CTX 210126_DrLies_KO3_CTX)
SUBDIRS=(fast5_pass fast5_fail)

echo "=== Pre-cleanup: scanning for truncated files (< 50% of their own folder's median size) ==="
for SAMPLE in "${SAMPLES[@]}"; do
    for SUBDIR in "${SUBDIRS[@]}"; do
        FOLDER="$DEST_DIR/$SAMPLE/$SUBDIR"
        if [ -d "$FOLDER" ]; then
            python3 - "$FOLDER" << 'PYEOF'
import os, sys, statistics
folder = sys.argv[1]
files = [f for f in os.listdir(folder) if f.endswith(".fast5")]
if len(files) < 2:
    sys.exit(0)  # not enough files to compute a meaningful median
sizes = {f: os.path.getsize(os.path.join(folder, f)) for f in files}
median = statistics.median(sizes.values())
threshold = median * 0.5
deleted_count = 0
deleted_bytes = 0
for f, size in sizes.items():
    if size < threshold:
        path = os.path.join(folder, f)
        print(f"  DELETING (truncated): {path} ({size} bytes, folder median {int(median)} bytes)")
        os.remove(path)
        deleted_count += 1
        deleted_bytes += size
if deleted_count:
    print(f"  {folder}: deleted {deleted_count} truncated files ({deleted_bytes/1e9:.2f} GB) for fresh re-download")
PYEOF
        fi
    done
done
echo "=== Pre-cleanup done ==="
echo ""

echo "=== Clearing cached index.html listing pages (forces fresh directory listing each run) ==="
find "$DEST_DIR" -name "index.html*" -delete
echo "=== Done ==="
echo ""

# Credentials came from Kiran's forwarded email (one-time, kero.hgc.jp).
# The 3.4TB download is complete and migrated (see 08_fast5_polya_prep/
# migrate_fast5_from_igfldb.sh) -- this script does not need to run again.
# Redacted below since this file is now going into git / the lab share
# folder. If this download ever needs to be repeated, get fresh
# credentials from Kiran and set them as environment variables instead of
# hardcoding them here:
#   export KERO_USER='...'
#   export KERO_PASS='...'
KERO_USER="${KERO_USER:?set KERO_USER before running}"
KERO_PASS="${KERO_PASS:?set KERO_PASS before running}"

for SAMPLE in "${SAMPLES[@]}"; do
    for SUBDIR in "${SUBDIRS[@]}"; do
        echo "=== Downloading $SAMPLE/$SUBDIR ==="
        wget -nH --cut-dirs=4 -r --no-parent -nc \
             -P "$DEST_DIR" \
             --http-user="$KERO_USER" --http-password="$KERO_PASS" \
             "https://kero.hgc.jp/PAGS/Joint/20K_7/download/${SAMPLE}/${SUBDIR}/"
    done
    echo "=== $SAMPLE done, running total: ==="
    du -sh "$DEST_DIR"
done

echo "=== Download attempt finished ==="
echo "Actual size downloaded:"
du -sh "$DEST_DIR"
echo "Actual file count:"
find "$DEST_DIR" -type f | wc -l
echo "Compare against the ~3.4TB stated in Japan's email -- if noticeably"
echo "short, re-submit this same script (sbatch download_fast5_japan.sh)."
