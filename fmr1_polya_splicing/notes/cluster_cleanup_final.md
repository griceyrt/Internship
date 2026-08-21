# Final cluster cleanup — 2026-08-21 (last internship day)

Run these on `cl5218comp1` (`ssh -A -J grangel@ssh.psmn.ens-lyon.fr grangel@allo-psmn.psmn.ens-lyon.fr` then `ssh -A cl5218comp1`).

## 1. Confirm the scripts/ reorg actually landed, then remove the backup

```bash
find ~/fmr1_polya_splicing/scripts -type f | sort
# should show 34 files across 10_ numbered subfolders (00-10), matching
# the local Mac copy exactly. The nanopolish_polya.sh job already ran
# successfully from scripts/10_polya_length/, which is itself proof the
# reorg + path patching worked -- safe to delete the backup:
rm -rf ~/fmr1_polya_splicing/scripts_old_flat_backup
```

## 2. Secure the fast5 download credentials

```bash
chmod 600 ~/fmr1_polya_splicing/scripts/08_fast5_polya_prep/download_fast5_japan.sh
```

## 3. Clean dead empty result stubs (same ones cleaned locally already)

```bash
find ~/fmr1_polya_splicing/results -maxdepth 1 -type d -empty
# rmdir whatever that lists -- these were leftover from an early
# folder-structure template, superseded by the real dated output folders
```

## 4. Decide on the raw .sra files (36G, your call — not forced)

The derived fastq.gz files were already fully verified against these
back on 2026-07-27 (read counts, length distributions all checked out),
so the .sra archives are redundant. Re-downloading has a small
re-acquisition cost if ever needed again.

```bash
du -sh ~/fmr1_polya_splicing/raw_data
# if you're comfortable deleting:
# rm -rf ~/fmr1_polya_splicing/raw_data/*
```

## 5. Check $HOME usage is actually back under control

```bash
df -h /home/grangel
```
Should be well below the 76%-full state from earlier in the project,
now that the BAMs/fastq/fast5 all live on `/scratch/Lake` instead.

## 6. Check the PSMN gateway for accumulated clutter one more time

Files single-hop-transferred straight to the gateway (rather than through
`allo-psmn`) silently land in a ~5MB stub invisible from any real node —
this has happened multiple times across this project and project 4.

```bash
ssh grangel@ssh.psmn.ens-lyon.fr "find ~ -type f; du -sh ~"
# delete anything found under a project subfolder (real copies live
# elsewhere by definition) -- but leave .ssh/ and .bash_history alone
```

## 7. Final sanity check — confirm fast5 is where it should be

```bash
du -sh /scratch/Lake/grangel/fmr1_polya_splicing/fast5
# should show ~3.4T (the du -sh display rounds down from the true
# 3.70TB, this is expected -- see main README)
ls /Xnfs/igfldb/grangel/ 2>&1
# should show "No such file or directory" -- confirms nothing was
# left behind on the periodically-wiped staging mount
```

## 8. Leave a note for whoever inherits this

Once the above is done, the cluster-side `~/fmr1_polya_splicing/` should
mirror the local Mac copy exactly (`scripts/` only — `plot_scripts/`,
`figures/`, `notes/`, `papers/` never needed to be on the cluster). Point
whoever continues this to the main `README.md` in the project root.
