# PSMN Cluster Workflow — Project 3 (fmr1_polya_splicing)

## 1. Connect (from your notes)

Three hops — each one only lets you go further, none of them run your job:

```
ssh -A grangel@ssh.psmn.ens-lyon.fr      # gateway, do nothing here
ssh -A allo-psmn.psmn.ens-lyon.fr        # entry node, can transfer files / basic commands
ssh -A cl5218comp1                       # login/build node -- this is where you submit jobs from
```

Then activate your environment:

```
source ~/.bashrc
conda activate biotools
```

`prefetch`, `fasterq-dump`, and `pigz` (used by the download script) all confirmed present in this `biotools` env (checked 2026-07-27). Note: this env does *not* have `parallel-fastq-dump` (the tool project 2's `BHARATH_get_sra_dump_fastq_run_salmon.sh` used) — the download script uses `fasterq-dump` + `pigz` instead, which is the modern official sra-tools dumper anyway.

## 2. Get the project folder onto the cluster

**Gotcha (hit 2026-07-27):** the three hops are three separate machines with three separate home directories — `ssh.psmn.ens-lyon.fr` does NOT share storage with `cl5218comp1`. A plain single-hop `rsync`/`scp` to the gateway lands the files somewhere useless. You have to chain through both intermediate hosts in one connection using SSH's `ProxyJump` (`-J`), landing directly on `cl5218comp1`.

Run this from your **Mac terminal**, not the cluster:

```
rsync -av --exclude 'raw_data' --exclude 'temp' \
    -e "ssh -A -J grangel@ssh.psmn.ens-lyon.fr,grangel@allo-psmn.psmn.ens-lyon.fr" \
    ~/Desktop/Internship/fmr1_polya_splicing/ \
    grangel@cl5218comp1:~/fmr1_polya_splicing/
```

(excludes `raw_data`/`temp` since those will be populated directly on the cluster, no need to push empty folders or, later, drag multi-GB files back and forth). If this ever needs to go the other way (pulling results back down to your Mac), same `-J` chain applies, just swap source/destination.

The leftover `/home/grangel/fmr1_polya_splicing` on the gateway host from the earlier single-hop attempt is harmless clutter — safe to ignore or `ssh -A grangel@ssh.psmn.ens-lyon.fr` in and `rm -rf` it later, not urgent.

## 3. Submit the download job

From `cl5218comp1`:

```
cd ~/fmr1_polya_splicing/scripts
sbatch download_gse188840.sh
```

## 4. Check on it

```
squeue -u grangel
sacct -j <jobid>
```

Log lands at `~/fmr1_polya_splicing/results/GSE188840_download_<date>.log` and also as `<jobid>_gse188840_download.log` in the scripts folder (SLURM's own stdout capture).

## 5. The one thing to actually look at when it finishes

Open the log and check the `du -sh` line for each of the 6 samples. This tells us whether SRA gave us fastq-only or fastq+raw signal:

- A few GB per sample → fastq-only. Nanopolish poly(A) estimation (part 2 of Kiran's ask) won't be possible from this SRA download alone — we'd need to ask Kiran/Lies where the original fast5/pod5 files are archived.
- Tens to 100+ GB per sample → raw signal is very likely included, and we can proceed straight into basecalling/Nanopolish.

Either way, this determines the next script, so don't kick off anything downstream until this is checked.
