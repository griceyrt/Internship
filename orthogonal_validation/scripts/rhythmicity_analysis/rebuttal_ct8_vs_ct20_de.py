#!/usr/bin/env python3
"""
Rebuttal point (d) -- WT CT8 vs CT20 differential expression. Reproduces
fig4_1_2024_12_04.R's method for WT-CT8 vs WT-CT20 (transcript-level, global
expressed-transcript filter across all 16 samples, plain DESeq2 Wald test).
Author: Gricey

Uses pydeseq2 since this sandbox has no R -- rerun rebuttal_ct8_vs_ct20_de.R
on the cluster for exact parity with DESeq2-in-R.

Run from: orthogonal_validation/
    python3 scripts/rebuttal_ct8_vs_ct20_de.py
"""
import glob
import pickle
import re
from pathlib import Path

import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

PROJECT_ROOT = Path(__file__).resolve().parent.parent
SALMON_DIR = PROJECT_ROOT / "data" / "salmon_2024-09-25"
GTF_PATH = (
    PROJECT_ROOT.parent / "boundary_analysis" / "data" / "transcriptome_productivity.gtf"
)
OUT_DIR = PROJECT_ROOT / "results" / "rhythmicity_analysis" / "rebuttal_ct8_vs_ct20"
OUT_DIR.mkdir(parents=True, exist_ok=True)

COUNT_THRESHOLD = 5
MIN_SAMPLES = 2
FDR_CUTOFF = 0.05
LOG2FC_CUTOFF = 1.5

CT8_SAMPLES = ["CT8_1", "CT8_2"]
CT20_SAMPLES = ["CT20_1", "CT20_2"]


def build_tx2gene(gtf_path: Path) -> dict:
    tx2gene = {}
    with open(gtf_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "transcript":
                continue
            attrs = fields[8]
            g = re.search(r'gene_id "([^"]+)"', attrs)
            t = re.search(r'transcript_id "([^"]+)"', attrs)
            if g and t:
                tx2gene[t.group(1)] = g.group(1)
    return tx2gene


def load_all_16_sample_counts() -> pd.DataFrame:
    sf_files = sorted(glob.glob(str(SALMON_DIR / "quant_*" / "quant_*.sf")))
    if len(sf_files) != 16:
        raise RuntimeError(f"expected 16 salmon quant.sf files, found {len(sf_files)}")
    counts = pd.DataFrame()
    for f in sf_files:
        sample = re.search(r"quant_([A-Za-z0-9]+_\d)\.sf$", f).group(1)
        df = pd.read_csv(f, sep="\t", index_col="Name")
        counts[sample] = df["NumReads"]
    return counts


def load_sample_counts(samples: list) -> pd.DataFrame:
    counts = pd.DataFrame()
    for s in samples:
        path = SALMON_DIR / f"quant_{s}" / f"quant_{s}.sf"
        df = pd.read_csv(path, sep="\t", index_col="Name")
        counts[s] = df["NumReads"]
    return counts


def main():
    print("Building tx2gene map from", GTF_PATH)
    tx2gene = build_tx2gene(GTF_PATH)
    print(f"  {len(tx2gene)} transcript -> gene mappings")

    print("\nComputing global expressed-transcript filter across all 16 samples...")
    counts_all = load_all_16_sample_counts()
    keep = (counts_all > COUNT_THRESHOLD).sum(axis=1) >= MIN_SAMPLES
    expressed_tx = counts_all.index[keep]
    print(f"  {len(expressed_tx)} / {len(counts_all)} transcripts pass the filter")
    pd.Series(expressed_tx, name="transcript_id").to_csv(
        OUT_DIR / "expressed_tx_list.csv", index=False
    )

    print("\nLoading CT8 + CT20 transcript counts...")
    tx_counts = load_sample_counts(CT8_SAMPLES + CT20_SAMPLES)
    tx_counts_f = tx_counts.loc[tx_counts.index.isin(expressed_tx)].round().astype(int)
    print(f"  {tx_counts_f.shape[0]} transcripts x {tx_counts_f.shape[1]} samples after filter")

    metadata = pd.DataFrame(
        {"condition": ["CT8"] * len(CT8_SAMPLES) + ["CT20"] * len(CT20_SAMPLES)},
        index=CT8_SAMPLES + CT20_SAMPLES,
    )

    print("\nRunning pydeseq2 (Wald test, condition CT20 vs CT8)...")
    dds = DeseqDataSet(counts=tx_counts_f.T, metadata=metadata, design="~condition", refit_cooks=True)
    dds.deseq2()
    stat = DeseqStats(dds, contrast=["condition", "CT20", "CT8"])
    stat.summary()
    res = stat.results_df
    res["gene_id"] = res.index.map(tx2gene)
    res.to_csv(OUT_DIR / "deseq2_tx_level_ct8_vs_ct20.csv")

    res_valid = res.dropna(subset=["padj"])
    sig = res_valid[(res_valid["log2FoldChange"].abs() > LOG2FC_CUTOFF) & (res_valid["padj"] < FDR_CUTOFF)]
    print(f"\nTranscripts with valid padj: {len(res_valid)}")
    print(f"Significant (|log2FC|>{LOG2FC_CUTOFF} & padj<{FDR_CUTOFF}): {len(sig)} transcripts, "
          f"{sig['gene_id'].nunique()} unique genes")
    print(f"\nSaved: {OUT_DIR / 'deseq2_tx_level_ct8_vs_ct20.csv'}")


if __name__ == "__main__":
    main()
