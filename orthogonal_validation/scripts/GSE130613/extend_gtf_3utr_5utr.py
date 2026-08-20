"""
Extends transcript boundaries in transcriptome_productivity.gtf to match
Asher's own MARS-seq processing (SI appendix, PNAS 2020): 1000bp toward the
5' edge, 100bp toward the 3' end, to better capture 3'-biased MARS-seq reads.
Author: Gricey

Extends the TERMINAL EXON of each transcript (so gffread's FASTA extraction
picks up the extra flanking bases) plus the parent transcript/gene lines.
CDS/UTR/codon lines untouched. Bounds capped at [1, chromosome_length] via
the genome .fai index.

Two-pass streaming (keeps memory low on the 900MB file): pass 1 records each
transcript's strand + min/max exon coordinates; pass 2 rewrites every line,
extending boundaries that match the transcript's terminal coordinate.

Usage:
    python3 extend_gtf_3utr_5utr.py <input.gtf> <genome.fa.fai> <output.gtf>
"""
import sys
import re

EXT_5PRIME = 1000
EXT_3PRIME = 100

TX_ID_RE = re.compile(r'transcript_id "([^"]+)"')


def get_chrom_lengths(fai_path):
    lengths = {}
    with open(fai_path) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            lengths[parts[0]] = int(parts[1])
    return lengths


def pass1_transcript_bounds(gtf_path):
    """Returns dict: transcript_id -> {strand, chrom, min_start, max_end}"""
    bounds = {}
    with open(gtf_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "exon":
                continue
            chrom, _, _, start, end, _, strand, _, attrs = fields
            m = TX_ID_RE.search(attrs)
            if not m:
                continue
            txid = m.group(1)
            start, end = int(start), int(end)
            if txid not in bounds:
                bounds[txid] = {"strand": strand, "chrom": chrom,
                                 "min_start": start, "max_end": end}
            else:
                b = bounds[txid]
                if start < b["min_start"]:
                    b["min_start"] = start
                if end > b["max_end"]:
                    b["max_end"] = end
    return bounds


def extended_coords(strand, is_min_start, is_max_end, start, end, chrom_len):
    """Given an exon's start/end, and whether it's the transcript's terminal
    min_start / max_end exon, return the (possibly extended) start/end."""
    new_start, new_end = start, end
    if strand == "+":
        # 5' edge = min_start (extend backward/left by 1000)
        if is_min_start:
            new_start = max(1, start - EXT_5PRIME)
        # 3' edge = max_end (extend forward/right by 100)
        if is_max_end:
            new_end = min(chrom_len, end + EXT_3PRIME)
    elif strand == "-":
        # 5' edge = max_end (extend forward/right by 1000, since - strand
        # transcription runs from high to low genomic coordinate)
        if is_max_end:
            new_end = min(chrom_len, end + EXT_5PRIME)
        # 3' edge = min_start (extend backward/left by 100)
        if is_min_start:
            new_start = max(1, start - EXT_3PRIME)
    return new_start, new_end


def main():
    if len(sys.argv) != 4:
        print(__doc__)
        sys.exit(1)
    gtf_in, fai_path, gtf_out = sys.argv[1], sys.argv[2], sys.argv[3]

    print("Pass 1: recording per-transcript strand + min/max exon bounds...")
    tx_bounds = pass1_transcript_bounds(gtf_in)
    print(f"  {len(tx_bounds)} transcripts found")

    chrom_lengths = get_chrom_lengths(fai_path)

    print("Pass 2: rewriting GTF with extended terminal exon/transcript/gene bounds...")
    n_exon_extended = 0
    n_lines = 0
    with open(gtf_in) as fin, open(gtf_out, "w") as fout:
        for line in fin:
            n_lines += 1
            if line.startswith("#"):
                fout.write(line)
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                fout.write(line)
                continue
            chrom, source, feature, start, end, score, strand, frame, attrs = fields
            start, end = int(start), int(end)

            if feature in ("exon", "transcript", "gene"):
                m = TX_ID_RE.search(attrs)
                # gene lines have no transcript_id and are left untouched --
                # gffread's FASTA extraction doesn't need them
                if feature in ("exon", "transcript") and m:
                    txid = m.group(1)
                    b = tx_bounds.get(txid)
                    if b:
                        chrom_len = chrom_lengths.get(chrom, 10**12)
                        is_min = (start == b["min_start"])
                        is_max = (end == b["max_end"])
                        new_start, new_end = extended_coords(
                            b["strand"], is_min, is_max, start, end, chrom_len)
                        if (new_start, new_end) != (start, end):
                            n_exon_extended += 1
                        start, end = new_start, new_end

            fout.write(f"{chrom}\t{source}\t{feature}\t{start}\t{end}\t{score}\t{strand}\t{frame}\t{attrs}\n")

    print(f"Done. {n_lines} lines processed, {n_exon_extended} exon/transcript boundary lines extended.")
    print(f"Output: {gtf_out}")


if __name__ == "__main__":
    main()
