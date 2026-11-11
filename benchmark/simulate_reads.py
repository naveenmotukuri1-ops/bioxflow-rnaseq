#!/usr/bin/env python3
"""
Simulate stranded paired-end RNA-seq reads with a known ground truth.

Reads are sampled from annotated transcripts (spliced, strand-aware) of any genome
FASTA + GTF. Fragment counts follow a negative binomial distribution, a chosen
fraction of genes is differentially expressed between two conditions with known
log2 fold changes, and reads carry random substitution errors. The library layout
is reverse-stranded (dUTP / Illumina TruSeq Stranded): read 1 is antisense.

Outputs: FASTQ files, samplesheet.csv, truth_gene_counts.tsv, truth_de.tsv

Requires only numpy.
"""
import argparse
import gzip
import os
import re
from collections import defaultdict

import numpy as np

COMP = np.zeros(256, dtype=np.uint8)
for a, b in zip(b"ACGTN", b"TGCAN"):
    COMP[a] = b
BASES = np.frombuffer(b"ACGT", dtype=np.uint8)


def open_text(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def read_fasta(path):
    seqs, name, chunks = {}, None, []
    with open_text(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(chunks).upper()
                name, chunks = line[1:].split()[0], []
            else:
                chunks.append(line.strip())
    if name is not None:
        seqs[name] = "".join(chunks).upper()
    return seqs


def read_transcripts(gtf, genome):
    """Return {transcript_id: (gene_id, uint8 sequence in transcript orientation)}."""
    exons = defaultdict(list)
    info = {}
    tid_re, gid_re = re.compile(r'transcript_id "([^"]+)"'), re.compile(r'gene_id "([^"]+)"')
    with open_text(gtf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "exon" or f[0] not in genome:
                continue
            tid, gid = tid_re.search(f[8]), gid_re.search(f[8])
            if not tid or not gid:
                continue
            exons[tid.group(1)].append((int(f[3]), int(f[4])))
            info[tid.group(1)] = (gid.group(1), f[0], f[6])
    out = {}
    for tid, ex in exons.items():
        gid, chrom, strand = info[tid]
        seq = "".join(genome[chrom][s - 1:e] for s, e in sorted(ex))
        arr = np.frombuffer(seq.encode(), dtype=np.uint8)
        arr = np.where(np.isin(arr, BASES), arr, ord("N")).astype(np.uint8)
        if strand == "-":
            arr = COMP[arr][::-1]
        out[tid] = (gid, np.ascontiguousarray(arr))
    return out


def write_fastq(path, reads, names, qual):
    with gzip.open(path, "wb", compresslevel=4) as fh:
        step = 50000
        for i in range(0, len(reads), step):
            block = [b"@%s\n%s\n+\n%s\n" % (n, r.tobytes(), qual) for n, r in zip(names[i:i + step], reads[i:i + step])]
            fh.write(b"".join(block))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--replicates", type=int, default=3)
    ap.add_argument("--depth", type=int, default=500000, help="mean read pairs per sample")
    ap.add_argument("--read-length", type=int, default=100)
    ap.add_argument("--fragment-mean", type=float, default=250)
    ap.add_argument("--fragment-sd", type=float, default=30)
    ap.add_argument("--de-fraction", type=float, default=0.15)
    ap.add_argument("--dispersion", type=float, default=0.05)
    ap.add_argument("--error-rate", type=float, default=0.005)
    ap.add_argument("--seed", type=int, default=1)
    args = ap.parse_args()

    rng = np.random.default_rng(args.seed)
    rl = args.read_length
    os.makedirs(args.outdir, exist_ok=True)

    genome = read_fasta(args.fasta)
    transcripts = {t: v for t, v in read_transcripts(args.gtf, genome).items() if len(v[1]) >= max(rl, 200)}
    by_gene = defaultdict(list)
    for tid, (gid, _) in transcripts.items():
        by_gene[gid].append(tid)
    genes = sorted(by_gene)
    n = len(genes)

    # expression model
    weight = rng.lognormal(mean=0.0, sigma=1.4, size=n)
    mu_control = args.depth * weight / weight.sum()
    is_de = np.zeros(n, dtype=bool)
    is_de[rng.choice(n, size=int(round(args.de_fraction * n)), replace=False)] = True
    lfc = np.where(is_de, rng.uniform(1.0, 3.0, size=n) * rng.choice([-1.0, 1.0], size=n), 0.0)
    mu = {"control": mu_control, "treated": mu_control * 2.0 ** lfc}
    isoform_props = {g: rng.dirichlet(np.ones(len(by_gene[g]))) for g in genes}

    qual = b"I" * rl
    truth = {}
    sheet = ["sample,fastq_1,fastq_2,condition"]
    for cond in ("control", "treated"):
        for rep in range(1, args.replicates + 1):
            sample = f"{cond}_rep{rep}"
            size = rng.uniform(0.8, 1.25)
            r = 1.0 / args.dispersion
            m = mu[cond] * size
            counts = rng.negative_binomial(r, r / (r + m))
            truth[sample] = counts

            r1_blocks, r2_blocks = [], []
            for gi, g in enumerate(genes):
                if counts[gi] == 0:
                    continue
                per_tx = rng.multinomial(counts[gi], isoform_props[g])
                for tid, c in zip(by_gene[g], per_tx):
                    if c == 0:
                        continue
                    seq = transcripts[tid][1]
                    L = len(seq)
                    frag = np.clip(np.rint(rng.normal(args.fragment_mean, args.fragment_sd, size=c)), rl, L).astype(np.int64)
                    start = np.floor(rng.random(c) * (L - frag + 1)).astype(np.int64)
                    idx = np.arange(rl)
                    sense = seq[start[:, None] + idx]                       # 5' end of the fragment
                    other = seq[(start + frag - rl)[:, None] + idx]         # 3' end of the fragment
                    r1_blocks.append(COMP[other][:, ::-1])                  # read 1: antisense
                    r2_blocks.append(sense)                                 # read 2: sense
            r1 = np.concatenate(r1_blocks)
            r2 = np.concatenate(r2_blocks)
            for reads in (r1, r2):
                err = rng.random(reads.shape) < args.error_rate
                reads[err] = BASES[rng.integers(0, 4, size=int(err.sum()))]
            order = rng.permutation(len(r1))   # mappers expect reads in random order
            r1, r2 = r1[order], r2[order]
            names = [b"%s.%d" % (sample.encode(), i + 1) for i in range(len(r1))]
            f1, f2 = f"{sample}_R1.fastq.gz", f"{sample}_R2.fastq.gz"
            write_fastq(os.path.join(args.outdir, f1), r1, names, qual)
            write_fastq(os.path.join(args.outdir, f2), r2, names, qual)
            sheet.append(f"{sample},{f1},{f2},{cond}")
            print(f"{sample}: {len(r1)} read pairs")

    with open(os.path.join(args.outdir, "samplesheet.csv"), "w") as fh:
        fh.write("\n".join(sheet) + "\n")
    samples = list(truth)
    with open(os.path.join(args.outdir, "truth_gene_counts.tsv"), "w") as fh:
        fh.write("gene_id\t" + "\t".join(samples) + "\n")
        for gi, g in enumerate(genes):
            fh.write(g + "\t" + "\t".join(str(truth[s][gi]) for s in samples) + "\n")
    with open(os.path.join(args.outdir, "truth_de.tsv"), "w") as fh:
        fh.write("gene_id\ttrue_log2fc\tis_de\texpected_control_count\n")
        for gi, g in enumerate(genes):
            fh.write(f"{g}\t{lfc[gi]:.4f}\t{int(is_de[gi])}\t{mu_control[gi]:.2f}\n")
    print(f"genes simulated: {n}; truly differential: {int(is_de.sum())}")


if __name__ == "__main__":
    main()
