#!/usr/bin/env python3
"""
Compare bioxflow/rnaseq results with the ground truth written by simulate_reads.py.

  evaluate.py --truth <simulated dir> --results <pipeline outdir> [--contrast treated_vs_control]

Reports, as Markdown:
  * agreement between estimated and true gene counts (Salmon and STAR, per sample)
  * agreement between estimated and true log2 fold changes
  * sensitivity and false discovery rate of the differential expression calls
"""
import argparse
import csv
import os

import numpy as np


def table(path):
    with open(path) as fh:
        rows = list(csv.reader(fh, delimiter="\t"))
    return rows[0], rows[1:]


def rank(x):
    order = np.argsort(x, kind="mergesort")
    r = np.empty(len(x))
    r[order] = np.arange(len(x))
    # average ranks for ties
    xs = x[order]
    i = 0
    while i < len(xs):
        j = i
        while j + 1 < len(xs) and xs[j + 1] == xs[i]:
            j += 1
        r[order[i:j + 1]] = (i + j) / 2.0
        i = j + 1
    return r


def spearman(a, b):
    return float(np.corrcoef(rank(a), rank(b))[0, 1])


def pearson(a, b):
    return float(np.corrcoef(a, b)[0, 1])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--truth", required=True)
    ap.add_argument("--results", required=True)
    ap.add_argument("--contrast", default="treated_vs_control")
    ap.add_argument("--alpha", type=float, default=0.05)
    args = ap.parse_args()

    head, rows = table(os.path.join(args.truth, "truth_gene_counts.tsv"))
    samples = head[1:]
    truth = {r[0]: np.array(r[1:], dtype=float) for r in rows}
    genes = sorted(truth)

    print("## Quantification accuracy (estimated vs true gene counts)\n")
    print("| Sample | Method | Genes | Spearman | Pearson (log2) | Mean abs. rel. error |")
    print("|---|---|---|---|---|---|")
    for method, fname in (("Salmon", "salmon.gene_counts.tsv"), ("STAR", "star.gene_counts.tsv")):
        path = os.path.join(args.results, "quantification", fname)
        if not os.path.exists(path):
            continue
        h, rs = table(path)
        est = {r[0]: dict(zip(h[2:], map(float, r[2:]))) for r in rs}
        for si, s in enumerate(samples):
            t = np.array([truth[g][si] for g in genes])
            e = np.array([est.get(g, {}).get(s, 0.0) for g in genes])
            expressed = t >= 10
            mard = np.mean(np.abs(e[expressed] - t[expressed]) / t[expressed])
            print(f"| {s} | {method} | {len(genes)} | {spearman(t, e):.4f} | "
                  f"{pearson(np.log2(t + 1), np.log2(e + 1)):.4f} | {100 * mard:.2f}% |")

    de_path = os.path.join(args.results, "differential_expression", f"{args.contrast}.deseq2_results.tsv")
    if not os.path.exists(de_path):
        print(f"\nNo differential expression table found at {de_path}")
        return
    h, rs = table(os.path.join(args.truth, "truth_de.tsv"))
    true_lfc = {r[0]: float(r[1]) for r in rs}
    true_de = {r[0] for r in rs if r[2] == "1"}
    h, rs = table(de_path)
    col = {name: i for i, name in enumerate(h)}
    res = {r[0]: r for r in rs}

    def num(v):
        return float("nan") if v in ("NA", "") else float(v)

    tested = [g for g in genes if g in res and not np.isnan(num(res[g][col["padj"]]))]
    t = np.array([true_lfc[g] for g in tested])
    mle = np.array([num(res[g][col["log2FoldChange"]]) for g in tested])
    shr = np.array([num(res[g][col["log2FoldChange_shrunk"]]) for g in tested])
    called = {g for g in tested if num(res[g][col["padj"]]) < args.alpha}
    tp = len(called & true_de)
    fp = len(called - true_de)
    wrong_sign = sum(1 for g in called & true_de
                     if np.sign(num(res[g][col["log2FoldChange"]])) != np.sign(true_lfc[g]))

    print("\n## Differential expression accuracy\n")
    print("| Metric | Value |")
    print("|---|---|")
    print(f"| Genes simulated | {len(genes)} |")
    print(f"| Genes tested after low-count filtering | {len(tested)} |")
    print(f"| Truly differential genes (all / tested) | {len(true_de)} / {len(true_de & set(tested))} |")
    print(f"| Genes called at padj < {args.alpha} | {len(called)} |")
    print(f"| True positives | {tp} |")
    print(f"| False positives | {fp} |")
    print(f"| Sensitivity (of all truly differential genes) | {100 * tp / max(1, len(true_de)):.1f}% |")
    print(f"| Observed false discovery rate | {100 * fp / max(1, len(called)):.1f}% |")
    print(f"| True positives with the wrong direction | {wrong_sign} |")
    print(f"| log2FC Pearson r, estimate vs truth | {pearson(t, mle):.4f} |")
    print(f"| log2FC Pearson r, shrunken estimate vs truth | {pearson(t, shr):.4f} |")
    print(f"| log2FC RMSE, estimate | {np.sqrt(np.mean((mle - t) ** 2)):.3f} |")
    print(f"| log2FC RMSE, shrunken estimate | {np.sqrt(np.mean((shr - t) ** 2)):.3f} |")


if __name__ == "__main__":
    main()
