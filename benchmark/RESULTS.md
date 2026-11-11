# Benchmark results

All numbers on this page were measured by running the pipeline; nothing is estimated.
Reads were produced by [`simulate_reads.py`](simulate_reads.py), so the true count of every
gene and the true set of differentially expressed genes are known, and scored with
[`evaluate.py`](evaluate.py).

## What this benchmark shows, and what it does not

* It shows that the pipeline recovers known expression levels and known fold changes from
  FASTQ files, end to end, and that the differential expression calls are well calibrated
  on data that follow a negative binomial model.
* It is a simulation. Real libraries have biases that the simulator does not model
  (degradation, GC and positional bias, intronic and intergenic reads, contamination), and
  the yeast annotation has almost no alternative splicing, so isoform ambiguity is barely
  tested. Performance on real human or mouse data has **not** yet been measured against an
  independent reference such as the SEQC/MAQC qPCR sets. Treat these numbers as a
  correctness check, not as a claim about every dataset.

## Bundled test dataset (`-profile test`)

*S. cerevisiae* chromosome I, 117 expressed genes, 2 conditions x 3 replicates, about
10,000 read pairs of 2 x 75 bp per sample, 23 genes truly differential (|log2FC| between 1
and 3). Run on macOS (Apple Silicon) with the `conda` profile.

## Quantification accuracy (estimated vs true gene counts)

| Sample | Method | Genes | Spearman | Pearson (log2) | Mean abs. rel. error |
|---|---|---|---|---|---|
| control_rep1 | Salmon | 117 | 0.9984 | 0.9870 | 0.86% |
| control_rep2 | Salmon | 117 | 0.9979 | 0.9852 | 0.77% |
| control_rep3 | Salmon | 117 | 0.9977 | 0.9841 | 1.16% |
| treated_rep1 | Salmon | 117 | 0.9953 | 0.9701 | 2.59% |
| treated_rep2 | Salmon | 117 | 0.9979 | 0.9761 | 2.68% |
| treated_rep3 | Salmon | 117 | 0.9988 | 0.9857 | 2.01% |
| control_rep1 | STAR | 117 | 0.9428 | 0.8919 | 5.51% |
| control_rep2 | STAR | 117 | 0.9560 | 0.9022 | 5.53% |
| control_rep3 | STAR | 117 | 0.9299 | 0.8669 | 5.98% |
| treated_rep1 | STAR | 117 | 0.9207 | 0.8904 | 6.82% |
| treated_rep2 | STAR | 117 | 0.9297 | 0.8946 | 6.52% |
| treated_rep3 | STAR | 117 | 0.9364 | 0.9092 | 6.44% |

## Differential expression accuracy

| Metric | Value |
|---|---|
| Genes simulated | 117 |
| Genes tested after low-count filtering | 101 |
| Truly differential genes (all / tested) | 23 / 19 |
| Genes called at padj < 0.05 | 18 |
| True positives | 18 |
| False positives | 0 |
| Sensitivity (of all truly differential genes) | 78.3% |
| Observed false discovery rate | 0.0% |
| True positives with the wrong direction | 1 |
| log2FC Pearson r, estimate vs truth | 0.8637 |
| log2FC Pearson r, shrunken estimate vs truth | 0.8822 |
| log2FC RMSE, estimate | 0.492 |
| log2FC RMSE, shrunken estimate | 0.441 |

Notes:

* Salmon counts track the truth almost exactly. STAR gene counts are lower for a few genes
  that overlap another gene, because STAR sets aside reads it cannot assign to a single
  gene. This is why Salmon is the primary quantification and STAR is the cross-check.
* One gene was called significant with the opposite sign to the truth. At this tiny depth
  Salmon's automatic library-type detection does not trigger and it quantifies without
  strand information, so a gene that overlaps a differential gene on the opposite strand
  can pick up its neighbour's reads.

## Full yeast genome

Reproduce with `bash benchmark/run_benchmark.sh <profile>`. Whole *S. cerevisiae* R64-1-1
genome (Ensembl release 112), 6,531 expressed genes, 2 conditions x 3 replicates, about
440,000 to 690,000 read pairs of 2 x 100 bp per sample, 980 genes truly differential
(|log2FC| between 1 and 3). Run on macOS (Apple Silicon, 8 GB RAM) with the `conda` profile.

## Quantification accuracy (estimated vs true gene counts)

| Sample | Method | Genes | Spearman | Pearson (log2) | Mean abs. rel. error |
|---|---|---|---|---|---|
| control_rep1 | Salmon | 6531 | 0.9768 | 0.9612 | 3.52% |
| control_rep2 | Salmon | 6531 | 0.9752 | 0.9584 | 4.37% |
| control_rep3 | Salmon | 6531 | 0.9741 | 0.9605 | 3.40% |
| treated_rep1 | Salmon | 6531 | 0.9757 | 0.9622 | 3.73% |
| treated_rep2 | Salmon | 6531 | 0.9760 | 0.9638 | 3.54% |
| treated_rep3 | Salmon | 6531 | 0.9772 | 0.9620 | 4.27% |
| control_rep1 | STAR | 6531 | 0.9225 | 0.8814 | 5.51% |
| control_rep2 | STAR | 6531 | 0.9229 | 0.8847 | 5.39% |
| control_rep3 | STAR | 6531 | 0.9223 | 0.8895 | 5.37% |
| treated_rep1 | STAR | 6531 | 0.9244 | 0.8909 | 5.49% |
| treated_rep2 | STAR | 6531 | 0.9239 | 0.8943 | 5.44% |
| treated_rep3 | STAR | 6531 | 0.9247 | 0.8896 | 5.51% |

## Differential expression accuracy

| Metric | Value |
|---|---|
| Genes simulated | 6531 |
| Genes tested after low-count filtering | 5229 |
| Truly differential genes (all / tested) | 980 / 818 |
| Genes called at padj < 0.05 | 752 |
| True positives | 728 |
| False positives | 24 |
| Sensitivity (of all truly differential genes) | 74.3% |
| Observed false discovery rate | 3.2% |
| True positives with the wrong direction | 0 |
| log2FC Pearson r, estimate vs truth | 0.8980 |
| log2FC Pearson r, shrunken estimate vs truth | 0.9202 |
| log2FC RMSE, estimate | 0.398 |
| log2FC RMSE, shrunken estimate | 0.327 |

Salmon-vs-STAR concordance (the pipeline's own quality gate) was 0.916 to 0.922 Spearman
for all six samples, both methods inferred a reverse-stranded library, and every sample
passed.

Notes:

* The observed false discovery rate (3.2%) is below the 5% that was asked for, so the
  calls are not over-confident on this data.
* Sensitivity is counted over all 980 truly differential genes, including 162 that were
  expressed too weakly at this depth to be tested at all. Among the 818 tested, 89% were found.
* No gene was called in the wrong direction. Shrunken fold changes were closer to the truth
  than the raw estimates (RMSE 0.33 vs 0.40).
* STAR counts agree less well with the truth than Salmon counts. The yeast annotation has
  many overlapping genes, and STAR sets aside reads that it cannot assign to one gene.
