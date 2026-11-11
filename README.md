# bioxflow/rnaseq

[![CI](https://github.com/naveenmotukuri1-ops/bioxflow-rnaseq/actions/workflows/ci.yml/badge.svg)](https://github.com/naveenmotukuri1-ops/bioxflow-rnaseq/actions/workflows/ci.yml)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A525.04-18b6a4.svg)](https://www.nextflow.io)
[![License: MIT](https://img.shields.io/badge/license-MIT-0e1b2b.svg)](LICENSE)

A bulk RNA-seq pipeline from [Bioxflow](https://www.bioxflow.uk) that goes from FASTQ files to
differential expression, and checks its own answer on the way.

**What is different about it:** every sample is quantified twice, by two methods that share
no code and use different evidence — Salmon (selective alignment to the transcriptome) and
STAR (spliced alignment to the genome). The pipeline compares the two results gene by gene
and flags any sample where they disagree. It also infers library strandedness with both
methods and reports a mismatch. You get a quality gate based on your own data instead of
having to trust a single tool.

## What it does

```
FASTQ ──► FastQC ──► fastp ──┬──► Salmon (decoy-aware, bias-corrected) ──► tximport ──┐
                             │                                                        ├──► concordance + strandedness check
                             └──► STAR 2-pass (ENCODE settings) ──► gene counts ──────┘
                                                                                      │
                                              DESeq2 (apeglm shrinkage), PCA, plots ◄─┘
                                                                                      │
                                                                         MultiQC report
```

1. **Read QC and trimming** — FastQC, then fastp with automatic adapter detection.
2. **Quantification, route 1** — Salmon selective alignment against a decoy-aware index (the
   whole genome is added as decoy so reads from unannotated loci are not forced onto
   transcripts), with sequence and GC bias correction and automatic library-type detection.
3. **Quantification, route 2** — STAR two-pass alignment with ENCODE long-RNA settings and
   per-gene read counts; sorted, indexed BAM files for genome browsers.
4. **Concordance gate** — per-sample correlation of the two gene-count vectors, the share of
   genes within 2-fold, and a strandedness cross-check. Samples below the threshold are
   marked `REVIEW`.
5. **Differential expression** — tximport gene-level summarisation with transcript-length
   offsets, DESeq2 Wald tests with Benjamini-Hochberg correction, apeglm-shrunken fold
   changes, optional batch term, custom contrasts, PCA, sample distances, volcano and MA plots.
6. **Report** — one MultiQC page with everything, plus plain TSV tables.

All tool versions are pinned, and containers or conda environments are provided for every step.

## Quick start

Install [Nextflow](https://www.nextflow.io) (25.04 or newer) and Docker, Singularity or conda.

Run the bundled test data first (about 5 minutes):

```bash
nextflow run naveenmotukuri1-ops/bioxflow-rnaseq -profile test,docker
```

Then your own data:

```bash
nextflow run naveenmotukuri1-ops/bioxflow-rnaseq -profile docker \
    --input samplesheet.csv --fasta genome.fa.gz --gtf genes.gtf.gz --outdir results
```

```csv
sample,fastq_1,fastq_2,condition
ctrl_1,ctrl_1_R1.fastq.gz,ctrl_1_R2.fastq.gz,control
ctrl_2,ctrl_2_R1.fastq.gz,ctrl_2_R2.fastq.gz,control
ctrl_3,ctrl_3_R1.fastq.gz,ctrl_3_R2.fastq.gz,control
drug_1,drug_1_R1.fastq.gz,drug_1_R2.fastq.gz,treated
drug_2,drug_2_R1.fastq.gz,drug_2_R2.fastq.gz,treated
drug_3,drug_3_R1.fastq.gz,drug_3_R2.fastq.gz,treated
```

On a laptop with about 16 GB of RAM, add `-profile laptop` to run Salmon only (STAR needs
about 32 GB for a human genome).

* [Usage and parameters](docs/usage.md)
* [Output files and how to read them](docs/output.md)

## How accurate is it?

The repository contains a read simulator ([`benchmark/simulate_reads.py`](benchmark/simulate_reads.py))
that produces RNA-seq reads for which the true answer is known, and a scoring script
([`benchmark/evaluate.py`](benchmark/evaluate.py)). Anyone can re-run the benchmark:

```bash
bash benchmark/run_benchmark.sh docker
```

The measured results, and what this benchmark does and does not show, are in
[benchmark/RESULTS.md](benchmark/RESULTS.md).

The bundled test dataset is a small version of the same simulation, and the CI run on
GitHub scores the pipeline against the truth on every change.

## Tools

| Step | Tool | Version |
|---|---|---|
| Read QC | FastQC | 0.12.1 |
| Trimming | fastp | 0.24.0 |
| Transcript extraction | gffread | 0.12.7 |
| Quantification | Salmon | 1.10.3 |
| Alignment and gene counts | STAR | 2.7.11b |
| BAM statistics | SAMtools | 1.21 |
| Gene-level summarisation | tximport | 1.34.0 |
| Differential expression | DESeq2 / apeglm | 1.46.0 / 1.28.0 |
| Report | MultiQC | 1.25.1 |

Please cite the tools you use: see [CITATIONS.md](CITATIONS.md).

## Licence and acknowledgements

bioxflow/rnaseq is released under the [MIT License](LICENSE).

The repository layout and several conventions follow the [nf-core](https://nf-co.re)
community standards (Ewels et al., *Nat Biotechnol* 2020), which are MIT licensed; the
nf-core licence is reproduced in [LICENSE-nf-core](LICENSE-nf-core). This is an independent
pipeline and is not an nf-core pipeline.

## Contributing

Issues and pull requests are welcome. See [www.bioxflow.uk](https://www.bioxflow.uk) for the
community and the other Bioxflow pipelines.
