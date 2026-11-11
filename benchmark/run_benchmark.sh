#!/usr/bin/env bash
# Ground-truth benchmark for bioxflow/rnaseq on the full S. cerevisiae genome.
#
#   bash benchmark/run_benchmark.sh [profile] [workdir]
#
# 1. downloads the Ensembl R64-1-1 yeast genome and annotation (release 112)
# 2. simulates 2 conditions x 3 replicates of stranded paired-end reads with known
#    expression levels and known differentially expressed genes
# 3. runs the pipeline
# 4. scores the results against the truth
#
# Needs Nextflow, python3 with numpy, and conda or docker. About 2 GB of disk.
set -euo pipefail

PROFILE="${1:-docker}"
WORK="${2:-benchmark_run}"
HERE="$(cd "$(dirname "$0")/.." && pwd)"
ENSEMBL="https://ftp.ensembl.org/pub/release-112"

mkdir -p "$WORK/ref" && cd "$WORK"
[ -s ref/genome.fa.gz ] || curl -fL -o ref/genome.fa.gz "$ENSEMBL/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz"
[ -s ref/genes.gtf.gz ] || curl -fL -o ref/genes.gtf.gz "$ENSEMBL/gtf/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.112.gtf.gz"

python3 "$HERE/benchmark/simulate_reads.py" --fasta ref/genome.fa.gz --gtf ref/genes.gtf.gz \
    --outdir sim --depth 500000 --read-length 100 --de-fraction 0.15 --seed 1

nextflow run "$HERE" -profile "$PROFILE" \
    --input sim/samplesheet.csv --fasta ref/genome.fa.gz --gtf ref/genes.gtf.gz \
    --outdir results --max_cpus 4 --max_memory 6.GB

python3 "$HERE/benchmark/evaluate.py" --truth sim --results results | tee benchmark_results.md
