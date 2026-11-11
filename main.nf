#!/usr/bin/env nextflow
/*
 * bioxflow/rnaseq
 * Bulk RNA-seq: dual-evidence quantification (Salmon + STAR), cross-method
 * concordance checks and DESeq2 differential expression.
 *
 * https://www.bioxflow.co.uk
 */

include { RNASEQ } from './workflows/rnaseq'

def helpMessage() {
    return """
    bioxflow/rnaseq v${workflow.manifest.version}

    Usage:
      nextflow run bioxflow/rnaseq -profile <conda|docker|singularity> \\
          --input samplesheet.csv --fasta genome.fa --gtf genes.gtf --outdir results

    Required:
      --input                 CSV with columns: sample,fastq_1,fastq_2,condition[,batch]
      --fasta                 Genome FASTA (optionally gzipped)
      --gtf                   Gene annotation GTF (optionally gzipped)

    Reference options:
      --transcript_fasta      Transcript FASTA (default: extracted from --fasta and --gtf)
      --gencode               Transcript FASTA / GTF are from GENCODE
      --salmon_index          Pre-built Salmon index directory
      --star_index            Pre-built STAR index directory
      --salmon_decoys         Build a decoy-aware Salmon index [default: true]

    Workflow options:
      --skip_alignment        Salmon only: no STAR alignment and no concordance check
      --skip_trimming         Do not trim reads with fastp
      --skip_fastqc           Do not run FastQC
      --skip_de               Do not run DESeq2
      --save_bam              Publish STAR BAM files [default: true]

    Differential expression:
      --reference_condition   Control condition [default: first condition in the samplesheet]
      --contrasts             CSV with columns: treatment,control [default: every condition vs reference]
      --alpha                 Adjusted p-value cut-off [default: 0.05]
      --lfc_threshold         Minimum absolute shrunken log2 fold change [default: 0]
      --min_count             Low-count filter [default: 10]

    Quality gate:
      --concordance_threshold Minimum Salmon-vs-STAR Spearman correlation per sample [default: 0.9]

    Resources:
      --max_cpus, --max_memory, --max_time    Upper limits for any single process
      -profile laptop                         Salmon-only settings for a 16 GB machine
      -profile test                           Run the bundled test dataset
    """.stripIndent()
}

workflow {
    main:
    if (params.help) {
        log.info(helpMessage())
    }
    else {
        log.info("bioxflow/rnaseq v${workflow.manifest.version}  |  input: ${params.input}  |  outdir: ${params.outdir}")
        RNASEQ()
    }
}
