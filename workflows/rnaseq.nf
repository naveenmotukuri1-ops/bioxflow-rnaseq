/*
 * bioxflow/rnaseq main workflow
 */

include { PREPARE_GENOME    } from '../subworkflows/local/prepare_genome'
include { CAT_FASTQ         } from '../modules/local/cat_fastq'
include { FASTQC            } from '../modules/local/fastqc'
include { FASTP             } from '../modules/local/fastp'
include { SALMON_QUANT      } from '../modules/local/salmon_quant'
include { STAR_ALIGN        } from '../modules/local/star_align'
include { SAMTOOLS_STATS    } from '../modules/local/samtools_stats'
include { QUANT_CONCORDANCE } from '../modules/local/quant_concordance'
include { DESEQ2            } from '../modules/local/deseq2'
include { MULTIQC           } from '../modules/local/multiqc'

// Relative FASTQ paths are looked up in the launch directory first, then next to the samplesheet
def resolveFastq(String path, sheet) {
    if (path.startsWith('/') || path.contains('://')) {
        return file(path, checkIfExists: true)
    }
    def local = file(path)
    if (local.exists()) {
        return local
    }
    return file("${sheet.parent}/${path}", checkIfExists: true)
}

// Read and validate the samplesheet. Rows sharing a sample name are treated as
// sequencing runs of the same library and are concatenated.
def parseSamplesheet(sheet) {
    def rows = sheet.splitCsv(header: true, strip: true)
    if (!rows) {
        error("Samplesheet ${sheet} is empty")
    }
    def missing = ['sample', 'fastq_1', 'condition'].findAll { col -> !rows[0].containsKey(col) }
    if (missing) {
        error("Samplesheet is missing required column(s): ${missing.join(', ')}. Expected: sample,fastq_1,fastq_2,condition[,batch]")
    }
    def samples = [:]
    rows.eachWithIndex { row, i ->
        def line = i + 2
        def id = row.sample
        if (!id || !(id ==~ /^[A-Za-z0-9][A-Za-z0-9_.-]*$/)) {
            error("Samplesheet line ${line}: sample name '${id}' must contain only letters, numbers, '_', '.' or '-'")
        }
        if (!row.condition || !(row.condition ==~ /^[A-Za-z][A-Za-z0-9_.]*$/)) {
            error("Samplesheet line ${line}: condition '${row.condition}' must start with a letter and contain only letters, numbers, '_' or '.'")
        }
        if (!row.fastq_1) {
            error("Samplesheet line ${line}: fastq_1 is empty")
        }
        def r1 = resolveFastq(row.fastq_1, sheet)
        def r2 = row.fastq_2 ? resolveFastq(row.fastq_2, sheet) : null
        def batch = row.batch ?: 'NA'
        def meta = [id: id, single_end: r2 == null, condition: row.condition, batch: batch]
        if (!samples.containsKey(id)) {
            samples[id] = [meta: meta, r1: [], r2: []]
        }
        else if (samples[id].meta != meta) {
            error("Samplesheet line ${line}: sample '${id}' appears with a different condition, batch or read layout than before")
        }
        samples[id].r1 << r1
        if (r2) {
            samples[id].r2 << r2
        }
    }
    return samples.values().collect { s -> [s.meta, s.r1, s.r2] }
}

workflow RNASEQ {
    main:
    if (!params.input) { error("Please provide a samplesheet with --input") }
    if (!params.fasta) { error("Please provide a genome FASTA with --fasta") }
    if (!params.gtf)   { error("Please provide a gene annotation with --gtf") }

    def sheet   = file(params.input, checkIfExists: true)
    def samples = parseSamplesheet(sheet)
    def conditions = samples.collect { s -> s[0].condition }.unique()
    def reference  = params.reference_condition ?: conditions[0]
    if (!(reference in conditions)) {
        error("--reference_condition '${reference}' is not a condition in the samplesheet (${conditions.join(', ')})")
    }
    log.info("Samples: ${samples.size()}  |  conditions: ${conditions.join(', ')}  |  reference: ${reference}")

    ch_samples_tsv = channel
        .fromList(samples)
        .map { meta, _r1, _r2 -> "${meta.id}\t${meta.condition}\t${meta.batch}\n" }
        .collectFile(name: 'samples.tsv', seed: "sample\tcondition\tbatch\n", sort: false)
        .first()

    //
    // Reference
    //
    PREPARE_GENOME()

    //
    // Reads: merge multi-run samples, QC, trim
    //
    ch_input = channel.fromList(samples).branch { _meta, r1, _r2 ->
        multiple: r1.size() > 1
        single: true
    }
    CAT_FASTQ(ch_input.multiple)
    ch_raw = ch_input.single
        .map { meta, r1, r2 -> [meta, r1 + r2] }
        .mix(CAT_FASTQ.out.reads)

    ch_multiqc = channel.empty()

    if (!params.skip_fastqc) {
        FASTQC(ch_raw)
        ch_multiqc = ch_multiqc.mix(FASTQC.out.zip.map { _meta, zip -> zip })
    }

    ch_reads = ch_raw
    if (!params.skip_trimming) {
        FASTP(ch_raw)
        ch_reads   = FASTP.out.reads
        ch_multiqc = ch_multiqc.mix(FASTP.out.json.map { _meta, json -> json })
    }

    //
    // Route 1: Salmon selective alignment against a decoy-aware transcriptome
    //
    SALMON_QUANT(ch_reads, PREPARE_GENOME.out.salmon_index)
    ch_salmon_dirs = SALMON_QUANT.out.results.map { _meta, dir -> dir }
    ch_multiqc     = ch_multiqc.mix(ch_salmon_dirs)

    //
    // Route 2: STAR genome alignment with gene counting
    //
    ch_star_counts = channel.value([])
    if (!params.skip_alignment) {
        STAR_ALIGN(ch_reads, PREPARE_GENOME.out.star_index)
        SAMTOOLS_STATS(STAR_ALIGN.out.bam)
        ch_star_counts = STAR_ALIGN.out.counts.map { _meta, tab -> tab }.collect()
        ch_multiqc = ch_multiqc
            .mix(STAR_ALIGN.out.log_final.map { _meta, log -> log })
            .mix(STAR_ALIGN.out.counts.map { _meta, tab -> tab })
            .mix(SAMTOOLS_STATS.out.flagstat.map { _meta, f -> f })
    }

    //
    // Gene-level summarisation and cross-method concordance
    //
    QUANT_CONCORDANCE(
        ch_samples_tsv,
        PREPARE_GENOME.out.tx2gene,
        ch_salmon_dirs.collect(),
        ch_star_counts,
        params.concordance_threshold
    )
    ch_multiqc = ch_multiqc.mix(QUANT_CONCORDANCE.out.multiqc.flatten())

    //
    // Differential expression
    //
    if (!params.skip_de) {
        DESEQ2(
            QUANT_CONCORDANCE.out.txi,
            ch_samples_tsv,
            PREPARE_GENOME.out.tx2gene,
            params.contrasts ? file(params.contrasts, checkIfExists: true) : [],
            reference,
            params.alpha,
            params.lfc_threshold,
            params.min_count
        )
        ch_multiqc = ch_multiqc.mix(DESEQ2.out.multiqc.flatten())
    }

    //
    // Software versions and MultiQC report
    //
    ch_versions = channel
        .topic('versions')
        .map { _process, tool, version -> "${tool}: \"${version.toString().trim()}\"" }
        .unique()
        .collectFile(name: 'bioxflow_mqc_versions.yml', newLine: true, sort: true, storeDir: "${params.outdir}/pipeline_info")

    MULTIQC(
        file("${projectDir}/assets/multiqc_config.yml", checkIfExists: true),
        ch_multiqc.mix(ch_versions).collect()
    )
}
