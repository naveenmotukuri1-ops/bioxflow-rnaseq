process CAT_FASTQ {
    tag "${meta.id}"
    label 'process_single'

    conda 'bioconda::samtools=1.21'
    container 'quay.io/biocontainers/samtools:1.21--h96c455f_1'

    input:
    tuple val(meta), path(reads_1, stageAs: 'r1/?/*'), path(reads_2, stageAs: 'r2/?/*')

    output:
    tuple val(meta), path("${meta.id}_*.merged.fastq.gz"), emit: reads

    script:
    def second = meta.single_end ? '' : "cat ${reads_2} > ${meta.id}_2.merged.fastq.gz"
    """
    cat ${reads_1} > ${meta.id}_1.merged.fastq.gz
    ${second}
    """
}
