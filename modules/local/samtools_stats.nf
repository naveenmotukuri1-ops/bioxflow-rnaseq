process SAMTOOLS_STATS {
    tag "${meta.id}"
    label 'process_low'

    conda 'bioconda::samtools=1.21'
    container 'quay.io/biocontainers/samtools:1.21--h96c455f_1'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path('*.bai'),      emit: bai
    tuple val(meta), path('*.flagstat'), emit: flagstat
    tuple val(meta), path('*.idxstats'), emit: idxstats
    tuple val("${task.process}"), val('samtools'), eval('samtools --version | head -n 1 | sed "s/samtools //"'), topic: versions

    script:
    """
    samtools index -@ ${task.cpus} ${bam}
    samtools flagstat -@ ${task.cpus} ${bam} > ${meta.id}.flagstat
    samtools idxstats ${bam} > ${meta.id}.idxstats
    """
}
