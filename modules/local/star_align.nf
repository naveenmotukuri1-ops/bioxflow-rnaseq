process STAR_ALIGN {
    tag "${meta.id}"
    label 'process_high'

    conda 'bioconda::star=2.7.11b'
    container 'quay.io/biocontainers/star:2.7.11b--h5814d7d_9'

    input:
    tuple val(meta), path(reads)
    path index

    output:
    tuple val(meta), path('*.Aligned.sortedByCoord.out.bam'), emit: bam
    tuple val(meta), path('*.ReadsPerGene.out.tab'),          emit: counts
    tuple val(meta), path('*.Log.final.out'),                 emit: log_final
    tuple val(meta), path('*.SJ.out.tab'),                    emit: junctions
    tuple val("${task.process}"), val('star'), eval('STAR --version'), topic: versions

    script:
    def args    = task.ext.args ?: ''
    def sort_gb = task.memory ? (task.memory.toBytes() * 0.5) as long : 2000000000
    """
    STAR --runMode alignReads \\
        --runThreadN ${task.cpus} \\
        --genomeDir ${index} \\
        --readFilesIn ${[reads].flatten().join(' ')} \\
        --readFilesCommand gunzip -c \\
        --outFileNamePrefix ${meta.id}. \\
        --outSAMtype BAM SortedByCoordinate \\
        --limitBAMsortRAM ${sort_gb} \\
        --quantMode GeneCounts \\
        --twopassMode Basic \\
        ${args}
    """
}
