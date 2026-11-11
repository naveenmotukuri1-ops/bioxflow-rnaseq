process GFFREAD {
    tag "${gtf.name}"
    label 'process_single'

    conda 'bioconda::gffread=0.12.7'
    container 'quay.io/biocontainers/gffread:0.12.7--h077b44d_6'

    input:
    path fasta
    path gtf

    output:
    path 'transcripts.fa', emit: transcripts
    tuple val("${task.process}"), val('gffread'), eval('gffread --version 2>&1'), topic: versions

    script:
    """
    gffread ${gtf} -g ${fasta} -w transcripts.fa
    """
}
