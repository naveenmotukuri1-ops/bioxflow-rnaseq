process FASTP {
    tag "${meta.id}"
    label 'process_low'

    conda 'bioconda::fastp=0.24.0'
    container 'quay.io/biocontainers/fastp:0.24.0--heae3180_1'

    input:
    tuple val(meta), path(reads, stageAs: 'input/*')

    output:
    tuple val(meta), path('*.trimmed.fastq.gz'), emit: reads
    tuple val(meta), path('*.fastp.json'),       emit: json
    tuple val(meta), path('*.fastp.html'),       emit: html
    tuple val(meta), path('*.fastp.log'),        emit: log
    tuple val("${task.process}"), val('fastp'), eval('fastp --version 2>&1 | sed "s/fastp //"'), topic: versions

    script:
    def args = task.ext.args ?: ''
    def r    = [reads].flatten()
    def io   = meta.single_end
        ? "--in1 ${r[0]} --out1 ${meta.id}.trimmed.fastq.gz"
        : "--in1 ${r[0]} --in2 ${r[1]} --out1 ${meta.id}_1.trimmed.fastq.gz --out2 ${meta.id}_2.trimmed.fastq.gz --detect_adapter_for_pe"
    """
    fastp ${io} \\
        --thread ${task.cpus} \\
        --json ${meta.id}.fastp.json \\
        --html ${meta.id}.fastp.html \\
        ${args} \\
        2> ${meta.id}.fastp.log
    """
}
