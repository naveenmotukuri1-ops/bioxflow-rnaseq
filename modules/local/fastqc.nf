process FASTQC {
    tag "${meta.id}"
    label 'process_low'

    conda 'bioconda::fastqc=0.12.1'
    container 'quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0'

    input:
    tuple val(meta), path(reads, stageAs: 'input/*')

    output:
    tuple val(meta), path('*.html'), emit: html
    tuple val(meta), path('*.zip'),  emit: zip
    tuple val("${task.process}"), val('fastqc'), eval('fastqc --version | sed "s/FastQC v//"'), topic: versions

    script:
    def args  = task.ext.args ?: ''
    def names = meta.single_end ? ["${meta.id}.fastq.gz"] : ["${meta.id}_1.fastq.gz", "${meta.id}_2.fastq.gz"]
    def links = [reads].flatten().withIndex().collect { f, i -> "ln -sf ${f} ${names[i]}" }.join('\n    ')
    """
    ${links}
    fastqc ${args} --threads ${task.cpus} ${names.join(' ')}
    """
}
