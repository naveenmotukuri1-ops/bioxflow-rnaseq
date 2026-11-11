process SALMON_QUANT {
    tag "${meta.id}"
    label 'process_medium'

    conda 'bioconda::salmon=1.10.3'
    container 'quay.io/biocontainers/salmon:1.10.3--h45fbf2d_5'

    input:
    tuple val(meta), path(reads)
    path index

    output:
    tuple val(meta), path("${meta.id}"), emit: results
    tuple val("${task.process}"), val('salmon'), eval('salmon --version | sed "s/salmon //"'), topic: versions

    script:
    def args = task.ext.args ?: ''
    def r    = [reads].flatten()
    // library type is detected automatically (-l A) and later cross-checked against STAR
    def io   = meta.single_end ? "-r ${r[0]}" : "-1 ${r[0]} -2 ${r[1]} --gcBias"
    """
    salmon quant \\
        --index ${index} \\
        --libType A \\
        ${io} \\
        --threads ${task.cpus} \\
        ${args} \\
        --output ${meta.id}
    """
}
