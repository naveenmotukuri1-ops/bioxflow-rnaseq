process STAR_GENOMEGENERATE {
    tag "${fasta.name}"
    label 'process_high'

    conda 'bioconda::star=2.7.11b'
    container 'quay.io/biocontainers/star:2.7.11b--h5814d7d_9'

    input:
    path fasta
    path gtf

    output:
    path 'star_index', emit: index
    tuple val("${task.process}"), val('star'), eval('STAR --version'), topic: versions

    script:
    def args = task.ext.args ?: ''
    """
    # scale the suffix-array index to the genome size (STAR manual, small genomes)
    GENOME_LEN=\$(grep -v '^>' ${fasta} | tr -d '\\n' | wc -c)
    SA_NBASES=\$(awk -v n="\$GENOME_LEN" 'BEGIN { v = int(log(n) / log(2) / 2 - 1); if (v > 14) v = 14; print v }')

    mkdir star_index
    STAR --runMode genomeGenerate \\
        --runThreadN ${task.cpus} \\
        --genomeDir star_index \\
        --genomeFastaFiles ${fasta} \\
        --sjdbGTFfile ${gtf} \\
        --sjdbOverhang 100 \\
        --genomeSAindexNbases \$SA_NBASES \\
        ${args}
    """
}
