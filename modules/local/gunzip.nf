process GUNZIP {
    tag "${archive.name}"
    label 'process_single'

    conda 'bioconda::samtools=1.21'
    container 'quay.io/biocontainers/samtools:1.21--h96c455f_1'

    input:
    path archive

    output:
    path "${archive.baseName}", emit: gunzip

    script:
    """
    gunzip -c ${archive} > ${archive.baseName}
    """
}
