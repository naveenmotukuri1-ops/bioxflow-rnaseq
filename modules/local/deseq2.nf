process DESEQ2 {
    label 'process_low'

    conda 'conda-forge::r-base=4.4 bioconda::bioconductor-deseq2=1.46.0 bioconda::bioconductor-tximport=1.34.0 bioconda::bioconductor-apeglm=1.28.0 conda-forge::r-ggplot2 conda-forge::r-pheatmap conda-forge::r-ggrepel'
    container 'community.wave.seqera.io/library/r-base_bioconductor-deseq2_bioconductor-tximport_bioconductor-apeglm_pruned:2015b8d56d444ba2'

    input:
    path txi
    path samples
    path tx2gene
    path contrasts
    val  reference
    val  alpha
    val  lfc_threshold
    val  min_count

    output:
    path 'deseq2/**',  emit: results
    path '*_mqc.*',    emit: multiqc
    tuple val("${task.process}"), val('deseq2'), eval("Rscript -e \"cat(as.character(packageVersion('DESeq2')))\""), topic: versions
    tuple val("${task.process}"), val('apeglm'), eval("Rscript -e \"cat(as.character(packageVersion('apeglm')))\""), topic: versions

    script:
    def con = contrasts ? "--contrasts ${contrasts}" : ''
    """
    bioxflow_deseq2.R \\
        --txi ${txi} \\
        --samples ${samples} \\
        --tx2gene ${tx2gene} \\
        --reference ${reference} \\
        --alpha ${alpha} \\
        --lfc_threshold ${lfc_threshold} \\
        --min_count ${min_count} \\
        ${con} \\
        --outdir deseq2
    """
}
