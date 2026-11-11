process TX2GENE {
    tag "${gtf.name}"
    label 'process_single'

    conda 'conda-forge::r-base=4.4 bioconda::bioconductor-deseq2=1.46.0 bioconda::bioconductor-tximport=1.34.0 bioconda::bioconductor-apeglm=1.28.0 conda-forge::r-ggplot2 conda-forge::r-pheatmap conda-forge::r-ggrepel'
    container 'community.wave.seqera.io/library/r-base_bioconductor-deseq2_bioconductor-tximport_bioconductor-apeglm_pruned:2015b8d56d444ba2'

    input:
    path gtf

    output:
    path 'tx2gene.tsv', emit: tsv

    script:
    """
    gtf2tx2gene.R ${gtf} tx2gene.tsv
    """
}
