process QUANT_CONCORDANCE {
    label 'process_low'

    conda 'conda-forge::r-base=4.4 bioconda::bioconductor-deseq2=1.46.0 bioconda::bioconductor-tximport=1.34.0 bioconda::bioconductor-apeglm=1.28.0 conda-forge::r-ggplot2 conda-forge::r-pheatmap conda-forge::r-ggrepel'
    container 'community.wave.seqera.io/library/r-base_bioconductor-deseq2_bioconductor-tximport_bioconductor-apeglm_pruned:2015b8d56d444ba2'

    input:
    path samples
    path tx2gene
    path salmon_dirs, stageAs: 'salmon/*'
    path star_counts, stageAs: 'star/*'
    val  threshold

    output:
    path 'txi.rds',                    emit: txi
    path 'salmon.*.tsv',               emit: salmon_tables
    path 'star.gene_counts.tsv',       emit: star_table, optional: true
    path 'strandedness.tsv',           emit: strandedness
    path 'concordance_summary.tsv',    emit: summary, optional: true
    path 'concordance_scatter.png',    emit: scatter, optional: true
    path '*_mqc.*',                    emit: multiqc
    tuple val("${task.process}"), val('tximport'), eval("Rscript -e \"cat(as.character(packageVersion('tximport')))\""), topic: versions

    script:
    """
    bioxflow_quant.R \\
        --samples ${samples} \\
        --tx2gene ${tx2gene} \\
        --salmon salmon \\
        --star star \\
        --threshold ${threshold}
    """
}
