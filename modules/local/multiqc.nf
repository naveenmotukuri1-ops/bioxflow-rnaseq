process MULTIQC {
    label 'process_low'

    conda 'bioconda::multiqc=1.25.1'
    container 'quay.io/biocontainers/multiqc:1.25.1--pyhdfd78af_0'

    input:
    path config
    path 'inputs/*'

    output:
    path 'multiqc_report.html', emit: report
    path 'multiqc_report_data', emit: data

    script:
    """
    multiqc --force --config ${config} --filename multiqc_report.html inputs
    """
}
