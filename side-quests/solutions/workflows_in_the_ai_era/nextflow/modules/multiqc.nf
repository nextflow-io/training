/*
 * MULTIQC - Aggregate QC metrics from FastQC, fastp, and Salmon
 */
process MULTIQC {
    container 'quay.io/biocontainers/multiqc:1.19--pyhdfd78af_0'

    input:
    path('*')

    output:
    path "multiqc_report.html", emit: report
    path "multiqc_data", emit: data

    script:
    """
    multiqc .
    """
}
