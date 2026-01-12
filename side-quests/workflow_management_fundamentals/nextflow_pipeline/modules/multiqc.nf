/*
 * MULTIQC - Aggregate QC reports from multiple tools
 *
 * This process collects outputs from FastQC, fastp, and Salmon
 * into a single interactive HTML report.
 */
process MULTIQC {
    container 'biocontainers/multiqc:1.25.1--pyhdfd78af_0'
    publishDir "${params.outdir}", mode: 'copy'

    cpus 1
    memory '2.GB'

    input:
    path('*')

    output:
    path "multiqc_report.html", emit: report
    path "multiqc_data", emit: data

    script:
    """
    multiqc . --filename multiqc_report
    """
}
