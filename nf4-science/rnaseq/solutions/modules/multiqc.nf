#!/usr/bin/env nextflow

process MULTIQC {
    publishDir "results/multiqc", mode: 'copy'

    input:
    path '*'
    val output_name

    output:
    path "${output_name}.html", emit: report
    path "${output_name}_data", emit: data

    script:
    """
    multiqc . -n ${output_name}.html
    """
}
