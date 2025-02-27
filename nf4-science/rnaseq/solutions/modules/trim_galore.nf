#!/usr/bin/env nextflow

process TRIM_GALORE {
    publishDir "results/trimming", mode: 'copy'

    input:
    path reads

    output:
    path "${reads.simpleName}_trimmed.fq.gz", emit: trimmed_reads
    path "${reads.simpleName}_trimming_report.txt", emit: trimming_reports
    path "${reads.simpleName}_trimmed_fastqc.{zip,html}", emit: fastqc_reports

    script:
    """
    trim_galore --fastqc $reads
    """
}
