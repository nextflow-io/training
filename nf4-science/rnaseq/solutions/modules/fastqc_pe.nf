#!/usr/bin/env nextflow

process FASTQC {
    publishDir "results/fastqc", mode: 'copy'

    input:
    path reads

    output:
    path "*_fastqc.zip", emit: zip
    path "*_fastqc.html", emit: html

    script:
    """
    fastqc ${reads}
    """
}
