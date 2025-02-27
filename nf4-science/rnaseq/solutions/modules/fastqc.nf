#!/usr/bin/env nextflow

process FASTQC {
    publishDir "results/fastqc", mode: 'copy'

    input:
    path reads

    output:
    path "${reads.simpleName}_fastqc.zip", emit: zip
    path "${reads.simpleName}_fastqc.html", emit: html

    script:
    """
    fastqc $reads
    """
}
