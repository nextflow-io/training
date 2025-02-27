#!/usr/bin/env nextflow

process HISAT2_ALIGN {
    publishDir "results/align", mode: 'copy'

    input:
    path reads
    path index
    path splice_sites

    output:
    path "${reads.simpleName}.bam", emit: bam
    path "${reads.simpleName}.hisat2.log", emit: log

    script:
    """
    hisat2 -x ${index.simpleName} -U $reads --known-splicesite-infile $splice_sites --new-summary --summary-file ${reads.simpleName}.hisat2.log | \
    samtools view -bS - > ${reads.simpleName}.bam
    """
}
