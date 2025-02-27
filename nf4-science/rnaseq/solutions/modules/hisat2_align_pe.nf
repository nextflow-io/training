#!/usr/bin/env nextflow

process HISAT2_ALIGN {
    publishDir "results/align", mode: 'copy'

    input:
    tuple path(read1), path(read2)
    path index
    path splice_sites

    output:
    path "${read1.simpleName}.bam", emit: bam
    path "${read1.simpleName}.hisat2.log", emit: log

    script:
    """
    hisat2 -x ${index.simpleName} -1 ${read1} -2 ${read2} --known-splicesite-infile $splice_sites --new-summary --summary-file ${read1.simpleName}.hisat2.log | \
    samtools view -bS - > ${read1.simpleName}.bam
    """
}
