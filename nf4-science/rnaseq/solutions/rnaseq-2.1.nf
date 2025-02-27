#!/usr/bin/env nextflow

// Module INCLUDE statements
include { FASTQC } from './modules/fastqc'

/*
 * Pipeline parameters
 */
params.hisat2_index = "path/to/hisat2/index"
params.splice_sites = "path/to/splice_sites.txt"

// Primary input
params.reads = "path/to/reads/*.fastq.gz"

workflow {

    // Create input channel from a file path
    read_ch = channel.fromPath(params.reads)

    // Initial quality control
    FASTQC(read_ch)

}
