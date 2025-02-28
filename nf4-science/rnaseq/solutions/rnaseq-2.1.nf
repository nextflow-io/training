#!/usr/bin/env nextflow

// Module INCLUDE statements
include { FASTQC } from './modules/fastqc.nf'

/*
 * Pipeline parameters
 */

// Primary input
params.reads = "data/reads/ENCSR000COQ1_1.fastq.gz"

workflow {

    // Create input channel from a file path
    read_ch = channel.fromPath(params.reads)

    // Initial quality control
    FASTQC(read_ch)

}
