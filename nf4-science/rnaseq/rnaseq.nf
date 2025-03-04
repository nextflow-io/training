#!/usr/bin/env nextflow

// Module INCLUDE statements
include { FASTQC } from './modules/fastqc.nf'
include {TRIM_GALORE} from './modules/trim_galore.nf'
include {HISAT2_ALIGN} from './modules/hisat2_align.nf'
/*
 * Pipeline parameters
 */
params.hisat2_index_zip = "data/genome_index.tar.gz"

// Primary input
params.reads = "data/reads/ENCSR000COQ1_1.fastq.gz"
workflow {

    // Create input channel
    read_ch = channel.fromPath(params.reads)

    // Call processes
    // initial fastqc
    FASTQC(read_ch)

    // adapter trimming and post trimming  QC
    TRIM_GALORE(read_ch)

    // align to a reference genome
    HISAT2_ALIGN(TRIM_GALORE.out.trimmed_reads, file (params.hisat2_index_zip))
}
