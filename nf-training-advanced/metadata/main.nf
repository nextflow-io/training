#!/usr/bin/env nextflow

workflow {
    Channel.fromFilePairs("data/reads/*/*_R{1,2}.fastq.gz")
    | view
}
