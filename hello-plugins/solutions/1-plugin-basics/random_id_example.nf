#!/usr/bin/env nextflow

include { randomString } from 'plugin/nf-hello'

workflow {
    // Generate random IDs for each sample
    Channel.of('sample_A', 'sample_B', 'sample_C')
        .map { sample -> "${sample}_${randomString(8)}" }
        .view()
}
