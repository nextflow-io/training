#!/usr/bin/env nextflow

/*
 * Example: Local function vs plugin function
 * This script uses a locally defined randomString function.
 * Compare with importing randomString from the nf-hello plugin.
 */

// Import function from plugin - no local definition needed
include { randomString } from 'plugin/nf-hello'

workflow {
    // Generate random IDs for each sample
    Channel.of('sample_A', 'sample_B', 'sample_C')
        | map { sample -> "${sample}_${randomString(8)}" }
        | view
}
