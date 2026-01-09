#!/usr/bin/env nextflow

// Local function - must be copied to every pipeline that needs it
def randomString(int length) {
    def chars = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789'
    def random = new Random()
    return (1..length).collect { chars[random.nextInt(chars.length())] }.join()
}

workflow {
    // Generate random IDs for each sample
    Channel.of('sample_A', 'sample_B', 'sample_C')
        | map { sample -> "${sample}_${randomString(8)}" }
        | view
}
