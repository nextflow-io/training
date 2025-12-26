#!/usr/bin/env nextflow

include { ANALYZE_READS } from './modules/analyze_reads.nf'

workflow {

    // Create a file object from a string path
    ch_files = channel.fromFilePairs('data/*_R{1,2}_001.fastq.gz')
    ch_samples = ch_files.map { id, files ->
        def (sample, replicate, type) = id.tokenize('_')
        [
            [
                id: sample,
                replicate: replicate.replace('rep', ''),
                type: type
            ],
            files
        ]
    }

    // Run the analysis
    ANALYZE_READS(ch_samples)
}
