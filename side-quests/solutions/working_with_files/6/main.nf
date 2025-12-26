#!/usr/bin/env nextflow

include { ANALYZE_READS } from './modules/analyze_reads.nf'

workflow {

    // Load files with channel.fromFilePairs
    ch_files = channel.fromFilePairs('data/*_R{1,2}_001.fastq.gz')
    ch_files.map { id,  files ->
        def (sample, replicate, type, readNum) = id.tokenize('_')
        [
            [
                id: sample,
                replicate: replicate.replace('rep', ''),
                type: type
            ],
            files
        ]
    }
        .set { ch_samples }

    // Run the analysis
    ANALYZE_READS(ch_samples)

}
