#!/usr/bin/env nextflow

include { ANALYZE_READS } from './modules/analyze_reads.nf'

workflow {
    main:
    // Load files with channel.fromFilePairs
    ch_files = channel.fromFilePairs('data/*_R{1,2}_001.fastq.gz')
    ch_samples = ch_files.map { id, files ->
            def (sample, replicate, type, readNum) = id.tokenize('_')
            tuple(
                [
                    id: sample,
                    replicate: replicate.replace('rep', ''),
                    type: type,
                ],
                files,
            )
        }

    // Run the analysis
    ANALYZE_READS(ch_samples)

    publish:
    analysis_results = ANALYZE_READS.out
}

output {
    analysis_results {
        path { meta, file -> "${meta.type}/${meta.id}/${meta.replicate}" }
    }
}
