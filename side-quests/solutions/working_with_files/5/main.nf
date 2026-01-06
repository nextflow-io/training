#!/usr/bin/env nextflow

workflow {

    // Load files with channel.fromFilePairs
    ch_files = channel.fromFilePairs('data/patientA_rep1_normal_R{1,2}_001.fastq.gz')
    ch_files.map { id, files ->
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
    .view()

}
