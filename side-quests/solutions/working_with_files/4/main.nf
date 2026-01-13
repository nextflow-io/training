#!/usr/bin/env nextflow

include { COUNT_LINES } from './modules/count_lines.nf'

workflow {

    // Load files with channel.fromPath
    ch_files = channel.fromPath('data/patientA_rep1_normal_R*_001.fastq.gz')
    ch_files.map { myFile ->
        def (patient, replicate, type, readNum) = myFile.simpleName.tokenize('_')
        [
            [
            id: patient,
            replicate: replicate.replace('rep', ''),
            type: type,
            readNum: readNum.replace('R', ''),
            ],
            myFile
        ]
    }
    .view()

    // Count the lines in the file
    COUNT_LINES(ch_files)
}
