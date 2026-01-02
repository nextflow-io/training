#!/usr/bin/env nextflow

include { COUNT_LINES } from './modules/count_lines.nf'

workflow {

    // Load files with channel.fromPath
    ch_files = channel.fromPath('data/patientA_rep1_normal_R*_001.fastq.gz')
    ch_files.view { myFile ->
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
    }

    // Count the lines in the file
    COUNT_LINES(ch_files)
}
