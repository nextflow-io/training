#!/usr/bin/env nextflow

include { COUNT_LINES } from './modules/count_lines.nf'

workflow {

    // Create a Path object from a string path
    myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

    // Print file attributes
    println "File object class: ${myFile.class}"
    println "File name: ${myFile.name}"
    println "Simple name: ${myFile.simpleName}"
    println "Extension: ${myFile.extension}"
    println "Parent directory: ${myFile.parent}"

    // Count the lines in the file
    COUNT_LINES(myFile)
}
