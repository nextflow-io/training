#!/usr/bin/env nextflow

workflow {
    main:
    // Create a Path object from a string path
    myFile = 'data/patientA_rep1_normal_R1_001.fastq.gz'

    println "${myFile} is of class ${myFile.class}"

    publish:
    analysis_results = channel.empty()
}

output {
    analysis_results {
    }
}
