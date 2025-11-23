#!/usr/bin/env nextflow

process PROCESS_FILES {
    input:
        val sample_name  // Process expects only 1 input channel

    output:
        path "${sample_name}_output.txt"

    script:
    """
    echo "Processing ${sample_name}" > ${sample_name}_output.txt
    """
}

workflow {

    // Create two separate channels
    samples_ch = channel.of('sample1', 'sample2', 'sample3')
    files_ch = channel.of('file1.txt', 'file2.txt', 'file3.txt')

    // ERROR: Passing 2 channels but process expects only 1
    PROCESS_FILES(samples_ch, files_ch)
}
