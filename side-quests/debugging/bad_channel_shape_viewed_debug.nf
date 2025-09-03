#!/usr/bin/env nextflow

process PROCESS_FILES {
    debug true

    input:
    val sample_name

    output:
    path "${sample_name}_output.txt"

    script:
    """
    echo "Sample name inside process is ${sample_name}"
    echo "Processing ${sample_name}" > ${sample_name}_output.txt
    """
}

workflow {

    input_ch = Channel.of(
            ['sample1', 'file1.txt'],
            ['sample2', 'file2.txt'],
            ['sample3', 'file3.txt'],
        )
        .view { "Channel content: ${it}" }
        .map { tuple -> tuple[0] }
        .view { "After mapping: ${it}" }

    PROCESS_FILES(input_ch)
}
