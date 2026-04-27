#!/usr/bin/env nextflow

process PROCESS_FILES {

    container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'

    time '1 ms'

    input:
    val sample_name

    output:
    path "${sample_name}_output.txt"

    script:
    """
    sleep 1  # Takes 1 second, but time limit is 1ms
    cowpy ${sample_name} > ${sample_name}_output.txt
    """
}

workflow {

    // Create input channel
    input_ch = channel.of('sample1', 'sample2', 'sample3')

    // Call the process with the input channel
    PROCESS_FILES(input_ch)
}
