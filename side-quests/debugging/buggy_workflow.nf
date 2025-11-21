#!/usr/bin/env nextflow

/*
 * Buggy workflow for debugging exercises
 * This workflow contains several intentional bugs for learning purposes
 */

// Parameters with missing validation
params.input = 'data/sample_data.csv'
params.output = 'results'

/*
 * Process with input/output mismatch
 */
process processFiles {
    publishDir "${params.output}/processed", mode: 'copy'

    input:
        tuple val(sample_id), path(input_file)

    output:
        path "${sample_id}_result.txt",

    script:
    """
    echo "Processing: ${sample}"
    cat ${input_file} > ${sample}_result.txt
    """

/*
 * Process with resource issues
 */
process heavyProcess {
    publishDir "${params.output}/heavy", mode: 'copy'

    time '1 ms'

    input:
        val sample_id

    output:
        path "${sample_id}_heavy.txt"

    script:
    """
    # Simulate heavy computation
    for i in {1..1000000}; do
        echo "Heavy computation $i for ${sample_id}"
    done > ${sample_id}.txt
    """
}

/*
 * Process with file handling issues
 */
process handleFiles {
    publishDir "${params.output}/files", mode: 'copy'

    input:
        path input_file

    output:
        path "processed_${input_file}"

    script:
    """
    if [ -f "${input_file}" ]; then
        cp ${input_file} processed_${input_file}
    fi
    """
}

/*
 * Main workflow with channel issues
 */
workflow {

    // Channel with incorrect usage
    input_ch = channel
        .fromPath(params.input)
        .splitCsv(header: true)
        .map { row -> row.sample_id }

    processed_ch = processFiles(input_ch)

    heavy_ch = heavyProcess(sample_ids)

    file_ch = channel.fromPath("*.txt")
    handleFiles(file_ch)
}
