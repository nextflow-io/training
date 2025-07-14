#!/usr/bin/env nextflow

/*
 * Buggy workflow for debugging exercises
 * This workflow contains several intentional bugs for learning purposes
 */

// Parameters with missing validation
params.input = 'data/sample_data.csv'
params.output = 'results'

// Channel with incorrect usage
input_ch = Channel
    .fromPath(params.input)
    .splitCsv(header: true)
    .map { row -> row.sample_id }  // BUG: Missing file_path from tuple

/*
 * Process with input/output mismatch
 */
process processFiles {
    publishDir "${params.output}/processed", mode: 'copy'

    input:
        tuple val(sample_id), path(input_file)  // BUG: Expects tuple but receives single value

    output:
        path "${sample_id}_result.txt"

    script:
    """
    echo "Processing: ${sample_id}"
    # BUG: Using undefined variable
    cat ${input_file} > ${sample_id}_result.txt
    """
}

/*
 * Process with resource issues
 */
process heavyProcess {
    // BUG: No resource requirements specified for resource-intensive task
    publishDir "${params.output}/heavy", mode: 'copy'

    input:
        val sample_id

    output:
        path "${sample_id}_heavy.txt"

    script:
    """
    # Simulate heavy computation
    for i in {1..1000000}; do
        echo "Heavy computation $i for ${sample_id}"
    done > ${sample_id}_heavy.txt
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
    # BUG: Incorrect file handling - input_file might not exist
    if [ -f "${input_file}" ]; then
        cp ${input_file} processed_${input_file}
    fi
    # BUG: No error handling if file doesn't exist
    """
}

/*
 * Main workflow with channel issues
 */
workflow {
    // BUG: Using channel incorrectly
    processed_ch = processFiles(input_ch)

    // BUG: Trying to use channel that doesn't exist
    heavy_ch = heavyProcess(sample_ids)

    // BUG: Incorrect channel usage
    file_ch = Channel.fromPath("*.txt")
    handleFiles(file_ch)
}
