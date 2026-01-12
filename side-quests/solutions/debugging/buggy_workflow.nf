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
    // Fixed: removed trailing comma
    path "${sample_id}_result.txt"

    script:
    """
    # Fixed: Used correct variable name
    echo "Processing: ${sample_id}"
    cat ${input_file} > ${sample_id}_result.txt
    """
}
// Fixed to add missing closing brace

/*
 * Process with resource issues
 */
process heavyProcess {
    publishDir "${params.output}/heavy", mode: 'copy'

    // Fixed: Increase execution time to avoid timeout
    time '100 s'

    input:
    val sample_id

    output:
    path "${sample_id}_heavy.txt"

    script:
    """
    # Simulate heavy computation
    for i in {1..10000}; do
        echo "Heavy computation \${i} for ${sample_id}" # Fixed: escaped variable
    done > ${sample_id}_heavy.txt # Fixed: output file name corrected
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
    """
}

/*
 * Main workflow with channel issues
 */
workflow {

    // Channel with incorrect usage
    input_ch = channel.fromPath(params.input)
        .splitCsv(header: true)
        .map { row -> [row.sample_id, file(row.fastq_path)] }
    // Fixed: made the map return a tuple with sample_id and fastq_path as a file

    processed_ch = processFiles(input_ch)

    heavy_ch = heavyProcess(input_ch.map { it[0] })
    // Fixed: using the first element of the tuple as input

    // BUG: Incorrect channel usage
    // Fixed: take output from second process
    handleFiles(heavyProcess.out)
}
