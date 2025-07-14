#!/usr/bin/env nextflow

/*
 * Debugged workflow - corrected version of buggy_workflow.nf
 * This workflow demonstrates proper debugging and error handling
 */

// Parameters with validation
params.input = 'data/sample_data.csv'
params.output = 'results'

// Validate required parameters
if (!params.input) {
    error "Please provide an input file with --input"
}

// Channel with proper structure
input_ch = Channel
    .fromPath(params.input)
    .splitCsv(header: true)
    .map { row -> tuple(row.sample_id, row.file_path) }  // FIXED: Include both values in tuple

/*
 * Process with correct input/output matching
 */
process processFiles {
    publishDir "${params.output}/processed", mode: 'copy'

    input:
        tuple val(sample_id), val(file_path)  // FIXED: Match the channel structure

    output:
        tuple val(sample_id), path("${sample_id}_result.txt")

    script:
    """
    echo "Processing: ${sample_id}"
    echo "File path: ${file_path}" > ${sample_id}_result.txt
    echo "Processed at: \$(date)" >> ${sample_id}_result.txt
    """
}

/*
 * Process with proper resource requirements
 */
process heavyProcess {
    // FIXED: Add appropriate resource requirements
    cpus 4
    memory '8.GB'
    time '2.h'

    publishDir "${params.output}/heavy", mode: 'copy'

    input:
        tuple val(sample_id), path(processed_file)

    output:
        path "${sample_id}_heavy.txt"

    script:
    """
    # Simulate heavy computation with proper resource allocation
    echo "Starting heavy computation for ${sample_id}"
    for i in {1..1000}; do
        echo "Heavy computation \$i for ${sample_id}"
    done > ${sample_id}_heavy.txt
    echo "Completed heavy computation for ${sample_id}"
    """
}

/*
 * Process with proper file handling and error checking
 */
process handleFiles {
    publishDir "${params.output}/files", mode: 'copy'

    input:
        tuple val(sample_id), path(input_file)

    output:
        path "processed_${sample_id}.txt"

    script:
    """
    # FIXED: Proper error handling and file validation
    if [ -f "${input_file}" ]; then
        cp ${input_file} processed_${sample_id}.txt
        echo "Successfully processed file for ${sample_id}"
    else
        echo "Error: Input file ${input_file} not found for ${sample_id}"
        exit 1
    fi
    """
}

/*
 * Main workflow with correct channel handling
 */
workflow {
    // FIXED: Use channels correctly
    processed_ch = processFiles(input_ch)

    // FIXED: Use existing channel
    heavy_ch = heavyProcess(processed_ch)

    // FIXED: Use proper channel and process combination
    file_ch = handleFiles(processed_ch)

    // Debug output
    processed_ch.view { sample_id, file ->
        "Processed: ${sample_id} -> ${file}"
    }

    heavy_ch.view { file ->
        "Heavy processing complete: ${file}"
    }

    file_ch.view { file ->
        "File handling complete: ${file}"
    }
}
