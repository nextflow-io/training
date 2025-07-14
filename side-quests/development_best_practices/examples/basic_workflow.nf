#!/usr/bin/env nextflow

/*
 * Basic Nextflow workflow for demonstrating development best practices
 * This workflow reads a CSV file and processes each row
 */

// Default parameters
params.input = 'data/sample_data.csv'
params.output = 'results'

// Input channel from CSV file
input_ch = Channel
    .fromPath(params.input)
    .splitCsv(header: true)
    .map { row -> tuple(row.sample_id, row.file_path) }

/*
 * Process each sample
 */
process processSample {
    publishDir "${params.output}/processed", mode: 'copy'

    input:
        tuple val(sample_id), val(file_path)

    output:
        tuple val(sample_id), path("${sample_id}_processed.txt")

    script:
    """
    echo "Processing sample: ${sample_id}"
    echo "Input file: ${file_path}" > ${sample_id}_processed.txt
    echo "Processed at: \$(date)" >> ${sample_id}_processed.txt
    """
}

/*
 * Main workflow
 */
workflow {
    // Process samples
    processed_ch = processSample(input_ch)

    // Print completion message
    processed_ch.view { sample_id, file ->
        "Completed processing: ${sample_id}"
    }
}
