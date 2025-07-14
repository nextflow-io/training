#!/usr/bin/env nextflow

/*
 * Formatted workflow - properly styled version of messy_workflow.nf
 * This workflow demonstrates proper formatting and coding standards
 */

// Parameters with proper spacing and quotes
params.input = 'data/sample_data.csv'
params.output = 'results'

// Input channel with clear formatting and descriptive name
samples_ch = Channel
    .fromPath(params.input)
    .splitCsv(header: true)
    .map { row -> tuple(row.sample_id, row.file_path) }

/*
 * Process samples with proper formatting
 * Input: Sample ID and file path
 * Output: Processed sample file
 */
process processSamples {
    publishDir "${params.output}/processed", mode: 'copy'

    input:
        tuple val(sample_id), val(file_path)

    output:
        tuple val(sample_id), path("${sample_id}_processed.txt")

    script:
    """
    echo "Processing sample: ${sample_id}"
    echo "Input file: ${file_path}" > ${sample_id}_processed.txt
    echo "Processing completed at: \$(date)" >> ${sample_id}_processed.txt
    """
}

/*
 * Final processing step with consistent formatting
 * Input: Sample ID and processed file
 * Output: Final processed file
 */
process finalProcessing {
    publishDir "${params.output}/final", mode: 'copy'

    input:
        tuple val(sample_id), path(processed_file)

    output:
        path "${sample_id}_final.txt"

    script:
    """
    echo "Final processing for sample: ${sample_id}"
    cat ${processed_file} > ${sample_id}_final.txt
    echo "Final processing completed at: \$(date)" >> ${sample_id}_final.txt
    """
}

/*
 * Main workflow with clear structure
 */
workflow {
    // Process samples
    processed_samples_ch = processSamples(samples_ch)

    // Final processing
    final_output_ch = finalProcessing(processed_samples_ch)

    // Display results
    final_output_ch.view { file ->
        "Final output created: ${file}"
    }
}
