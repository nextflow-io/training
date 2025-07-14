#!/usr/bin/env nextflow

/*
 * Advanced workflow demonstrating comprehensive best practices
 * This workflow showcases modular design, error handling, and proper documentation
 */

// Parameter validation and defaults
params.input = null
params.output = 'results'
params.skip_qc = false
params.threads = 2

// Validate required parameters
if (!params.input) {
    error "Please provide an input file with --input"
}

// Check if input file exists
if (!file(params.input).exists()) {
    error "Input file does not exist: ${params.input}"
}

// Display parameters
log.info """\
Advanced Nextflow Workflow
========================
Input file    : ${params.input}
Output dir    : ${params.output}
Skip QC       : ${params.skip_qc}
Threads       : ${params.threads}
"""

// Input channel with validation
input_ch = Channel
    .fromPath(params.input, checkIfExists: true)
    .splitCsv(header: true)
    .map { row ->
        if (!row.sample_id || !row.file_path) {
            error "Invalid CSV format: missing sample_id or file_path"
        }
        tuple(row.sample_id, row.file_path, row.type ?: 'unknown')
    }

/*
 * Quality control process
 * Performs initial quality checks on input data
 */
process qualityControl {
    tag "${sample_id}"
    publishDir "${params.output}/qc", mode: 'copy'

    cpus params.threads
    memory '4.GB'
    time '1.h'

    errorStrategy 'retry'
    maxRetries 2

    input:
        tuple val(sample_id), val(file_path), val(type)

    output:
        tuple val(sample_id), path("${sample_id}_qc_report.txt"), val(type)

    when:
        !params.skip_qc

    script:
    """
    echo "Quality control for sample: ${sample_id}"
    echo "Sample ID: ${sample_id}" > ${sample_id}_qc_report.txt
    echo "File path: ${file_path}" >> ${sample_id}_qc_report.txt
    echo "Type: ${type}" >> ${sample_id}_qc_report.txt
    echo "QC completed at: \$(date)" >> ${sample_id}_qc_report.txt

    # Simulate quality control checks
    if [ "${type}" == "paired" ]; then
        echo "Performing paired-end QC checks..." >> ${sample_id}_qc_report.txt
    else
        echo "Performing single-end QC checks..." >> ${sample_id}_qc_report.txt
    fi

    echo "QC Status: PASS" >> ${sample_id}_qc_report.txt
    """

    stub:
    """
    touch ${sample_id}_qc_report.txt
    """
}

/*
 * Data processing process with resource optimization
 * Processes data based on type with appropriate resources
 */
process processData {
    tag "${sample_id}"
    publishDir "${params.output}/processed", mode: 'copy'

    cpus { type == 'paired' ? params.threads * 2 : params.threads }
    memory { type == 'paired' ? '8.GB' : '4.GB' }
    time '2.h'

    errorStrategy { task.attempt <= 2 ? 'retry' : 'ignore' }
    maxRetries 3

    input:
        tuple val(sample_id), path(qc_report), val(type)

    output:
        tuple val(sample_id), path("${sample_id}_processed.txt"), val(type)

    script:
    """
    echo "Processing data for sample: ${sample_id}"
    echo "Sample ID: ${sample_id}" > ${sample_id}_processed.txt
    echo "Type: ${type}" >> ${sample_id}_processed.txt
    echo "Processing started at: \$(date)" >> ${sample_id}_processed.txt

    # Include QC information
    echo "QC Results:" >> ${sample_id}_processed.txt
    cat ${qc_report} >> ${sample_id}_processed.txt

    # Type-specific processing
    if [ "${type}" == "paired" ]; then
        echo "Performing paired-end processing with ${task.cpus} CPUs..." >> ${sample_id}_processed.txt
        # Simulate paired-end processing
        for i in {1..10}; do
            echo "Paired processing step \$i" >> ${sample_id}_processed.txt
        done
    else
        echo "Performing single-end processing with ${task.cpus} CPUs..." >> ${sample_id}_processed.txt
        # Simulate single-end processing
        for i in {1..5}; do
            echo "Single processing step \$i" >> ${sample_id}_processed.txt
        done
    fi

    echo "Processing completed at: \$(date)" >> ${sample_id}_processed.txt
    """
}

/*
 * Summary generation process
 * Creates a comprehensive summary of all processed samples
 */
process generateSummary {
    publishDir "${params.output}/summary", mode: 'copy'

    cpus 1
    memory '2.GB'
    time '30.m'

    input:
        path processed_files

    output:
        path "processing_summary.txt"

    script:
    """
    echo "Processing Summary Report" > processing_summary.txt
    echo "=========================" >> processing_summary.txt
    echo "Generated at: \$(date)" >> processing_summary.txt
    echo "" >> processing_summary.txt

    echo "Total samples processed: \$(ls ${processed_files} | wc -l)" >> processing_summary.txt
    echo "" >> processing_summary.txt

    echo "Individual sample results:" >> processing_summary.txt
    echo "-------------------------" >> processing_summary.txt

    for file in ${processed_files}; do
        echo "File: \$file" >> processing_summary.txt
        echo "Sample details:" >> processing_summary.txt
        head -3 \$file >> processing_summary.txt
        echo "" >> processing_summary.txt
    done
    """
}

/*
 * Main workflow
 */
workflow {
    // Quality control (conditional)
    if (params.skip_qc) {
        // Skip QC and create dummy channel
        qc_ch = input_ch.map { sample_id, file_path, type ->
            tuple(sample_id, file("${workflow.workDir}/dummy_qc.txt"), type)
        }
    } else {
        qc_ch = qualityControl(input_ch)
    }

    // Process data
    processed_ch = processData(qc_ch)

    // Generate summary
    summary_ch = generateSummary(processed_ch.map { it[1] }.collect())

    // Display progress
    processed_ch.view { sample_id, file, type ->
        "✓ Processed ${sample_id} (${type}): ${file}"
    }

    summary_ch.view { file ->
        "✓ Summary generated: ${file}"
    }
}

/*
 * Workflow completion hook
 */
workflow.onComplete {
    log.info """\
Pipeline completed!
===================
Success       : ${workflow.success}
Exit status   : ${workflow.exitStatus}
Duration      : ${workflow.duration}
Work dir      : ${workflow.workDir}
Results dir   : ${params.output}
"""
}

/*
 * Error handling
 */
workflow.onError {
    log.error """\
Pipeline failed!
================
Error message : ${workflow.errorMessage}
Exit status   : ${workflow.exitStatus}
Work dir      : ${workflow.workDir}
"""
}
