#!/usr/bin/env nextflow

include { FASTQC } from './modules/fastqc.nf'

/*
 * Workflow parameters
 */
params.input = 'data/sample_data.csv'
params.output_dir = 'results'

/*
 * Main workflow demonstrating module usage and navigation
 */
workflow {
    // Create input channel from CSV file
    ch_input = Channel.fromPath(params.input, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            return [row.sample_id, file(row.fastq_path)]
        }
        .view { "Processing sample: ${it[0]} -> ${it[1]}" }

    // Run FastQC process - Ctrl-click on FASTQC to navigate to module
    FASTQC(ch_input)

    // View outputs
    FASTQC.out.html.view { "FastQC report: ${it}" }
    FASTQC.out.zip.view { "FastQC data: ${it}" }
}
