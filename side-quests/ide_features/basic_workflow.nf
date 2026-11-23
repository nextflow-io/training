#!/usr/bin/env nextflow

process FASTQC {
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("*.html"), emit: html
    tuple val(sample_id), path("*.zip"), emit: zip

    script:
    def args = task.ext.args ?: ''
    """
    fastqc \\
        ${args} \\
        --threads ${task.cpus} \\
        ${reads}
    """
}

/*
 * Workflow parameters
 */
params{
    input: Path = 'data/sample_data.csv'
    output_dir: String = 'results'
}

/*
 * Main workflow demonstrating module usage and navigation
 */
workflow {
    main:
    // Create input channel from CSV file
    ch_input = channel.fromPath(params.input, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            return [row.sample_id, file(row.fastq_path)]
        }
        .view { sample_id, fastq -> "Processing sample: ${sample_id} -> ${fastq}" }

    // Run FastQC process - Ctrl-click on FASTQC to navigate to module
    FASTQC(ch_input)

    // View outputs
    FASTQC.out.html.view { report -> "FastQC report: ${report}" }
    FASTQC.out.zip.view { data -> "FastQC data: ${data}" }

    publish:
    fastqc_html = FASTQC.out.html
    fastqc_zip = FASTQC.out.zip
}

output {
    fastqc_html {
        path 'fastqc'
    }
    fastqc_zip {
        path 'fastqc'
    }
}
