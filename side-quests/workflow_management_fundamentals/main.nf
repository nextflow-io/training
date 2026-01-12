#!/usr/bin/env nextflow

// Include process definitions
include { FASTQC } from './modules/fastqc'

// Pipeline parameters
params.samples = 'data/samples.csv'
params.outdir = 'results'

// Main workflow
workflow {

    // Read sample sheet and create channel
    ch_samples = Channel
        .fromPath(params.samples)
        .splitCsv(header: true)
        .map { row ->
            def meta = [id: row.sample_id, organism: row.organism]
            def reads = [file(row.read1), file(row.read2)]
            return [meta, reads]
        }

    // Run FastQC on all samples
    FASTQC(ch_samples)
}
