#!/usr/bin/env nextflow

/*
 * Bacterial Genome Analysis Pipeline - Nextflow Version
 *
 * This workflow demonstrates the benefits of workflow management:
 * - Automatic parallelization based on data dependencies
 * - Built-in resume/caching after failures
 * - Per-process container isolation
 * - Declarative resource management
 * - Automatic provenance tracking
 */

// Include process definitions
include { FASTQC } from './modules/fastqc'
include { FASTP } from './modules/fastp'
include { SPADES } from './modules/spades'
include { QUAST } from './modules/quast'

// Pipeline parameters
params.samples = '../data/samples.csv'
params.outdir = '../results'

// Main workflow
workflow {

    // Read sample sheet and create channel
    // Each sample becomes an independent item that can be processed in parallel
    ch_samples = Channel
        .fromPath(params.samples)
        .splitCsv(header: true)
        .map { row ->
            def meta = [id: row.sample_id, organism: row.organism]
            def reads = [file(row.read1), file(row.read2)]
            return [meta, reads]
        }

    // Run the pipeline
    // Nextflow automatically determines what can run in parallel

    // FastQC and FASTP both need only the raw reads - they run in PARALLEL
    FASTQC(ch_samples)
    FASTP(ch_samples)

    // SPADES needs FASTP output - it WAITS for FASTP, but runs samples in PARALLEL
    SPADES(FASTP.out.reads)

    // QUAST needs SPADES output - it WAITS for SPADES, but runs samples in PARALLEL
    QUAST(SPADES.out.assembly)
}
