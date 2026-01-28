#!/usr/bin/env nextflow

/*
 * RNA-seq Analysis Pipeline - Nextflow Version
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
include { UNTAR; SALMON_QUANT } from './modules/salmon'

// Pipeline parameters
params.samples = '../data/samples.csv'
params.salmon_index = 'https://raw.githubusercontent.com/nf-core/test-datasets/rnaseq/reference/salmon.tar.gz'
params.outdir = '../results'

// Main workflow
workflow {

    // Read sample sheet and create channel
    // Each sample becomes an independent item that can be processed in parallel
    ch_samples = Channel
        .fromPath(params.samples)
        .splitCsv(header: true)
        .map { row ->
            def meta = [id: row.sample]
            def reads = [file(row.fastq_1), file(row.fastq_2)]
            return [meta, reads]
        }

    // Download and extract pre-built salmon index
    ch_salmon_index = Channel.fromPath(params.salmon_index)
    UNTAR(ch_salmon_index)

    // Run the pipeline
    // Nextflow automatically determines what can run in parallel

    // FastQC and FASTP both need only the raw reads - they run in PARALLEL
    FASTQC(ch_samples)
    FASTP(ch_samples)

    // SALMON needs FASTP output - it WAITS for FASTP, but runs samples in PARALLEL
    // .first() converts the index to a value channel so it's reused for all samples
    SALMON_QUANT(FASTP.out.reads, UNTAR.out.index.first())
}
