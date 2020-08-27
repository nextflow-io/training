nextflow.preview.dsl=2

/*
 * define main params
 */
params.reads = "data/ggal/gut_{1,2}.fq"
params.trans = "data/ggal/transcriptome_*.fa"
params.multiqc = "data/multiqc"
params.outdir = "results"

/*
 * include pipeline logic
 */
include{ RNASEQ } from './modules/flow.nf'

/*
 * For the sake of the example we run each transcript file against all reads
 */
workflow {
  // basic channels
  trans = Channel .fromPath( params.trans )
  reads = Channel .fromFilePairs( params.reads, checkExists: true )

  // create a multi-channel combining all transcripts against all read pairs
  trans
    .combine(reads)
    .dump(tag:'combination')
    .multiMap { row -> trans: row[0]; reads: tuple(row[1],row[2]) }
    .set { trans_and_reads }

  // call pipeline invocation, note: one argument holding two channels
  RNASEQ(trans_and_reads)
}