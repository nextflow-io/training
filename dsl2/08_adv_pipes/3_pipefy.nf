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
 * Same as before using piping operators
 */
workflow {
  // basic channels
  trans = Channel .fromPath( params.trans )
  reads = Channel .fromFilePairs( params.reads, checkExists: true )

  // create a multi-channel combining all transcripts against all read pairs
  trans \
    | combine(reads) \
    | multiMap { row -> trans: row[0]; reads: tuple(row[1],row[2]) } \
    | RNASEQ
}