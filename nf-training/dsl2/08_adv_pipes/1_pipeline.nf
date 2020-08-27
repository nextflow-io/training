nextflow.preview.dsl=2

/*
 * define main params
 */
params.reads = "data/ggal/gut_{1,2}.fq"
params.trans = "data/ggal/transcript.fa"
params.multiqc = "data/multiqc"
params.outdir = "results"

/*
 * include pipeline logic
 */
include{ RNASEQ } from './modules/flow.nf'

/*
 * define params and run it
 */
workflow {
  trans = file(params.trans)
  reads = Channel .fromFilePairs( params.reads, checkExists: true )

  RNASEQ(trans, reads)
}