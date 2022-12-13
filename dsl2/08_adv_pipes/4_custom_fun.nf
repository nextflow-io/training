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
 * custom function to split the transcript and reads channel
 */
def split_trans_and_reads(combined) {
  combined.multiMap { row ->
            trans: row[0];
            reads: tuple(row[1],row[2])
            }
}

/*
 * main flow
 */
workflow {

  // basic channels
  trans = Channel .fromPath( params.trans )
  reads = Channel .fromFilePairs( params.reads, checkExists: true )

  // assemble the flow logic
  trans \
    | combine(reads) \
    | split_trans_and_reads \
    | RNASEQ
}