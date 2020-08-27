nextflow.preview.dsl=2

/*
 * define main params
 */
params.reads = "data/ggal/gut_{1,2}.fq"
params.trans = "data/ggal/transcriptome_*.fa"
params.multiqc = "data/multiqc"
params.outdir = "results"

/*
 * include pipeline and helper components
 */
include{ RNASEQ } from './modules/flow.nf'
include{ rna_inputs } from './modules/helpers.nf'

/*
 * main flow
 */
workflow {
    rna_inputs(params.trans, params.reads) | RNASEQ
}
