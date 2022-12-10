nextflow.preview.dsl=2

include{ INDEX; QUANT; FASTQC; MULTIQC } from './tasks.nf'

workflow RNASEQ {
  take:
    transcriptome
    read_pairs_ch

  main:
    INDEX( transcriptome )
    QUANT( INDEX.out, read_pairs_ch )
    FASTQC( read_pairs_ch )
    MULTIQC(
            QUANT.out.mix(FASTQC.out).collect(),
            file(params.multiqc) )
}