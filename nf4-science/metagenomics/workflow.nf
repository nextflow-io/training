/*
 * required tasks
 */
include { BOWTIE2                     }     from './modules/bowtie2.nf'
include { KRAKEN2                     }    from './modules/kraken2.nf'
include { BRACKEN                     }     from './modules/bracken.nf'
include { K_REPORT_TO_KRONA         }     from './modules/kReport2Krona.nf'
include { KT_IMPORT_TEXT             }     from './modules/ktImportText.nf'
include { KRAKEN_BIOM               }   from './modules/kraken_biom.nf'
include { KNIT_PHYLOSEQ             }   from './modules/knit_phyloseq.nf'

/*
 * workflow
 */

workflow kraken2Flow {
    // required inputs
    take:
        bowtie2_index
        kraken2_db
        reads_ch
    // workflow implementation
    main:
        BOWTIE2(reads_ch, bowtie2_index)
        KRAKEN2(BOWTIE2.out, kraken2_db)
        BRACKEN(KRAKEN2.out, kraken2_db)
        K_REPORT_TO_KRONA(BRACKEN.out)
        KT_IMPORT_TEXT(K_REPORT_TO_KRONA.out)
        if(params.sheet_csv){
            KRAKEN_BIOM(BRACKEN.out.collect())
            KNIT_PHYLOSEQ(KRAKEN_BIOM.out)
        }
}
