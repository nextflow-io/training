#!/usr/bin/env nextflow

include { COWPY } from './modules/cowpy.nf'
include { IDENTIFY_LANGUAGE } from './modules/langid.nf'

workflow  {
    main:
    ch_datasheet = channel.fromPath("./data/datasheet.csv")
        .splitCsv(header: true)
        .map { row ->
            [[id: row.id, character: row.character], row.recording]
        }

    COWPY(ch_datasheet)

    // Run langid to identify the language of each greeting
    IDENTIFY_LANGUAGE(ch_datasheet)
    IDENTIFY_LANGUAGE.out
        .map { meta, file, lang_id ->
            [meta + [lang: lang_id], file]
        }
        .view()

    publish:
    cowpy_art = COWPY.out
}

output {
    cowpy_art {
    }
}
