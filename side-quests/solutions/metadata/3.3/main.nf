#!/usr/bin/env nextflow

include { IDENTIFY_LANGUAGE } from './modules/langid.nf'
include { COWPY } from './modules/cowpy.nf'

workflow {

    ch_datasheet = channel.fromPath("./data/datasheet.csv")
        .splitCsv(header: true)
        .map { row ->
            [[id: row.id, character: row.character], row.recording]
        }

    // Run langid to identify the language of each greeting
    IDENTIFY_LANGUAGE(ch_datasheet)

    // Reorganize and enrich the metadata
    IDENTIFY_LANGUAGE.out
        .map { meta, file, lang_id ->
            [meta + [lang: lang_id], file]
        }
        .map { meta, file ->

            def lang_group = "unknown"
            if (meta.lang.equals("de") || meta.lang.equals('en')) {
                lang_group = "germanic"
            }
            else if (meta.lang in ["fr", "es", "it"]) {
                lang_group = "romance"
            }

            [meta + [lang_group: lang_group], file]
        }
        .set { ch_languages }

    // Run cowpy to generate ASCII art
    COWPY(ch_languages)
}
