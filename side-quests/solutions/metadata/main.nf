/*
    * Use langid to predict the language of each input file
    */
process IDENTIFY_LANGUAGE {

    container 'community.wave.seqera.io/library/pip_langid:b2269f456a5629ff'

    input:
    tuple val(meta), path(file)

    output:
    tuple val(meta), path(file), stdout

    script:
    """
    langid < ${file} -l en,de,fr,es,it | sed -E "s/.*\\('([a-z]+)'.*/\\1/" | tr -d '\\n'
    """
}

/*
    * Generate ASCII art with cowpy
*/
process COWPY {

    publishDir "results/${meta.lang_group}", mode: 'copy'

    container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'

    input:
    tuple val(meta), path(input_file)

    output:
    tuple val(meta), path("cowpy-${input_file}")

    script:
    """
    cat ${input_file} | cowpy -c ${meta.character} > cowpy-${input_file}
    """
}

workflow {

    ch_samplesheet = Channel.fromPath("./data/samplesheet.csv")
        .splitCsv(header: true)
        .map { row ->
            [[id: row.id, character: row.character], row.recording]
        }

    ch_prediction = IDENTIFY_LANGUAGE(ch_samplesheet)
    ch_languages = ch_prediction
        .map { meta, file, lang ->
            [meta + [lang: lang], file]
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
    COWPY(ch_languages)
}
