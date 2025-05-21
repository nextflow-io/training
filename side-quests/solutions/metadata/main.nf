/*
 * Use langid to predict the language of each input file
 */
process IDENTIFY_LANGUAGE {

    container 'community.wave.seqera.io/library/pip_langid:b2269f456a5629ff'

    input:
    tuple val(meta), path(file)

    output:
    tuple val(meta), stdout

    script:
    """
    langid < ${greeting} -l en,de,fr,es,it | sed -E "s/.*\\('([a-z]+)'.*/\\1/" | tr -d '\\n'
    """

}

/*
 * Generate ASCII art with cowpy
*/
process COWPY {

    publishDir "results/${meta.lang}", mode: 'copy'

    container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'

    input:
    tuple val(meta), path(file)

    output:
    tuple val(meta), path("cowpy-${file}")

    script:
    """
    cat $file | cowpy -c ${meta.character} > cowpy-${file}
    """

}

workflow  {

    ch_samplesheet = Channel.fromPath("./data/samplesheet.csv")
                    .splitCsv(header: true)
                    .map { row ->
                        [ [id:row.id, character:row.character], row.recording ]
                    }

    ch_prediction = IDENTIFY_LANGUAGE(ch_samplesheet)

    ch_languages = ch_samplesheet.join(ch_prediction)
                                .map { meta, file, lang ->
                                    [ meta + [lang:lang], file ]
                                }
                                .map{ meta, file ->
                                    def lang_group = (meta.lang.equals('de') || meta.lang.equals('en')) ? 'germanic' : 'romanic'
                                    [ meta + [lang_group:lang_group], file ]
                                }

    romanic_languages = ch_languages.filter { meta, file ->
                                        meta.lang_group == 'romanic'
                                    }

    COWPY(romanic_languages)

}
