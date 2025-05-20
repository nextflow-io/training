/*
 * Use echo to print 'Hello World!' to a file
 */
process IDENTIFY_LANGUAGE {
    publishDir 'results', mode: 'copy'
    tag "${meta.id}"

    container 'community.wave.seqera.io/library/pip_langid:b2269f456a5629ff'

    input:
    tuple val(meta), path(greeting)

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
    tag "${meta.id}"

    publishDir "results/${meta.lang}", mode: 'copy'

    container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'

    input:
    tuple val(meta), path(input_file)
    val character

    output:
    tuple val(meta), path("cowpy-${input_file}")

    script:
    """
    cat $input_file | cowpy -c "$character" > cowpy-${input_file}
    """

}

workflow  {

    files = Channel.fromPath("./data/*.txt").map { file -> [ [id:file.getName()], file] }

    ch_prediction = IDENTIFY_LANGUAGE(files)

    ch_language_groups = files.join(ch_prediction)
                        //Uses meta map to join
                        .map { meta, file, lang ->
                            [ meta + [lang:lang], file ]
                        }
                        //generates a new key in the meta map
                        .map{ meta, file ->
                            def lang_group = (meta.lang.equals('de') || meta.lang.equals('en')) ? 'germanic' : 'romanic'
                            [ meta + [lang_group:lang_group], file ]
                        }

    //uses meta map field to filter the channel to only include languages in the romanic group
    romanic_languages = ch_language_groups.filter { meta, file ->
                                                    meta.lang_group == 'romanic'
                                                }

    COWPY(romanic_languages, params.character)

}
