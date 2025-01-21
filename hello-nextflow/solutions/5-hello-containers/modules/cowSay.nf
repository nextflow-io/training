#!/usr/bin/env nextflow

// Generate ASCII art with cowsay
process cowSay {

    publishDir 'results', mode: 'copy'

    container 'community.wave.seqera.io/library/pip_cowsay:131d6a1b707a8e65'

    input:
        path input_file
        val character

    output:
        path "cowsay-${input_file}"

    script:
    """
    cowsay -c "$character" -t "\$(cat $input_file)" > cowsay-${input_file}
    """
}
