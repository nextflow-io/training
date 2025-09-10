#!/usr/bin/env nextflow

// Generate ASCII art with cowpy
process cowpy {

    container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'

    publishDir 'results', mode: 'copy'

    input:
        path input_file
        val character

    output:
        path "cowpy-${input_file}"

    script:
    """
    cat $input_file | cowpy -c "$character" > cowpy-${input_file}
    """

}
