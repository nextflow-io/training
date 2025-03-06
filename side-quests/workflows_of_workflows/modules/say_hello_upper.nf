#!/usr/bin/env nextflow

/*
 * Use a text replacement tool to convert the greeting to uppercase
 */
process SAY_HELLO_UPPER {
    publishDir 'results', mode: 'copy'

    tag "converting ${input_file}"

    input:
        path input_file

    output:
        path "UPPER-${input_file}"

    script:
    """
    cat '${input_file}' | tr '[a-z]' '[A-Z]' > 'UPPER-${input_file}'
    """
}
