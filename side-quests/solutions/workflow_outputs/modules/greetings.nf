#!/usr/bin/env nextflow

/*
 * Create a greeting file from input text
 * Note: No publishDir - outputs managed by workflow output block
 */
process SAY_HELLO {

    input:
        tuple val(greeting), val(language)

    output:
        tuple val(greeting), val(language), path("${greeting}-output.txt")

    script:
    """
    echo '$greeting' > '${greeting}-output.txt'
    """
}

/*
 * Convert greeting to uppercase
 * Note: No publishDir - outputs managed by workflow output block
 */
process CONVERT_TO_UPPER {

    input:
        tuple val(greeting), val(language), path(input_file)

    output:
        tuple val(greeting), val(language), path("UPPER-${input_file}")

    script:
    """
    cat '$input_file' | tr '[a-z]' '[A-Z]' > 'UPPER-${input_file}'
    """
}
