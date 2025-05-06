#!/usr/bin/env nextflow

process convertToUpper {

    publishDir 'results', mode: 'copy'

    input:
    path input_file

    output:
    path "UPPER-${input_file}"

    script: 
    """
    cat $input_file | tr '[a-z]' '[A-Z]'> UPPER-${input_file}
    """
}