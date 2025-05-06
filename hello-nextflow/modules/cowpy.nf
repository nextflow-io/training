#!/usr/bin/env nextflow

// Generate ASCII art using cowsay
process cowpy {

    publishDir 'results', mode: 'copy'
    container 'community.wave.seqera.io/library/pip_cowpy:85c5ef8d359fd9e6'
    conda 'conda-forge::cowpy='

    input:
    path input_file
    val character

    output:
    path "cowpy-${input_file}"

    script:
    """
    cat $input_file | cowpy -c "$character" > "cowpy-${input_file}"
    """
    
}