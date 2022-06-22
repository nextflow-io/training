#!/usr/bin/env nextflow

params.greeting = 'Hello world!'
greeting_ch = Channel.from(params.greeting)

process SPLITLETTERS {
    input:
    val x

    output:
    file 'chunk_*'

    """
    printf '$x' | split -b 6 - chunk_
    """
}

process CONVERTTOUPPER {
    input:
    file y

    output:
    stdout

    """
    cat $y | tr '[a-z]' '[A-Z]' 
    """
}

workflow {
    letters_ch = SPLITLETTERS(greeting_ch)
    results_ch = CONVERTTOUPPER(letters_ch.flatten())
    results_ch.view{ it }
}

