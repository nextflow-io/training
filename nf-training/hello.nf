#!/usr/bin/env nextflow

params.greeting  = 'Hello world!'
greeting_ch = Channel.from(params.greeting)

process splitLetters {

    input:
    val x

    output:
    file 'chunk_*'

    """
    printf '$x' | split -b 6 - chunk_
    """
}

process convertToUpper {

    input:
    file y

    output:
    stdout

    """
    cat $y | tr '[a-z]' '[A-Z]' 
    """
}

workflow{
    letters_ch = splitLetters(greeting_ch)
    results_ch = convertToUpper(letters_ch.flatten())
    results_ch.view{ it }
}


