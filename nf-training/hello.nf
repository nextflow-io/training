#!/usr/bin/env nextflow

params.greeting  = 'Hello world!'
greeting_ch = Channel.from(params.greeting)

process splitLetters {

    input:
    val x from greeting_ch

    output:
    file 'chunk_*' into letters_ch

    """
    printf '$x' | split -b 6 - chunk_
    """
}

process convertToUpper {

    input:
    file y from letters_ch.flatten()

    output:
    stdout into result_ch

    """
    cat $y | tr '[a-z]' '[A-Z]' 
    """
}

result_ch.view{ it }
