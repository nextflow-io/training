#!/usr/bin/env nextflow

nextflow.enable.dsl=2

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

params.greeting  = 'Hello world!'
greeting_ch = Channel.from(params.greeting)

workflow {

    letters_ch = splitLetters(greeting_ch)
    uppercase_ch = convertToUpper( letters_ch.flatten() )
    uppercase_ch.view{ it.trim() }

}
