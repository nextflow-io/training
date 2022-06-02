#!/usr/bin/env nextflow

params.greeting  = 'Hello world!'
greeting_ch = Channel.from(params.greeting)

include { splitLetters   } from './modules.nf'
include { convertToUpper } from './modules.nf'

workflow{
    letters_ch = splitLetters(greeting_ch)
    results_ch = convertToUpper(letters_ch.flatten())
    results_ch.view{ it }
}

