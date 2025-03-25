#!/usr/bin/env nextflow

include { SPLITLETTERS } from './modules.hi.nf'
include { CONVERTTOUPPER } from './modules.hi.nf'
include { COUNTLETTERS } from './modules.hi.nf'

//params.greeting = 'Hello world!'
//greeting_ch = Channel.of(params.greeting)
greeting_ch = Channel.fromPath('nf-training/saudacoes.txt').splitText()

workflow {
    letters_ch = SPLITLETTERS(greeting_ch)
    count_ch = COUNTLETTERS(greeting_ch)
    results_ch = CONVERTTOUPPER(letters_ch.flatten())
    results_ch.view{ it }
    //count_ch.view{ it }
}
