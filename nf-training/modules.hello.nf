#!/usr/bin/env nextflow

params.greeting  = 'Hello world!'

include { SPLITLETTERS   } from './modules.nf'
include { CONVERTTOUPPER } from './modules.nf'

workflow{
    greeting_ch = Channel.of(params.greeting)
    letters_ch = SPLITLETTERS(greeting_ch)
    results_ch = CONVERTTOUPPER(letters_ch.flatten())
    results_ch.view{ it }
}
