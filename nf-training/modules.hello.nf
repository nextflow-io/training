#!/usr/bin/env nextflow

params.greeting  = 'Hello world!'
greeting_ch = channel.of(params.greeting)

include { SPLITLETTERS   } from './modules.nf'
include { CONVERTTOUPPER } from './modules.nf'

workflow{
    letters_ch = SPLITLETTERS(greeting_ch)
    results_ch = CONVERTTOUPPER(letters_ch.flatten())
    results_ch.view{ it }
}
