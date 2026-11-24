#!/usr/bin/env nextflow

// Include modules
include { sayHello } from './modules/sayHello.nf'
include { translateGreeting } from './modules/translateGreeting.nf'
include { collectGreetings } from './modules/collectGreetings.nf'
include { cowpy } from './modules/cowpy.nf'

/*
 * Pipeline parameters
 */
params {
    input: Path = 'data/names.csv'
    language: String = 'es'
    batch: String = 'batch'
    character: String = 'turkey'
}

workflow {

    main:
    // create a channel for inputs from a CSV file
    name_ch = channel.fromPath(params.input)
                        .splitCsv()
                        .map { row -> row[0] }
    // emit a greeting
    greetings_ch = sayHello(name_ch)
    // translate the greeting
    translated_ch = translateGreeting(greetings_ch, params.language)
    // collect all the greetings into one file
    collected_out = collectGreetings(translated_ch.collect(), params.batch, params.language)
    // generate ASCII art of the greetings with cowpy
    cowpy_ch = cowpy(collected_out.outfile, params.character)

    publish:
    greetings = greetings_ch
    translated = translated_ch
    collected = collected_out.outfile
    batch_report = collected_out.report
    cowpy_art = cowpy_ch
}

output {
    greetings {
        path "${params.batch}/intermediates"
        mode 'copy'
    }
    translated {
        path "${params.batch}/intermediates"
        mode 'copy'
    }
    collected {
        path "${params.batch}/intermediates"
        mode 'copy'
    }
    batch_report {
        path "${params.batch}"
        mode 'copy'
    }
    cowpy_art {
        path "${params.batch}"
        mode 'copy'
    }
}
