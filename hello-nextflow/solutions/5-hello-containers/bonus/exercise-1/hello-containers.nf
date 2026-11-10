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
    input: Path
    language: String
    batch: String
    character: String
}

workflow {

    main:
    // create a channel for inputs from a CSV file
    name_ch = channel.fromPath(params.input)
                        .splitCsv()
                        .map { row -> row[0] }
    // emit a greeting
    sayHello(name_ch)
    // translate the greeting
    translateGreeting(sayHello.out, params.language)
    // collect all the greetings into one file
    collectGreetings(translateGreeting.out.collect(), params.batch)
    // generate ASCII art of the greetings with cowpy
    cowpy(collectGreetings.out.outfile, params.character)

    publish:
    greetings = sayHello.out
    translated = translateGreeting.out
    collected = collectGreetings.out.outfile
    batch_report = collectGreetings.out.report
    cowpy_art = cowpy.out
}

output {
    greetings {
        path "${params.batch}/intermediates"
    }
    translated {
        path "${params.batch}/intermediates"
    }
    collected {
        path "${params.batch}/intermediates"
    }
    batch_report {
        path "${params.batch}"
    }
    cowpy_art {
        path "${params.batch}"
        mode 'copy'
    }
}
