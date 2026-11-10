#!/usr/bin/env nextflow

// Include modules
include { sayHello } from './challenge-modules/sayHello.nf'
include { translateGreeting } from './challenge-modules/translateGreeting.nf'
include { collectGreetings } from './challenge-modules/collectGreetings.nf'
include { cowpy } from './challenge-modules/cowpy.nf'

/*
 * Pipeline parameters
 */
params {
    input: Path
    language: String
    batch: String
}

workflow {

    main:
    // create a channel of [name, character] pairs from a CSV file
    person_ch = channel.fromPath(params.input)
                        .splitCsv()
                        .map { row -> [row[0], row[2]] }
    // emit a greeting
    sayHello(person_ch)
    // translate the greeting
    translateGreeting(sayHello.out, params.language)
    // collect all the greetings into one file
    collectGreetings(translateGreeting.out.map { greeting_file, _character -> greeting_file }.collect(), params.batch, params.language)
    // generate ASCII art of each greeting with its own character
    cowpy(translateGreeting.out)

    publish:
    greetings = sayHello.out.map { greeting_file, _character -> greeting_file }
    translated = translateGreeting.out.map { greeting_file, _character -> greeting_file }
    collected = collectGreetings.out.outfile
    batch_report = collectGreetings.out.report
    cowpy_art = cowpy.out
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
