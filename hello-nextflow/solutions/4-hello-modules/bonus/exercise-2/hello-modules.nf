#!/usr/bin/env nextflow

// Include modules
include { sayHello } from './modules/sayHello/main.nf'
include { translateGreeting } from './modules/translateGreeting/main.nf'
include { collectGreetings } from './modules/collectGreetings/main.nf'

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

    publish:
    greetings = sayHello.out
    translated = translateGreeting.out
    collected = collectGreetings.out.outfile
    batch_report = collectGreetings.out.report
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
}
