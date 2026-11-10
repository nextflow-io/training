#!/usr/bin/env nextflow

// Include modules
include { sayHello } from './modules/sayHello.nf'
include { translateGreeting } from './modules/translateGreeting.nf'
include { collectGreetings } from './modules/collectGreetings.nf'
include { translateGreeting as translateSecond } from './modules/translateGreeting.nf'
include { collectGreetings as collectSecond } from './modules/collectGreetings.nf'
include { cowpy } from './modules/cowpy.nf'

/*
 * Pipeline parameters
 */
params {
    input: Path
    language: String
    second_language: String
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
    collectGreetings(translateGreeting.out.collect(), params.batch, params.language)
    // translate the greeting into the second language
    translateSecond(sayHello.out, params.second_language)
    // collect the second-language greetings into one file
    collectSecond(translateSecond.out.collect(), params.batch, params.second_language)
    // generate ASCII art of the greetings with cowpy
    cowpy(collectGreetings.out.outfile, params.character)

    publish:
    greetings = sayHello.out
    translated = translateGreeting.out
    collected = collectGreetings.out.outfile
    batch_report = collectGreetings.out.report
    cowpy_art = cowpy.out
    translated_second = translateSecond.out
    collected_second = collectSecond.out.outfile
    batch_report_second = collectSecond.out.report
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
    translated_second {
        path "${params.batch}/intermediates"
        mode 'copy'
    }
    collected_second {
        path "${params.batch}/intermediates"
        mode 'copy'
    }
    batch_report_second {
        path "${params.batch}"
        mode 'copy'
    }
}
