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
    greetings_ch = sayHello(name_ch)
    // translate the greeting
    translated_ch = translateGreeting(greetings_ch, params.language)
    // collect all the greetings into one file
    collected_out = collectGreetings(translated_ch.collect(), params.batch, params.language)
    // translate the greeting into the second language
    translated_second_ch = translateSecond(greetings_ch, params.second_language)
    // collect the second-language greetings into one file
    collected_second_out = collectSecond(translated_second_ch.collect(), params.batch, params.second_language)
    // generate ASCII art of the greetings with cowpy
    cowpy_ch = cowpy(collected_out.outfile, params.character)

    publish:
    greetings = greetings_ch
    translated = translated_ch
    collected = collected_out.outfile
    batch_report = collected_out.report
    cowpy_art = cowpy_ch
    translated_second = translated_second_ch
    collected_second = collected_second_out.outfile
    batch_report_second = collected_second_out.report
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
