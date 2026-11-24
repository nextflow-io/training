#!/usr/bin/env nextflow

// Include the greeting workflow and the cowpy module
include { GREET } from './workflows/greet.nf'
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
    // greet, translate and collect
    greet_out = GREET(name_ch, params.language, params.batch)
    // generate ASCII art of the greetings with cowpy
    cowpy_ch = cowpy(greet_out.collected, params.character)

    publish:
    greetings = greet_out.greetings
    translated = greet_out.translated
    collected = greet_out.collected
    batch_report = greet_out.report
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
