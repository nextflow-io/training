#!/usr/bin/env nextflow

include { sayHello } from './modules/sayHello.nf'
include { convertToUpper } from './modules/convertToUpper.nf'
include { collectGreetings } from './modules/collectGreetings.nf'
include { cowpy } from './modules/cowpy.nf'

/*
 * Pipeline parameters
 */
params {
    input: Path
    batch: String
    character: String
}

workflow {

    main:
    // create a channel for inputs from a CSV file
    greeting_ch = channel.fromPath(params.input)
                        .splitCsv()
                        .map { line -> line[0] }
    sayHello(greeting_ch)
    convertToUpper(sayHello.out)
    collectGreetings(convertToUpper.out.collect(), params.batch)
    cowpy(collectGreetings.out.outfile, params.character)

    publish:
    first_output = sayHello.out
    uppercased = convertToUpper.out
    collected = collectGreetings.out.outfile
    batch_report = collectGreetings.out.report
    cowpy_art = cowpy.out
}

output {
    first_output {
        path { sayHello.name }
    }
    uppercased {
        path { convertToUpper.name }
    }
    collected {
        path { collectGreetings.name }
    }
    batch_report {
        path { collectGreetings.name }
    }
    cowpy_art {
        path { cowpy.name }
    }
}
