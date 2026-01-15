#!/usr/bin/env nextflow

// Include modules
include { sayHello } from './modules/sayHello.nf'
include { convertToUpper } from './modules/convertToUpper.nf'
include { collectGreetings } from './modules/collectGreetings.nf'
include { cowpy } from './modules/cowpy.nf'

/*
 * Pipeline parameters
 */
params {
    input: Path = 'data/greetings.csv'
    batch: String = 'test-batch'
    character: String = 'turkey'
}

workflow {

    main:
    // create a channel for inputs from a CSV file
    greeting_ch = channel.fromPath(params.input)
                        .splitCsv()
                        .map { line -> line[0] }
    // emit a greeting
    sayHello(greeting_ch)
    // convert the greeting to uppercase
    convertToUpper(sayHello.out)
    // collect all the greetings into one file
    collectGreetings(convertToUpper.out.collect(), params.batch)
    // generate ASCII art of the greetings with cowpy
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
        path 'hello_containers/intermediates'
        mode 'copy'
    }
    uppercased {
        path 'hello_containers/intermediates'
        mode 'copy'
    }
    collected {
        path 'hello_containers/intermediates'
        mode 'copy'
    }
    batch_report {
        path 'hello_containers'
        mode 'copy'
    }
    cowpy_art {
        path 'hello_containers'
        mode 'copy'
    }
}
