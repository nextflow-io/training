#!/usr/bin/env nextflow

include { sayHello } from './modules/sayHello.nf'

/*
 * Pipeline parameters
 */
params {
    input: Path
}

workflow {

    main:
    // create a channel for inputs from a CSV file
    greeting_ch = channel.fromPath(params.input)
                        .splitCsv()
                        .map { line -> line[0] }
    // emit a greeting
    sayHello(greeting_ch)

    publish:
    first_output = sayHello.out
}

output {
    first_output {
        path '2-inputs'
        mode 'copy'
    }
}
