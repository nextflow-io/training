#!/usr/bin/env nextflow

include { sayHello } from './modules/sayHello.nf'

/*
 * Pipeline parameters
 */
params {
    input: String
}

workflow {

    main:
    // emit a greeting
    sayHello(params.input)

    publish:
    first_output = sayHello.out
}

output {
    first_output {
        path '1-hello'
        mode 'copy'
    }
}
