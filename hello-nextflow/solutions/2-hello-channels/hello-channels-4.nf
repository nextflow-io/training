#!/usr/bin/env nextflow

/*
 * Use echo to print 'Hello World!' to a file
 */
process sayHello {

    input:
    val greeting

    output:
    path "${greeting}-output.txt"

    script:
    """
    echo '${greeting}' > '${greeting}-output.txt'
    """
}

/*
 * Pipeline parameters
 */
params {
    input: Path = 'greetings.csv'
}

workflow {

    main:
    greetings_array = ['Hello', 'Bonjour', 'Holà']

    // create a channel for inputs from a CSV file
    greeting_ch = channel.fromPath(params.greeting)
        .view { item -> "Before splitCsv: ${item}" }
        .splitCsv()
        .view { item -> "After splitCsv: ${item}" }
        .map { line -> line[0] }
        .view { item -> "After map: ${item}" }

    // emit a greeting
    sayHello(greeting_ch)

    publish:
    first_output = sayHello.out
}

output {
    first_output {
        path 'hello_channels'
        mode 'copy'
    }
}
