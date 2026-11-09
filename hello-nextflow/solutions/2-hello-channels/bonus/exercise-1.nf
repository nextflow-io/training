#!/usr/bin/env nextflow

/*
 * Use echo to print a greeting to a file
 */
process sayHello {

    input:
    val name

    output:
    path "${name}-output.txt"

    script:
    """
    echo 'Hello, ${name}!' > '${name}-output.txt'
    """
}

/*
 * Pipeline parameters
 */
params {
    input: Path = 'data/names-header.csv'
}

workflow {

    main:
    // create a channel for inputs from a CSV file
    name_ch = channel.fromPath(params.input)
                        .view { csv -> "Before splitCsv: $csv" }
                        .splitCsv(header: true)
                        .view { row -> "After splitCsv: $row" }
                        .map { row -> row.name }
                        .view { name -> "After map: $name" }
    // emit a greeting
    sayHello(name_ch)

    publish:
    first_output = sayHello.out
}

output {
    first_output {
        path '.'
    }
}
