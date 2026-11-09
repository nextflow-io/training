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
    input: Path = 'data/names.csv'
}

workflow {

    main:
    // create a channel for inputs from a CSV file
    name_ch = channel.fromPath(params.input)
                        .view { csv -> "Before splitCsv: $csv" }
                        .splitCsv()
                        .view { row -> "After splitCsv: $row" }
                        .filter { row -> row[1] == 'french' }
                        .map { row -> row[0] }
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
