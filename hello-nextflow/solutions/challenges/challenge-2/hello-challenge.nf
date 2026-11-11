#!/usr/bin/env nextflow

/*
 * Use echo to print a greeting to a file
 */
process sayHello {

    tag "${name}"

    input:
    val name

    output:
    path "${name}.txt"

    script:
    """
    echo 'Hello, ${name}!' > '${name}.txt'
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
                        .filter { row -> row[1] == 'fr' }
                        .map { row -> row[0] }
                        .view { name -> "After map: $name" }
    // emit a greeting
    sayHello(name_ch)

    publish:
    greetings = sayHello.out
}

output {
    greetings {
        path 'hello'
        mode 'copy'
    }
}
