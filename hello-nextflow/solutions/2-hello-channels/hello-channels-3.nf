#!/usr/bin/env nextflow

/*
 * Use echo to print 'Hello World!' to a file
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
    name: String = 'World'
}

workflow {

    main:
    // declare an array of input names
    names_array = ['Mateo','Andrea','Celine']
    // create a channel for inputs
    name_ch = channel.of(names_array)
                        .view { names -> "Before flatten: $names" }
                        .flatten()
                        .view { name -> "After flatten: $name" }
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
