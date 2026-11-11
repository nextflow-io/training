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
    name: String = 'Friend'
}

workflow {

    main:
    // declare an array of input names
    names_array = ['Mateo','Raphael','Celine']
    // create a channel for inputs
    name_ch = channel.of(names_array)
                        .view { names -> "Before flatten: $names" }
                        .flatten()
                        .view { name -> "After flatten: $name" }
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
