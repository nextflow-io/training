#!/usr/bin/env nextflow

/*
 * Use echo to print 'Hello World!' to a file
 */
process sayHello {

    publishDir 'results', mode: 'copy'

    input:
        val greeting

    output:
        path "$greeting-output.txt" //requires double quotes to eval the variable

    script:
    """
    echo '$greeting' > '$greeting-output.txt'
    """
}

/*
 * Pipeline parameters
 */
params.greeting = 'greetings.csv'   // 'Holà mundo!'

workflow {    //note nextflow does parallelization of processes - processing order of dataflow logic can be random

    greeting_ch = channel.fromPath(params.greeting) // read the file from path
                         .view { csv -> "before splitCSV $csv"}  // .view operator to see the content of the channel
                         .splitCsv()
                         .map( item -> item[0] ) // map operator to extract the first column of the CSV file
                         .view { csv -> "after splitCSV $csv"}  // .view operator to see the content of the channel

    // greetings_array = ['Hello, world!', 'Bonjour le monde!', 'Hola mundo!', 'Hallo Welt!', 'Ciao!']
    // create a channel from the array - the channel can include multiple operators chained behind each other
    // greeting_ch = Channel.of(greetings_array)
                        //.view { greeting -> "before flatten: $greeting"} // .view operator to see the content of the channel
                        // .flatten()   // collapses array into series of values (character strings here)
                        //.view { greeting -> "after flatten: $greeting"} // .view operator to see the content of the channel}


    // emit a greeting
    sayHello(greeting_ch)
    //sayHello(params.greeting)
}
