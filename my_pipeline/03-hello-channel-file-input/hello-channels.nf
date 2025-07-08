#!/usr/bin/env nextflow

/*
 * Use echo to print 'Hello World!' to a file
 */
process sayHello {
    
    publishDir 'results', mode: 'copy'

    input:
    val greeting

    output:
        path "$greeting-output.txt"

    script:
    """
    echo '$greeting' > '$greeting-output.txt'
    """
}

/* * Passing default parameters to the workflow */
params.greeting = 'greetings.csv'

workflow {
    // create a channel from the CSV file
    greeting_ch = channel.fromPath(params.greeting)
                         .view{ csv -> "Before splitCsv: $csv" }
                         .splitCsv()
                         .view{ csv -> "After splitCsv: $csv" }
                         .map{ item -> item[0] }
                         .view{ csv -> "After map: $csv" }

    // emit a greeting
    // nf created the channels automatically
    sayHello(greeting_ch)
}