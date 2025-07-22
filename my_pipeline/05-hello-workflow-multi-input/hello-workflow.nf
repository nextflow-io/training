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

process convertToUppercase {
    
    publishDir 'results', mode: 'copy'

    input:
    path input_file

    output:
        path "UPPER-${input_file}-output.txt"

    script:
    """
    cat $input_file | tr '[a-z]' '[A-Z]' > "UPPER-${input_file}-output.txt"
    """
}

process collectGreetings {
    
    publishDir 'results', mode: 'copy'

    input:
    path input_file
    val batch_name 

    output:
        path "COLLECTED-${batch_name}-output.txt"

    script:
    """
    cat $input_file > "COLLECTED-${batch_name}-output.txt"
    """
}

/* * Passing default parameters to the workflow */
params.greeting = 'greetings.csv'
params.batch_name = 'test-batch'

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

    // Convert the output of sayHello to uppercase
    convertToUppercase(sayHello.out)

    // Collect all greetings into a single file
    collectGreetings(sayHello.out.collect(), params.batch_name)
    convertToUppercase.out.view{ greetings -> "Before collectGreetings: $greetings" }
    convertToUppercase.out.collect()
        .view{ greetings -> "After collectGreetings: $greetings" }


}