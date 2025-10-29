#!/usr/bin/env nextflow

/*
 * Use echo to print 'Hello World!' to a file
 */
process sayHello {

    publishDir 'results', mode: 'copy'

    input:
        val greeting

    output:
        path "${greeting}-output.txt"

    script:
    """
    echo '$greeting' > '$greeting-output.txt'
    """
}

process toUpper {
    publishDir 'results', mode: 'copy'

    input:
    path input_file

    output:
    path "UPPER-${input_file}.txt"

    script:
    """
    cat $input_file | tr '[a-z]' '[A-Z]' > UPPER-${input_file}.txt
    """
}

process collectAllGret {
    publishDir 'results', mode: 'copy'

    input:
    path input_file
    val batch_name

    output:
    path "Collected-${batch_name}-output.txt", emit: outfile
    val greeting_count, emit: count

    script:
    greeting_count = input_file.size()
    """
    cat $input_file > Collected-${batch_name}-output.txt
    """
}

/*
 * Pipeline parameters
 */
params.greeting = 'greetings.csv'
params.batch_name = 'test-batch'

workflow {

    // create a channel for inputs from a CSV file
    greeting_ch = Channel.fromPath(params.greeting)
                        .splitCsv()
                        .map { line -> line[0] }

    // emit a greeting
    sayHello(greeting_ch)

    toUpper(sayHello.out)

    collectAllGret(toUpper.out.collect(), params.batch_name)
    collectAllGret.out.count.view { num -> "There were $num Greetings!"}
}
