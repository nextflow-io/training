#!/usr/bin/env nextflow

nextflow.preview.types = true

/*
 * Pipeline parameters
 */
params {
    greeting: Path = 'greetings.csv'
    batch: String = 'test-batch'
}

/*
 * Use echo to print 'Hello World!' to a file
 */
process sayHello {

    publishDir 'results', mode: 'copy'

    input:
    greeting: String

    output:
    file "${greeting}-output.txt"

    script:
    """
    echo '$greeting' > '$greeting-output.txt'
    """
}

/*
 * Use a text replacement tool to convert the greeting to uppercase
 */
process convertToUpper {

    publishDir 'results', mode: 'copy'

    input:
    input_file: Path

    output:
    file "UPPER-${input_file}"

    script:
    """
    cat '$input_file' | tr '[a-z]' '[A-Z]' > 'UPPER-${input_file}'
    """
}

/*
 * Collect uppercase greetings into a single output file
 */
process collectGreetings {

    publishDir 'results', mode: 'copy'

    input:
    input_files: List<Path>
    batch_name: String

    output:
    outfile: Path = file("COLLECTED-${batch_name}-output.txt")
    count: Integer = input_files.size()

    script:
    """
    cat ${input_files} > 'COLLECTED-${batch_name}-output.txt'
    """
}

workflow {

    // create a channel for inputs from a CSV file
    greeting_ch = channel.fromPath(params.greeting)
                        .splitCsv()
                        .map { line -> line[0] }

    // emit a greeting
    sayHello(greeting_ch)

    // convert the greeting to uppercase
    convertToUpper(sayHello.out)

    // collect all the greetings into one file
    collectGreetings(convertToUpper.out.collect(), params.batch)

    // emit a message about the size of the batch
    collectGreetings.out.count.view { "There were $it greetings in this batch" }
}
