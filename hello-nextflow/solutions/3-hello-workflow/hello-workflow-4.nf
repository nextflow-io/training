#!/usr/bin/env nextflow

/*
 * Use echo to print 'Hello World!' to a file
 */
process sayHello {

    input:
    val greeting

    output:
    path "${greeting}-output.txt"

    script:
    """
    echo '${greeting}' > '${greeting}-output.txt'
    """
}

/*
 * Use a text replacement tool to convert the greeting to uppercase
 */
process convertToUpper {

    input:
    path input_file

    output:
    path "UPPER-${input_file}"

    script:
    """
    cat '${input_file}' | tr '[a-z]' '[A-Z]' > 'UPPER-${input_file}'
    """
}

/*
 * Collect uppercase greetings into a single output file
 */
process collectGreetings {

    input:
    path input_files
    val batch_name

    output:
    path "COLLECTED-${batch_name}-output.txt", emit: outfile
    val count_greetings, emit: count

    script:
    count_greetings = input_files.size()
    """
    cat ${input_files} > 'COLLECTED-${batch_name}-output.txt'
    """
}

/*
 * Pipeline parameters
 */
params {
    input: Path = 'greetings.csv'
    batch: String = 'test-batch'
}

workflow {

    main:
    // create a channel for inputs from a CSV file
    greeting_ch = channel.fromPath(params.input)
                        .splitCsv()
                        .map { line -> line[0] }
    // emit a greeting
    sayHello(greeting_ch)
    // convert the greeting to uppercase
    convertToUpper(sayHello.out)
    // collect all the greetings into one file
    collectGreetings(convertToUpper.out.collect(), params.batch)

    // emit a message about the size of the batch
    collectGreetings.out.count.view { num_greetings -> "There were ${num_greetings} greetings in this batch" }

    publish:
    first_output = sayHello.out
    upper_output = convertToUpper.out
    collected_output = collectGreetings.out
}

output {
    first_output {
        path 'hello_workflow'
        mode 'copy'
    }
    upper_output {
        path 'hello_workflow'
        mode 'copy'
    }
    collected_output {
        path 'hello_workflow'
        mode 'copy'
    }
}
