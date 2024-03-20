/*
 * Pipeline parameters
 */
params.output_file = 'output.txt'

/*
 * Use echo to print 'Hello World!' to standard out
 */
process sayHello {

    output: 
        path params.output_file
    
    """
    echo 'Hello World!' > $params.output_file
    """
}

workflow {

    // emit a greeting
    sayHello()
}