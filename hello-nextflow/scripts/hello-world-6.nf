/*
 * Pipeline parameters
 */
params.output_file = 'output.txt'

/*
 * Use echo to print 'Hello World!' to standard out
 */
process sayHello {
    input:
        val greeting  

    output: 
        path params.output_file
    
    """
    echo '$greeting' > $params.output_file
    """
}

workflow {

    // create a channel for inputs
    greeting_ch = Channel.of('Hello world!')

    // emit a greeting
    sayHello(greeting_ch)
}