/*
 * Use echo to print 'Hello World!' to standard out
 */
process sayHello {

    output: 
        path 'output.txt'
    
    """
    echo 'Hello World!' > output.txt
    """
}

workflow {

    // emit a greeting
    sayHello()
}