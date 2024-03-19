params.output_file = 'hello-output.txt'

workflow {
    sayHello()
}

process sayHello {

    output: 
        path params.output_file
    
    """
    echo 'Hello World!' > $params.output_file
    """
}