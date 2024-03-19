params.output_file = 'hello-output.txt'
greeting_ch = Channel.of('Hello world!')

workflow {
    sayHello(greeting_ch)
}

process sayHello {
    input:
        val greeting  

    output: 
        path params.output_file
    
    """
    echo '$greeting' > $params.output_file
    """
}