workflow {
    sayHello()
}

process sayHello {

    output: 
        path 'hello-output.txt'
    
    """
    echo 'Hello World!' > hello-output.txt
    """
}