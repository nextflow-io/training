workflow {
    sayHello()
}

process sayHello {

    output: 
        stdout
    
    """
    echo 'Hello World!'
    """
}