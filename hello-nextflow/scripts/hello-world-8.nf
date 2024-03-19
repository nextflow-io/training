params.output_file = 'hello.txt'
greeting_ch = Channel.of('hi','why','no')

workflow {
    sayHello(greeting_ch)
    convertToUpper(sayHello.out)
}

process sayHello {
    input:
        val greeting  

    output: 
        path "${greeting}-${params.output_file}"
    
    """
    echo '$greeting' > '$greeting-$params.output_file'
    """
}

process convertToUpper {
    input:
        path input_file

    output:
        path "UPPER-${input_file}"

    """
    cat $input_file | tr '[a-z]' '[A-Z]' > 'UPPER-${input_file}'
    """
}

