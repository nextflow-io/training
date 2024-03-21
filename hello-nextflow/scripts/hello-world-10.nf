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
        path "${greeting}-${params.output_file}"
    
    """
    echo '$greeting' > '$greeting-$params.output_file'
    """
}

/*
 * Use a text replace utility to convert the greeting to uppercase
 */
process convertToUpper {
    input:
        path input_file

    output:
        path "UPPER-${input_file}"

    """
    cat $input_file | tr '[a-z]' '[A-Z]' > 'UPPER-${input_file}'
    """
}

workflow {

    // create a channel for inputs from a file
    greeting_ch = Channel.fromPath(params.input_file).splitText() { it.trim() }

    // emit a greeting
    sayHello(greeting_ch)

    // convert the greeting to uppercase
    convertToUpper(sayHello.out)
}
