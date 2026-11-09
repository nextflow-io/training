/*
 * Use echo to print a greeting to a file
 */
process sayHello {

    input:
    tuple val(name), val(character)

    output:
    tuple path("${name}-output.txt"), val(character)

    script:
    """
    echo 'Hello, ${name}!' > '${name}-output.txt'
    """
}
