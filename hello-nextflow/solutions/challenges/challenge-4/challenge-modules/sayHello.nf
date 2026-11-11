/*
 * Use echo to print a greeting to a file
 */
process sayHello {

    tag "${name}"

    input:
    tuple val(name), val(character)

    output:
    tuple path("${name}.txt"), val(character)

    script:
    """
    echo 'Hello, ${name}!' > '${name}.txt'
    """
}
