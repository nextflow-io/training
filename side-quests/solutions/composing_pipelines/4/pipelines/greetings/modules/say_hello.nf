nextflow.enable.types = true

/*
 * Use echo to print 'Hello World!' to a file
 */
process SAY_HELLO {
    tag "greeting ${name}"

    input:
    name: String

    output:
    record(name: name, file: file("${name}-output.txt"))

    script:
    """
    echo 'Hello, ${name}!' > "${name}-output.txt"
    """
}
