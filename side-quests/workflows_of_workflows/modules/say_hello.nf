/*
 * Use echo to print 'Hello World!' to a file
 */
process SAY_HELLO {
    tag "greeting ${name}"

    input:
    val name

    output:
    path "${name}-output.txt"

    script:
    """
    echo 'Hello, ${name}!' > "${name}-output.txt"
    """
}
