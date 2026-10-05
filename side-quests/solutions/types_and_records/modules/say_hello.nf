nextflow.enable.types = true

/*
 * Write a personalised greeting to a file
 */
process SAY_HELLO {
    tag "${id}"

    input:
    record(id: String, name: String, greeting: String)

    output:
    record(id: id, greeting_file: file("${id}.txt"))

    script:
    """
    echo '${greeting}, ${name}!' > ${id}.txt
    """
}
