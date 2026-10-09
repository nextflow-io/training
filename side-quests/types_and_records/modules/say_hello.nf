/*
 * Write a personalised greeting to a file
 */
process SAY_HELLO {
    tag "${meta.id}"

    input:
    val meta

    output:
    tuple val(meta), path("${meta.id}.txt")

    script:
    """
    echo '${meta.greeting}, ${meta.name}!' > ${meta.id}.txt
    """
}
