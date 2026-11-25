/*
 * Convert a greeting to uppercase
 */
process SHOUT {
    tag "${meta.id}"

    input:
    tuple val(meta), path(greeting_file)

    output:
    tuple val(meta), path("${meta.id}-shouted.txt")

    script:
    """
    tr '[:lower:]' '[:upper:]' < ${greeting_file} > ${meta.id}-shouted.txt
    """
}
