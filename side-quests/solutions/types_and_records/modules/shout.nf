nextflow.enable.types = true

/*
 * Convert a greeting to uppercase
 */
process SHOUT {
    tag "${id}"

    input:
    record(id: String, greeting_file: Path)

    output:
    record(id: id, shouted: file("${id}-shouted.txt"))

    script:
    """
    tr '[:lower:]' '[:upper:]' < ${greeting_file} > ${id}-shouted.txt
    """
}
