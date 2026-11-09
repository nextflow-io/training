nextflow.enable.types = true

include { Greeting } from '../types'

/*
 * Use a text replacement tool to convert the greeting to uppercase
 */
process SAY_HELLO_UPPER {
    tag "converting ${greeting.file.name}"

    input:
    greeting: Greeting

    output:
    record(name: greeting.name, file: file("UPPER-${greeting.file.name}"))

    script:
    """
    cat ${greeting.file} | tr '[a-z]' '[A-Z]' > UPPER-${greeting.file.name}
    """
}
