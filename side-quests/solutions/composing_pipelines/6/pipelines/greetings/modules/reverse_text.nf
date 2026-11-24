nextflow.enable.types = true

include { Greeting } from '../types'

/*
 * Use a text manipulation tool to reverse the text in a file
 */
process REVERSE_TEXT {
    tag "reversing ${greeting.file.name}"

    input:
    greeting: Greeting

    output:
    record(name: greeting.name, file: file("REVERSED-${greeting.file.name}"))

    script:
    """
    cat ${greeting.file} | rev > REVERSED-${greeting.file.name}
    """
}
