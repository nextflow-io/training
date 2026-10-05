nextflow.enable.types = true

include { Person } from '../types.nf'

/*
 * Write a personalised greeting to a file
 */
process SAY_HELLO {
    tag "${person.id}"

    input:
    person: Person

    output:
    person + record(greeting_file: file("${person.id}.txt"))

    script:
    """
    echo '${person.greeting}, ${person.name}!' > ${person.id}.txt
    """
}
