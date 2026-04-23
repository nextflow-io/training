#!/usr/bin/env nextflow

params.input = 'greetings.csv'

/**
 * Generate a random alphanumeric string
 */
def randomString(int length) {
    def chars = ('a'..'z') + ('A'..'Z') + ('0'..'9')
    def random = new Random()
    return (1..length).collect { chars[random.nextInt(chars.size())] }.join()
}

process SAY_HELLO {
    input:
        val greeting
    output:
        stdout
    script:
    """
    echo '$greeting'
    """
}

workflow {
    greeting_ch = channel.fromPath(params.input)
        .splitCsv(header: true)
        .map { row -> "${row.greeting}_${randomString(8)}" }
    SAY_HELLO(greeting_ch)
    SAY_HELLO.out.view { result -> "Output: ${result.trim()}" }
}
