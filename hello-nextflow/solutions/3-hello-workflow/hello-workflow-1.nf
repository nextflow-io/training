#!/usr/bin/env nextflow

/*
 * Use echo to print a greeting to a file
 */
process sayHello {

    input:
    val name

    output:
    path "${name}-output.txt"

    script:
    """
    echo 'Hello, ${name}!' > '${name}-output.txt'
    """
}

/*
 * Use Python to translate the greeting into Spanish
 */
process translateGreeting {

    input:
    path input_file

    output:
    path "spanish-${input_file}"

    script:
    """
    #!/usr/bin/env python3

    with open('${input_file}') as greeting_file:
        greeting = greeting_file.read()

    with open('spanish-${input_file}', 'w') as translated_file:
        translated_file.write(greeting.replace('Hello', 'Hola'))
    """
}

/*
 * Pipeline parameters
 */
params {
    input: Path = 'data/names.csv'
}

workflow {

    main:
    // create a channel for inputs from a CSV file
    name_ch = channel.fromPath(params.input)
                        .splitCsv()
                        .map { row -> row[0] }
    // emit a greeting
    sayHello(name_ch)
    // translate the greeting
    translateGreeting(sayHello.out)

    publish:
    first_output = sayHello.out
    translated = translateGreeting.out
}

output {
    first_output {
        path '.'
    }
    translated {
        path '.'
    }
}
