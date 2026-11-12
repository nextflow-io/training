#!/usr/bin/env nextflow

/*
 * Use echo to print a greeting to a file
 */
process sayHello {

    tag "${name}"

    input:
    val name

    output:
    path "${name}.txt"

    script:
    """
    echo 'Hello, ${name}!' > '${name}.txt'
    """
}

/*
 * Use Python to translate the greeting into Spanish
 */
process translateGreeting {

    tag "${input_file}"

    input:
    path input_file

    output:
    path "es-${input_file}"

    script:
    """
    #!/usr/bin/env python3

    with open('${input_file}') as greeting_file:
        greeting = greeting_file.read()

    with open('es-${input_file}', 'w') as translated_file:
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
    greetings_ch = sayHello(name_ch)
    // translate the greeting
    translated_ch = translateGreeting(greetings_ch)

    publish:
    greetings = greetings_ch
    translated = translated_ch
}

output {
    greetings {
        path 'hello'
        mode 'copy'
    }
    translated {
        path 'hello'
        mode 'copy'
    }
}
