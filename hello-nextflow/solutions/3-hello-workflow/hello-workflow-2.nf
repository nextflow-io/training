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
 * Use Python to translate the greeting into another language
 */
process translateGreeting {

    input:
    path input_file
    val language

    output:
    path "${language}-${input_file}"

    script:
    """
    #!/usr/bin/env python3

    translations = {
        'french': 'Bonjour',
        'spanish': 'Hola',
        'italian': 'Ciao',
        'dutch': 'Hallo',
        'arabic': 'Marhaban',
        'japanese': 'Konnichiwa',
        'english': 'Hello',
    }

    with open('${input_file}') as greeting_file:
        greeting = greeting_file.read()

    with open('${language}-${input_file}', 'w') as translated_file:
        translated_file.write(greeting.replace('Hello', translations['${language}']))
    """
}

/*
 * Pipeline parameters
 */
params {
    input: Path = 'data/names.csv'
    language: String = 'spanish'
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
    translateGreeting(sayHello.out, params.language)

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
