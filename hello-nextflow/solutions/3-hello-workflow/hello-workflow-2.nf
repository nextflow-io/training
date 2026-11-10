#!/usr/bin/env nextflow

/*
 * Use echo to print a greeting to a file
 */
process sayHello {

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
    import sys

    translations = {
        'es': 'Hola',        # Spanish
        'fr': 'Bonjour',     # French
        'it': 'Ciao',        # Italian
        'nl': 'Hallo',       # Dutch
        'ar': 'Marhaban',    # Arabic
        'ja': 'Konnichiwa',  # Japanese
        'en': 'Hello',       # English
    }

    if '${language}' not in translations:
        sys.exit('Unsupported language code: ${language}')

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
    language: String = 'es'
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
    greetings = sayHello.out
    translated = translateGreeting.out
}

output {
    greetings {
        path '.'
    }
    translated {
        path '.'
    }
}
