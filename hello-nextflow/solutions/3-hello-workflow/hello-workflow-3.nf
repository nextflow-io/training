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
 * Use Python to translate the greeting into another language
 */
process translateGreeting {

    tag "${language}-${input_file}"

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
 * Collect translated greetings into a single output file
 */
process collectGreetings {

    tag "${batch_name}"

    input:
    path input_files
    val batch_name
    val language

    output:
    path "${batch_name}-${language}-greetings.txt"

    script:
    """
    cat ${input_files} > '${batch_name}-${language}-greetings.txt'
    """
}

/*
 * Pipeline parameters
 */
params {
    input: Path = 'data/names.csv'
    language: String = 'es'
    batch: String = 'batch'
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
    translated_ch = translateGreeting(greetings_ch, params.language)
    // collect all the greetings into one file
    collected_ch = collectGreetings(translated_ch.collect(), params.batch, params.language)

    publish:
    greetings = greetings_ch
    translated = translated_ch
    collected = collected_ch
}

output {
    greetings {
        path "${params.batch}/intermediates"
        mode 'copy'
    }
    translated {
        path "${params.batch}/intermediates"
        mode 'copy'
    }
    collected {
        path "${params.batch}/intermediates"
        mode 'copy'
    }
}
