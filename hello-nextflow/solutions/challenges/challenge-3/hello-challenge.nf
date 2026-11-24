#!/usr/bin/env nextflow

/*
 * Use echo to print a greeting to a file
 */
process sayHello {

    tag "${name}"

    input:
    tuple val(name), val(language)

    output:
    tuple path("${name}.txt"), val(language)

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
    tuple path(input_file), val(language)

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
    path "${batch_name}-${language}-greetings.txt", emit: outfile
    path "${batch_name}-${language}-report.txt", emit: report

    script:
    count_greetings = input_files.size()
    """
    cat ${input_files} > '${batch_name}-${language}-greetings.txt'
    echo 'Number of greetings processed: ${count_greetings}' > '${batch_name}-${language}-report.txt'
    """
}

/*
 * Pipeline parameters
 */
params {
    input: Path = 'data/names.csv'
    batch: String = 'batch'
}

workflow {

    main:
    // create a channel of [name, language] pairs from a CSV file
    person_ch = channel.fromPath(params.input)
                        .splitCsv()
                        .map { row -> [row[0], row[1]] }
    // emit a greeting
    greetings_ch = sayHello(person_ch)
    // translate each greeting into its own language
    translated_ch = translateGreeting(greetings_ch)
    // collect all the greetings into one file
    collected_out = collectGreetings(translated_ch.collect(), params.batch, 'native')

    publish:
    greetings = greetings_ch.map { greeting_file, _language -> greeting_file }
    translated = translated_ch
    collected = collected_out.outfile
    batch_report = collected_out.report
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
    batch_report {
        path "${params.batch}"
        mode 'copy'
    }
}
