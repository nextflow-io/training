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
 * Collect translated greetings into a single output file
 */
process collectGreetings {

    input:
    path input_files
    val batch_name

    output:
    path "COLLECTED-${batch_name}-output.txt", emit: outfile
    path "${batch_name}-report.txt", emit: report

    script:
    count_greetings = input_files.size()
    """
    cat ${input_files} | sort > 'COLLECTED-${batch_name}-output.txt'
    echo 'Number of greetings processed: ${count_greetings}' > '${batch_name}-report.txt'
    """
}

/*
 * Pipeline parameters
 */
params {
    input: Path = 'data/names.csv'
    language: String = 'french'
    batch: String = 'batch'
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
    // collect all the greetings into one file
    collectGreetings(translateGreeting.out.collect(), params.batch)

    publish:
    first_output = sayHello.out
    translated = translateGreeting.out
    collected = collectGreetings.out.outfile
    batch_report = collectGreetings.out.report
}

output {
    first_output {
        path "${params.batch}/intermediates"
    }
    translated {
        path "${params.batch}/intermediates"
    }
    collected {
        path "${params.batch}/intermediates"
    }
    batch_report {
        path "${params.batch}"
    }
}
