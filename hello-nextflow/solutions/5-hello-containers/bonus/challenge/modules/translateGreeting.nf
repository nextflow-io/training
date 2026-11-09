/*
 * Use Python to translate the greeting into another language
 */
process translateGreeting {

    input:
    tuple path(input_file), val(character)
    val language

    output:
    tuple path("${language}-${input_file}"), val(character)

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
