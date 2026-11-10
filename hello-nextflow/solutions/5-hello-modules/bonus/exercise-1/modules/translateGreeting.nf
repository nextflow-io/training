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
