process splitLetters {

    output:
    file 'chunk_*' into letters

    '''
    printf 'Hola' | split -b 1 - chunk_
    '''
}

letters
    .println { "File: ${it.name} => ${it.text}" }