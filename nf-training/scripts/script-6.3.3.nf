process splitLetters {

    output:
    file 'chunk_*' into letters

    '''
    printf 'Hola' | split -b 1 - chunk_
    '''
}

letters
    .flatMap()
    .view { "File: ${it.name} => ${it.text}" }