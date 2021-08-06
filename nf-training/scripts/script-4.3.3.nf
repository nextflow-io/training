nextflow.enable.dsl=2

process splitLetters {

    output:
    file 'chunk_*'

    '''
    printf 'Hola' | split -b 1 - chunk_
    '''
}



workflow{
    letters = splitLetters()	
    letters
        .flatMap()
        .view { "File: ${it.name} => ${it.text}" }

}


