process SEPARARLETRAS {
    output:
    path 'pedaco_*'

    script:
    """
    printf 'Hola' | split -b 1 - pedaco_
    """
}

workflow {
    letters = SEPARARLETRAS()
    letters
        //.flatMap()
        .view { "Arquivo: ${it.name} => ${it.text}" }
}
