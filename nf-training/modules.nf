process SPLITLETTERS {
    input:
    val x

    output:
    file 'chunk_*'

    """
    printf '$x' | split -b 6 - chunk_
    """
}

process CONVERTTOUPPER {
    input:
    file y

    output:
    stdout

    """
    cat $y | tr '[a-z]' '[A-Z]' 
    """
}
