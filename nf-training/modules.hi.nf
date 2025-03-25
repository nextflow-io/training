process SPLITLETTERS {
    input:
    val x

    output:
    path 'chunk_*'

    script:
    """
    printf '$x' | split -b 6 - chunk_
    """
}

process CONVERTTOUPPER {
    input:
    path y

    output:
    stdout

    script:
    """
    cat $y | tr '[a-z]' '[A-Z]'
    """
}

process COUNTLETTERS {
    input:
    val z

    output:
    path 'count.txt'

    script:
    """
    printf '$z' | wc -m > count.txt
    """
}
