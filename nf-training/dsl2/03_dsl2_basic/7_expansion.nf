nextflow.preview.dsl=2

process foo {
  output:
    path '*.txt'
    path '*.bam'

  script:
    """
    touch file.txt
    touch file.bam
    """
}

 process bar {
  input:
    path x
    path y

  script:
    """
    echo this_command $x $y
    """
}

process baz {
  input:
    path x
    path y
    path z

  script:
    """
    echo that_command $x $y $z
    """
}

/*
 * composing with matching cardinality
 */
workflow test1 {
    foo()
    bar(foo.out)
}

/*
 * expanding multiple output
 */
workflow test2 {
    foo()
    baz(foo.out, '/some/file.csv')
}

workflow test3 {
    foo()
    baz('/some/file.csv', foo.out)
}


/*
 * accessing single elements
 */
process another {
  output:
    path '*.txt'
    path '*.csv'
    path '*.bam'
  script:
    """
    touch file.txt
    touch file.csv
    touch file.bam
    """
}

workflow test4 {
    another()
    bar(another.out[0], another.out[2])
}