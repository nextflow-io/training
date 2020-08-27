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
 * piping with matching cardinality
 */
workflow flow1 {
    foo | bar
}

/*
 * piping partial inputs
 *
 * currying ie. partial matching is not working
 */
workflow flow2 {
    foo | baz('/some/file.csv')
}

