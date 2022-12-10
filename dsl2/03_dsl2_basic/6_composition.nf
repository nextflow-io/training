nextflow.preview.dsl=2

process foo {
  output:
    path 'foo.txt'
  script:
    """
    echo your_command > foo.txt
    """
}

 process bar {
  input:
    path x
   output:
    path 'bar.txt'
  script:
    """
    echo another_command $x > bar.txt
    """
}

workflow test1 {
    foo()
    bar(foo.out)
}


workflow test2 {
    bar(foo())
}