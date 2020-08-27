
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
  output:
    tuple path('file1.txt'), path('file2.txt')
  script:
    """
    touch file1.txt
    touch file2.txt
    """
}

/*
 * the process output is implicitly returned by process invocation
 * any channel operator can be applied to it
 */
 
workflow flow1 {
  foo().view { "foo out: $it.name" }
}


/* 
 * the process output can be assigned to a variable 
 * and it can used as any other channel
 * 
 * note: DSL2 allow to read the same channel more than one time 
 * without the need to fork it explicitly 
 */

workflow flow2 {
  some_ch = foo()
  some_ch.view { "foo out 1: $it.name" }
  some_ch.view { "foo out 2: $it.name" }
}


/*
 * once a process has been invoked, the `out` implicit variable can be
 * used using the process name as scope to access the process output
 */

workflow flow3 {
  foo()

  foo.out.view { "foo out 1: $it.name" }
  foo.out.view { "foo out 2: $it.name" }
}
