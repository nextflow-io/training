
nextflow.preview.dsl=2

process foo {
  output:
    tuple path('file1.txt'), path('file2.txt')
  script:
    """
    touch file1.txt
    touch file2.txt
    """
}

process bar {
  output:
    path('file1.txt')
    path('file2.txt')
  script:
    """
    touch file1.txt
    touch file2.txt
    """
}

/*
 * out is a single channel, therefore operators can be applied
 */
workflow flow1 {
  foo()
  foo.out.view { "Out files: $it.name" }
  
  // ^^ wait, why `it.name` works since in the above line considering `it` is a tuple?
  // because groovy automatically applies the spread operator the invoked
  // property (ie `name`) is not defined in the tuple (ie list) object,
  // therefore it's the same as: `it*.name`
  // which is the same as: `it.collect { file -> file.name }`
  // see https://groovy-lang.org/operators.html#_spread_operator
}

/*
 * this fails because `bar` out is a multi channel
 */
workflow flow2 {
  bar()
  bar.out.view()
}

/*
 * components of a multi channel out can be accessed using array like syntax
 */
workflow flow3 {
  bar()
  bar.out[0].view { "First: $it.name" }
  bar.out[1].view { "Second: $it.name" }
}

/*
 * multi assignment syntax can be used with multi channel outputs
 */
workflow flow4 {
 def (ch1, ch2) = bar()
 ch1.view { "First: $it.name" }
 ch2.view { "Second: $it.name" }
}

/*
 * as before using the `out` property
 */
workflow flow5 {
 bar()
 def (ch1, ch2) = bar.out
 ch1.view { "First: $it.name" }
 ch2.view { "Second: $it.name" }
}
