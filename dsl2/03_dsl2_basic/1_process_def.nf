/*
 * Since DSL represents a major enhancement of Nextflow syntax with some breaking changes,
 * it needs to be enabled explicitly starting your script with the following directive
 *
 * This allows the use of the new features without breaking any existing scripts. Alternatively it's possible to
 * enable DLS2 using the command line option `-dsl2`.
 */
 
nextflow.enable.dsl=2


process foo {
  output:
    path 'foo.txt'
  script:
    """
    echo your_command > foo.txt
    """
}

workflow {
  foo()
}

// note: the same process can be invoked only once in the same context

// you can use `NXF_IGNORE_WARN_DSL2` to suppress warning message
