nextflow.enable.dsl=2

/*
 * A self-contained pipeline script can be used as
 * module reusing the workflow logic as a composable component
 */

include{ foo; bar } from '../some/module.nf'

workflow flow2 {
  take:
    data
  main:
    foo(data)
    bar(foo.out)
}

