nextflow.enable.dsl=2

/*
 * A self-contained pipeline script can be used as
 * module reusing the workflow logic as a composable component
 */

include{ foo; bar } from '../some/module.nf'


workflow flow1 {
  data = file('data/prot/BB11001.fa')
  foo(data)
  bar(foo.out)
}

workflow {
  flow1()
}
