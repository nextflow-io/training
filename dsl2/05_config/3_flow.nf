nextflow.preview.dsl=2

include{ foo; foo as bar } from './mod.nf'

/*
 * When a process is nested into named workflow the config is resolved
 * 1) against the extended name eg. `flow1:foo`
 * 2) the process name
 * 3) the generic config
 */

workflow flow1 {
  foo()
  bar()
}

/*
 * Same applies for nested workflow
 */

workflow flow2 {
  flow1()
}

workflow {
  flow1()
  flow2()
}