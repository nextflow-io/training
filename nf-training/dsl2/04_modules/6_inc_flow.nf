nextflow.enable.dsl=2

/*
 * A workflow component can be included and used as any other module
 */
include{ flow1 } from './flow/pipe1'
include{ foo } from './other/module'

workflow {
  foo()
  flow1()
}
