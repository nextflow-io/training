nextflow.enable.dsl=2

/*
 * Workflow takes input as any other process
 */
 
include{ flow2 } from './flow/pipe2'

workflow {
  x = file('data/prot/BB11001.fa')
  flow2(x)
}
