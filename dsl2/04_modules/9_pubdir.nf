nextflow.enable.dsl=2

/*
 * the usual `publishDir` works well, the plan is to
 * replace with a more decoupled mechanism in the future
 */
 
params.outdir = 'my_results'

include{ foo;bar } from './some/pub.nf'

workflow {
  x = file('data/prot/BB11001.fa')
  foo(x)
  bar(foo.out)
}
