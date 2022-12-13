nextflow.enable.dsl=2

/*
 * MULTIPLE INCLUSIONS
 *
 * multiple processes can be imported using
 * curly brackets notation
 */
 
include{ foo; bar } from './some/module.nf'

workflow {
    data = file('data/prot/BB11001.fa')
    foo(data)
    bar(foo.out)
    bar.out.view { it.text }
}
