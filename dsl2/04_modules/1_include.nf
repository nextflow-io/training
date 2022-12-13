nextflow.enable.dsl=2

/*
 * INCLUDE
 *
 * A process can be included from a separate script
 * 
 * note: module path must start with `./` or `../` or `/` 
 */

include { foo } from './some/module.nf'

workflow {
    data = file('data/prot/BB11002.fa')
    foo(data)
}

