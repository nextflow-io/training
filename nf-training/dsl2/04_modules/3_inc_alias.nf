nextflow.enable.dsl=2

/*
 * INCLUDE ALIAS
 *
 * The same process can be imported (and used) more than one time
 * using a name alias.
 */

include{ foo; foo as bar } from './other/module.nf'

workflow {
    foo()
    bar()
}
