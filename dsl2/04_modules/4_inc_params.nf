nextflow.enable.dsl=2

/*
 * MODULE PARAMS
 *
 * Included modules inherit the params from the invoking context
 *
 * NOTE: parameters are resolved when the inclusion is declared
 * therefore it's suggested to keep module inclusions *after* the
 * params definition block
 *
 */

params.alpha = 'Hello'

include{ foo } from './other/module.nf'

workflow {
    foo()
}
