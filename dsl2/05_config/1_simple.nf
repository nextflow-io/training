nextflow.preview.dsl=2

include{ foo } from './mod.nf'

/*
 * For included processes the usual config rules applied.
 *
 * It tries to match the config for process named `foo` (and label)
 * then fallback on the generic process config
 */

workflow {
  foo()
}