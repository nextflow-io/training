nextflow.enable.dsl=2

/*
 * MODIFY SPECIFIC MODULE PARAMS
 *
 * It is possible to specific module params with the options:
 *
 * - params
 * - addParams
 */

params.alpha = 'Hello'

// in this case the module will only receive the params `omega`
// from the invoking environment

include{ foo } from './other/module.nf' params(omega:1)


// the `addParams` option extends the invoking environment
// adding/overriding the specified parameters. In this case
// the resulting params are [alpha: 'Hello', omega: 9]

include{ foo as bar } from './other/module.nf' addParams(omega:9)

workflow {
  foo()
  bar()
}
