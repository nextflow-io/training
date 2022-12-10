nextflow.preview.dsl=2

include{ foo as bar; foo as baz } from './mod.nf'

/*
 * Process imported with an alias, it tries
 * 1) to match first the alias name eg `bar`,
 * 2) then fallback with the process original eg `foo`
 * 3) finally it used the process generic settings
 */

workflow {
  bar()
  baz()
}