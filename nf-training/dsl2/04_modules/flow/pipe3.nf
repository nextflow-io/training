nextflow.enable.dsl=2

/*
 * A self-contained pipeline script can be used as
 * module reusing the workflow logic as a composable component
 */

include{ foo; bar } from '../some/module.nf'

workflow flow3 {
  // the `take` block can list as many input channel as needed
  take:
    fa_file
    bam_file

  // once used `take` we need to use `main` annotation to identify
  // where the workflow logic starts
  main:
    foo(fa_file)
    bar(bam_file)

  // the `emit` block allows the definition of the channel that
  // need to be exposed as the workflow output
  emit:
    foo.out
  // alternatively the assignment can be used to give a name to the emitted channel
    bar_data = bar.out
}

