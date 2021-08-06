nextflow.enable.dsl=2

/*
 * Workflow can produce outputs as a process.
 * See script './flow/pipe3' 
 */
 
include{ flow3 } from './flow/pipe3'

process baz {
  input:
    path x
    path y
  script:
    """
    echo your_command $x $y
    """
}


workflow flow1 {
  x = file('data/ggal/transcript.fa')
  y = file('data/ggal/sample.bam')

  // invoke the include workflow 
  flow3(x,y)

  // access result channel by name
  flow3.out.bar_data.view { "bar: $it.name" }
}



workflow flow2 {
  x = file('data/ggal/transcript.fa')
  y = file('data/ggal/sample.bam')

  // workflow execution
  flow3(x,y)

  // workflow output can be composoed as/with any other process
  baz(flow3.out)
}
