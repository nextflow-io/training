/*
 * Workflow script can contain one ore more named `workflow` definition.
 * Use the `-entry` command line option to select the flow to run.
 */

nextflow.preview.dsl=2

process foo {
  input:
    path x
  output:
    path 'foo.txt'
  script:
    """
    echo your_command $x > foo.txt
    """
}


workflow flow1 {
  data = file("data/ggal/liver_1.fq")
  foo(data)
}

workflow flow2 {
  data = file("data/ggal/lung_1.fq")
  foo(data)
}

workflow {
  data = file("data/ggal/gut_1.fq")
  foo(data)
}
