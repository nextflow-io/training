nextflow.preview.dsl=2

/*
 * Nextflow parameters can also be specified using a yaml or json file
 *
 * This spares using low level APIs, the resulting values can be accessed
 * through the usual `params` implicit variable
 *
 */

params.samples = [
    ['foo', '/data/pair1.fq', '/data/pair2.fq']
]

process foo {
  echo true
  input:
  tuple val(sampleId), path(p1), path(p2)

  """
  echo your_command $sampleId $p1 $p2
  """
}


workflow {
  Channel.fromList(params.samples) | foo
}

