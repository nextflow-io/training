nextflow.enable.dsl=2

methods = ['prot','dna', 'rna']

process foo {
  input:
  val x

  output:
  val x

  """
  echo $x > file
  """
}

workflow{

  receiver_ch = foo(Channel.from(methods))
  receiver_ch.view { "Received: $it" }

}

