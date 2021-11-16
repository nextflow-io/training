methods = ['prot','dna', 'rna']

process foo {
  input:
  val x from methods

  output:
  val x into receiver

  """
  echo $x > file
  """
}

receiver.view { "Received: $it" }