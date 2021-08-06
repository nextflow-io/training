nextflow.enable.dsl=2

sequences = Channel.fromPath('data/prots/*.tfa')
methods = ['regular', 'expresso', 'psicoffee']

process alignSequences {
  echo true
  input:
  path seq
  each mode

  """
  echo t_coffee -in $seq -mode $mode
  """
}


workflow{
	
	alignSequences(sequences,methods)

}