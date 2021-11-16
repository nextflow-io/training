sequences = Channel.fromPath('data/prots/*.tfa')
methods = ['regular', 'expresso', 'psicoffee']

process alignSequences {
  input:
  path seq from sequences
  each mode from methods

  """
  echo t_coffee -in $seq -mode $mode
  """
}
