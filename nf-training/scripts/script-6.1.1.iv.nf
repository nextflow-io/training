params.data = 'le monde'

process baz {
  shell:
  '''
  X='Bonjour'
  echo $X !{params.data}
  '''
}