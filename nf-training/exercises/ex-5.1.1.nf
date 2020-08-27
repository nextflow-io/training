params.data = 'World'

process foo {
  script:
  """
  echo Hello $params.data
  """
}

process bar {
  script:
  '''
  echo $PATH | tr : '\\n'
  '''
}

process baz {
  shell:
  '''
  X='Bonjour' 
  echo $X !{params.data}
  '''
}
