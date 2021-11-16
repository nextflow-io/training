process bar {
  script:
  '''
  echo $PATH | tr : '\\n'
  '''
}