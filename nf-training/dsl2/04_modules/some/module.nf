process foo {
  input:
    path x
  output:
    path 'file.txt'

  """
   rev $x > file.txt
  """
}

process bar {
  input:
    path x
  output:
    path 'lower.txt'

  """
   cat $x | tr '[A-Z]' '[a-z]' > lower.txt
  """
}
