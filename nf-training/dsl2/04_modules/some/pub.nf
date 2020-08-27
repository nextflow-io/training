params.outdir = 'results'

process foo {
  publishDir "$params.outdir/$task.process"

  input:
    path x
  output:
    path 'file.txt'

  """
   rev $x > file.txt
  """
}

process bar {
  publishDir "$params.outdir/$task.process"

  input:
    path x
  output:
    path 'lower.txt'

  """
   cat $x | tr '[A-Z]' '[a-z]' > lower.txt
  """
}
