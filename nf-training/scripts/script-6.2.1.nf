num = Channel.from( 1, 2, 3 )

process basicExample {
  input:
  val x from num

  """
  echo process job $x
  """
}