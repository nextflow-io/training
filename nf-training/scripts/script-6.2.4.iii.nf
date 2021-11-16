process bar {
  echo true
  input:
  val x from Channel.value(1)
  val y from Channel.from('a','b','c')
  script:
   """
   echo $x and $y
   """
}