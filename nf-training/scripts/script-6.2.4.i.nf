process foo {
  echo true
  input:
  val x from Channel.from(1,2,3)
  val y from Channel.from('a','b','c')
  script:
   """
   echo $x and $y
   """
}