ch1 = Channel.from(1,2,3)
ch2 = Channel.from('a','b','c')

process foo {
  echo true
  input:
  val x from ch1
  val y from ch2
  script:
   """
   echo $x and $y
   """
}