nextflow.enable.dsl=2

process bar {
  echo true
  input:
  val x
  val y
  script:
   """
   echo $x and $y
   """
}

input1 = Channel.value(1)
input2 = Channel.from('a','b','c')


workflow{

  bar (input1, input2)

}
