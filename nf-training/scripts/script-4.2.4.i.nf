nextflow.enable.dsl=2

process foo {
  echo true
  input:
  val x
  val y
  script:
   """
   echo $x and $y
   """
}

input1 = Channel.from(1,2,3)
input2 = Channel.from('a','b','c')

workflow{

  foo (input1,input2)

}
