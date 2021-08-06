nextflow.enable.dsl=2

num = Channel.from( 1, 2, 3 )

process basicExample {
  input:
  val x 

  """
  echo process job $x
  """
}

workflow {
   
    myrun = basicExample(num)

}
