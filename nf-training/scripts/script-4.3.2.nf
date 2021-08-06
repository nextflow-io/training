nextflow.enable.dsl=2

process randomNum {

    output:
    file 'result.txt'

    '''
    echo $RANDOM > result.txt
    '''
}


workflow{

  receiver_ch = randomNum()
  receiver_ch.view { "Received: " + it.text }

}

