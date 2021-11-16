process randomNum {

    output:
    file 'result.txt' into numbers

    '''
    echo $RANDOM > result.txt
    '''
}

numbers.view { "Received: " + it.text }