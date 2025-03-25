canal1 = Channel.of(1, 2, 3)
canal2 = Channel.of(1) // usa só o primeiro elemento do que vai somar
//canal2 = Channel.value(1) // usa todos os elementos do canal maior

process SUM {
    input:
    val x
    val y

    output:
    stdout

    script:
    """
    echo \$(($x+$y))
    """
}

workflow {
    //SUM(canal1, canal2).view() // consome só 1 vez
    SUM(canal1, canal2.first()).view() // consome todas as vezes
}
