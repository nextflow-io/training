nextflow.enable.dsl=2

Channel.from(1,2,3,4)
    .map { it -> it * it }
    .view()
