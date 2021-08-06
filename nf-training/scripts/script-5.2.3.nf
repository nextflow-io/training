nextflow.enable.dsl=2

Channel
    .from( 'a', 'b', 'c' )
    .set{ abc }

abc
    .view{ "Foo emits: $it" }

abc
    .view{ "Bar emits: $it" }
