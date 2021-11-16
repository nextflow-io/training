Channel
    .from( 'hello', 'world' )
    .map { it -> it.reverse() }
    .view()