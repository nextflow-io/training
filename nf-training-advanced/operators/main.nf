
workflow {
    channel.of( 1, 2, 3, 4, 5 )
        .map { num -> num * num }
        .view()
}
