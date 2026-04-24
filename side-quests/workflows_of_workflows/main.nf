workflow {
    main:
    names = channel.of('Alice', 'Bob', 'Charlie')

    publish:
    greetings = channel.empty()
}

output {
    greetings {
        path 'greetings'
    }
}
