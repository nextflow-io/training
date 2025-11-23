#!/usr/bin/env nextflow

params.greeting = 'Hello world!'
greeting_ch = channel.of(params.greeting)

process SPLITLETTERS {
    input:
    val x

    output:
    path 'chunk_*'

    script:
    """
    #!/usr/bin/env python
    x="$x"
    for i, word in enumerate(x.split()):
        with open(f"chunk_{i}", "w") as f:
            f.write(word)
    """
}

process CONVERTTOUPPER {
    input:
    path y

    output:
    stdout

    script:
    """
    #!/usr/bin/env python
    with open("$y") as f:
        print(f.read().upper(), end="")
    """
}

workflow {
    letters_ch = SPLITLETTERS(greeting_ch)
    results_ch = CONVERTTOUPPER(letters_ch.flatten())
    results_ch.view{ it }
}
