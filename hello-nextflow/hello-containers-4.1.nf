#!/usr/bin/env nextflow

/*
 * Pipeline parameters
 */
params.input_file = "containers/data/pioneers.csv"
// params.character can be any of 'beavis', 'cheese', 'cow', 'daemon', 'dragon', 'fox', 'ghostbusters', 'kitty',
// 'meow', 'miki', 'milk', 'octopus', 'pig', 'stegosaurus', 'stimpy', 'trex', 'turkey', 'turtle', 'tux'
params.character = "cow"

/*
 * Use echo to print 'Hello World!' to standard out
 */
process sayHello {

    publishDir 'containers/results', mode: 'copy'

    input:
        val greeting  

    output: 
        path "output-*.txt"
    
    script:
        // Replace the spaces in the greeting with hyphens for the output filename
        def safe_greeting = greeting.tokenize(' ').join('-')
        """
        echo '$greeting' > 'output-${safe_greeting}.txt'
        """
}

/*
 * Use a text replace utility to convert the greeting to uppercase
 */
process cowSay {

    publishDir 'containers/results', mode: 'copy'
    container 'community.wave.seqera.io/library/pip_cowsay:131d6a1b707a8e65'
    
    input:
        path input_file

    output:
        path "cowsay-*"

    """
    cowsay -c "$params.character" -t "\$(cat $input_file)" > cowsay-${input_file}
    """
}

process quote {

    publishDir 'containers/results', mode: 'copy'
    container 'community.wave.seqera.io/library/pip_quote:25b3982790125217'

    input:
        val author

    output:
        path "quote-*.txt"

    script:
        // Replace the spaces in the author with hyphens for the output filename
        def safe_author = author.tokenize(' ').join('-')
        """
        # The stdout and stderr of the quote command are swapped so we redirect them to the correct files
        # Pending fix: https://github.com/maxhumber/quote/issues/12
        # quote "$author" > quote-${safe_author}.txt
        quote "$author" > /dev/null 2> quote-${safe_author}.txt || true
        """
}

workflow {

    // create a channel for inputs from a CSV file
    input_ch = Channel.fromPath(params.input_file)
        .splitCsv()
        .flatten()

    getQuote(input_ch)

    // cowSay the text
    cowSay(getQuote.out)
}
