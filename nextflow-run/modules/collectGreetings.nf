/*
 * Collect uppercase greetings into a single output file
 */
process collectGreetings {

    publishDir 'results', mode: 'copy'

    input:
        path input_files

    output:
        path "COLLECTED-output.txt"

    script:
    """
    cat ${input_files} > 'COLLECTED-output.txt'
    """
}
