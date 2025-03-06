/*
 * Use a text manipulation tool to reverse the text in a file
 */
process REVERSE_TEXT {
    publishDir 'results', mode: 'copy'

    tag "reversing ${input_file}"

    input:
        path input_file

    output:
        path "REVERSED-${input_file}"

    script:
    """
    cat '${input_file}' | rev > 'REVERSED-${input_file}'
    """
}
