/*
 * Collect translated greetings into a single output file
 */
process collectGreetings {

    tag "${batch_name}"

    input:
    path input_files
    val batch_name
    val language

    output:
    path "${batch_name}-${language}-greetings.txt", emit: outfile
    path "${batch_name}-${language}-report.txt", emit: report

    script:
    count_greetings = input_files.size()
    """
    cat ${input_files} > '${batch_name}-${language}-greetings.txt'
    echo 'Number of greetings processed: ${count_greetings}' > '${batch_name}-${language}-report.txt'
    """
}
