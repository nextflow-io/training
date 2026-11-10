/*
 * Collect translated greetings into a single output file
 */
process collectGreetings {

    input:
    path input_files
    val batch_name

    output:
    path "${batch_name}-greetings.txt", emit: outfile
    path "${batch_name}-report.txt", emit: report

    script:
    count_greetings = input_files.size()
    """
    cat ${input_files} > '${batch_name}-greetings.txt'
    echo 'Number of greetings processed: ${count_greetings}' > '${batch_name}-report.txt'
    """
}
