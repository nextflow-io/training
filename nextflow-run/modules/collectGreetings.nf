/*
 * Collect uppercase greetings into a single output file
 */
process collectGreetings {

    input:
    path input_files, arity: '1..*'
    val batch_name

    output:
    path "COLLECTED-${batch_name}-output.txt", emit: outfile
    path "${batch_name}-report.txt", emit: report

    script:
    count_greetings = input_files.size()
    """
    cat ${input_files} > 'COLLECTED-${batch_name}-output.txt'
    echo 'Number of greetings processed: ${count_greetings}' > '${batch_name}-report.txt'
    """
}
