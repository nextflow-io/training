/*
 * Collect all greetings into a single file
 */
process COLLECT_GREETINGS {

    input:
    path input_files
    val batch_name

    output:
    path "COLLECTED-${batch_name}.txt"

    script:
    """
    cat ${input_files} > COLLECTED-${batch_name}.txt
    """
}
