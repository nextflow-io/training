nextflow.enable.types = true

/*
 * Collect all greetings into a single file
 */
process COLLECT_GREETINGS {

    input:
    input_files: Bag<Path>
    batch_name: String

    output:
    file("COLLECTED-${batch_name}.txt")

    script:
    """
    cat ${input_files.join(' ')} > COLLECTED-${batch_name}.txt
    """
}
