nextflow.preview.types = true

/*
 * Collect uppercase greetings into a single output file
 */
process collectGreetings {

    publishDir 'results', mode: 'copy'

    input:
    input_files: List<Path>
    batch_name: String

    output:
    outfile: Path = file("COLLECTED-${batch_name}-output.txt")
    count: Integer = input_files.size()

    script:
    """
    cat ${input_files} > 'COLLECTED-${batch_name}-output.txt'
    """
}
