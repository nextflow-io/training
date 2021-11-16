/*
 * define the `index` process that creates a binary index
 * given the transcriptome file
 */
process index {

    input:
    path transcriptome from params.transcriptome_file

    output:
    path 'index' into index_ch

    script:
    """
    salmon index --threads $task.cpus -t $transcriptome -i index
    """
}