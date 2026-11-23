// Generate ASCII art with cowpy (https://github.com/jeffbuttars/cowpy)
process COWPY {

    container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
    conda 'conda-forge::cowpy==1.1.5'

    input:
        tuple val(meta), path(input_file)

    output:
        tuple val(meta), path("${prefix}.txt")                                    , emit: cowpy_output
        tuple val("${task.process}"), val('cowpy'), val("1.1.5"), topic: versions , emit: versions_cowpy

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    cat $input_file | cowpy $args > ${prefix}.txt
    """
}
