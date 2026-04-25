// Generate ASCII art with cowpy
process COWPY {

    container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'

    input:
    tuple val(meta), path(input_file)

    output:
    tuple val(meta), path("${meta.lang}-${input_file}")

    script:
    """
    cat ${input_file} | cowpy -c ${meta.character} > ${meta.lang}-${input_file}
    """
}
