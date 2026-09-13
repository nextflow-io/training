// Generate ASCII art with cowpy
process COWPY {

    container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'

    input:
    path recording
    val character

    output:
    path "cowpy-${recording}"

    script:
    """
    cat ${recording} | cowpy -c ${character} > cowpy-${recording}
    """
}
