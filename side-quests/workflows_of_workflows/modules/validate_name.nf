process VALIDATE_NAME {
    tag "validating ${name}"

    input:
    val name

    output:
    val name

    script:
    """
    if [[ "${name}" =~ [^a-zA-Z] ]]; then
        echo "Error: Name '${name}' contains invalid characters" >&2
        exit 1
    fi
    variable=$name
    if [[ "\${#variable}" -lt 2 ]]; then
        echo "Error: Name '${name}' is too short" >&2
        exit 1
    fi
    """
}
