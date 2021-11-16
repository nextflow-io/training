process convertToUpper {

    input:
    file y from letters.flatten()

    output:
    stdout into result

    """
    rev $y
    """
}