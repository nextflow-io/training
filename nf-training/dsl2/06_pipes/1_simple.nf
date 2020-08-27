nextflow.preview.dsl=2

process foo {
    input: val data
    output: val result
    exec:
    result = "$data world"
}

/*
 * Nextflow processes and operators can be composed using the `|` pipe operator
 */

workflow {
   Channel.of('Hello','Hola','Ciao') | foo | map { it.toUpperCase() } | view
}


/*
 * same using a concise formatting
 */

workflow pretty {
   Channel.of('Hello','Hola','Ciao') \
   | foo \
   | map { it.toUpperCase() } \
   | view
}