nextflow.preview.dsl=2

process foo {
  input:
    val data
  output:
    val result
  exec:
    result = "$data world"
}

process bar {
  input:
    val data
  output:
    val result
  exec:
    result = data.toUpperCase()
}

process baz {
  input:
    val x
    val y
  output:
    val z
  exec:
    z = x + y
}

/*
 * Two, or more processes, can take the same inputs
 * using the `&` operator. The result is a multi-channel output \
 */
workflow flow1 {
   Channel.of('Hello') | map { it.reverse() } | (foo & bar) | baz | view
}

/*
 * The `mix` operator can be used to "multiplex" multiple outs to one channel
 */
workflow flow2 {
   Channel.of('Hello') | map { it.reverse() } | (foo & bar) | mix | view
}

/*
 * alternative formatting
 */

workflow pretty {
   Channel.of('Hello') \
     | map { it.reverse() } \
     | (foo & bar) \
     | mix \
     | view
}