nextflow.preview.dsl=2

/*
 * the usual operators can be used to parse many files
 */

include{ parseJsonFile } from './modules/parsers.nf'

process foo {
  input:
  tuple val(meta), path(data_file)

  """
  echo your_command $meta.region_id $data_file
  """
}

/*
 * show how:
 * - create a library of parse helpers
 * - parse multiple json files
 * - pipe json metadata to NF process
 */
workflow {
    Channel.fromPath('data/meta/regions*.json') \
      | flatMap { parseJsonFile(it) } \
      | map { entry -> tuple(entry,"/some/data/${entry.patient_id}.txt") } \
      | foo
}