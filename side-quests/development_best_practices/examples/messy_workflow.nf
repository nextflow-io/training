#!/usr/bin/env nextflow

// Messy workflow with poor formatting and style
params.input='data/sample_data.csv'
params.output='results'

ch1=Channel.fromPath(params.input).splitCsv(header: true).map{row->tuple(row.sample_id,row.file_path)}

process proc1{
publishDir "${params.output}/proc1",mode:'copy'
input:
tuple val(x),val(y)
output:
tuple val(x),path("${x}_out.txt")
script:
"""
echo "processing: ${x}"
echo "file: ${y}">>${x}_out.txt
"""
}

process proc2 {
  publishDir "${params.output}/proc2", mode: 'copy'

  input:
    tuple val(sample), path(file)

  output:
    path "${sample}_final.txt"

  script:
  """
  cat ${file} > ${sample}_final.txt
  echo "final processing" >> ${sample}_final.txt
  """
}

workflow{
ch2=proc1(ch1)
ch3=proc2(ch2)
ch3.view{file->"Output: ${file}"}
}
