nextflow.preview.dsl=2

process foo {
  output:
    tuple path('file1.txt'), path('file2.txt'), emit: my_txt_files
  script:
    """
    touch file1.txt
    touch file2.txt
    """
}

process bar {
  output:
    path 'file1.csv', emit: csv_files
    path 'file2.bam', emit: bam_files
  script:
    """
    touch file1.csv
    touch file2.bam
    """
}

workflow flow1 {
    foo()
    foo.out.my_txt_files.view { "txt file: $it.name"}
}

workflow flow2 {
    bar()
    bar.out.csv_files.view { "cvs file: $it.name" }
    bar.out.bam_files.view { "cvs file: $it.name" }
}