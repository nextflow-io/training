nextflow.preview.dsl=2

/*
 * The operator splitCsv allow to parse csv files
 */


workflow flow1 {
    Channel
      .fromPath("data/meta/regions.tsv", checkIfExists:true)
      // use `sep` option to parse TAB separated files
      .splitCsv(sep:'\t')
      // row is a list object 
      .view()
}



