nextflow.preview.dsl=2

/*
 * The operator splitCsv allow to parse csv files
 */


workflow flow1 {
    Channel
      .fromPath("data/meta/patients_1.csv")
      .splitCsv()
      // row is a list object 
      .view { row -> "${row[0]} +++ ${row[2]}" }
}

/*
 * Using the header to create dictionary for each row
 */
workflow flow2 {
    Channel
      .fromPath("data/meta/patients_1.csv")
      // the option `header` specify to interpret first line as column names
      .splitCsv(header:true)
      // now each column can be accessed by name
      .view { row -> "${row.patient_id} +++ ${row.s3_dir}" }
}


/*
 * process many csv at once 
 */
workflow flow3 {
    Channel
      .fromPath("data/meta/patients_*.csv") // <-- just use a pattern
      .splitCsv(header:true)
      .view { row -> "${row.patient_id} +++ ${row.s3_dir}" }
}


