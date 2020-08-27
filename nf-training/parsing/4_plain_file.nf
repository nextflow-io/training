nextflow.preview.dsl=2

/*
 * Splitter operators `splitText`, `splitCsv`, `splitFast` can also
 * be applied to regular file path
 */


workflow flow1 {

  def f = file('data/meta/random.txt')
  def lines = f.splitText()
  def count=0
  for( String row : lines ) {
    log.info "${count++} ${row.toUpperCase()}"
  }

}


workflow flow2 {

  def f = file('data/meta/patients_1.csv')
  def lines = f.splitCsv()
  for( List row : lines ) {
    log.info "${row[0]} -- ${row[2]}"
  }

}