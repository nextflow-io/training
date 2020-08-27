/*
 * custom function to split the transcript and reads channel
 */
def split_trans_and_reads(combined) {
  combined.multiMap { row ->
       trans: row[0];
       reads: tuple(row[1],row[2])
       }
}

/*
 * custom function that capture the pipeline input logic
 */
def rna_inputs(String pathTranscript, String pathReads) {

  def trans = Channel .fromPath( pathTranscript )
  def reads = Channel .fromFilePairs( pathReads, checkExists: true )

  trans \
    | combine(reads) \
    | split_trans_and_reads

}
