Channel
    .fromFilePairs( params.reads, checkIfExists: true )
    .set { read_pairs_ch }