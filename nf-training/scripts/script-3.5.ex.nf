Channel
    .fromFilePairs( params.reads, checkIfExists: true )
    .into { read_pairs_ch; read_pairs2_ch }