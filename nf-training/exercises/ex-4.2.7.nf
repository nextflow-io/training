Channel.fromPath('data/ggal/*.fq')
        .map{ [ it.name.tokenize('_')[0], it] }
        .groupTuple()
        .println()
