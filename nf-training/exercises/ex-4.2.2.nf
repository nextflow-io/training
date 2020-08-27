Channel.fromPath('data/ggal/*.fq')
        .map { file -> [ file.name, file ] }
        .println { name, file -> "> file: $name" }
        