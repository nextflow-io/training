Channel.fromPath('data/ggal/*.fq')
        .map { file -> [ file.name, file ] }
        .view { name, file -> "> $name : $file" }