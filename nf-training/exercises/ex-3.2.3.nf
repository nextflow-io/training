params.path = 'data/ggal/*.fq'
Channel.fromPath(params.path)
        .println { file -> "$file.name \tdirectory: $file.parent" }
