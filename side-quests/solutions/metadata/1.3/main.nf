#!/usr/bin/env nextflow

include { COWPY } from './modules/cowpy.nf'

workflow  {
    main:
    ch_datasheet = channel.fromPath("./data/datasheet.csv")
        .splitCsv(header: true)
        .multiMap { row ->
            file: row.recording
            character: row.character
        }

    COWPY(ch_datasheet.file, ch_datasheet.character)

    publish:
    cowpy_art = COWPY.out
}

output {
    cowpy_art {
    }
}
