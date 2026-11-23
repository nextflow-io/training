#!/usr/bin/env nextflow

include { COWPY } from './modules/cowpy.nf'

workflow  {
    main:
    ch_datasheet = channel.fromPath("./data/datasheet.csv")
        .splitCsv(header: true)
        .map { row ->
            [row.id, row.character, row.recording]
        }

    COWPY(ch_datasheet)

    publish:
    cowpy_art = COWPY.out
}

output {
    cowpy_art {
    }
}
