#!/usr/bin/env nextflow

ch = Channel.of("1","2","3")
ch.view()
