nextflow.enable.types = true

include { Greeting } from '../types'

process TIMESTAMP_GREETING {
    tag "adding timestamp to greeting"

    input:
    greeting: Greeting

    output:
    record(name: greeting.name, file: file("timestamped_${greeting.file.baseName}.txt"))

    script:
    """
    echo "[\$(date '+%Y-%m-%d %H:%M:%S')] \$(cat ${greeting.file})" > timestamped_${greeting.file.baseName}.txt
    """
}
