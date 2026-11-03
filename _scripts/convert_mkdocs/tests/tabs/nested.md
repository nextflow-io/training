There are two ways to fix this:

=== "Option 1: Extract first element"

    === "After"

        ```nextflow {2} title="bad_channel_shape.nf" showLineNumbers startLineNumber=16
        workflow {
            PROCESS_FILES(input_ch.map { sample_id, file -> sample_id })
        }
        ```

    === "Before"

        ```nextflow {2} title="bad_channel_shape.nf" showLineNumbers startLineNumber=16
        workflow {
            PROCESS_FILES(input_ch)
        }
        ```

=== "Option 2: Accept tuple in process"

    === "After"

        ```nextflow {3} title="bad_channel_shape.nf" showLineNumbers
        process PROCESS_FILES {
            input:
                tuple val(sample_name), val(file_name)
        }
        ```

    === "Before"

        ```nextflow {3} title="bad_channel_shape.nf" showLineNumbers
        process PROCESS_FILES {
            input:
                val sample_name
        }
        ```

Both options work.
