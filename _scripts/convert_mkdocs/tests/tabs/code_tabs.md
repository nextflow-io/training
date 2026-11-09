=== "Process definition"

    ```nextflow title="modules/cowpy.nf" showLineNumbers startLineNumber=6
    input:
    tuple val(meta), val(character), path(recording)
    ```

=== "Workflow call"

    ```nextflow title="main.nf"
    COWPY(ch_languages.map { meta, file -> [meta, meta.character, file] })
    ```

A lone code tab becomes a one-item `<Tabs>`, since code tabs need two blocks:

=== "Single tab"

    ```bash
    nextflow run main.nf
    ```
