??? exercise "Fix the workflow"

    Run it first:

    ```bash
    nextflow run buggy.nf
    ```

    ??? failure "Command output"

        ```console
        ERROR ~ Script compilation failed
        ```

=== "Docker"

    ```bash
    nextflow run main.nf -profile docker
    ```

    ??? success "Command output"

        ```console
        Done
        ```
