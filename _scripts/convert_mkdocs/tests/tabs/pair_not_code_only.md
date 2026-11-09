An After/Before pair with prose in a tab cannot be a walkthrough:

=== "After"

    The fixed process:

    ```nextflow
    process A {}
    ```

=== "Before"

    ```nextflow
    process B {}
    ```

Nor can one where a tab holds two blocks:

=== "Before"

    ```nextflow title="main.nf"
    workflow {}
    ```

    ```nextflow title="nextflow.config"
    docker.enabled = false
    ```

=== "After"

    ```nextflow title="main.nf"
    workflow { A() }
    ```

    ```nextflow title="nextflow.config"
    docker.enabled = true
    ```

Three code tabs with After/Before labels are not a pair:

=== "Before"

    ```nextflow
    a
    ```

=== "After"

    ```nextflow
    b
    ```

=== "Original"

    ```nextflow
    c
    ```
