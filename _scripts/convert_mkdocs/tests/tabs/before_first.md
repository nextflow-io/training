Fill in the module, then check it against the "After" tab.

=== "Before"

    ```nextflow title="modules/index.nf" showLineNumbers
    process INDEX {
    }
    ```

===+ "After"

    ```nextflow {2} title="modules/index.nf" showLineNumbers
    process INDEX {
        container 'samtools'
    }
    ```
