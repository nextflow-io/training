Make the change:

=== "After"

    ```nextflow {3} title="hello.nf" showLineNumbers startLineNumber=4
    process SAY_HELLO {

        publishDir 'results'
    }
    ```

=== "Before"

    ```nextflow title="hello.nf" showLineNumbers startLineNumber=4
    process SAY_HELLO {
    }
    ```

Run the workflow.
