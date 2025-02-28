The directory has the file content like:

```console title="Directory contents"
nf-test
├── greetings.csv
├── main.nf
└── nextflow.config
```

For a detailed description of the files, see the [warmup from Hello Nextflow](../hello_nextflow/00_orientation.md).
The workflow we'll be testing is part of the workflow built in [Hello Workflow](../hello_nextflow/03_hello_workflow.md), and is composed of two processes: `sayHello` and `convertToUpper`:

```bash title="Workflow code"
/*
 * Pipeline parameters
 */
params.input_file = "greetings.csv"

/*
 * Use echo to print 'Hello World!' to standard out
 */
process sayHello {

    publishDir 'results', mode: 'copy'

    input:
        val greeting

    output:
        path "${greeting}-output.txt"

    script:
    """
    echo '$greeting' > '$greeting-output.txt'
    """
}

/*
 * Use a text replace utility to convert the greeting to uppercase
 */
process convertToUpper {

    publishDir 'results', mode: 'copy'

    input:
        path input_file

    output:
        path "UPPER-${input_file}"

    script:
    """
    cat '$input_file' | tr '[a-z]' '[A-Z]' > UPPER-${input_file}
    """
}

workflow {

    // create a channel for inputs from a CSV file
    greeting_ch = Channel.fromPath(params.input_file).splitCsv().flatten()

    // emit a greeting
    sayHello(greeting_ch)

    // convert the greeting to uppercase
    convertToUpper(sayHello.out)
}
```

We're going to assume an understanding of this workflow, but if you're not sure, you can refer back to [Hello Workflow](../hello_nextflow/03_hello_workflow.md).

### 0.1. Run the workflow

Let's run the workflow to make sure it's working as expected.

```bash
nextflow run main.nf
```

```console title="Result of running the workflow"
 N E X T F L O W   ~  version 24.10.2

Launching `main.nf` [soggy_linnaeus] DSL2 - revision: bbf79d5c31

executor >  local (6)
[f7/c3be66] sayHello (3)       | 3 of 3 ✔
[cd/e15303] convertToUpper (3) | 3 of 3 ✔
```
