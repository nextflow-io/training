# Exercise 6

Move into the exercise 6 directory and execute the `hello-gatk.nf` script.

```bash
cd /workspace/gitpod/troubleshoot/exercise6
nextflow run hello-gatk.nf
```

!!! warning "Error message"

    ```
    ERROR ~ No such file or directory: /workspace/gitpod/troubleshoot/exercise6/data/samplesheet.csv

    -- Check '.nextflow.log' file for details
    ```

??? Solution

    The error message suggests that the path to the `samplesheet.csv` is incorrect.

    Searching the file directory, indeed, `/workspace/gitpod/troubleshoot/exercise6/data/samplesheet.csv` does not exist.

    The complete path should be `/workspace/gitpod/troubleshoot/data/samplesheet.csv`.

    By examining how the path to `samplesheet.csv` is derived it can be seen that there is no obvious path error.

    ```console title="hello-gatk.nf" linenums="5"
    // Execution environment setup
    params.projectDir = "/workspace/gitpod/troubleshoot"
    $projectDir = params.projectDir

    // Primary input
    params.reads_bam = "${projectDir}/data/samplesheet.csv"
    ```

    However, it can be noted that there is a `$` before `projectDir` on line 7 which would suggest it is a variable.

    `$projectDir` was once an implicit variable, and its usage here is causing the wrong path.

    By removing the `$` the proper path should be resolved.

    ```console title="hello-gatk.nf" linenums="7"
    baseDir = params.baseDir
    ```
