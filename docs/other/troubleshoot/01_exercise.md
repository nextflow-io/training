# Exercise 1

Move into the exercise 1 directory and execute the `hello-gatk.nf` script.

```bash
cd /workspaces/training/troubleshoot/exercise1
nextflow run hello-gatk.nf
```

!!! warning "Error message"

    Nextflow stops the workflow execution and reports the failing task when a process execution exits with a non-zero exit status:

    ```
    ERROR ~ Error executing process > 'SAMTOOLS_INDEX (2)'

    Caused by:
      Process `SAMTOOLS_INDEX (2)` terminated with an error exit status (127)

    Command executed:

      samtools index 'reads_father.bam'

    Command exit status:
      127

    Command output:
      (empty)

    Command error:
      .command.sh: line 2: samtools: command not found

    Work dir:
      /workspaces/training/troubleshoot/exercise1/work/46/e01f274e2ef1735164061d62c51169

    Tip: when you have fixed the problem you can continue the execution adding the option `-resume` to the run command line

    -- Check '.nextflow.log' file for details
    ```

    **Caused by**: A description of the error cause

    **Command executed**: The command executed

    **Command exit status**: The command exit status

    **Command output**: The command standard output, when available

    **Command error**: The command standard error

    **Work dir**: The command work directory

??? Solution

    The command error `.command.sh: line 2: samtools: command not found` suggests that `samtools` was not available.

    Every process in the `hello-gatk.nf` pipeline contains a container directive with a relevant container image. However, in this case, docker was not enabled and the `samtools` tooling was not available locally.

    Enabling docker in the `nextflow.config` file will activate the docker containers for each process and resolve this error.

    ```console title="nextflow.config" linenums="1"
    docker.enabled = true
    ```
