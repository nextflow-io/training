
# Script 1

```bash
cd /workspace/gitpod/troubleshoot/exercise1
```

```bash
nextflow run hello-gatk.nf
```

## Error message

Nextflow stops the workflow execution and reports the failing task when a process execution exits with a non-zero exit status:

```
ERROR ~ Error executing process > 'SAMTOOLS_INDEX (2)'

Caused by: # (1)!
  Process `SAMTOOLS_INDEX (2)` terminated with an error exit status (127)

Command executed: # (2)!

  samtools index 'reads_father.bam'

Command exit status: # (3)!
  127

Command output: # (4)!
  (empty)

Command error: # (5)!
  .command.sh: line 2: samtools: command not found

Work dir: # (6)!
  /workspace/gitpod/troubleshoot/exercise1/work/46/e01f274e2ef1735164061d62c51169

Tip: when you have fixed the problem you can continue the execution adding the option `-resume` to the run command line

 -- Check '.nextflow.log' file for details
```

1. A description of the error cause
2. The command executed
3. The command exit status
4. The command standard output, when available
5. The command standard error
6. The command work directory

Carefully review all error data as it can provide information valuable for debugging.

If this is not enough, `cd` into the task work directory.

It contains all the files to replicate the issue in an isolated manner.

The task execution directory contains these files:

-   `.command.sh`: The command script.
-   `.command.run`: The command wrapped used to run the task.
-   `.command.out`: The complete task standard output.
-   `.command.err`: The complete task standard error.
-   `.command.log`: The wrapper execution output.
-   `.command.begin`: Sentinel file created as soon as the task is launched.
-   `.exitcode`: A file containing the task exit code.
-   Task input files (symlinks)
-   Task output files

Verify that the `.command.sh` file contains the expected command and all variables are correctly resolved.

??? Solution

    The command error `.command.sh: line 2: samtools: command not found` suggests that `samtools` was not available.

    Every process in the `hello-gatk.nf` pipeline contains a container directive with a relevant container. However, in this case, docker was not enabled and the `samtools` tooling was not available locally. Enabling docker in the `nextflow.config` file will activate the docker containers for each process and resolve this error.

    ```console title="nextflow.config"
    docker.enabled = true
    ```