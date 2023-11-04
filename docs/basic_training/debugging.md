---
title: Troubleshooting
description: Error handling and troubleshooting
---

# Error handling and troubleshooting

## Execution error debugging

When a process execution exits with a non-zero exit status, Nextflow stops the workflow execution and reports the failing task:

!!! info ""

    Click the :material-plus-circle: icons in the code for explanations.

```bash
ERROR ~ Error executing process > 'INDEX'

Caused by: # (1)!
  Process `INDEX` terminated with an error exit status (127)

Command executed: # (2)!

  salmon index --threads 1 -t transcriptome.fa -i index

Command exit status: # (3)!
  127

Command output: # (4)!
  (empty)

Command error: # (5)!
  .command.sh: line 2: salmon: command not found

Work dir: # (6)!
  /Users/pditommaso/work/0b/b59f362980defd7376ee0a75b41f62
```

1. A description of the error cause
2. The command executed
3. The command exit status
4. The command standard output, when available
5. The command standard error
6. The command work directory

Carefully review all error data as it can provide information valuable for debugging.

If this is not enough, `cd` into the task work directory. It contains all the files to replicate the issue in an isolated manner.

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

Also verify the existence of the `.exitcode` and `.command.begin` files, which if absent, suggest the task was never executed by the subsystem (e.g. the batch scheduler). If the `.command.begin` file exists, the task was launched but was likely killed abruptly.

You can replicate the failing execution using the command `bash .command.run` to verify the cause of the error.

## Explore exactly what is run

If you want to see the full commands executed when running `bash .command.run`,
you can edit the `.command.run` file and add the following line somewhere near
the top, but under the `#!/bin/bash` part:

```bash
set -x
```

This will make sure that all executed commands are printed out as they are
executed.

Note that you might want to do this for `.command.sh` also, since this is the
actual command generated from the task definition in the pipeline code.

You can also explore the `nxf_launch()` bash function inside `.command.run` in
a text editor, to make sure that it looks as expected. This can be useful
especially in a new environment, to make sure that containers are initialized
correctly.

Finally, if you want to explore the temporary folder where commands are being
executed, you need to comment out the line that deletes this folder after the
command is finished or interrupted. This line you can find in the
`nxf_onexit()` function in `.command.run`.

There, you can comment out (Add a `#` character in the beginning of the line)
for the line that looks like so:

```bash
rm -rf $NXF_SCRATCH || true
```

... so that it instead becomes:

```bash
#rm -rf $NXF_SCRATCH || true
```

And then you can run `bash .command.run` again. Make note of a folder named
something like `/tmp/nxf.XXXXXXXXX` in the output from the command, as this is
the temporary execution folder. Now you can `cd` into this folder, since it is
not removed after the command is finished. There you can check for example that
the input files are correctly linked and accessible.

## Ignore errors

There are cases in which a process error may be expected and it should not stop the overall workflow execution.

To handle this use case, set the process `errorStrategy` to `ignore`:

```groovy linenums="1"
process FOO {
    errorStrategy 'ignore'

    script:
    """
    your_command --this --that
    """
}
```

If you want to ignore any error you can set the same directive in the config file as a default setting:

```groovy
process.errorStrategy = 'ignore'
```

## Automatic error fail-over

In rare cases, errors may be caused by transient conditions. In this situation, an effective strategy is re-executing the failing task.

```groovy linenums="1"
process FOO {
    errorStrategy 'retry'

    script:
    """
    your_command --this --that
    """
}
```

Using the `retry` error strategy the task is re-executed a second time if it returns a non-zero exit status before stopping the complete workflow execution.

The directive [maxRetries](https://www.nextflow.io/docs/latest/process.html#maxretries) can be used to set the number of attempts the task can be re-executed before declaring it failed with an error condition.

## Retry with backoff

There are cases in which the required execution resources may be temporarily unavailable (e.g. network congestion). In these cases simply re-executing the same task will likely result in an identical error. A retry with an exponential backoff delay can better recover these error conditions.

```groovy linenums="1"
process FOO {
    errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
    maxRetries 5

    script:
    '''
    your_command --here
    '''
}
```

## Dynamic resources allocation

It’s a very common scenario that different instances of the same process may have very different needs in terms of computing resources. In such situations requesting, for example, an amount of memory too low will cause some tasks to fail. Instead, using a higher limit that fits all the tasks in your execution could significantly decrease the execution priority of your job in a scheduling system.

To handle this use case, you can use a `retry` error strategy and increase the computing resources allocated by the task at each successive _attempt_.

```groovy linenums="1"
process FOO {
    cpus 4
    memory { 2.GB * task.attempt } // (1)!
    time { 1.hour * task.attempt } // (2)!
    errorStrategy { task.exitStatus == 140 ? 'retry' : 'terminate' } // (3)!
    maxRetries 3 // (4)!

    script:
    """
    your_command --cpus $task.cpus --mem $task.memory
    """
}
```

1. The memory is defined in a dynamic manner, the first attempt is 2 GB, the second 4 GB, and so on.
2. The wall execution time is set dynamically as well, the first execution attempt is set to 1 hour, the second 2 hours, and so on.
3. If the task returns an exit status equal to `140` it will set the error strategy to `retry` otherwise it will terminate the execution.
4. It will retry the process execution up to three times.
