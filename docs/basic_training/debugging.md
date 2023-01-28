---
title: Troubleshooting
description: Error handling and troubleshooting
---

# Error handling and troubleshooting

## Execution error debugging

When a process execution exits with a non-zero exit status, Nextflow stops the workflow execution and reports the failing task:

    ERROR ~ Error executing process > 'index'

    Caused by:
      Process `index` terminated with an error exit status (127)

    Command executed:

      salmon index --threads 1 -t transcriptome.fa -i index

    Command exit status:
      127

    Command output:
      (empty)

    Command error:
      .command.sh: line 2: salmon: command not found

    Work dir:
      /Users/pditommaso/work/0b/b59f362980defd7376ee0a75b41f62

-   A description of the error cause

-   The command executed

-   The command exit status

-   The command standard output, when available

-   The command standard error

-   The command work directory

Carefully review all error data as it can provide information valuable for debugging.

If this is not enough, `cd` into the task work directory. It contains all the files to replicate the issue in an isolated manner.

The task execution directory contains these files:

-   `.command.sh`: The command script.

-   `.command.run`: The command wrapped used to run the job.

-   `.command.out`: The complete job standard output.

-   `.command.err`: The complete job standard error.

-   `.command.log`: The wrapper execution output.

-   `.command.begin`: Sentinel file created as soon as the job is launched.

-   `.exitcode`: A file containing the task exit code.

-   Task input files (symlinks)

-   Task output files

Verify that the `.command.sh` file contains the expected command and all variables are correctly resolved.

Also verify the existence of the `.exitcode` and `.command.begin` files, which if absent, suggest the task was never executed by the subsystem (e.g. the batch scheduler). If the `.command.begin` file exists, the job was launched but was likely killed abruptly.

You can replicate the failing execution using the command `bash .command.run` to verify the cause of the error.

## Ignore errors

There are cases in which a process error may be expected and it should not stop the overall workflow execution.

To handle this use case, set the process `errorStrategy` to `ignore`:

    process foo {
      errorStrategy 'ignore'
      script:
      """
        your_command --this --that
      """
    }

If you want to ignore any error you can set the same directive in the config file as a default setting:

    process.errorStrategy = 'ignore'

## Automatic error fail-over

In rare cases, errors may be caused by transient conditions. In this situation, an effective strategy is re-executing the failing task.

    process foo {
      errorStrategy 'retry'
      script:
      """
        your_command --this --that
      """
    }

Using the `retry` error strategy the task is re-executed a second time if it returns a non-zero exit status before stopping the complete workflow execution.

The directive [maxRetries](https://www.nextflow.io/docs/latest/process.html#maxretries) can be used to set the number of attempts the task can be re-executed before declaring it failed with an error condition.

## Retry with backoff

There are cases in which the required execution resources may be temporarily unavailable (e.g. network congestion). In these cases simply re-executing the same task will likely result in an identical error. A retry with an exponential backoff delay can better recover these error conditions.

    process foo {
      errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
      maxRetries 5
      script:
      '''
      your_command --here
      '''
    }

## Dynamic resources allocation

It’s a very common scenario that different instances of the same process may have very different needs in terms of computing resources. In such situations requesting, for example, an amount of memory too low will cause some tasks to fail. Instead, using a higher limit that fits all the tasks in your execution could significantly decrease the execution priority of your jobs.

To handle this use case, you can use a `retry` error strategy and increase the computing resources allocated by the job at each successive _attempt_.

    process foo {
      cpus 4
      memory { 2.GB * task.attempt }
      time { 1.hour * task.attempt }
      errorStrategy { task.exitStatus == 140 ? 'retry' : 'terminate' }
      maxRetries 3

      script:
      """
        your_command --cpus $task.cpus --mem $task.memory
      """
    }

-   The memory is defined in a dynamic manner, the first attempt is 2 GB, the second 4 GB, and so on.

-   The wall execution time is set dynamically as well, the first execution attempt is set to 1 hour, the second 2 hours, and so on.

-   If the task returns an exit status equal to `140` it will set the error strategy to `retry` otherwise it will terminate the execution.

-   It will retry the process execution up to three times.
