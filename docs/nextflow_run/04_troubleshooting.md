---
title: Troubleshooting
description: Diagnose and fix common issues when running Nextflow pipelines
---

# Part 4: Troubleshooting

When running pipelines, things don't always go smoothly.
This section covers common issues you'll encounter and how to diagnose and resolve them.

---

## 1. Understanding task failures

When a Nextflow task fails, understanding what went wrong is the first step to fixing it.

### 1.1. Reading the error output

When a task fails, Nextflow displays an error message that includes the task work directory:

```console
ERROR ~ Error executing process > 'SAYHELLO (1)'

Caused by:
  Process `SAYHELLO (1)` terminated with an error exit status (1)

Command executed:

  echo 'Hello, World!'  > output.txt

Command exit status:
  1

Command output:
  (empty)

Command error:
  .command.sh: line 2: unexpected EOF while looking for matching `''

Work dir:
  /workspace/gitpod/nf-training/work/a3/7be2fa1e...
```

The key pieces of information here are:

- **Exit status**: Non-zero means failure (1 is a general error, 127 means command not found, 137 means out of memory)
- **Command error**: The actual error message from the failed command
- **Work dir**: Where to find all the task files for debugging

### 1.2. Examining the work directory

Navigate to the task's work directory to examine what happened:

```bash
cd /workspace/gitpod/nf-training/work/a3/7be2fa1e...
ls -la
```

The most useful files for troubleshooting are:

| File           | Purpose                                              |
| -------------- | ---------------------------------------------------- |
| `.command.sh`  | The exact script that was executed                   |
| `.command.err` | Standard error output (error messages)               |
| `.command.out` | Standard output                                      |
| `.command.log` | Combined stdout and stderr                           |
| `.command.run` | The wrapper script Nextflow uses to execute the task |
| `.exitcode`    | The exit code of the command                         |

Start with `.command.err` - it usually contains the actual error message:

```bash
cat .command.err
```

If that's empty, check `.command.log`:

```bash
cat .command.log
```

### 1.3. Re-running a failed command manually

Sometimes you need to re-run the command manually to debug it.
From inside the task's work directory:

```bash
bash .command.sh
```

This lets you see the error in real-time and experiment with fixes.

---

## 2. Common error messages and solutions

### 2.1. Command not found

```console
Command error:
  .command.sh: line 3: cowsay: command not found
```

**Cause**: The required software isn't installed or isn't in the PATH.

**Solutions**:

1. Check if a container is configured for the process
2. Enable Docker/Singularity in your config: `docker.enabled = true`
3. If using Conda, ensure `conda.enabled = true` is set
4. For local execution, install the missing software

### 2.2. File not found

```console
Command error:
  cat: input.txt: No such file or directory
```

**Cause**: An expected input file wasn't staged into the work directory.

**Solutions**:

1. Check that input paths in your samplesheet or params are correct
2. Verify the file exists: `ls -la /path/to/input.txt`
3. Check if the file was produced by an upstream process that failed
4. Look for typos in file paths (case-sensitive on Linux)

### 2.3. Out of memory (exit code 137)

```console
Command exit status:
  137
```

**Cause**: The task was killed because it exceeded its memory allocation.

**Solutions**:

1. Increase memory allocation in your config:

   ```groovy title="nextflow.config"
   process {
       withName: 'MEMORY_INTENSIVE_PROCESS' {
           memory = '16.GB'
       }
   }
   ```

2. Check if your system has enough memory available
3. Use `-with-report` to see actual memory usage and tune allocations

### 2.4. Permission denied

```console
Command error:
  bash: ./script.sh: Permission denied
```

**Cause**: The script or file lacks execute permissions, or you don't have write access to a directory.

**Solutions**:

1. Check file permissions: `ls -la script.sh`
2. Ensure output directories are writable
3. If using containers, check that volume mounts have correct permissions

### 2.5. Container image not found

```console
ERROR ~ Error executing process > 'PROCESS_NAME'

Caused by:
  Failed to pull image 'docker://nonexistent/image:latest'
```

**Cause**: The specified container image doesn't exist or can't be accessed.

**Solutions**:

1. Verify the image name and tag are correct
2. Check if the image is in a private registry (may need authentication)
3. Try pulling the image manually: `docker pull image:tag`
4. Check your internet connection

---

## 3. Debugging techniques

### 3.1. Use verbose logging

When console output isn't enough, disable ANSI logging to see complete output:

```bash
nextflow run main.nf -ansi-log false
```

This prints all log messages line by line instead of updating in place, making it easier to capture and review.

### 3.2. Inspect resolved configuration

To see the complete configuration Nextflow is using (after merging all config files):

```bash
nextflow config -flat
```

Or to see it as a structured file:

```bash
nextflow config
```

This helps identify configuration issues like incorrect parameter values or missing settings.

### 3.3. Use the execution report

Generate a detailed execution report to understand resource usage:

```bash
nextflow run main.nf -with-report report.html
```

Open `report.html` in a browser to see:

- Which tasks took the longest
- Memory and CPU usage per task
- Timeline of execution
- Tasks that were cached vs. executed

### 3.4. Check the Nextflow log

The `.nextflow.log` file in your execution directory contains detailed debugging information:

```bash
cat .nextflow.log | grep -i error
```

For even more detail, enable trace logging:

```bash
NF_DEBUG=true nextflow run main.nf
```

---

## 4. Pipeline hangs or runs slowly

### 4.1. Check for stuck tasks

If a pipeline seems to hang, look at the console output.
Tasks showing `0 of N` for a long time may indicate:

- A process waiting for upstream data that will never arrive
- A container that's slow to pull
- A scheduler queue that's full

### 4.2. Identify bottlenecks

Run with `-with-report` and examine the timeline view.
Look for:

- Single long-running tasks blocking others
- Many tasks queued waiting for resources
- Excessive time in the "submit" state

### 4.3. Local executor limits

By default, the local executor uses all available CPUs.
If you're running out of memory, limit parallelism:

```groovy title="nextflow.config"
executor {
    name = 'local'
    cpus = 4
    memory = '8.GB'
}
```

---

## 5. Resume isn't working

### 5.1. Why tasks re-run unexpectedly

Tasks re-run instead of using cached results when:

1. **Input files changed**: Even a timestamp change can invalidate cache
2. **Process code changed**: Any modification to the process script
3. **Container changed**: Different image or tag
4. **Work directory was deleted**: Cache depends on work directory contents
5. **Different parameters**: Any change to parameters used by the process

### 5.2. Check what changed

Use `nextflow log` to compare runs:

```bash
nextflow log
```

Then examine specific runs:

```bash
nextflow log <run_name> -f hash,name,status
```

Compare hashes between runs to identify what changed.

### 5.3. Preserve work directories

If you're cleaning up work directories, be selective:

```bash
# Keep work directories from the last run
nextflow clean -but last -n  # dry run first
nextflow clean -but last
```

---

## Takeaway

You now know how to diagnose task failures by examining work directory files, recognize common error messages and their solutions, use debugging techniques like verbose logging and execution reports, and troubleshoot performance and caching issues.

### What's next?

Congratulations on completing the Nextflow Run course!
You now have the skills to run, configure, and troubleshoot Nextflow pipelines.

Head to the [Next steps](next_steps.md) page to learn about where to go from here.

---

## Quiz

<quiz>
A task fails with exit code 137. What is the most likely cause?
- [ ] The command wasn't found
- [ ] A syntax error in the script
- [x] The task ran out of memory
- [ ] A file permission issue

Learn more: [2.3. Out of memory (exit code 137)](#23-out-of-memory-exit-code-137)
</quiz>

<quiz>
Where should you look first when a task fails to find the error message?
- [ ] `.command.sh`
- [x] `.command.err`
- [ ] `.command.run`
- [ ] `.exitcode`

Learn more: [1.2. Examining the work directory](#12-examining-the-work-directory)
</quiz>

<quiz>
What command shows the complete configuration Nextflow will use?
- [ ] `nextflow run --show-config`
- [ ] `nextflow params`
- [x] `nextflow config`
- [ ] `nextflow info`

Learn more: [3.2. Inspect resolved configuration](#32-inspect-resolved-configuration)
</quiz>

<quiz>
Which of the following will NOT cause a task to re-run instead of using cached results?
- [ ] Modifying the process script
- [ ] Changing the container image
- [x] Running from a different terminal window
- [ ] Changing an input file's timestamp

Learn more: [5.1. Why tasks re-run unexpectedly](#51-why-tasks-re-run-unexpectedly)
</quiz>

<quiz>
You get "command not found" for a tool that should be in a container. What should you check first?
- [ ] That the tool is installed on your system
- [x] That `docker.enabled = true` is set in your config
- [ ] That the work directory exists
- [ ] That the input files are correct

Learn more: [2.1. Command not found](#21-command-not-found)
</quiz>
