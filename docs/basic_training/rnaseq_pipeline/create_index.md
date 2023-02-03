---
title: Transcriptome index
---

# Create a transcriptome index file

Nextflow allows the execution of any command or script by using a `process` definition.

A `process` is defined by providing three main declarations: the process [`input`](https://www.nextflow.io/docs/latest/process.html#inputs), [`output`](https://www.nextflow.io/docs/latest/process.html#outputs) and command [`script`](https://www.nextflow.io/docs/latest/process.html#script).

To add a transcriptome `INDEX` processing step, try adding the following code blocks to your `script1.nf`. Alternatively, these code blocks have already been added to `script2.nf`.

```groovy
/*
 * define the `index` process that creates a binary index
 * given the transcriptome file
 */
process INDEX {
  input:
  path transcriptome

  output:
  path 'salmon_index'

  script:
  """
  salmon index --threads $task.cpus -t $transcriptome -i salmon_index
  """
}
```

Additionally, add a workflow scope containing an input channel definition and the index process:

```groovy
workflow {
    index_ch = INDEX(params.transcriptome_file)
}
```

Here, the `params.transcriptome_file` parameter is used as the input for the `INDEX` process. The `INDEX` process (using the `salmon` tool) creates `salmon_index`, an indexed transcriptome that is passed as an output to the `index_ch` channel.

!!! info

    The `input` declaration defines a `transcriptome` path variable which is used in the `script` as a reference (using the dollar symbol) in the Salmon command line.

!!! warning

    Resource requirements such as CPUs and memory limits can change with different workflow executions and platforms. Nextflow can use `$task.cpus` as a variable for the number of CPUs. See [process directives documentation](https://www.nextflow.io/docs/latest/process.html#directives) for more details.

Run it by using the command:

```bash
nextflow run script2.nf
```

The execution will fail because `salmon` is not installed in your environment.

Add the command line option `-with-docker` to launch the execution through a Docker container, as shown below:

```bash
nextflow run script2.nf -with-docker
```

This time the execution will work because it uses the Docker container `nextflow/rnaseq-nf` that is defined in the `nextflow.config` file in your current directory. If you are running this script locally then you will need to download docker to your machine, log in and activate docker, and allow the script to download the container containing the run scripts. You can learn more about docker [here](https://www.nextflow.io/docs/latest/docker.html).

To avoid adding `-with-docker` each time you execute the script, add the following line to the `nextflow.config` file:

```groovy
docker.enabled = true
```

## :material-progress-question: Exercises

!!! exercise

    Enable the Docker execution by default by adding the above setting in the `nextflow.config` file.

!!! exercise

    Print the output of the `index_ch` channel by using the [view](https://www.nextflow.io/docs/latest/operator.html#view) operator.

    ??? result

        Add the following to the end of your workflow block in your script file

        ```groovy
        index_ch.view()
        ```

!!! exercise

    If you have more CPUs available, try changing your script to request more resources for this process. For example, see the [directive docs](https://www.nextflow.io/docs/latest/process.html#cpus). `$task.cpus` is already specified in this script, so setting the number of CPUs as a directive will tell Nextflow to run this job.

    ??? result

        Add `cpus 2` to the top of the index process:

        ```groovy
        process INDEX {
            cpus 2
            input:
            ...
        ```

        Then check it worked by looking at the script executed in the work directory. Look for the hexadecimal (e.g. `work/7f/f285b80022d9f61e82cd7f90436aa4/`), Then `cat` the `.command.sh` file.

!!! exercise "Bonus Exercise"

    Use the command `tree work` to see how Nextflow organizes the process work directory. Check [here](https://www.tecmint.com/linux-tree-command-examples/) if you need to download `tree`.

    ??? result

        It should look something like this:

        ```
        work
        ├── 17
        │   └── 263d3517b457de4525513ae5e34ea8
        │       ├── index
        │       │   ├── complete_ref_lens.bin
        │       │   ├── ctable.bin
        │       │   ├── ctg_offsets.bin
        │       │   ├── duplicate_clusters.tsv
        │       │   ├── eqtable.bin
        │       │   ├── info.json
        │       │   ├── mphf.bin
        │       │   ├── pos.bin
        │       │   ├── pre_indexing.log
        │       │   ├── rank.bin
        │       │   ├── refAccumLengths.bin
        │       │   ├── ref_indexing.log
        │       │   ├── reflengths.bin
        │       │   ├── refseq.bin
        │       │   ├── seq.bin
        │       │   └── versionInfo.json
        │       └── transcriptome.fa -> /workspace/Gitpod_test/data/ggal/transcriptome.fa
        ├── 7f
        ```

## :material-check-all: Summary

In this step you have learned:

1. How to define a process executing a custom command
2. How process inputs are declared
3. How process outputs are declared
4. How to print the content of a channel
5. How to access the number of available CPUs
