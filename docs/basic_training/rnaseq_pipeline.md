---
description: Basic Nextflow Training Workshop
---

# Simple RNA-Seq workflow

To demonstrate a real-world biomedical scenario, we will implement a proof of concept RNA-Seq workflow which:

1. Indexes a transcriptome file
2. Performs quality controls
3. Performs quantification
4. Creates a MultiQC report

This will be done using a series of seven scripts, each of which builds on the previous to create a complete workflow. You can find these in the tutorial folder (`script1.nf` - `script7.nf`).

## Define the workflow parameters

Parameters are inputs and options that can be changed when the workflow is run.

The script `script1.nf` defines the workflow input parameters.

```groovy
params.reads = "$projectDir/data/ggal/gut_{1,2}.fq"
params.transcriptome_file = "$projectDir/data/ggal/transcriptome.fa"
params.multiqc = "$projectDir/multiqc"

println "reads: $params.reads"
```

Run it by using the following command:

```bash
nextflow run script1.nf
```

Try to specify a different input parameter in your execution command, for example:

```bash
nextflow run script1.nf --reads '/workspace/gitpod/nf-training/data/ggal/lung_{1,2}.fq'
```

### :material-progress-question: Exercises

!!! exercise

    Modify the `script1.nf` by adding a fourth parameter named `outdir` and set it to a default path that will be used as the workflow output directory.

    ??? result

        ```groovy
        params.reads = "$projectDir/data/ggal/gut_{1,2}.fq"
        params.transcriptome_file = "$projectDir/data/ggal/transcriptome.fa"
        params.multiqc = "$projectDir/multiqc"
        params.outdir = "results"
        ```

!!! exercise

    Modify `script1.nf` to print all of the workflow parameters by using a single `log.info` command as a [multiline string](https://www.nextflow.io/docs/latest/script.html#multi-line-strings) statement.

    !!! tip ""

        :material-lightbulb: See an example [here](https://github.com/nextflow-io/rnaseq-nf/blob/3b5b49f/main.nf#L41-L48).


    ??? result

        Add the following to your script file:

        ```groovy
        log.info """\
            R N A S E Q - N F   P I P E L I N E
            ===================================
            transcriptome: ${params.transcriptome_file}
            reads        : ${params.reads}
            outdir       : ${params.outdir}
            """
            .stripIndent(true)
        ```

### :material-check-all: Summary

In this step you have learned:

1. How to define parameters in your workflow script
2. How to pass parameters by using the command line
3. The use of `$var` and `${var}` variable placeholders
4. How to use multiline strings
5. How to use `log.info` to print information and save it in the log execution file

## Create a transcriptome index file

Nextflow allows the execution of any command or script by using a `process` definition.

A `process` is defined by providing three main declarations: the process [`input`](https://www.nextflow.io/docs/latest/process.html#inputs), [`output`](https://www.nextflow.io/docs/latest/process.html#outputs) and command [`script`](https://www.nextflow.io/docs/latest/process.html#script).

To add a transcriptome `INDEX` processing step, try adding the following code blocks to your `script1.nf`. Alternatively, these code blocks have already been added to `script2.nf`.

```groovy
/*
 * define the INDEX process that creates a binary index
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

This time the execution will work because it uses the Docker container `nextflow/rnaseq-nf` that is defined in the `nextflow.config` file in your current directory. If you are running this script locally then you will need to download Docker to your machine, log in and activate Docker, and allow the script to download the container containing the run scripts. You can learn more about Docker [here](https://www.nextflow.io/docs/latest/docker.html).

To avoid adding `-with-docker` each time you execute the script, add the following line to the `nextflow.config` file:

```groovy
docker.enabled = true
```

### :material-progress-question: Exercises

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

    If you have more CPUs available, try changing your script to request more resources for this process. For example, see the [directive docs](https://www.nextflow.io/docs/latest/process.html#cpus). `$task.cpus` is already specified in this script, so setting the number of CPUs as a directive will tell Nextflow how to execute this process, in terms of number of CPUs.

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

### :material-check-all: Summary

In this step you have learned:

1. How to define a process executing a custom command
2. How process inputs are declared
3. How process outputs are declared
4. How to print the content of a channel
5. How to access the number of available CPUs

## Collect read files by pairs

This step shows how to match **read** files into pairs, so they can be mapped by **Salmon**.

Edit the script `script3.nf` by adding the following statement as the last line of the file:

```groovy
read_pairs_ch.view()
```

Save it and execute it with the following command:

```bash
nextflow run script3.nf
```

It will print something similar to this:

```bash
[gut, [/.../data/ggal/gut_1.fq, /.../data/ggal/gut_2.fq]]
```

The above example shows how the `read_pairs_ch` channel emits tuples composed of two elements, where the first is the read pair prefix and the second is a list representing the actual files.

Try it again specifying different read files by using a glob pattern:

```bash
nextflow run script3.nf --reads 'data/ggal/*_{1,2}.fq'
```

!!! warning

    File paths that include one or more wildcards ie. `*`, `?`, etc., MUST be wrapped in single-quoted characters to avoid Bash expanding the glob.

### :material-progress-question: Exercises

!!! exercise

    Use the [set](https://www.nextflow.io/docs/latest/operator.html#set) operator in place of `=` assignment to define the `read_pairs_ch` channel.

    ??? result

        ```groovy
        Channel
            .fromFilePairs(params.reads)
            .set { read_pairs_ch }
        ```

!!! exercise

    Use the `checkIfExists` option for the [fromFilePairs](https://www.nextflow.io/docs/latest/channel.html#fromfilepairs) channel factory to check if the specified path contains file pairs.

    ??? result

        ```groovy
        Channel
            .fromFilePairs(params.reads, checkIfExists: true)
            .set { read_pairs_ch }
        ```

### :material-check-all: Summary

In this step you have learned:

1. How to use `fromFilePairs` to handle read pair files
2. How to use the `checkIfExists` option to check for the existence of input files
3. How to use the `set` operator to define a new channel variable

!!! info

    The declaration of a channel can be before the workflow scope or within it. As long as it is upstream of the process that requires the specific channel.

## Perform expression quantification

`script4.nf` adds a gene expression `QUANTIFICATION` process and a call to it within the workflow scope. Quantification requires the index transcriptome and RNA-Seq read pair fastq files.

In the workflow scope, note how the `index_ch` channel is assigned as output in the `INDEX` process.

Next, note that the first input channel for the `QUANTIFICATION` process is the previously declared `index_ch`, which contains the `path` to the `salmon_index`.

Also, note that the second input channel for the `QUANTIFICATION` process, is the `read_pair_ch` we just created. This being a `tuple` composed of two elements (a value: `sample_id` and a list of paths to the fastq reads: `reads`) in order to match the structure of the items emitted by the `fromFilePairs` channel factory.

Execute it by using the following command:

```bash
nextflow run script4.nf -resume
```

You will see the execution of the `QUANTIFICATION` process.

When using the `-resume` option, any step that has already been processed is skipped.

Try to execute the same script again with more read files, as shown below:

```bash
nextflow run script4.nf -resume --reads 'data/ggal/*_{1,2}.fq'
```

You will notice that the `QUANTIFICATION` process is executed multiple times.

Nextflow parallelizes the execution of your workflow simply by providing multiple sets of input data to your script.

!!! tip

    It may be useful to apply optional settings to a specific process using [directives](https://www.nextflow.io/docs/latest/process.html#directives) by specifying them in the process body.

### :material-progress-question: Exercises

!!! exercise

    Add a [tag](https://www.nextflow.io/docs/latest/process.html#tag) directive to the `QUANTIFICATION` process to provide a more readable execution log.

    ??? result

        Add the following before the input declaration:

        ```groovy
        tag "Salmon on $sample_id"
        ```

!!! exercise

    Add a [publishDir](https://www.nextflow.io/docs/latest/process.html#publishdir) directive to the `QUANTIFICATION` process to store the process results in a directory of your choice.

    ??? result

        Add the following before the `input` declaration in the `QUANTIFICATION` process:

        ```groovy
        publishDir params.outdir, mode: 'copy'
        ```

### :material-check-all: Summary

In this step you have learned:

1. How to connect two processes together by using the channel declarations
2. How to `resume` the script execution and skip cached steps
3. How to use the `tag` directive to provide a more readable execution output
4. How to use the `publishDir` directive to store a process results in a path of your choice

## Quality control

Next, we implement a `FASTQC` quality control step for your input reads (using the label `fastqc`). The inputs are the same as the read pairs used in the `QUANTIFICATION` step.

You can run it by using the following command:

```bash
nextflow run script5.nf -resume
```

Nextflow DSL2 knows to split the `reads_pair_ch` into two identical channels as they are required twice as an input for both of the `FASTQC` and the `QUANTIFICATION` process.

## MultiQC report

This step collects the outputs from the `QUANTIFICATION` and `FASTQC` processes to create a final report using the [MultiQC](http://multiqc.info/) tool.

Execute the next script with the following command:

```bash
nextflow run script6.nf -resume --reads 'data/ggal/*_{1,2}.fq'
```

It creates the final report in the `results` folder in the current `work` directory.

In this script, note the use of the [mix](https://www.nextflow.io/docs/latest/operator.html#mix) and [collect](https://www.nextflow.io/docs/latest/operator.html#collect) operators chained together to gather the outputs of the `QUANTIFICATION` and `FASTQC` processes as a single input. [Operators](https://www.nextflow.io/docs/latest/operator.html) can be used to combine and transform channels.

```groovy
MULTIQC(quant_ch.mix(fastqc_ch).collect())
```

We only want one task of MultiQC to be executed to produce one report. Therefore, we use the `mix` channel operator to combine the two channels followed by the `collect` operator, to return the complete channel contents as a single element.

### :material-check-all: Summary

In this step you have learned:

1. How to collect many outputs to a single input with the `collect` operator
2. How to `mix` two channels into a single channel
3. How to chain two or more operators together

## Handle completion event

This step shows how to execute an action when the workflow completes the execution.

Note that Nextflow processes define the execution of **asynchronous** tasks i.e. they are not executed one after another as if they were written in the workflow script in a common **imperative** programming language.

The script uses the `workflow.onComplete` event handler to print a confirmation message when the script completes.

Try to run it by using the following command:

```bash
nextflow run script7.nf -resume --reads 'data/ggal/*_{1,2}.fq'
```

## Email notifications

Send a notification email when the workflow execution completes using the `-N <email address>` command-line option.

Note: this requires the configuration of a SMTP server in the nextflow config file. Below is an example `nextflow.config` file showing the settings you would have to configure:

```groovy
mail {
    from = 'info@nextflow.io'
    smtp.host = 'email-smtp.eu-west-1.amazonaws.com'
    smtp.port = 587
    smtp.user = "xxxxx"
    smtp.password = "yyyyy"
    smtp.auth = true
    smtp.starttls.enable = true
    smtp.starttls.required = true
}
```

See [mail documentation](https://www.nextflow.io/docs/latest/mail.html#mail-configuration) for details.

## Custom scripts

Real-world workflows use a lot of custom user scripts (BASH, R, Python, etc.). Nextflow allows you to consistently use and manage these scripts. Simply put them in a directory named `bin` in the workflow project root. They will be automatically added to the workflow execution `PATH`.

For example, create a file named `fastqc.sh` with the following content:

```bash
#!/bin/bash
set -e
set -u

sample_id=${1}
reads=${2}

mkdir fastqc_${sample_id}_logs
fastqc -o fastqc_${sample_id}_logs -f fastq -q ${reads}
```

Save it, give execute permission, and move it into the `bin` directory as shown below:

```bash
chmod +x fastqc.sh
mkdir -p bin
mv fastqc.sh bin
```

Then, open the `script7.nf` file and replace the `FASTQC` process’ script with the following code:

```groovy
script:
"""
fastqc.sh "$sample_id" "$reads"
"""
```

Run it as before:

```bash
nextflow run script7.nf -resume --reads 'data/ggal/*_{1,2}.fq'
```

### :material-check-all: Summary

In this step you have learned:

1. How to write or use existing custom scripts in your Nextflow workflow.
2. How to avoid the use of absolute paths by having your scripts in the `bin/` folder.

## Metrics and reports

Nextflow can produce multiple reports and charts providing several runtime metrics and execution information.

Run the [rnaseq-nf](https://github.com/nextflow-io/rnaseq-nf) workflow previously introduced as shown below:

```bash
nextflow run rnaseq-nf -with-docker -with-report -with-trace -with-timeline -with-dag dag.png
```

The `-with-docker` option launches each task of the execution as a Docker container run command.

The `-with-report` option enables the creation of the workflow execution report. Open the file `report.html` with a browser to see the report created with the above command.

The `-with-trace` option enables the creation of a tab separated value (TSV) file containing runtime information for each executed task. Check the `trace.txt` for an example.

The `-with-timeline` option enables the creation of the workflow timeline report showing how processes were executed over time. This may be useful to identify the most time consuming tasks and bottlenecks. See an example at [this link](https://www.nextflow.io/docs/latest/tracing.html#timeline-report).

Finally, the `-with-dag` option enables the rendering of the workflow execution direct acyclic graph representation. Note: This feature requires the installation of [Graphviz](http://www.graphviz.org/) on your computer. See [here](https://www.nextflow.io/docs/latest/tracing.html#dag-visualisation) for further details. Then try running:

```bash
open dag.png
```

!!! warning

    Run time metrics may be incomplete for runs with short running tasks, as in the case of this tutorial.

!!! info

    You view the HTML files by right-clicking on the file name in the left side-bar and choosing the **Preview** menu item.

## Run a project from GitHub

Nextflow allows the execution of a workflow project directly from a GitHub repository (or similar services, e.g., BitBucket and GitLab).

This simplifies the sharing and deployment of complex projects and tracking changes in a consistent manner.

The following GitHub repository hosts a complete version of the workflow introduced in this tutorial: <https://github.com/nextflow-io/rnaseq-nf>

You can run it by specifying the project name and launching each task of the execution as a Docker container run command:

```bash
nextflow run nextflow-io/rnaseq-nf -with-docker
```

It automatically downloads the container and stores it in the `$HOME/.nextflow` folder.

Use the command `info` to show the project information:

```bash
nextflow info nextflow-io/rnaseq-nf
```

Nextflow allows the execution of a specific revision of your project by using the `-r` command line option. For example:

```bash
nextflow run nextflow-io/rnaseq-nf -r v2.1 -with-docker
```

Revision are defined by using Git tags or branches defined in the project repository.

Tags enable precise control of the changes in your project files and dependencies over time.

## More resources

-   [Nextflow documentation](http://docs.nextflow.io) - The Nextflow docs home.
-   [Nextflow patterns](https://github.com/nextflow-io/patterns) - A collection of Nextflow implementation patterns.
-   [CalliNGS-NF](https://github.com/CRG-CNAG/CalliNGS-NF) - A Variant calling workflow implementing GATK best practices.
-   [nf-core](http://nf-co.re/) - A community collection of production ready genomic workflows.
