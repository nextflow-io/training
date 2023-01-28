# Processes

In Nextflow, a `process` is the basic computing primitive to execute foreign functions (i.e., custom scripts or tools).

The `process` definition starts with the keyword `process`, followed by the process name and finally the process body delimited by curly brackets.

The `process` name is commonly written in upper case by convention.

A basic `process`, only using the `script` definition block, looks like the following:

    process SAYHELLO {

      script:
      """
      echo 'Hello world!'
      """
    }

In more complex examples, the process body can contain up to **five** definition blocks:

1. **Directives** are initial declarations that define optional settings

2. **Input** defines the expected input file(s) and the channel from where to find them

3. **Output** defines the expected output file(s) and the channel to send the data to

4. **When** is an optional clause statement to allow conditional processes

5. **Script** is a string statement that defines the command to be executed by the process

The full process syntax is defined as follows:

    process < name > {

      [ directives ]        //

      input:                //
      < process inputs >

      output:               //
      < process outputs >

      when:                 //
      < condition >

      [script|shell|exec]:  //
      """
      < user script to be executed >
      """
    }

-   Zero, one or more process directives

-   Zero, one or more process inputs

-   Zero, one or more process outputs

-   An optional boolean conditional to trigger the process execution

-   The command to be executed

## Script

The `script` block is a string statement that defines the command to be executed by the process.

A process can execute only one `script` block. It must be the last statement when the process contains input and output declarations.

The `script` block can be a single or a multi-line string. The latter simplifies the writing of non-trivial scripts composed of multiple commands spanning over multiple lines. For example:

    process EXAMPLE {

      script:
      """
      echo 'Hello world!\nHola mundo!\nCiao mondo!\nHallo Welt!' > file
      cat file | head -n 1 | head -c 5 > chunk_1.txt
      gzip -c chunk_1.txt  > chunk_archive.gz
      """
    }

    workflow {
      EXAMPLE()
    }

By default, the `process` command is interpreted as a **Bash** script. However, any other scripting language can be used by simply starting the script with the corresponding [Shebang](<https://en.wikipedia.org/wiki/Shebang_(Unix)>) declaration. For example:

    process PYSTUFF {

      script:
      """
      #!/usr/bin/env python

      x = 'Hello'
      y = 'world!'
      print ("%s - %s" % (x,y))
      """
    }

    workflow {
      PYSTUFF()
    }

Multiple programming languages can be used within the same workflow script. However, for large chunks of code it is better to save them into separate files and invoke them from the process script. One can store the specific scripts in the `./bin/` folder.

### Script parameters

Script parameters (`params`) can be defined dynamically using variable values. For example:

    params.data = 'World'

    process FOO {

      script:
      """
      echo Hello $params.data
      """
    }

    workflow {
      FOO()
    }

A process script can contain any string format supported by the Groovy programming language. This allows us to use string interpolation as in the script above or multiline strings. Refer to [String interpolation](#groovy.adoc#_string_interpolation) for more information.

Since Nextflow uses the same Bash syntax for variable substitutions in strings, Bash environment variables need to be escaped using the `\` character.

    process FOO {

      script:
      """
      echo "The current directory is \$PWD"
      """
    }

    workflow {
      FOO()
    }

It can be tricky to write a script uses many Bash variables. One possible alternative is to use a script string delimited by single-quote characters

    process BAR {

      script:
      """
      echo $PATH | tr : '\\n'
      """
    }

    workflow {
      BAR()
    }

However, this blocks the usage of Nextflow variables in the command script.

Another alternative is to use a `shell` statement instead of `script` and use a different syntax for Nextflow variables, e.g., `!{..}`. This allows the use of both Nextflow and Bash variables in the same script.

    params.data = 'le monde'

    process BAZ {

      shell:
      '''
      X='Bonjour'
      echo $X !{params.data}
      '''
    }

    workflow {
      BAZ()
    }

### Conditional script

The process script can also be defined in a completely dynamic manner using an `if` statement or any other expression for evaluating a string value. For example:

    params.compress = 'gzip'
    params.file2compress = "$baseDir/data/ggal/transcriptome.fa"

    process FOO {

      input:
      path file

      script:
      if( params.compress == 'gzip' )
        """
        gzip -c $file > ${file}.gz
        """
      else if( params.compress == 'bzip2' )
        """
        bzip2 -c $file > ${file}.bz2
        """
      else
        throw new IllegalArgumentException("Unknown aligner $params.compress")
    }

    workflow {
      FOO(params.file2compress)
    }

## Inputs

Nextflow processes are isolated from each other but can communicate between themselves by sending values through channels.

Inputs implicitly determine the dependencies and the parallel execution of the process. The process execution is fired each time _new_ data is ready to be consumed from the input channel:

![](img/channel-process.png)

The `input` block defines which channels the `process` is expecting to receive data from. You can only define one `input` block at a time, and it must contain one or more input declarations.

The `input` block follows the syntax shown below:

    input:
      <input qualifier> <input name>

### Input values

The `val` qualifier allows you to receive data of any type as input. It can be accessed in the process script by using the specified input name, as shown in the following example:

    num = Channel.of( 1, 2, 3 )

    process BASICEXAMPLE {
      debug true

      input:
      val x

      script:
      """
      echo process job $x
      """
    }

    workflow {
      myrun = BASICEXAMPLE(num)
    }

In the above example the process is executed three times, each time a value is received from the channel `num` and used to process the script. Thus, it results in an output similar to the one shown below:

    process job 3
    process job 1
    process job 2

The channel guarantees that items are delivered in the same order as they have been sent - but - since the process is executed in a parallel manner, there is no guarantee that they are processed in the same order as they are received.

### Input files

The `path` qualifier allows the handling of file values in the process execution context. This means that Nextflow will stage it in the process execution directory, and it can be accessed in the script by using the name specified in the input declaration.

    reads = Channel.fromPath( 'data/ggal/*.fq' )

    process FOO {
      debug true

      input:
      path 'sample.fastq'

      script:
      """
      ls sample.fastq
      """
    }

    workflow {
      result = FOO(reads)
    }

The input file name can also be defined using a variable reference as shown below:

    reads = Channel.fromPath( 'data/ggal/*.fq' )

    process FOO {
      debug true

      input:
      path sample

      script:
      """
      ls  $sample
      """
    }

    workflow {
      result = FOO(reads)
    }

The same syntax is also able to handle more than one input file in the same execution and only requires changing the channel composition.

    reads = Channel.fromPath( 'data/ggal/*.fq' )

    process FOO {
      debug true

      input:
      path sample

      script:
      """
      ls -lh $sample
      """
    }

    workflow {
      FOO(reads.collect())
    }

In the past, the `file` qualifier was used for files, but the `path` qualifier should be preferred over file to handle process input files when using Nextflow 19.10.0 or later. When a process declares an input file, the corresponding channel elements must be **file** objects created with the file helper function from the file specific channel factories (e.g., `Channel.fromPath` or `Channel.fromFilePairs`).

Write a script that creates a channel containing all read files matching the pattern `data/ggal/*_1.fq` followed by a process that concatenates them into a single file and prints the first 20 lines.

    params.reads = "$baseDir/data/ggal/*_1.fq"

    Channel
      .fromPath( params.reads )
      .set { read_ch }

    process CONCATENATE {
      tag "Concat all files"

      input:
      path '*'

      output:
      path 'top_10_lines'

      script:
      """
      cat * > concatenated.txt
      head -n 20 concatenated.txt > top_10_lines
      """
    }

    workflow {
      concat_ch = CONCATENATE(read_ch.collect())
      concat_ch.view()
    }

### Combine input channels

A key feature of processes is the ability to handle inputs from multiple channels. However, it’s important to understand how channel contents and their semantics affect the execution of a process.

Consider the following example:

    ch1 = Channel.of(1,2,3)
    ch2 = Channel.of('a','b','c')

    process FOO {
      debug true

      input:
      val x
      val y

      script:
       """
       echo $x and $y
       """
    }

    workflow {
      FOO(ch1, ch2)
    }

Both channels emit three values, therefore the process is executed three times, each time with a different pair:

-   (1, a)

-   (2, b)

-   (3, c)

What is happening is that the process waits until there’s a complete input configuration, i.e., it receives an input value from all the channels declared as input.

When this condition is verified, it consumes the input values coming from the respective channels, spawns a task execution, then repeats the same logic until one or more channels have no more content.

This means channel values are consumed serially one after another and the first empty channel causes the process execution to stop, even if there are other values in other channels.

**So what happens when channels do not have the same cardinality (i.e., they emit a different number of elements)?**

For example:

    input1 = Channel.of(1,2)
    input2 = Channel.of('a','b','c','d')

    process FOO {
      debug true

      input:
      val x
      val y

      script:
       """
       echo $x and $y
       """
    }

    workflow {
      FOO(input1, input2)
    }

In the above example, the process is only executed twice because the process stops when a channel has no more data to be processed.

However, what happens if you replace value x with a `value` channel?

Compare the previous example with the following one :

    input1 = Channel.value(1)
    input2 = Channel.of('a','b','c')

    process BAR {
      debug true

      input:
      val x
      val y

      script:
       """
       echo $x and $y
       """
    }

    workflow {
      BAR(input1, input2)
    }

    1 and b
    1 and a
    1 and c

This is because _value_ channels can be consumed multiple times and do not affect process termination.

Write a process that is executed for each read file matching the pattern `data/ggal/*_1.fq` and use the same `data/ggal/transcriptome.fa` in each execution.

    params.reads = "$baseDir/data/ggal/*_1.fq"
    params.transcriptome_file = "$baseDir/data/ggal/transcriptome.fa"

    Channel
        .fromPath( params.reads )
        .set { read_ch }

    process COMMAND {
      tag "Run_command"

      input:
      path reads
      path transcriptome

      output:
      path result

      script:
      """
      echo your_command $reads $transcriptome > result
      """
    }

    workflow {
      concat_ch = COMMAND(read_ch, params.transcriptome_file)
      concat_ch.view()
    }

### Input repeaters

The `each` qualifier allows you to repeat the execution of a process for each item in a collection every time new data is received. For example:

    sequences = Channel.fromPath('data/prots/*.tfa')
    methods = ['regular', 'espresso', 'psicoffee']

    process ALIGNSEQUENCES {
      debug true

      input:
      path seq
      each mode

      script:
      """
      echo t_coffee -in $seq -mode $mode
      """
    }

    workflow {
      ALIGNSEQUENCES(sequences, methods)
    }

In the above example, every time a file of sequences is received as an input by the process, it executes three tasks, each running a different alignment method set as a `mode` variable. This is useful when you need to repeat the same task for a given set of parameters.

Extend the previous example so a task is executed for each read file matching the pattern `data/ggal/*_1.fq` and repeat the same task with both `salmon` and `kallisto`.

    params.reads = "$baseDir/data/ggal/*_1.fq"
    params.transcriptome_file = "$baseDir/data/ggal/transcriptome.fa"
    methods= ['salmon', 'kallisto']

    Channel
        .fromPath( params.reads )
        .set { read_ch }

    process COMMAND {
      tag "Run_command"

      input:
      path reads
      path transcriptome
      each mode

      output:
      path result

      script:
      """
      echo $mode $reads $transcriptome > result
      """
    }

    workflow {
      concat_ch = COMMAND(read_ch , params.transcriptome_file, methods)
      concat_ch
          .view { "To run : ${it.text}" }
    }

## Outputs

The _output_ declaration block defines the channels used by the process to send out the results produced.

Only one output block, that can contain one or more output declaration, can be defined. The output block follows the syntax shown below:

    output:
      <output qualifier> <output name> , emit: <output channel>

### Output values

The `val` qualifier specifies a defined _value_ in the script context. Values are frequently defined in the _input_ and/or _output_ declaration blocks, as shown in the following example:

    methods = ['prot','dna', 'rna']

    process FOO {

      input:
      val x

      output:
      val x

      script:
      """
      echo $x > file
      """
    }

    workflow {
      receiver_ch = FOO(Channel.of(methods))
      receiver_ch.view { "Received: $it" }
    }

### Output files

The `path` qualifier specifies one or more files produced by the process into the specified channel as an output.

    process RANDOMNUM {

        output:
        path 'result.txt'

        script:
        """
        echo $RANDOM > result.txt
        """
    }


    workflow {
      receiver_ch = RANDOMNUM()
      receiver_ch.view { "Received: " + it.text }
    }

In the above example the process `RANDOMNUM` creates a file named `result.txt` containing a random number.

Since a file parameter using the same name is declared in the output block, the file is sent over the `receiver_ch` channel when the task is complete. A downstream `process` declaring the same channel as _input_ will be able to receive it.

### Multiple output files

When an output file name contains a wildcard character (`*` or `?`) it is interpreted as a [glob](http://docs.oracle.com/javase/tutorial/essential/io/fileOps.html#glob) path matcher. This allows us to _capture_ multiple files into a list object and output them as a sole emission. For example:

    process SPLITLETTERS {

        output:
        path 'chunk_*'

        """
        printf 'Hola' | split -b 1 - chunk_
        """
    }

    workflow {
        letters = SPLITLETTERS()
        letters
            .flatMap()
            .view { "File: ${it.name} => ${it.text}" }
    }

Prints the following:

    File: chunk_aa => H
    File: chunk_ab => o
    File: chunk_ac => l
    File: chunk_ad => a

Some caveats on glob pattern behavior:

-   Input files are not included in the list of possible matches

-   Glob pattern matches both files and directory paths

-   When a two stars pattern \`\`\*\*\`\` is used to recourse across directories, only file paths are matched i.e., directories are not included in the result list.

Remove the `flatMap` operator and see out the output change. The documentation for the `flatMap` operator is available at [this link](https://www.nextflow.io/docs/latest/operator.html#flatmap).

    File: [chunk_aa, chunk_ab, chunk_ac, chunk_ad] => [H, o, l, a]

### Dynamic output file names

When an output file name needs to be expressed dynamically, it is possible to define it using a dynamic string that references values defined in the input declaration block or in the script global context. For example:

    species = ['cat','dog', 'sloth']
    sequences = ['AGATAG','ATGCTCT', 'ATCCCAA']

    Channel.fromList(species)
           .set { species_ch }

    process ALIGN {

      input:
      val x
      val seq

      output:
      path "${x}.aln"

      script:
      """
      echo align -in $seq > ${x}.aln
      """
    }

    workflow {
      genomes = ALIGN( species_ch, sequences )
      genomes.view()
    }

In the above example, each time the process is executed an alignment file is produced whose name depends on the actual value of the `x` input.

### Composite inputs and outputs

So far we have seen how to declare multiple input and output channels that can handle one value at a time. However, Nextflow can also handle a _tuple_ of values.

The input and output declarations for tuples must be declared with a `tuple` qualifier followed by the definition of each element in the tuple.

    reads_ch = Channel.fromFilePairs('data/ggal/*_{1,2}.fq')

    process FOO {

      input:
        tuple val(sample_id), path(sample_id)

      output:
        tuple val(sample_id), path('sample.bam')

      script:
      """
        echo your_command_here --reads $sample_id > sample.bam
      """
    }

    workflow {
      bam_ch = FOO(reads_ch)
      bam_ch.view()
    }

In previous versions of Nextflow `tuple` was called `set` but it was used the same way with the same semantic.

Modify the script of the previous exercise so that the _bam_ file is named as the given `sample_id`.

    reads_ch = Channel.fromFilePairs('data/ggal/*_{1,2}.fq')

    process FOO {

      input:
        tuple val(sample_id), path(sample_files)

      output:
        tuple val(sample_id), path("${sample_id}.bam")

      script:
      """
        echo your_command_here --reads $sample_id > ${sample_id}.bam
      """
    }

    workflow {
      bam_ch = FOO(reads_ch)
      bam_ch.view()
    }

## When

The `when` declaration allows you to define a condition that must be verified in order to execute the process. This can be any expression that evaluates a boolean value.

It is useful to enable/disable the process execution depending on the state of various inputs and parameters. For example:

    params.dbtype = 'nr'
    params.prot = 'data/prots/*.tfa'
    proteins = Channel.fromPath(params.prot)

    process FIND {
      debug true

      input:
      path fasta
      val type

      when:
      fasta.name =~ /^BB11.*/ && type == 'nr'

      script:
      """
      echo blastp -query $fasta -db nr
      """
    }

    workflow {
      result = FIND(proteins, params.dbtype)
    }

## Directives

Directive declarations allow the definition of optional settings that affect the execution of the current process without affecting the _semantic_ of the task itself.

They must be entered at the top of the process body, before any other declaration blocks (i.e., `input`, `output`, etc.).

Directives are commonly used to define the amount of computing resources to be used or other meta directives that allow the definition of extra configuration of logging information. For example:

    process FOO {
      cpus 2
      memory 1.GB
      container 'image/name'

      script:
      """
      echo your_command --this --that
      """
    }

The complete list of directives is available [at this link](https://www.nextflow.io/docs/latest/process.html#directives).

<table>
<caption>Commonly used directives</caption>
<colgroup>
<col style="width: 100%" />
</colgroup>
<tbody>
<tr class="odd">
<td style="text-align: left;"><p>Name</p></td>
</tr>
<tr class="even">
<td style="text-align: left;"><p>Description</p></td>
</tr>
<tr class="odd">
<td style="text-align: left;"><p><a href="https://www.nextflow.io/docs/latest/process.html#cpus">cpus</a></p></td>
</tr>
<tr class="even">
<td style="text-align: left;"><p>Allows you to define the number of (logical) CPUs required by the process’ task.</p></td>
</tr>
<tr class="odd">
<td style="text-align: left;"><p><a href="https://www.nextflow.io/docs/latest/process.html#time">time</a></p></td>
</tr>
<tr class="even">
<td style="text-align: left;"><p>Allows you to define how long a process is allowed to run (e.g., time <em>1h</em>: 1 hour, <em>1s</em> 1 second, <em>1m</em> 1 minute, <em>1d</em> 1 day).</p></td>
</tr>
<tr class="odd">
<td style="text-align: left;"><p><a href="https://www.nextflow.io/docs/latest/process.html#memory">memory</a></p></td>
</tr>
<tr class="even">
<td style="text-align: left;"><p>Allows you to define how much memory the process is allowed to use (e.g., <em>2 GB</em> is 2 GB). Can also use B, KB,MB,GB and TB.</p></td>
</tr>
<tr class="odd">
<td style="text-align: left;"><p><a href="https://www.nextflow.io/docs/latest/process.html#disk">disk</a></p></td>
</tr>
<tr class="even">
<td style="text-align: left;"><p>Allows you to define how much local disk storage the process is allowed to use.</p></td>
</tr>
<tr class="odd">
<td style="text-align: left;"><p><a href="https://www.nextflow.io/docs/latest/process.html#tag">tag</a></p></td>
</tr>
<tr class="even">
<td style="text-align: left;"><p>Allows you to associate each process execution with a custom label to make it easier to identify them in the log file or the trace execution report.</p></td>
</tr>
</tbody>
</table>

Commonly used directives

## Organize outputs

### PublishDir directive

Given each process is being executed in separate temporary `work/` folder (e.g., work/f1/850698…; work/g3/239712…; etc.), we may want to save important, non-intermediary, and/or final files in a results folder.

Remember to delete the work folder from time to time to clear your intermediate files and stop them from filling your computer!

To store our workflow result files, we need to explicitly mark them using the directive [publishDir](https://www.nextflow.io/docs/latest/process.html#publishdir) in the process that’s creating the files. For example:

    params.outdir = 'my-results'
    params.prot = 'data/prots/*.tfa'
    proteins = Channel.fromPath(params.prot)


    process BLASTSEQ {
        publishDir "$params.outdir/bam_files", mode: 'copy'

        input:
        path fasta

        output:
        path ('*.txt')

        script:
        """
        echo blastp $fasta > ${fasta}_result.txt
        """
    }

    workflow {
      blast_ch = BLASTSEQ(proteins)
      blast_ch.view()
    }

The above example will copy all blast script files created by the `BLASTSEQ` task into the directory path `my-results`.

The publish directory can be local or remote. For example, output files could be stored using an [AWS S3 bucket](https://aws.amazon.com/s3/) by using the `s3://` prefix in the target path.

### Manage semantic sub-directories

You can use more than one `publishDir` to keep different outputs in separate directories. For example:

    params.reads = 'data/reads/*_{1,2}.fq.gz'
    params.outdir = 'my-results'

    samples_ch = Channel.fromFilePairs(params.reads, flat: true)

    process FOO {
      publishDir "$params.outdir/$sampleId/", pattern: '*.fq'
      publishDir "$params.outdir/$sampleId/counts", pattern: "*_counts.txt"
      publishDir "$params.outdir/$sampleId/outlooks", pattern: '*_outlook.txt'

      input:
        tuple val(sampleId), path('sample1.fq.gz'), path('sample2.fq.gz')

      output:
        path "*"

      script:
      """
        < sample1.fq.gz zcat > sample1.fq
        < sample2.fq.gz zcat > sample2.fq

        awk '{s++}END{print s/4}' sample1.fq > sample1_counts.txt
        awk '{s++}END{print s/4}' sample2.fq > sample2_counts.txt

        head -n 50 sample1.fq > sample1_outlook.txt
        head -n 50 sample2.fq > sample2_outlook.txt
      """
    }

    workflow {
      out_channel = FOO(samples_ch)
    }

The above example will create an output structure in the directory `my-results`, that contains a separate sub-directory for each given sample ID, each containing the folders `counts` and `outlooks`.
