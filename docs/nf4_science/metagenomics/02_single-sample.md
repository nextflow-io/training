# Part 2: Single-sample analysis

Now that we have a clear overview of what we want to achieve, we can start writing code!
We are going to do this in four steps:

1. Create the modules that are going to perform each task
2. Create a `workflow.nf` containing the core logic of the pipeline
3. Create a `main.nf` entrypoint workflow
4. Create a `nextflow.config` configuration file

Then we will test that we can run the pipeline and that it produces the expected outputs.

---

## 1. Create the modules

You can think of modules as building blocks that we can string together to form an assembly line.
Each module corresponds to a step in the pipeline.

### 1.1. Remove host sequence with Bowtie2

As noted in the method overview, the objective of this module is to clean the reads by aligning them against a reference genome and remove any reads that map to that reference.

Let's create the `bowtie2.nf` file inside the `modules/` folder and write the following code:

```groovy title="modules/bowtie2.nf" linenums="1"
process BOWTIE2 {
    tag "${sample_id}"
    publishDir "$params.outdir/${sample_id}", pattern: "*.sam", mode:'copy'
    container "community.wave.seqera.io/library/bowtie2:2.5.4--d51920539234bea7"

    input:
    tuple val(sample_id), path(reads)
    path bowtie2_index

    output:
    tuple val("${sample_id}"), path("${sample_id}.1"), path("${sample_id}.2"), path("${sample_id}.sam")

    script:
    """
    export BOWTIE2_INDEXES=/workspaces/training/nf4-science/metagenomics/data/yeast
    bowtie2 -x $bowtie2_index -1 ${reads[0]} -2 ${reads[1]} -p 2 -S ${sample_id}.sam --un-conc-gz ${sample_id}
    """
}
```

Let's take a moment to break down what we are seeing here.

First, let's look at the 'housekeeping' parts of the process, which you should be familiar with if you worked through the beginner Nextflow training.

- The process name is `BOWTIE2`, this is important when creating the workflow file.
- The `tag` directive is used to indicate which sample is being processed at a determined moment.
  This will be useful when running the pipeline.
- `publishDir` points out to the directory where the output is stored.
  In this case we are taking the desired output path (supplied as a parameter) and creating subdirectories based on the sample names where each `.sam` file will be output.
- `container` indicates the docker container with which the process will be run.
  We have retrieved all the containers from [Seqera Containers](https://seqera.io/containers/).

Now let's have a look at the interesting bits!

#### 1.1.1. Process inputs

The `input` for this process will be the `sample id`, the paired-end `reads`, as well as the path to the `bowtie index`.

We'll go over these in more detail when we create the `main.nf` file.

#### 1.1.2. Process outputs

This process will produce a tuple containing the `sample id`, the path to the cleaned `reads` and the path to the `*.sam` file.

The SAM file contains the information about the alignment of the sequences against the reference indexed genome.
To learn more about this format, please see the [SAM format specification](https://samtools.github.io/hts-specs/SAMv1.pdf).

#### 1.1.3. Command script

The `script` block contains two commands:

1. An environment variable is exported pointing to the directory where the indexed genome is stored.
   This is required by Bowtie2.
   If you plan to use another reference genome, you will need to adapt this.
   **[TODO: parameterize this]**

2. The Bowtie2 command that includes again the path to the indexed genome.
   The `-1` and `-2` arguments capture the path to forward and reverse reads, respectively.The `-S` argument sets the desired format output.
   The `--un-conc-gz` argument is used to write out any paired-end reads that fail to align concordantly to the reference genome.
   To learn more about this parameter, see the [Bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml) manual.

### 1.2. Apply taxonomic classification with Kraken2

Now, let's create the module for our taxonomic classifier.

Let's create `kraken2.nf` file inside the `modules/` folder and write this code:

```groovy title="modules/kraken2.nf" linenums="1"
process KRAKEN2 {
	tag "${sample_id}"
	publishDir "$params.outdir/${sample_id}", mode:'copy'
	container "community.wave.seqera.io/library/kraken2:2.14--83aa57048e304f01"

	input:
	tuple val(sample_id), path(reads_1), path(reads_2), path(sam)
	path kraken2_db

	output:
	tuple val("${sample_id}"), path("${sample_id}.k2report"), path("${sample_id}.kraken2")

	script:
	"""
	kraken2 --db $kraken2_db --threads 2 \
	--report ${sample_id}.k2report \
	--report-minimizer-data \
	--minimum-hit-groups 2 \
	--gzip-compressed \
	--paired \
	${reads_1} ${reads_2} > ${sample_id}.kraken2
	"""
}
```

At first glance, you can see that it follows the same structure as the previous process.
The directives `tag`, `publishDir` and `container` play the same role.

#### 1.2.1. Process inputs

The `input` in this case is a tuple containing the `sample id`, the cleaned reads and the `*.sam` file.
This one is declared just to maintain the correspondence between the output from Bowtie2 and the Kraken2 input.

Also, the path to the Kraken database is declared.
We'll go over this in more detail when create the `main.nf` file.

#### 1.2.2. Process outputs

The output from this process will be a tuple containing the `sample id`, the path to the `.k2report` file, as well as the path to the `.kraken2` file.

#### 1.2.3. Command script

This is a standard Kraken2 command that specifies the following arguments:

- the path to the database
- the number of threads to use
- the path to the reports Kraken2 generates
- the minimum number of 'hit groups' needed to make a classification call
- the flag to report the minimizers and distinct minimizer count
- the parameters that indicate that the received reads are paired-end and compressed in a _.gz_ format.

You are strongly encouraged to check out both the [protocol](https://www.nature.com/articles/s41596-022-00738-y) that documents this methodology, as well as the [source publication](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1891-0) to understand how Kraken2 works and the type of files it generates.
This [article](https://homolog.us/blogs/bioinfo/2017/10/25/intro-minimizer/) about minimizers may also be useful.

### 1.3. Apply Bayesian re-estimation of species abundance with Bracken

Let's create the file `bracken.nf` inside `modules/` and write the following code:

```groovy title="modules/bracken.nf" linenums="1"
process BRACKEN {
	tag "${sample_id}"
	publishDir "$params.outdir/${sample_id}", mode:'copy'
	container "community.wave.seqera.io/library/bracken:3.1--22a4e66ce04c5e01"

	input:
	tuple val(sample_id), path(k2report), path(kraken2)
	path kraken2_db

	output:
	tuple val("${sample_id}"), path("${sample_id}.breport"), path("${sample_id}.bracken")

	script:
	"""
	bracken -d $kraken2_db \
	-i ${k2report} -r 250 -l S -t 10 \
	-o ${sample_id}.bracken \
	-w ${sample_id}.breport
	"""
}
```

Once again, this module follows the same structure as the previous process.
The directives `tag`, `publishDir` and `container` play the same role.

#### 1.2.1. Process inputs

The declared `input` for this Bracken process matches the output generated by Kraken2.
It expects to receive a tuple containing the `sample id` and the paths to `.k2report` and `kraken2` reports.

#### 1.2.2. Process outputs

Similarly, the `output` establishes that this process will generate a tuple with the `sample id` and the two types of `reports` Bracken produces.
You can find more information about these reports in the [GitHub repository for Bracken](https://github.com/jenniferlu717/Bracken).

#### 1.2.3. Command script

The Bracken command specifies the following arguments:

- the path to the database
- the input represented by the Kraken2 report
- the length of the reads
- the taxonomic level at which we want to re-estimate the abundance (species in this case)
- the required number of reads before abundance estimation to perform re-estimation
- the format for the output reports

### 1.4. Generate plots with Krona

Krona plots are interactive pie charts that are widely used in the metagenomics field.

To generate them, we will need two distinct processes; one to extract species abundance metrics from the Bracken report, and another to render the Krona plot itself based on that file.

#### 1.4.1. Extract the species abundance metrics

Let's create the file `kReport2Krona.nf` inside `modules` and write the following code:

```groovy title="modules/kReport2Krona.nf" linenums="1"
process K_REPORT_TO_KRONA {
	tag "${sample_id}"
	publishDir "$params.outdir/${sample_id}", mode:'copy'
	container "community.wave.seqera.io/library/krakentools:1.2--db94e0b19cfa397b"

	input:
	tuple val(sample_id), path(b_report), path(bracken)

	output:
	tuple val("${sample_id}"), path("${sample_id}.b.krona.txt")

	script:
	"""
	kreport2krona.py -r ${b_report} \
	-o ${sample_id}.b.krona.txt \
	--no-intermediate-ranks
	"""
}
```

As before, this process takes as input the exact output from Bracken in a tuple form, and it is going to generate another tuple containing `sample id` and the `*.txt` file with the abundance metrics.

The command will run a script called `kreport2krona.py` that is provided as part of [KrakenTools](https://github.com/jenniferlu717/KrakenTools) and uses the Bracken report as input to produce a plain text file containing the species abundance information.

Note that the command includes an additional flag to leave out non-traditional ranks.
For more information about that option, see the KrakenTools documentation referenced above.

#### 1.4.2. Extract the species abundance metrics

Let's create the file `ktImportText.nf` inside `modules/` and write the following code:

```groovy title="modules/ktImportText.nf" linenums="1"
process KT_IMPORT_TEXT {
	tag "${sample_id}"
	publishDir "$params.outdir/${sample_id}", mode:'copy'
	container "community.wave.seqera.io/library/krona:2.8.1--2f750080982f027e"

	input:
	tuple val(sample_id), path(krona_txt)

	output:
	path "${sample_id}.krona.html"

	script:
	"""
	ktImportText ${krona_txt} \
	-o ${sample_id}.krona.html
	"""
}
```

This process takes as input the plain `*.txt` in a tuple with the sample ID, and uses Krona to render the Krona plot in a self-contained `*.html` file, which is the final output of our pipeline at this stage.

The Krona plot can be opened in a standard modern browser, just like the one you are using to follow this tutorial.

This completes the implementation of the modules we need for the single-sample use case, so we can now move on to creating the files that will describe the workflow and control its execution.

---

## 2. Create the `workflow.nf`

There are several strategies for building Nextflow workflows.
Here we are going to use a composable workflow structure as described in the [`Workflows of Workflows`](../../side_quests/workflows_of_workflows.md) Side Quest, which uses an entrypoint workflow in a `main.nf` file, an embedded `workflow.nf` file containing the core logic of the workflow, and the `take` syntax to declare inputs.

Let's start by creating the `workflow.nf` file that will the core logic of the workflow.
Note that we create this file in our current directory, NOT in the `modules/` directory.

We start by adding the following code to import modules:

```groovy title="workflow.nf" linenums="1"
/*
 * required tasks
 */
include { BOWTIE2							} from './modules/bowtie2.nf'
include { KRAKEN2							}	from './modules/kraken2.nf'
include { BRACKEN							} from './modules/bracken.nf'
include { K_REPORT_TO_KRONA		} from './modules/kReport2Krona.nf'
include { KT_IMPORT_TEXT			} from './modules/ktImportText.nf'
```

Here we list all the modules that we wish to import into the workflow using the standard `include` syntax.
This consists of two parts: invoking the name of the process inside the curly brackets, and pointing to the relative path where the files are located.
Note that the process names should match exactly how they are written in the module files.

Now, let's write the workflow itself.
We need to declare the primary input for the workflow and invoke the processes on the appropriate inputs.

```groovy title="workflow.nf" linenums="10"
/*
 * workflow
 */

workflow kraken2Flow {
	// required inputs
	take:
		bowtie2_index
		kraken2_db
		reads_ch
	// workflow implementation
	main:
		BOWTIE2(reads_ch, bowtie2_index)
		KRAKEN2(BOWTIE2.out, kraken2_db)
		BRACKEN(KRAKEN2.out, kraken2_db)
		K_REPORT_TO_KRONA(BRACKEN.out)
		KT_IMPORT_TEXT(K_REPORT_TO_KRONA.out)
```

In this declaration you see that we need three primary inputs for the pipeline: the indexed reference genome for Bowtie2, the indexed database for Kraken2 and Bracken, and the paths to the reads; when creating the `nextflow.config`, you will see how to specify these paths.

In the block `main`, you see the names of the processes with one or more parameters inside the parenthesis depicting the exact data flow:

1. `BOWTIE2(reads_ch, bowtie2_index)` will be the first executed process using the reads and the indexed genome; the other processes must wait until this one is finished.
2. Once `BOWTIE2` has completed the task, `KRAKEN2` will take the output, along with the database path to perform the specified task; the remaining process are on hold until Kraken2 is finished.
3. `BRACKEN`, in turn, will take `KRAKEN2` output to run the abundance re-estimation using the same database; `K_REPORT_TO_KRONA` and `KT_IMPORT_TEXT` can not be run until `BRACKEN` is finished with its task.
4. Finally, `K_REPORT_TO_KRONA`, and subsequently `KT_IMPORT_TEXT`, will be run to generate our final goal file which is the Krona plot.

---

## 3. Create the `main.nf` entrypoint workflow

We are getting closer to running the pipeline!
Let's create the `main.nf` file.

Within this file, we create a banner with the pipeline name to be shown when the execution starts using `log.info`.

`log.info` is a method call on a logging object, and can be used to write informational messages to the Nextflow log during a pipeline's execution.
You can learn more about this [here](https://carpentries-incubator.github.io/workflows-nextflow/11-Simple_Rna-Seq_pipeline.html#define-the-pipeline-parameters).

```groovy title="main.nf" linenums="1"
#!/usr/bin/env nextflow

log.info """\
	__________________________________________________________________________________________________________________________________________________
	__________________________________________________________________________________________________________________________________________________
	>=>   >=>                       >=>                                         >=> >=>>=>                                >=>
	>=>  >=>                        >=>                           >=>>=>       >=>  >>   >=>                              >=>
	>=> >=>     >> >==>    >=> >=>  >=>  >=>   >==>    >==>>==>  >>   >=>     >=>   >>    >=> >> >==>    >=> >=>     >==> >=>  >=>   >==>    >==>>==>
	>>=>>        >=>     >=>   >=>  >=> >=>  >>   >=>   >=>  >=>     >=>     >=>    >==>>=>    >=>     >=>   >=>   >=>    >=> >=>  >>   >=>   >=>  >=>
	>=>  >=>     >=>    >=>    >=>  >=>=>    >>===>>=>  >=>  >=>    >=>     >=>     >>    >=>  >=>    >=>    >=>  >=>     >=>=>    >>===>>=>  >=>  >=>
	>=>   >=>    >=>     >=>   >=>  >=> >=>  >>         >=>  >=>  >=>      >=>      >>     >>  >=>     >=>   >=>   >=>    >=> >=>  >>         >=>  >=>
	>=>     >=> >==>      >==>>>==> >=>  >=>  >====>   >==>  >=> >======> >=>       >===>>=>  >==>      >==>>>==>    >==> >=>  >=>  >====>   >==>  >=>
	__________________________________________________________________________________________________________________________________________________
	__________________________________________________________________________________________________________________________________________________
"""
.stripIndent()
```

Then we add an `include` statement to import the `kraken2Flow` workflow from the './workflow.nf' file, as well as a `workflow` block that sets up an input channel and invokes the `kraken2Flow` workflow:

```groovy title="main.nf" linenums="18"
include { kraken2Flow } from './workflow.nf'

workflow {

	reads_ch = Channel .fromFilePairs( params.reads, checkIfExists:true )
	kraken2Flow( params.bowtie2_index, params.kraken2_db, reads_ch )
}
```

1. Creating a channel for the paired-end reads using the channel factory [`fromFilePairs`](https://nextflow.io/docs/latest/reference/channel.html#fromfilepairs).

2. Running the workflow using the reference indexed genome, the Kraken2 database and the channel created for the reads.

---

## 4. Create the `nextflow.config` configuration file

Finally, we create the file `nextflow.config`, where we'll set up the configuration of our pipeline.

This is where we provide default input parameters for the pipeline and enable the use of Docker containers.

```groovy title="nextflow.config" linenums="1"
/*
 * pipeline input parameters
 */

params {
    reads                                 = null
    outdir                                = "/workspaces/training/nf4-science/metagenomics/output"
    bowtie2_index                         = "/workspaces/training/nf4-science/metagenomics/data/yeast/yeast"
    kraken2_db                            = "/workspaces/training/nf4-science/metagenomics/data/viral_db"
}

// Enable using docker as the container engine to run the pipeline
docker.enabled = true
```

!!!tip

    Please note that we are using absolute paths for the parameters. Even though relative paths are preferred for system-independence of the pipelines, often some tools or scripts will fail if parameters are not provided as absolute paths.

---

## 5. Run the pipeline

That's it, we are all set to run the pipeline!

Let's just pick one of the samples provided (you can choose any of them) and run the following command:

```bash
nextflow run main.nf --reads 'data/samples/ERR2143768/ERR2143768_{1,2}.fastq'
```

On the output of the command line, you will see:

```console title="Output"
 N E X T F L O W   ~  version 24.10.4

Launching `main.nf` [determined_lorenz] DSL2 - revision: 8f65b983e6

        __________________________________________________________________________________________________________________________________________________
        __________________________________________________________________________________________________________________________________________________
        >=>   >=>                       >=>                                         >=> >=>>=>                                >=>
        >=>  >=>                        >=>                           >=>>=>       >=>  >>   >=>                              >=>
        >=> >=>     >> >==>    >=> >=>  >=>  >=>   >==>    >==>>==>  >>   >=>     >=>   >>    >=> >> >==>    >=> >=>     >==> >=>  >=>   >==>    >==>>==>
        >>=>>        >=>     >=>   >=>  >=> >=>  >>   >=>   >=>  >=>     >=>     >=>    >==>>=>    >=>     >=>   >=>   >=>    >=> >=>  >>   >=>   >=>  >=>
        >=>  >=>     >=>    >=>    >=>  >=>=>    >>===>>=>  >=>  >=>    >=>     >=>     >>    >=>  >=>    >=>    >=>  >=>     >=>=>    >>===>>=>  >=>  >=>
        >=>   >=>    >=>     >=>   >=>  >=> >=>  >>         >=>  >=>  >=>      >=>      >>     >>  >=>     >=>   >=>   >=>    >=> >=>  >>         >=>  >=>
        >=>     >=> >==>      >==>>>==> >=>  >=>  >====>   >==>  >=> >======> >=>       >===>>=>  >==>      >==>>>==>    >==> >=>  >=>  >====>   >==>  >=>
        __________________________________________________________________________________________________________________________________________________
        __________________________________________________________________________________________________________________________________________________

executor >  local (5)
[fe/4b8409] process > kraken2Flow:BOWTIE2 (ERR2143768)           [100%] 1 of 1 ✔
[14/c3d787] process > kraken2Flow:KRAKEN2 (ERR2143768)           [100%] 1 of 1 ✔
[4c/5d2db3] process > kraken2Flow:BRACKEN (ERR2143768)           [100%] 1 of 1 ✔
[e4/c305af] process > kraken2Flow:K_REPORT_TO_KRONA (ERR2143768) [100%] 1 of 1 ✔
[39/08b32c] process > kraken2Flow:KT_IMPORT_TEXT (ERR2143768)    [100%] 1 of 1 ✔
```

If that worked for you, it's finally time to analyze the results!

You can find all the output files in the `output/` directory by running the `tree` command:

```bash
tree output
```

You should see the following:

```console title="Output contents"
output/
└── ERR2143758
    ├── ERR2143768.b.krona.txt
    ├── ERR2143768.bracken
    ├── ERR2143768.breport
    ├── ERR2143768.k2report
    ├── ERR2143768.kraken2
    ├── ERR2143768.krona.html
    └── ERR2143768.sam
```

Feel free to explore each of the files to understand each process and how data were handled.

The file we are most interested in is the `*.html` file containing the Krona plot.
You can either download it and open it in your browser, or install the [preview](https://marketplace.visualstudio.com/items?itemName=ms-vscode.live-server) extension for Visual Studio.

Keep in mind that we used the viral database, so you will only see information about the viruses contained in the sample.
To perform the analysis with bacteria and beyond you have to download a different database or build your own.

---

### Takeaway

Congratulations! You just have performed a taxonomic annotation and abundance estimation of your metagenomics sample.

### What's next?

Learn how to modify the workflow to process multiple samples in parallel, create a report with different diversity measurements, execute customized scripts, and use Nextflow operators and conditionals to control the workflow.
