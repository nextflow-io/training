# Part 2: Process parallelization, operators, conditional execution and custom scripts

In this part, we are going to rely on same pipeline structure we built in Part 1 to extend for: 

1. Multi-sample analysis
2. Use of a Nextflow operator 
3. Control the execution of the workflow according to the input
4. Include a process that runs a customized script

---

## 1. Multi-sample input

With our shining brand-new pipeline, we are at this moment able to analyze each sample individually by running the workflow multiple times. Nonetheless, one of the most powerful capabilities by Nextflow is its native parallel execution according to the available resources the executor finds. You can think of this as a sort of "integrated _for_ loop" that will process all the samples in parallel in a single run without the need of re-running the pipeline.

To achieve this purpose, there are two possibilities:

* The use of wildcards in the input (this can be tricky and requires to take into account particular folder structures).
* Create a file that points out to the sample files regardless of their location in the file system. 

In this course, we will target the second input option, albeit you are welcome to explore how you can use the first option by checking the [Nextflow documentation](https://www.nextflow.io/docs/latest/working-with-files.html).

To move forward, let's create then the file `samplesheet.csv` inside the folder **data**:

```csv title="data/samplesheet.csv" linenums="1"
sample_id,fastq_1,fastq_2
ERR2143768,/workspaces/training/nf4-science/metagenomics/data/samples/ERR2143768/ERR2143768_1.fastq,/workspaces/training/nf4-science/metagenomics/data/samples/ERR2143768/ERR2143768_2.fastq
ERR2143769,/workspaces/training/nf4-science/metagenomics/data/samples/ERR2143769/ERR2143769_1.fastq,/workspaces/training/nf4-science/metagenomics/data/samples/ERR2143769/ERR2143769_2.fastq
ERR2143770,/workspaces/training/nf4-science/metagenomics/data/samples/ERR2143770/ERR2143770_1.fastq,/workspaces/training/nf4-science/metagenomics/data/samples/ERR2143770/ERR2143770_2.fastq
ERR2143774,/workspaces/training/nf4-science/metagenomics/data/samples/ERR2143774/ERR2143774_1.fastq,/workspaces/training/nf4-science/metagenomics/data/samples/ERR2143774/ERR2143774_2.fastq
```

Here, we have provided the `sample id` and the absolute paths to both forward and reverse reads per sample. Please notice that the files are not required to be stored in the directory; however, it is recommend to maintain a consistent folder structure.

Now, we can not use this file as input in the current state of the pipeline given that it expects only a path to create a paire-end channel. Let's include then an additional parameter in the `nextflow.config` file (notice that it would go inside the parameter block, keeping the same structure):

```groovy title="nextflow.config" linenums="10"
    sheet_csv                             = null
```

We initialize this parameter as `null` since it can be used or not. Now, we need to modify the `main.nf` file to state how the input should be handled depending of the type of input:

```groovy title="main.nf" linenums="22"
	    if(params.reads){
		    reads_ch = Channel .fromFilePairs( params.reads, checkIfExists:true )
		} else {
		    reads_ch = Channel.fromPath( params.sheet_csv )
				   .splitCsv(header:true)
				   .map { row-> tuple(row.sample_id, [file(row.fastq_1), file(row.fastq_2)]) }
		}
```

This modified declaration states that if we use the parameter `--reads` when we invoke the `nextflow run main.nf`, the _reads_ channel will be created using only the path to paired-end files. Otherwise, we must include the parameter `--sheet_csv` with the corresponding file containing the sample information. Being so, it is necessary to use one of the two forms of input; if we use both at the same time, the `--reads` will predominate or if none of them is indicated, the pipeline will fail. Do not worry now for the way in which channel is created using the `.csv` file, this declaration is quite stantard and you can just copy and paste for other pipelines in which you would like to use it; however, you can learn more about this [here](https://nextflow-io.github.io/patterns/process-per-csv-record/).

Now, we would be ready to re-run the pipeline to process all the samples in a single call. Notwithstanding, the inclusion of additional samples has the advantage that we can expand the analysis to estimate β-diversity and compare them to extract important insights. 

---

# 2. Additional processes

## 2.1 Kraken-biom

Let's create a new module that is going to handle the Bracken output to produce a Biological Observation Matrix (BIOM) file that concatenates the species abundance in each sample. The `kraken_biom.nf` file will be located in the **modules** directory:

```groovy title="modules/kraken_biom.nf" linenums="1"
process KRAKEN_BIOM {
	  tag "merge_samples"
    publishDir "$params.outdir", mode:'copy'
    container "community.wave.seqera.io/library/kraken-biom:1.2.0--f040ab91c9691136"

    input:
    val "files"

    output:
    path "merged.biom"

    script:
    """
    list=(${files.join(' ')})
    extracted=\$(echo "\${list[@]}" | tr ' ' '\n' | awk 'NR % 3 == 2')
    kraken-biom \${extracted} --fmt json -o merged.biom
    """
}
```

This process will _collect_ each output from the Bracken files to build a single `*.biom` file that contains the abundance species data of all the samples. In the `script` statement we find three tasks to execute, the first two lines are for variable manipulation required to handle the type of input this process receives (more about this when modifying `workflow.nf` below), and the second line executes the kraken-biom command that is available thanks to specified container.

### 2.1.1 Operator _collect()_ and conditional execution

Nextflow provides a high number of operators that smooth data handling and orchestrates the workflow to do exactly what we want. In this case, the process `KRAKEN_BIOM` requires all the files produced by Bracken belonging to each sample, which means that `KRAKEN_BIOM` can not be triggered until all Bracken processes are finished. For this task, the operator _collect()_ comes really handy, and therefore let's include it in our `workflow.nf`... but wait! Let's recall that `KRAKEN_BIOM` and the following `KNIT_PHYLOSEQ` are only triggered if the execution is aiming at processing more than one sample. Being so, we will include these processes and modify the workflow execution to add the conditional statement in the `workflow.nf`:

```groovy title="workflow.nf" linenums="9"
include { KRAKEN_BIOM               }   from './modules/kraken_biom.nf'
```

```groovy title="workflow.nf" linenums="29"
        if(params.sheet_csv){
		    KRAKEN_BIOM(BRACKEN.out.collect())
		}
```

Here, you can see that we have added the operator _collect()_ to capture all the output files from `BRACKEN`, and this is happening only if we are using as input `--sheet_csv`. This operator is going to return a list of the elements specified in the output of the process (`BRACKEN`), and, for instance, we are interested in each "second" (indices 1,4,7...) element of the list to run the _kraken-biom_ command; this is the reason why within the `script` statement in `kraken_biom.nf` we have incluced two codelines to obtain the paths to these files. If this is not entirely clear, please check the [Nextflow documentation](https://www.nextflow.io/docs/latest/reference/operator.html#collect).    

## 2.2 Phyloseq

### 2.2.1 Including a customized script

We are at the last step of the pipeline execution, and now we need to process the `*.biom` file by transforming it into a Phyloseq object, which is easier to use, more intuitive to understand, and is equipped with multiple tools and methods to plot. Another amazing feature by Nextflow is the possibility to run the so-called _Scripts à la carte_, which means that a process does not necessarily requires an external tool to execute, and hence you can develop your own analysis with customized scripts, i.e., R or Python. Here, we will run an R script inside the module `knit_phyloseq.nf` to create and process the Phyloseq object taking as input the ouput from `kraken_biom.nf`:

```groovy title="modules/kraken_biom.nf" linenums="1"
process KNIT_PHYLOSEQ {
	tag "knit_phyloseq"
    publishDir "$params.outdir", mode:'copy'
    container "community.wave.seqera.io/library/bioconductor-phyloseq_knit_r-base_r-ggplot2_r-rmdformats:6efceb52eb05eb44"

    input:
    path "merged"

    output:
    stdout

    script:
    def report = params.report
    def outdir = params.outdir
    """
    biom_path=\$(realpath ${merged})
    outreport=\$(realpath ${outdir})
    Rscript -e "rmarkdown::render('${report}', params=list(args='\${biom_path}'),output_file='\${outreport}/report.html')"
    """
}
```

As you can see, we are declaring some variables both in Nextflow and bash to able to call the script. This is a special case since this type of scripts can be stored in the **bin** directory for Nextflow to find them directly. Nevertheless, as we are not "running the script" directly but we are calling `Rscript` to render a final `*.html` report, Nextflow is not able to automatically find the customized script nor detect when report is rendered. As a result the ouput from this process is just a standard/command-line ouput, and we have to include an additional parameter in the `nextflow.config` file:

```groovy title="nextflow.config" linenums="11"
    report                             = "/workspaces/training/nf4-science/metagenomics/bin/report.Rmd"
```

In addition, please notice the `container` used for the `KNIT_PHYLOSEQ`, which is combination of multiple packages required to render the `*.html` report. This is possible thanks to an awesome tool called [Seqera Containers](https://seqera.io/containers/), which is able to build almost any container (for docker or singularity!) by just "merging" different PyPI or Conda packages; please give it a try and be amazed by Seqera Containers.

Also, we have to include this new process within `workflow.nf`:

```groovy title="workflow.nf" linenums="10"
include { KNIT_PHYLOSEQ             }   from './modules/knit_phyloseq.nf'
```

We need to call it as well inside the conditional execution if multi-sample is being handled:

```groovy title="workflow.nf" linenums="31"
        KNIT_PHYLOSEQ(KRAKEN_BIOM.out)
```

---

## 3. Execution

Now, we are completely set to run the analysis for as many samples as we would like, and we will obtain a final report depicting different metrics regarding taxonomic abundance, network analysis, and α and β-diversity. Let's execute:

```bash
nextflow run main.nf --sheet_csv 'data/samplesheet.csv'
```

On the output of the command line, you will see:

```console title="Output"
 N E X T F L O W   ~  version 24.10.4

Launching `main.nf` [stoic_miescher] DSL2 - revision: 8f65b983e6

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

executor >  local (22)
[4e/914152] kraken2Flow:BOWTIE2 (ERR2143774)           [100%] 4 of 4 ✔
[bf/7fcac7] kraken2Flow:KRAKEN2 (ERR2143774)           [100%] 4 of 4 ✔
[f5/aa12aa] kraken2Flow:BRACKEN (ERR2143774)           [100%] 4 of 4 ✔
[e9/84eb9d] kraken2Flow:K_REPORT_TO_KRONA (ERR2143774) [100%] 4 of 4 ✔
[59/456551] kraken2Flow:KT_IMPORT_TEXT (ERR2143768)    [100%] 4 of 4 ✔
[da/7b9f45] kraken2Flow:KRAKEN_BIOM (merge_samples)    [100%] 1 of 1 ✔
[d0/deccc9] kraken2Flow:KNIT_PHYLOSEQ (knit_phyloseq)  [100%] 1 of 1 ✔
```

Keep in mind that since the execution is in parallel, the order in which the samples are processed is random and the order in which `sample ids` appear will differ among executions. Also, during while the pipeline is running you will see that `KRAKEN_BIOM`, and hence `KNIT_PHYLOSEQ`, will not be triggered until all the samples are processed by the previous processes.

Finally, inside the **output** directory, you will see multiple folders with the exact `sample ids`, and within these all the output files, including the files to visualize the Krona plots. Likewise, in the **output** folder you will see the file `report.html` which is ready to be opened and explored. It's your time to analyze it!

---

### Takeaway

You just learnt how control workflow execution by including conditionals and operators, process multiple samples simultaneously and running a customized script to perform a metagenomics data analysis at read level. 

### What's next?

Great! You are well equipped now to start developing your ownn pipelines.
