# Part 2: Process parallelization, operator and conditional execution

In this part, we are going to rely on same pipeline structure we built in Part 1 to extend for: 

1. Multi-sample analysis
2. Use of a Nextflow operator 
3. Control the execution of the workflow according to the input

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
