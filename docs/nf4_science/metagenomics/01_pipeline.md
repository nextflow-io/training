# Part 1: Pipeline development and single sample run

In the field of metagenomics data analysis, there is an endless universe of pipelines or methodologies you can follow to explore and characterize your samples. 
For this course, we propose to wrap with Nextflow the methodology published by [Jennifer Lu et al. (2022)](https://www.nature.com/articles/s41596-022-00738-y). 

## 1. Workflow

The workflow is designed as follows:

<p align="center">
    <img src="src/workflow_kraken.png" alt="MAGFlow" width="100%">
</p>

As you can see from the picture, the pipeline will undergo as follows:

1. The input is **FASTQ** files from one or multiple samples. For this course, we will be using only paired-end reads recovered from an oligotrophic, phosphorus-deficient pond in Cuatro Ciénegas, Mexico [(Okie et al.,2020)](https://elifesciences.org/articles/49816); the BioProject accesion number is [PRJEB22811](https://www.ncbi.nlm.nih.gov/bioproject/PRJEB22811).
2. **Host removal** with [**Bowtie2**](https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml) by aligning the reads against an indexed reference genome. In this case we are using the index of _O. zativa_ given the storage limitations we have in this codespace, although you can use any organism of your interest by building [your own index](https://www.metagenomics.wiki/tools/bowtie2/index) or downloading a [precomputed](https://benlangmead.github.io/aws-indexes/bowtie) one.
3. **Taxonomic classification** with [**Kraken2**](https://ccb.jhu.edu/software/kraken2/). This tools relies on a indexed database that can be [downloaded](https://benlangmead.github.io/aws-indexes/k2) or you can build your customized version following specific [instructions](https://avilpage.com/2024/07/mastering-kraken2-build-custom-db.html). We will be using the Viral database given the storage limitations (I know, it's a bit annoying), therefore this methodology is pointing at "viral metagenomics"; however, by just switching to any other database you can analyze your samples to annotate bacteria, archaea and more.
4. **Bayesian species abundance re-estimation** with [**Bracken**](https://ccb.jhu.edu/software/bracken/index.shtml?t=manual). This software is designed to compute species abundance using Kraken classification results as stated in the reference paper, and it also uses some files contained in the dabatase folders such as the kmer distribution files (don't worry about this now, you can learn about how the method works afterwards).
5. [**Krona**](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-12-385) plots are generated using the Bracken output to visualize interactively the relative abundance of each annotated species.
6. If multiple samples are used as input to the pipeline, the bracken reports will be concatenated and converted into a Biological Observation Matrix [(**BIOM**)](https://biom-format.org/) file.
7. Finally, the BIOM file will be first converted in a [**Phyloseq**](https://joey711.github.io/phyloseq/index.html) object, and this object will be further processed to generate absolute plots, estimate both α and β-diversity and perform a network analysis; this information will be presented in a final `report.html`. If you can to follow the code to generate the plots and metrics, you can check this Phyloseq [tutorial](https://vaulot.github.io/tutorials/Phyloseq_tutorial.html). 

!!!tip

    If you feel a bit overwhelmed by the theoretical background of the methodology, we strongly encourage you to check this [Carpentries](https://carpentries-lab.github.io/metagenomics-analysis/) lesson first, where each step will be explained step by step using interesting examples.

## 2. Modules

Once we have a clear overview of what we want to achieve, we can start developing the modules that are going to perform each task; you can think of them as building blocks that we can stack up to construct a big tower (this is valid for this pipeline as it is linear).

### 2.1 Bowtie2

As previously, the objective with this module is to align the reads against a reference genome. Let's create then the `bowtie2.nf` file inside the **modules** folder to write the following code:

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
    export BOWTIE2_INDEXES=/workspaces/training/nf4-science/metagenomics/data/oryza
    bowtie2 -x $bowtie2_index -1 ${reads[0]} -2 ${reads[1]} -p 2 -S ${sample_id}.sam --un-conc-gz ${sample_id}
 	  """
}
```
Let's take a moment to break down what we are seeing here:
- The process name is `BOWTIE2`, this is important when creating the workflow file.
- The `tag` directive is used to indicate which sample is being processed at a determined moment. This will be useful when running the pipeline.
- `publishDir` points out to the directory where the ouput is stored. In this case we are taking the path from the parameters, and within it subfolders with the sample namse will be created to store each _.sam_ file, creating a copy of such files.
- `container` indicates the docker container on which the process will be run. More information about this can be found in the part 1 of the [RNASeq course](../rnaseq). 
