# Part 1: Method overview

In the field of metagenomics data analysis, there is an endless universe of pipelines or methodologies you can follow to explore and characterize your samples.
We recommend this comprehensive [review](https://www.sciencedirect.com/science/article/pii/S2001037021004931) for you to explore the different existing approaches.
For this course, we propose to wrap with Nextflow the protocol published by [Jennifer Lu et al. (2022)](https://www.nature.com/articles/s41596-022-00738-y).

The example dataset we will use to demonstrate the analysis consists of only paired-end reads recovered from an oligotrophic, phosphorus-deficient pond in Cuatro Ciénegas, Mexico ([Okie et al.,2020](https://elifesciences.org/articles/49816)) in FASTQ format.
The BioProject accession number is [PRJEB22811](https://www.ncbi.nlm.nih.gov/bioproject/PRJEB22811).

---

## 1. Workflow design

Our goal is to develop a workflow that takes **FASTQ** files from one or multiple samples as input and applies the following processing steps: host removal, taxonomic classification, Bayesian re-estimation of species abundance, and generation of plots and metrics.

<div markdown class="metagenomics">

![Metagenomics](../../assets/img/workflow_kraken.png)

</div>

To perform these steps, we will use the following tools:

1. **Host removal** with [**Bowtie2**](https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml) by aligning the reads against an indexed reference genome.
   Here, we are using the indexed genome of yeast, but you can use any organism you are interested in by building [your own index](https://www.metagenomics.wiki/tools/bowtie2/index) or downloading a [precomputed one](https://benlangmead.github.io/aws-indexes/bowtie).
2. **Taxonomic classification** with [**Kraken2**](https://ccb.jhu.edu/software/kraken2/).
   This tools relies on a indexed database that can be [downloaded](https://benlangmead.github.io/aws-indexes/k2).
   Alternatively, you can build your customized version following [these instructions](https://avilpage.com/2024/07/mastering-kraken2-build-custom-db.html).
   Here, we will use the Viral database, therefore this methodology is labeled as "viral metagenomics".
   However, you can annotate bacteria, archaea and more simply by switching to another database.
3. **Bayesian re-estimation of species abundance** with [**Bracken**](https://ccb.jhu.edu/software/bracken/index.shtml?t=manual).
   This software is designed to compute species abundance using Kraken classification results as described in the reference paper.
   It also uses some files contained in the dabatase folders such as the kmer distribution files.
   This is a fairly complex analysis, but you don't need to know the details in order to follow this tutorial; you can learn about how the method works afterwards.
4. **Plot generation** with [**Krona**](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-12-385) from the Bracken output.
   This will allow us to visualize interactively the relative abundance of each annotated species.
5. (Multi-sample) **Concatenation with kraken-biom.**
   If multiple samples are provided, the Bracken reports will be concatenated and converted into a [Biological Observation Matrix (BIOM)](https://biom-format.org/) file.
6. (Multi-sample) **Generation of final report with Phyloseq**
   The BIOM file will be converted to a [Phyloseq](https://joey711.github.io/phyloseq/index.html) object, and this object will be further processed to generate absolute plots, estimate both α and β-diversity and perform a network analysis.
   This information will be presented in a final `report.html`.
   To learn more about the code used to generate the plots and metrics, check out this Phyloseq [tutorial](https://vaulot.github.io/tutorials/Phyloseq_tutorial.html).

!!!tip

    If you feel a bit overwhelmed by the theoretical background of the methodology, we strongly encourage you to check this [Carpentries](https://carpentries-lab.github.io/metagenomics-analysis/) lesson first, where the concepts are explained step by step using interesting examples.

---

## 2. [TODO: add optional manual testing of the various tools via containers]

---

### Takeaway

You understand the underlying method and the overall design of the workflow.

**[TODO: You have tested all the individual commands interactively in the relevant containers.]
**

### What's next?

Learn how to wrap those same commands into a multi-step workflow that uses containers to execute the work.
