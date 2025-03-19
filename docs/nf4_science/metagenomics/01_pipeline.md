# Part 1: Pipeline development and single sample run

In the field of metagenomics data analysis, there is an endless universe of pipelines or methodologies you can follow to explore and characterize your samples. 
For this course, we propose to wrap with Nextflow the methodology published by [Jennifer Lu et al. (2022)](https://www.nature.com/articles/s41596-022-00738-y). 
The workflow is designed as follows:

<p align="center">
    <img src="src/workflow_kraken.png" alt="MAGFlow" width="100%">
</p>

As you can see from the picture, the pipeline will undergo as follows:
- Host removal with [**Bowtie2**](https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml) by aligning the reads against an indexed reference genome. In this case we are using the index of _O. zativa_ given the storage limitations we have in this codespace, although you can use any organism of your interest by building your own index or downloading a [precomputed](https://benlangmead.github.io/aws-indexes/bowtie) one.
- Taxonomic classification with [**Kraken2**](https://ccb.jhu.edu/software/kraken2/). This tools relies on a indexed database that can be [downloaded](https://benlangmead.github.io/aws-indexes/k2) or you can build your customized version following specific [instructions](https://avilpage.com/2024/07/mastering-kraken2-build-custom-db.html). We will be using the Viral database given the storage limitations (I know, it's a bit annoying), therefore this methodology is pointing at "viral metagenomics"; however, by just switching to any other database you can analyze your samples to annotate bacteria, archaea and more.
- Bayesian species abundance re-estimation with [**Bracken**](https://ccb.jhu.edu/software/bracken/index.shtml?t=manual). This software is designed to compute species abundance using Kraken classification results as stated in the reference paper, and it also uses some files contained in the dabatase folders such as the kmer distribution files (don't worry about this now, you can learn about how the method works afterwards).
- [**Krona**](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-12-385) plots are generated using the Bracken output to visualize interactively the relative abundance of each annotated species.
- If multiple samples are used as input to the pipeline, the bracken reports will be concatenated and converted into a Biological Observation Matrix [(**BIOM**)](https://biom-format.org/) file.
- Finally, the BIOM file will be first converted in a [**Phyloseq**](https://joey711.github.io/phyloseq/index.html) object, and this object will be further processed to generate absolute plots, estimate both α and β-diversity and perform a network analysis; this information will be presented in a final `report.html`.

!!!tip

    If you feel a bit overwhelmed by the theoretical background of the methodology, we strongly encourage you to check this [Carpentries](https://carpentries-lab.github.io/metagenomics-analysis/) lesson first, where each step will be explained step by step using interesting examples.
