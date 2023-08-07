# Data Description

The input data used to test the pipeline implementation is described below. For the purpose of this project, only a subset of the original data is used for most of the data types.

## Genome assembly

`genome.fa`

The human genome assembly hg19 (GRCh37) from [GenBank](https://www.ncbi.nlm.nih.gov/assembly/GCA_000001405.1), chromosome 22 only.

## RNA-seq reads

`ENCSR000COQ[12]_[12].fastq.gz`

The RNA-seq data comes from the human GM12878 cell line from whole cell, cytosol and nucleus extraction (see table below).

The libraries are stranded PE76 Illumina GAIIx RNA-Seq from rRNA-depleted Poly-A+ long RNA (`> 200` nucleotides in size).

Only reads mapped to the [22q11^](http://genome-euro.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr22%3A14700001-25900000&hgsid=221945779_QucOFSFGagd1cn9uVki0TFjrxSBU) locus of the human genome (`chr22:16000000-18000000`) are used.

| ENCODE ID                                                             | Cellular fraction | Replicate ID | File names                                             |                                                        |
| --------------------------------------------------------------------- | ----------------- | ------------ | ------------------------------------------------------ | ------------------------------------------------------ |
| [ENCSR000COQ](https://www.encodeproject.org/experiments/ENCSR000COQ/) | Whole Cell        | 1<br>2       | `ENCSR000COQ1_1.fastq.gz`<br>`ENCSR000COQ2_1.fastq.gz` | `ENCSR000COQ1_2.fastq.gz`<br>`ENCSR000COQ2_2.fastq.gz` |
| [ENCSR000CPO](https://www.encodeproject.org/experiments/ENCSR000CPO/) | Nuclear           | 1<br>2       | `ENCSR000CPO1_1.fastq.gz`<br>`ENCSR000CPO2_1.fastq.gz` | `ENCSR000CPO1_2.fastq.gz`<br>`ENCSR000CPO2_2.fastq.gz` |
| [ENCSR000COR](https://www.encodeproject.org/experiments/ENCSR000COR/) | Cytosolic         | 1<br>2       | `ENCSR000COR1_1.fastq.gz`<br>`ENCSR000COR2_1.fastq.gz` | `ENCSR000COR1_2.fastq.gz`<br>`ENCSR000COR2_2.fastq.gz` |

## "Known" variants

`known_variants.vcf.gz`

Known variants come from high confident variant calls for GM12878 from the [Illumina Platinum Genomes](https://www.illumina.com/platinumgenomes.html) project. These variant calls were obtained by taking into account pedigree information and the concordance of calls across different methods.

We’re using the subset from chromosome 22 only.

## Blacklisted regions

`blacklist.bed`

Blacklisted regions are regions of the genomes with anomalous coverage. We use regions for the hg19 assembly, taken from the [ENCODE project portal](https://www.encodeproject.org/annotations/ENCSR636HFF/). These regions were identified with DNAse and ChiP-seq samples over ~60 human tissues/cell types, and had a very high ratio of multi-mapping to unique-mapping reads and high variance in mappability.
