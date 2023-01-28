The input data used to test the pipeline implementation is described below. For the purpose of this project, only a subset of the original data is used for most of the data types.

Genome assembly  
`genome.fa`

The human genome assembly <span class="crg">hg19 (GRCh37)</span> from [GenBank](https://www.ncbi.nlm.nih.gov/assembly/GCA_000001405.1), chromosome 22 only.

RNA-seq reads  
`ENCSR000COQ[12]_[12].fastq.gz`

The RNA-seq data comes from the human <span class="crg">GM12878</span> cell line from whole cell, cytosol and nucleous extraction (see table below).

The libraries are <span class="crg">stranded PE76 Illumina GAIIx</span> RNA-Seq from <span class="crg">rRNA-depleted Poly-A+</span> long RNA (`> 200` nucleotides in size).

Only reads mapped to the [22q11^](http://genome-euro.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr22%3A14700001-25900000&hgsid=221945779_QucOFSFGagd1cn9uVki0TFjrxSBU) locus of the human genome (`chr22:16000000-18000000`) are used.

<table>
<colgroup>
<col style="width: 25%" />
<col style="width: 25%" />
<col style="width: 25%" />
<col style="width: 25%" />
</colgroup>
<tbody>
<tr class="odd">
<td style="text-align: left;"><p>ENCODE ID</p></td>
<td style="text-align: left;"><p>Cellular fraction</p></td>
<td style="text-align: left;"><p>replicate ID</p></td>
<td style="text-align: left;"><p>file names</p></td>
</tr>
<tr class="even">
<td style="text-align: left;"><p><a href="https://www.encodeproject.org/experiments/ENCSR000COQ/">ENCSR000COQ</a></p></td>
<td style="text-align: left;"><p>Whole Cell</p></td>
<td style="text-align: left;"><p>1</p></td>
<td style="text-align: left;"><pre><code>ENCSR000COQ1_1.fastq.gz
ENCSR000COQ1_2.fastq.gz</code></pre></td>
</tr>
<tr class="odd">
<td style="text-align: left;"><p>2</p></td>
<td style="text-align: left;"><pre><code>ENCSR000COQ2_1.fastq.gz
ENCSR000COQ2_2.fastq.gz</code></pre></td>
<td></td>
<td></td>
</tr>
<tr class="even">
<td style="text-align: left;"><p><a href="https://www.encodeproject.org/experiments/ENCSR000CPO/">ENCSR000CPO</a></p></td>
<td style="text-align: left;"><p>Nuclear</p></td>
<td style="text-align: left;"><p>1</p></td>
<td style="text-align: left;"><pre><code>ENCSR000CPO1_1.fastq.gz
ENCSR000CPO1_2.fastq.gz</code></pre></td>
</tr>
<tr class="odd">
<td style="text-align: left;"><p>2</p></td>
<td style="text-align: left;"><pre><code>ENCSR000CPO2_1.fastq.gz
ENCSR000CPO2_2.fastq.gz</code></pre></td>
<td></td>
<td></td>
</tr>
<tr class="even">
<td style="text-align: left;"><p><a href="https://www.encodeproject.org/experiments/ENCSR000COR/">ENCSR000COR</a></p></td>
<td style="text-align: left;"><p>Cytosolic</p></td>
<td style="text-align: left;"><p>1</p></td>
<td style="text-align: left;"><pre><code>ENCSR000COR1_1.fastq.gz
ENCSR000COR1_1.fastq.gz</code></pre></td>
</tr>
<tr class="odd">
<td style="text-align: left;"><p>2</p></td>
<td style="text-align: left;"><pre><code>ENCSR000COR2_1.fastq.gz
ENCSR000COR2_1.fastq.gz</code></pre></td>
<td></td>
<td></td>
</tr>
</tbody>
</table>

"Known" variants  
`known_variants.vcf.gz`

Known variants come from high confident variant calls for <span class="crg">GM12878</span> from the [Illumina Platinum Genomes](https://www.illumina.com/platinumgenomes.html) project. These variant calls were obtained by taking into account pedigree information and the concordance of calls across different methods.

We’re using the subset from chromosome 22 only.

Blacklisted regions  
`blacklist.bed`

Blacklisted regions are regions of the genomes with anomalous coverage. We use regions for the <span class="crg">hg19</span> assembly, taken from the [ENCODE project portal](https://www.encodeproject.org/annotations/ENCSR636HFF/). These regions were identified with DNAse and ChiP-seq samples over ~60 human tissues/cell types, and had a very high ratio of multi-mapping to unique-mapping reads and high variance in mappability.
