For each processed sample the pipeline stores results into a folder named after the sample identifier. These folders are created in the directory specified as a parameter in `params.results`.

Result files for this workshop can be found in the folder `results` within the current folder. There you should see a directory called `ENCSR000COQ/` containing the following files:

Variant calls  
`final.vcf`

This file contains all somatic variants (SNVs) called from RNAseq data. You will see variants that pass all filters, with the `PASS` keyword in the <span class="red">7th</span> field of the vcf file (`filter status`), and also those that did not pass one or more filters.

`commonSNPs.diff.sites_in_files`

Tab-separated file with comparison between variants obtained from RNAseq and "known" variants from DNA.

The file is sorted by genomic position and contains 8 fields:

<table>
<colgroup>
<col style="width: 33%" />
<col style="width: 33%" />
<col style="width: 33%" />
</colgroup>
<tbody>
<tr class="odd">
<td style="text-align: left;"><p><span class="red">1</span></p></td>
<td style="text-align: center;"><pre><code> CHROM   </code></pre></td>
<td style="text-align: left;"><p>chromosome name;</p></td>
</tr>
<tr class="even">
<td style="text-align: left;"><p><span class="red">2</span></p></td>
<td style="text-align: center;"><pre><code> POS1    </code></pre></td>
<td style="text-align: left;"><p>position of the SNV in file #1 (RNAseq data);</p></td>
</tr>
<tr class="odd">
<td style="text-align: left;"><p><span class="red">3</span></p></td>
<td style="text-align: center;"><pre><code> POS2    </code></pre></td>
<td style="text-align: left;"><p>position of SNV in file #2 (DNA "known" variants);</p></td>
</tr>
<tr class="even">
<td style="text-align: left;"><p><span class="red">4</span></p></td>
<td style="text-align: center;"><pre><code> IN_FILE </code></pre></td>
<td style="text-align: left;"><p>flag whether SNV is present in the file #1 <em>1</em>, in the file #2 <em>2</em>, or in both files <em>B</em>;</p></td>
</tr>
<tr class="odd">
<td style="text-align: left;"><p><span class="red">5</span></p></td>
<td style="text-align: center;"><pre><code> REF1    </code></pre></td>
<td style="text-align: left;"><p>reference sequence in the file 1;</p></td>
</tr>
<tr class="even">
<td style="text-align: left;"><p><span class="red">6</span></p></td>
<td style="text-align: center;"><pre><code> REF2    </code></pre></td>
<td style="text-align: left;"><p>reference sequence in the file 2;</p></td>
</tr>
<tr class="odd">
<td style="text-align: left;"><p><span class="red">7</span></p></td>
<td style="text-align: center;"><pre><code> ALT1    </code></pre></td>
<td style="text-align: left;"><p>alternative sequence in the file 1;</p></td>
</tr>
<tr class="even">
<td style="text-align: left;"><p><span class="red">8</span></p></td>
<td style="text-align: center;"><pre><code> ALT2    </code></pre></td>
<td style="text-align: left;"><p>alternative sequence in the file 2</p></td>
</tr>
</tbody>
</table>

`known_snps.vcf`

Variants that are common to RNAseq and "known" variants from DNA.

Allele specific expression quantification  
`ASE.tsv`

Tab-separated file with allele counts at common SNVs positions (only SNVs from the file `known_snps.vcf`)

The file is sorted by coordinates and contains 13 fields:

<table>
<colgroup>
<col style="width: 33%" />
<col style="width: 33%" />
<col style="width: 33%" />
</colgroup>
<tbody>
<tr class="odd">
<td style="text-align: left;"><p><span class="red">1</span></p></td>
<td style="text-align: center;"><pre><code> contig        </code></pre></td>
<td style="text-align: left;"><p>contig, scaffold or chromosome name of the variant</p></td>
</tr>
<tr class="even">
<td style="text-align: left;"><p><span class="red">2</span></p></td>
<td style="text-align: center;"><pre><code> position      </code></pre></td>
<td style="text-align: left;"><p>position of the variant</p></td>
</tr>
<tr class="odd">
<td style="text-align: left;"><p><span class="red">3</span></p></td>
<td style="text-align: center;"><pre><code> variant ID    </code></pre></td>
<td style="text-align: left;"><p>variant ID in the dbSNP</p></td>
</tr>
<tr class="even">
<td style="text-align: left;"><p><span class="red">4</span></p></td>
<td style="text-align: center;"><pre><code> refAllele     </code></pre></td>
<td style="text-align: left;"><p>reference allele sequence</p></td>
</tr>
<tr class="odd">
<td style="text-align: left;"><p><span class="red">5</span></p></td>
<td style="text-align: center;"><pre><code> altAllele     </code></pre></td>
<td style="text-align: left;"><p>alternate allele sequence</p></td>
</tr>
<tr class="even">
<td style="text-align: left;"><p><span class="red">6</span></p></td>
<td style="text-align: center;"><pre><code> refCount      </code></pre></td>
<td style="text-align: left;"><p>number of reads that support the reference allele</p></td>
</tr>
<tr class="odd">
<td style="text-align: left;"><p><span class="red">7</span></p></td>
<td style="text-align: center;"><pre><code> altCount      </code></pre></td>
<td style="text-align: left;"><p>number of reads that support the alternate allele</p></td>
</tr>
<tr class="even">
<td style="text-align: left;"><p><span class="red">8</span></p></td>
<td style="text-align: center;"><pre><code> totalCount    </code></pre></td>
<td style="text-align: left;"><p>total number of reads at the site that support both reference and alternate allele and any other alleles present at the site</p></td>
</tr>
<tr class="odd">
<td style="text-align: left;"><p><span class="red">9</span></p></td>
<td style="text-align: center;"><pre><code> lowMAPQDepth  </code></pre></td>
<td style="text-align: left;"><p>number of reads that have low mapping quality</p></td>
</tr>
<tr class="even">
<td style="text-align: left;"><p><span class="red">10</span></p></td>
<td style="text-align: center;"><pre><code> lowBaseQDepth </code></pre></td>
<td style="text-align: left;"><p>number of reads that have low base quality</p></td>
</tr>
<tr class="odd">
<td style="text-align: left;"><p><span class="red">11</span></p></td>
<td style="text-align: center;"><pre><code> rawDepth      </code></pre></td>
<td style="text-align: left;"><p>total number of reads at the site that support both reference and alternate allele and any other alleles present at the site</p></td>
</tr>
<tr class="even">
<td style="text-align: left;"><p><span class="red">12</span></p></td>
<td style="text-align: center;"><pre><code> otherBases    </code></pre></td>
<td style="text-align: left;"><p>number of reads that support bases other than reference and alternate bases</p></td>
</tr>
<tr class="odd">
<td style="text-align: left;"><p><span class="red">13</span></p></td>
<td style="text-align: center;"><pre><code> improperPairs </code></pre></td>
<td style="text-align: left;"><p>number of reads that have malformed pairs</p></td>
</tr>
</tbody>
</table>

Allele frequency histogram  
`AF.histogram.pdf`

This file contains a histogram plot of allele frequency for SNVs common to RNA-seq and "known" variants from DNA.
