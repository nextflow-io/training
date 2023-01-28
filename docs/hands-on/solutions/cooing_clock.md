<div class="formalpara-title">

**Solution**

</div>

``` nextflow
/*
 * Process 1B: Create a FASTA genome sequence dictionary with Picard for GATK
 */

process '1B_prepare_genome_picard' {

  input:
    path genome from params.genome 

  output:
    path "${genome.baseName}.dict" into genome_dict_ch 

  script:
  """
  PICARD=`which picard.jar`
  java -jar \$PICARD CreateSequenceDictionary R= $genome O= ${genome.baseName}.dict
  """
}
```

- Take as input the `genome` file from the `params.genome` parameter

- Give as output the file `${genome.baseName}.dict` and adds it to the channel `genome_dict_ch`
