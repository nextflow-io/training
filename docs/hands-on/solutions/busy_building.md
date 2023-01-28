<div class="formalpara-title">

**Solution**

</div>

``` nextflow
/*
 * Process 1A: Create a FASTA genome index with samtools
 */

process '1A_prepare_genome_samtools' {

  input:
    path genome from params.genome 

  output:
    path "${genome}.fai" into genome_index_ch

  script:
  """
  samtools faidx ${genome}
  """
}
```

- The solution is to use the **`params.genome`** parameter defined at the beginning of the script.
