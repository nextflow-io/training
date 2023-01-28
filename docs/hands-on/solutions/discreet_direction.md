<div class="formalpara-title">

**Solution**

</div>

``` nextflow
/*
 * Process 1C: Create the genome index file for STAR
 */

process '1C_prepare_star_genome_index' {

  input:
    path genome from params.genome 

  output:
    path 'genome_dir' into genome_dir_ch 

  script:
  """
  mkdir genome_dir 

  STAR --runMode genomeGenerate \
       --genomeDir genome_dir \
       --genomeFastaFiles ${genome} \
       --runThreadN ${task.cpus}
  """
}
```

- Take as input the `genome` file from the `params.genome` parameter.

- The `output` is a `file`\* called `genome_dir` and is added `into` a channel called `genome_dir_ch`. You can call the channel whatever you wish.

- Creates the output directory that will contain the resulting STAR genome index.

<div class="note">

\* The file in this case is a directory however it makes no difference.

</div>
