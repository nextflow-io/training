<div class="formalpara-title">

**Solution**

</div>

``` nextflow
/*
 * Process 1D: Create a file containing the filtered and recoded set of variants
 */

process '1D_prepare_vcf_file' {

  input:
    path variantsFile from params.variants 
    path blacklisted from params.blacklist 

  output:
    tuple path("${variantsFile.baseName}.filtered.recode.vcf.gz"), \
          path("${variantsFile.baseName}.filtered.recode.vcf.gz.tbi") into prepared_vcf_ch 

  script:
  """
  vcftools --gzvcf $variantsFile -c \
           --exclude-bed ${blacklisted} \
           --recode | bgzip -c \
           > ${variantsFile.baseName}.filtered.recode.vcf.gz 

  tabix ${variantsFile.baseName}.filtered.recode.vcf.gz 
  """
}
```

- Take as input the variants file, assigning the name `${variantsFile}`.

- Take as input the blacklisted file, assigning the name `${blacklisted}`.

- Out a tuple (or set) of two files into the `prepared_vcf_ch` channel.

- Defines the name of the first output file.

- Generates the secound output file (with `.tbi` suffix).
