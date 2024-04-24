# Script 3

## Error message

```
WARN: Input tuple does not match tuple declaration in process `GATK_HAPLOTYPECALLER` -- offending value: [NA12878, /workspace/gitpod/troubleshoot/exercise3/work/aa/e0741cf287ce3a0c23f4dc12604b8c/reads_mother.bam, /workspace/gitpod/troubleshoot/exercise3/work/aa/e0741cf287ce3a0c23f4dc12604b8c/reads_mother.bam.bai]
ERROR ~ Error executing process > 'GATK_HAPLOTYPECALLER (1)'

Caused by:
  Not a valid path value: 'NA12878'


Tip: view the complete command output by changing to the process work dir and entering the command `cat .command.out`

 -- Check '.nextflow.log' file for details
```

```
tuple val(id), path(input_bam), path(input_bam_index)
```

??? Solution

    The warning can be broken down to help identify the cause: 

    I can be immediately seen that the `GATK_HAPLOTYPECALLER` has an offending value `[NA12878, /workspace/gitpod/troubleshoot/exercise3/work/aa/e0741cf287ce3a0c23f4dc12604b8c/reads_mother.bam, /workspace/gitpod/troubleshoot/exercise3/work/aa/e0741cf287ce3a0c23f4dc12604b8c/reads_mother.bam.bai]`. Specifically, the value 'NA12878' is not a valid path.

    There is only one tuple for the `GATK_HAPLOTYPECALLER` process:

    ```console title="hello.gatk.nf" linenums="48"
    tuple path(input_bam), path(input_bam_index)
    ```

    ```
    /*
    * Call variants with GATK HapolotypeCaller in GVCF mode
    */
    process GATK_HAPLOTYPECALLER {

        container "broadinstitute/gatk:4.5.0.0"

        input:
            tuple val(id), path(input_bam), path(input_bam_index)
            path ref_fasta
            path ref_index
            path ref_dict
            path interval_list

        output:
            tuple val(id), path("${input_bam}.g.vcf"), path("${input_bam}.g.vcf.idx")

        """
        gatk HaplotypeCaller \
            -R ${ref_fasta} \
            -I ${input_bam} \
            -O ${input_bam}.g.vcf \
            -L ${interval_list} \
            -ERC GVCF
        """
    }
    ```