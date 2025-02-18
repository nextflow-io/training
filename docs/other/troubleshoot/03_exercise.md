# Exercise 3

Move into the exercise 3 directory and execute the `hello-gatk.nf` script.

```bash
cd /workspaces/training/troubleshoot/exercise3
nextflow run hello-gatk.nf
```

!!! warning "Error message"

    ```console
    WARN: Input tuple does not match tuple declaration in process `GATK_HAPLOTYPECALLER` -- offending value: [NA12878, /workspaces/training/troubleshoot/exercise3/work/aa/e0741cf287ce3a0c23f4dc12604b8c/reads_mother.bam, /workspaces/training/troubleshoot/exercise3/work/aa/e0741cf287ce3a0c23f4dc12604b8c/reads_mother.bam.bai]
    ERROR ~ Error executing process > 'GATK_HAPLOTYPECALLER (1)'

    Caused by:
    Not a valid path value: 'NA12878'


    Tip: view the complete command output by changing to the process work dir and entering the command `cat .command.out`

    -- Check '.nextflow.log' file for details
    ```

??? Solution

    This warning can be broken down to help identify the cause:

    You can immediately see that the `GATK_HAPLOTYPECALLER` has an offending value `[NA12878, /workspaces/training/troubleshoot/exercise3/work/aa/e0741cf287ce3a0c23f4dc12604b8c/reads_mother.bam, /workspaces/training/troubleshoot/exercise3/work/aa/e0741cf287ce3a0c23f4dc12604b8c/reads_mother.bam.bai]`.

    _Your offending value may be NA12878, NA12877, or NA12882_

    In this example, the value 'NA12878' is not a valid path.

    There is only one tuple for the `GATK_HAPLOTYPECALLER` process:

    ```console title="hello.gatk.nf" linenums="48"
    tuple path(input_bam), path(input_bam_index)
    ```

    It can be seen that this tuple is made up of two paths and offending value has one value and two paths.

    It can be concluded that the tuple has a missing value in the first position.

    To resolve this error a `val` must be added to this input.

    The name of this `val` might be used in the script block and must be considered.

    While all of the variables in the script block are accounted for, the `id` value in the output has not been specified.

    It can be concluded that the input val should be `val(id)` to match the output.

    ```console title="hello-gatk.nf" linenums="40"
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
