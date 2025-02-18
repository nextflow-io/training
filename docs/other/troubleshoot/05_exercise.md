# Exercise 5

Move into the exercise 5 directory and execute the `hello-gatk.nf` script.

```bash
cd /workspaces/training/troubleshoot/exercise5
nextflow run hello-gatk.nf
```

!!! warning "Error message"

    ```console
    ERROR ~ Error executing process > 'GATK_JOINTGENOTYPING (1)'

    Caused by:
    Process `GATK_JOINTGENOTYPING (1)` terminated with an error exit status (2)

    Command executed:

    gatk GenomicsDBImport         --sample-name-map family_trio_map.tsv         --genomicsdb-workspace-path family_trio_gdb         -L intervals.list

    gatk GenotypeGVCFs         -R ref.fasta         -V gendb://family_trio_gdb         -O family_trio.joint.vcf         -L intervals.list

    Command exit status:
    2

    Command output:
    (empty)

    Command error:
    16:33:09.066 INFO  NativeLibraryLoader - Loading libgkl_compression.so from jar:file:/gatk/gatk-package-4.5.0.0-local.jar!/com/intel/gkl/native/libgkl_compression.so
    16:33:09.374 INFO  GenomicsDBImport - ------------------------------------------------------------
    16:33:09.380 INFO  GenomicsDBImport - The Genome Analysis Toolkit (GATK) v4.5.0.0
    16:33:09.380 INFO  GenomicsDBImport - For support and documentation go to https://software.broadinstitute.org/gatk/
    16:33:09.381 INFO  GenomicsDBImport - Executing as root@5f94bf7c2fe0 on Linux v6.1.75-060175-generic amd64
    16:33:09.381 INFO  GenomicsDBImport - Java runtime: OpenJDK 64-Bit Server VM v17.0.9+9-Ubuntu-122.04
    16:33:09.382 INFO  GenomicsDBImport - Start Date/Time: April 24, 2024 at 4:33:08 PM GMT
    16:33:09.383 INFO  GenomicsDBImport - ------------------------------------------------------------
    16:33:09.383 INFO  GenomicsDBImport - ------------------------------------------------------------
    16:33:09.385 INFO  GenomicsDBImport - HTSJDK Version: 4.1.0
    16:33:09.386 INFO  GenomicsDBImport - Picard Version: 3.1.1
    16:33:09.387 INFO  GenomicsDBImport - Built for Spark Version: 3.5.0
    16:33:09.388 INFO  GenomicsDBImport - HTSJDK Defaults.COMPRESSION_LEVEL : 2
    16:33:09.388 INFO  GenomicsDBImport - HTSJDK Defaults.USE_ASYNC_IO_READ_FOR_SAMTOOLS : false
    16:33:09.388 INFO  GenomicsDBImport - HTSJDK Defaults.USE_ASYNC_IO_WRITE_FOR_SAMTOOLS : true
    16:33:09.389 INFO  GenomicsDBImport - HTSJDK Defaults.USE_ASYNC_IO_WRITE_FOR_TRIBBLE : false
    16:33:09.389 INFO  GenomicsDBImport - Deflater: IntelDeflater
    16:33:09.389 INFO  GenomicsDBImport - Inflater: IntelInflater
    16:33:09.389 INFO  GenomicsDBImport - GCS max retries/reopens: 20
    16:33:09.389 INFO  GenomicsDBImport - Requester pays: disabled
    16:33:09.391 INFO  GenomicsDBImport - Initializing engine
    16:33:09.397 INFO  GenomicsDBImport - Shutting down engine
    [April 24, 2024 at 4:33:09 PM GMT] org.broadinstitute.hellbender.tools.genomicsdb.GenomicsDBImport done. Elapsed time: 0.01 minutes.
    Runtime.totalMemory()=201326592
    ***********************************************************************

    A USER ERROR has occurred: Bad input: Sample name map file must have 2 or 3 fields per line in the format:
    Sample        File
    or:
    Sample        File    Index
    but found line: "NA12877/t/workspaces/training/troubleshoot/exercise5/work/db/0d4c6d0c8bad080cec4a1e09217159/reads_father.bam.g.vcf/t/workspaces/training/troubleshoot/exercise5/work/db/0d4c6d0c8bad080cec4a1e09217159/reads_father.bam.g.vcf.idx/nNA12882/t/workspaces/training/troubleshoot/exercise5/work/a7/84fe9e9038bed36bdbf620847b47c2/reads_son.bam.g.vcf/t/workspaces/training/troubleshoot/exercise5/work/a7/84fe9e9038bed36bdbf620847b47c2/reads_son.bam.g.vcf.idx/nNA12878/t/workspaces/training/troubleshoot/exercise5/work/59/a653a3855cf6f6c70bc7b42aff0a0d/reads_mother.bam.g.vcf/t/workspaces/training/troubleshoot/exercise5/work/59/a653a3855cf6f6c70bc7b42aff0a0d/reads_mother.bam.g.vcf.idx/n" with 1 fields

    ***********************************************************************
    Set the system property GATK_STACKTRACE_ON_USER_EXCEPTION (--java-options '-DGATK_STACKTRACE_ON_USER_EXCEPTION=true') to print the stack trace.
    Using GATK jar /gatk/gatk-package-4.5.0.0-local.jar
    Running:
        java -Dsamjdk.use_async_io_read_samtools=false -Dsamjdk.use_async_io_write_samtools=true -Dsamjdk.use_async_io_write_tribble=false -Dsamjdk.compression_level=2 -jar /gatk/gatk-package-4.5.0.0-local.jar GenomicsDBImport --sample-name-map family_trio_map.tsv --genomicsdb-workspace-path family_trio_gdb -L intervals.list

    Work dir:
    /workspaces/training/troubleshoot/exercise5/work/a3/fa7fab1e9d1a1e0ec2d65e87b9cebf

    Tip: view the complete command output by changing to the process work dir and entering the command `cat .command.out`

    -- Check '.nextflow.log' file for details
    ```

??? Solution

    The error message from this execution can be broken down to give clues.

    ```
    Command exit status:
    2
    ```

    Exit code `2` suggests invalid usage of a shell built-in command. Examples of built-in commands include alias, echo, and printf. Alternatively, it could mean that you are trying to access a file or directory that doesn't exist or requires permissions.

    ```
    A USER ERROR has occurred: Bad input: Sample name map file must have 2 or 3 fields per line in the format
    ```

    The user error `Bad input: Sample name map` suggests that the sample name map is the problematic.

    ```console title="hello-gatk.nf" linenums="87"
    gatk GenomicsDBImport \
        --sample-name-map ${sample_map} \
        --genomicsdb-workspace-path ${cohort_name}_gdb \
        -L ${interval_list}

    gatk GenotypeGVCFs \
        -R ${ref_fasta} \
        -V gendb://${cohort_name}_gdb \
        -O ${cohort_name}.joint.vcf \
        -L ${interval_list}
    ```

    The [GATK GenomicsDBImport documentation](https://gatk.broadinstitute.org/hc/en-us/articles/360036883491-GenomicsDBImport) can be used to help us interstand the GenomicsDBImport specific error.

    | Argument name     | Default value   | Summary                            |
    | ----------------- | ----------------| -----------------------------------|
    | --sample-name-map | null            | Path to file containing a mapping of sample name to file uri in tab delimited format. If this is specified then the header from the first sample will be treated as the merged header rather than merging the headers, and the sample names will be taken from this file. This may be used to rename input samples. This is a performance optimization that relaxes the normal checks for consistent headers. Using vcfs with incompatible headers may result in silent data corruption. |

    Although the documentation does not explicitly describe the error, it highlights the importance of the structure of the map, something that is also suggested in the error message.

    By reviewing the code where the `sample_map` channel is created, it can be seen that the map has `/` rather than `\` for the tab delimiter.

    ```console title="hello-gatk.nf" linenums="119"
    // Create a sample map of the output GVCFs
    sample_map = GATK_HAPLOTYPECALLER.out.collectFile(){ id, gvcf, idx ->
            ["${params.cohort_name}_map.tsv", "${id}/t${gvcf}/t${idx}/n"]
    }
    ```

    By correcting this the error will be resolved.

    ```console title="hello-gatk.nf" linenums="119"
    // Create a sample map of the output GVCFs
    sample_map = GATK_HAPLOTYPECALLER.out.collectFile(){ id, gvcf, idx ->
            ["${params.cohort_name}_map.tsv", "${id}\t${gvcf}\t${idx}\n"]
    }
    ```
