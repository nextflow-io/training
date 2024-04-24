# Script 2

## Error message

`GATK_HAPLOTYPECALLER` only runs once

```
executor >  local (5)
[d2/839095] process > SAMTOOLS_INDEX (1)       [100%] 3 of 3 ✔
[1a/2dafee] process > GATK_HAPLOTYPECALLER (1) [100%] 1 of 1 ✔
[59/278fc4] process > GATK_JOINTGENOTYPING (1) [100%] 1 of 1 ✔
```

Commonly a channel is the wrong type.

## Solution

??? Solution

    ```
        // Call variants from the indexed BAM file
    GATK_HAPLOTYPECALLER(
        SAMTOOLS_INDEX.out,
        params.genome_reference,
        params.genome_reference_index,
        params.genome_reference_dict,
        params.calling_intervals
    )
    ```