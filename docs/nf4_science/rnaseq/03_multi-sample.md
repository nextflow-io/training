# Part 3: Multi-sample paired-end implementation of RNAseq processing

TODO: short blurb

- Make the workflow accept multiple input samples and parallelize execution
- Add comprehensive QC report generation
- Switch to paired-end RNAseq data

---

## 1. Make the workflow accept multiple input samples and parallelize execution

TODO

## 1.1. Modify the input channel factory to accept a CSV file of file paths as input

TODO

## 1.2. Run the workflow to test that it works

TODO

---

## 2. Aggregate pre-processing QC metrics into a single MultiQC report

TODO

### 2.1. Describe the MultiQC process

Let's write a process, which we'll call `MULTIQC`, that collect QC metrics with MultiQC in a generic way.

```groovy title="modules/trim_galore.nf" linenums="1"
#!/usr/bin/env nextflow

process MULTIQC {

    container "community.wave.seqera.io/library/pip_multiqc:ad8f247edb55897c"
    publishDir "results/multiqc", mode: 'copy'

    input:
    path '*'
    val output_name

    output:
    path "${output_name}.html", emit: report
    path "${output_name}_data", emit: data

    script:
    """
    multiqc . -n ${output_name}.html
    """
}
```

### 2.2. Import the module into the workflow file

TODO

### 2.3. Call the process on the outputs of the previous steps

TODO

### 2.4. Set a default value for the `report_id` parameter

TODO

### 2.5. Run the workflow to test that it works

TODO

---

<!-- Bonus if there's time -->

## 3. Switch to paired-end RNAseq data

TODO

## 3.1. Modify channel factory to grab the second file of reads

TODO

## 3.2. Generalize the FASTQC process

TODO

## 3.3. Adapt the TRIM_GALORE process to expect paired-end reads

TODO

## 3.4. Adapt the HISAT2 process to expect paired-end reads

TODO

## 3.5. Run the workflow to test that it works

TODO

---

### Takeaway

You know how to adapt a single-sample workflow to parallelize processing of multiple samples, generate a comprehensive QC report and switch to using paired-end read data if needed.

_Making the workflow accept either data type on the fly is out of scope for this training._

### What's next?

Congratulations, you've completed the Nextflow For RNAseq mini-course! Celebrate your success and take a well deserved break!

Next, we ask you to complete a very short survey about your experience with this training course, then we'll take you to a page with links to further training resources and helpful links.
