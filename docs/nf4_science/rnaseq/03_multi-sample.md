# Part 3: Multi-sample paired-end implementation of RNAseq processing

TODO: short blurb

- Make the workflow accept multiple input samples and parallelize execution
- Switch to paired-end RNAseq data

---

## 1. Make the workflow accept multiple input samples and parallelize execution

TODO

## 1.1. Modify the input channel factory to accept a CSV file of file paths as input

TODO

## 1.2. Apply the `collect()` operator to the MULTIQC process to aggregate QC reports across all samples

TODO

## 1.3. Update the default value for the `report_id` parameter

TODO

## 1.4. Run the workflow to test that it works

TODO

---

## 2. Switch to paired-end RNAseq data

TODO

## 2.1. Modify channel factory to grab the second file of reads

TODO

## 2.2. Generalize the FASTQC process

TODO

## 2.3. Adapt the TRIM_GALORE process to expect paired-end reads

TODO

## 2.4. Adapt the HISAT2 process to expect paired-end reads

TODO

## 2.5. Run the workflow to test that it works

TODO

---

### Takeaway

You know how to adapt a single-sample workflow to parallelize processing of multiple samples, and to use paired-end read data.

### What's next?

Celebrate your success and take a well deserved break!

Next, we ask you to complete a very short survey about your experience with this training course, then we'll take you to a page with links to further training resources and helpful links.
