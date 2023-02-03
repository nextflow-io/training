# Quality control

Next, we implement a `FASTQC` quality control step for your input reads (using the label `fastqc`). The inputs are the same as the read pairs used in the `QUANTIFICATION` step.

You can run it by using the following command:

```bash
nextflow run script5.nf -resume
```

Nextflow DSL2 knows to split the `reads_pair_ch` into two identical channels as they are required twice as an input for both of the `FASTQC` and the `QUANTIFICATION` process.

## MultiQC report

This step collects the outputs from the `QUANTIFICATION` and `FASTQC` processes to create a final report using the [MultiQC](http://multiqc.info/) tool.

Execute the next script with the following command:

```bash
nextflow run script6.nf -resume --reads 'data/ggal/*_{1,2}.fq'
```

It creates the final report in the `results` folder in the current `work` directory.

In this script, note the use of the [mix](https://www.nextflow.io/docs/latest/operator.html#mix) and [collect](https://www.nextflow.io/docs/latest/operator.html#collect) operators chained together to gather the outputs of the `QUANTIFICATION` and `FASTQC` processes as a single input. [Operators](https://www.nextflow.io/docs/latest/operator.html) can be used to combine and transform channels.

```groovy
MULTIQC(quant_ch.mix(fastqc_ch).collect())
```

We only want one task of MultiQC to be executed to produce one report. Therefore, we use the `mix` channel operator to combine the two channels followed by the `collect` operator, to return the complete channel contents as a single element.

## :material-check-all: Summary

In this step you have learned:

1. How to collect many outputs to a single input with the `collect` operator
2. How to `mix` two channels into a single channel
3. How to chain two or more operators together
