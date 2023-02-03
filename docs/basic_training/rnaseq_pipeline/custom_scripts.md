# Custom scripts

Real-world pipelines use a lot of custom user scripts (BASH, R, Python, etc.) Nextflow allows you to consistently use and manage these scripts. Simply put them in a directory named `bin` in the pipeline project root. They will be automatically added to the pipeline execution `PATH`.

For example, create a file named `fastqc.sh` with the following content:

```bash
!/bin/bash
set -e
set -u

sample_id=${1}
reads=${2}

mkdir fastqc_${sample_id}_logs
fastqc -o fastqc_${sample_id}_logs -f fastq -q ${reads}
```

Save it, give execute permission, and move it into the `bin` directory as shown below:

```bash
chmod +x fastqc.sh
mkdir -p bin
mv fastqc.sh bin
```

Then, open the `script7.nf` file and replace the `FASTQC` process’ script with the following code:

```groovy
script:
"""
fastqc.sh "$sample_id" "$reads"
"""
```

Run it as before:

```bash
nextflow run script7.nf -resume --reads 'data/ggal/*_{1,2}.fq'
```

## :material-check-all: Summary

In this step you have learned:

1. How to write or use existing custom scripts in your Nextflow pipeline.
2. How to avoid the use of absolute paths by having your scripts in the `bin/` folder.
