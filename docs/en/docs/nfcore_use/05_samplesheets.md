# Part 5: Prepare a samplesheet

In Parts 1 through 3, you ran pipelines with their test profiles, which supply a ready-made input samplesheet.
To run a pipeline on your own data, you need to write that samplesheet yourself.

Most nf-core pipelines take their main input as a samplesheet: a CSV file in which each row describes a sample (or one sequencing run of a sample) and each column provides a piece of information the pipeline needs, such as a sample name or the path to a FASTQ file.
Each pipeline defines its own set of columns, so the first step is always to find out what the pipeline expects.

In this part, you write a samplesheet for `nf-core/rnaseq` from a set of FASTQ files, then check it without running the full pipeline.

---

## 1. Find out what the pipeline expects

Every nf-core pipeline documents its samplesheet format in two places: the usage documentation, written for humans, and the input schema, which the pipeline uses to check the samplesheet at launch.

### 1.1. Read the usage documentation

Open the [nf-core/rnaseq usage documentation](https://nf-co.re/rnaseq/usage) in your browser and find the **Samplesheet input** section.

It tells you that the samplesheet must be a comma-separated file with the following columns:

| Column         | Description                                                                      |
| -------------- | -------------------------------------------------------------------------------- |
| `sample`       | Sample name. Rows with the same sample name are merged into a single sample.     |
| `fastq_1`      | Path to the FASTQ file for read 1.                                               |
| `fastq_2`      | Path to the FASTQ file for read 2. Leave empty for single-end data.              |
| `strandedness` | Library strandedness: `forward`, `reverse`, `unstranded`, or `auto` to infer it. |

The documentation also explains two behaviors that matter when you prepare real data:

- **Multiple rows per sample**: if a sample was sequenced more than once (for example on several lanes), you list each run on its own row with the same sample name, and the pipeline concatenates the reads before processing.
- **Strandedness inference**: if you don't know the strandedness of your library, `auto` makes the pipeline infer it from a subset of the reads.

### 1.2. Examine the input schema

The usage documentation is written by the pipeline developers, so it can fall behind the code.
The input schema is what the pipeline actually enforces.

The input schema lives in the pipeline repository at `assets/schema_input.json`.
You can view it on GitHub at [github.com/nf-core/rnaseq/blob/master/assets/schema_input.json](https://github.com/nf-core/rnaseq/blob/master/assets/schema_input.json).

??? abstract "Input schema (excerpt)"

    ```json title="assets/schema_input.json"
    {
        "items": {
            "type": "object",
            "properties": {
                "sample": {
                    "type": ["string", "integer"],
                    "pattern": "^\\S+$",
                    "errorMessage": "Sample name must be provided and cannot contain spaces",
                    "meta": ["id"]
                },
                "fastq_1": {
                    "type": "string",
                    "format": "file-path",
                    "exists": true,
                    "pattern": "^([\\S\\s]*\\/)?[^\\s\\/]+\\.f(ast)?q(\\.gz)?$",
                    "errorMessage": "FastQ file for reads 1 must be provided, cannot contain spaces and must have extension '.fq', '.fastq', '.fq.gz' or '.fastq.gz'"
                },
                "fastq_2": {
                    "type": "string",
                    "format": "file-path",
                    "exists": true,
                    "pattern": "^([\\S\\s]*\\/)?[^\\s\\/]+\\.f(ast)?q(\\.gz)?$",
                    "errorMessage": "FastQ file for reads 2 cannot contain spaces and must have extension '.fq', '.fastq', '.fq.gz' or '.fastq.gz'"
                },
                "strandedness": {
                    "type": "string",
                    "errorMessage": "Strandedness must be provided and be one of 'auto', 'forward', 'reverse' or 'unstranded'",
                    "enum": ["forward", "reverse", "unstranded", "auto"],
                    "meta": ["strandedness"]
                }
            },
            "required": ["sample", "fastq_1", "strandedness"]
        }
    }
    ```

The schema states the same rules as the documentation, in a form the pipeline can check automatically:

- `required` lists the columns that must have a value in every row: `sample`, `fastq_1` and `strandedness`.
  `fastq_2` is optional, which is how the pipeline supports single-end data.
- `pattern` constrains the format of a value: sample names cannot contain spaces, and FASTQ paths must end in `.fq`, `.fastq`, `.fq.gz` or `.fastq.gz`.
- `exists: true` means the file at that path must exist when the pipeline launches.
- `enum` lists the only accepted values for `strandedness`.

The full schema also defines a few optional columns (such as `seq_platform` and `seq_center`) that you can ignore for now.

### 1.3. Examine the test samplesheet

The samplesheet used by the `test` profile in Part 3 is a good example of these rules in practice.
Its location is set by the `input` parameter in the pipeline's `conf/test.config` file.

??? abstract "Test samplesheet (excerpt; paths shortened)"

    ```csv title="samplesheet_test.csv"
    sample,fastq_1,fastq_2,strandedness
    WT_REP1,https://.../SRR6357070_1.fastq.gz,https://.../SRR6357070_2.fastq.gz,auto
    WT_REP1,https://.../SRR6357071_1.fastq.gz,https://.../SRR6357071_2.fastq.gz,auto
    WT_REP2,https://.../SRR6357072_1.fastq.gz,https://.../SRR6357072_2.fastq.gz,reverse
    RAP1_UNINDUCED_REP1,https://.../SRR6357073_1.fastq.gz,,reverse
    RAP1_UNINDUCED_REP2,https://.../SRR6357074_1.fastq.gz,,reverse
    RAP1_UNINDUCED_REP2,https://.../SRR6357075_1.fastq.gz,,reverse
    ```

This excerpt shows several features of the format:

- `WT_REP1` and `RAP1_UNINDUCED_REP2` each appear on two rows, so each will be processed as a single sample with the reads of both runs combined.
- The `RAP1_UNINDUCED` samples are single-end: their `fastq_2` field is empty.
- File paths can be URLs.
  Nextflow downloads remote files automatically, and the same applies to cloud storage paths such as `s3://` or `gs://`.

### Takeaway

You know how to find the samplesheet requirements of an nf-core pipeline, both in its usage documentation and in its input schema.

### What's next?

Write a samplesheet for a set of FASTQ files.

---

## 2. Write a samplesheet for your data

The `data/reads/` directory contains FASTQ files for a small RNA-seq experiment, named the way a sequencing facility typically delivers them.

### 2.1. Survey the data files

List the files:

```bash
tree data/reads
```

```console title="Directory contents"
data/reads
├── KO_1_S3_L001_R1_001.fastq.gz
├── KO_1_S3_L001_R2_001.fastq.gz
├── KO_2_S4_L001_R1_001.fastq.gz
├── KO_2_S4_L001_R2_001.fastq.gz
├── WT_1_S1_L001_R1_001.fastq.gz
├── WT_1_S1_L001_R2_001.fastq.gz
├── WT_1_S1_L002_R1_001.fastq.gz
├── WT_1_S1_L002_R2_001.fastq.gz
├── WT_2_S2_L001_R1_001.fastq.gz
└── WT_2_S2_L001_R2_001.fastq.gz
```

The filenames follow the Illumina convention `<sample>_S<number>_L<lane>_R<read>_001.fastq.gz`, which gives you everything you need to fill in the samplesheet:

- There are four samples: two wild-type (`WT_1`, `WT_2`) and two knockout (`KO_1`, `KO_2`).
- Every sample has an `R1` and an `R2` file, so the data is paired-end.
- `WT_1` was sequenced on two lanes (`L001` and `L002`), so it needs two rows.

The strandedness of the library is not recorded in the filenames.
When you don't have that information from the sequencing facility, use `auto`.

### 2.2. Write the samplesheet

Create a file named `samplesheet.csv` in your working directory and fill it in based on the information above.

??? exercise

    Write `samplesheet.csv` with a header row and one row per sequencing run.
    Use paths relative to your working directory, for example `data/reads/WT_2_S2_L001_R1_001.fastq.gz`.

??? solution

    ```csv title="samplesheet.csv"
    --8<-- "nfcore-use/solutions/samplesheet.csv"
    ```

    `WT_1` has two rows, one per lane, and the pipeline merges them into a single sample.

Relative paths are resolved from the directory where you launch the pipeline.
If you plan to launch from a different directory, or share the samplesheet with collaborators, use absolute paths instead.

??? tip "Generating a samplesheet with a script"

    Writing a samplesheet by hand is fine for a handful of files, but error-prone for a real experiment with dozens of samples.
    A short shell loop can generate the rows from the filenames:

    ```bash
    echo "sample,fastq_1,fastq_2,strandedness" > samplesheet.csv
    for r1 in data/reads/*_R1_001.fastq.gz; do
        r2=${r1/_R1_/_R2_}
        sample=$(basename $r1 | sed 's/_S[0-9]*_L[0-9]*_R1_001.fastq.gz//')
        echo "$sample,$r1,$r2,auto" >> samplesheet.csv
    done
    ```

    Always check the result before using it, since this kind of script depends on every file following the same naming convention.

### 2.3. Check the samplesheet without running the pipeline

The full pipeline takes a while to run, so it is useful to check the samplesheet on its own first.
The `-preview` option runs the pipeline's launch steps, including input validation, without executing any tasks.

Run the pipeline with your samplesheet and the `-preview` option:

```bash
nextflow run nf-core/rnaseq -profile test,docker --input samplesheet.csv --outdir results -preview
```

The `test` profile provides the reference files, and `--input` replaces its samplesheet with yours.

??? success "Command output"

    ```console
    ...
    -[nf-core/rnaseq] Pipeline completed successfully-
    ```

The pipeline reports success without running any tasks, which means your samplesheet passed validation.
If you made a mistake, the pipeline stops with a list of errors; Part 6 covers how to interpret and fix them.

??? tip "Running the full pipeline on this data"

    The files only contain 2,000 read pairs each, which is below the pipeline's default minimum of 10,000 reads per sample after trimming.
    Samples below that threshold are skipped, so lower it with a params file:

    ```yaml title="small_reads.yml"
    min_trimmed_reads: 1000
    ```

    Then remove `-preview` and add `-c laptop.config` as in Part 3:

    ```bash
    nextflow run nf-core/rnaseq -profile test,docker --input samplesheet.csv --outdir results -c laptop.config -params-file small_reads.yml
    ```

    The parameter must be set in a params file rather than as `--min_trimmed_reads 1000` on the command line, because command-line values are passed as strings and the pipeline expects an integer.

### Takeaway

You can write a samplesheet for a set of FASTQ files, including samples sequenced over multiple lanes, and check it with `-preview` before committing to a full run.

### What's next?

Generate a samplesheet for data from a public archive.

---

## 3. Generate a samplesheet for public data

Analyses often reuse public data, for example to add samples to a study or to reproduce a published result.
The `nf-core/fetchngs` pipeline downloads sequencing data and metadata from public archives (SRA, ENA, DDBJ and GEO), given a list of accession numbers.
It can also write a samplesheet formatted for a downstream nf-core pipeline, including `nf-core/rnaseq`.

### 3.1. Download data with fetchngs

The file `data/sra_ids.csv` lists four sequencing runs, one per line:

```csv title="data/sra_ids.csv"
--8<-- "nfcore-use/data/sra_ids.csv"
```

These runs come from a drug-screening study ([GEO series GSE206599](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE206599)) in which human cells were treated with the MEK inhibitor trametinib or with DMSO as a control.
They were chosen because their files are only a few megabytes each.

Run `nf-core/fetchngs` with these accessions, and request a samplesheet for `nf-core/rnaseq`:

```bash
nextflow run nf-core/fetchngs -profile docker --input data/sra_ids.csv --outdir fetchngs-results --nf_core_pipeline rnaseq
```

The pipeline downloads the FASTQ files from the ENA, so it needs access to the internet.
It finishes in about a minute, then prints a warning:

```console title="Command output (excerpt)"
WARN: =============================================================================
  Please double-check the samplesheet that has been auto-created by the pipeline.
  Public databases don't reliably hold information such as strandedness
  information, controls etc

  All of the sample metadata obtained from the ENA has been appended
  as additional columns to help you manually curate the samplesheet before
  running nf-core/other pipelines.
===================================================================================
```

The results contain the FASTQ files, the metadata for each run, and the samplesheet:

```console title="Directory contents (excerpt)"
fetchngs-results
├── fastq
│   ├── SRX15815511_SRR19771055_1.fastq.gz
│   ├── SRX15815511_SRR19771055_2.fastq.gz
│   └── ...
├── metadata
├── pipeline_info
└── samplesheet
    └── samplesheet.csv
```

### 3.2. Examine the generated samplesheet

Open `fetchngs-results/samplesheet/samplesheet.csv`.
Here are its first four columns and one of the metadata columns, `sample_title`, for the first two rows:

| `sample`      | `fastq_1`                                                   | `fastq_2`                                                   | `strandedness` | `sample_title`                            |
| ------------- | ----------------------------------------------------------- | ----------------------------------------------------------- | -------------- | ----------------------------------------- |
| `SRX15815511` | `fetchngs-results/fastq/SRX15815511_SRR19771055_1.fastq.gz` | `fetchngs-results/fastq/SRX15815511_SRR19771055_2.fastq.gz` | `auto`         | `22s000068_Trametinib_DMSO_countsPerGene` |
| `SRX15815525` | `fetchngs-results/fastq/SRX15815525_SRR19771041_1.fastq.gz` | `fetchngs-results/fastq/SRX15815525_SRR19771041_2.fastq.gz` | `auto`         | `22s000075_Trametinib_DMSO_countsPerGene` |

The first four columns are the ones `nf-core/rnaseq` requires:

- **`sample`** is the experiment accession, which identifies the sample but says nothing about it.
  The experimental condition only appears in metadata columns such as `sample_title`.
- **`fastq_1` and `fastq_2`** point to the downloaded files, relative to the directory where you ran `nf-core/fetchngs`.
- **`strandedness`** is set to `auto`, because the archives do not record it reliably.
  You can set a different value with the `--nf_core_rnaseq_strandedness` parameter.

The remaining 31 columns contain metadata from the ENA, such as the instrument model, read counts and study title.

The samplesheet is valid as it is.
You can check it with `-preview`:

```bash
nextflow run nf-core/rnaseq -profile test,docker --input fetchngs-results/samplesheet/samplesheet.csv --outdir results -preview
```

```console title="Command output (excerpt)"
WARN: Found the following unidentified headers in fetchngs-results/samplesheet/samplesheet.csv: - run_accession - experiment_accession - ...
```

The pipeline passes validation, with a warning listing the columns it does not use.

### 3.3. Curate the samplesheet

The generated samplesheet works, but the accession numbers make the pipeline's outputs hard to interpret: every report and output file will be named after them.
Curating the samplesheet means giving the samples meaningful names and keeping only the columns the pipeline uses.

??? exercise

    Create a file named `samplesheet_public.csv` with the four columns required by `nf-core/rnaseq`.
    Name each sample after its treatment, using the `sample_title` column to find out which is which.
    Check the result with `-preview`.

??? solution

    ```csv title="samplesheet_public.csv"
    --8<-- "nfcore-use/solutions/samplesheet_public.csv"
    ```

    The `-preview` run now completes without the warning about unidentified headers.

Keep the generated samplesheet alongside the curated one: its metadata columns record where each sample came from.

??? exercise "Apply it to the pipeline you chose in Part 4"

    Find the samplesheet requirements of the pipeline you chose in [Part 4](./04_choose_pipeline.md), in its usage documentation and in its `assets/schema_input.json` file on GitHub.
    Then write the header and one example row of a samplesheet for your scenario, using placeholder paths for the data files.

    Note any column whose value your scenario does not provide, and where you would get that information.

### Takeaway

`nf-core/fetchngs` downloads public data and writes a samplesheet for a downstream pipeline.
The generated samplesheet is a starting point: check the values the archives cannot provide reliably, and rename the samples before running the analysis.

### What's next?

Learn to interpret and fix the errors the pipeline reports when a samplesheet is invalid.

---

## Summary

In this part you learned to:

- Find a pipeline's samplesheet requirements in its usage documentation and input schema
- Write a samplesheet for paired-end data, including samples sequenced on multiple lanes
- Check a samplesheet with `-preview` without running the pipeline
- Generate a samplesheet for public data with `nf-core/fetchngs`, and curate it before use
