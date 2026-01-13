# Part 3: Organizing inputs

In Part 2, we ran molkart with multiple parameters on the command line.
Now we'll learn two better approaches for managing inputs: **parameter files** and **samplesheets**.

## 1. Using parameter files

### 1.1. The problem with long command lines

Recall our command from Part 2:

```bash
nextflow run ./molkart \
  --input 'data/samplesheet.csv' \
  --mindagap_tilesize 90 \
  --mindagap_boxsize 7 \
  --mindagap_loopnum 100 \
  --clahe_pyramid_tile 368 \
  --segmentation_method "cellpose" \
  --outdir results
```

This works, but it's hard to reproduce, share, or modify.
What if you need to run the same analysis again next month?
What if a collaborator wants to use your exact settings?

### 1.2. Solution: Use a parameter file

Create a file called `params.yaml`:

```yaml title="params.yaml"
input: "data/samplesheet.csv"
outdir: "results"
mindagap_tilesize: 90
mindagap_boxsize: 7
mindagap_loopnum: 100
clahe_pyramid_tile: 368
segmentation_method: "cellpose"
```

Now your command becomes:

```bash
nextflow run ./molkart -params-file params.yaml -resume
```

That's it! The parameter file documents your exact configuration and makes it easy to rerun or share.

### 1.3. Overriding parameters

You can still override specific parameters from the command line:

```bash
nextflow run ./molkart -params-file params.yaml --segmentation_method "stardist" --outdir stardist_results -resume
```

The above line changes the `segmentation_method` to be `stardist` and the `--outdir` name to be `stardist_results` instead of the params in the `params.yaml` file.
Additionally, you can see that the `-resume` flag allowed us to reuse the pre-processing results from the previous run, saving time.
You can use this pattern to quickly test different variations of the pipeline.

### Takeaway

Parameter files make your analyses reproducible and easy to share.
Use them for any real analysis work.

### What's next?

Learn how samplesheets organize information about multiple samples.

---

## 2. The samplesheet pattern

### 2.1. What is a samplesheet?

A samplesheet is a CSV file that describes your input samples.
Each row is a sample, and columns specify the files and metadata for that sample.

Let's look at the samplesheet we've been using:

```bash
cat data/samplesheet.csv
```

```csv title="samplesheet.csv"
sample,nuclear_image,spot_table,membrane_image
mem_only,https://raw.githubusercontent.com/nf-core/test-datasets/molkart/test_data/input_data/nuclear.tiff,https://raw.githubusercontent.com/nf-core/test-datasets/molkart/test_data/input_data/spots.txt,https://raw.githubusercontent.com/nf-core/test-datasets/molkart/test_data/input_data/membrane.tiff
```

The columns are:

- `sample`: Unique sample identifier
- `nuclear_image`: Nuclear staining image (TIFF)
- `spot_table`: Transcript spots (TXT)
- `membrane_image`: Membrane staining image (TIFF, optional)

### 2.2. File paths

Samplesheets accept multiple path types:

- **URLs**: Nextflow downloads automatically (as shown above)
- **Local paths**: `data/nuclear.tiff` or `/absolute/path/to/nuclear.tiff`
- **Cloud storage**: `s3://bucket/nuclear.tiff`, `gs://bucket/nuclear.tiff`, `az://container/nuclear.tiff`

You can mix path types in the same samplesheet.

### 2.3. Creating your own samplesheet

First, let's download the test data files locally:

```bash
cd /workspaces/training/nf4-science/imaging/data
curl -O https://raw.githubusercontent.com/nf-core/test-datasets/molkart/test_data/input_data/nuclear.tiff
curl -O https://raw.githubusercontent.com/nf-core/test-datasets/molkart/test_data/input_data/spots.txt
curl -O https://raw.githubusercontent.com/nf-core/test-datasets/molkart/test_data/input_data/membrane.tiff
cd ..
```

Now let's modify the samplesheet to reference these local files:

```csv title="samplesheet.csv"
sample,nuclear_image,spot_table,membrane_image
mem_only,data/nuclear.tiff,data/spots.txt,data/membrane.tiff
```

!!! Warning

    Notice the paths in the samplesheet are relative to where you **run** Nextflow, not where the samplesheet is located.

Finally, let's execute nf-core/molkart one more time with the samplesheet with local file paths:

`nextflow run ./molkart -params-file params.yaml -resume`

As you can see, Nextflow executes this run similar to when the files were downloaded from Github. This is one of the great features of Nextflow, it stages the data properly for you, regardless of where it is located.

### Takeaway

Samplesheets organize multi-sample datasets in a way that allows you to explicitly define your metadata along with the file paths.
Most nf-core pipelines use this pattern.

### What's next?

Now that we've covered inputs, let's explore how to configure Nextflow pipelines for different computing environments.
