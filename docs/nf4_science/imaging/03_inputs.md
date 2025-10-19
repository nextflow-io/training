# Part 3: Understanding Inputs

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
input: 'data/samplesheet.csv'
outdir: 'results'
mindagap_tilesize: 90
mindagap_boxsize: 7
mindagap_loopnum: 100
clahe_pyramid_tile: 368
segmentation_method: "cellpose"
```

Now your command becomes:

```bash
nextflow run ./molkart -params-file params.yaml
```

That's it! The parameter file documents your exact configuration and makes it easy to rerun or share.

### 1.3. Overriding parameters

You can still override specific parameters from the command line:

```bash
nextflow run ./molkart -params-file params.yaml --outdir test_results
```

This is useful for quick tests without editing the file.


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

Let's modify the samplesheet to reference the local files we downloaded earlier.:

```csv title="samplesheet.csv"
sample,nuclear_image,spot_table,membrane_image
mem_only,data/nuclear.tiff,data/spots.txt,data/membrane.tiff
```

!!! Warning

    Notice the paths in the samplesheet are relative to where you **run** Nextflow, not where the samplesheet is located.

### Takeaway

Samplesheets organize multi-sample datasets in a way that allows you to explicitly define your metadata along with the file paths.
Most nf-core pipelines use this pattern.

### What's next?

Learn how to create training data for custom cell segmentation models.

---

## 3. Training custom segmentation models

### 3.1. Why custom models?

The molkart pipeline uses pre-trained machine learning models (like Cellpose) to segment cells.
These pre-trained models work well for common cell types, but your cells might look different:

- Unusual cell morphology
- Different staining patterns
- Non-standard imaging conditions

nf-core/molkart contains functionality for creating training subsets to help you train a **custom cellpose model** on your specific data.

### 3.2. The training workflow

Here's the process for training a custom model:

1. **Generate training images**: Use molkart's `create_training_subset` to create small cropped images
2. **Annotate in Cellpose**: Load crops into Cellpose GUI and manually outline cells
3. **Train custom model**: Let Cellpose learn from your annotations
4. **Use custom model**: Run molkart with your trained model

For time purposes, and to showcase how parameter files help manage inputs, we'll only focus on Step 1 here.

### 3.3. Creating training crops

Create a parameter file for generating training images:

```yaml title="params-training-subset.yaml"
input: 'data/samplesheet.csv'
outdir: 'results-training-crops'

# Preprocessing (same as before)
mindagap_tilesize: 90
mindagap_boxsize: 7
mindagap_loopnum: 100
clahe_pyramid_tile: 368

# Training subset parameters
create_training_subset: true
crop_size_x: 30
crop_size_y: 30
crop_amount: 4
```

Run it:

```bash
nextflow run ./molkart -params-file params-training-subset.yaml -resume
```

This creates small image crops in your output directory that you can load into Cellpose for annotation and training.

### 3.5. Next steps (outside this course)

If you were to continue, the next steps would be:

1. Install Cellpose on your local machine
2. Load the cropped images into the Cellpose GUI
3. Manually draw cell outlines on ~50-100 cells
4. Train a custom model using Cellpose's training function
5. Run molkart with `--cellpose_model /path/to/your/custom/model`

Importantly you can continue using nextflow's `-resume` functionality in between executions so that all the pre-processing steps are re-used between executions.

### Takeaway

When running Nextflow pipelines the combination of parameter files, samplesheets, and resume functionality makes it easy to manage complex iterative research patterns with a minimal amount of overhead.

### What's next?

Now that we've covered inputs, let's explore how to configure Nextflow pipelines for different computing environments.
