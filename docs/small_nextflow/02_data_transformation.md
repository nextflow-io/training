# Part 2: Data Transformation & Analysis

In this part, we'll build a multi-step workflow that classifies images using machine learning, manages computational resources, and groups results intelligently.

---

## Classification

Let's get to the fun part - the cat sorting!
We have a little classification script - `classify.py` that I've provided in the `.stuff` directory.
In your research sometimes you have small accessory scripts that are useful for your pipelines.
We're using a python script here in this workshop example, but this pattern will hold for scripts written in perl, ruby, R, python, closurescript, or any of the other interpreted languages.

### Set up the classification script

Let's pull the file out into a new `bin` directory:

```bash
mkdir -p bin
cp .stuff/classify.py bin/
```

The script requires some dependencies.
Again, we'll do this the slow/painful way one time before we demonstrate how to use containers to encapsulate the software dependencies.

We'll grab one more file from our `.stuff` directory - a `pyproject.toml` file which is a way of describing software dependencies for Python projects.
This is unrelated to Nextflow, but an example of one of the (many) ways in which different languages and frameworks might install software.

You can install the dependencies and activate the environment with:

```bash
cp .stuff/pyproject.toml .
uv sync
source .venv/bin/activate
```

which you can run with:

```bash
bin/classify.py --help
```

```console title="Output"
usage: classify.py [-h] [--model-path MODEL_PATH] [--labels LABELS [LABELS ...]] [--json] image

Classify a single image using MetaCLIP

positional arguments:
  image                 Path to the image file to classify

options:
  -h, --help            show this help message and exit
  --model-path MODEL_PATH
                        Path to MetaCLIP model weights (default: data/models/b32_400m.pt)
  --labels LABELS [LABELS ...]
                        Labels for classification (default: ["good cat", "bad cat"])
  --json                Output result as JSON to stdout
  --architecture {ViT-B-32-quickgelu,ViT-B-16-quickgelu,ViT-L-14-quickgelu,ViT-H-14-quickgelu}
                        Model architecture (auto-detected from filename if not specified)
```

### 7.2. Download the classification model

The script takes images, a model, and a set of labels and classifies each of the images according to the labels.
To run the script outside of Nextflow, we'll need to download one of the models.
Do so with:

```bash
mkdir -p data/models
(cd data/models && wget https://dl.fbaipublicfiles.com/MMPT/metaclip/b32_400m.pt)
```

### 7.3. Create the Classify process

Now let's create a `Classify` process that will take two channels - one channel of images and one channel that supplies the model:

```groovy title="Process definition" linenums="1"
process Classify {
    input:
        tuple val(meta), path(img)
        path(model)
    output: tuple val(meta), stdout
    script: "classify.py --model-path $model ${img}"
}
```

Note here that we're calling the `classify.py` script directly, even though we can't do that from the command line (we had to provide the relative or absolute path).
This is because Nextflow automatically adds the `bin` directory (relative to the main.nf) to the `$PATH` for all Nextflow tasks.
This is a very convenient way to bundle accessory scripts and snippets with your workflow.

### 7.4. Understanding queue vs. value channels

Processes can have multiple channels as input or as output.
A process will continue to emit tasks as long as it can pull an item from each of the input channels.
We could create a new channel for the model, and define a sensible default:

```groovy title="Workflow with model channel" linenums="1"
params.model = "${projectDir}/data/models/b32_400m.pt"

workflow {
    images = channel.fromPath("data/pics/*.{png,gif,jpg}")
    | map { img -> [[id: img.baseName], img] }

    model_channel = channel.fromPath(params.model)
    Classify(images, model_channel)
}
```

What happens when you run the workflow?
Given what we know about channels, what might be happening?

**Answer:** The Classify process only spawns a single task.
This is because after pulling the model path from the second input channel on the first iteration, the channel is empty, so no more Classify tasks can be submitted for execution.

There are two types of channel in Nextflow - **queue channels** and **value channels**.
Queue channels are exhaustible - they have a set number of items in the channel and each process can only take each item in the channel once.
The second type of channel is a value channel, which is a channel of only a single item.
This item is emitted without exhaustion.

### 7.5. Using value channels

There are some operators which will always return a value channel.
Examples are `first`, `collect`, `count`, etc.

We could also create a value channel using the `channel.value` factory:

```groovy title="Using channel.value" linenums="1"
params.model = "${projectDir}/data/models/b32_400m.pt"

workflow {
    images = channel.fromPath("data/pics/*.{png,gif,jpg}")
    | map { img -> [[id: img.baseName], img] }

    model_channel = channel.value(file(params.model))
    Classify(images, model_channel)
}
```

Note here that we're wrapping the params.model value (a String) in the `file()` function, which turns an ordinary String into an object that Nextflow can use as a path.
We've not needed to use this until now because the `channel.fromPath` factory necessarily returns paths, so it automatically does this conversion for us.

### 7.6. Implicit value channels

An even simpler solution is to provide the path object directly when calling the process.
Any non-channel object will automatically be converted into a value channel for you:

```groovy title="main.nf" hl_lines="8" linenums="1"
#!/usr/bin/env nextflow

params.model = "${projectDir}/data/models/b32_400m.pt"

workflow {
    images = channel.fromPath("data/pics/*.{png,gif,jpg}")
    | map { img -> [[id: img.baseName], img] }

    Resize(images, params.width)

    Classify(images, file(params.model))
    | view
}
```

Add the Classify process definition to your workflow and run it:

```bash
nextflow run main.nf
```

You might find that the process errors out with a 137 exit code.
This generally means that we've run out of RAM because we're running too many of these classification jobs at the same time.
Let's talk about how we tell Nextflow that a particular process requires more resources.

### Takeaway

Understanding queue channels vs. value channels is crucial for controlling how data flows through multi-input processes.

### What's next?

Let's learn how to manage computational resources for our processes.

---

## 2. Resources

Our processes are currently composed of the `input:`, `output:`, and `script:` blocks.
In addition to these blocks, processes can use "process directives" which are optional annotations which modify the behaviour of the processes.
There are many directives ([documentation](https://www.nextflow.io/docs/latest/reference/process.html#directives)), but we can introduce the concept with two important process directives - `memory` and `cpus`.

### 8.1. Understanding executors

So far, we've been using the local executor to run Nextflow - running on the local machine.
There are many other executors targeting different backends, from HPC executors like SLURM and PBS to cloud executors like AWS Batch, Google Batch, and Azure Batch.
There are more than a dozen supported executors ([documentation](https://www.nextflow.io/docs/latest/executor.html)).

Each of these have a concept of the resources a particular task will require - resources such as cpus, memory, gpus, disk, etc.

### 8.2. Resource defaults and management

If not otherwise specified, the defaults are to request 1 cpu, 1 GB of RAM and 0 GPUs for each task.

When using the local executor, Nextflow scans the machine it is running on and determines how many cpus and how much RAM the system has.
It will ensure that (given the resources specified or defaults applied) the running tasks never exceed the available limits.
If the system has 16 GB of RAM, for example, and a particular process requires 6 GB of ram, Nextflow will ensure that _at most_ 2 of those tasks are running at any one time.
As a task finishes, Nextflow begins the next task in line.

### 8.3. Add resource directives

Update your Classify process to request more memory:

```groovy title="Process with memory directive" hl_lines="2" linenums="1"
process Classify {
    memory '13 GB'

    input:
        tuple val(meta), path(img)
        path(model)
    output: tuple val(meta), stdout
    script: "classify.py --model-path $model ${img}"
}
```

Now run the workflow again:

```bash
nextflow run main.nf
```

### Takeaway

Process directives like `memory` and `cpus` communicate resource requirements to Nextflow executors, enabling proper scheduling and preventing resource exhaustion.

### What's next?

Let's learn how to combine related data using the join and groupTuple operators.

---

## 3. Grouping

Now we want to combine our classification results with our resized images.
We can use the `join` operator, which finds pairs of items (one from each channel) that share a key.
By default, the `join` operator will use the first element of each item in the channel as the key.
In our case, that first item was the image metadata, which occupies the first position in both the Classify process output and the Resize process output.

### 9.1. Join classification results with images

Update your workflow to join the channels:

```groovy title="Workflow with join" hl_lines="12-14" linenums="1"
#!/usr/bin/env nextflow

params.model = "${projectDir}/data/models/b32_400m.pt"

workflow {
    images = channel.fromPath("data/pics/*.{png,gif,jpg}")
    | map { img -> [[id: img.baseName], img] }

    Resize(images, params.width)

    Classify(images, file(params.model))
    | join(Resize.out)
    | view
}
```

This produces a channel like:

```
[metadata, label, img]
[metadata, label, img]
[metadata, label, img]
[metadata, label, img]
```

### 9.2. Group items by label

In order to make a picture of just the good cats and a second picture of just the bad cats, we'll need to group the items in the channel based on the label.
We can do this with the `groupTuple` operator.
Normally the groupTuple expects that the grouping key will be the first element in each item in the channel.
In our case, it is the second item, i.e. index "1" if the first item is index "0".
To ask Nextflow to group on the item with index 1, we add a `by: 1` argument to the operator:

```groovy title="Workflow with grouping" hl_lines="13-15" linenums="1"
#!/usr/bin/env nextflow

params.model = "${projectDir}/data/models/b32_400m.pt"

workflow {
    images = channel.fromPath("data/pics/*.{png,gif,jpg}")
    | map { img -> [[id: img.baseName], img] }

    Resize(images, params.width)

    Classify(images, file(params.model))
    | join(Resize.out)
    | groupTuple(by: 1)
    | view
}
```

This produces a channel of the form:

```
[metadatas, label, images]
[metadatas, label, images]
```

### Takeaway

The `join` and `groupTuple` operators allow you to match related items and collect them by common attributes.

### What's next?

Let's create visual collages for each group of classified images.

---

## 4. Collage

Let's create a `Collage` process that takes this channel and produces a collage of all of the images for each label.
The script block here is a little involved, but it uses ImageMagick's montage command to arrange images into a grid.

### 10.1. Create the Collage process

```groovy title="Collage process" linenums="1"
process Collage {
    input: tuple val(metadatas), val(label), path("inputs/*.png")
    output: tuple val(label), path("collage.png")
    script:
    """
    magick montage inputs/* \\
        -geometry +10+10 \\
        -background black \\
        +polaroid \\
        -background '#ffbe76' \\
        collage_nolabel.png
    magick montage \\
        -pointsize 48 \\
        -label '$label' \\
        -geometry +0+0 \\
        -background "#f0932b" \\
        collage_nolabel.png collage.png
    """
}
```

### 10.2. Connect to the workflow

We can then hook this into our channel chain:

```groovy title="Workflow with Collage" hl_lines="13-15" linenums="1"
#!/usr/bin/env nextflow

params.model = "${projectDir}/data/models/b32_400m.pt"

workflow {
    images = channel.fromPath("data/pics/*.{png,gif,jpg}")
    | map { img -> [[id: img.baseName], img] }

    Resize(images, params.width)

    Classify(images, file(params.model))
    | join(Resize.out)
    | groupTuple(by: 1)
    | Collage
    | view
}
```

### 10.3. Optimize with resized images

Those collage tasks are taking a little too long, but that might be because we're collaging the original full-sized images and not our resized images.
Because the `images` channel and the output channel from the `Resize` process both have the same shape, we can simply replace them in the workflow:

```groovy title="Optimized workflow" hl_lines="12" linenums="1"
#!/usr/bin/env nextflow

params.model = "${projectDir}/data/models/b32_400m.pt"

workflow {
    images = channel.fromPath("data/pics/*.{png,gif,jpg}")
    | map { img -> [[id: img.baseName], img] }

    Resize(images, params.width)

    Classify(images, file(params.model))
    | join(Resize.out)
    | groupTuple(by: 1)
    | Collage
    | view
}
```

### 10.4. Combine all collages

For our final process, let's combine these two collages together into a single final image.
We'll create a process that takes a collection of images (we don't care what they are called) and produces a final `collage_all.png` image:

```groovy title="CombineImages process" linenums="1"
process CombineImages {
    input: path "in.*.png"
    output: path "collage_all.png"
    script:
    """
    magick montage \\
        -geometry +10+10 \\
        -quality 05 \\
        -background '#ffbe76' \\
        -border 5 \\
        -bordercolor '#f0932b' \\
        in.*.png \\
        collage_all.png
    """
}
```

### 10.5. Transform the channel

The channel coming from the Collage process looks like:

```
[label, collageImage]
[label, collageImage]
```

but we need it to look like:

```
[collageImage, collageImage]
```

So we'll drop the labels and collect all images:

```groovy title="Final workflow" hl_lines="14-18" linenums="1"
#!/usr/bin/env nextflow

params.model = "${projectDir}/data/models/b32_400m.pt"

workflow {
    images = channel.fromPath("data/pics/*.{png,gif,jpg}")
    | map { img -> [[id: img.baseName], img] }

    Resize(images, params.width)

    Classify(images, file(params.model))
    | join(Resize.out)
    | groupTuple(by: 1)
    | Collage
    | map { _label, img -> img }
    | collect
    | CombineImages
    | view
}
```

The `collect` operator takes all the items in a channel and then emits them as a single "wide" collection.

Run the complete workflow:

```bash
nextflow run main.nf
```

### 10.6. Scaling up without code changes

One of Nextflow's key strengths is automatic scalability.
Let's see this in action by adding more data to our analysis!

While your workflow is still running (or right after it completes), open a new terminal and add more cat images:

```bash
# Add 20 more cats to our dataset
.stuff/cat_me.sh --count 20 --prefix data/pics
```

This brings our total from 4 cats to 24 cats.
Now run the workflow again with `-resume`:

```bash
nextflow run main.nf -resume
```

Notice what happens in the output:

- Tasks for the original 4 images show as **[cached]** in gray
- Only the 20 new images are processed through Resize and Classify
- The groupTuple, Collage, and CombineImages steps run again (because their inputs changed)
- The final collage now includes all 24 cats

**You didn't change a single line of code** - the workflow automatically:

- Detected the new input files via the glob pattern `data/pics/*.{png,gif,jpg}`
- Processed only the new images that hadn't been seen before
- Reused cached results for the original 4 images
- Scaled the grouping and collage operations to handle more data

This is the power of Nextflow's declarative approach: you describe **what** you want to do, and Nextflow figures out **how** to do it efficiently, whether you have 4 files or 4,000 files.

!!! tip "Scalability in practice"

    This same pattern works at any scale:

    - **Local development**: Test with 4 samples
    - **Pilot study**: Scale to 24 samples with no code changes
    - **Production**: Process thousands of samples with the same workflow
    - **HPC/Cloud**: Nextflow automatically distributes tasks across available resources

### Takeaway

You can chain together multiple processes and operators to build sophisticated multi-step workflows that transform and aggregate data.
Nextflow automatically scales your workflow as your data grows, without requiring any code changes.

### What's next?

Now that we have a working workflow, let's learn how to publish the results in an organized way.
