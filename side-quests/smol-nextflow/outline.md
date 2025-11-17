# Outline

## Introduction

Nextflow offers you a way to iterate over a collection of files, so let's grab some files to iterate over. We're going to write a workflow which produces a gallery of good and bad cats. First things we're going to need are some cats.

We're starting with a (nearly) empty directory. There is a `.stuff` directory that contains some bits and pieces to help us along during the workshop, but you can imagine that we're essentially starting from scratch.

The first thing we're going to need is some data. I've created a small script that pull from the Cat-As-A-Service API to give us some random cats.

```bash
.stuff/cat_me.sh --help
```

To start, lets grab 4 cats. By default, the script will save the images to `./data/pics`:

```bash
.stuff/cat_me.sh --count 4 --prefix data/pics
```

This will generate some example data. It will look something like this:

```
data
└── pics
    ├── 5n4MTAC6ld0bVeCe.jpg
    ├── 5n4MTAC6ld0bVeCe.txt
    ├── IOSNx33kkgkiPfaP.jpg
    ├── IOSNx33kkgkiPfaP.txt
    ├── IRuDgdPJZFA39dyf.jpg
    ├── IRuDgdPJZFA39dyf.txt
    ├── uq5KqqiF0qpgTQVA.jpg
    └── uq5KqqiF0qpgTQVA.txt
```

Now lets iterate over those images in Nextflow. To start, we'll just create a channel of those images. We're not gong to do anything with them, but to make sure that everything is working, we connect the channel to the `view` operator which takes the things in the channel (files in our case) and prints a String representation of those things to the command line:

```nextflow
workflow {
    channel.fromPath("data/pics/*.{png,gif,jpg}")
    | view
}
```

TODO: Claude: Brief explanation what the {png,gif,jpg} syntax is doing (multiple glob pattern)

## Channels

Let's try actually doing something with those images. We'll start with something simple - resizing images. We'll need to download some software to do so. Eventually, we'll talk about containers and reproducibility, but just to start, let's download the software to our machine:

```bash
sudo apt update
sudo apt-get install -Y imagemagick
```

We'll create a new "process" for our resize operation. You can think of these processes as templates, or a classes. We'll connect the process to a channel and fire of a new "task" for each thing in the channel.

In our Resize process definition (indeed in most process definitions), there are three blocks - input, output, and script. We connect these processes by channels and these blocks describe what we expect to get, what we expect to emit, and the work we want to do in-between.

```nextflow
workflow {
    images = channel.fromPath("data/pics/*.{png,gif,jpg}")

    images
    | Resize
    | view
}

process Resize {
    input: path(img)
    output: path("resized-*")
    script: "magick ${img} -resize ${width}x resized-${img.baseName}.png"
}
```

### Input

The input block describes what we expect to take from the input channel. The "things" in the channel can be have a type. The most common "types" are

- `val` (values Strings, integers, Maps, complex obects),
- `path` (paths like directories, files), or
- `tuple` (a collection of values and paths).

We'll see more of those later, but for the moment, the channel we created in the `workflow` block is a channel of files, so in our process definition, we'll say "I'm gong to supply a channel of paths to this process, and as our process takes thinsg from that channel to spawn a new process, we'll call the thing in channel `img`.

### Output

The output block describes what we want to emit into the process' output channel. Again, we can describe the "type" of thing emitted - `val`, `path`, `tuple` and others. For now, we'll promise to produce a file (or directory) that matches the glob pattern `reseized-*`.

### Script

The script block describes what work we want to do on each of the things in the channel - how we're going to transform each of the things we pull from the input channel into the files or values we promised to emit into the outupt channel.

By default, the script block will be rendered into a bash script, but you can use any interpreted language that makes sese to you - python, ruby, R, zsh, closure, whatever. In this introductory workshop, we'll stick with the default bash.

We run the "convert" command from imagemagick which performs many types of manipulation. In our case, we'll use the `-resize` argument to resize the image to a width of 400 pixels. We also supply an output filename. You'll notice that we use the `${img}` varliable twice in our script block. This `${img}` is the varibale we defined in the input block. For each iteration of our process (each task), the variable will be the path to our individual image.

For example, if the "thing" in the channel is the image `kitten.jpg`, then when Nextflow creates a new Resize task for this file, it will "render" our script block into bash, replacing the `${img}` variables with the path to produce this valid bash:

```bash
magick kitten.jpg -resize 400x resized-kitten.jpg
```

Now let's run our workflow! We'll iterate over all of the images in `data/pics` (relative to our current location) and produce a channel of resized pictures that we then pipe into the `view` operator to print the channel contents to stdout.

## Investigate .command files

TODO: Explain that each tasks is run in a separate directory.This is to ensure independence of each of the tasks - they can't interfere with each other.

TODO: Show the contents of one of the task work directories.

## Harmonization

One of the nice features of the `convert` utility is that it will also do file format conversion for us. It will infer the format from the extension of the final argument. For example, if we execute

```bash
magick kitten.jpg -resize 400x resized-kitten.png
```

The `convert` utility will both resize the image and convert the jpg to png format. Let's say we want to ensure that downstream in our workflow, we'd like to ensure all images are in the png format. How might we modify our `script` block to replace the extension or pull out the file basename so that we can append the `.png` extension?

If you're a bash wizard, you might know that if you have a variable `$myFile` with the path to our file, you can replace the extension with this arcane incantation:

```bash
file=kitten.jpg
magick "$file" -resize 400x "${file%.*}.png"
```

Or perhaps you use the `basename` utility:

```bash
magick "$file" -resize 400x "$(basename "$file" .${file##*.}).png"
```

I love bash, but it's easy to forget this syntax or mistype it. Fortunately for us, when inside the script block the `img` variable is not a bash variable - it's a Nextflow variable, and Nextflow provides some convenience methods for operating on those path objects. The full list is available in the [Nextflow stdlib documentation](https://www.nextflow.io/docs/latest/reference/stdlib-types.html#stdlib-types-path), but one handy method is `baseName`
We can simply call `${img.baseName}` to retun the file base name. For example:

```nextflow
process Resize {
    input: path(img)
    output: path("resized-*")
    script: "magick ${img} -resize 400x resized-${img.baseName}.png"
}
```

## Parameters

What if we want to make our workflow a little more flexible. Let's pull out the width and expose it as a parameter to the user.

```nextflow
process Resize {
    input: path(img)
    output: path("resized-*")
    script: "magick $img -resize ${width}x resized-${img.baseName}.png"
}
```

Now we can run wth

```bash
nextflow run . --width 300
```

This is great, but it's considered best practice (and we'll see why in a bit) to make the inputs to a process explicit. We can do this by adding a second channel as input:

```nextflow
workflow {
    images = channel.fromPath("data/pics/*.{png,gif,jpg}")

    Resize(images, params.width)
    | view
}

process Resize {
    input:
      path(img)
      val(width)
    output: path("resized-*")
    script: "magick $img -resize ${width}x resized-${img.baseName}.png"
}
```

The params object still works in the same way:

```bash
nextflow run . --width 500
```

## Extracting an ID

Great, but I'd like a way of retaining the original IDs.

TODO: Claude: Explain the `map` operator and explain closures (3 sentences)

```nextflow
workflow {
    images = channel.fromPath("data/pics/*.{png,gif,jpg}")
    | map { img -> [img.baseName, img] }
    | view
}
```

Better - we have the id extracted as a String. What if we want to add other metadata later? Let's turn it into a Map

```nextflow
workflow {
    images = channel.fromPath("data/pics/*.{png,gif,jpg}")
    | map { img -> [[id: img.baseName], img] }
}
```

Now we've change the "shape" of the items in the channel, so we'll update the downstream process:

```nextflow
process Resize {
    input:
        tuple val(meta), path(img)
        val(width)
    output: path("resized-*")
    script: "magick $img -resize ${width}x resized-${img.baseName}.png"
}
```

Run the workflow and view the output:

```bash
nextflow run .
```

## Classification

Let's get the fun part - the cat sorting. We have a little classification script - `classify.py` that I've provided in the `.stuff` directory. In your research sometimes you have small accessory scripts that are useful for your pipeliens. We're using a python script here in this workshop example, but this patttern will hold for scripts written in perl, ruby, R, python, closurescript, or any of the other interpreted languages.

Let's pull the file out into a new `bin` directory:

```
mkdir -p bin
cp .stuff/classify.py bin/
```

The script requires some dependencies. Again, we'll do this the slow/painful way one time before we demonstrate how to use containers to encapsulate the software dependencies.

We'll grab one more file from our `.stuff` directory - a pyproject.toml file which is a way of describing softare dependencies for Python projects. This is unrelated to Nextflow, but an example of one of the (many) ways in which different languages and frameworks might install software.

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

````
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
                        Model architecture (auto-detected from filename if not specified)```
````

The script takes images, a model, and a set of labels and classifies each of the images according to the labels. To run the script outside of Nextflow, we'll need to download one of the models. Do so with:

```bash
mkdir -p data/models
(cd data/models && wget https://dl.fbaipublicfiles.com/MMPT/metaclip/b32_400m.pt)
```

Now let's create a `Classify` process that will take two channels - one channel of images and one channel that supplies the model:

```nextflow
process Classify {
    input:
        tuple val(meta), path(img)
        path(model)
    output: tuple val(meta), stdout
    script: "classify.py --model-path $model ${img}"
}
```

Note here that we're calling the `classify.py` script directly, even though we can't do that from the command line (we had to provide the relative or absolute path). This is because Nextflow automatically adds the `bin` directory (relative to the main.nf) to the `$PATH` for all Nextflow tasks. This is a very convenient way to bundle accessory scripts and snippets with your workflow.

Processes can have multiple channels as input or as output. A process will continue to emit tasks as long as it can pull an item from each of the input channels. We could create a new channel for the model, and define a sensible default:

```nextflow
params.model = "${projectDir}/data/models/b32_400m.pt"

workflow {
    images = channel.fromPath("data/pics/*.{png,gif,jpg}")
    | map { img -> [[id: img.baseName], img] }

    model_channel = channel.fromPath(params.model)
    // rest of the workflow
}
```

... which would return a channel with a single item. Try supplying this channel as input to our Classify process:

```nextflow
params.model = "${projectDir}/data/models/b32_400m.pt"

workflow {
    images = channel.fromPath("data/pics/*.{png,gif,jpg}")
    | map { img -> [[id: img.baseName], img] }

    model_channel = channel.fromPath(params.model)
    Classify(images, model_channel)
    // rest of the workflow
}
```

What happens when you run the workflow? Given what we know about the channels, what might be happening?

Answers: The Classify process only spawns a single task. This is because after pulling the model path from the second input channel on the first iteration, the channel is empty, so no more Classify tasks can be submitted for execution.

There are two types of channel in Nextflow - queue channels and value channels. Queue channels are exhaustible - they have a set number of items in the channel and each processes can only take each item in the channel once. The second type of channel is a value channel, which is a channel of only a single item. This item is emiited without exhaustion.

There are some operators which will alwys return a value channel. Examples are `first`, `collect`, `count`, etc.

We could also create a value channel using the `channel.value` factory:

```nextflow
params.model = "${projectDir}/data/models/b32_400m.pt"

workflow {
    images = channel.fromPath("data/pics/*.{png,gif,jpg}")
    | map { img -> [[id: img.baseName], img] }

    model_channel = channel.value(file(params.model))
    Classify(images, model_channel)
    // rest of the workflow
}
```

Note here that we're wrapping the params.model value (a String) in the `file()` function, which turns an ordinary String into an object that Nextflow can use as a path. We've not needed to use this until now because the `channel.fromPath` factory necessarily returns paths, so it automatically does this conversion for us.

An even simpler solution is to provide the path object directly when calling the process. Any non-channel object will automatically be converted into a value channel for you.

```nextflow
params.model = "${projectDir}/data/models/b32_400m.pt"

workflow {
    images = channel.fromPath("data/pics/*.{png,gif,jpg}")
    | map { img -> [[id: img.baseName], img] }

    Classify(images, file(params.model))
    // rest of the workflow
}
```

Now we'd like to take our channel of images and pass them each through this classification program.

We'd like to combine each of our images with each of our models (only one at the moment). Do do these cross-product operations, we can use the `combine` operator:

```nextflow
images
| combine(models)
| view
```

This will produce a channel that "combines" our images and the reference file. The channel now looks something like:

```
[imageMetaData, image, modelMetaData, model]
[imageMetaData, image, modelMetaData, model]
[imageMetaData, image, modelMetaData, model]
[imageMetaData, image, modelMetaData, model]
```

We might as well combine these two metadata Maps with a `map` operation:

```nextflow
images
| combine(models)
| map { imgMeta, img, modelMeta, model -> [imgMeta + modelMeta, model, img]}
| view
```

That now gives us a channel that looks like:

```
[metaData, image, model]
[metaData, image, model]
[metaData, image, model]
[metaData, image, model]
```

Given that shape, let's create a `Classify` process:

```nextflow
process Classify {
    input:
        tuple val(meta), path(img)
        path(model)
    output: tuple val(meta), stdout
    script: "classify.py --model-path $model ${img}"
}
```

Now we can run this:

```bash
nextflow run .
```

You might find that the process errors out with a 137 exit code. This generally means that we've run out of RAM because we're running too many of these classification jobs at the same time. Let's talk about how we tell Nextflow that a particular process requires more resources.

## Resources

Our processes are currently composed of the `input:`, `output:`, and `script:` blocks. In addition to these blocks, processes can use "process directives" which are optional annotations which modify the behaviour of the processes. There are many directives ([documentation](https://www.nextflow.io/docs/latest/reference/process.html#directives)), but we can introduce the concept with two important process directives - `memory` and `cpus`.

So far, we've been using the local executor to run Nextflow - running on the local machine. There are many other executors targetting different backends, from HPC executors like SLURM and PBS to cloud executors like {AWS,GCP,Azure} Batch. There are more than a dozen supported executors ([documentation](https://www.nextflow.io/docs/latest/executor.html)).

Each of these have a concept of the resources a particular task will require - resources such as cpus, memory, gpus, disk, etc.

If not otherwise specified, the defaults are to request 1 cpu, 1 GB of RAM and 0 GPUs for each task.

When using the local executor, Nextflow scans the machine it is running on and determines how many cpus and how much RAM the system has. It will ensure that (given the resources specified or defaults applied) the running tasks never exceed the available limits. If the system has 16 GB of RAM, for example, and a particular process requires 6 GB of ram, Nextflow will ensure that _at most_ 2 of those tasks are running at any one time. As a task finishes, Nextflow begins the next task in line.

## Grouping

Which we can now join up to our channel chain using the `join` operator, which finds pairs of items (one from each channel) that share a key. By default, the `join` operator will use the first element of each item in the channel as the key. In our case, that first item was the image metadata, which occupies the first position in both the Classify process output and the images channel.

```nextflow
Classify.out
| join(images)
| view
```

which produces a channel like:

```
[metadata, label, img]
[metadata, label, img]
[metadata, label, img]
[metadata, label, img]
```

In order to make a picture of just the good cats and a second picture of just the bad cats, we'll need to group the items in the channel based on the label. We can do this with the `groupTuple` operator. Normally the groupTuple expects that the grouping key will be the first element in each item in the channel. In our case, it is the second item, i.e. index "1" if the first item is index "0". To ask Nextflow to group on the item with index 1, we add a `by: 1` argument to the operator:

```nextflow
Classify.out
| join(images)
| groupTuple(by: 1)
| view
```

which produces a channenl of the form:

```
[metadatas, label, images]
[metadatas, label, images]
```

## Collage

Let's create a `collage` process that takes this channel and produces a collage of all of the images for each label. The script block here is a little involved.

```nextflow
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

We can then hoook this into our channel chain:

```nextflow
Classify.out
| join(images)
| groupTuple(by: 1)
| Collage
| view
```

Those collage tasks are taking a little too long, but that might be because we're collaging the original full sized images and not our resized images. Because the `images` channel and the output channel from the `Resize` process both have the same shape, we can simply replace them in the workflow:

```nextflow
Classify.out
| join(Resize.out)
| groupTuple(by: 1)
| Collage
| view
```

For our final process, let's combine these two collages together into a single final image. We'll create a process that takes collection of images (we don't care what they are called) and produces a final `collage_all.png` image.

```nextflow
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

The channel coming from the Collage process looks like

```
[label, collageImage]
[label, collageImage]
```

but we need it to look like:

```
[collageImage, collageImage]
```

So we'll drop the labels:

```
Classify.out
| join(Resize.out)
| groupTuple(by: 1)
| Collage
| map { _label, img -> img }
| view
```

to give us a channel that looks like:

```
collageImage
collageImage
```

and then pass that to `collect` which takes all the items in a channel and then emits them as a single "wide" collection:

```
Classify.out
| join(Resize.out)
| groupTuple(by: 1)
| Collage
| map { _label, img -> img }
| collect
| view
```

We can now pass this to our new CombineImages process:

```
Classify.out
| join(Resize.out)
| groupTuple(by: 1)
| Collage
| map { _label, img -> img }
| collect
| CombineImages
| view
```

## Workflow Outputs

Great! We have a workflow that (arguably cruely) collects our cat into "good" and "bad" groupings! Unfortunately, the final output file is still deep in this work directory in a hostile-looking hash-addressed directory. We'd like to define some final workflow outputs that should be published somewhere safe, outside of the work directory.

To define the workflow outputs, we'll need to define a `publish:` block in the workflow. We'll also need to put the existing workflow in a `main:` block as shown below:

```nextflow
workflow {
    main:
        images = channel.fromPath("data/pics/*.{png,gif,jpg}")
        | map { img -> [[id: img.baseName], img] }

        Resize(images, params.width)

        Classify(images, file(params.model))

        Classify.out
        | join(Resize.out)
        | groupTuple(by: 1)
        | Collage
        | map { _label, img -> img }
        | collect
        | CombineImages

    publish:
        // Something to go here
}
```

### Publishing the Collage

In the `publish:` block, we define channels that we'd like to publish.

```nextflow
workflow {
    main:
        // workflow here
    publish:
        collage = CombineImages.out
}

output {
    collage {}
}
```

Now when we run, the final collage will be symlinked into `results/collage_all.png`

We can control the publication mechanism my adding arguments:

```nextflow
output {
    collage {
        mode 'copy'
    }
}
```

... which will now cause Nextflow to copy the output file rather than symlink

### Publishing the classifications

The more intereting outputs might be those with more metadata associated with them. For example, we might want to record the classification for each image ID. To publsh metadata-rich outputs, we'll first create a channel that is composed of Maps, e.g.

```nextflow
workflow {
    main:
        // workflow here

        Classify.out
        | join(Resize.out)
        | map { meta, label, image -> meta + [label:label, image:image] }
        | set { classifiedMaps }

    publish:
        collage = CombineImages.out
        classification = classifiedMaps
}

output {
    collage {
        mode 'copy'
    }
    classified {
        mode 'copy'
    }
}
```

This will cause the resized images to also be published in the `results` directory, but it's looking a bit cluttered now.

```
results
├── collage_all.png
├── resized-4skdDxHm4yDsSJIr.png
├── resized-4y6Hyu0uzVZcEx89.png
├── resized-6Nb0ipGrHDHqCEmZ.png
└── resized-wfMCf1lHc9YPw455.png
```

Let's bring a little bit of order:

```nextflow
output {
    collage {
        mode 'copy'
    }
    classified {
        mode 'copy'
        path { sample -> "images/${sample.label.replaceAll(/\s+/, '_')}" }
    }
}
```

```
results
├── collage_all.png
└── images
    ├── good_cat
    │   ├── resized-4skdDxHm4yDsSJIr.png
    │   └── resized-wfDuHvt6VIn2tM8T.png
    └── bad_cat
        ├── resized-6Nb0ipGrHDHqCEmZ.png
        └── resized-wfMCf1lHc9YPw455.png
```

Now we sanitized the label names so that they'd be in more sensibly named directories (no spaces, etc), but this risks corrupting that metadata. Lets ask Nextflow to publish a more digestible samplesheet or "index" of the published outputs, that includes the real, unsanitized labels:

```nextflow
output {
    collage {
        mode 'copy'
    }
    classified {
        mode 'copy'
        path { sample -> "images/${sample.label.replaceAll(/\s+/, '_')}" }
        index {
            header true
            path 'images/cats.csv'
        }
    }
}
```

which produces a csv at results/images/cats.csv.

For more structured data, you can also choose yaml or json:

```nextflow
output {
    collage {
        mode 'copy'
    }
    classified {
        mode 'copy'
        path { sample -> "images/${sample.label.replaceAll(/\s+/, '_')}" }
        index {
            header true
            path 'images/cats.json'
        }
    }
}
```

## Filesystem Independence

Nextflow speasks many different communication protocols, allowing you to seamlessly move from using data on a local or shared filesystem, to http/https://, to object storage protocols like s3://, az://, gcp:// or even older ftp:// protocols. You can provide support for new protocols yourself via Nextflow's plugin system.

For example, our current workflow uses the

```bash
nextflow run . --model https://dl.fbaipublicfiles.com/MMPT/metaclip/b32_400m.pt
```

## Containerization

All of our Nextflow tasks are currently using the software installed on the host operating system. This practice can quickly become a problem for you for a number of reasons:

- As the workflow grows, and the number of software dependency stacks will also likely grow, and it becomes increasingly likely that the installation of one piece of software accidentally updates a depencency of another piece of software. Similarly, you many end up with incompatible software dependency stacks.
- The analysis becomes tied to a very specific machine or infrastructure, difficult to repdocuce in exactly the same way elsewhere (by yourself or by a colleague)
- Managing software is a thankless and boring task.

Nextflow provides the opportunity to run each task in an isolated software environment, and can do so via a variety of technologies, including

- conda
- containers (docker, apptainer/singularity, charliecloud, sarus, shifter, and podman)
- spack

Let's improve the reproducibility and portability of our workflow

TODO: Claude: include a brief (3-4 sentence) explanation of containerization (with a focus on Docker)

You'll remember that we manually installed software two different ways:

- imagemagick (via `apt-get install`), and
- python packages (via `uv sync`)

We could use a single container for all of the steps in the workflow, but this might limit the reusability of the containers, and upgrading one piece of software for one task would mean changing the container for all of the tasks. Most researchers prefer (and Nextflow supports) defining container per-process.

To replace the imagemagick we installed via apt-get, we'll use the public container 'minidocks/imagemagick:7'

We've already talked about the `memory` and `cpus` process directives, but another useful directive is the `container` directive. We'll use this to add the container to our `Resize`, `Classify`, `Collage`, and `CombineImages` processes, e.g.

```nextflow
process Resize {
    container 'minidocks/imagemagick:7'
    //...

process Classify {
    container 'minidocks/imagemagick:7'
    //...

process Collage {
    container 'minidocks/imagemagick:7'
    //...

process CombineImages {
    container 'minidocks/imagemagick:7'
    //...
```

Our `classify.py` process includes three specific python packages (torch, pillow, and openclip-torch) at specific versions. It's unlikely that there is an existing container that provides these specific packages. We could opt to build our own

There are a number of ways of building containers, but we'll use the [Seqera Containers](https://seqera.io/containers/) web interface. You can add multiple packages

![Creating a new container using Seqera Containers](./img/seqera-container-python-00.png)

## Version Control

TOOD: Create git repository at project root
TODO: Commit current state, create git tag, and create branch, and then change and re-commit.
TODO: Change directories and then run using revision argument, pointing to branch, tag, and then specific commmit.

## Cloud Executors

## Extension Exercise 1

Our team is interested in which cat is the custest cat and which cat is the ugliest cat. Can you extend the workflow to identify (for each label) which picture scores the highest?

Hint 1: You can use the --json flat on the `classify.py` script
Hint 2: You can parse a json file in a closure by using the JsonSlurper class, part of the standard library. It will return a standard

```nextflow
| map { meta, jsonFile -> new groovy.json.JsonSlurper().parseText(jsonFile.text) }
```

Hint 3: You can use the `min` and `max` operators to return a channel containing the minimum or maximum item, and you can pass a closure to those operators to describe how the elements in the channel should be compared ([docs](https://www.nextflow.io/docs/latest/reference/operator.html#min))

## Extension Exercise 2

We've decided that "bad" and "good" are too cruel a classification system for the cats. Can you modify the workflow to add a `--labels` parameter. The parameter should take a comma-separated list of labels and use those labels in preference to the default "good cat" and "bad cat". E.g.

```
nextflow run . --labels 'red cat','orange cat','black cat'
```

## Last TODOs

TODO: explain difference between path(img) and path("inputs/\*.png")
TODO: Add in resouces directive memory, cpus, etc.

<!--
This block is not part of the cource material and is safe to ignore
# Notes on wrapping and unwrapping metadata
I've heard many variants of a common request which might be phrased as "I'd like to call a subworkflow on a channel of "things".
Let's say we have a channel of samples and a second channel of references. The subworkflow expects a channel of samples, but only a single reference. The subworkflow produces a single vcf file. -->
