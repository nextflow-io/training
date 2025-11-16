# Outline

What are the Nextflow concepts being introduced at each stage?

## Introduction

Nextflow offers you a way to iterate over a collection of files, so let's grab some files to iterate over. We're going to write a workflow which produces a gallery of cute and not-so-cute cats. First things we're going to need are some cats.
I've created a small script that pull from the Cat-As-A-Service API to give us some random cats.

```bash
bin/cat_me.sh --help
```

To start, lets grab 4 cats. By default, the script will save the images to `./data/pics`:

```bash
bin/cat_me.sh --count 4 --prefix data/pics
```

Now lets iterate over those images in Nextflow. To start, we'll just create a channel of those images. We're not gong to do anything with them, but to make sure that everything is working, we connect the channel to the `view` operator which takes the things in the channel (files in our case) and prints a String representation of those things to the command line:

```nextflow
workflow {
    channel.fromPath("data/pics/*.{png,gif,jpg}")
    | view
}
```

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
    script: "convert ${img} -resize 400x resized-${img}"
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
convert kitten.jpg -resize 400x resized-kitten.jpg
```

Now let's run our workflow! We'll iterate over all of the images in `data/pics` (relative to our current location) and produce a channel of resized pictures that we then pipe into the `view` operator to print the channel contents to stdout.

## Investigate .command files

TODO: Explain that each tasks is run in a separate directory.This is to ensure independence of each of the tasks - they can't interfere with each other.

TODO: Show the contents of one of the task work directories.

## Harmonization
One of the nice features of the `convert` utility is that it will also do file format conversion for us. It will infer the format from the extension of the final argument. For example, if we execute

```bash
convert kitten.jpg -resize 400x resized-kitten.png
```

The `convert` utility will both resize the image and convert the jpg to png format. Let's say we want to ensure that downstream in our workflow, we'd like to ensure all images are in the png format. How might we modify our `script` block to replace the extension or pull out the file basename so that we can append the `.png` extension?

If you're a bash wizard, you might know that if you have a variable `$myFile` with the path to our file, you can replace the extension with this arcane incantation:

```bash
file=kitten.jpg
convert "$file" -resize 400x "${file%.*}.png"
```

Or perhaps you use the `basename` utility:

```bash
convert "$file" -resize 400x "$(basename "$file" .${file##*.}).png"
```

I can never remember these syntax. Forgunately for me, when inside the script block the `img` variable is not a bash variable - it's a Nextflow variable, and Nextflow provides some convenience methods for operating on those path objects.
We can simply call `${img.baseName}` to retun the file base name. For example:

```nextflow
process Resize {
    input: path(img)
    output: path("resized-*")
    script: "convert ${img} -resize 400x resized-${img.baseName}.png"
}
```

## Parameters
What if we want to make our workflow a little more flexible. Let's pull out the width and expose it as a parameter to the user.

```nextflow
process Resize {
    input: path(img)
    output: path("resized-*")
    script: "convert $img -resize ${width}x resized-${img.baseName}.png"
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
    script: "convert $img -resize ${width}x resized-${img.baseName}.png"
}
```

The params object still works in the same way:

```bash
nextflow run . --width 500
```

## Extracting an ID
Great, but I'd like a way of retaining the original IDs.

TODO: Explain the map operator and explain closures

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
    output: path("resized/*")
    script: "convert $img -resize ${width}x resized-${img.baseName}.png"
}
```

This is better, but now that we have that critical metadata being passed through the channels in our "meta" object, we can stop encoding the id in the filenames. We can also just pass the "meta" object through to the output channel so that the newly resized image stays locked with the meta.

```nextflow
process Resize {
    input:
        tuple val(meta), path(img)
        val(width)
    output: tuple val(meta), path("resized.png")
    script: "convert $img -resize ${width}x resized.png"
}
```

Run the workflow and view the output:

```bash
nextflow run .
```

## Classification

We have a little classification script - bin/classify.py, which you can run with:

```bash
bin/classify.py --help
```

```
usage: classify.py [-h] [--image-dir IMAGE_DIR] [--model-path MODEL_PATH] [--labels LABELS [LABELS ...]] [--json] [--architecture {ViT-B-32-quickgelu,ViT-B-16-quickgelu,ViT-L-14-quickgelu,ViT-H-14-quickgelu}]

Classify images using MetaCLIP

options:
  -h, --help            show this help message and exit
  --image-dir IMAGE_DIR
                        Directory containing images to classify (default: data/pics)
  --model-path MODEL_PATH
                        Path to MetaCLIP model weights (default: data/models/b32_400m.pt)
  --labels LABELS [LABELS ...]
                        Labels for classification (default: ["cute cat", "ugly cat"])
  --json                Output results as JSONL (one JSON object per line) to stdout
  --architecture {ViT-B-32-quickgelu,ViT-B-16-quickgelu,ViT-L-14-quickgelu,ViT-H-14-quickgelu}
                        Model architecture (auto-detected from filename if not specified)
```

The script takes images, a model, and a set of labels and classifies each of the images according to the labels. To run the script outside of Nextflow, we'll need to download one of the models. Do so with:

```bash
cd data/models
wget https://dl.fbaipublicfiles.com/MMPT/metaclip/b32_400m.pt
cd ../../
```

