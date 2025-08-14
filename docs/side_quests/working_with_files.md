# Working with Files

Bioinformatics workflows often involve processing large numbers of files. Nextflow provides powerful tools to handle files efficiently, helping you organize and process your data with minimal code.

In this side quest, we'll explore how Nextflow handles files, from basic file operations to more advanced techniques for working with file collections. You'll learn how to extract metadata from filenames - a common requirement in bioinformatics pipelines.

By the end of this side quest, you'll be able to:

- Create Path objects from file path strings using Nextflow's `file()` method
- Access file attributes such as name, extension, and parent directory
- Handle both local and remote files transparently using URIs
- Use channels to automate file handling with `Channel.fromPath()` and `Channel.fromFilePairs()`
- Extract and structure metadata from filenames using string manipulation
- Group related files using pattern matching and glob expressions
- Integrate file operations into Nextflow processes with proper input handling
- Organize process outputs using metadata-driven directory structures

---

## 0. Warmup

### 0.1. Prerequisites

Before taking on this side quest you should:

- Complete the [Hello Nextflow](../hello_nextflow/README.md) tutorial
- Understand basic Nextflow concepts (processes, channels, operators)
- Basic understanding of using sample-specific data (meta data): [Working with sample-specific data](./metadata.md)

### 0.2. Starting Point

Let's move into the project directory:

```bash title="Navigate to project directory"
cd side-quests/working_with_files
```

You can set VSCode to focus on this directory:

```bash title="Open directory in VSCode"
code .
```

You'll find a simple workflow file (`file_operations.nf`) and a data directory containing some example files.

```console title="Directory contents"
> tree
.
├── count_lines.nf
├── data
│   ├── patientA_rep1_normal_R1_001.fastq.gz
│   ├── patientA_rep1_normal_R2_001.fastq.gz
│   ├── patientA_rep1_tumor_R1_001.fastq.gz
│   ├── patientA_rep1_tumor_R2_001.fastq.gz
│   ├── patientA_rep2_normal_R1_001.fastq.gz
│   ├── patientA_rep2_normal_R2_001.fastq.gz
│   ├── patientA_rep2_tumor_R1_001.fastq.gz
│   ├── patientA_rep2_tumor_R2_001.fastq.gz
│   ├── patientB_rep1_normal_R1_001.fastq.gz
│   ├── patientB_rep1_normal_R2_001.fastq.gz
│   ├── patientB_rep1_tumor_R1_001.fastq.gz
│   ├── patientB_rep1_tumor_R2_001.fastq.gz
│   ├── patientC_rep1_normal_R1_001.fastq.gz
│   ├── patientC_rep1_normal_R2_001.fastq.gz
│   ├── patientC_rep1_tumor_R1_001.fastq.gz
│   └── patientC_rep1_tumor_R2_001.fastq.gz
└── file_operations.nf

2 directories, 18 files
```

This directory contains paired-end sequencing data for three patients (A, B, C), with a typical `_R1_` and `_R2_` naming convention for what are known as 'forward' and 'reverse reads' (don't worry if you don't know what this means, it's not important for this session). Each patient has normal and tumor tissue types, and patient A has two replicates.

### 0.3. Running the Workflow

Take a look at the workflow file `file_operations.nf`:

```groovy title="file_operations.nf" linenums="1"
workflow {
    // Create a Path object from a string path
    myFile = 'data/patientA_rep1_normal_R1_001.fastq.gz'
    println "${myFile} is of class ${myFile.class}"
}
```

We have a mini-workflow that refers to a single file path in it's workflow, then prints it to the console, along with its class.

!!! note

    **What is `.class`?**

    In Groovy (the language Nextflow uses), `.class` tells us what type of object we're working with. It's like asking "what kind of thing is this?" - whether it's a string, a number, a file, or something else. This will help us see the difference between a plain string and a Path object in the next sections.

Run the workflow:

```bash title="Run the workflow"
nextflow run file_operations.nf
```

```console title="Starting Output"
 N E X T F L O W   ~  version 25.04.3

Launching `file_operations.nf` [romantic_chandrasekhar] DSL2 - revision: 5a4a89bc3a

data/patientA_rep1_normal_R1_001.fastq.gz is of class java.lang.String
```

We printed the string path exactly as we wrote it. This is just text output - Nextflow hasn't done anything special with it yet. We've also confirmed that, so far, to Nextflow this is only a string (of class `java.lang.String`), we haven't yet told Nextflow about its file nature.

## 1. Basic File Operations

### 1.1. Creating Path Objects

We can tell Nextflow how to handle files by creating [Path objects](https://www.nextflow.io/docs/latest/reference/stdlib-types.html#path) from path strings. In our workflow, we have a string path `data/patientA_rep1_normal_R1_001.fastq.gz`, and we covert that to a Path object using the `file()` method, which provides access to file properties and operations.

Edit the `file_operations.nf` to wrap the string with `file()` as follows:

=== "After"

    ```groovy title="file_operations.nf" linenums="2" hl_lines="2"
    // Create a Path object from a string path
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')
        println "${myFile} is of class ${myFile.class}"
    ```

=== "Before"

    ```groovy title="file_operations.nf" linenums="2" hl_lines="2"
    // Create a Path object from a string path
        myFile = 'data/patientA_rep1_normal_R1_001.fastq.gz'
        println "${myFile} is of class ${myFile.class}"
    ```

Run the workflow:

```bash title="Test Path object creation"
nextflow run file_operations.nf
```

```console title="Path object output"
 N E X T F L O W   ~  version 25.04.3

Launching `file_operations.nf` [kickass_coulomb] DSL2 - revision: 5af44b1b59

/workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz is of class class sun.nio.fs.UnixPath
```

Now we see the full absolute path instead of the relative path we wrote. Nextflow has converted our string into a Path object and resolved it to the actual file location on the system. The file path will now be absolute, like `/workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz`, and notice that the Path object class is `sun.nio.fs.UnixPath`, this is Nextflow's way of representing local files. As we'll see later, remote files will have different class names (like `nextflow.file.http.XPath` for HTTP files), but they all work exactly the same way and can be used identically in your workflows.

!!! note

    **The key difference:**

    - **Path string**: Just text that Nextflow treats as characters
    - **Path object**: A smart file reference that Nextflow can work with

    Think of it like this: a path string is like writing an address on paper, while a Path object is like having a GPS device that knows how to navigate and can tell you details about the journey.

### 1.2. File Attributes

Why is this helpful? Well now Nextflow understands that `myFile` is a Path object and not just a string, we can access the various attributes of the Path object.

Let's update our workflow to print out the file attributes:

=== "After"

    ```groovy title="file_operations.nf" linenums="2" hl_lines="5-9"
        // Create a Path object from a string path
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        // Print file attributes
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
    ```

=== "Before"

    ```groovy title="file_operations.nf" linenums="2" hl_lines="3"
        // Create a file object from a string path
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')
        println "${myFile} is of class ${myFile.class}"
    ```

Run the workflow:

```bash title="Test file attributes"
nextflow run file_operations.nf
```

You'll see various file attributes printed to the console:

```console title="File Attributes Output"
 N E X T F L O W   ~  version 25.04.3

Launching `file_operations.nf` [ecstatic_ampere] DSL2 - revision: f3fa3dcb48

File object class: sun.nio.fs.UnixPath
File name: patientA_rep1_normal_R1_001.fastq.gz
Simple name: patientA_rep1_normal_R1_001
Extension: gz
Parent directory: /workspaces/training/side-quests/working_with_files/data
```

### 1.3. Why proper file handling matters

The difference between strings and Path objects becomes critical when you start building actual workflows with processes. Let's side-step for a moment to take a look at a workflow where this has been done wrong.

`count_lines.nf` contains a process that takes a `val` input and tries to treat it as a file:

```groovy title="count_lines.nf" linenums="1" hl_lines="5 16"
process COUNT_LINES {
    debug true

    input:
    val fastq_file

    script:
    """
    set -o pipefail
    echo "Processing file: $fastq_file"
    gzip -dc $fastq_file | wc -l
    """
}

workflow {
    myFile = 'data/patientA_rep1_normal_R1_001.fastq.gz'
    COUNT_LINES(myFile)
}
```

!!! note

    **About `debug = true`:**

    The `debug = true` directive in the process definition causes Nextflow to print the output from your script (like the line count "40") directly in the execution log. Without this, you would only see the process execution status but not the actual output from your script.

    For more information on debugging Nextflow processes, see the [Debugging Nextflow Workflows](./debugging.md) side quest.

Run this workflow to see the error:

```bash title="Test val input with string error"
nextflow run count_lines.nf
```

You'll get an error like this:

```console title="Val input with string error"
 N E X T F L O W   ~  version 25.04.3

Launching `count_lines.nf` [goofy_koch] DSL2 - revision: 4d9e909d80

executor >  local (1)
[7f/c22b7f] COUNT_LINES [  0%] 0 of 1
ERROR ~ Error executing process > 'COUNT_LINES'

Caused by:
  Process `COUNT_LINES` terminated with an error exit status (1)


Command executed:

executor >  local (1)
[7f/c22b7f] COUNT_LINES [  0%] 0 of 1 ✘
WARN: Got an interrupted exception while taking agent result | java.lang.InterruptedException
ERROR ~ Error executing process > 'COUNT_LINES'

Caused by:
  Process `COUNT_LINES` terminated with an error exit status (1)


Command executed:

  set -o pipefail
  echo "Processing file: data/patientA_rep1_normal_R1_001.fastq.gz"
  gzip -dc data/patientA_rep1_normal_R1_001.fastq.gz | wc -l

Command exit status:
  1

Command output:
  Processing file: data/patientA_rep1_normal_R1_001.fastq.gz
  0

Command error:
  Processing file: data/patientA_rep1_normal_R1_001.fastq.gz
  gzip: data/patientA_rep1_normal_R1_001.fastq.gz: No such file or directory
  0

Work dir:
  /workspaces/training/side-quests/working_with_files/work/7f/c22b7f6f86c81f14d53de15584fdd5

Tip: you can replicate the issue by changing to the process work dir and entering the command `bash .command.run`

 -- Check '.nextflow.log' file for details
```

The process failed because the file `data/patientA_rep1_normal_R1_001.fastq.gz` doesn't exist in the process working directory. When you use `val` input, Nextflow passes the string value through to your script, but it doesn't stage the actual file. The process tries to use the string as a file path, but the file isn't there.

Now let's fix this by changing the process to use a `path` input:

=== "After"

    ```groovy title="file_operations.nf" linenums="1" hl_lines="5 16"
    process COUNT_LINES {
        debug true

        input:
        path fastq_file

        script:
        """
        set -o pipefail
        echo "Processing file: $fastq_file"
        gzip -dc $fastq_file | wc -l
        """
    }

    workflow {
        myFile = 'data/patientA_rep1_normal_R1_001.fastq.gz'
        COUNT_LINES(myFile)
    }
    ```

=== "Before"

    ```groovy title="file_operations.nf" linenums="1" hl_lines="5 16"
    process COUNT_LINES {
        debug true

        input:
        val fastq_file

        script:
        """
        set -o pipefail
        echo "Processing file: $fastq_file"
        gzip -dc $fastq_file | wc -l
        """
    }

    workflow {
        myFile = 'data/patientA_rep1_normal_R1_001.fastq.gz'
        COUNT_LINES(myFile)
    }
    ```

Run this updated version and you'll get a different error:

```console title="Path input with string error"
 N E X T F L O W   ~  version 25.04.3

Launching `file_operations.nf` [mighty_poitras] DSL2 - revision: e996edfc53

[-        ] COUNT_LINES -
ERROR ~ Error executing process > 'COUNT_LINES'

Caused by:
  Not a valid path value: 'data/patientA_rep1_normal_R1_001.fastq.gz'

Tip: when you have fixed the problem you can continue the execution adding the option `-resume` to the run command line
```

**What this error means:**

This is much better! Nextflow immediately detected the problem and failed before even starting the process. When you use `path` input, Nextflow validates that you're passing actual file references, not just strings. It's telling you that `'data/patientA_rep1_normal_R1_001.fastq.gz'` is not a valid path value because it's a string, not a Path object.

Now let's fix this properly by using the `file()` method to create a Path object:

=== "After"

    ```groovy title="file_operations.nf" linenums="1" hl_lines="5 16"
    process COUNT_LINES {
        debug true

        input:
        path fastq_file

        script:
        """
        set -o pipefail
        echo "Processing file: $fastq_file"
        gzip -dc $fastq_file | wc -l
        """
    }

    workflow {
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')
        COUNT_LINES(myFile)
    }
    ```

=== "Before"

    ```groovy title="file_operations.nf" linenums="1" hl_lines="5 16"
    process COUNT_LINES {
        debug true

        input:
        path fastq_file

        script:
        """
        set -o pipefail
        echo "Processing file: $fastq_file"
        gzip -dc $fastq_file | wc -l
        """
    }

    workflow {
        myFile = 'data/patientA_rep1_normal_R1_001.fastq.gz'
        COUNT_LINES(myFile)
    }
    ```

Now when you run this, it should work correctly! The file will be staged into the process working directory and the `wc -l` command will succeed.

```bash title="Test successful execution"
nextflow run count_lines.nf
```

You should see output like this:

```console title="Successful execution"
 N E X T F L O W   ~  version 25.04.3

Launching `count_lines.nf` [astonishing_tesla] DSL2 - revision: ee38f96485

executor >  local (1)
[8a/3f23d5] COUNT_LINES [100%] 1 of 1 ✔
Processing file: patientA_rep1_normal_R1_001.fastq.gz
      40
```

The process successfully:

- Staged the file into the working directory
- Decompressed the .gz file
- Counted the lines (40 lines in this case)
- Completed without errors

### Takeaway

- Path strings vs Path objects: Strings are just text, Path objects are smart file references
- The `file()` method converts a string path into a Path object that Nextflow can work with
- You can access file properties like `name`, `simpleName`, `extension`, and `parent` [using file attributes](https://www.nextflow.io/docs/latest/working-with-files.html#getting-file-attributes)
- Using Path objects instead of strings allows Nextflow to properly manage files in your workflow
- Process Input Outcomes: Proper file handling requires Path objects, not strings, to ensure files are correctly staged and accessible in processes.

---

## 2. Using Remote Files

One of the key features of Nextflow is the ability to transparently switch between local files (on the same machine) to remote files accessible over the internet. To do this, all you need to do as a user is switch the file path you supply to the workflow from a normal file path (e.g. `/path/to/data`) to a file path with a remote protocol at the start. Importantly, you should **never** have to change the workflow logic to accomodate files coming from different locations.

For example, replacing `/path/to/data` with `s3://path/to/data` in your inputs will switch to using the S3 protocol. Many different protocols are supported:

- HTTP(S)/FTP (http://, https://, ftp://)
- Amazon S3 (s3://)
- Azure Blob Storage (az://)
- Google Cloud Storage (gs://)

To use these, simply replace the string and Nextflow will handle authentication and staging the files to the right place, downloading or uploading and all other file operations you would expect. We call this string the Uniform Resource Identifier (URI).

The key strength of this is we can switch between environments without changing any pipeline logic. For example, you can develop with a small, local test set before moving to a remote set of data by changing the URI.

### 2.1. Using a file from the internet

!!! warning

    Accessing remote data requires an internet connection!

In your workflow, you can replace the string path with an HTTPS one to download this file from the internet. We are going to swap the relative path of the FASTQ files with the remote one. This is the same data as we have been previously using.

Open `file_operations.nf` again and make changes like this:

=== "After"

    ```groovy title="file_operations.nf" linenums="2" hl_lines="2"
    // Using a remote file from the internet
        myFile = file('https://github.com/nextflow-io/training/blob/bb187e3bfdf4eec2c53b3b08d2b60fdd7003b763/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz')

        // Print file attributes
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
    ```

=== "Before"

    ```groovy title="file_operations.nf" linenums="2" hl_lines="2"
        // Create a file object from a string path
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        // Print file attributes
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
    ```

!!! note

    HTTPS remote data does not accept globs because HTTPS cannot list multiple files, however other storage protocols such as blob storage (`s3://`, `az://`, `gs://`) can.

Run the workflow and it will automatically pull the data from the internet:

```bash title="Test remote file access"
nextflow run file_operations.nf
```

```console title="Remote file output"
 N E X T F L O W   ~  version 25.04.3

Launching `file_operations.nf` [insane_swartz] DSL2 - revision: fff18abe6d

File object class: class nextflow.file.http.XPath
File name: patientA_rep1_normal_R1_001.fastq.gz
Simple name: patientA_rep1_normal_R1_001
Extension: gz
Parent directory: /nextflow-io/training/blob/bb187e3bfdf4eec2c53b3b08d2b60fdd7003b763/side-quests/working_with_files/data
```

In this example very little has changed! This shows how easy it is to switch between local and remote data using Nextflow.

To see the difference, check the working directory located at the hash value of the process. The file will be downloaded to a staging directory located within the work directory, then symlinked into the relevant process directory. If the file was used again, Nextflow would only download it once.

In this way, you can replace local with remote data without changing any pipeline logic.

#### Cloud Storage with Glob Patterns

While HTTP doesn't support globs, cloud storage protocols do. Here's how you could use glob patterns with cloud storage:

```groovy title="Cloud storage examples (not runnable in this environment)"
// S3 with glob patterns - would match multiple files
ch_s3_files = Channel.fromPath('s3://my-bucket/data/*.fastq.gz')

// Azure Blob Storage with glob patterns
ch_azure_files = Channel.fromPath('az://container/data/patient*_R{1,2}.fastq.gz')

// Google Cloud Storage with glob patterns
ch_gcs_files = Channel.fromPath('gs://bucket/data/sample_*.fastq.gz')
```

These examples show the power of Nextflow's unified file handling - the same code works whether files are local or in the cloud, as long as the protocol supports the operations you need.

### Takeaway

- Remote data is accessed using a URI (HTTP, FTP, S3, Azure, Google Cloud)
- Nextflow will automatically download and stage the data to the right place
- Do not write logic to download or upload remote files!
- Local and remote files produce different object types but work identically
- **Important**: HTTP/HTTPS only work with single files (no glob patterns)
- Cloud storage (S3, Azure, GCS) supports both single files and glob patterns
- You can seamlessly switch between local and remote data sources without changing code logic (as long as the protocol supports your required operations)

!!! note

    **Note on Object Types**: Notice that local files produce `sun.nio.fs.UnixPath` objects while remote files produce `nextflow.file.http.XPath` objects. Despite these different class names, both work exactly the same way and can be used identically in your workflows. This is a key feature of Nextflow - you can seamlessly switch between local and remote data sources without changing your code logic.

### 2.2. Switching Back to Local Files

For the remainder of this side quest, we'll use local files in our examples. This allows us to demonstrate powerful features like glob patterns and batch processing that aren't available with HTTP URLs. Remember: the same concepts apply to cloud storage (S3, Azure, GCS) where glob patterns are fully supported.

Let's update our workflow to use local files again:

=== "After"

    ```groovy title="file_operations.nf" linenums="2" hl_lines="2"
        // Create a file object from a string path
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        // Print file attributes
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
    ```

=== "Before"

    ```groovy title="file_operations.nf" linenums="2" hl_lines="5-8"
        // Create a Path object from a string path
        myFile = file('https://github.com/nextflow-io/training/blob/bb187e3bfdf4eec2c53b3b08d2b60fdd7003b763/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz')

        // Print file attributes
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
    ```

---

## 3. Reading files using the `fromPath()` channel factory

The `file()` method is useful for simple file operations, and we can combine that with [`Channel.of()`](https://www.nextflow.io/docs/latest/reference/channel.html#of) to build channels from files like:

```groovy title="Channel.of() with file()"
    ch_fastq = Channel.of([file('data/patientA_rep1_normal_R1_001.fastq.gz')])
```

But we have a much more convenient tool called [`Channel.fromPath()`](https://www.nextflow.io/docs/latest/reference/channel.html#frompath) which generates a channel from static file strings as well as glob patterns.

### 3.1. Reading Files with Channel.fromPath

Update your `file_operations.nf` file:

=== "After"

    ```groovy title="file_operations.nf" linenums="2" hl_lines="2"
        // Reading files with Channel.fromPath
        ch_fastq = Channel.fromPath('data/patientA_rep1_normal_R1_001.fastq.gz')
        ch_fastq.view { "Found file: $it of type ${it.class}" }

        // // Print file attributes
        // Comment these out for now, we'll come back to them!
        // println "File object class: ${myFile.class}"
        // println "File name: ${myFile.name}"
        // println "Simple name: ${myFile.simpleName}"
        // println "Extension: ${myFile.extension}"
        // println "Parent directory: ${myFile.parent}"
    ```

=== "Before"

    ```groovy title="file_operations.nf" linenums="2" hl_lines="5-8"
    // Create a Path object from a string path
        // Create a file object from a string path
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        // Print file attributes
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
    ```

Run the workflow:

```bash title="Test Channel.fromPath"
nextflow run file_operations.nf
```

You'll see each file path being emitted as a separate element in the channel:

```console title="Channel.fromPath Output"
 N E X T F L O W   ~  version 25.04.3

Launching `file_operations.nf` [grave_meucci] DSL2 - revision: b09964a583

Found file: /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz of type class sun.nio.fs.UnixPath
```

Note how Nextflow has grabbed the file we specified and turned it into a `Path` type object, in exactly the same way that `file()` would have done. `Channel.fromPath()` is just a convenient way of creating a new channel populated by a list of files.

### 3.2. Viewing Channel Contents

In our first version, we use `.view()` to print the file name. Let's update our workflow to print out the file attributes:

=== "After"

    ```groovy title="file_operations.nf" linenums="2" hl_lines="3-9"
        // Reading files with Channel.fromPath
        ch_fastq = Channel.fromPath('data/patientA_rep1_normal_R1_001.fastq.gz')
        ch_fastq.view { myFile ->
            println "File object class: ${myFile.class}"
            println "File name: ${myFile.name}"
            println "Simple name: ${myFile.simpleName}"
            println "Extension: ${myFile.extension}"
            println "Parent directory: ${myFile.parent}"
        }
    ```

=== "Before"

    ```groovy title="file_operations.nf" linenums="2" hl_lines="3"
        // Reading files with Channel.fromPath
        ch_fastq = Channel.fromPath('data/patientA_rep1_normal_R1_001.fastq.gz')
        ch_fastq.view { myFile -> "Found file: $myFile" }

        // // Print file attributes
        // Comment these out for now, we'll come back to them!
        // println "File name: ${myFile.name}"
        // println "Simple name: ${myFile.simpleName}"
        // println "Extension: ${myFile.extension}"
        // println "Parent directory: ${myFile.parent}"
    ```

Run the workflow:

```bash title="Test file attributes with Channel.fromPath"
nextflow run file_operations.nf
```

```console title="Channel.fromPath Output"
 N E X T F L O W   ~  version 25.04.3

Launching `file_operations.nf` [furious_swanson] DSL2 - revision: c35c34950d

File object class: sun.nio.fs.UnixPath
File name: patientA_rep1_normal_R1_001.fastq.gz
Simple name: patientA_rep1_normal_R1_001
Extension: gz
Parent directory: /workspaces/training/side-quests/working_with_files/data
```

### 3.3. Using a glob to match multiple files

`Channel.fromPath()` can take a glob pattern as an argument, which will match all files in the directory that match the pattern. Let's grab both of the pair of FASTQs associated with this patient.

A glob pattern is a pattern that matches one or more characters in a string. The `*` wildcard is the most common glob pattern, which will match any character in it's place. To do this, we replace the full path with a `*` wildcard, which will match any character in it's place. In this case, we will replace the read number from `R1` to `R*`.

=== "After"

    ```groovy title="file_operations.nf" linenums="3"
        ch_fastq = Channel.fromPath('data/patientA_rep1_normal_R*_001.fastq.gz')
    ```

=== "Before"

    ```groovy title="file_operations.nf" linenums="3"
        ch_fastq = Channel.fromPath('data/patientA_rep1_normal_R1_001.fastq.gz')
    ```

Run the workflow:

```bash title="Test glob pattern matching"
nextflow run file_operations.nf
```

```console title="Channel.fromPath Glob Output"
 N E X T F L O W   ~  version 25.04.3

Launching `file_operations.nf` [boring_sammet] DSL2 - revision: d2aa789c9a

File object class: sun.nio.fs.UnixPath
File name: patientA_rep1_normal_R1_001.fastq.gz
Simple name: patientA_rep1_normal_R1_001
Extension: gz
Parent directory: /workspaces/training/side-quests/working_with_files/data
File object class: sun.nio.fs.UnixPath
File name: patientA_rep1_normal_R2_001.fastq.gz
Simple name: patientA_rep1_normal_R2_001
Extension: gz
Parent directory: /workspaces/training/side-quests/working_with_files/data
```

Using this method, we could grab as many or as few files as we want just by changing the glob pattern. If we made it more generous, we could grab all the files in the `data` directory, but we'll come back to that later.

### Takeaway

- `Channel.fromPath()` creates a channel with files matching a pattern
- Each file is emitted as a separate element in the channel
- We can use a glob pattern to match multiple files
- Files are automatically converted to Path objects with full attributes
- The `.view()` method allows inspection of channel contents

---

## 4. Extracting Patient Metadata from Filenames

One of the most common tasks in bioinformatics workflows is extracting metadata from filenames. This is usually feasible when working with sequencing data, where filenames often contain information about the sample, condition, replicate, and read number.

This isn't ideal - metadata should never be embedded in filenames, but it's a common reality. We want to extract that metadata in a standardised manner so we can use it later.

Let's explore how to extract metadata from our FASTQ filenames using Nextflow's powerful data transformation capabilities.

### 4.1. Basic Metadata Extraction

First, let's modify our workflow to extract metadata from the filenames.

First we will grab the simpleName of the file, which includes the metadata, and return with the file. Then, we will separate out the metadata by underscores using tokenize. Finally, we will use string handling to remove additional text like "rep" which aren't required right now.

=== "After"

    ```groovy title="file_operations.nf" linenums="4"
        ch_fastq.map { myFile ->
            [ myFile.simpleName, myFile ]
        }
        .view()
    ```

=== "Before"

    ```groovy title="file_operations.nf" linenums="4"
        ch_fastq.view {
            println "File object class: ${myFile.class}"
            println "File name: ${it.name}"
            println "Simple name: ${it.simpleName}"
            println "Extension: ${it.extension}"
            println "Parent directory: ${it.parent}"
        }
    ```

```bash title="Test filename metadata extraction"
nextflow run file_operations.nf
```

```console title="Sample Metadata Output"
 N E X T F L O W   ~  version 25.04.3

Launching `file_operations.nf` [suspicious_mahavira] DSL2 - revision: ae8edc4e48

[patientA_rep1_normal_R2_001, /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R2_001.fastq.gz]
[patientA_rep1_normal_R1_001, /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz]
```

Note how we have separated the patient `simpleName`, which includes the metadata, from the `file` object. This is useful if we want to use the patient metadata in a later process.

### 4.2. Extracting Metadata from Filenames

Our metadata is embedded in the filename, but it's not in a standard format. We need to split up the filename into it's components which are separated by underscores.

Groovy includes a method called `tokenize()` which is perfect for this task.

=== "After"

    ```groovy title="file_operations.nf" linenums="4" hl_lines="2"
        ch_fastq.map { myFile ->
            [ myFile.simpleName.tokenize('_'), myFile ]
        }
    ```

=== "Before"

    ```groovy title="file_operations.nf" linenums="4" hl_lines="2"
        ch_fastq.map { myFile ->
            [ myFile.simpleName, myFile ]
        }
    ```

Once we run this, we should see the patient metadata as a list of strings, and the Path object as the second element in the tuple.

```bash title="Test filename tokenization"
nextflow run file_operations.nf
```

```console title="Sample Tokenize Output"
 N E X T F L O W   ~  version 25.04.3

Launching `file_operations.nf` [gigantic_gauss] DSL2 - revision: a39baabb57

[[patientA, rep1, normal, R1, 001], /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz]
[[patientA, rep1, normal, R2, 001], /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R2_001.fastq.gz]
```

Success! We've broken down our patient information from a single string into a list of strings. We can now handle each part of the patient information separately.

### 4.3. Using a map to organise the data

Our meta data is just a flat list at the moment. It's easy to use but hard to read. What is the item at index 3? Can you tell without checking?

A [map](https://www.baeldung.com/groovy-maps) is Groovy's version of a key-value store. Every item has a key and a value and we can refer to each key to get the value. This will make our code much easier to read, i.e. we go from this:

```groovy
data = [patientA, 1, normal, R1]

println data[3]
```

to this:

```groovy
data = [id: patientA, replicate: 1, type: normal, readNum: 1]

println data.readNum
```

Let's convert our flat list into a map now.

=== "After"

    ```groovy title="file_operations.nf" linenums="4" hl_lines="2-11"
        ch_fastq.map { myFile ->
            def (patient, replicate, type, readNum) = myFile.simpleName.tokenize('_')
            [
              [
                id: patient,
                replicate: replicate.replace('rep', ''),
                type: type,
                readNum: readNum.replace('rep', ''),
              ],
              myFile
            ]
        }
    ```

=== "Before"

    ```groovy title="file_operations.nf" linenums="4" hl_lines="2"
        ch_fastq.map { myFile ->
            [ myFile.simpleName.tokenize('_'), myFile ]
        }
    ```

Notice that we're simplifying a couple of the meta data items as we go (e.g. `readNum.replace('rep', '')`).

Now re-run the workflow:

```bash title="Test metadata map structure"
nextflow run file_operations.nf
```

```console title="Map Output"
 N E X T F L O W   ~  version 25.04.3

Launching `file_operations.nf` [infallible_swartz] DSL2 - revision: 7f4e68c0cb

[[id:patientA, replicate:rep1, type:normal, readNum:R2], /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R2_001.fastq.gz]
[[id:patientA, replicate:rep1, type:normal, readNum:R1], /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz]
```

We have converted our flat list into a map, and now we can refer to each bit of sample data by name instead of by index. This makes our code easier to read and more maintainable.

### Takeaway

- We can handle filenames in Nextflow with the power of a full programming language
- We can treat the filenames as strings to extract relevant information
- Use of methods like `tokenize()` and `replace()` allows us to manipulate strings in the filename
- The `.map()` operation transforms channel elements while preserving structure
- Structured metadata (maps) makes code more readable and maintainable than positional lists

Next up, we will look at how to handle paired-end reads.

---

## 5. Simplifying with Channel.fromFilePairs

Nextflow provides a specialized channel factory method for working with paired files: `Channel.fromFilePairs()`. This method automatically groups files that share a common prefix. This is particularly useful for paired-end sequencing data, where you have two files (e.g., R1 and R2) for each sample.

### 5.1. Basic Usage of fromFilePairs

Complete your `file_operations.nf` file with the following (deleting the map operation):

=== "After"

    ```groovy title="file_operations.nf" linenums="3" hl_lines="1"
        ch_fastq = Channel.fromFilePairs('data/patientA_rep1_normal_R{1,2}_001.fastq.gz')
            .view()
    ```

=== "Before"

    ```groovy title="file_operations.nf" linenums="3" hl_lines="1"
        ch_fastq = Channel.fromPath('data/patientA_rep1_normal_R*_001.fastq.gz')
        ch_fastq.map { myFile ->
            def (sample, replicate, type, readNum) = myFile.simpleName.tokenize('_')
            [
                [
                    id: sample,
                    replicate: replicate.replace('rep', ''),
                    type: type,
                    readNum: readNum,
                ],
                myFile
            ]
        }
        .view()
    ```

Run the workflow:

```bash title="Test Channel.fromFilePairs"
nextflow run file_operations.nf
```

The output will show the paired files grouped together:

```console title="Channel.fromFilePairs Output"
 N E X T F L O W   ~  version 25.04.3

Launching `file_operations.nf` [chaotic_cuvier] DSL2 - revision: 472265a440

[patientA_rep1_normal_R, [/workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz, /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R2_001.fastq.gz]]
```

Note the difference in data structure. Rather than being a list of results, we have a single result in the format `id, [ fastq1, fastq2 ]`. Nextflow has done the hard work of extracting the patient name by examining the shared prefix and using it as a patient id.

### 5.2. Extract metadata from file pairs

We still need the metadata. Our `map` operation from before won't work because it doesn't match the data structure, but we can modify it to work. We already have access to the patient name in the `id` variable, so we can use that to extract the metadata without grabbing the `simpleName` from the Path object like before.

=== "After"

    ```groovy title="file_operations.nf" linenums="3" hl_lines="2-13"
        ch_fastq = Channel.fromFilePairs('data/patientA_rep1_normal_R{1,2}_001.fastq.gz')
        ch_fastq.map { id, fastqs ->
            def (sample, replicate, type, readNum) = id.tokenize('_')
            [
                [
                    id: sample,
                    replicate: replicate.replace('rep', ''),
                    type: type
                ],
                fastqs
            ]
        }
        .view()
    ```

=== "Before"

    ```groovy title="file_operations.nf" linenums="2" hl_lines="3-11"
        ch_fastq = Channel.fromFilePairs('data/patientA_rep1_normal_R{1,2}_001.fastq.gz')
            .view()
    ```

```bash title="Test file pairs metadata extraction"
nextflow run file_operations.nf
```

```console title="File Pairs Output parsed"

 N E X T F L O W   ~  version 25.04.3

Launching `file_operations.nf` [prickly_stonebraker] DSL2 - revision: f62ab10a3f

[[id:patientA, replicate:1, type:normal], [/workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz, /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R2_001.fastq.gz]]
```

Notice, that this time we don't have a `readNum`. Because read 1 and read 2 are kept together, we do not need to track this in the meta data.

Well done! We have grabbed the metadata from the filenames and used them as values in the tuple.

### Takeaway

- [`Channel.fromFilePairs()` automatically finds and pairs related files](https://www.nextflow.io/docs/latest/reference/channel.html#fromfilepairs)
- This simplifies handling paired-end reads in your pipeline
- Paired files can be grouped as `[id, [file1, file2]]` tuples
- Metadata extraction can be done from the paired file ID rather than individual files

---

## 6. Using File Operations in Processes

Now let's put all this together in a simple process to reinforce how to use file operations in a Nextflow process.

### 6.1. Create the process

We'll keep it simple and make a process called `ANALYZE_READS` that takes in a tuple of metadata and a pair of fastq files and analyses them. We could imagine this is an alignment, or variant calling or any other step.

Add the following to the top of your `file_operations.nf` file:

=== "After"

    ```groovy title="file_operations.nf - process example" linenums="1" hl_lines="1-24"
    process ANALYZE_READS {
        tag "${meta.id}"

        publishDir "results/${meta.id}", mode: 'copy'

        input:
        tuple val(meta), path(fastqs)

        output:
        tuple val(meta.id), path("${meta.id}_stats.txt")

        script:
        """
        echo "Sample metadata: ${meta.id}" > ${meta.id}_stats.txt
        echo "Replicate: ${meta.replicate}" >> ${meta.id}_stats.txt
        echo "Type: ${meta.type}" >> ${meta.id}_stats.txt
        echo "Read 1: ${fastqs[0]}" >> ${meta.id}_stats.txt
        echo "Read 2: ${fastqs[1]}" >> ${meta.id}_stats.txt
        echo "Read 1 size: \$(gunzip -dc ${fastqs[0]} | wc -l | awk '{print \$1/4}') reads" >> ${meta.id}_stats.txt
        echo "Read 2 size: \$(gunzip -dc ${fastqs[1]} | wc -l | awk '{print \$1/4}') reads" >> ${meta.id}_stats.txt
        """
    }

    workflow {
    ```

=== "Before"

    ```groovy title="file_operations.nf" linenums="1" hl_lines="1"
    workflow {
    ```

!!! note

    We are calling our map '`meta`'. For a more in-depth introduction to meta maps, see [Working with metadata](./metadata.md).

### 6.2. Implement the process in the workflow

Then implement the process in the workflow:

=== "After"

    ```groovy title="file_operations.nf" linenums="26" hl_lines="2 13"
        ch_fastq = Channel.fromFilePairs('data/patientA_rep1_normal_R{1,2}_001.fastq.gz')
        ch_samples = ch_fastq.map { id, fastqs ->
            def (sample, replicate, type, readNum) = id.tokenize('_')
            [
                [
                    id: sample,
                    replicate: replicate.replace('rep', ''),
                    type: type
                ],
                fastqs
            ]
        }
        ANALYZE_READS(ch_samples)
    }
    ```

=== "Before"

    ```groovy title="file_operations.nf" linenums="26" hl_lines="2 13"
        ch_fastq = Channel.fromFilePairs('data/patientA_rep1_normal_R{1,2}_001.fastq.gz')
        ch_fastq.map { id, fastqs ->
            def (sample, replicate, type, readNum) = id.tokenize('_')
            [
                [
                    id: sample,
                    replicate: replicate.replace('rep', ''),
                    type: type
                ],
                fastqs
            ]
        }
        .view()
    }
    ```

```bash title="Test ANALYZE_READS process"
nextflow run file_operations.nf
```

```console title="ANALYZE_READS Output"
 N E X T F L O W   ~  version 25.04.3

Launching `./file_operations.nf` [shrivelled_cori] DSL2 - revision: b546a31769

executor >  local (1)
[b5/110360] process > ANALYZE_READS (patientA) [100%] 1 of 1 ✔
```

We should see the following files in the `results/patientA` directory:

```console title="Results Directory"
> tree results/patientA
results/patientA
└── patientA_stats.txt
```

The process took our inputs and created a new file with the patient metadata. Based on what you learned in hello-nextflow, what occurred in the working directory?

### 6.3. Include many more patients

Remember Channel.fromPath() accepts a _glob_ as input, which means it can accept any number of files that match the pattern. Therefore if we want to include all the patients we can just modify the input string to include more patients.

=== "After"

    ```groovy title="file_operations.nf" linenums="26"
        ch_fastq = Channel.fromFilePairs('data/*_R{1,2}_001.fastq.gz')
    ```

=== "Before"

    ```groovy title="file_operations.nf" linenums="26"
        ch_fastq = Channel.fromFilePairs('data/patientA_rep1_normal_R{1,2}_001.fastq.gz')
    ```

Run the pipeline now and see all the results:

```bash title="Test processing multiple samples"
nextflow run file_operations.nf
```

```console title="ANALYZE_READS Multiple Samples"
 N E X T F L O W   ~  version 25.04.3

Launching `./file_operations.nf` [big_stonebraker] DSL2 - revision: f7f9b8a76c

executor >  local (8)
[d5/441891] process > ANALYZE_READS (patientC) [100%] 8 of 8 ✔
```

Check the results directory now:

```console title="Results Directory"
> tree results
results
├── patientA
│   └── patientA_stats.txt
├── patientB
│   └── patientB_stats.txt
└── patientC
    └── patientC_stats.txt
```

See how we have analyzed all the patients in one go!

Wait, we have a problem. We have 2 replicates for patientA, but only 1 output file! We are overwriting the output file each time.

### 6.4. Make the published files unique

Since we have access to the patient metadata, we can use it to make the output files unique.

=== "After"

    ```groovy title="file_operations.nf" linenums="4"
        publishDir "results/${meta.type}/${meta.id}/${meta.replicate}", mode: 'copy'
    ```

=== "Before"

    ```groovy title="file_operations.nf" linenums="4"
        publishDir "results/${id}", mode: 'copy'
    ```

We have grabbed the metadata from the patients and used it to construct an output directory for each patient.

Run the pipeline now and see all the results. Remove the results directory first to give yourself a clean workspace:

```bash title="Test unique published files"
rm -r results
nextflow run file_operations.nf
```

```console title="Results Directory"
 N E X T F L O W   ~  version 25.04.3

Launching `./file_operations.nf` [insane_swartz] DSL2 - revision: fff18abe6d

executor >  local (8)
[e3/449081] process > ANALYZE_READS (patientC) [100%] 8 of 8 ✔
```

Check the results directory now:

```console title="Results Directory"
> tree results
results/
├── normal
│   ├── patientA
│   │   ├── 1
│   │   │   └── patientA_stats.txt
│   │   └── 2
│   │       └── patientA_stats.txt
│   ├── patientB
│   │   └── 1
│   │       └── patientB_stats.txt
│   └── patientC
│       └── 1
│           └── patientC_stats.txt
└── tumor
    ├── patientA
    │   ├── 1
    │   │   └── patientA_stats.txt
    │   └── 2
    │       └── patientA_stats.txt
    ├── patientB
    │   └── 1
    │       └── patientB_stats.txt
    └── patientC
        └── 1
            └── patientC_stats.txt
```

See how using patient metadata as values gives us powerful flexibility in our pipeline. By propagating metadata alongside our data in tuples, we can:

1. Create organized output directories based on patient attributes
2. Make decisions in processes based on patient properties
3. Split, join, and recombine data based on metadata values

This pattern of keeping metadata explicit and attached to the data (rather than encoded in filenames) is a core best practice in Nextflow that enables building robust, maintainable bioinformatics workflows. Learn more about this in [Working with metadata](./metadata.md)

### Takeaway

- The `publishDir` directive can organize outputs based on metadata values
- Metadata in tuples enables structured organization of results
- This approach creates maintainable workflows with clear data provenance
- Processes can take tuples of metadata and files as input
- The `tag` directive provides process identification in execution logs
- Workflow structure separates channel creation from process execution

## Summary

In this side quest, you've learned how to work with files in Nextflow, from basic operations to more advanced techniques for handling collections of files. Here's a summary of what we covered:

1. **Basic File Operations**: We created Path objects with `file()` and accessed file attributes like name, extension, and parent directory, learning the difference between strings and Path objects.

2. **Using Remote Files**: We learned how to transparently switch between local and remote files using URIs, demonstrating Nextflow's ability to handle files from various sources without changing workflow logic.

3. **Reading files using the `fromPath()` channel factory**: We created channels from file patterns with `Channel.fromPath()` and viewed their file attributes, including object types.

4. **Extracting Patient Metadata from Filenames**: We used `tokenize()` and `replace()` to extract and structure metadata from filenames, converting them to organized maps.

5. **Simplifying with Channel.fromFilePairs**: We used `Channel.fromFilePairs()` to automatically pair related files and extract metadata from paired file IDs.

6. **Using File Operations in Processes**: We integrated file operations into Nextflow processes with proper input handling, using `publishDir` to organize outputs based on metadata.

These techniques will help you build more efficient and maintainable workflows, especially when working with large numbers of files with complex naming conventions.

### Key Concepts

- **Path Object Creation**

  ```groovy
  // Create a Path object from a string path
  myFile = file('path/to/file.txt')
  ```

- **File Attributes**

  ```groovy
  // Get file attributes
  println myFile.name       // file.txt
  println myFile.baseName   // file
  println myFile.extension  // txt
  println myFile.parent     // path/to
  ```

- **Channel Creation from Files**

  ```groovy
  // Create a channel from a file pattern
  ch_fastq = Channel.fromPath('data/*.fastq.gz')

  // Create a channel from paired files
  ch_pairs = Channel.fromFilePairs('data/*_R{1,2}_001.fastq.gz')
  ```

- **Extracting Metadata**

  ```groovy
  // Extract metadata with tokenize
  def name = file.name.tokenize('_')
  def patientId = name[0]
  def replicate = name[1].replace('rep', '')
  def type = name[2]
  def readNum = name[3].replace('R', '')
  ```

- **Using remote files**

  ```groovy
  // Use a local file
  myFile = file('path/to/file.txt')

  // Use a file on FTP
  myFile = file('ftp://path/to/file.txt')

  // Use a file on HTTPS
  myFile = file('https://path/to/file.txt')

  // Use a file on S3
  myFile = file('s3://path/to/file.txt')

  // Use a file on Azure Blob Storage
  myFile = file('az://path/to/file.txt')

  // Use a file on Google Cloud Storage
  myFile = file('gs://path/to/file.txt')
  ```

## Resources

- [Nextflow Documentation: Working with Files](https://www.nextflow.io/docs/latest/working-with-files.html)
- [Channel.fromPath](https://www.nextflow.io/docs/latest/channel.html#frompath)
- [Channel.fromFilePairs](https://www.nextflow.io/docs/latest/channel.html#fromfilepairs)
