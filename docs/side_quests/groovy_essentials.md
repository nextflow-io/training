# Groovy Essentials for Nextflow Developers

Nextflow is built on Apache Groovy, a powerful dynamic language that runs on the Java Virtual Machine. This foundation gives Nextflow its flexibility and expressiveness, but it also creates a common source of confusion for developers.

**Here's the challenge:** Most Nextflow tutorials focus on the workflow orchestration - channels, processes, and data flow - but when you need to manipulate data, parse filenames, implement conditional logic, or handle errors gracefully, you're actually writing Groovy code. Many developers don't realize when they've crossed this boundary.

**Why does this matter?** The difference between a brittle workflow that breaks on unexpected input and a robust pipeline that adapts gracefully often comes down to understanding and leveraging Groovy's powerful features within your Nextflow workflows.

**The common struggle:** Most Nextflow developers can write basic workflows, but they hit walls when they need to:
- Process messy, real-world data with missing fields or inconsistent formats
- Extract metadata from complex file naming schemes
- Route samples through different analysis strategies based on their characteristics
- Handle errors gracefully instead of crashing on invalid input
- Build reusable, maintainable code that doesn't repeat the same patterns everywhere

Understanding where Nextflow ends and Groovy begins is crucial for effective workflow development. Nextflow provides channels, processes, and workflow orchestration, while Groovy handles data manipulation, string processing, conditional logic, and general programming tasks within your workflow scripts.

This side quest will bridge that gap by taking you on a hands-on journey from basic concepts to production-ready mastery. We'll transform a simple CSV-reading workflow into a sophisticated, production-ready bioinformatics pipeline that can handle any dataset thrown at it. Starting with a basic workflow that processes sample metadata, we'll evolve it step-by-step through realistic challenges you'll face in production:

- **Messy data?** We'll add robust parsing and null-safe operators, learning to distinguish between Nextflow and Groovy constructs
- **Complex file naming schemes?** We'll master regex patterns and string manipulation for bioinformatics file names
- **Need intelligent sample routing?** We'll implement conditional logic and strategy selection, transforming file collections into command-line arguments
- **Worried about failures?** We'll add comprehensive error handling and validation patterns
- **Code getting repetitive?** We'll learn functional programming with closures and composition, mastering essential Groovy operators like safe navigation and Elvis
- **Processing thousands of samples?** We'll leverage powerful collection operations for file path manipulations

---

## 0. Warmup

### 0.1. Prerequisites

Before taking on this side quest you should:

- Complete the [Hello Nextflow](../hello_nextflow/README.md) tutorial or have equivalent experience
- Understand basic Nextflow concepts (processes, channels, workflows)
- Have basic familiarity with Groovy syntax (variables, maps, lists)

This tutorial will explain Groovy concepts as we encounter them, so you don't need extensive prior Groovy knowledge. We'll start with fundamental concepts and build up to advanced patterns.

### 0.2. Starting Point

Let's move into the project directory and explore our working materials.

```bash title="Navigate to project directory"
cd side-quests/groovy_essentials
```

You'll find a `data` directory with sample files and a main workflow file that we'll evolve throughout this tutorial.

```console title="Directory contents"
> tree
.
├── data
│   ├── metadata
│   │   └── analysis_parameters.yaml
│   ├── samples.csv
│   └── sequences
│       ├── sample_001.fastq
│       ├── sample_002.fastq
│       └── sample_003.fastq
├── main.nf
├── nextflow.config
├── README.md
└── templates
    └── analysis_script.sh

5 directories, 9 files
```

Our sample CSV contains information about biological samples that need different processing based on their characteristics:

```console title="samples.csv"
sample_id,organism,tissue_type,sequencing_depth,file_path,quality_score
SAMPLE_001,human,liver,30000000,data/sequences/sample_001.fastq,38.5
SAMPLE_002,mouse,brain,25000000,data/sequences/sample_002.fastq,35.2
SAMPLE_003,human,kidney,45000000,data/sequences/sample_003.fastq,42.1
```

We'll use this realistic dataset to explore practical Groovy techniques that you'll encounter in real bioinformatics workflows.

---

## 1. Nextflow vs Groovy: Understanding the Boundaries

### 1.1. Identifying What's What

One of the most common sources of confusion for Nextflow developers is understanding when they're working with Nextflow constructs versus Groovy language features. Let's build a workflow step by step to see how they work together.

#### Step 1: Basic Nextflow Workflow

Start with a simple workflow that just reads the CSV file (we've already done this for you in `main.nf`):

```groovy title="main.nf" linenums="1"
workflow {
    ch_samples = Channel.fromPath("./data/samples.csv")
        .splitCsv(header: true)
        .view()
}
```

The `workflow` block defines our pipeline structure, while `Channel.fromPath()` creates a channel from a file path. The `.splitCsv()` operator processes the CSV file and converts each row into a map data structure.

Run this workflow to see the raw CSV data:

```bash title="Test basic workflow"
nextflow run main.nf
```

You should see output like:
```console title="Raw CSV data"
[id:sample_001, organism:human, tissue_type:liver, sequencing_depth:30000000, file_path:data/sequences/sample_001.fastq, quality_score:38.5]
[id:sample_002, organism:mouse, tissue_type:brain, sequencing_depth:25000000, file_path:data/sequences/sample_002.fastq, quality_score:35.2]
[id:sample_003, organism:human, tissue_type:kidney, sequencing_depth:45000000, file_path:data/sequences/sample_003.fastq, quality_score:42.1]
```

#### Step 2: Adding the Map Operator

Now we're going to use some Groovy code to transform the data, using the `.map()` operator you will probably already be familiar with. This operator takes a 'closure' where we can write Groovy code to transform each item.

!!! note

    A **closure** is a block of code that can be passed around and executed later. Think of it as a function that you define inline. In Groovy, closures are written with curly braces `{ }` and can take parameters. They're fundamental to how Nextflow operators work and if you've been writing Nextflow for a while, you may already have been using them without realizing it!

Here's what that map operation looks like:

=== "After"

    ```groovy title="main.nf" linenums="2" hl_lines="3-6"
        ch_samples = Channel.fromPath("./data/samplesheet.csv")
            .splitCsv(header: true)
            .map { row ->
                return row
            }
            .view()
    ```

=== "Before"

    ```groovy title="main.nf" linenums="2" hl_lines="3"
        ch_samples = Channel.fromPath("./data/samplesheet.csv")
            .splitCsv(header: true)
            .view()
    ```

The closure we've added is `{ row -> return row }`. We've named the parameter `row`, but it could be called anything, so you could write `.map { item -> ... }` or `.map { sample -> ... }` and it would work exactly the same way. It's also possible not to name the parameter at all and just use the implicit variable `it`, like `.map { return it }`, but naming it makes the code clearer.

When Nextflow processes each item in the channel, it passes that item to your closure as the parameter you named. So if your channel contains CSV rows, `row` will hold one complete row at a time.

Apply this change and run the workflow:

```bash title="Test map operator"
nextflow run main.nf
```

You'll see the same output as before, because we're simply returning the input unchanged. This confirms that the map operator is working correctly. Now let's start transforming the data.

#### Step 3: Creating a Map Data Structure

Now we're going to write **pure Groovy code** inside our closure. Everything from this point forward is Groovy syntax and methods, not Nextflow operators.

=== "After"

    ```groovy title="main.nf" linenums="2" hl_lines="4-12"
        ch_samples = Channel.fromPath("./data/samplesheet.csv")
            .splitCsv(header: true)
            .map { row ->
                // This is all Groovy code now!
                def sample_meta = [
                    id: row.sample_id.toLowerCase(),
                    organism: row.organism,
                    tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                    depth: row.sequencing_depth.toInteger(),
                    quality: row.quality_score.toDouble()
                ]
                return sample_meta
            }
            .view()
    ```

=== "Before"

    ```groovy title="main.nf" linenums="2" hl_lines="4"
        ch_samples = Channel.fromPath("./data/samplesheet.csv")
            .splitCsv(header: true)
            .map { row ->
                return row
            }
            .view()
    ```

Notice how we've left Nextflow syntax behind and are now writing pure Groovy code. A map is a key-value data structure similar to dictionaries in Python, objects in JavaScript, or hashes in Ruby. It lets us store related pieces of information together. In this map, we're storing the sample ID, organism, tissue type, sequencing depth, and quality score.

We use Groovy's string manipulation methods like `.toLowerCase()` and `.replaceAll()` to clean up our data, and type conversion methods like `.toInteger()` and `.toDouble()` to convert string data from the CSV into the appropriate numeric types.

Apply this change and run the workflow:

```bash title="Test map data structure"
nextflow run main.nf
```

You should see the refined map output like:

```console title="Transformed metadata"
[id:sample_001, organism:human, tissue:liver, depth:30000000, quality:38.5]
[id:sample_002, organism:mouse, tissue:brain, depth:25000000, quality:35.2]
[id:sample_003, organism:human, tissue:kidney, depth:45000000, quality:42.1]
```

#### Step 4: Adding Conditional Logic

Now let's add more Groovy logic - this time using a ternary operator to make decisions based on data values.

Make the following change:

=== "After"

    ```groovy title="main.nf" linenums="2" hl_lines="11-12"
        ch_samples = Channel.fromPath("./data/samplesheet.csv")
            .splitCsv(header: true)
            .map { row ->
                def sample_meta = [
                    id: row.sample_id.toLowerCase(),
                    organism: row.organism,
                    tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                    depth: row.sequencing_depth.toInteger(),
                    quality: row.quality_score.toDouble()
                ]
                def priority = sample_meta.quality > 40 ? 'high' : 'normal'
                return sample_meta + [priority: priority]
            }
            .view()
    ```

=== "Before"

    ```groovy title="main.nf" linenums="2" hl_lines="11"
        ch_samples = Channel.fromPath("./data/samplesheet.csv")
            .splitCsv(header: true)
            .map { row ->
                def sample_meta = [
                    id: row.sample_id.toLowerCase(),
                    organism: row.organism,
                    tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                    depth: row.sequencing_depth.toInteger(),
                    quality: row.quality_score.toDouble()
                ]
                return sample_meta
            }
            .view()
    ```

The ternary operator is a shorthand for an if/else statement that follows the pattern `condition ? value_if_true : value_if_false`. This line means: "If the quality is greater than 40, use 'high', otherwise use 'normal'". Its cousin, the **Elvis operator** (`?:`), provides default values when something is null or empty - we'll explore that pattern later in this tutorial.

The map addition operator `+` creates a **new map** rather than modifying the existing one. This line creates a new map that contains all the key-value pairs from `sample_meta` plus the new `priority` key.

!!! Note

    Using the addition operator `+` creates a new map rather than modifying the existing one, which is a useful practice to adopt. Never directly modify maps passed into closures, as this can lead to unexpected behavior in Nextflow. This is especially important because in Nextflow workflows, the same data often flows through multiple channel operations or gets processed by different processes simultaneously. When multiple operations reference the same map object, modifying it in-place can cause unpredictable side effects - one operation might change data that another operation is still using. By creating new maps instead of modifying existing ones, you ensure that each operation works with its own clean copy of the data, making your workflows more predictable and easier to debug.

Run the modified workflow:

```bash title="Test conditional logic"
nextflow run main.nf
```

You should see output like:

```console title="Metadata with priority"
[id:sample_001, organism:human, tissue:liver, depth:30000000, quality:38.5, priority:normal]
[id:sample_002, organism:mouse, tissue:brain, depth:25000000, quality:35.2, priority:normal]
[id:sample_003, organism:human, tissue:kidney, depth:45000000, quality:42.1, priority:high]
```

We've successfully added conditional logic to enrich our metadata with a priority level based on quality scores.

#### Step 5: Combining Maps and Returning Results

So far, we've only been returning what Nextflow community calls the 'meta map', and we've been ignoring the files those metadata relate to. But if you're writing Nextflow workflows, you probably want to do something with those files.

Let's output a channel structure comprising a tuple of 2 elements: the enriched metadata map and the corresponding file path. This is a common pattern in Nextflow for passing data to processes.

=== "After"

    ```groovy title="main.nf" linenums="2" hl_lines="12"
        ch_samples = Channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row ->
                def sample_meta = [
                    id: row.sample_id.toLowerCase(),
                    organism: row.organism,
                    tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                    depth: row.sequencing_depth.toInteger(),
                    quality: row.quality_score.toDouble()
                ]
                def priority = sample_meta.quality > 40 ? 'high' : 'normal'
                return [sample_meta + [priority: priority], file(row.file_path)]
            }
            .view()
    ```

=== "Before"

    ```groovy title="main.nf" linenums="2" hl_lines="12"
        ch_samples = Channel.fromPath("./data/samplesheet.csv")
            .splitCsv(header: true)
            .map { row ->
                def sample_meta = [
                    id: row.sample_id.toLowerCase(),
                    organism: row.organism,
                    tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                    depth: row.sequencing_depth.toInteger(),
                    quality: row.quality_score.toDouble()
                ]
                def priority = sample_meta.quality > 40 ? 'high' : 'normal'
                return sample_meta + [priority: priority]
            }
            .view()
    ```

Apply this change and run the workflow:

```bash title="Test complete workflow"
nextflow run main.nf
```

You should see output like:
```console title="Complete workflow output"
[[id:sample_001, organism:human, tissue:liver, depth:30000000, quality:38.5, priority:normal], /workspaces/training/side-quests/groovy_essentials/data/sequences/sample_001.fastq]
[[id:sample_002, organism:mouse, tissue:brain, depth:25000000, quality:35.2, priority:normal], /workspaces/training/side-quests/groovy_essentials/data/sequences/sample_002.fastq]
[[id:sample_003, organism:human, tissue:kidney, depth:45000000, quality:42.1, priority:high], /workspaces/training/side-quests/groovy_essentials/data/sequences/sample_003.fastq]
```

This `[meta, file]` tuple structure is a common pattern in Nextflow for passing both metadata and associated files to processes.

!!! note

    **Maps and Metadata**: Maps are fundamental to working with metadata in Nextflow. For a more detailed explanation of working with metadata maps, see the [Working with metadata](./metadata.md) side quest.

Our workflow demonstrates the core pattern: **Nextflow constructs** (`workflow`, `Channel.fromPath()`, `.splitCsv()`, `.map()`, `.view()`) orchestrate data flow, while **basic Groovy constructs** (maps `[key: value]`, string methods, type conversions, ternary operators) handle the data processing logic inside the `.map()` closure.

### 1.2. Distinguishing Nextflow operators from Groovy functions

So far, so good, we can distinguish between Nextflow constructs and basic Groovy constructs. But what about when the syntax overlaps?

A perfect example of this confusion is the `collect` operation, which exists in both contexts but does completely different things. Groovy's `collect` transforms each element, while Nextflow's `collect` gathers all channel elements into a single-item channel.

Let's demonstrate this with some sample data, starting by refreshing ourselves on what the Nextflow `collect()` operator does. Check out `collect.nf`:

```groovy title="collect.nf" linenums="1"
def sample_ids = ['sample_001', 'sample_002', 'sample_003']

// Nextflow collect() - groups multiple channel emissions into one
ch_input = Channel.fromList(sample_ids)
ch_input.view { "Individual channel item: ${it}" }
ch_collected = ch_input.collect()
ch_collected.view { "Nextflow collect() result: ${it} (${it.size()} items grouped into 1)" }
```

We're using the `fromList()` channel factory to create a channel that emits each sample ID as a separate item, and we use `view()` to print each item as it flows through the channel.Then we apply Nextflow's `collect()` operator to gather all items into a single list and use a second `view()` to print the collected result which appears as a single item containing a list of all sample IDs. We've changed the structure of the channel, but we haven't changed the data itself.

Run the workflow to confirm this:

```bash title="Test collect operations"
nextflow run collect.nf
```

```console title="Different collect behaviors"
 N E X T F L O W   ~  version 25.04.6

Launching `collect.nf` [loving_mendel] DSL2 - revision: e8d054a46e

Individual channel item: sample_001
Individual channel item: sample_002
Individual channel item: sample_003
Nextflow collect() result: [sample_001, sample_002, sample_003] (3 items grouped into 1)
```

Now let's see Groovy's `collect` method in action. Modify `collect.nf` to apply Groovy's `collect` method to the original list of sample IDs:

=== "After"

    ```groovy title="main.nf" linenums="1" hl_lines="9-13"
        def sample_ids = ['sample_001', 'sample_002', 'sample_003']

        // Nextflow collect() - groups multiple channel emissions into one
        ch_input = Channel.fromList(sample_ids)
        ch_input.view { "Individual channel item: ${it}" }
        ch_collected = ch_input.collect()
        ch_collected.view { "Nextflow collect() result: ${it} (${it.size()} items grouped into 1)" }

        // Groovy collect - transforms each element, preserves structure
        def formatted_ids = sample_ids.collect { id ->
            id.toUpperCase().replace('SAMPLE_', 'SPECIMEN_')
        }
        println "Groovy collect result: ${formatted_ids} (${sample_ids.size()} items transformed into ${formatted_ids.size()})"
    ```

=== "Before"

    ```groovy title="main.nf" linenums="1"
        def sample_ids = ['sample_001', 'sample_002', 'sample_003']

        // Nextflow collect() - groups multiple channel emissions into one
        ch_input = Channel.fromList(sample_ids)
        ch_input.view { "Individual channel item: ${it}" }
        ch_collected = ch_input.collect()
        ch_collected.view { "Nextflow collect() result: ${it} (${it.size()} items grouped into 1)" }
    ```

Run the modified workflow:

```bash title="Test Groovy collect"
nextflow run collect.nf
```

```console title="Groovy collect results" hl_lines="9"
 N E X T F L O W   ~  version 25.04.6

Launching `collect.nf` [cheeky_stonebraker] DSL2 - revision: 2d5039fb47

Groovy collect result: [SPECIMEN_001, SPECIMEN_002, SPECIMEN_003] (3 items transformed into 3)
Individual channel item: sample_001
Individual channel item: sample_002
Individual channel item: sample_003
Nextflow collect() result: [sample_001, sample_002, sample_003] (3 items grouped into 1)
```

This time, have NOT changed the structure of the data, we still have 3 items in the list, but we HAVE transformed each item using Groovy's `collect` method to produce a new list with modified values. This is sort of like using the `map` operator in Nextflow, but it's pure Groovy code operating on a standard Groovy list.

`collect` is an extreme case we're using here to make a point. The key lesson is that when you're writing workflows always distinguish between **Groovy constructs** (data structures) and **Nextflow constructs** (channels/workflows). Operations can share names but behave completely differently.

### Takeaway

In this section, you've learned:

- **It takes both Nextflow and Groovy**: Nextflow provides the workflow structure and data flow, while Groovy provides the data manipulation and logic
- **Distinguishing Nextflow from Groovy**: How to identify which language construct you're using given the context
- **Context matters**: The same operation name can have completely different behaviors

Understanding these boundaries is essential for debugging, documentation, and writing maintainable workflows.

Next we'll dive deeper into Groovy's powerful string processing capabilities, which are essential for handling real-world data.

---

## 2. Advanced String Processing for Bioinformatics

The difference between a brittle workflow that breaks on unexpected input and a robust pipeline that adapts gracefully often comes down to mastering Groovy's string processing capabilities. Let's transform our pipeline to handle the messy realities of real-world bioinformatics data.

### 2.1. Pattern Matching and Regular Expressions

Many bioinformatics workflows encounter files with complex naming conventions that encode important metadata. Let's see how Groovy's pattern matching can extract this information automatically.

We're going to return to our `main.nf` workflow and add some pattern matching logic to extract additional sample information from file names. The FASTQ files in our dataset follow Illumina-style naming conventions with names like `SAMPLE_001_S1_L001_R1_001.fastq.gz`. These might look cryptic, but they actually encode useful metadata like sample ID, lane number, and read direction. We're going to use Groovy's regex capabilities to parse these names.

Make the following change to your existing `main.nf` workflow:

=== "After"

    ```groovy title="main.nf" linenums="2" hl_lines="9-19,21"
            .map { row ->
                def sample_meta = [
                    id: row.sample_id.toLowerCase(),
                    organism: row.organism,
                    tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                    depth: row.sequencing_depth.toInteger(),
                    quality: row.quality_score.toDouble()
                ]
                def fastq_path = file(row.file_path)

                def m = (fastq_path.name =~ /^(.+)_S(\d+)_L(\d{3})_(R[12])_(\d{3})\.fastq(?:\.gz)?$/)
                def file_meta = m ? [
                    sample_num: m[0][2].toInteger(),
                    lane: m[0][3],
                    read: m[0][4],
                    chunk: m[0][5]
                ] : [:]

                def priority = sample_meta.quality > 40 ? 'high' : 'normal'
                return [sample_meta + file_meta + [priority: priority], fastq_path]
            }
    ```

=== "Before"

    ```groovy title="main.nf" linenums="2"
        ch_samples = Channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row ->
                def sample_meta = [
                    id: row.sample_id.toLowerCase(),
                    organism: row.organism,
                    tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                    depth: row.sequencing_depth.toInteger(),
                    quality: row.quality_score.toDouble()
                ]
                def priority = sample_meta.quality > 40 ? 'high' : 'normal'
                return [sample_meta + [priority: priority], file(row.file_path)]
            }
            .view()
    ```

This demonstrates key Groovy string processing concepts:

1. **Regular expression literals** using `~/pattern/` syntax - this creates a regex pattern without needing to escape backslashes
2. **Pattern matching** with the `=~` operator - this attempts to match a string against a regex pattern
3. **Matcher objects** that capture groups with `[0][1]`, `[0][2]`, etc. - `[0]` refers to the entire match, `[1]`, `[2]`, etc. refer to captured groups in parentheses

The regex pattern `^(.+)_S(\d+)_L(\d{3})_(R[12])_(\d{3})\.fastq(?:\.gz)?$` is designed to match the Illumina-style naming convention and capture key components, namely the sample number, lane number, read direction, and chunk number.

Run the modified workflow:

```bash title="Test pattern matching"
nextflow run main.nf
```

You should see output with metadata enriched from the file names, like

```console title="Metadata with file parsing"
 N E X T F L O W   ~  version 25.04.6

Launching `main.nf` [clever_pauling] DSL2 - revision: 605d2058b4

[[id:sample_001, organism:human, tissue:liver, depth:30000000, quality:38.5, sample_num:1, lane:001, read:R1, chunk:001, priority:normal], /Users/jonathan.manning/projects/training/side-quests/groovy_essentials/data/sequences/SAMPLE_001_S1_L001_R1_001.fastq.gz]
[[id:sample_002, organism:mouse, tissue:brain, depth:25000000, quality:35.2, sample_num:2, lane:001, read:R1, chunk:001, priority:normal], /Users/jonathan.manning/projects/training/side-quests/groovy_essentials/data/sequences/SAMPLE_002_S2_L001_R1_001.fastq.gz]
[[id:sample_003, organism:human, tissue:kidney, depth:45000000, quality:42.1, sample_num:3, lane:001, read:R1, chunk:001, priority:high], /Users/jonathan.manning/projects/training/side-quests/groovy_essentials/data/sequences/SAMPLE_003_S3_L001_R1_001.fastq.gz]
```

### 2.2. Creating Reusable closures

You may have noticed that the content of our map operation is getting quite long and complex. To keep our workflow maintainable, it's a good idea to break out complex logic into reusable functions or closures.

To do that we simply define a closure using the assignment operator `=` and the `{}` syntax, within the `workflow{}`. Then we can call that closure by name inside our map operation using standard function call syntax (not the curly braces).

Make that change like so:

=== "After"

    ```groovy title="main.nf" linenums="1" hl_lines="3-23,27"
        workflow {

            separateMetadata = { row ->
                def sample_meta = [
                    id: row.sample_id.toLowerCase(),
                    organism: row.organism,
                    tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                    depth: row.sequencing_depth.toInteger(),
                    quality: row.quality_score.toDouble()
                ]
                def fastq_path = file(row.file_path)

                def m = (fastq_path.name =~ /^(.+)_S(\d+)_L(\d{3})_(R[12])_(\d{3})\.fastq(?:\.gz)?$/)
                def file_meta = m ? [
                    sample_num: m[0][2].toInteger(),
                    lane: m[0][3],
                    read: m[0][4],
                    chunk: m[0][5]
                ] : [:]

                def priority = sample_meta.quality > 40 ? 'high' : 'normal'
                return [sample_meta + file_meta + [priority: priority], fastq_path]
            }

            ch_samples = Channel.fromPath("./data/samples.csv")
                .splitCsv(header: true)
                    .map(separateMetadata)
                    .view()
        }
    ```

=== "Before"

    ```groovy title="main.nf" linenums="1" hl_lines="4-24"
        workflow {
            ch_samples = Channel.fromPath("./data/samples.csv")
                .splitCsv(header: true)
                    .map { row ->
                        def sample_meta = [
                            id: row.sample_id.toLowerCase(),
                            organism: row.organism,
                            tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                            depth: row.sequencing_depth.toInteger(),
                            quality: row.quality_score.toDouble()
                        ]
                        def fastq_path = file(row.file_path)

                        def m = (fastq_path.name =~ /^(.+)_S(\d+)_L(\d{3})_(R[12])_(\d{3})\.fastq(?:\.gz)?$/)
                        def file_meta = m ? [
                            sample_num: m[0][2].toInteger(),
                            lane: m[0][3],
                            read: m[0][4],
                            chunk: m[0][5]
                        ] : [:]

                        def priority = sample_meta.quality > 40 ? 'high' : 'normal'
                        return [sample_meta + file_meta + [priority: priority], fastq_path]
                    }
                    .view()
        }
    ```

By doing this we've reduced the actual workflow logic down to something really trivial:

```groovy title="minimal workflow"
            ch_samples = Channel.fromPath("./data/samples.csv")
                .splitCsv(header: true)
                    .map(separateMetadata)
                    .view()
```

... which makes the logic much easier to read and understand at a glance. The closure `separateMetadata` encapsulates all the complex logic for parsing and enriching metadata, making it reusable and testable.

You can run that to make sure it still works:

```bash title="Test reusable closure"
nextflow run main.nf
```

```console title="Closure results"
 N E X T F L O W   ~  version 25.04.6

Launching `main.nf` [tender_archimedes] DSL2 - revision: 8bfb9b2485

[[id:sample_001, organism:human, tissue:liver, depth:30000000, quality:38.5, sample_num:1, lane:001, read:R1, chunk:001, priority:normal], /Users/jonathan.manning/projects/training/side-quests/groovy_essentials/data/sequences/SAMPLE_001_S1_L001_R1_001.fastq.gz]
[[id:sample_002, organism:mouse, tissue:brain, depth:25000000, quality:35.2, sample_num:2, lane:001, read:R1, chunk:001, priority:normal], /Users/jonathan.manning/projects/training/side-quests/groovy_essentials/data/sequences/SAMPLE_002_S2_L001_R1_001.fastq.gz]
[[id:sample_003, organism:human, tissue:kidney, depth:45000000, quality:42.1, sample_num:3, lane:001, read:R1, chunk:001, priority:high], /Users/jonathan.manning/projects/training/side-quests/groovy_essentials/data/sequences/SAMPLE_003_S3_L001_R1_001.fastq.gz]
```

### 2.3. Dynamic Script Logic in Processes

Another place you'll find it very useful to break out your Groovy toolbox is in process script blocks. You can use Groovy logic to make your scripts dynamic and adaptable to different input conditions.

To illustrate what we mean, let's add some processes to our existing `main.nf` workflow that demonstrate common patterns for dynamic script generation. Open `modules/fastp.nf` and take a look:

```groovtitle="modules/fastp.nf" linenums="1"
process FASTP {
    container 'community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(sample_id), path("*_trimmed*.fastq.gz"), emit: reads

    script:
    """
    fastp \\
        --in1 ${reads[0]} \\
        --in2 ${reads[1]} \\
        --out1 ${meta.id}_trimmed_R1.fastq.gz \\
        --out2 ${meta.id}_trimmed_R2.fastq.gz \\
        --json ${meta.id}.fastp.json \\
        --html ${meta.id}.fastp.html \\
        --thread $task.cpus
    """
}
```

The process takes FASTQ files as input and runs the `fastp` tool to trim adapters and filter low-quality reads. Unfortunately, the person who wrote this process didn't allow for the single-end reads we have in our example dataset. Let's add it to our workflow and see what happens:

=== "After"

    ```groovy title="main.nf" linenums="1" hl_lines="1,31"
        include { FASTP } from './modules/fastp.nf'

        workflow {

            separateMetadata = { row ->
                def sample_meta = [
                    id: row.sample_id.toLowerCase(),
                    organism: row.organism,
                    tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                    depth: row.sequencing_depth.toInteger(),
                    quality: row.quality_score.toDouble()
                ]
                def fastq_path = file(row.file_path)

                def m = (fastq_path.name =~ /^(.+)_S(\d+)_L(\d{3})_(R[12])_(\d{3})\.fastq(?:\.gz)?$/)
                def file_meta = m ? [
                    sample_num: m[0][2].toInteger(),
                    lane: m[0][3],
                    read: m[0][4],
                    chunk: m[0][5]
                ] : [:]

                def priority = sample_meta.quality > 40 ? 'high' : 'normal'
                return [sample_meta + file_meta + [priority: priority], fastq_path]
            }

            ch_samples = Channel.fromPath("./data/samples.csv")
                .splitCsv(header: true)
                    .map(separateMetadata)

            ch_fastp = FASTP(ch_samples)
        }
    ```

=== "Before"

    ```groovy title="main.nf" linenums="1" hl_lines="4-24"
        workflow {

            separateMetadata = { row ->
                def sample_meta = [
                    id: row.sample_id.toLowerCase(),
                    organism: row.organism,
                    tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                    depth: row.sequencing_depth.toInteger(),
                    quality: row.quality_score.toDouble()
                ]
                def fastq_path = file(row.file_path)

                def m = (fastq_path.name =~ /^(.+)_S(\d+)_L(\d{3})_(R[12])_(\d{3})\.fastq(?:\.gz)?$/)
                def file_meta = m ? [
                    sample_num: m[0][2].toInteger(),
                    lane: m[0][3],
                    read: m[0][4],
                    chunk: m[0][5]
                ] : [:]

                def priority = sample_meta.quality > 40 ? 'high' : 'normal'
                return [sample_meta + file_meta + [priority: priority], fastq_path]
            }

            ch_samples = Channel.fromPath("./data/samples.csv")
                .splitCsv(header: true)
                    .map(separateMetadata)
                    .view()
        }
    ```

Run this modified workflow:

```bash title="Test fastp process"
nextflow run main.nf
```

You'll see a long error trace with some content like:

```console title="Process error"
ERROR ~ Error executing process > 'FASTP (3)'

Caused by:
  Process `FASTP (3)` terminated with an error exit status (255)


Command executed:

  fastp \
      --in1 SAMPLE_003_S3_L001_R1_001.fastq \
      --in2 null \
      --out1 sample_003_trimmed_R1.fastq.gz \
      --out2 sample_003_trimmed_R2.fastq.gz \
      --json sample_003.fastp.json \
      --html sample_003.fastp.html \
      --thread 2

Command exit status:
  255

Command output:
  (empty)
```

You can see that the process is trying to run `fastp` with a `null` value for the second input file, which is causing it to fail. This is because our dataset contains single-end reads, but the process is hardcoded to expect paired-end reads (two input files at a time).

Let's fix this by adding some Groovy logic to the `script:` block of the `FASTP` process to handle both single-end and paired-end reads dynamically. We'll use an if/else statement to check how many read files are are present and adjust the command accordingly.

=== "After"

    ```groovy title="main.nf" linenums="11" hl_lines="3,5,15"
        script:
        // Simple single-end vs paired-end detection
        def is_single = reads instanceof List ? reads.size() == 1 : true

        if (is_single) {
            def input_file = reads instanceof List ? reads[0] : reads
            """
            fastp \\
                --in1 ${input_file} \\
                --out1 ${meta.id}_trimmed.fastq.gz \\
                --json ${meta.id}.fastp.json \\
                --html ${meta.id}.fastp.html \\
                --thread $task.cpus
            """
        } else {
            """
            fastp \\
                --in1 ${reads[0]} \\
                --in2 ${reads[1]} \\
                --out1 ${meta.id}_trimmed_R1.fastq.gz \\
                --out2 ${meta.id}_trimmed_R2.fastq.gz \\
                --json ${meta.id}.fastp.json \\
                --html ${meta.id}.fastp.html \\
                --thread $task.cpus
            """
        }
    ```

=== "Before"

    ```groovy title="main.nf" linenums="11"
            script:
            """
            fastp \\
                --in1 ${reads[0]} \\
                --in2 ${reads[1]} \\
                --out1 ${meta.id}_trimmed_R1.fastq.gz \\
                --out2 ${meta.id}_trimmed_R2.fastq.gz \\
                --json ${meta.id}.fastp.json \\
                --html ${meta.id}.fastp.html \\
                --thread $task.cpus
            """
        }
    ```

Now the workflow can handle both single-end and paired-end reads gracefully. The Groovy logic checks the number of input files and constructs the appropriate command for `fastp`. Let's see if it works:

```bash title="Test dynamic fastp"
nextflow run main.nf
```

```console title="Successful run"
 N E X T F L O W   ~  version 25.04.6

Launching `main.nf` [adoring_rosalind] DSL2 - revision: 04b1cd93e9

executor >  local (3)
[31/a8ad4d] process > FASTP (3) [100%] 3 of 3 ✔
```

Looks good! If we check the actual commands that were run (customise for your task hash):

```console title="Check commands executed"
cat work/31/a8ad4d95749e685a6d842d3007957f/.command.sh
```

We can see that Nextflow correctly picked the right command for single-end reads:

```bash title=".command.sh"
#!/bin/bash -ue
fastp \
    --in1 SAMPLE_003_S3_L001_R1_001.fastq \
    --out1 sample_003_trimmed.fastq.gz \
    --json sample_003.fastp.json \
    --html sample_003.fastp.html \
    --thread 2
```

Another common one can be seen in [the Nextflow for Science Genomics module](../nf4science/genomics/02_joint_calling.md). In that module, the GATK process being called can take multiple input files, but each must be prefixed with `-V` to form a correct command line. The process uses Groovy logic to transform a collection of input files (`all_gvcfs`) into the correct command arguments:

```groovy title="command line manipulation for GATK" linenums="1"
    script:
    def gvcfs_line = all_gvcfs.collect { gvcf -> "-V ${gvcf}" }.join(' ')
    """
    gatk GenomicsDBImport \
        ${gvcfs_line} \
        -L ${interval_list} \
        --genomicsdb-workspace-path ${cohort_name}_gdb
    """
```

These patterns of using Groovy logic in process script blocks are extremely powerful and can be applied in many scenarios - from handling variable input types to building complex command-line arguments from file collections, making your processes truly adaptable to the diverse requirements of real-world data.

### Takeaway

In this section, you've learned:

- **Regular expressions for file parsing**: Using Groovy's `=~` operator and regex patterns to extract metadata from complex bioinformatics file naming conventions
- **Reusable closures**: Extracting complex logic into named closures that can be passed to channel operators, making workflows more readable and maintainable
- **Dynamic script generation**: Using Groovy conditional logic within process script blocks to adapt commands based on input characteristics (like single-end vs paired-end reads)
- **Command-line argument construction**: Transforming file collections into properly formatted command arguments using `collect()` and `join()` methods

These string processing patterns are essential for handling the diverse file formats and naming conventions you'll encounter in real-world bioinformatics workflows.


---

## 3. Conditional Logic and Process Control

Earlier on, we discussed how to use the `.map()` operator to use snippets of Groovy code to transform data flowing through channels. The counterpart to that is using Groovy to not just transform data, but to control which processes get executed based on the data itself. This is essential for building flexible workflows that can adapt to different sample types and analysis requirements.

Nextflow has several [operators](https://www.nextflow.io/docs/latest/reference/operator.html) that control process flow, including, many of which take closures as arguments, meanint their content is evaluated at run time, allowing us to use Groovy logic to drive workflow decisions based on channel content.

For example, let's pretend that our sequencing samples need to be trimmed with FASTP only if they're human samples with a coverage above a certain threshold. Mouse samples or low-coverage samples should be run with Trimgalore instead (this is a contrived example, but it illustrates the point).

Add a new process for Trimgalore in `modules/trimgalore.nf`:

=== "After"

    ```groovy title="main.nf" linenums="1" hl_lines="2"
    include { FASTP } from './modules/fastp.nf'
    include { TRIMGALORE } from './modules/trimgalore.nf'
    ```

=== "Before"

    ```groovy title="main.nf" linenums="1"
    include { FASTP } from './modules/fastp.nf'
    ```

... and then modify your `main.nf` workflow to branch samples based on their metadata and route them through the appropriate trimming process, like this:

=== "After"

    ```groovy title="branched workflow" linenums="28" hl_lines="5-12"
            ch_samples = Channel.fromPath("./data/samples.csv")
                .splitCsv(header: true)
                    .map(separateMetadata)

            trim_branches = ch_samples
                .branch { meta, reads ->
                    fastp: meta.organism == 'human' && meta.depth >= 30000000
                    trimgalore: true
                }

            ch_fastp = FASTP(trim_branches.fastp)
            ch_trimgalore = TRIMGALORE(trim_branches.trimgalore)

=== "Before"

    ```groovy title="branched workflow" linenums="28" hl_lines="5-12"
            ch_samples = Channel.fromPath("./data/samples.csv")
                .splitCsv(header: true)
                    .map(separateMetadata)

            ch_fastp = FASTP(ch_samples)


Run this modified workflow:

```bash title="Test conditional trimming"
nextflow run main.nf
```

```console title="Conditional trimming results"
 N E X T F L O W   ~  version 25.04.6

Launching `main.nf` [boring_koch] DSL2 - revision: 68a6bc7bd8

executor >  local (3)
[3d/bb1e90] process > FASTP (2)      [100%] 2 of 2 ✔
[4c/455334] process > TRIMGALORE (1) [100%] 1 of 1 ✔
```

Here, we've used small but mighty Groovy expressions inside the `.branch{}` operator to route samples based on their metadata. Human samples with high coverage go through `FASTP`, while all other samples go through `TRIMGALORE`. Combined with other closure-taking operators such as `.filter{}`, this allows us to build complex conditional workflows that adapt to the data itself.

### Takeaway

In this section, you've learned to use Groovy logic to control workflow execution with using the closure interfaces of Nextflow opearators.

Our pipeline now intelligently routes samples through appropriate processes, but production workflows need to handle invalid data gracefully. Let's add validation and error handling to make our pipeline robust enough for real-world use.

---

## 4. Error Handling and Validation Patterns

### 4.1. Basic Input Validation

Before our pipeline processes samples through complex conditional logic, we should validate that the input data meets our requirements. Let's create validation functions that check sample metadata and provide useful error messages:

=== "After"

    ```groovy title="main.nf" linenums="330" hl_lines="1-25"

    // Simple validation function
    def validateSample(Map sample) {
        def errors = []

        // Check required fields
        if (!sample.sample_id) {
            errors << "Missing sample_id"
        }

        if (!sample.organism) {
            errors << "Missing organism"
        }

        // Validate organism
        def valid_organisms = ['human', 'mouse', 'rat']
        if (sample.organism && !valid_organisms.contains(sample.organism.toLowerCase())) {
            errors << "Invalid organism: ${sample.organism}"
        }

        // Check sequencing depth is numeric
        if (sample.sequencing_depth) {
            try {
                def depth = sample.sequencing_depth as Integer
                if (depth < 1000000) {
                    errors << "Sequencing depth too low: ${depth}"
                }
            } catch (NumberFormatException e) {
                errors << "Invalid sequencing depth: ${sample.sequencing_depth}"
            }
        }

        return errors
    }

    // Test validation
    def test_samples = [
        [sample_id: 'SAMPLE_001', organism: 'human', sequencing_depth: '30000000'],
        [sample_id: '', organism: 'alien', sequencing_depth: 'invalid'],
        [sample_id: 'SAMPLE_003', organism: 'mouse', sequencing_depth: '500000']
    ]

    test_samples.each { sample ->
        def errors = validateSample(sample)
        if (errors) {
            println "Sample ${sample.sample_id}: ${errors.join(', ')}"
        } else {
            println "Sample ${sample.sample_id}: Valid"
        }
    }
    ```

=== "Before"

    ```groovy title="main.nf" linenums="330"
    ```

### 4.2. Try-Catch Error Handling

Let's implement simple try-catch patterns for handling errors:

=== "After"

    ```groovy title="main.nf" linenums="370" hl_lines="1-25"

    // Process sample with error handling
    def processSample(Map sample) {
        try {
            // Validate first
            def errors = validateSample(sample)
            if (errors) {
                throw new RuntimeException("Validation failed: ${errors.join(', ')}")
            }

            // Simulate processing
            def result = [
                id: sample.sample_id,
                organism: sample.organism,
                processed: true
            ]

            println "✓ Successfully processed ${sample.sample_id}"
            return result

        } catch (Exception e) {
            println "✗ Error processing ${sample.sample_id}: ${e.message}"

            // Return partial result
            return [
                id: sample.sample_id ?: 'unknown',
                organism: sample.organism ?: 'unknown',
                processed: false,
                error: e.message
            ]
        }
    }

    // Test error handling
    test_samples.each { sample ->
        def result = processSample(sample)
        println "Result for ${result.id}: processed = ${result.processed}"
    }
    ```

=== "Before"

    ```groovy title="main.nf" linenums="370"
    ```

### 4.3. Setting Defaults and Validation

Let's create a simple function that provides defaults and validates configuration:

=== "After"

    ```groovy title="main.nf" linenums="400" hl_lines="1-25"

    // Simple configuration with defaults
    def getConfig(Map user_params) {
        // Set defaults
        def defaults = [
            quality_threshold: 30,
            max_cpus: 4,
            output_dir: './results'
        ]

        // Merge user params with defaults
        def config = defaults + user_params

        // Simple validation
        if (config.quality_threshold < 0 || config.quality_threshold > 40) {
            println "Warning: Quality threshold ${config.quality_threshold} out of range, using default"
            config.quality_threshold = defaults.quality_threshold
        }

        if (config.max_cpus < 1) {
            println "Warning: Invalid CPU count ${config.max_cpus}, using default"
            config.max_cpus = defaults.max_cpus
        }

        return config
    }

    // Test configuration
    def test_configs = [
        [:], // Empty - should get defaults
        [quality_threshold: 35, max_cpus: 8], // Valid values
        [quality_threshold: -5, max_cpus: 0] // Invalid values
    ]

    test_configs.each { user_config ->
        def config = getConfig(user_config)
        println "Input: ${user_config} -> Output: ${config}"
    }
    ```

=== "Before"

    ```groovy title="main.nf" linenums="400"
    ```

### Takeaway

In this section, you've learned:

- **Basic validation functions** that check required fields and data types
- **Try-catch error handling** for graceful failure handling
- **Configuration with defaults** using map merging and validation

These patterns help you write workflows that handle invalid input gracefully and provide useful feedback to users.

Before diving into advanced closures, let's master some essential Groovy language features that make code more concise and null-safe. These operators and patterns are used throughout production Nextflow workflows and will make your code more robust and readable.

---

## 5. Essential Groovy Operators and Patterns

With our pipeline now handling complex conditional logic, we need to make it more robust against missing or malformed data. Bioinformatics workflows often deal with incomplete metadata, optional configuration parameters, and varying input formats. Let's enhance our pipeline with essential Groovy operators that handle these challenges gracefully.

### 5.1. Safe Navigation and Elvis Operators in Workflows

!!! note

    **Safe Navigation (`?.`) and Elvis (`?:`) Operators**: These are essential for null-safe programming. Safe navigation returns null instead of throwing an exception if the object is null, while the Elvis operator provides a default value if the left side is null, empty, or false.

The safe navigation operator (`?.`) and Elvis operator (`?:`) are essential for null-safe programming when processing real-world biological data:

- **Safe navigation (`?.`)** - returns null instead of throwing an exception if the object is null
- **Elvis operator (`?:`)** - provides a default value if the left side is null, empty, or false

=== "After"

    ```groovy title="main.nf" linenums="320" hl_lines="1-25"

    workflow {
        ch_samples = Channel.fromPath(params.input)
            .splitCsv(header: true)
            .map { row ->
                // Safe navigation prevents crashes on missing fields
                def sample_id = row.sample_id?.toLowerCase() ?: 'unknown_sample'
                def organism = row.organism?.toLowerCase() ?: 'unknown'

                // Elvis operator provides defaults
                def quality = (row.quality_score as Double) ?: 30.0
                def depth = (row.sequencing_depth as Integer) ?: 1_000_000

                // Chain operators for conditional defaults
                def reference = row.reference ?: (organism == 'human' ? 'GRCh38' : 'custom')

                // Groovy Truth - empty strings and nulls are false
                def priority = row.priority ?: (quality > 40 ? 'high' : 'normal')

                return [
                    id: sample_id,
                    organism: organism,
                    quality: quality,
                    depth: depth,
                    reference: reference,
                    priority: priority
                ]
            }
            .view { meta ->
                "Sample: ${meta.id} (${meta.organism}) - Quality: ${meta.quality}, Priority: ${meta.priority}"
            }
    }
    ```

=== "Before"

    ```groovy title="main.nf" linenums="320"
    ```

### 5.2. String Patterns and Multi-line Templates

Groovy provides powerful string features for parsing filenames and generating dynamic commands:

=== "After"

    ```groovy title="main.nf" linenums="370" hl_lines="1-30"

    workflow {
        // Demonstrate slashy strings for regex (no need to escape backslashes)
        def parseFilename = { filename ->
                    // Slashy string - compare to regular string: "^(\\w+)_(\\w+)_(\\d+)\\.fastq$"
        // Slashy strings don't require escaping backslashes, making regex patterns much cleaner
        def pattern = /^(\w+)_(\w+)_(\d+)\.fastq$/
            def matcher = filename =~ pattern

            if (matcher) {
                return [
                    organism: matcher[0][1].toLowerCase(),
                    tissue: matcher[0][2].toLowerCase(),
                    sample_id: matcher[0][3]
                ]
            } else {
                return [organism: 'unknown', tissue: 'unknown', sample_id: 'unknown']
            }
        }

        // Multi-line strings with interpolation for command generation
        def generateCommand = { meta ->
            def depth_category = meta.depth > 10_000_000 ? 'high' : 'standard'
            def db_path = meta.organism == 'human' ? '/db/human' : '/db/other'

            // Multi-line string with variable interpolation
            """
            echo "Processing ${meta.organism} sample: ${meta.sample_id}"
            analysis_tool \\
                --sample ${meta.sample_id} \\
                --depth-category ${depth_category} \\
                --database ${db_path} \\
                --threads ${params.max_cpus ?: 4}
            """
        }

        // Test the patterns
        ch_files = Channel.of('Human_Liver_001.fastq', 'Mouse_Brain_002.fastq')
            .map { filename ->
                def parsed = parseFilename(filename)
                def command = generateCommand([sample_id: parsed.sample_id, organism: parsed.organism, depth: 15_000_000])
                return [parsed, command]
            }
            .view { parsed, command -> "Parsed: ${parsed}, Command: ${command.split('\n')[0]}..." }
    }
    ```

=== "Before"

    ```groovy title="main.nf" linenums="370"
    ```

### 5.3. Combining Operators for Robust Data Handling

Let's combine these operators in a realistic workflow scenario:

=== "After"

    ```groovy title="main.nf" linenums="420" hl_lines="1-20"

    workflow {
        ch_samples = Channel.fromPath(params.input)
            .splitCsv(header: true)
            .map { row ->
                // Combine safe navigation and Elvis operators
                def meta = [
                    id: row.sample_id?.toLowerCase() ?: 'unknown',
                    organism: row.organism ?: 'unknown',
                    quality: (row.quality_score as Double) ?: 30.0,
                    files: row.file_path ? [file(row.file_path)] : []
                ]

                // Use Groovy Truth for validation
                if (meta.files && meta.id != 'unknown') {
                    return [meta, meta.files]
                } else {
                    log.info "Skipping sample with missing data: ${meta.id}"
                    return null
                }
            }
            .filter { it != null }  // Remove invalid samples using Groovy Truth
            .view { meta, files ->
                "Valid sample: ${meta.id} (${meta.organism}) - Quality: ${meta.quality}"
            }
    }
    ```

=== "Before"

    ```groovy title="main.nf" linenums="420"
    ```

### Takeaway

In this section, you've learned:

- **Safe navigation operator** (`?.`) for null-safe property access
- **Elvis operator** (`?:`) for default values and null coalescing
!!! note

    **Groovy Truth**: In Groovy, null, empty strings, empty collections, and zero are all considered "false" in boolean contexts. This is different from many other languages and is essential to understand for proper conditional logic.

- **Groovy Truth** - how null, empty strings, and empty collections evaluate to false - in Groovy, null, empty strings, empty collections, and zero are all considered "false" in boolean contexts
- **Slashy strings** (`/pattern/`) for regex patterns without escaping
- **Multi-line string interpolation** for command templates
- **Numeric literals with underscores** for improved readability

These patterns make your code more resilient to missing data and easier to read, which is essential when processing diverse bioinformatics datasets.

---

## 6. Advanced Closures and Functional Programming

Our pipeline now handles missing data gracefully and processes complex input formats robustly. But as our workflow grows more sophisticated, we start seeing repeated patterns in our data transformation code. Instead of copy-pasting similar closures across different processes or workflows, let's learn how to create reusable, composable functions that make our code cleaner and more maintainable.

### 6.1. Named Closures for Reusability

!!! note

    **Closures**: A closure is a block of code that can be assigned to a variable and executed later. Think of it as a function that can be passed around and reused. They're fundamental to Groovy's functional programming capabilities.

So far we've used anonymous closures defined inline within channel operations. When you find yourself repeating the same transformation logic across multiple processes or workflows, named closures can eliminate duplication and improve readability:

A **closure** is a block of code that can be assigned to a variable and executed later. Think of it as a function that can be passed around and reused.

=== "After"

    ```groovy title="main.nf" linenums="350" hl_lines="1-30"

    // Define reusable closures for common transformations
    def extractSampleInfo = { row ->
        [
            id: row.sample_id.toLowerCase(),
            organism: row.organism,
            quality: row.quality_score.toDouble(),
            depth: row.sequencing_depth.toInteger()
        ]
    }

    def addPriority = { meta ->
        meta + [priority: meta.quality > 40 ? 'high' : 'normal']
    }

    def formatForDisplay = { meta, file_path ->
        "Sample: ${meta.id} (${meta.organism}) - Quality: ${meta.quality}, Priority: ${meta.priority}"
    }

    workflow {
        // Use named closures in channel operations
        ch_samples = Channel.fromPath(params.input)
            .splitCsv(header: true)
            .map(extractSampleInfo)        // Named closure
            .map(addPriority)              // Named closure
            .map { meta -> [meta, file("./data/sequences/${meta.id}.fastq")] }
            .view(formatForDisplay)        // Named closure

        // Reuse the same closures elsewhere
        ch_filtered = ch_samples
            .filter { meta, file -> meta.quality > 30 }
            .map { meta, file -> addPriority(meta) }  // Reuse closure
            .view(formatForDisplay)                    // Reuse closure
    }
    ```

=== "Before"

    ```groovy title="main.nf" linenums="350"
    ```

### 6.2. Function Composition

Groovy closures can be composed together using the `>>` operator, allowing you to build complex transformations from simple, reusable pieces:

**Function composition** means chaining functions together so the output of one becomes the input of the next. The `>>` operator creates a new closure that applies multiple transformations in sequence.

=== "After"

    ```groovy title="main.nf" linenums="390" hl_lines="1-25"

    // Simple transformation closures
    def normalizeId = { meta ->
        meta + [id: meta.id.toLowerCase().replaceAll(/[^a-z0-9_]/, '_')]
    }

    def addQualityCategory = { meta ->
        def category = meta.quality > 40 ? 'excellent' :
                      meta.quality > 30 ? 'good' :
                      meta.quality > 20 ? 'acceptable' : 'poor'
        meta + [quality_category: category]
    }

    def addProcessingFlags = { meta ->
        meta + [
            needs_extra_qc: meta.quality < 30,
            high_priority: meta.organism == 'human' && meta.quality > 35
        ]
    }

    // Compose transformations using >> operator
    def enrichSample = normalizeId >> addQualityCategory >> addProcessingFlags

    workflow {
        Channel.fromPath(params.input)
            .splitCsv(header: true)
            .map(extractSampleInfo)
            .map(enrichSample)          // Apply composed transformation
            .view { meta ->
                "Processed: ${meta.id} (${meta.quality_category}) - Extra QC: ${meta.needs_extra_qc}"
            }
    }
    ```

=== "Before"

    ```groovy title="main.nf" linenums="390"
    ```

### 6.3. Currying for Specialized Functions

Currying allows you to create specialized versions of general-purpose closures by fixing some of their parameters:

**Currying** is a technique where you take a function with multiple parameters and create a new function with some of those parameters "fixed" or "pre-filled". This creates specialized versions of general-purpose functions.

=== "After"

    ```groovy title="main.nf" linenums="430" hl_lines="1-20"

    // General-purpose filtering closure
    def qualityFilter = { threshold, meta -> meta.quality >= threshold }

    // Create specialized filters using currying
    def highQualityFilter = qualityFilter.curry(40)
    def standardQualityFilter = qualityFilter.curry(30)

    workflow {
        ch_samples = Channel.fromPath(params.input)
            .splitCsv(header: true)
            .map(extractSampleInfo)

        // Use the specialized filters in different channel operations
        ch_high_quality = ch_samples.filter(highQualityFilter)
        ch_standard_quality = ch_samples.filter(standardQualityFilter)

        // Both channels can be processed differently
        ch_high_quality.view { meta -> "High quality: ${meta.id}" }
        ch_standard_quality.view { meta -> "Standard quality: ${meta.id}" }
    }
    ```

=== "Before"

    ```groovy title="main.nf" linenums="430"
    ```

### 6.4. Closures Accessing External Variables

Closures can access and modify variables from their defining scope, which is useful for collecting statistics:

=== "After"

    ```groovy title="main.nf" linenums="480" hl_lines="1-20"

    workflow {
        // Variable in the workflow scope
        def sample_count = 0
        def human_samples = 0

        // Closure that accesses and modifies external variables
        def countSamples = { meta ->
            sample_count++  // Modifies external variable
            if (meta.organism == 'human') {
                human_samples++  // Modifies another external variable
            }
            return meta  // Pass data through unchanged
        }

        Channel.fromPath(params.input)
            .splitCsv(header: true)
            .map(extractSampleInfo)
            .map(countSamples)          // Closure modifies external variables
            .collect()                  // Wait for all samples to be processed
            .view {
                "Processing complete: ${sample_count} total samples, ${human_samples} human samples"
            }
    }
    ```

=== "Before"

    ```groovy title="main.nf" linenums="480"
    ```

### Takeaway

In this section, you've learned:

- **Named closures** for eliminating code duplication and improving readability
- **Function composition** with `>>` operator to build complex transformations
- **Currying** to create specialized versions of general-purpose closures
- **Variable scope access** in closures for collecting statistics and generating reports

These advanced patterns help you write more maintainable, reusable workflows that follow functional programming principles while remaining easy to understand and debug.

With our pipeline now capable of intelligent routing, robust error handling, and advanced functional programming patterns, we're ready for the final enhancement. As your workflows scale to process hundreds or thousands of samples, you'll need sophisticated data processing capabilities that can organize, filter, and analyze large collections efficiently.

The functional programming patterns we just learned work beautifully with Groovy's powerful collection methods. Instead of writing loops and conditional logic, you can chain together expressive operations that clearly describe what you want to accomplish.

---

## 7. Collection Operations and File Path Manipulations

### 7.1. Common Collection Methods in Channel Operations

When processing large datasets, channel operations often need to organize and analyze sample collections. Groovy's collection methods integrate seamlessly with Nextflow channels to provide powerful data processing capabilities:

Groovy provides many built-in methods for working with collections (lists, maps, etc.) that make data processing much more expressive than traditional loops.

=== "After"

    ```groovy title="main.nf" linenums="500" hl_lines="1-40"

    // Sample data with mixed quality and organisms
    def samples = [
        [id: 'sample_001', organism: 'human', quality: 42, files: ['data1.txt', 'data2.txt']],
        [id: 'sample_002', organism: 'mouse', quality: 28, files: ['data3.txt']],
        [id: 'sample_003', organism: 'human', quality: 35, files: ['data4.txt', 'data5.txt', 'data6.txt']],
        [id: 'sample_004', organism: 'rat', quality: 45, files: ['data7.txt']],
        [id: 'sample_005', organism: 'human', quality: 30, files: ['data8.txt', 'data9.txt']]
    ]

    // findAll - filter collections based on conditions
    def high_quality_samples = samples.findAll { it.quality > 40 }
    println "High quality samples: ${high_quality_samples.collect { it.id }.join(', ')}"

    // groupBy - group samples by organism
    def samples_by_organism = samples.groupBy { it.organism }
    println "Grouping by organism:"
    samples_by_organism.each { organism, sample_list ->
        println "  ${organism}: ${sample_list.size()} samples"
    }

    // unique - get unique organisms
    def organisms = samples.collect { it.organism }.unique()
    println "Unique organisms: ${organisms.join(', ')}"

    // flatten - flatten nested file lists
    def all_files = samples.collect { it.files }.flatten()
    println "All files: ${all_files.take(5).join(', ')}... (${all_files.size()} total)"

    // sort - sort samples by quality
    def sorted_by_quality = samples.sort { it.quality }
    println "Quality range: ${sorted_by_quality.first().quality} to ${sorted_by_quality.last().quality}"

    // reverse - reverse the order
    def reverse_quality = samples.sort { it.quality }.reverse()
    println "Highest quality first: ${reverse_quality.collect { "${it.id}(${it.quality})" }.join(', ')}"

    // count - count items matching condition
    def human_samples = samples.count { it.organism == 'human' }
    println "Human samples: ${human_samples} out of ${samples.size()}"

    // any/every - check conditions across collection
    def has_high_quality = samples.any { it.quality > 40 }
    def all_have_files = samples.every { it.files.size() > 0 }
    println "Has high quality samples: ${has_high_quality}"
    println "All samples have files: ${all_have_files}"
    ```

=== "Before"

    ```groovy title="main.nf" linenums="500"
    ```

### 7.2. File Path Manipulations

Working with file paths is essential in bioinformatics workflows. Groovy provides many useful methods for extracting information from file paths:

=== "After"

    ```groovy title="main.nf" linenums="550" hl_lines="1-30"

    // File path manipulation examples
    def sample_files = [
        '/path/to/data/patient_001_R1.fastq.gz',
        '/path/to/data/patient_001_R2.fastq.gz',
        '/path/to/results/patient_002_analysis.bam',
        '/path/to/configs/experiment_setup.json'
    ]

    sample_files.each { file_path ->
        def f = file(file_path)  // Create Nextflow file object

        println "\nFile: ${file_path}"
        println "  Name: ${f.getName()}"                    // Just filename
        println "  BaseName: ${f.getBaseName()}"            // Filename without extension
        println "  Extension: ${f.getExtension()}"          // File extension
        println "  Parent: ${f.getParent()}"                // Parent directory
        println "  Parent name: ${f.getParent().getName()}" // Just parent directory name

        // Extract sample ID from filename
        def matcher = f.getName() =~ /^(patient_\d+)/
        if (matcher) {
            println "  Sample ID: ${matcher[0][1]}"
        }
    }

    // Group files by sample ID using path manipulation
    def files_by_sample = sample_files
        .findAll { it.contains('patient') }  // Only patient files
        .groupBy { file_path ->
            def filename = file(file_path).getName()
            def matcher = filename =~ /^(patient_\d+)/
            return matcher ? matcher[0][1] : 'unknown'
        }

    println "\nFiles grouped by sample:"
    files_by_sample.each { sample_id, files ->
        println "  ${sample_id}: ${files.size()} files"
    }
    ```

=== "Before"

    ```groovy title="main.nf" linenums="550"
    ```

### 7.3. The Spread Operator

The spread operator (`*.`) is a powerful Groovy feature for calling methods on all elements in a collection:

The **spread operator** (`*.`) is a shorthand way to call the same method on every element in a collection. It's equivalent to using `.collect { it.methodName() }` but more concise.

=== "After"

    ```groovy title="main.nf" linenums="590" hl_lines="1-20"

    // Spread operator examples
    def file_paths = [
        '/data/sample1.fastq',
        '/data/sample2.fastq',
        '/results/output1.bam',
        '/results/output2.bam'
    ]

    // Convert to file objects
    def files = file_paths.collect { file(it) }

    // Using spread operator - equivalent to files.collect { it.getName() }
    def filenames = files*.getName()
    println "Filenames: ${filenames.join(', ')}"

    // Get all parent directories
    def parent_dirs = files*.getParent()*.getName()
    println "Parent directories: ${parent_dirs.unique().join(', ')}"

    // Get all extensions
    def extensions = files*.getExtension().unique()
    println "File types: ${extensions.join(', ')}"
    ```

=== "Before"

    ```groovy title="main.nf" linenums="590"
    ```

### Takeaway

In this section, you've learned:

- **Collection filtering** with `findAll` and conditional logic
- **Grouping and organizing** data with `groupBy` and `sort`
- **File path manipulation** using Nextflow's file object methods
- **Spread operator** (`*.`) for concise collection operations

These patterns help you process and organize complex datasets efficiently, which is essential for handling real-world bioinformatics data.

---

## Summary

Throughout this side quest, you've built a comprehensive sample processing pipeline that evolved from basic metadata handling to a sophisticated, production-ready workflow. Each section built upon the previous, demonstrating how Groovy transforms simple Nextflow workflows into powerful data processing systems.

Here's how we progressively enhanced our pipeline:

1. **Nextflow vs Groovy Boundaries**: You learned to distinguish between workflow orchestration (Nextflow) and programming logic (Groovy), including the crucial differences between constructs like `collect`.

2. **String Processing**: You learned regular expressions, parsing functions, and file collection transformation for building dynamic command-line arguments.

3. **Conditional Logic**: You added intelligent routing that automatically selects analysis strategies based on sample characteristics like organism, quality scores, and sequencing depth.

4. **Error Handling**: You made the pipeline robust by adding validation functions, try-catch error handling, and configuration management with sensible defaults.

5. **Essential Groovy Operators**: You mastered safe navigation (`?.`), Elvis (`?:`), Groovy Truth, slashy strings, and other key language features that make code more resilient and readable.

6. **Advanced Closures**: You learned functional programming techniques including named closures, function composition, currying, and closures with variable scope access for building reusable, maintainable code.

7. **Collection Operations**: You added sophisticated data processing capabilities using Groovy collection methods like `findAll`, `groupBy`, `unique`, `flatten`, and the spread operator to handle large-scale sample processing.

### Key Benefits

- **Clearer code**: Understanding when to use Nextflow and Groovy helps you write more organized workflows
- **Better error handling**: Basic validation and try-catch patterns help your workflows handle problems gracefully
- **Flexible processing**: Conditional logic lets your workflows process different sample types appropriately
- **Configuration management**: Using defaults and simple validation makes your workflows easier to use

### From Simple to Sophisticated

The pipeline journey you completed demonstrates the evolution from basic data processing to production-ready bioinformatics workflows:

1. **Started simple**: Basic CSV processing and metadata extraction with clear Nextflow vs Groovy boundaries
2. **Added intelligence**: Dynamic file name parsing with regex patterns and conditional routing based on sample characteristics
3. **Made it robust**: Null-safe operators, validation, error handling, and graceful failure management
4. **Made it maintainable**: Advanced closure patterns, function composition, and reusable components that eliminate code duplication
5. **Scaled it efficiently**: Collection operations for processing hundreds of samples with powerful data filtering and organization

This progression mirrors the real-world evolution of bioinformatics pipelines - from research prototypes handling a few samples to production systems processing thousands of samples across laboratories and institutions. Every challenge you solved and pattern you learned reflects actual problems developers face when scaling Nextflow workflows.

### Next Steps

With these Groovy fundamentals mastered, you're ready to:

- Write cleaner workflows with proper separation between Nextflow and Groovy logic
- Transform file collections into properly formatted command-line arguments
- Handle different file naming conventions and input formats gracefully
- Build reusable, maintainable code using advanced closure patterns and functional programming
- Process and organize complex datasets using collection operations
- Add basic validation and error handling to make your workflows more user-friendly

Continue practicing these patterns in your own workflows, and refer to the [Groovy documentation](http://groovy-lang.org/documentation.html) when you need to explore more advanced features.

### Key Concepts Reference

- **Language Boundaries**
  ```groovy title="Nextflow vs Groovy examples"
  // Nextflow: workflow orchestration
  Channel.fromPath('*.fastq').splitCsv(header: true)

  // Groovy: data processing
  sample_data.collect { it.toUpperCase() }
  ```

- **String Processing**
  ```groovy title="String processing examples"
  // Pattern matching
  filename =~ ~/^(\w+)_(\w+)_(\d+)\.fastq$/

  // Function with conditional return
  def parseSample(filename) {
      def matcher = filename =~ pattern
      return matcher ? [valid: true, data: matcher[0]] : [valid: false]
  }

  // File collection to command arguments (in process script block)
  script:
  def file_args = input_files.collect { file -> "--input ${file}" }.join(' ')
  """
  analysis_tool ${file_args} --output results.txt
  """
  ```

- **Error Handling**
  ```groovy title="Error handling patterns"
  try {
      def errors = validateSample(sample)
      if (errors) throw new RuntimeException("Invalid: ${errors.join(', ')}")
  } catch (Exception e) {
      println "Error: ${e.message}"
  }
  ```

- **Essential Groovy Operators**
  ```groovy title="Essential operators examples"
  // Safe navigation and Elvis operators
  def id = data?.sample?.id ?: 'unknown'
  if (sample.files) println "Has files"  // Groovy Truth

  // Slashy strings for regex
  def pattern = /^\w+_R[12]\.fastq$/
  def script = """
  echo "Processing ${sample.id}"
  analysis --depth ${depth ?: 1_000_000}
  """
  ```

- **Advanced Closures**
  ```groovy title="Advanced closure patterns"
  // Named closures and composition
  def enrichData = normalizeId >> addQualityCategory >> addFlags
  def processor = generalFunction.curry(fixedParam)

  // Closures with scope access
  def collectStats = { data -> stats.count++; return data }
  ```

- **Collection Operations**
  ```groovy title="Collection operations examples"
  // Filter, group, and organize data
  def high_quality = samples.findAll { it.quality > 40 }
  def by_organism = samples.groupBy { it.organism }
  def file_names = files*.getName()  // Spread operator
  def all_files = nested_lists.flatten()
  ```

## Resources

- [Groovy Documentation](http://groovy-lang.org/documentation.html)
- [Nextflow Operators](https://www.nextflow.io/docs/latest/operator.html)
- [Regular Expressions in Groovy](https://groovy-lang.org/syntax.html#_regular_expression_operators)
- [JSON Processing](https://groovy-lang.org/json.html)
- [XML Processing](https://groovy-lang.org/processing-xml.html)
