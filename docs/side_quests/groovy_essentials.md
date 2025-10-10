# Groovy Essentials for Nextflow Developers

Nextflow is built on Groovy, a powerful dynamic language that runs on the Java Virtual Machine. You can write a lot of Nextflow without ever feeling like you've learned Groovy - many workflows use only basic syntax for variables, maps, and lists. Most Nextflow tutorials focus on workflow orchestration (channels, processes, and data flow), and you can go surprisingly far with just that.

However, when you need to manipulate data, parse complex filenames, implement conditional logic, or build robust production workflows, you're writing Groovy code - and knowing a few key Groovy concepts can dramatically improve your ability to solve real-world problems efficiently. Understanding where Nextflow ends and Groovy begins helps you write clearer, more maintainable workflows.

This side quest takes you on a hands-on journey from basic concepts to production-ready patterns. We'll transform a simple CSV-reading workflow into a sophisticated bioinformatics pipeline, evolving it step-by-step through realistic challenges:

- **Understanding boundaries:** Distinguish between Nextflow operators and Groovy methods, and master when to use each
- **Data manipulation:** Extract, transform, and subset maps and collections using Groovy's powerful operators
- **String processing:** Parse complex file naming schemes with regex patterns and master variable interpolation
- **Reusable functions:** Extract complex logic into named functions for cleaner, more maintainable workflows
- **Dynamic logic:** Build processes that adapt to different input types and use closures for dynamic resource allocation
- **Conditional routing:** Intelligently route samples through different processes based on their metadata characteristics
- **Safe operations:** Handle missing data gracefully with null-safe operators and validate inputs with clear error messages
- **Configuration-based handlers:** Use workflow event handlers for logging, notifications, and lifecycle management

---

## 0. Warmup

### 0.1. Prerequisites

Before taking on this side quest you should:

- Complete the [Hello Nextflow](../hello_nextflow/README.md) tutorial or have equivalent experience
- Understand basic Nextflow concepts (processes, channels, workflows)
- Have basic familiarity with common programming constructs used in Groovy syntax (variables, maps, lists)

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
├── collect.nf
├── data
│   ├── samples.csv
│   └── sequences
│       ├── SAMPLE_001_S1_L001_R1_001.fastq
│       ├── SAMPLE_002_S2_L001_R1_001.fastq
│       └── SAMPLE_003_S3_L001_R1_001.fastq
├── main.nf
├── modules
│   ├── fastp.nf
│   ├── generate_report.nf
│   └── trimgalore.nf
└── nextflow.config

4 directories, 10 files
```

Our sample CSV contains information about biological samples that need different processing based on their characteristics:

```console title="samples.csv"
sample_id,organism,tissue_type,sequencing_depth,file_path,quality_score
SAMPLE_001,human,liver,30000000,data/sequences/SAMPLE_001_S1_L001_R1_001.fastq,38.5
SAMPLE_002,mouse,brain,25000000,data/sequences/SAMPLE_002_S2_L001_R1_001.fastq,35.2
SAMPLE_003,human,kidney,45000000,data/sequences/SAMPLE_003_S3_L001_R1_001.fastq,42.1
```

We'll use this realistic dataset to explore practical Groovy techniques that you'll encounter in real bioinformatics workflows.

---

## 1. Nextflow vs Groovy: Understanding the Boundaries

### 1.1. Identifying What's What

One of the most common sources of confusion for Nextflow developers is understanding when they're working with Nextflow constructs versus Groovy language features. Let's build a workflow step by step to see a common example of how they work together.

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
Launching `main.nf` [marvelous_tuckerman] DSL2 - revision: 6113e05c17

[sample_id:SAMPLE_001, organism:human, tissue_type:liver, sequencing_depth:30000000, file_path:data/sequences/SAMPLE_001_S1_L001_R1_001.fastq, quality_score:38.5]
[sample_id:SAMPLE_002, organism:mouse, tissue_type:brain, sequencing_depth:25000000, file_path:data/sequences/SAMPLE_002_S2_L001_R1_001.fastq, quality_score:35.2]
[sample_id:SAMPLE_003, organism:human, tissue_type:kidney, sequencing_depth:45000000, file_path:data/sequences/SAMPLE_003_S3_L001_R1_001.fastq, quality_score:42.1]
```

#### Step 2: Adding the Map Operator

Now we're going to use some Groovy code to transform the data, using the `.map()` operator you will probably already be familiar with. This operator takes a 'closure' where we can write Groovy code to transform each item.

!!! note

    A **closure** is a block of code that can be passed around and executed later. Think of it as a function that you define inline. In Groovy, closures are written with curly braces `{ }` and can take parameters. They're fundamental to how Nextflow operators work and if you've been writing Nextflow for a while, you may already have been using them without realizing it!

Here's what that map operation looks like:

=== "After"

    ```groovy title="main.nf" linenums="2" hl_lines="3-6"
        ch_samples = Channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row ->
                return row
            }
            .view()
    ```

=== "Before"

    ```groovy title="main.nf" linenums="2" hl_lines="3"
        ch_samples = Channel.fromPath("./data/samples.csv")
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

Now we're going to write **pure Groovy code** inside our closure. Everything from this point forward in this section is Groovy syntax and methods, not Nextflow operators.

=== "After"

    ```groovy title="main.nf" linenums="2" hl_lines="4-12"
        ch_samples = Channel.fromPath("./data/samples.csv")
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
        ch_samples = Channel.fromPath("./data/samples.csv")
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
                return sample_meta + [priority: priority]
            }
            .view()
    ```

=== "Before"

    ```groovy title="main.nf" linenums="2" hl_lines="11"
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
                return sample_meta
            }
            .view()
    ```

The ternary operator is a shorthand for an if/else statement that follows the pattern `condition ? value_if_true : value_if_false`. This line means: "If the quality is greater than 40, use 'high', otherwise use 'normal'". Its cousin, the **Elvis operator** (`?:`), provides default values when something is null or empty - we'll explore that pattern later in this tutorial.

The map addition operator `+` creates a **new map** rather than modifying the existing one. This line creates a new map that contains all the key-value pairs from `sample_meta` plus the new `priority` key.

!!! Note

    Never modify maps passed into closures - always create new ones using `+` (for example). In Nextflow, the same data often flows through multiple operations simultaneously. Modifying a map in-place can cause unpredictable side effects when other operations reference that same object. Creating new maps ensures each operation has its own clean copy.

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

#### Step 4.5: Subsetting Maps with `.subMap()`

While the `+` operator adds keys to a map, sometimes you need to do the opposite - extract only specific keys. Groovy's `.subMap()` method is perfect for this.

Let's add a line to create a simplified version of our metadata that only contains identification fields:

=== "After"

    ```groovy title="main.nf" linenums="2" hl_lines="12-15"
        ch_samples = Channel.fromPath("./data/samples.csv")
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
                def id_only = sample_meta.subMap(['id', 'organism', 'tissue'])
                println "ID fields only: ${id_only}"

                def priority = sample_meta.quality > 40 ? 'high' : 'normal'
                return sample_meta + [priority: priority]
            }
            .view()
    ```

=== "Before"

    ```groovy title="main.nf" linenums="2" hl_lines="12"
        ch_samples = Channel.fromPath("./data/samples.csv")
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
                def priority = sample_meta.quality > 40 ? 'high' : 'normal'
                return sample_meta + [priority: priority]
            }
            .view()
    ```

Run the modified workflow:

```bash title="Test subMap"
nextflow run main.nf
```

You should see output showing both the full metadata displayed by the `view()` operation and the extracted subset we printed with `println`:

```console title="SubMap results"
 N E X T F L O W   ~  version 25.04.6

Launching `main.nf` [peaceful_cori] DSL2 - revision: 4cc4a8340f

ID fields only: [id:sample_001, organism:human, tissue:liver]
ID fields only: [id:sample_002, organism:mouse, tissue:brain]
ID fields only: [id:sample_003, organism:human, tissue:kidney]
[id:sample_001, organism:human, tissue:liver, depth:30000000, quality:38.5, priority:normal]
[id:sample_002, organism:mouse, tissue:brain, depth:25000000, quality:35.2, priority:normal]
[id:sample_003, organism:human, tissue:kidney, depth:45000000, quality:42.1, priority:high]
```

The `.subMap()` method takes a list of keys and returns a new map containing only those keys. If a key doesn't exist in the original map, it's simply not included in the result.

This is particularly useful when you need to create different metadata versions for different processes - some might need full metadata while others need only minimal identification fields.

Now remove those println statements to restore your workflow to its previous state, as we don't need them going forward.

!!! tip "Map Operations Summary"

    - **Add keys**: `map1 + [new_key: value]` - Creates new map with additional keys
    - **Extract keys**: `map1.subMap(['key1', 'key2'])` - Creates new map with only specified keys
    - **Both operations create new maps** - Original maps remain unchanged

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
                return [sample_meta + [priority: priority], file(row.file_path) ]
            }
            .view()
    ```

=== "Before"

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
[[id:sample_001, organism:human, tissue:liver, depth:30000000, quality:38.5, priority:normal], /workspaces/training/side-quests/groovy_essentials/data/sequences/SAMPLE_001_S1_L001_R1_001.fastq]
[[id:sample_002, organism:mouse, tissue:brain, depth:25000000, quality:35.2, priority:normal], /workspaces/training/side-quests/groovy_essentials/data/sequences/SAMPLE_002_S2_L001_R1_001.fastq]
[[id:sample_003, organism:human, tissue:kidney, depth:45000000, quality:42.1, priority:high], /workspaces/training/side-quests/groovy_essentials/data/sequences/SAMPLE_003_S3_L001_R1_001.fastq]
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

We are:

- Defining a (Groovy) list
- Using the `fromList()` channel factory to create a channel that emits each sample ID as a separate item
- Using `view()` to print each item as it flows through the channel
- Applying Nextflow's `collect()` operator to gather all items into a single list
- Using a second `view()` to print the collected result which appears as a single item containing a list of all sample IDs

We've changed the structure of the channel, but we haven't changed the data itself.

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

`view()` returns an output for every channel emission, so we know that this single output contains all 3 original items grouped into one list.

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

In this new snippet we:

- Define a new variable `formatted_ids` that uses Groovy's `collect` method to transform each sample ID in the original list
- Print the result using `println`

Run the modified workflow:

```bash title="Test Groovy collect"
nextflow run collect.nf
```

```console title="Groovy collect results" hl_lines="5"
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

### 1.3. The Spread Operator (`*.`) - Shorthand for Property Extraction

Related to Groovy's `collect` is the spread operator (`*.`), which provides a concise way to extract properties from collections. It's essentially syntactic sugar for a common `collect` pattern.

Let's add a demonstration to our `collect.nf` file:

=== "After"

    ```groovy title="collect.nf" linenums="1" hl_lines="15-18"
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

    // Spread operator - concise property access
    def sample_data = [[id: 's1', quality: 38.5], [id: 's2', quality: 42.1], [id: 's3', quality: 35.2]]
    def all_ids = sample_data*.id
    println "Spread operator result: ${all_ids}"
    ```

=== "Before"

    ```groovy title="collect.nf" linenums="1"
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

Run the updated workflow:

```bash title="Test spread operator"
nextflow run collect.nf
```

You should see output like:

```console title="Spread operator output" hl_lines="6"
 N E X T F L O W   ~  version 25.04.6

Launching `collect.nf` [cranky_galileo] DSL2 - revision: 5f3c8b2a91

Groovy collect result: [SPECIMEN_001, SPECIMEN_002, SPECIMEN_003] (3 items transformed into 3)
Spread operator result: [s1, s2, s3]
Individual channel item: sample_001
Individual channel item: sample_002
Individual channel item: sample_003
Nextflow collect() result: [sample_001, sample_002, sample_003] (3 items grouped into 1)
```

The spread operator `*.` is a shorthand for a common collect pattern:

```groovy
// These are equivalent:
def ids = samples*.id
def ids = samples.collect { it.id }

// Also works with method calls:
def names = files*.getName()
def names = files.collect { it.getName() }
```

The spread operator is particularly useful when you need to extract a single property from a list of objects - it's more readable than writing out the full `collect` closure.

!!! tip "When to Use Groovy's Spread vs Collect"

    - **Use spread (`*.`)** for simple property access: `samples*.id`, `files*.name`
    - **Use collect** for transformations or complex logic: `samples.collect { it.id.toUpperCase() }`, `samples.collect { [it.id, it.quality > 40] }`

### Takeaway

In this section, you've learned:

- **It takes both Nextflow and Groovy**: Nextflow provides the workflow structure and data flow, while Groovy provides the data manipulation and logic
- **Distinguishing Nextflow from Groovy**: How to identify which language construct you're using given the context
- **Context matters**: The same operation name can have completely different behaviors

Understanding these boundaries is essential for debugging, documentation, and writing maintainable workflows.

Next we'll dive deeper into Groovy's powerful string processing capabilities, which are essential for handling real-world data.

---

## 2. String Processing and Dynamic Script Generation

The difference between a brittle workflow that breaks on unexpected input and a robust pipeline that adapts gracefully often comes down to mastering Groovy's string processing capabilities. In this section, we'll explore how to parse complex file names, generate process scripts dynamically based on input characteristics, and properly interpolate variables in different contexts.

### 2.1. Pattern Matching and Regular Expressions

Many bioinformatics workflows encounter files with complex naming conventions that encode important metadata. Let's see how Groovy's pattern matching can extract this information automatically.

We're going to return to our `main.nf` workflow and add some pattern matching logic to extract additional sample information from file names. The FASTQ files in our dataset follow Illumina-style naming conventions with names like `SAMPLE_001_S1_L001_R1_001.fastq.gz`. These might look cryptic, but they actually encode useful metadata like sample ID, lane number, and read direction. We're going to use Groovy's regex capabilities to parse these names.

Make the following change to your existing `main.nf` workflow:

=== "After"

    ```groovy title="main.nf" linenums="4" hl_lines="10-21"
            .map { row ->
                // This is all Groovy code now!
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

    ```groovy title="main.nf" linenums="4" "hl_lines="11"
            .map { row ->
                // This is all Groovy code now!
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

[[id:sample_001, organism:human, tissue:liver, depth:30000000, quality:38.5, sample_num:1, lane:001, read:R1, chunk:001, priority:normal], /workspaces/training/side-quests/groovy_essentials/data/sequences/SAMPLE_001_S1_L001_R1_001.fastq]
[[id:sample_002, organism:mouse, tissue:brain, depth:25000000, quality:35.2, sample_num:2, lane:001, read:R1, chunk:001, priority:normal], /workspaces/training/side-quests/groovy_essentials/data/sequences/SAMPLE_002_S2_L001_R1_001.fastq]
[[id:sample_003, organism:human, tissue:kidney, depth:45000000, quality:42.1, sample_num:3, lane:001, read:R1, chunk:001, priority:high], /workspaces/training/side-quests/groovy_essentials/data/sequences/SAMPLE_003_S3_L001_R1_001.fastq]
```

### 2.2. Dynamic Script Generation in Processes

Process script blocks are essentially multi-line strings that get passed to the shell. You can use Groovy logic to dynamically generate different script strings based on input characteristics, making your processes adaptable to different input conditions.

To illustrate what we mean, let's add some processes to our existing `main.nf` workflow that demonstrate common patterns for dynamic script generation. Open `modules/fastp.nf` and take a look:

```groovy title="modules/fastp.nf" linenums="1"
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

First, include the module at the very first line of your `main.nf` workflow:

```groovy title="main.nf" linenums="1"
include { FASTP } from './modules/fastp.nf'
```

Then modify the `workflow` block to connect the `ch_samples` channel to the `FASTP` process:

=== "After"

    ```groovy title="main.nf" linenums="25" hl_lines="7"
    workflow {

        ch_samples = Channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map{ row -> separateMetadata(row) }

        ch_fastp = FASTP(ch_samples)
    }
    ```

=== "Before"

    ```groovy title="main.nf" linenums="25" hl_lines="6"
    workflow {

        ch_samples = Channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map{ row -> separateMetadata(row) }
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

    ```groovy title="main.nf" linenums="10" hl_lines="3 5 15"
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

    ```groovy title="main.nf" linenums="10"
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

Another common usage of dynamic script logic can be seen in [the Nextflow for Science Genomics module](../../nf4science/genomics/02_joint_calling). In that module, the GATK process being called can take multiple input files, but each must be prefixed with `-V` to form a correct command line. The process uses Groovy logic to transform a collection of input files (`all_gvcfs`) into the correct command arguments:

```groovy title="command line manipulation for GATK" linenums="1" hl_lines="2 5"
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

### 2.3. Variable Interpolation: Groovy, Bash, and Shell Variables

When writing process scripts, you're actually working with three different types of variables, and using the wrong syntax is a common source of errors. Let's add a process that creates a processing report to demonstrate the differences.

Take a look a the module file `modules/generate_report.nf`:

```groovy title="modules/generate_report.nf" linenums="1"
process GENERATE_REPORT {

    publishDir 'results/reports', mode: 'copy'

    input:
    tuple val(meta), path(reads)

    output:
    path "${meta.id}_report.txt"

    script:
    """
    echo "Processing ${reads}" > ${meta.id}_report.txt
    echo "Sample: ${meta.id}" >> ${meta.id}_report.txt
    """
}
```

This process writes a simple report with the sample ID and filename. Now let's run it to see what happens when we need to mix different types of variables.

Include the process in your `main.nf` and add it to the workflow:

=== "After"

    ```groovy title="main.nf" linenums="1" hl_lines="2 12"
    include { FASTP } from './modules/fastp.nf'
    include { GENERATE_REPORT } from './modules/generate_report.nf'

    // ... separateMetadata function ...

    workflow {
        ch_samples = Channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map{ row -> separateMetadata(row) }

        ch_fastp = FASTP(ch_samples)
        GENERATE_REPORT(ch_samples)
    }
    ```

=== "Before"

    ```groovy title="main.nf" linenums="1" hl_lines="1 10"
    include { FASTP } from './modules/fastp.nf'

    // ... separateMetadata function ...

    workflow {
        ch_samples = Channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map{ row -> separateMetadata(row) }

        ch_fastp = FASTP(ch_samples)
    }
    ```

Now run the workflow and check the generated reports in `results/reports/`. They should contain basic information about each sample.

But what if we want to add information about when and where the processing occurred? Let's modify the process to use **shell** variables and a bit of command substitution to include the current user, hostname, and date in the report:

=== "After"

    ```groovy title="modules/generate_report.nf" linenums="10" hl_lines="5-7"
    script:
    """
    echo "Processing ${reads}" > ${meta.id}_report.txt
    echo "Sample: ${meta.id}" >> ${meta.id}_report.txt
    echo "Processed by: ${USER}" >> ${meta.id}_report.txt
    echo "Hostname: $(hostname)" >> ${meta.id}_report.txt
    echo "Date: $(date)" >> ${meta.id}_report.txt
    """
    ```

=== "Before"

    ```groovy title="modules/generate_report.nf" linenums="10"
    script:
    """
    echo "Processing ${reads}" > ${meta.id}_report.txt
    echo "Sample: ${meta.id}" >> ${meta.id}_report.txt
    """
    ```

If you run this, you'll notice an error or unexpected behavior - Nextflow tries to interpret `$(hostname)` as a Groovy variable that doesn't exist:

```console title="Error with shell variables"
unknown recognition error type: groovyjarjarantlr4.v4.runtime.LexerNoViableAltException
ERROR ~ Module compilation error
- file : /workspaces/training/side-quests/groovy_essentials/modules/generate_report.nf
- cause: token recognition error at: '(' @ line 16, column 22.
       echo "Hostname: $(hostname)" >> ${meta.id}_report.txt
                        ^

1 error
```

 We need to escape it so Bash can handle it instead.

Fix this by escaping the shell variables and command substitutions with a backslash (`\`):

=== "After - Fixed"

    ```groovy title="modules/generate_report.nf" linenums="10" hl_lines="5-7"
    script:
    """
    echo "Processing ${reads}" > ${meta.id}_report.txt
    echo "Sample: ${meta.id}" >> ${meta.id}_report.txt
    echo "Processed by: \${USER}" >> ${meta.id}_report.txt
    echo "Hostname: \$(hostname)" >> ${meta.id}_report.txt
    echo "Date: \$(date)" >> ${meta.id}_report.txt
    """
    ```

=== "Before - Broken"

    ```groovy title="modules/generate_report.nf" linenums="10"
    script:
    """
    echo "Processing ${reads}" > ${meta.id}_report.txt
    echo "Sample: ${meta.id}" >> ${meta.id}_report.txt
    echo "Processed by: ${USER}" >> ${meta.id}_report.txt
    echo "Hostname: $(hostname)" >> ${meta.id}_report.txt
    echo "Date: $(date)" >> ${meta.id}_report.txt
    """
    ```

Now it works! The backslash (`\`) tells Nextflow "don't interpret this, pass it through to Bash."

### Takeaway

In this section, you've learned:

- **Regular expressions for file parsing**: Using Groovy's `=~` operator and regex patterns to extract metadata from complex bioinformatics file naming conventions
- **Dynamic script generation**: Using Groovy conditional logic to generate different script strings based on input characteristics (like single-end vs paired-end reads)
- **Variable interpolation**: Understanding the difference between Nextflow/Groovy variables (`${var}`), shell environment variables (`\${var}`), and shell command substitution (`\$(cmd)`)

These string processing and generation patterns are essential for handling the diverse file formats and naming conventions you'll encounter in real-world bioinformatics workflows.

---

## 3. Creating Reusable Functions

As your workflow logic becomes more complex, keeping everything inline in channel operators or process definitions can make your code hard to read and maintain. Groovy functions let you extract complex logic into named, reusable components that can be called from anywhere in your workflow.

You may have noticed that the content of our map operation is getting quite long and complex. To keep our workflow maintainable, it's a good idea to break out complex logic into reusable functions.

To illustrate what that looks like with our existing workflow, make the modification below, using `def` to define a reusable function called `separateMetadata`:

=== "After"

    ```groovy title="main.nf" linenums="1" hl_lines="4-24 29"
    include { FASTP } from './modules/fastp.nf'
    include { GENERATE_REPORT } from './modules/generate_report.nf'

    def separateMetadata(row) {
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

    workflow {
        ch_samples = Channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map{ row -> separateMetadata(row) }

        ch_fastp = FASTP(ch_samples)
        GENERATE_REPORT(ch_samples)
    }
    ```

=== "Before"

    ```groovy title="main.nf" linenums="1" hl_lines="7-27"
    include { FASTP } from './modules/fastp.nf'
    include { GENERATE_REPORT } from './modules/generate_report.nf'

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

        ch_fastp = FASTP(ch_samples)
        GENERATE_REPORT(ch_samples)
    }
    ```

By extracting this logic into a function, we've reduced the actual workflow logic down to something much cleaner:

```groovy title="minimal workflow"
    ch_samples = Channel.fromPath("./data/samples.csv")
        .splitCsv(header: true)
        .map{ row -> separateMetadata(row) }

    ch_fastp = FASTP(ch_samples)
    GENERATE_REPORT(ch_samples)
```

This makes the workflow logic much easier to read and understand at a glance. The function `separateMetadata` encapsulates all the complex logic for parsing and enriching metadata, making it reusable and testable.

Run the workflow to make sure it still works:

```bash title="Test reusable function"
nextflow run main.nf
```

```console title="Function results"
 N E X T F L O W   ~  version 25.04.6

Launching `main.nf` [admiring_panini] DSL2 - revision: 8cc832e32f

executor >  local (6)
[8c/2e3f91] process > FASTP (3)           [100%] 3 of 3 ✔
[7a/1b4c92] process > GENERATE_REPORT (3) [100%] 3 of 3 ✔
```

The output should show both processes completing successfully. The workflow is now much cleaner and easier to maintain, with all the complex metadata processing logic encapsulated in the `separateMetadata` function.

### Takeaway

In this section, you've learned:

- **Extracting functions**: Moving complex logic from inline closures into named functions
- **Function scope**: Functions defined at the script level can be called from anywhere in your workflow
- **Cleaner workflows**: Using functions makes your workflow blocks more concise and readable

Next, we'll explore how to use Groovy closures in process directives for dynamic resource allocation.

---

## 4. Dynamic Resource Directives with Closures

So far we've used Groovy in the `script` block of processes. But Groovy closures are also incredibly useful in process directives, especially for dynamic resource allocation. Let's add resource directives to our FASTP process that adapt based on the sample characteristics.

Currently, our FASTP process uses default resources. Let's make it smarter by allocating more CPUs for high-depth samples. Edit `modules/fastp.nf` to include a dynamic `cpus` directive and a static `memory` directive:

=== "After"

    ```groovy title="modules/fastp.nf" linenums="1" hl_lines="4-5"
    process FASTP {
        container 'community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690'

        cpus { meta.depth > 40000000 ? 4 : 2 }
        memory '2 GB'

        input:
        tuple val(meta), path(reads)
    ```

=== "Before"

    ```groovy title="modules/fastp.nf" linenums="1"
    process FASTP {
        container 'community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690'

        input:
        tuple val(meta), path(reads)
    ```

The closure `{ meta.depth > 40000000 ? 4 : 2 }` is evaluated for each task, allowing per-sample resource allocation. High-depth samples (>40M reads) get 4 CPUs, while others get 2 CPUs.

!!! note "Accessing Input Variables in Directives"

    The closure can access any input variables (like `meta` here) because Nextflow evaluates these closures in the context of each task execution.

Run the workflow again:

```bash title="Test resource allocation"
nextflow run main.nf -no-ansi-log
```

We're using the `-no-ansi-log` option to make it easier to see the task hashes.

```console title="Resource allocation output"
N E X T F L O W  ~  version 25.04.6
Launching `main.nf` [fervent_albattani] DSL2 - revision: fa8f249759
[bd/ff3d41] Submitted process > FASTP (2)
[a4/a3aab2] Submitted process > FASTP (1)
[48/6db0c9] Submitted process > FASTP (3)
[ec/83439d] Submitted process > GENERATE_REPORT (3)
[bd/15d7cc] Submitted process > GENERATE_REPORT (2)
[42/699357] Submitted process > GENERATE_REPORT (1)
```

You can check the exact `docker` command that was run to see the CPU allocation for any given task:

```console title="Check docker command"
cat work/48/6db0c9e9d8aa65e4bb4936cd3bd59e/.command.run | grep "docker run"
```

You should see something like:

```bash title="docker command"
    docker run -i --cpu-shares 4096 --memory 2048m -e "NXF_TASK_WORKDIR" -v /workspaces/training/side-quests/groovy_essentials:/workspaces/training/side-quests/groovy_essentials -w "$NXF_TASK_WORKDIR" --name $NXF_BOXID community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690 /bin/bash -ue /workspaces/training/side-quests/groovy_essentials/work/48/6db0c9e9d8aa65e4bb4936cd3bd59e/.command.sh
```

In this example we've chosen an example that requested 4 CPUs (`--cpu-shares 4096`), because it was a high-depth sample, but you should see different CPU allocations depending on the sample depth. Try this for the other tasks as well.

Another powerful pattern is using `task.attempt` for retry strategies. To show why this is useful, we're going to start by reducing the memory allocation to FASTP to less than it needs. Change the `memory` directive in `modules/fastp.nf` to `1.GB`:

=== "After"

    ```groovy title="modules/fastp.nf" linenums="1" hl_lines="4-5"
    process FASTP {
        container 'community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690'

        cpus { meta.depth > 40000000 ? 4 : 2 }
        memory '1 GB'

        input:
        tuple val(meta), path(reads)
    ```

=== "Before"

    ```groovy title="modules/fastp.nf" linenums="1"
    process FASTP {
        container 'community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690'

        cpus { meta.depth > 40000000 ? 4 : 2 }
        memory '2 GB'

        input:
        tuple val(meta), path(reads)
    ```

... and run the workflow again:

```bash title="Test insufficient memory"
nextflow run main.nf
```

You'll see an error indicating that the process was killed for exceeding memory limits:

```console title="Memory error output" hl_lines="2 11"
Command exit status:
  137

Command output:
  (empty)

Command error:
  Detecting adapter sequence for read1...
  No adapter detected for read1

  .command.sh: line 7:   101 Killed                  fastp --in1 SAMPLE_002_S2_L001_R1_001.fastq --out1 sample_002_trimmed.fastq.gz --json sample_002.fastp.json --html sample_002.fastp.html --thread 2
```

This is a very common scenario in real-world workflows - sometimes you just don't know how much memory a task will need until you run it. To make our workflow more robust, we can implement a retry strategy that increases memory allocation on each attempt, once again using a Groovy closure. Modify the `memory` directive to multiply the base memory by `task.attempt`, and add `errorStrategy 'retry'` and `maxRetries 2` directives:

=== "After"

    ```groovy title="modules/fastp.nf" linenums="1" hl_lines="5-7"
    process FASTP {
        container 'community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690'

        cpus { meta.depth > 40000000 ? 4 : 2 }
        memory { '1 GB' * task.attempt }
        errorStrategy 'retry'
        maxRetries 2

        input:
        tuple val(meta), path(reads)
    ```

=== "Before"

    ```groovy title="modules/fastp.nf" linenums="1" hl_lines="5"
    process FASTP {
        container 'community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690'

        cpus { meta.depth > 40000000 ? 4 : 2 }
        memory '2 GB'

        input:
        tuple val(meta), path(reads)
    ```

Now if the process fails due to insufficient memory, Nextflow will retry with more memory:

- First attempt: 1 GB (task.attempt = 1)
- Second attempt: 2 GB (task.attempt = 2)

... and so on, up to the `maxRetries` limit.

### Takeaway

Dynamic directives with Groovy closures let you:

- Allocate resources based on input characteristics
- Implement automatic retry strategies with increasing resources
- Combine multiple factors (metadata, attempt number, priorities)
- Use Groovy logic for complex resource calculations

This makes your workflows both more efficient (not over-allocating) and more robust (automatic retry with more resources).

---

## 5. Conditional Logic and Process Control

Earlier on, we discussed how to use the `.map()` operator to use snippets of Groovy code to transform data flowing through channels. The counterpart to that is using Groovy to not just transform data, but to control which processes get executed based on the data itself. This is essential for building flexible workflows that can adapt to different sample types and analysis requirements.

Nextflow has several [operators](https://www.nextflow.io/docs/latest/reference/operator.html) that control process flow, including, many of which take closures as arguments, meanint their content is evaluated at run time, allowing us to use Groovy logic to drive workflow decisions based on channel content.

### 5.1. Routing with `.branch()`

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

    ```groovy title="main.nf" linenums="28" hl_lines="5-12"
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
    ```

=== "Before"

    ```groovy title="main.nf" linenums="28"
        ch_samples = Channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map(separateMetadata)

        ch_fastp = FASTP(ch_samples)
    ```

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

Here, we've used small but mighty Groovy expressions inside the `.branch{}` operator to route samples based on their metadata. Human samples with high coverage go through `FASTP`, while all other samples go through `TRIMGALORE`.

### 5.2. Using `.filter()` with Groovy Truth

Another powerful pattern for controlling workflow execution is the `.filter()` operator, which uses a closure to determine which items should continue down the pipeline. Let's add a validation step to filter out samples that don't meet our quality requirements.

Add the following before the branch operation:

```groovy title="main.nf - Adding filter" hl_lines="5-9"
    ch_samples = Channel.fromPath("./data/samples.csv")
        .splitCsv(header: true)
            .map(separateMetadata)

    // Filter out invalid or low-quality samples
    ch_valid_samples = ch_samples
        .filter { meta, reads ->
            meta.id && meta.organism && meta.depth >= 1000000
        }

    trim_branches = ch_valid_samples
        .branch { meta, reads ->
            fastp: meta.organism == 'human' && meta.depth >= 30000000
            trimgalore: true
        }
```

This filter uses **Groovy Truth** - Groovy's way of evaluating expressions in boolean contexts:

- `null`, empty strings, empty collections, and zero are all "false"
- Non-null values, non-empty strings, and non-zero numbers are "true"

So `meta.id && meta.organism` checks that both fields exist and are non-empty, while `meta.depth >= 1000000` ensures we have sufficient sequencing depth.

!!! note "Groovy Truth in Practice"

    The expression `meta.id && meta.organism` is more concise than writing:
    ```groovy
    meta.id != null && meta.id != '' && meta.organism != null && meta.organism != ''
    ```

    This makes filtering logic much cleaner and easier to read.

You could also combine `.filter()` with more complex Groovy logic:

```groovy title="Complex filtering examples"
// Filter using safe navigation and Elvis operators
ch_samples
    .filter { meta, reads ->
        (meta.quality ?: 0) > 30 && meta.organism?.toLowerCase() in ['human', 'mouse']
    }

// Filter using regular expressions
ch_samples
    .filter { meta, reads ->
        meta.id =~ /^SAMPLE_\d+$/ && reads.exists()
    }

// Filter using multiple conditions with Groovy Truth
ch_samples
    .filter { meta, reads ->
        meta.files          // Non-empty file list
        && meta.paired      // Boolean flag is true
        && !meta.failed     // Negative check
    }
```

### Takeaway

In this section, you've learned to use Groovy logic to control workflow execution using the closure interfaces of Nextflow operators like `.branch{}` and `.filter{}`, leveraging Groovy Truth to write concise conditional expressions.

Our pipeline now intelligently routes samples through appropriate processes, but production workflows need to handle invalid data gracefully. Let's make our workflow robust against missing or null values.

---

## 6. Safe Navigation and Elvis Operators

Our `separateMetadata` function currently assumes all CSV fields are present and valid. But what happens with incomplete data? Let's find out.

### 6.1. The Problem: Null Pointer Crashes

Add a row with missing data to your `data/samples.csv`:

```csv
SAMPLE_004,,unknown_tissue,20000000,data/sequences/SAMPLE_004_S4_L001_R1_001.fastq,
```

Notice the empty organism field and missing quality_score. Now try running the workflow:

```bash
nextflow run main.nf
```

It crashes with a NullPointerException! This is where Groovy's safe operators save the day.

### 6.2. Safe Navigation Operator (`?.`)

The safe navigation operator (`?.`) returns null instead of throwing an exception. Update your `separateMetadata` function:

=== "After"

    ```groovy title="main.nf" linenums="4" hl_lines="6-8"
    def separateMetadata(row) {
        def sample_meta = [
            id: row.sample_id.toLowerCase(),
            organism: row.organism?.toLowerCase(),
            tissue: row.tissue_type?.replaceAll('_', ' ')?.toLowerCase(),
            depth: row.sequencing_depth.toInteger(),
            quality: row.quality_score?.toDouble()
        ]
        // ... rest unchanged
    ```

=== "Before"

    ```groovy title="main.nf" linenums="4"
    def separateMetadata(row) {
        def sample_meta = [
            id: row.sample_id.toLowerCase(),
            organism: row.organism,
            tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
            depth: row.sequencing_depth.toInteger(),
            quality: row.quality_score.toDouble()
        ]
        // ... rest unchanged
    ```

Run again:

```bash
nextflow run main.nf
```

No crash! But SAMPLE_004 now has `null` values which could cause problems downstream.

### 6.3. Elvis Operator (`?:`) for Defaults

The Elvis operator (`?:`) provides default values. Update again:

=== "After"

    ```groovy title="main.nf" linenums="4" hl_lines="6-8"
    def separateMetadata(row) {
        def sample_meta = [
            id: row.sample_id.toLowerCase(),
            organism: row.organism?.toLowerCase() ?: 'unknown',
            tissue: row.tissue_type?.replaceAll('_', ' ')?.toLowerCase() ?: 'unknown',
            depth: row.sequencing_depth.toInteger(),
            quality: row.quality_score?.toDouble() ?: 0.0
        ]
        // ... rest unchanged
    ```

=== "Before"

    ```groovy title="main.nf" linenums="4"
    def separateMetadata(row) {
        def sample_meta = [
            id: row.sample_id.toLowerCase(),
            organism: row.organism?.toLowerCase(),
            tissue: row.tissue_type?.replaceAll('_', ' ')?.toLowerCase(),
            depth: row.sequencing_depth.toInteger(),
            quality: row.quality_score?.toDouble()
        ]
        // ... rest unchanged
    ```

Run once more:

```bash
nextflow run main.nf
```

Perfect! SAMPLE_004 now has safe defaults: 'unknown' for organism/tissue, 0.0 for quality.

### 6.4. Filtering with Safe Operators

Now let's filter out samples with missing data. Update your workflow:

=== "After"

    ```groovy title="main.nf" linenums="28" hl_lines="4-7"
    workflow {
        ch_samples = Channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map(separateMetadata)
            .filter { meta, reads ->
                meta.organism != 'unknown' && (meta.quality ?: 0) > 0
            }

        // ... rest of workflow
    ```

=== "Before"

    ```groovy title="main.nf" linenums="28"
    workflow {
        ch_samples = Channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map(separateMetadata)

        // ... rest of workflow
    ```

Run the workflow:

```bash
nextflow run main.nf
```

SAMPLE_004 is now filtered out! Only valid samples proceed.

### Takeaway

- **Safe navigation (`?.`)**: Prevents crashes on null values - returns null instead of throwing exception
- **Elvis operator (`?:`)**: Provides defaults - `value ?: 'default'`
- **Combining**: `value?.method() ?: 'default'` is the common pattern
- **In filters**: Use to handle missing data: `(meta.quality ?: 0) > threshold`

These operators make workflows resilient to incomplete data - essential for real-world bioinformatics.

---

## 7. Validation with `error()` and `log.warn`

Sometimes you need to stop the workflow immediately if input parameters are invalid. Nextflow provides `error()` for this. Let's add validation to our workflow.

Create a validation function before your workflow block:

=== "After"

    ```groovy title="main.nf" linenums="1" hl_lines="3-18"
    include { FASTP } from './modules/fastp.nf'
    include { TRIMGALORE } from './modules/trimgalore.nf'
    include { GENERATE_REPORT } from './modules/generate_report.nf'

    def validateInputs() {
        // Check CSV file exists
        if (!file(params.input ?: './data/samples.csv').exists()) {
            error("Input CSV file not found: ${params.input ?: './data/samples.csv'}")
        }

        // Warn if output directory already exists
        if (file(params.outdir ?: 'results').exists()) {
            log.warn "Output directory already exists: ${params.outdir ?: 'results'}"
        }

        // Check for required genome parameter
        if (params.run_gatk && !params.genome) {
            error("Genome reference required when running GATK. Please provide --genome")
        }
    }

    // ... separateMetadata function ...

    workflow {
        validateInputs()  // Call validation first

        ch_samples = Channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map(separateMetadata)
        // ... rest of workflow
    }
    ```

=== "Before"

    ```groovy title="main.nf" linenums="1"
    include { FASTP } from './modules/fastp.nf'
    include { TRIMGALORE } from './modules/trimgalore.nf'
    include { GENERATE_REPORT } from './modules/generate_report.nf'

    // ... separateMetadata function ...

    workflow {
        ch_samples = Channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map(separateMetadata)
        // ... rest of workflow
    }
    ```

Now try running without the CSV file:

```bash
mv data/samples.csv data/samples.csv.bak
nextflow run main.nf
```

The workflow stops immediately with a clear error message instead of failing mysteriously later!

You can also add validation within the `separateMetadata` function:

```groovy title="main.nf - Validation in function"
def separateMetadata(row) {
    // Validate required fields
    if (!row.sample_id) {
        error("Missing sample_id in CSV row: ${row}")
    }

    def sample_meta = [
        id: row.sample_id.toLowerCase(),
        organism: row.organism?.toLowerCase() ?: 'unknown',
        // ... rest of fields
    ]

    // Validate data makes sense
    if (sample_meta.depth < 1000000) {
        log.warn "Low sequencing depth for ${sample_meta.id}: ${sample_meta.depth}"
    }

    return [sample_meta, file(row.file_path)]
}
```

### Takeaway

- **`error()`**: Stops workflow immediately with clear message
- **`log.warn`**: Issues warnings without stopping workflow
- **Early validation**: Check inputs before processing to fail fast with helpful errors
- **Validation functions**: Create reusable validation logic that can be called at workflow start

Proper validation makes workflows more robust and user-friendly by catching problems early with clear error messages.

---

## 8. Groovy in Configuration: Workflow Event Handlers

Up until now, we've been writing Groovy code in our workflow scripts and process definitions. But there's one more important place where Groovy is essential: workflow event handlers in your `nextflow.config` file.

Event handlers are Groovy closures that run at specific points in your workflow's lifecycle. They're perfect for adding logging, notifications, or cleanup operations without cluttering your main workflow code.

### 8.1. The `onComplete` Handler

The most commonly used event handler is `onComplete`, which runs when your workflow finishes (whether it succeeded or failed). Let's add one to summarize our pipeline results.

Your `nextflow.config` file already has Docker enabled. Add an event handler after the existing configuration:

=== "After"

    ```groovy title="nextflow.config" linenums="1" hl_lines="5-15"
    // Nextflow configuration for Groovy Essentials side quest

    docker.enabled = true

    workflow.onComplete = {
        println ""
        println "Pipeline execution summary:"
        println "=========================="
        println "Completed at: ${workflow.complete}"
        println "Duration    : ${workflow.duration}"
        println "Success     : ${workflow.success}"
        println "workDir     : ${workflow.workDir}"
        println "exit status : ${workflow.exitStatus}"
        println ""
    }
    ```

=== "Before"

    ```groovy title="nextflow.config" linenums="1"
    // Nextflow configuration for Groovy Essentials side quest

    docker.enabled = true
    ```

This is a Groovy closure being assigned to `workflow.onComplete`. Inside, you have access to the `workflow` object which provides useful properties about the execution.

Run your workflow and you'll see this summary appear at the end!

Let's make it more useful by adding conditional logic:

=== "After"

    ```groovy title="nextflow.config" linenums="5" hl_lines="11-18"
    workflow.onComplete = {
        println ""
        println "Pipeline execution summary:"
        println "=========================="
        println "Completed at: ${workflow.complete}"
        println "Duration    : ${workflow.duration}"
        println "Success     : ${workflow.success}"
        println "workDir     : ${workflow.workDir}"
        println "exit status : ${workflow.exitStatus}"
        println ""

        if (workflow.success) {
            println "✅ Pipeline completed successfully!"
            println "Results are in: ${params.outdir ?: 'results'}"
        } else {
            println "❌ Pipeline failed!"
            println "Error: ${workflow.errorMessage}"
        }
    }
    ```

=== "Before"

    ```groovy title="nextflow.config" linenums="5"
    workflow.onComplete = {
        println ""
        println "Pipeline execution summary:"
        println "=========================="
        println "Completed at: ${workflow.complete}"
        println "Duration    : ${workflow.duration}"
        println "Success     : ${workflow.success}"
        println "workDir     : ${workflow.workDir}"
        println "exit status : ${workflow.exitStatus}"
        println ""
    }
    ```

You can also write the summary to a file using Groovy file operations:

```groovy title="nextflow.config - Writing summary to file"
workflow.onComplete = {
    def summary = """
    Pipeline Execution Summary
    ===========================
    Completed: ${workflow.complete}
    Duration : ${workflow.duration}
    Success  : ${workflow.success}
    Command  : ${workflow.commandLine}
    """

    println summary

    // Write to a log file
    def log_file = file("${workflow.launchDir}/pipeline_summary.txt")
    log_file.text = summary
}
```

### 8.2. Other Useful Event Handlers

Besides `onComplete`, there are other event handlers you can use:

**`onStart`** - Runs when the workflow begins:

```groovy title="nextflow.config - onStart handler"
workflow.onStart = {
    println "="* 50
    println "Starting pipeline: ${workflow.runName}"
    println "Project directory: ${workflow.projectDir}"
    println "Launch directory: ${workflow.launchDir}"
    println "Work directory: ${workflow.workDir}"
    println "="* 50
}
```

**`onError`** - Runs only if the workflow fails:

```groovy title="nextflow.config - onError handler"
workflow.onError = {
    println "="* 50
    println "Pipeline execution failed!"
    println "Error message: ${workflow.errorMessage}"
    println "="* 50

    // Write detailed error log
    def error_file = file("${workflow.launchDir}/error.log")
    error_file.text = """
    Workflow Error Report
    =====================
    Time: ${new Date()}
    Error: ${workflow.errorMessage}
    Error report: ${workflow.errorReport ?: 'No detailed report available'}
    """

    println "Error details written to: ${error_file}"
}
```

You can use multiple handlers together:

```groovy title="nextflow.config - Combined handlers"
workflow.onStart = {
    println "Starting ${workflow.runName} at ${workflow.start}"
}

workflow.onError = {
    println "Workflow failed: ${workflow.errorMessage}"
}

workflow.onComplete = {
    def duration_mins = workflow.duration.toMinutes().round(2)
    def status = workflow.success ? "SUCCESS ✅" : "FAILED ❌"

    println """
    Pipeline finished: ${status}
    Duration: ${duration_mins} minutes
    """
}
```

### Takeaway

In this section, you've learned:

- **Event handler closures**: Groovy closures in `nextflow.config` that run at different lifecycle points
- **`onComplete` handler**: For execution summaries and result reporting
- **`onStart` handler**: For logging pipeline initialization
- **`onError` handler**: For error handling and logging failures
- **Workflow object properties**: Accessing `workflow.success`, `workflow.duration`, `workflow.errorMessage`, etc.

Event handlers are pure Groovy code running in your config file, demonstrating that Nextflow configuration is actually a Groovy script with access to the full language.

---

## Summary

Throughout this side quest, you've built a comprehensive sample processing pipeline that evolved from basic metadata handling to a sophisticated, production-ready workflow. Each section built upon the previous, demonstrating how Groovy transforms simple Nextflow workflows into powerful data processing systems.

Here's how we progressively enhanced our pipeline:

1. **Nextflow vs Groovy Boundaries**: You learned to distinguish between workflow orchestration (Nextflow) and programming logic (Groovy), including the crucial differences between constructs like `collect`.

2. **Advanced String Processing**: You mastered regular expressions for parsing file names, dynamic script generation in processes, and variable interpolation (Groovy vs Bash vs Shell).

3. **Creating Reusable Functions**: You learned to extract complex logic into named functions that can be called from channel operators, making workflows more readable and maintainable.

4. **Dynamic Resource Directives with Closures**: You explored using Groovy closures in process directives for adaptive resource allocation based on input characteristics.

5. **Conditional Logic and Process Control**: You added intelligent routing using `.branch()` and `.filter()` operators, leveraging Groovy Truth for concise conditional expressions.

6. **Safe Navigation and Elvis Operators**: You made the pipeline robust against missing data using `?.` for null-safe property access and `?:` for providing default values.

7. **Validation with error() and log.warn**: You learned to validate inputs early and fail fast with clear error messages.

8. **Groovy in Configuration**: You learned to use workflow event handlers (`onComplete`, `onStart`, `onError`) for logging, notifications, and lifecycle management.

### Key Benefits

- **Clearer code**: Understanding when to use Nextflow and Groovy helps you write more organized workflows
- **Robust handling**: Safe navigation and Elvis operators make workflows resilient to missing data
- **Flexible processing**: Conditional logic lets your workflows process different sample types appropriately
- **Adaptive resources**: Dynamic directives optimize resource usage based on input characteristics

### From Simple to Sophisticated

The pipeline journey you completed demonstrates the evolution from basic data processing to production-ready bioinformatics workflows:

1. **Started simple**: Basic CSV processing and metadata extraction with clear Nextflow vs Groovy boundaries
2. **Added intelligence**: Dynamic file name parsing with regex patterns, variable interpolation mastery, and dynamic script generation based on input types
3. **Made it maintainable**: Extracted complex logic into reusable functions for cleaner, more testable code
4. **Made it efficient**: Dynamic resource allocation with closures in directives and retry strategies
5. **Added routing**: Conditional logic to route samples through appropriate processes based on their characteristics
6. **Made it robust**: Safe navigation and Elvis operators for handling missing data gracefully, plus validation for early error detection
7. **Added observability**: Workflow event handlers for logging, notifications, and lifecycle management

This progression mirrors the real-world evolution of bioinformatics pipelines - from research prototypes handling a few samples to production systems processing thousands of samples across laboratories and institutions. Every challenge you solved and pattern you learned reflects actual problems developers face when scaling Nextflow workflows.

### Next Steps

With these Groovy fundamentals mastered, you're ready to:

- Write cleaner workflows with proper separation between Nextflow and Groovy logic
- Master variable interpolation to avoid common pitfalls with Groovy, Bash, and shell variables
- Use dynamic resource directives for efficient, adaptive workflows
- Transform file collections into properly formatted command-line arguments
- Handle different file naming conventions and input formats gracefully using regex and string processing
- Build reusable, maintainable code using advanced closure patterns and functional programming
- Process and organize complex datasets using collection operations
- Add validation, error handling, and logging to make your workflows production-ready
- Implement workflow lifecycle management with event handlers

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
