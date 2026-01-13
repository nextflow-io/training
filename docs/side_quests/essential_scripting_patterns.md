# Essential Nextflow Scripting Patterns

Nextflow is a programming language that runs on the Java Virtual Machine. While Nextflow is built on [Groovy](http://groovy-lang.org/) and shares much of its syntax, Nextflow is more than just "Groovy with extensions" -- it is a standalone language with a fully-specified [syntax](https://nextflow.io/docs/latest/reference/syntax.html) and [standard library](https://nextflow.io/docs/latest/reference/stdlib.html).

You can write a lot of Nextflow without venturing beyond basic syntax for variables, maps, and lists. Most Nextflow tutorials focus on workflow orchestration (channels, processes, and data flow), and you can go surprisingly far with just that.

However, when you need to manipulate data, parse complex filenames, implement conditional logic, or build robust production workflows, it helps to think about two distinct aspects of your code: **dataflow** (channels, operators, processes, and workflows) and **scripting** (the code inside closures, functions, and process scripts). While this distinction is somewhat arbitrary—it's all Nextflow code—it provides a useful mental model for understanding when you're orchestrating your pipeline versus when you're manipulating data. Mastering both dramatically improves your ability to write clear, maintainable workflows.

### Learning goals

This side quest takes you on a hands-on journey from basic concepts to production-ready patterns.
We'll transform a simple CSV-reading workflow into a sophisticated bioinformatics pipeline, evolving it step-by-step through realistic challenges:

- **Understanding boundaries:** Distinguish between dataflow operations and scripting, and understand how they work together
- **Data manipulation:** Extract, transform, and subset maps and collections using powerful operators
- **String processing:** Parse complex file naming schemes with regex patterns and master variable interpolation
- **Reusable functions:** Extract complex logic into named functions for cleaner, more maintainable workflows
- **Dynamic logic:** Build processes that adapt to different input types and use closures for dynamic resource allocation
- **Conditional routing:** Intelligently route samples through different processes based on their metadata characteristics
- **Safe operations:** Handle missing data gracefully with null-safe operators and validate inputs with clear error messages
- **Configuration-based handlers:** Use workflow event handlers for logging, notifications, and lifecycle management

### Prerequisites

Before taking on this side quest, you should:

- Have completed the [Hello Nextflow](../hello_nextflow/README.md) tutorial or equivalent beginner's course.
- Be comfortable using basic Nextflow concepts and mechanisms (processes, channels, operators, working with files, meta data)
- Have basic familiarity with common programming constructs (variables, maps, lists)

This tutorial will explain programming concepts as we encounter them, so you don't need extensive programming experience.
We'll start with fundamental concepts and build up to advanced patterns.

---

## 0. Get started

#### Open the training codespace

If you haven't yet done so, make sure to open the training environment as described in the [Environment Setup](../envsetup/index.md).

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### Move into the project directory

Let's move into the directory where the files for this tutorial are located.

```bash
cd side-quests/essential_scripting_patterns
```

#### Review the materials

You'll find a main workflow file and a `data` directory containing example data files.

```console title="Directory contents"
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
```

Our sample CSV contains information about biological samples that need different processing based on their characteristics:

```console title="samples.csv"
sample_id,organism,tissue_type,sequencing_depth,file_path,quality_score
SAMPLE_001,human,liver,30000000,data/sequences/SAMPLE_001_S1_L001_R1_001.fastq,38.5
SAMPLE_002,mouse,brain,25000000,data/sequences/SAMPLE_002_S2_L001_R1_001.fastq,35.2
SAMPLE_003,human,kidney,45000000,data/sequences/SAMPLE_003_S3_L001_R1_001.fastq,42.1
```

We'll use this realistic dataset to explore practical programming techniques that you'll encounter in real bioinformatics workflows.

<!-- TODO: Can we make this more domain-agnostic? -->

<!-- TODO: add an assignment statement? #### Review the assignment -->

#### Readiness checklist

Think you're ready to dive in?

- [ ] I understand the goal of this course and its prerequisites
- [ ] My codespace is up and running
- [ ] I've set my working directory appropriately
<!-- - [ ] I understand the assignment -->

If you can check all the boxes, you're good to go.

---

## 1. Dataflow vs Scripting: Understanding the Boundaries

### 1.1. Identifying What's What

When writing Nextflow workflows, it's important to distinguish between **dataflow** (how data moves through channels and processes) and **scripting** (the code that manipulates data and makes decisions). Let's build a workflow demonstrating how they work together.

#### 1.1.1. Basic Nextflow Workflow

Start with a simple workflow that just reads the CSV file (we've already done this for you in `main.nf`):

```groovy title="main.nf" linenums="1"
workflow {
    ch_samples = channel.fromPath("./data/samples.csv")
        .splitCsv(header: true)
        .view()
}
```

The `workflow` block defines our pipeline structure, while `channel.fromPath()` creates a channel from a file path. The `.splitCsv()` operator processes the CSV file and converts each row into a map data structure.

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

#### 1.1.2. Adding the Map Operator

Now we're going to add scripting to transform the data, using the `.map()` operator you will probably already be familiar with. This operator takes a 'closure' where we can write code to transform each item.

!!! note

    A **closure** is a block of code that can be passed around and executed later. Think of it as a function that you define inline. Closures are written with curly braces `{ }` and can take parameters. They're fundamental to how Nextflow operators work and if you've been writing Nextflow for a while, you may already have been using them without realizing it!

Here's what that map operation looks like:

=== "After"

    ```groovy title="main.nf" linenums="2" hl_lines="3-6"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row ->
                return row
            }
            .view()
    ```

=== "Before"

    ```groovy title="main.nf" linenums="2" hl_lines="3"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .view()
    ```

This is our first **closure** - an anonymous function you can pass as an argument (similar to lambdas in Python or arrow functions in JavaScript). Closures are essential for working with Nextflow operators.

The closure `{ row -> return row }` takes a parameter `row` (could be any name: `item`, `sample`, etc.).

When the `.map()` operator processes each channel item, it passes that item to your closure. Here, `row` holds one CSV row at a time.

Apply this change and run the workflow:

```bash title="Test map operator"
nextflow run main.nf
```

You'll see the same output as before, because we're simply returning the input unchanged. This confirms that the map operator is working correctly. Now let's start transforming the data.

#### 1.1.3. Creating a Map Data Structure

Now we're going to write **scripting** logic inside our closure to transform each row of data. This is where we process individual data items rather than orchestrating data flow.

=== "After"

    ```groovy title="main.nf" linenums="2" hl_lines="4-12"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row ->
                // Scripting for data transformation
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
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row ->
                return row
            }
            .view()
    ```

The `sample_meta` map is a key-value data structure (like dictionaries in Python, objects in JavaScript, or hashes in Ruby) storing related information: sample ID, organism, tissue type, sequencing depth, and quality score.

We use string manipulation methods like `.toLowerCase()` and `.replaceAll()` to clean up our data, and type conversion methods like `.toInteger()` and `.toDouble()` to convert string data from the CSV into the appropriate numeric types.

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

#### 1.1.4. Adding Conditional Logic

Now let's add more scripting - this time using a ternary operator to make decisions based on data values.

Make the following change:

=== "After"

    ```groovy title="main.nf" linenums="2" hl_lines="11-12"
        ch_samples = channel.fromPath("./data/samples.csv")
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
        ch_samples = channel.fromPath("./data/samples.csv")
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

#### 1.1.5. Subsetting Maps with `.subMap()`

While the `+` operator adds keys to a map, sometimes you need to do the opposite - extract only specific keys. The `.subMap()` method is perfect for this.

Let's add a line to create a simplified version of our metadata that only contains identification fields:

=== "After"

    ```groovy title="main.nf" linenums="2" hl_lines="12-15"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row ->
                // Scripting for data transformation
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
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row ->
                // Scripting for data transformation
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
 N E X T F L O W   ~  version 25.04.3

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

#### 1.1.6. Combining Maps and Returning Results

So far, we've only been returning what the Nextflow community calls the 'meta map', and we've been ignoring the files those metadata relate to. But if you're writing Nextflow workflows, you probably want to do something with those files.

Let's output a channel structure comprising a tuple of 2 elements: the enriched metadata map and the corresponding file path. This is a common pattern in Nextflow for passing data to processes.

=== "After"

    ```groovy title="main.nf" linenums="2" hl_lines="12"
        ch_samples = channel.fromPath("./data/samples.csv")
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
                return tuple( sample_meta + [priority: priority], file(row.file_path) )
            }
            .view()
    ```

=== "Before"

    ```groovy title="main.nf" linenums="2" hl_lines="12"
        ch_samples = channel.fromPath("./data/samples.csv")
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
[[id:sample_001, organism:human, tissue:liver, depth:30000000, quality:38.5, priority:normal], /workspaces/training/side-quests/essential_scripting_patterns/data/sequences/SAMPLE_001_S1_L001_R1_001.fastq]
[[id:sample_002, organism:mouse, tissue:brain, depth:25000000, quality:35.2, priority:normal], /workspaces/training/side-quests/essential_scripting_patterns/data/sequences/SAMPLE_002_S2_L001_R1_001.fastq]
[[id:sample_003, organism:human, tissue:kidney, depth:45000000, quality:42.1, priority:high], /workspaces/training/side-quests/essential_scripting_patterns/data/sequences/SAMPLE_003_S3_L001_R1_001.fastq]
```

This `[meta, file]` tuple structure is a common pattern in Nextflow for passing both metadata and associated files to processes.

!!! note

    **Maps and Metadata**: Maps are fundamental to working with metadata in Nextflow. For a more detailed explanation of working with metadata maps, see the [Working with metadata](./metadata.md) side quest.

Our workflow demonstrates the core pattern: **dataflow operations** (`workflow`, `channel.fromPath()`, `.splitCsv()`, `.map()`, `.view()`) orchestrate how data moves through the pipeline, while **scripting** (maps `[key: value]`, string methods, type conversions, ternary operators) inside the `.map()` closure handles the transformation of individual data items.

### 1.2. Understanding Different Types: Channel vs List

So far, so good, we can distinguish between dataflow operations and scripting. But what about when the same method name exists in both contexts?

A perfect example is the `collect` method, which exists for both channel types and List types in the Nextflow standard library. The `collect()` method on a List transforms each element, while the `collect()` operator on a channel gathers all channel emissions into a single-item channel.

Let's demonstrate this with some sample data, starting by refreshing ourselves on what the channel `collect()` operator does. Check out `collect.nf`:

```groovy title="collect.nf" linenums="1"
def sample_ids = ['sample_001', 'sample_002', 'sample_003']

// channel.collect() - groups multiple channel emissions into one
ch_input = channel.fromList(sample_ids)
ch_input.view { sample -> "Individual channel item: ${sample}" }
ch_collected = ch_input.collect()
ch_collected.view { list -> "channel.collect() result: ${list} (${list.size()} items grouped into 1)" }
```

Steps:

- Define a List of sample IDs
- Create a channel with `fromList()` that emits each sample ID separately
- Print each item with `view()` as it flows through
- Gather all items into a single list with the channel's `collect()` operator
- Print the collected result (single item containing all sample IDs) with a second `view()`

We've changed the structure of the channel, but we haven't changed the data itself.

Run the workflow to confirm this:

```bash title="Test collect operations"
nextflow run collect.nf
```

```console title="Different collect behaviors"
 N E X T F L O W   ~  version 25.04.3

Launching `collect.nf` [loving_mendel] DSL2 - revision: e8d054a46e

Individual channel item: sample_001
Individual channel item: sample_002
Individual channel item: sample_003
channel.collect() result: [sample_001, sample_002, sample_003] (3 items grouped into 1)
```

`view()` returns an output for every channel emission, so we know that this single output contains all 3 original items grouped into one list.

Now let's see the `collect` method on a List in action. Modify `collect.nf` to apply the List's `collect` method to the original list of sample IDs:

=== "After"

    ```groovy title="main.nf" linenums="1" hl_lines="9-13"
    def sample_ids = ['sample_001', 'sample_002', 'sample_003']

    // channel.collect() - groups multiple channel emissions into one
    ch_input = channel.fromList(sample_ids)
    ch_input.view { sample -> "Individual channel item: ${sample}" }
    ch_collected = ch_input.collect()
    ch_collected.view { list -> "channel.collect() result: ${list} (${list.size()} items grouped into 1)" }

    // List.collect() - transforms each element, preserves structure
    def formatted_ids = sample_ids.collect { id ->
        id.toUpperCase().replace('SAMPLE_', 'SPECIMEN_')
    }
    println "List.collect() result: ${formatted_ids} (${sample_ids.size()} items transformed into ${formatted_ids.size()})"
    ```

=== "Before"

    ```groovy title="main.nf" linenums="1"
    def sample_ids = ['sample_001', 'sample_002', 'sample_003']

    // channel.collect() - groups multiple channel emissions into one
    ch_input = channel.fromList(sample_ids)
    ch_input.view { sample -> "Individual channel item: ${sample}" }
    ch_collected = ch_input.collect()
    ch_collected.view { list -> "channel.collect() result: ${list} (${list.size()} items grouped into 1)" }
    ```

In this new snippet we:

- Define a new variable `formatted_ids` that uses the List's `collect` method to transform each sample ID in the original list
- Print the result using `println`

Run the modified workflow:

```bash title="Test List collect"
nextflow run collect.nf
```

```console title="List collect results" hl_lines="5"
 N E X T F L O W   ~  version 25.04.3

Launching `collect.nf` [cheeky_stonebraker] DSL2 - revision: 2d5039fb47

List.collect() result: [SPECIMEN_001, SPECIMEN_002, SPECIMEN_003] (3 items transformed into 3)
Individual channel item: sample_001
Individual channel item: sample_002
Individual channel item: sample_003
channel.collect() result: [sample_001, sample_002, sample_003] (3 items grouped into 1)
```

This time, we have NOT changed the structure of the data, we still have 3 items in the list, but we HAVE transformed each item using the List's `collect` method to produce a new list with modified values. This is similar to using the `map` operator on a channel, but it's operating on a List data structure rather than a channel.

`collect` is an extreme case we're using here to make a point. The key lesson is that when you're writing workflows, always distinguish between **data structures** (Lists, Maps, etc.) and **channels** (dataflow constructs). Operations can share names but behave completely differently depending on the type they're called on.

### 1.3. The Spread Operator (`*.`) - Shorthand for Property Extraction

Related to the List's `collect` method is the spread operator (`*.`), which provides a concise way to extract properties from collections. It's essentially syntactic sugar for a common `collect` pattern.

Let's add a demonstration to our `collect.nf` file:

=== "After"

    ```groovy title="collect.nf" linenums="1" hl_lines="15-18"
    def sample_ids = ['sample_001', 'sample_002', 'sample_003']

    // channel.collect() - groups multiple channel emissions into one
    ch_input = channel.fromList(sample_ids)
    ch_input.view { sample -> "Individual channel item: ${sample}" }
    ch_collected = ch_input.collect()
    ch_collected.view { list -> "channel.collect() result: ${list} (${list.size()} items grouped into 1)" }

    // List.collect() - transforms each element, preserves structure
    def formatted_ids = sample_ids.collect { id ->
        id.toUpperCase().replace('SAMPLE_', 'SPECIMEN_')
    }
    println "List.collect() result: ${formatted_ids} (${sample_ids.size()} items transformed into ${formatted_ids.size()})"

    // Spread operator - concise property access
    def sample_data = [[id: 's1', quality: 38.5], [id: 's2', quality: 42.1], [id: 's3', quality: 35.2]]
    def all_ids = sample_data*.id
    println "Spread operator result: ${all_ids}"
    ```

=== "Before"

    ```groovy title="collect.nf" linenums="1"
    def sample_ids = ['sample_001', 'sample_002', 'sample_003']

    // channel.collect() - groups multiple channel emissions into one
    ch_input = channel.fromList(sample_ids)
    ch_input.view { sample -> "Individual channel item: ${sample}" }
    ch_collected = ch_input.collect()
    ch_collected.view { list -> "channel.collect() result: ${list} (${list.size()} items grouped into 1)" }

    // List.collect() - transforms each element, preserves structure
    def formatted_ids = sample_ids.collect { id ->
        id.toUpperCase().replace('SAMPLE_', 'SPECIMEN_')
    }
    println "List.collect() result: ${formatted_ids} (${sample_ids.size()} items transformed into ${formatted_ids.size()})"
    ```

Run the updated workflow:

```bash title="Test spread operator"
nextflow run collect.nf
```

You should see output like:

```console title="Spread operator output" hl_lines="6"
 N E X T F L O W   ~  version 25.04.3

Launching `collect.nf` [cranky_galileo] DSL2 - revision: 5f3c8b2a91

List.collect() result: [SPECIMEN_001, SPECIMEN_002, SPECIMEN_003] (3 items transformed into 3)
Spread operator result: [s1, s2, s3]
Individual channel item: sample_001
Individual channel item: sample_002
Individual channel item: sample_003
channel.collect() result: [sample_001, sample_002, sample_003] (3 items grouped into 1)
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

!!! tip "When to Use Spread vs Collect"

    - **Use spread (`*.`)** for simple property access: `samples*.id`, `files*.name`
    - **Use collect** for transformations or complex logic: `samples.collect { it.id.toUpperCase() }`, `samples.collect { [it.id, it.quality > 40] }`

### Takeaway

In this section, you've learned:

- **Dataflow vs scripting**: Channel operators orchestrate how data flows through your pipeline, while scripting transforms individual data items
- **Understanding types**: The same method name (like `collect`) can behave differently depending on the type it's called on (Channel vs List)
- **Context matters**: Always be aware of whether you're working with channels (dataflow) or data structures (scripting)

Understanding these boundaries is essential for debugging, documentation, and writing maintainable workflows.

Next we'll dive deeper into string processing capabilities, which are essential for handling real-world data.

---

## 2. String Processing and Dynamic Script Generation

Mastering string processing separates brittle workflows from robust pipelines. This section covers parsing complex file names, dynamic script generation, and variable interpolation.

### 2.1. Pattern Matching and Regular Expressions

Bioinformatics files often have complex naming conventions encoding metadata. Let's extract this automatically using pattern matching with regular expressions.

We're going to return to our `main.nf` workflow and add some pattern matching logic to extract additional sample information from file names. The FASTQ files in our dataset follow Illumina-style naming conventions with names like `SAMPLE_001_S1_L001_R1_001.fastq.gz`. These might look cryptic, but they actually encode useful metadata like sample ID, lane number, and read direction. We're going to use regex capabilities to parse these names.

Make the following change to your existing `main.nf` workflow:

=== "After"

    ```groovy title="main.nf" linenums="4" hl_lines="10-21"
            .map { row ->
                // Scripting for data transformation
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
                return tuple(sample_meta + file_meta + [priority: priority], fastq_path)
            }
    ```

=== "Before"

    ```groovy title="main.nf" linenums="4" hl_lines="10-11"
            .map { row ->
                // Scripting for data transformation
                def sample_meta = [
                    id: row.sample_id.toLowerCase(),
                    organism: row.organism,
                    tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                    depth: row.sequencing_depth.toInteger(),
                    quality: row.quality_score.toDouble()
                ]
                def priority = sample_meta.quality > 40 ? 'high' : 'normal'
                return tuple(sample_meta + [priority: priority], file(row.file_path))
            }
    ```

This demonstrates key **string processing concepts**:

1. **Regular expression literals** using `~/pattern/` syntax - this creates a regex pattern without needing to escape backslashes
2. **Pattern matching** with the `=~` operator - this attempts to match a string against a regex pattern
3. **Matcher objects** that capture groups with `[0][1]`, `[0][2]`, etc. - `[0]` refers to the entire match, `[1]`, `[2]`, etc. refer to captured groups in parentheses

Let's break down the regex pattern `^(.+)_S(\d+)_L(\d{3})_(R[12])_(\d{3})\.fastq(?:\.gz)?$`:

| Pattern             | Matches                                | Captures                           |
| ------------------- | -------------------------------------- | ---------------------------------- |
| `^(.+)`             | Sample name from start                 | Group 1: sample name               |
| `_S(\d+)`           | Sample number `_S1`, `_S2`, etc.       | Group 2: sample number             |
| `_L(\d{3})`         | Lane number `_L001`                    | Group 3: lane (3 digits)           |
| `_(R[12])`          | Read direction `_R1` or `_R2`          | Group 4: read direction            |
| `_(\d{3})`          | Chunk number `_001`                    | Group 5: chunk (3 digits)          |
| `\.fastq(?:\.gz)?$` | File extension `.fastq` or `.fastq.gz` | Not captured (?: is non-capturing) |

This parses Illumina-style naming conventions to extract metadata automatically.

Run the modified workflow:

```bash title="Test pattern matching"
nextflow run main.nf
```

You should see output with metadata enriched from the file names, like

```console title="Metadata with file parsing"
 N E X T F L O W   ~  version 25.04.3

Launching `main.nf` [clever_pauling] DSL2 - revision: 605d2058b4

[[id:sample_001, organism:human, tissue:liver, depth:30000000, quality:38.5, sample_num:1, lane:001, read:R1, chunk:001, priority:normal], /workspaces/training/side-quests/essential_scripting_patterns/data/sequences/SAMPLE_001_S1_L001_R1_001.fastq]
[[id:sample_002, organism:mouse, tissue:brain, depth:25000000, quality:35.2, sample_num:2, lane:001, read:R1, chunk:001, priority:normal], /workspaces/training/side-quests/essential_scripting_patterns/data/sequences/SAMPLE_002_S2_L001_R1_001.fastq]
[[id:sample_003, organism:human, tissue:kidney, depth:45000000, quality:42.1, sample_num:3, lane:001, read:R1, chunk:001, priority:high], /workspaces/training/side-quests/essential_scripting_patterns/data/sequences/SAMPLE_003_S3_L001_R1_001.fastq]
```

### 2.2. Dynamic Script Generation in Processes

Process script blocks are essentially multi-line strings that get passed to the shell. You can use **conditional logic** (if/else, ternary operators) to dynamically generate different script strings based on input characteristics. This is essential for handling diverse input types—like single-end vs paired-end sequencing reads—without duplicating process definitions.

Let's add a process to our workflow that demonstrates this pattern. Open `modules/fastp.nf` and take a look:

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

    ```groovy title="main.nf" linenums="25" hl_lines="27"
    workflow {

        ch_samples = channel.fromPath("./data/samples.csv")
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
                return tuple(sample_meta + file_meta + [priority: priority], fastq_path)
            }

        ch_fastp = FASTP(ch_samples)
    }
    ```

=== "Before"

    ```groovy title="main.nf" linenums="25" hl_lines="26"
    workflow {

        ch_samples = channel.fromPath("./data/samples.csv")
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
                return [sample_meta + file_meta + [priority: priority], file(row.file_path)]
            }
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

Fix this by adding conditional logic to the `FASTP` process `script:` block. An if/else statement checks read file count and adjusts the command accordingly.

=== "After"

    ```groovy title="main.nf" linenums="10" hl_lines="3-27"
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

    ```groovy title="main.nf" linenums="10" hl_lines="2-11"
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

Now the workflow can handle both single-end and paired-end reads gracefully. The conditional logic checks the number of input files and constructs the appropriate command for `fastp`. Let's see if it works:

```bash title="Test dynamic fastp"
nextflow run main.nf
```

```console title="Successful run"
 N E X T F L O W   ~  version 25.04.3

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

Another common usage of dynamic script logic can be seen in [the Nextflow for Science Genomics module](../../nf4science/genomics/02_joint_calling). In that module, the GATK process being called can take multiple input files, but each must be prefixed with `-V` to form a correct command line. The process uses scripting to transform a collection of input files (`all_gvcfs`) into the correct command arguments:

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

These patterns of using scripting in process script blocks are extremely powerful and can be applied in many scenarios - from handling variable input types to building complex command-line arguments from file collections, making your processes truly adaptable to the diverse requirements of real-world data.

### 2.3. Variable Interpolation: Nextflow and Shell Variables

Process scripts mix Nextflow variables, shell variables, and command substitutions, each with different interpolation syntax. Using the wrong syntax causes errors. Let's explore these with a process that creates a processing report.

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

    ```groovy title="main.nf" linenums="1" hl_lines="2 30"
    include { FASTP } from './modules/fastp.nf'
    include { GENERATE_REPORT } from './modules/generate_report.nf'

    workflow {
        ch_samples = channel.fromPath("./data/samples.csv")
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
                return tuple(sample_meta + file_meta + [priority: priority], fastq_path)
            }

        ch_fastp = FASTP(ch_samples)
        GENERATE_REPORT(ch_samples)
    }
    ```

=== "Before"

    ```groovy title="main.nf" linenums="1" hl_lines="1 10-29"
    include { FASTP } from './modules/fastp.nf'

    workflow {
        ch_samples = channel.fromPath("./data/samples.csv")
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
                return tuple(sample_meta + file_meta + [priority: priority], fastq_path)
            }

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

If you run this, you'll notice an error or unexpected behavior - Nextflow tries to interpret `$(hostname)` as a Nextflow variable that doesn't exist:

```console title="Error with shell variables"
unknown recognition error type: groovyjarjarantlr4.v4.runtime.LexerNoViableAltException
ERROR ~ Module compilation error
- file : /workspaces/training/side-quests/essential_scripting_patterns/modules/generate_report.nf
- cause: token recognition error at: '(' @ line 16, column 22.
       echo "Hostname: $(hostname)" >> ${meta.id}_report.txt
                        ^

1 error
```

We need to escape it so Bash can handle it instead.

Fix this by escaping the shell variables and command substitutions with a backslash (`\`):

=== "After"

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

=== "Before"

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

In this section, you've learned **string processing** techniques:

- **Regular expressions for file parsing**: Using the `=~` operator and regex patterns (`~/pattern/`) to extract metadata from complex file naming conventions
- **Dynamic script generation**: Using conditional logic (if/else, ternary operators) to generate different script strings based on input characteristics
- **Variable interpolation**: Understanding when Nextflow interprets strings vs when the shell does
  - `${var}` - Nextflow variables (interpolated by Nextflow at workflow compile time)
  - `\${var}` - Shell environment variables (escaped, passed to bash at runtime)
  - `\$(cmd)` - Shell command substitution (escaped, executed by bash at runtime)

These string processing and generation patterns are essential for handling the diverse file formats and naming conventions you'll encounter in real-world bioinformatics workflows.

---

## 3. Creating Reusable Functions

Complex workflow logic inline in channel operators or process definitions reduces readability and maintainability. **Functions** let you extract this logic into named, reusable components.

Our map operation has grown long and complex. Let's extract it into a reusable function using the `def` keyword.

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
        return tuple(sample_meta + file_meta + [priority: priority], fastq_path)
    }

    workflow {
        ch_samples = channel.fromPath("./data/samples.csv")
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
        ch_samples = channel.fromPath("./data/samples.csv")
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
                return tuple(sample_meta + file_meta + [priority: priority], fastq_path)
            }

        ch_fastp = FASTP(ch_samples)
        GENERATE_REPORT(ch_samples)
    }
    ```

By extracting this logic into a function, we've reduced the actual workflow logic down to something much cleaner:

```groovy title="minimal workflow"
    ch_samples = channel.fromPath("./data/samples.csv")
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
 N E X T F L O W   ~  version 25.04.3

Launching `main.nf` [admiring_panini] DSL2 - revision: 8cc832e32f

executor >  local (6)
[8c/2e3f91] process > FASTP (3)           [100%] 3 of 3 ✔
[7a/1b4c92] process > GENERATE_REPORT (3) [100%] 3 of 3 ✔
```

The output should show both processes completing successfully. The workflow is now much cleaner and easier to maintain, with all the complex metadata processing logic encapsulated in the `separateMetadata` function.

### Takeaway

In this section, you've learned **function creation**:

- **Defining functions with `def`**: The keyword for creating named functions (like `def` in Python or `function` in JavaScript)
- **Function scope**: Functions defined at the script level are accessible throughout your Nextflow workflow
- **Return values**: Functions automatically return the last expression, or use explicit `return`
- **Cleaner code**: Extracting complex logic into functions is a fundamental software engineering practice in any language

Next, we'll explore how to use closures in process directives for dynamic resource allocation.

---

## 4. Dynamic Resource Directives with Closures

So far we've used scripting in the `script` block of processes. But **closures** (introduced in Section 1.1) are also incredibly useful in process directives, especially for dynamic resource allocation. Let's add resource directives to our FASTP process that adapt based on the sample characteristics.

Currently, our FASTP process uses default resources. Let's make it smarter by allocating more CPUs for high-depth samples. Edit `modules/fastp.nf` to include a dynamic `cpus` directive and a static `memory` directive:

=== "After"

    ```groovy title="modules/fastp.nf" linenums="1" hl_lines="4-5"
    process FASTP {
        container 'community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690'

        cpus { meta.depth > 40000000 ? 2 : 1 }
        memory 2.GB

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

The closure `{ meta.depth > 40000000 ? 2 : 1 }` uses the **ternary operator** (covered in Section 1.1) and is evaluated for each task, allowing per-sample resource allocation. High-depth samples (>40M reads) get 2 CPUs, while others get 1 CPU.

!!! note "Accessing Input Variables in Directives"

    The closure can access any input variables (like `meta` here) because Nextflow evaluates these closures in the context of each task execution.

Run the workflow again:

```bash title="Test resource allocation"
nextflow run main.nf -ansi-log false
```

We're using the `-ansi-log false` option to make it easier to see the task hashes.

```console title="Resource allocation output"
N E X T F L O W  ~  version 25.04.3
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
    docker run -i --cpu-shares 4096 --memory 2048m -e "NXF_TASK_WORKDIR" -v /workspaces/training/side-quests/essential_scripting_patterns:/workspaces/training/side-quests/essential_scripting_patterns -w "$NXF_TASK_WORKDIR" --name $NXF_BOXID community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690 /bin/bash -ue /workspaces/training/side-quests/essential_scripting_patterns/work/48/6db0c9e9d8aa65e4bb4936cd3bd59e/.command.sh
```

In this example we've chosen an example that requested 2 CPUs (`--cpu-shares 2048`), because it was a high-depth sample, but you should see different CPU allocations depending on the sample depth. Try this for the other tasks as well.

Another powerful pattern is using `task.attempt` for retry strategies. To show why this is useful, we're going to start by reducing the memory allocation to FASTP to less than it needs. Change the `memory` directive in `modules/fastp.nf` to `1.GB`:

=== "After"

    ```groovy title="modules/fastp.nf" linenums="1" hl_lines="5"
    process FASTP {
        container 'community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690'

        cpus { meta.depth > 40000000 ? 4 : 2 }
        memory 1.GB

        input:
        tuple val(meta), path(reads)
    ```

=== "Before"

    ```groovy title="modules/fastp.nf" linenums="1" hl_lines="5"
    process FASTP {
        container 'community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690'

        cpus { meta.depth > 40000000 ? 4 : 2 }
        memory 2.GB

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
        memory { 1.GB * task.attempt }
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
        memory 2.GB

        input:
        tuple val(meta), path(reads)
    ```

Now if the process fails due to insufficient memory, Nextflow will retry with more memory:

- First attempt: 1 GB (task.attempt = 1)
- Second attempt: 2.GB (task.attempt = 2)

... and so on, up to the `maxRetries` limit.

### Takeaway

Dynamic directives with closures let you:

- Allocate resources based on input characteristics
- Implement automatic retry strategies with increasing resources
- Combine multiple factors (metadata, attempt number, priorities)
- Use conditional logic for complex resource calculations

This makes your workflows both more efficient (not over-allocating) and more robust (automatic retry with more resources).

---

## 5. Conditional Logic and Process Control

Previously, we used `.map()` with scripting to transform channel data. Now we'll use conditional logic to control which processes execute based on data—essential for flexible workflows adapting to different sample types.

Nextflow's [dataflow operators](https://www.nextflow.io/docs/latest/reference/operator.html) take closures evaluated at runtime, enabling conditional logic to drive workflow decisions based on channel content.

### 5.1. Routing with `.branch()`

For example, let's pretend that our sequencing samples need to be trimmed with FASTP only if they're human samples with a coverage above a certain threshold. Mouse samples or low-coverage samples should be run with Trimgalore instead (this is a contrived example, but it illustrates the point).

We've provided a simple Trimgalore process in `modules/trimgalore.nf`, take a look if you like, but the details aren't important for this exercise. The key point is that we want to route samples based on their metadata.

Include the new from in `modules/trimgalore.nf`:

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
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map(separateMetadata)

        trim_branches = ch_samples
            .branch { meta, reads ->
                fastp: meta.organism == 'human' && meta.depth >= 30000000
                trimgalore: true
            }

        ch_fastp = FASTP(trim_branches.fastp)
        ch_trimgalore = TRIMGALORE(trim_branches.trimgalore)
        GENERATE_REPORT(ch_samples)
    ```

=== "Before"

    ```groovy title="main.nf" linenums="28" hl_lines="5"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map(separateMetadata)

        ch_fastp = FASTP(ch_samples)
        GENERATE_REPORT(ch_samples)
    ```

Run this modified workflow:

```bash title="Test conditional trimming"
nextflow run main.nf
```

```console title="Conditional trimming results"
 N E X T F L O W   ~  version 25.04.3

Launching `main.nf` [adoring_galileo] DSL2 - revision: c9e83aaef1

executor >  local (6)
[1d/0747ac] process > FASTP (2)           [100%] 2 of 2 ✔
[cc/c44caf] process > TRIMGALORE (1)      [100%] 1 of 1 ✔
[34/bd5a9f] process > GENERATE_REPORT (1) [100%] 3 of 3 ✔
```

Here, we've used small but mighty conditional expressions inside the `.branch{}` operator to route samples based on their metadata. Human samples with high coverage go through `FASTP`, while all other samples go through `TRIMGALORE`.

### 5.2. Using `.filter()` with Truthiness

Another powerful pattern for controlling workflow execution is the `.filter()` operator, which uses a closure to determine which items should continue down the pipeline. Inside the filter closure, you'll write **boolean expressions** that decide which items pass through.

Nextflow (like many dynamic languages) has a concept of **"truthiness"** that determines what values evaluate to `true` or `false` in boolean contexts:

- **Truthy**: Non-null values, non-empty strings, non-zero numbers, non-empty collections
- **Falsy**: `null`, empty strings `""`, zero `0`, empty collections `[]` or `[:]`, `false`

This means `meta.id` alone (without explicit `!= null`) checks if the ID exists and isn't empty. Let's use this to filter out samples that don't meet our quality requirements.

Add the following before the branch operation:

=== "After"

    ```groovy title="main.nf" linenums="28" hl_lines="5-11"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map(separateMetadata)

        // Filter out invalid or low-quality samples
        ch_valid_samples = ch_samples
            .filter { meta, reads ->
                meta.id && meta.organism && meta.depth >= 25000000
            }

        trim_branches = ch_valid_samples
            .branch { meta, reads ->
                fastp: meta.organism == 'human' && meta.depth >= 30000000
                trimgalore: true
            }
    ```

=== "Before"

    ```groovy title="main.nf" linenums="28" hl_lines="5"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map(separateMetadata)

        trim_branches = ch_samples
            .branch { meta, reads ->
                fastp: meta.organism == 'human' && meta.depth >= 30000000
                trimgalore: true
            }
    ```

Run the workflow again:

```bash title="Test filtering samples"
nextflow run main.nf
```

Because we've chosen a filter that excludes some samples, you should see fewer tasks executed:

```console title="Filtered samples results"
N E X T F L O W  ~  version 25.04.3
Launching `main.nf` [lonely_williams] DSL2 - revision: d0b3f121ec
[94/b48eac] Submitted process > FASTP (2)
[2c/d2b28f] Submitted process > GENERATE_REPORT (2)
[65/2e3be4] Submitted process > GENERATE_REPORT (1)
[94/b48eac] NOTE: Process `FASTP (2)` terminated with an error exit status (137) -- Execution is retried (1)
[3e/0d8664] Submitted process > TRIMGALORE (1)
[6a/9137b0] Submitted process > FASTP (1)
[6a/9137b0] NOTE: Process `FASTP (1)` terminated with an error exit status (137) -- Execution is retried (1)
[83/577ac0] Submitted process > GENERATE_REPORT (3)
[a2/5117de] Re-submitted process > FASTP (1)
[1f/a1a4ca] Re-submitted process > FASTP (2)
```

The filter expression `meta.id && meta.organism && meta.depth >= 25000000` combines truthiness with explicit comparisons:

- `meta.id && meta.organism` checks that both fields exist and are non-empty (using truthiness)
- `meta.depth >= 25000000` ensures sufficient sequencing depth with an explicit comparison

!!! note "Truthiness in Practice"

    The expression `meta.id && meta.organism` is more concise than writing:
    ```groovy
    meta.id != null && meta.id != '' && meta.organism != null && meta.organism != ''
    ```

    This makes filtering logic much cleaner and easier to read.

### Takeaway

In this section, you've learned to use conditional logic to control workflow execution using the closure interfaces of Nextflow operators like `.branch{}` and `.filter{}`, leveraging truthiness to write concise conditional expressions.

Our pipeline now intelligently routes samples through appropriate processes, but production workflows need to handle invalid data gracefully. Let's make our workflow robust against missing or null values.

---

## 6. Safe Navigation and Elvis Operators

Our `separateMetadata` function currently assumes all CSV fields are present and valid. But what happens with incomplete data? Let's find out.

### 6.1. The Problem: Accessing Properties That Don't Exist

Let's say we want to add support for optional sequencing run information. In some labs, samples might have an additional field for the sequencing run ID or batch number, but our current CSV doesn't have this column. Let's try to access it anyway.

Modify the `separateMetadata` function to include a run_id field:

=== "After"

    ```groovy title="main.nf" linenums="5" hl_lines="9"
    def separateMetadata(row) {
        def sample_meta = [
            id: row.sample_id.toLowerCase(),
            organism: row.organism,
            tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
            depth: row.sequencing_depth.toInteger(),
            quality: row.quality_score.toDouble()
        ]
        def run_id = row.run_id.toUpperCase()
    ```

=== "Before"

    ```groovy title="main.nf" linenums="5"
    def separateMetadata(row) {
        def sample_meta = [
            id: row.sample_id.toLowerCase(),
            organism: row.organism,
            tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
            depth: row.sequencing_depth.toInteger(),
            quality: row.quality_score.toDouble()
        ]
    ```

Now run the workflow:

```bash
nextflow run main.nf
```

It crashes with a NullPointerException:

```console title="Null pointer error"
 N E X T F L O W   ~  version 25.04.3

Launching `main.nf` [trusting_torvalds] DSL2 - revision: b56fbfbce2

ERROR ~ Cannot invoke method toUpperCase() on null object

 -- Check script 'main.nf' at line: 13 or see '.nextflow.log' file for more details
```

The problem is that `row.run_id` returns `null` because the `run_id` column doesn't exist in our CSV. When we try to call `.toUpperCase()` on `null`, it crashes. This is where the safe navigation operator saves the day.

### 6.2. Safe Navigation Operator (`?.`)

The safe navigation operator (`?.`) returns `null` instead of throwing an exception when called on a `null` value. If the object before `?.` is `null`, the entire expression evaluates to `null` without executing the method.

Update the function to use safe navigation:

=== "After"

    ```groovy title="main.nf" linenums="4" hl_lines="9"
    def separateMetadata(row) {
        def sample_meta = [
            id: row.sample_id.toLowerCase(),
            organism: row.organism,
            tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
            depth: row.sequencing_depth.toInteger(),
            quality: row.quality_score.toDouble()
        ]
        def run_id = row.run_id?.toUpperCase()
    ```

=== "Before"

    ```groovy title="main.nf" linenums="4" hl_lines="9"
    def separateMetadata(row) {
        def sample_meta = [
            id: row.sample_id.toLowerCase(),
            organism: row.organism,
            tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
            depth: row.sequencing_depth.toInteger(),
            quality: row.quality_score.toDouble()
        ]
        def run_id = row.run_id.toUpperCase()
    ```

Run again:

```bash
nextflow run main.nf
```

No crash! The workflow now handles the missing field gracefully. When `row.run_id` is `null`, the `?.` operator prevents the `.toUpperCase()` call, and `run_id` becomes `null` instead of causing an exception.

### 6.3. Elvis Operator (`?:`) for Defaults

The Elvis operator (`?:`) provides default values when the left side is "falsy" (as explained previously). It's named after Elvis Presley because `?:` looks like his famous hair and eyes when viewed sideways!

Now that we're using safe navigation, `run_id` will be `null` for samples without that field. Let's use the Elvis operator to provide a default value and add it to our `sample_meta` map:

=== "After"

    ```groovy title="main.nf" linenums="5" hl_lines="9-10"
    def separateMetadata(row) {
        def sample_meta = [
            id: row.sample_id.toLowerCase(),
            organism: row.organism,
            tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
            depth: row.sequencing_depth.toInteger(),
            quality: row.quality_score.toDouble()
        ]
        def run_id = row.run_id?.toUpperCase() ?: 'UNSPECIFIED'
        sample_meta.run = run_id
    ```

=== "Before"

    ```groovy title="main.nf" linenums="5" hl_lines="9"
    def separateMetadata(row) {
        def sample_meta = [
            id: row.sample_id.toLowerCase(),
            organism: row.organism,
            tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
            depth: row.sequencing_depth.toInteger(),
            quality: row.quality_score.toDouble()
        ]
        def run_id = row.run_id?.toUpperCase()
    ```

Also add a `view()` operator in the workflow to see the results:

=== "After"

    ```groovy title="main.nf" linenums="30" hl_lines="4"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map{ row -> separateMetadata(row) }
            .view()
    ```

=== "Before"

    ```groovy title="main.nf" linenums="30"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map{ row -> separateMetadata(row) }
    ```

and run the workflow:

```bash
nextflow run main.nf
```

You'll see output like this:

```console title="View output with run field"
[[id:sample_001, organism:human, tissue:liver, depth:30000000, quality:38.5, run:UNSPECIFIED, sample_num:1, lane:001, read:R1, chunk:001, priority:normal], /workspaces/training/side-quests/essential_scripting_patterns/data/sequences/SAMPLE_001_S1_L001_R1_001.fastq]
[[id:sample_002, organism:mouse, tissue:brain, depth:25000000, quality:35.2, run:UNSPECIFIED, sample_num:2, lane:001, read:R1, chunk:001, priority:normal], /workspaces/training/side-quests/essential_scripting_patterns/data/sequences/SAMPLE_002_S2_L001_R1_001.fastq]
[[id:sample_003, organism:human, tissue:kidney, depth:45000000, quality:42.1, run:UNSPECIFIED, sample_num:3, lane:001, read:R1, chunk:001, priority:high], /workspaces/training/side-quests/essential_scripting_patterns/data/sequences/SAMPLE_003_S3_L001_R1_001.fastq]
```

Perfect! Now all samples have a `run` field with either their actual run ID (in uppercase) or the default value 'UNSPECIFIED'. The combination of `?.` and `?:` provides both safety (no crashes) and sensible defaults.

Take out the `.view()` operator now that we've confirmed it works.

!!! tip "Combining Safe Navigation and Elvis"

    The pattern `value?.method() ?: 'default'` is common in production workflows:

    - `value?.method()` - Safely calls method, returns `null` if `value` is `null`
    - `?: 'default'` - Provides fallback if result is `null`

    This pattern handles missing/incomplete data gracefully.

Use these operators consistently in functions, operator closures (`.map{}`, `.filter{}`), process scripts, and config files. They prevent crashes when handling real-world data.

### Takeaway

- **Safe navigation (`?.`)**: Prevents crashes on null values - returns null instead of throwing exception
- **Elvis operator (`?:`)**: Provides defaults - `value ?: 'default'`
- **Combining**: `value?.method() ?: 'default'` is the common pattern

These operators make workflows resilient to incomplete data - essential for real-world work.

---

## 7. Validation with `error()` and `log.warn`

Sometimes you need to stop the workflow immediately if input parameters are invalid. In Nextflow, you can use built-in functions like `error()` and `log.warn`, as well as standard programming constructs like `if` statements and boolean logic, to implement validation logic. Let's add validation to our workflow.

Create a validation function before your workflow block, call it from the workflow, and change the channel creation to use a parameter for the CSV file path. If the parameter is missing or the file doesn't exist, call `error()` to stop execution with a clear message.

=== "After"

    ```groovy title="main.nf" linenums="1" hl_lines="5-20 23-24"
    include { FASTP } from './modules/fastp.nf'
    include { TRIMGALORE } from './modules/trimgalore.nf'
    include { GENERATE_REPORT } from './modules/generate_report.nf'

    def validateInputs() {
        // Check input parameter is provided
        if (!params.input) {
            error("Input CSV file path not provided. Please specify --input <file.csv>")
        }

        // Check CSV file exists
        if (!file(params.input).exists()) {
            error("Input CSV file not found: ${params.input}")
        }
    }
    ...
    workflow {
        validateInputs()
        ch_samples = channel.fromPath(params.input)
    ```

=== "Before"

    ```groovy title="main.nf" linenums="1"
    include { FASTP } from './modules/fastp.nf'
    include { TRIMGALORE } from './modules/trimgalore.nf'
    include { GENERATE_REPORT } from './modules/generate_report.nf'

    ...
    workflow {
        ch_samples = channel.fromPath("./data/samples.csv")
    ```

Now try running without the CSV file:

```bash
nextflow run main.nf
```

The workflow stops immediately with a clear error message instead of failing mysteriously later!

```console title="Validation error output"
 N E X T F L O W   ~  version 25.04.3

Launching `main.nf` [confident_coulomb] DSL2 - revision: 07059399ed

WARN: Access to undefined parameter `input` -- Initialise it to a default value eg. `params.input = some_value`
Input CSV file path not provided. Please specify --input <file.csv>
```

Now run it with a non-existent file:

```bash
nextflow run main.nf --input ./data/nonexistent.csv
```

Observe the error:

```console title="File not found error output"
 N E X T F L O W   ~  version 25.04.3

Launching `main.nf` [cranky_gates] DSL2 - revision: 26839ae3eb

Input CSV file not found: ./data/nonexistent.csv
```

Finally, run it with the correct file:

```bash
nextflow run main.nf --input ./data/samples.csv
```

This time it runs successfully.

You can also add validation within the `separateMetadata` function. Let's use the non-fatal `log.warn` to issue warnings for samples with low sequencing depth, but still allow the workflow to continue:

=== "After"

    ```groovy title="main.nf" linenums="1" hl_lines="3-6"
        def priority = sample_meta.quality > 40 ? 'high' : 'normal'

        // Validate data makes sense
        if (sample_meta.depth < 30000000) {
            log.warn "Low sequencing depth for ${sample_meta.id}: ${sample_meta.depth}"
        }

        return tuple(sample_meta + file_meta + [priority: priority], fastq_path)
    }
    ```

=== "Before"

    ```groovy title="main.nf" linenums="1"
        def priority = sample_meta.quality > 40 ? 'high' : 'normal'

        return tuple(sample_meta + file_meta + [priority: priority], fastq_path)
    }
    ```

Run the workflow again with the original CSV:

```bash
nextflow run main.nf --input ./data/samples.csv
```

... and you'll see a warning about low sequencing depth for one of the samples:

```console title="Warning output"
 N E X T F L O W   ~  version 25.04.3

Launching `main.nf` [awesome_goldwasser] DSL2 - revision: a31662a7c1

executor >  local (5)
[ce/df5eeb] process > FASTP (2)           [100%] 2 of 2 ✔
[-        ] process > TRIMGALORE          -
[d1/7d2b4b] process > GENERATE_REPORT (3) [100%] 3 of 3 ✔
WARN: Low sequencing depth for sample_002: 25000000
```

### Takeaway

- **`error()`**: Stops workflow immediately with clear message
- **`log.warn`**: Issues warnings without stopping workflow
- **Early validation**: Check inputs before processing to fail fast with helpful errors
- **Validation functions**: Create reusable validation logic that can be called at workflow start

Proper validation makes workflows more robust and user-friendly by catching problems early with clear error messages.

---

## 8. Workflow Event Handlers

Up until now, we've been writing code in our workflow scripts and process definitions. But there's one more important feature you should know about: workflow event handlers.

Event handlers are closures that run at specific points in your workflow's lifecycle. They're perfect for adding logging, notifications, or cleanup operations. These handlers should be defined in your workflow script alongside your workflow definition.

### 8.1. The `onComplete` Handler

The most commonly used event handler is `onComplete`, which runs when your workflow finishes (whether it succeeded or failed). Let's add one to summarize our pipeline results.

Add the event handler to your `main.nf` file, inside your workflow definition:

=== "After"

    ```groovy title="main.nf" linenums="66" hl_lines="5-16"
        ch_fastp = FASTP(trim_branches.fastp)
        ch_trimgalore = TRIMGALORE(trim_branches.trimgalore)
        GENERATE_REPORT(ch_samples)

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
    }
    ```

=== "Before"

    ```groovy title="main.nf" linenums="66" hl_lines="4"
        ch_fastp = FASTP(trim_branches.fastp)
        ch_trimgalore = TRIMGALORE(trim_branches.trimgalore)
        GENERATE_REPORT(ch_samples)
    }
    ```

This closure runs when the workflow completes. Inside, you have access to the `workflow` object which provides useful properties about the execution.

Run your workflow and you'll see this summary appear at the end!

```bash title="Run with onComplete handler"
nextflow run main.nf --input ./data/samples.csv -ansi-log false
```

```console title="onComplete output"
N E X T F L O W  ~  version 25.04.3
Launching `main.nf` [marvelous_boltzmann] DSL2 - revision: a31662a7c1
WARN: Low sequencing depth for sample_002: 25000000
[9b/d48e40] Submitted process > FASTP (2)
[6a/73867a] Submitted process > GENERATE_REPORT (2)
[79/ad0ac5] Submitted process > GENERATE_REPORT (1)
[f3/bda6cb] Submitted process > FASTP (1)
[34/d5b52f] Submitted process > GENERATE_REPORT (3)

Pipeline execution summary:
==========================
Completed at: 2025-10-10T12:14:24.885384+01:00
Duration    : 2.9s
Success     : true
workDir     : /workspaces/training/side-quests/essential_scripting_patterns/work
exit status : 0
```

Let's make it more useful by adding conditional logic:

=== "After"

    ```groovy title="main.nf" linenums="66" hl_lines="5-22"
        ch_fastp = FASTP(trim_branches.fastp)
        ch_trimgalore = TRIMGALORE(trim_branches.trimgalore)
        GENERATE_REPORT(ch_samples)

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
            } else {
                println "❌ Pipeline failed!"
                println "Error: ${workflow.errorMessage}"
            }
        }
    }
    ```

=== "Before"

    ```groovy title="main.nf" linenums="66" hl_lines="5-16"
        ch_fastp = FASTP(trim_branches.fastp)
        ch_trimgalore = TRIMGALORE(trim_branches.trimgalore)
        GENERATE_REPORT(ch_samples)

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
    }
    ```

Now we get an even more informative summary, including a success/failure message and the output directory if specified:

```console title="Enhanced onComplete output"
N E X T F L O W  ~  version 25.04.3
Launching `main.nf` [boring_linnaeus] DSL2 - revision: a31662a7c1
WARN: Low sequencing depth for sample_002: 25000000
[e5/242efc] Submitted process > FASTP (2)
[3b/74047c] Submitted process > GENERATE_REPORT (3)
[8a/7a57e6] Submitted process > GENERATE_REPORT (1)
[a8/b1a31f] Submitted process > GENERATE_REPORT (2)
[40/648429] Submitted process > FASTP (1)

Pipeline execution summary:
==========================
Completed at: 2025-10-10T12:16:00.522569+01:00
Duration    : 3.6s
Success     : true
workDir     : /workspaces/training/side-quests/essential_scripting_patterns/work
exit status : 0

✅ Pipeline completed successfully!
```

You can also write the summary to a file using file operations:

```groovy title="main.nf - Writing summary to file"
workflow {
    // ... your workflow code ...

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
}
```

### 8.2. The `onError` Handler

Besides `onComplete`, there is one other event handler you can use: `onError`, which runs only if the workflow fails:

```groovy title="main.nf - onError handler"
workflow {
    // ... your workflow code ...

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
}
```

You can use multiple handlers together in your workflow script:

```groovy title="main.nf - Combined handlers"
workflow {
    // ... your workflow code ...

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
}
```

### Takeaway

In this section, you've learned:

- **Event handler closures**: Closures in your workflow script that run at different lifecycle points
- **`onComplete` handler**: For execution summaries and result reporting
- **`onError` handler**: For error handling and logging failures
- **Workflow object properties**: Accessing `workflow.success`, `workflow.duration`, `workflow.errorMessage`, etc.

Event handlers show how you can use the full power of the Nextflow language within your workflow scripts to add sophisticated logging and notification capabilities.

---

## Summary

Congratulations, you made it!

Throughout this side quest, you've built a comprehensive sample processing pipeline that evolved from basic metadata handling to a sophisticated, production-ready workflow.
Each section built upon the previous one, demonstrating how programming constructs transform simple workflows into powerful data processing systems, with the following benefits:

- **Clearer code**: Understanding dataflow vs scripting helps you write more organized workflows
- **Robust handling**: Safe navigation and Elvis operators make workflows resilient to missing data
- **Flexible processing**: Conditional logic lets your workflows process different sample types appropriately
- **Adaptive resources**: Dynamic directives optimize resource usage based on input characteristics

This progression mirrors the real-world evolution of bioinformatics pipelines, from research prototypes handling a few samples to production systems processing thousands of samples across laboratories and institutions.
Every challenge you solved and pattern you learned reflects actual problems developers face when scaling Nextflow workflows.

Applying these patterns in your own work will enable you to build robust, production-ready workflows.

### Key patterns

1.  **Dataflow vs Scripting:** You learned to distinguish between dataflow operations (channel orchestration) and scripting (code that manipulates data), including the crucial differences between operations on different types like `collect` on Channel vs List.

    - Dataflow: channel orchestration

    ```groovy
    channel.fromPath('*.fastq').splitCsv(header: true)
    ```

    - Scripting: data processing on collections

    ```groovy
    sample_data.collect { it.toUpperCase() }
    ```

2.  **Advanced String Processing**: You mastered regular expressions for parsing file names, dynamic script generation in processes, and variable interpolation (Nextflow vs Bash vs Shell).

    - Pattern matching

    ```groovy
    filename =~ ~/^(\w+)_(\w+)_(\d+)\.fastq$/
    ```

    - Function with conditional return

    ```groovy
    def parseSample(filename) {
        def matcher = filename =~ pattern
        return matcher ? [valid: true, data: matcher[0]] : [valid: false]
    }
    ```

    - File collection to command arguments (in process script block)

    ```groovy
    script:
    def file_args = input_files.collect { file -> "--input ${file}" }.join(' ')
    """
    analysis_tool ${file_args} --output results.txt
    """
    ```

3.  **Creating Reusable Functions**: You learned to extract complex logic into named functions that can be called from channel operators, making workflows more readable and maintainable.

    - Define a named function

    ```groovy
    def separateMetadata(row) {
        def sample_meta = [ /* code hidden for brevity */ ]
        def fastq_path = file(row.file_path)
        def m = (fastq_path.name =~ /^(.+)_S(\d+)_L(\d{3})_(R[12])_(\d{3})\.fastq(?:\.gz)?$/)
        def file_meta = m ? [ /* code hidden for brevity */ ] : [:]
        def priority = sample_meta.quality > 40 ? 'high' : 'normal'

        return tuple(sample_meta + file_meta + [priority: priority], fastq_path)
    }
    ```

    - Call the named function in a workflow

    ```groovy
    workflow {
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map{ row -> separateMetadata(row) }

        ch_fastp = FASTP(ch_samples)
    }
    ```

4.  **Dynamic Resource Directives with Closures**: You explored using closures in process directives for adaptive resource allocation based on input characteristics.

    - Named closures and composition

    ```groovy
    def enrichData = normalizeId >> addQualityCategory >> addFlags
    def processor = generalFunction.curry(fixedParam)
    ```

    - Closures with scope access

    ```groovy
    def collectStats = { data -> stats.count++; return data }
    ```

5.  **Conditional Logic and Process Control**: You added intelligent routing using `.branch()` and `.filter()` operators, leveraging truthiness for concise conditional expressions.

    - Use `.branch()` to route data through different workflow branches

    ```groovy
    trim_branches = ch_samples
    .branch { meta, reads ->
        fastp: meta.organism == 'human' && meta.depth >= 30000000
        trimgalore: true
    }

    ch_fastp = FASTP(trim_branches.fastp)
    ch_trimgalore = TRIMGALORE(trim_branches.trimgalore)
    ```

    - Boolean evaluation with Groovy Truth

    ```groovy
    if (sample.files) println "Has files"
    ```

    - Use `filter()` to subset data with 'truthiness'

    ```groovy
    ch_valid_samples = ch_samples
        .filter { meta, reads ->
            meta.id && meta.organism && meta.depth >= 25000000
        }
    ```

6.  **Safe Navigation and Elvis Operators**: You made the pipeline robust against missing data using `?.` for null-safe property access and `?:` for providing default values.

    ```groovy
    def id = data?.sample?.id ?: 'unknown'
    ```

7.  **Validation with error() and log.warn**: You learned to validate inputs early and fail fast with clear error messages.

    ```groovy
    try {
        def errors = validateSample(sample)
        if (errors) throw new RuntimeException("Invalid: ${errors.join(', ')}")
    } catch (Exception e) {
        println "Error: ${e.message}"
    }
    ```

8.  **Configuration Event Handlers**: You learned to use workflow event handlers (`onComplete` and `onError`) for logging, notifications, and lifecycle management.

    - Using `onComplete` to log and notify

    ```groovy
    workflow.onComplete = {
        println "Success     : ${workflow.success}"
        println "exit status : ${workflow.exitStatus}"

        if (workflow.success) {
            println "✅ Pipeline completed successfully!"
        } else {
            println "❌ Pipeline failed!"
            println "Error: ${workflow.errorMessage}"
        }
    }
    ```

    - Using `onError` to take action specifically in case of failure

    ```groovy
    workflow.onError = {
        // Write detailed error log
        def error_file = file("${workflow.launchDir}/error.log")
        error_file.text = """
        Time: ${new Date()}
        Error: ${workflow.errorMessage}
        Error report: ${workflow.errorReport ?: 'No detailed report available'}
        """

        println "Error details written to: ${error_file}"
    }
    ```

### Additional resources

- [Nextflow Language Reference](https://nextflow.io/docs/latest/reference/syntax.html)
- [Nextflow Operators](https://www.nextflow.io/docs/latest/operator.html)
- [Nextflow Script Syntax](https://www.nextflow.io/docs/latest/script.html)
- [Nextflow Standard Library](https://nextflow.io/docs/latest/reference/stdlib.html)

Be sure to check out these resources when you need to explore more advanced features.

You'll benefit from practicing and expanding your skills in order to:

- Write cleaner workflows with proper separation between dataflow and scripting
- Master variable interpolation to avoid common pitfalls with Nextflow, Bash, and shell variables
- Use dynamic resource directives for efficient, adaptive workflows
- Transform file collections into properly formatted command-line arguments
- Handle different file naming conventions and input formats gracefully using regex and string processing
- Build reusable, maintainable code using advanced closure patterns and functional programming
- Process and organize complex datasets using collection operations
- Add validation, error handling, and logging to make your workflows production-ready
- Implement workflow lifecycle management with event handlers

---

## What's next?

Return to the [menu of Side Quests](./index.md) or click the button in the bottom right of the page to move on to the next topic in the list.
