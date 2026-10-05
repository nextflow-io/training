# Types and Records

When you're developing a pipeline, most of the data flowing through your channels has no declared shape.
A meta map is just a bag of keys, a tuple is just a list of values in some order, and nothing tells Nextflow which keys or positions a process expects.
That works fine until something changes: a collaborator renames a column in a samplesheet, or a module reads `meta.sample` when the map contains `meta.id`.
The pipeline usually keeps running, quietly substitutes `null`, and you only find out when someone reads the results.

Nextflow's **static typing** lets you declare what your data looks like, and **records** give your metadata named, checkable fields in place of maps and positional tuples.
Together they turn many of these silent mistakes into errors that point straight at the problem, often before the pipeline even runs.

### Learning goals

In this side quest, you'll take a small, working, untyped pipeline, watch it produce wrong results without complaint, and migrate it step by step to static types and records.

By the end of this side quest, you'll be able to:

- Enable static typing and convert processes to typed inputs and outputs
- Replace meta maps and tuples with records, and merge records with `+`
- Use destructured record inputs so a process takes only the fields it needs
- Declare record types and give named workflows typed `take:` and `emit:` blocks
- Explain the difference between a `Channel` and a `Value`
- Use `nextflow lint` to catch mistakes before you run the pipeline

### Prerequisites

Before taking on this side quest, you should:

- Have completed the [Hello Nextflow](../../hello_nextflow/index.md) tutorial or equivalent beginner's course.
- Be comfortable using basic Nextflow concepts and mechanisms (processes, channels, operators, modules)
- Be familiar with meta maps, as covered in the [Metadata and Meta Maps](../metadata/index.md) side quest

---

## 0. Get started

#### Open the training codespace

If you haven't yet done so, make sure to open the training environment as described in the [Environment Setup](../../envsetup/index.md).

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### Move into the project directory

Move into the directory where the files for this tutorial are located.

```bash
cd side-quests/types_and_records
```

You can set VSCode to focus on this directory:

```bash
code .
```

The editor opens with the project directory in focus.

#### Review the materials

You'll find a main workflow file, three process modules, and a `data` directory with two samplesheets:

```console title="Directory contents"
├── data
│   ├── people.csv               # The samplesheet the pipeline was written for
│   └── people_v2.csv            # A collaborator's updated samplesheet
├── main.nf                      # The (untyped) pipeline you will migrate
├── modules
│   ├── collect_greetings.nf     # Concatenates all greetings into one file
│   ├── say_hello.nf             # Writes a greeting for one person
│   └── shout.nf                 # Converts a greeting to uppercase
└── nextflow.config
```

The samplesheet lists a few people, the language they speak, and how to greet them:

```csv title="data/people.csv"
id,name,language,greeting
alice,alice,en,Hello
bruno,BRUNO,fr,Bonjour
chiara,Chiara,it,Ciao
dieter,dieter,de,Hallo
```

The pipeline reads each row into a meta map, writes a greeting file for each person, makes an uppercase copy of each greeting, and collects all the uppercase greetings into a single file.
The processes are deliberately tiny (`echo`, `tr` and `cat`), so you can focus on how the data flows between them.

#### Review the assignment

Your challenge is to migrate this pipeline to static types and records without changing what it produces, and to make it fail loudly when its input doesn't have the shape it expects.

Along the way, you will:

1. **Watch the untyped pipeline fail silently** on a samplesheet with a renamed column
2. **Enable static typing** and convert a process to typed inputs and outputs
3. **Replace meta maps and tuples with records**, so the same bad samplesheet fails immediately
4. **Wrap the processing in a typed workflow** with a clear, declared interface
5. **Let the type checker catch mistakes** before you run anything

#### Readiness checklist

Think you're ready to dive in?

- [ ] I understand the goal of this course and its prerequisites
- [ ] My codespace is up and running
- [ ] I've set my working directory appropriately
- [ ] I understand the assignment

If you can check all the boxes, you're good to go.

---

## 1. Meet the pipeline and its blind spot

Before changing anything, run the pipeline as it is, then see what happens when its input changes under it.

### 1.1. Run the starter pipeline

Open `main.nf` and take a look at the code:

```groovy title="main.nf" linenums="1"
#!/usr/bin/env nextflow

include { SAY_HELLO } from './modules/say_hello.nf'
include { SHOUT } from './modules/shout.nf'
include { COLLECT_GREETINGS } from './modules/collect_greetings.nf'

params {
    input: Path = 'data/people.csv'
    batch: String = 'batch'
}

workflow {
    main:
    // Read the CSV into a channel of meta maps, one per person
    people = channel.fromPath(params.input)
        .splitCsv(header: true)
        .map { row ->
            [id: row.id, name: row.name, language: row.language, greeting: row.greeting]
        }

    greetings = SAY_HELLO(people)
    shouted = SHOUT(greetings)
    collected = COLLECT_GREETINGS(shouted.map { _meta, file -> file }.collect(), params.batch)

    publish:
    greetings = greetings
    shouted = shouted
    collected = collected
}

output {
    greetings {
        path 'greetings'
    }
    shouted {
        path 'shouted'
    }
    collected {
    }
}
```

This is the same structure you've seen in Hello Nextflow and the metadata side quest: a `params` block, a CSV read with `splitCsv`, a meta map per row, three processes, and a `publish:` section with a matching `output {}` block.

Now look at `modules/say_hello.nf`:

```groovy title="modules/say_hello.nf" linenums="1"
/*
 * Write a personalised greeting to a file
 */
process SAY_HELLO {
    tag "${meta.id}"

    input:
    val meta

    output:
    tuple val(meta), path("${meta.id}.txt")

    script:
    """
    echo '${meta.greeting}, ${meta.name}!' > ${meta.id}.txt
    """
}
```

The process takes the whole meta map and reaches into it for `meta.greeting` and `meta.name`.
Its output is a `[meta, file]` tuple, which `SHOUT` (in `modules/shout.nf`) receives as `tuple val(meta), path(greeting_file)`.

Run the pipeline:

```bash
nextflow run main.nf
```

??? success "Command output"

    ```console
     N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [curious_curie] revision: e3a24f8e11

    executor >  local (9)
    [fb/398978] SAY_HELLO (chiara) | 4 of 4 ✔
    [4b/abf3c0] SHOUT (chiara)     | 4 of 4 ✔
    [7b/5cd51f] COLLECT_GREETINGS  | 1 of 1 ✔

    Outputs:

      /workspaces/training/side-quests/types_and_records/results

      greetings:
        - [{id: chiara, name: Chiara, language: it, greeting: Ciao}, greetings/chiara.txt]
        - [{id: alice, name: alice, language: en, greeting: Hello}, greetings/alice.txt]
        - [{id: bruno, name: BRUNO, language: fr, greeting: Bonjour}, greetings/bruno.txt]
        - [{id: dieter, name: dieter, language: de, greeting: Hallo}, greetings/dieter.txt]

      shouted:
        - [{id: bruno, name: BRUNO, language: fr, greeting: Bonjour}, shouted/bruno-shouted.txt]
        - [{id: alice, name: alice, language: en, greeting: Hello}, shouted/alice-shouted.txt]
        - [{id: dieter, name: dieter, language: de, greeting: Hallo}, shouted/dieter-shouted.txt]
        - [{id: chiara, name: Chiara, language: it, greeting: Ciao}, shouted/chiara-shouted.txt]

      collected: COLLECTED-batch.txt
    ```

Check two of the greetings it wrote:

```bash
cat results/greetings/alice.txt results/greetings/bruno.txt
```

??? abstract "File contents"

    ```console
    Hello, alice!
    Bonjour, BRUNO!
    ```

The pipeline works, although the names come out exactly as they were typed in the samplesheet, in whatever case someone happened to use.
You'll tidy that up later.

### 1.2. Run it on a collaborator's samplesheet

A collaborator has exported the same people from a different system and sent you `data/people_v2.csv`.
Take a look at it:

```csv title="data/people_v2.csv"
id,name,language,salutation
alice,alice,en,Hello
bruno,BRUNO,fr,Bonjour
chiara,Chiara,it,Ciao
dieter,dieter,de,Hallo
```

The other system calls the greeting column `salutation`.
It's an easy difference to miss, so run the pipeline on the new file and see what happens:

```bash
nextflow run main.nf --input data/people_v2.csv
```

??? success "Command output"

    ```console
     N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [cranky_goldwasser] revision: e3a24f8e11

    executor >  local (9)
    [4c/066b44] SAY_HELLO (chiara) | 4 of 4 ✔
    [88/0eb8b7] SHOUT (bruno)      | 4 of 4 ✔
    [a6/b3b523] COLLECT_GREETINGS  | 1 of 1 ✔

    Outputs:

      /workspaces/training/side-quests/types_and_records/results

      greetings:
        - [{id: chiara, name: Chiara, language: it, greeting: null}, greetings/chiara.txt]
        - [{id: alice, name: alice, language: en, greeting: null}, greetings/alice.txt]
        - [{id: bruno, name: BRUNO, language: fr, greeting: null}, greetings/bruno.txt]
        - [{id: dieter, name: dieter, language: de, greeting: null}, greetings/dieter.txt]

      shouted:
        - [{id: alice, name: alice, language: en, greeting: null}, shouted/alice-shouted.txt]
        - [{id: chiara, name: Chiara, language: it, greeting: null}, shouted/chiara-shouted.txt]
        - [{id: dieter, name: dieter, language: de, greeting: null}, shouted/dieter-shouted.txt]
        - [{id: bruno, name: BRUNO, language: fr, greeting: null}, shouted/bruno-shouted.txt]

      collected: COLLECTED-batch.txt
    ```

Every task succeeded and every output was published.
Now look at what was actually written:

```bash
cat results/greetings/alice.txt
```

??? abstract "File contents"

    ```console
    null, alice!
    ```

The greeting is `null`.
The row has no `greeting` column, so `row.greeting` quietly returned `null`, the meta map stored it, and `SAY_HELLO` interpolated it into the command as the text `null`.
Nothing in the pipeline knew that `greeting` was required, so nothing complained.

Positional tuples have the same weakness.
`SHOUT` only works because its `tuple val(meta), path(greeting_file)` input lists its elements in the same order as the tuple that `SAY_HELLO` emits, and `#!groovy shouted.map { _meta, file -> file }` relies on that order too.
Swap two elements and you get a confusing failure, or worse, a wrong result.

What's missing is a way to tell Nextflow what shape the data should have, so it can check that shape for you.

### Takeaway

In this section, you've seen how an untyped pipeline fails:

- **Meta maps accept anything**: a missing key gives `null`, not an error
- **Tuples rely on position**: every producer and consumer must agree on the order of elements
- **The failure surfaces late**: the run looks successful, and the problem only shows up in the results

Static typing is how you give Nextflow that knowledge.
In the next section, you'll switch it on and convert the first process.

---

## 2. Enable static typing

Static typing in Nextflow is opt-in, one script at a time, through a feature flag.
In this section you'll turn it on for `main.nf` and for one module, fix the problems it reveals, and meet `nextflow lint`.

### 2.1. Turn on static typing in the main script

Add the feature flag near the top of `main.nf`:

=== "After"

    ```groovy title="main.nf" linenums="1" hl_lines="3"
    #!/usr/bin/env nextflow

    nextflow.enable.types = true

    include { SAY_HELLO } from './modules/say_hello.nf'
    include { SHOUT } from './modules/shout.nf'
    include { COLLECT_GREETINGS } from './modules/collect_greetings.nf'
    ```

=== "Before"

    ```groovy title="main.nf" linenums="1"
    #!/usr/bin/env nextflow

    include { SAY_HELLO } from './modules/say_hello.nf'
    include { SHOUT } from './modules/shout.nf'
    include { COLLECT_GREETINGS } from './modules/collect_greetings.nf'
    ```

The flag applies only to the script that declares it.
The three modules are separate scripts, so they are still untyped.

Run the pipeline:

```bash
nextflow run main.nf
```

??? failure "Command output"

    ```console
     N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [naughty_joliot] revision: 9e571fe973

    WARN: Static typing is a preview feature -- syntax and behavior may change in future releases

    Missing process or function fromPath([data/people.csv])

     -- Check script 'main.nf' at line: 17 or see '.nextflow.log' file for more details
    ```

Two things happened.
First, Nextflow printed a warning that static typing is a preview feature.
You'll see that warning on every run from now on; it's expected on this version of Nextflow.

Second, the pipeline failed before running a single task.
With static typing on, Nextflow checks the arguments of functions such as `channel.fromPath()` against their declared types.
`channel.fromPath()` expects a file pattern (a `String`), but the `params` block declares `input` as a `Path`, so there is no matching version of the function to call.
This is static typing doing its job: the parameter is already a file, so there is nothing for `fromPath` to resolve.

### 2.2. Read the samplesheet the typed way

The idiomatic way to read a samplesheet in a typed script is to put the `Path` into a channel with `channel.of()`, and then split each file into rows inside `flatMap`:

=== "After"

    ```groovy title="main.nf" linenums="17" hl_lines="1 2"
        people = channel.of(params.input)
            .flatMap { csv -> csv.splitCsv(header: true) }
    ```

=== "Before"

    ```groovy title="main.nf" linenums="17" hl_lines="1 2"
        people = channel.fromPath(params.input)
            .splitCsv(header: true)
    ```

Here, `splitCsv` is called as a method on the CSV file itself, which returns a list of rows, and `flatMap` emits each row as a separate channel element.
The result is the same channel of rows as before, but every step now has a type that Nextflow can follow.
Avoid the `splitCsv` _operator_ in typed code, because newer versions of Nextflow warn that it is discouraged with static typing.

Run the pipeline again:

```bash
nextflow run main.nf
```

??? success "Command output"

    ```console
     N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [nauseous_cray] revision: a900104c7f

    WARN: Static typing is a preview feature -- syntax and behavior may change in future releases

    executor >  local (9)
    [7f/5bcf30] SAY_HELLO (bruno) | 4 of 4 ✔
    [e9/82ff1f] SHOUT (bruno)     | 4 of 4 ✔
    [46/529836] COLLECT_GREETINGS | 1 of 1 ✔

    Outputs:

      /workspaces/training/side-quests/types_and_records/results

      greetings:
        - [{id: bruno, name: BRUNO, language: fr, greeting: Bonjour}, greetings/bruno.txt]
        - [{id: chiara, name: Chiara, language: it, greeting: Ciao}, greetings/chiara.txt]
        - [{id: alice, name: alice, language: en, greeting: Hello}, greetings/alice.txt]
        - [{id: dieter, name: dieter, language: de, greeting: Hallo}, greetings/dieter.txt]

      shouted:
        - [{id: dieter, name: dieter, language: de, greeting: Hallo}, shouted/dieter-shouted.txt]
        - [{id: chiara, name: Chiara, language: it, greeting: Ciao}, shouted/chiara-shouted.txt]
        - [{id: alice, name: alice, language: en, greeting: Hello}, shouted/alice-shouted.txt]
        - [{id: bruno, name: BRUNO, language: fr, greeting: Bonjour}, shouted/bruno-shouted.txt]

      collected: COLLECTED-batch.txt
    ```

The pipeline runs all the way through on the original samplesheet, so the files in `results/` hold the correct greetings again.

Notice that `main.nf` is now typed, but the three processes it includes are still written with the legacy `val`/`path`/`tuple` syntax.
A typed script can include and call untyped modules, which means you can migrate a real pipeline gradually, one file at a time, instead of all at once.

### 2.3. Convert a process to typed inputs and outputs

Start with the simplest process: `COLLECT_GREETINGS`, which takes a list of files and a batch name and produces one file.

#### 2.3.1. Turn on static typing in the module

Add the same feature flag at the top of `modules/collect_greetings.nf`:

=== "After"

    ```groovy title="modules/collect_greetings.nf" linenums="1" hl_lines="1"
    nextflow.enable.types = true

    /*
     * Collect all greetings into a single file
    ```

=== "Before"

    ```groovy title="modules/collect_greetings.nf" linenums="1"
    /*
     * Collect all greetings into a single file
    ```

Run the pipeline:

```bash
nextflow run main.nf -resume
```

??? failure "Command output"

    ```console
     N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [condescending_mccarthy] revision: a900104c7f

    Error modules/collect_greetings.nf:9:5: Invalid input declaration in typed process
    │   9 |     path input_files
    ╰     |     ^^^^^^^^^^^^^^^^

    ERROR ~ Script compilation failed

     -- Check '.nextflow.log' file for details
    ```

This time the script doesn't even compile.
Once a script enables static typing, every process in it must use the typed syntax, and `path input_files` is a legacy input declaration.
Typed and legacy processes can live in the same pipeline, but not in the same file.

#### 2.3.2. Rewrite the inputs and outputs

In a typed process, each input is a name followed by a type, and each output is a value:

=== "After"

    ```groovy title="modules/collect_greetings.nf" linenums="8" hl_lines="2 3 6 10"
        input:
        input_files: Bag<Path>
        batch_name: String

        output:
        file("COLLECTED-${batch_name}.txt")

        script:
        """
        cat ${input_files.join(' ')} > COLLECTED-${batch_name}.txt
        """
    ```

=== "Before"

    ```groovy title="modules/collect_greetings.nf" linenums="8" hl_lines="2 3 6 10"
        input:
        path input_files
        val batch_name

        output:
        path "COLLECTED-${batch_name}.txt"

        script:
        """
        cat ${input_files} > COLLECTED-${batch_name}.txt
        """
    ```

Here's what changed:

- `val batch_name` becomes `batch_name: String`, an input with a declared type.
- `path input_files` becomes `input_files: Bag<Path>`. `Path` inputs are staged into the task directory just like `path` inputs, and `Bag<Path>` says the input is an unordered collection of files, which is exactly what `collect()` produces.
- The output is now a plain value: `#!groovy file("COLLECTED-${batch_name}.txt")` returns the file of that name from the task directory. A process with a single output doesn't need to name it.
- `input_files` is now a real collection rather than a staged file list, so the script joins the file names with spaces explicitly.

#### 2.3.3. Run the pipeline

Run the pipeline to check that the converted process still produces the same result:

```bash
nextflow run main.nf -resume
```

??? success "Command output"

    ```console
     N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [hopeful_varahamihira] revision: a900104c7f

    WARN: Static typing is a preview feature -- syntax and behavior may change in future releases

    executor >  local (1)
    [2c/5bfadb] SAY_HELLO (dieter) | 4 of 4, cached: 4 ✔
    [66/1e7c57] SHOUT (alice)      | 4 of 4, cached: 4 ✔
    [73/1c891d] COLLECT_GREETINGS  | 1 of 1 ✔

    Outputs:

      /workspaces/training/side-quests/types_and_records/results

      greetings:
        - [{id: chiara, name: Chiara, language: it, greeting: Ciao}, greetings/chiara.txt]
        - [{id: dieter, name: dieter, language: de, greeting: Hallo}, greetings/dieter.txt]
        - [{id: bruno, name: BRUNO, language: fr, greeting: Bonjour}, greetings/bruno.txt]
        - [{id: alice, name: alice, language: en, greeting: Hello}, greetings/alice.txt]

      shouted:
        - [{id: dieter, name: dieter, language: de, greeting: Hallo}, shouted/dieter-shouted.txt]
        - [{id: alice, name: alice, language: en, greeting: Hello}, shouted/alice-shouted.txt]
        - [{id: bruno, name: BRUNO, language: fr, greeting: Bonjour}, shouted/bruno-shouted.txt]
        - [{id: chiara, name: Chiara, language: it, greeting: Ciao}, shouted/chiara-shouted.txt]

      collected: COLLECTED-batch.txt
    ```

`COLLECT_GREETINGS` ran again because its definition changed, while `SAY_HELLO` and `SHOUT` stayed cached.

### 2.4. Check the code with `nextflow lint`

Running the pipeline is a slow way to find out whether your code is valid.
`nextflow lint` parses and checks every script and config file in a directory without running anything:

```bash
nextflow lint .
```

??? success "Command output"

    ```console
    Linting Nextflow code..
    Linting: modules/shout.nf
    Linting: modules/say_hello.nf
    Linting: modules/collect_greetings.nf
    Linting: nextflow.config
    Linting: main.nf
    Nextflow linting complete!
     ✅ 5 files had no errors
    ```

All files pass.
If you run `nextflow lint` on the module from step 2.3.1, before you rewrote it, it reports the same `Invalid input declaration` error you got from `nextflow run`, along with every other line that would need to change.
Get into the habit of running it after every change: it's fast, and as your code gains more type information, it can catch more mistakes.

### Takeaway

In this section, you've learned:

- **How to enable static typing**: `nextflow.enable.types = true`, in each script that uses typed code
- **How to migrate incrementally**: a typed script can include untyped modules, but each file is either typed or untyped
- **How to write a typed process**: inputs are `name: Type`, and outputs are values such as `file(...)`
- **How to check code without running it**: `nextflow lint`

You now have a typed main script and one typed process.
But the process that actually caused the `null` greeting, `SAY_HELLO`, still takes an untyped meta map.
Typing it as `meta: Map` wouldn't help much, because a `Map` can still hold any keys or none.
What you need is a data structure whose fields Nextflow knows about, and that's what records are.

---

## 3. Replace meta maps and tuples with records

A **record** is a set of named fields, like a meta map, but designed for static typing.
Records can hold files and values side by side, so they replace both the meta map and the `[meta, file]` tuple, and processes can declare exactly which fields they need.

In this section you'll switch the pipeline from meta maps to records, and then rerun the collaborator's samplesheet to see what changes.

### 3.1. Create records instead of meta maps

Start where the data is created.
Replace the map literal with a call to `record()`, and add a `view()` so you can see the result:

=== "After"

    ```groovy title="main.nf" linenums="16" hl_lines="1 5 8"
        // Read the CSV into a channel of records, one per person
        people = channel.of(params.input)
            .flatMap { csv -> csv.splitCsv(header: true) }
            .map { row ->
                record(id: row.id, name: row.name, language: row.language, greeting: row.greeting)
            }

        people.view()

        greetings = SAY_HELLO(people)
    ```

=== "Before"

    ```groovy title="main.nf" linenums="16" hl_lines="1 5"
        // Read the CSV into a channel of meta maps, one per person
        people = channel.of(params.input)
            .flatMap { csv -> csv.splitCsv(header: true) }
            .map { row ->
                [id: row.id, name: row.name, language: row.language, greeting: row.greeting]
            }

        greetings = SAY_HELLO(people)
    ```

The `record()` function takes the same `name: value` pairs as a map literal.

Run the pipeline:

```bash
nextflow run main.nf -resume
```

??? success "Command output"

    ```console
     N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [serene_dijkstra] revision: a649ea99fd

    WARN: Static typing is a preview feature -- syntax and behavior may change in future releases

    [id:alice, name:alice, language:en, greeting:Hello]
    [id:bruno, name:BRUNO, language:fr, greeting:Bonjour]
    [id:chiara, name:Chiara, language:it, greeting:Ciao]
    [id:dieter, name:dieter, language:de, greeting:Hallo]

    [2c/5bfadb] SAY_HELLO (dieter) | 4 of 4, cached: 4 ✔
    [06/2790fd] SHOUT (chiara)     | 4 of 4, cached: 4 ✔
    [73/1c891d] COLLECT_GREETINGS  | 1 of 1, cached: 1 ✔

    Outputs:

      /workspaces/training/side-quests/types_and_records/results

      greetings:
        - [{id: bruno, name: BRUNO, language: fr, greeting: Bonjour}, greetings/bruno.txt]
        - [{id: chiara, name: Chiara, language: it, greeting: Ciao}, greetings/chiara.txt]
        - [{id: dieter, name: dieter, language: de, greeting: Hallo}, greetings/dieter.txt]
        - [{id: alice, name: alice, language: en, greeting: Hello}, greetings/alice.txt]

      shouted:
        - [{id: bruno, name: BRUNO, language: fr, greeting: Bonjour}, shouted/bruno-shouted.txt]
        - [{id: chiara, name: Chiara, language: it, greeting: Ciao}, shouted/chiara-shouted.txt]
        - [{id: dieter, name: dieter, language: de, greeting: Hallo}, shouted/dieter-shouted.txt]
        - [{id: alice, name: alice, language: en, greeting: Hello}, shouted/alice-shouted.txt]

      collected: COLLECTED-batch.txt
    ```

A record prints just like a map, and you read its fields the same way, with dot notation (`person.name`).
That's why the untyped `SAY_HELLO` module keeps working unchanged: `meta.greeting` and `meta.name` resolve against the record's fields.

### 3.2. Merge records with `+`

The names in section 1.1 came out with inconsistent capitalisation (`alice`, `BRUNO`).
Fixing that introduces one of the most useful things you can do with records: merging them with `+`.

#### 3.2.1. Normalise the names

Add a second `map` that replaces each person's `name` with a tidied-up version:

=== "After"

    ```groovy title="main.nf" linenums="17" hl_lines="6 7 8"
        people = channel.of(params.input)
            .flatMap { csv -> csv.splitCsv(header: true) }
            .map { row ->
                record(id: row.id, name: row.name, language: row.language, greeting: row.greeting)
            }
            .map { person ->
                person + record(name: person.name.toLowerCase().capitalize())
            }
    ```

=== "Before"

    ```groovy title="main.nf" linenums="17"
        people = channel.of(params.input)
            .flatMap { csv -> csv.splitCsv(header: true) }
            .map { row ->
                record(id: row.id, name: row.name, language: row.language, greeting: row.greeting)
            }
    ```

The expression `#!groovy person + record(name: ...)` combines two records into a new one:

- Fields that appear on only one side are copied into the result.
- Fields that appear on both sides take the value from the **right-hand** record.

So the result has every field of `person`, with `name` replaced by the lowercased, capitalised version.
Records are never modified in place: `+` always builds a new record and leaves `person` untouched.

The same rule lets you add fields.
For example, `#!groovy record(a: 1, b: 'x') + record(b: 'y', c: true)` produces a record with `a: 1`, `b: 'y'` and `c: true`.
Think of the left side as the defaults and the right side as the overrides.

#### 3.2.2. Run the pipeline

Run the pipeline to see the merged records:

```bash
nextflow run main.nf -resume
```

??? success "Command output"

    ```console
     N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [sick_curie] revision: a0e9e62298

    WARN: Static typing is a preview feature -- syntax and behavior may change in future releases

    [id:alice, name:Alice, language:en, greeting:Hello]
    [id:bruno, name:Bruno, language:fr, greeting:Bonjour]
    [id:chiara, name:Chiara, language:it, greeting:Ciao]
    [id:dieter, name:Dieter, language:de, greeting:Hallo]

    executor >  local (7)
    [89/3233e6] SAY_HELLO (dieter) | 4 of 4, cached: 1 ✔
    [60/5347c4] SHOUT (dieter)     | 4 of 4, cached: 1 ✔
    [85/a37a48] COLLECT_GREETINGS  | 1 of 1 ✔

    Outputs:

      /workspaces/training/side-quests/types_and_records/results

      greetings:
        - [{id: chiara, name: Chiara, language: it, greeting: Ciao}, greetings/chiara.txt]
        - [{id: bruno, name: Bruno, language: fr, greeting: Bonjour}, greetings/bruno.txt]
        - [{id: alice, name: Alice, language: en, greeting: Hello}, greetings/alice.txt]
        - [{id: dieter, name: Dieter, language: de, greeting: Hallo}, greetings/dieter.txt]

      shouted:
        - [{id: chiara, name: Chiara, language: it, greeting: Ciao}, shouted/chiara-shouted.txt]
        - [{id: alice, name: Alice, language: en, greeting: Hello}, shouted/alice-shouted.txt]
        - [{id: bruno, name: Bruno, language: fr, greeting: Bonjour}, shouted/bruno-shouted.txt]
        - [{id: dieter, name: Dieter, language: de, greeting: Hallo}, shouted/dieter-shouted.txt]

      collected: COLLECTED-batch.txt
    ```

Every name is now capitalised the same way, and `name` is still in the same place in each record: the merged value replaced the old one rather than being added alongside it.

`SAY_HELLO` and `SHOUT` ran again for three of the four people, because their records changed.
`Chiara` was already capitalised, so her record is identical to before and her tasks came from the cache.

### 3.3. Give `SAY_HELLO` a record input

Now for the process that wrote `null, alice!`.
Instead of accepting any value, `SAY_HELLO` will declare exactly which fields it needs, with their types.

Replace the contents of `modules/say_hello.nf`:

=== "After"

    ```groovy title="modules/say_hello.nf" linenums="1" hl_lines="1 7 10 13 17"
    nextflow.enable.types = true

    /*
     * Write a personalised greeting to a file
     */
    process SAY_HELLO {
        tag "${id}"

        input:
        record(id: String, name: String, greeting: String)

        output:
        record(id: id, greeting_file: file("${id}.txt"))

        script:
        """
        echo '${greeting}, ${name}!' > ${id}.txt
        """
    }
    ```

=== "Before"

    ```groovy title="modules/say_hello.nf" linenums="1" hl_lines="5 8 11 15"
    /*
     * Write a personalised greeting to a file
     */
    process SAY_HELLO {
        tag "${meta.id}"

        input:
        val meta

        output:
        tuple val(meta), path("${meta.id}.txt")

        script:
        """
        echo '${meta.greeting}, ${meta.name}!' > ${meta.id}.txt
        """
    }
    ```

The input is a _destructured_ record: `record(id: String, name: String, greeting: String)` says "this process takes a record, and from it I need `id`, `name` and `greeting`, all strings".
Each field becomes a variable in the process, so the script refers to `#!groovy ${greeting}` and `#!groovy ${name}` directly.

Destructuring is how a process selects a subset of a larger record.
The records in `people` also carry a `language` field, which `SAY_HELLO` doesn't need.
It simply ignores it, and the process would work with any record that has these three fields.

The output is a record too.
`#!groovy record(id: id, greeting_file: file("${id}.txt"))` names each output field, so downstream code asks for `greeting_file` by name instead of remembering that the file is the second element of a tuple.

### 3.4. Give `SHOUT` a record input

`SHOUT` gets the same treatment.
It only needs the person's `id` and the greeting file, so that's all it asks for.

Replace the contents of `modules/shout.nf`:

=== "After"

    ```groovy title="modules/shout.nf" linenums="1" hl_lines="1 7 10 13 17"
    nextflow.enable.types = true

    /*
     * Convert a greeting to uppercase
     */
    process SHOUT {
        tag "${id}"

        input:
        record(id: String, greeting_file: Path)

        output:
        record(id: id, shouted: file("${id}-shouted.txt"))

        script:
        """
        tr '[:lower:]' '[:upper:]' < ${greeting_file} > ${id}-shouted.txt
        """
    }
    ```

=== "Before"

    ```groovy title="modules/shout.nf" linenums="1" hl_lines="5 8 11 15"
    /*
     * Convert a greeting to uppercase
     */
    process SHOUT {
        tag "${meta.id}"

        input:
        tuple val(meta), path(greeting_file)

        output:
        tuple val(meta), path("${meta.id}-shouted.txt")

        script:
        """
        tr '[:lower:]' '[:upper:]' < ${greeting_file} > ${meta.id}-shouted.txt
        """
    }
    ```

Because `greeting_file` is declared as a `Path`, Nextflow stages it into the task directory, exactly as `path(greeting_file)` did in the tuple.
The difference is that `SHOUT` now matches its input by field name, so it no longer cares where in the record the file sits.

### 3.5. Update the workflow

Back in `main.nf`, remove the `view()` and update the one place that unpacked a tuple by position:

=== "After"

    ```groovy title="main.nf" linenums="26" hl_lines="3"
        greetings = SAY_HELLO(people)
        shouted = SHOUT(greetings)
        collected = COLLECT_GREETINGS(shouted.map { r -> r.shouted }.collect(), params.batch)
    ```

=== "Before"

    ```groovy title="main.nf" linenums="26" hl_lines="1 5"
        people.view()

        greetings = SAY_HELLO(people)
        shouted = SHOUT(greetings)
        collected = COLLECT_GREETINGS(shouted.map { _meta, file -> file }.collect(), params.batch)
    ```

`#!groovy shouted.map { r -> r.shouted }` picks the `shouted` file out of each record by name.
The `output {}` block doesn't need to change: when you publish a channel of records, Nextflow publishes every file field they contain.

Run the pipeline:

```bash
nextflow run main.nf -resume
```

??? success "Command output"

    ```console
     N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [elegant_kilby] revision: 80bd6b48ff

    WARN: Static typing is a preview feature -- syntax and behavior may change in future releases

    executor >  local (9)
    [95/13cf02] SAY_HELLO (alice)  | 4 of 4 ✔
    [ad/bbcd55] SHOUT (alice)      | 4 of 4 ✔
    [88/390ce7] COLLECT_GREETINGS  | 1 of 1 ✔

    Outputs:

      /workspaces/training/side-quests/types_and_records/results

      greetings:
        - {id: alice, greeting_file: greetings/alice.txt}
        - {id: dieter, greeting_file: greetings/dieter.txt}
        - {id: chiara, greeting_file: greetings/chiara.txt}
        - {id: bruno, greeting_file: greetings/bruno.txt}

      shouted:
        - {id: dieter, shouted: shouted/dieter-shouted.txt}
        - {id: chiara, shouted: shouted/chiara-shouted.txt}
        - {id: bruno, shouted: shouted/bruno-shouted.txt}
        - {id: alice, shouted: shouted/alice-shouted.txt}

      collected: COLLECTED-batch.txt
    ```

The published outputs are now listed as records, with each file field next to the `id` it belongs to.
Check the greetings:

```bash
cat results/greetings/alice.txt results/greetings/bruno.txt
```

??? abstract "File contents"

    ```console
    Hello, Alice!
    Bonjour, Bruno!
    ```

The pipeline produces the same files as before, with tidier names.

### 3.6. Rerun the collaborator's samplesheet

Now for the real test.
Run the pipeline on `data/people_v2.csv` again, the samplesheet that silently produced `null, alice!` in section 1.2:

```bash
nextflow run main.nf --input data/people_v2.csv
```

??? failure "Command output"

    ```console
     N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [astonishing_tesla] revision: 80bd6b48ff

    WARN: Static typing is a preview feature -- syntax and behavior may change in future releases

    [-        ] SAY_HELLO         -
    [-        ] SHOUT             -
    [-        ] COLLECT_GREETINGS -
    ERROR ~ Error executing process > 'SAY_HELLO (1)'

    Caused by:
      [SAY_HELLO (1)] input at index 0 cannot be null -- append `?` to the type annotation to mark it as nullable

    Tip: you can replicate the issue by changing to the process work dir and entering the command `bash .command.run`

     -- Check '.nextflow.log' file for details
    ```

This time the pipeline stops at the first `SAY_HELLO` task, before any greeting is written.
The error says that the task's input (the record at index 0) cannot be `null`.
`SAY_HELLO` destructures `greeting: String` from that record, the record has no greeting, and a typed input refuses `null` unless its type is marked as nullable (for example `String?`).
Instead of a successful run with wrong results, you get a failure that names the process and tells you which input is missing a value, before any bad output is published.

To make the pipeline work with the new samplesheet, you would map the `salutation` column to the `greeting` field when building the records.
The important change is that the mismatch can no longer slip through unnoticed.

<!-- TODO(26.10): confirm Channel<T> params ship in 26.10 -->

!!! tip

    From Nextflow 26.10, a parameter can be declared as a channel of records (for example `input: Channel<Person>`, using the record types you'll meet in the next section).
    Nextflow then loads the samplesheet itself and checks every row against the record type, so a missing column is reported when the samplesheet is loaded, naming the missing field.

### Takeaway

In this section, you've learned:

- **How to create records**: `record(name: value, ...)` builds a record, and you read fields with dot notation
- **How to merge records**: `+` combines two records, and the right-hand side wins when both have the same field
- **How processes consume records**: a destructured input such as `record(id: String, greeting_file: Path)` takes just the fields a process needs, stages `Path` fields, and rejects `null` values
- **How processes produce records**: `record(id: id, file: file(...))` names every output field

The data now carries its own field names, and the processes declare what they need.
But the workflow itself is still a loose sequence of channels in the entry workflow.
In the next section, you'll give that logic a name and a typed interface, so that anyone calling it can see exactly what goes in and what comes out.

---

## 4. Give the workflow a typed interface

Your collaborator has seen the tidied-up greetings and wants to call the same logic from their own pipeline, on their own samplesheet.
Right now, everything happens in the unnamed entry workflow, so they would have to read it line by line to find out what it needs: a channel of what, with which fields, and what comes out the other end?

The answer is to move the logic into a named workflow, as in the [Workflows of Workflows](../workflows_of_workflows/index.md) side quest.
This time, though, the workflow's inputs and outputs will have types, so its interface answers those questions itself, and Nextflow can check every call against it.

### 4.1. Declare a record type and a typed workflow

A typed `take:` block needs a type for each input.
The input here is a channel of people, so you first need a name for "a person".

#### 4.1.1. Add a record type and the `GREET` workflow

Add a `record Person` declaration and a new `GREET` workflow between the `params` block and the entry workflow:

=== "After"

    ```groovy title="main.nf" linenums="9" hl_lines="6-22"
    params {
        input: Path = 'data/people.csv'
        batch: String = 'batch'
    }

    record Person {
        id: String
        name: String
        language: String
        greeting: String
    }

    workflow GREET {
        take:
        people: Channel<Person>

        main:
        greetings = SAY_HELLO(people)

        emit:
        greetings
    }

    workflow {
    ```

=== "Before"

    ```groovy title="main.nf" linenums="9"
    params {
        input: Path = 'data/people.csv'
        batch: String = 'batch'
    }

    workflow {
    ```

`#!groovy record Person { ... }` declares a **record type**: a name for a set of fields and their types.
A record created with `record()` doesn't need to be declared up front.
A record type is what you use when you want to _refer_ to that shape somewhere, such as in a workflow input.

Record types describe the minimum a record must contain.
Any record with these four fields counts as a `Person`, even if it carries extra fields.

The `take:` block declares one input, `people: Channel<Person>`: a channel whose elements are `Person` records.
The `emit:` block has a single output, so it doesn't need a name; the workflow returns the channel itself.

#### 4.1.2. Call the workflow

In the entry workflow, call `GREET` in place of `SAY_HELLO`:

=== "After"

    ```groovy title="main.nf" linenums="44" hl_lines="1"
        greetings = GREET(people)
        shouted = SHOUT(greetings)
        collected = COLLECT_GREETINGS(shouted.map { r -> r.shouted }.collect(), params.batch)
    ```

=== "Before"

    ```groovy title="main.nf" linenums="44" hl_lines="1"
        greetings = SAY_HELLO(people)
        shouted = SHOUT(greetings)
        collected = COLLECT_GREETINGS(shouted.map { r -> r.shouted }.collect(), params.batch)
    ```

The call returns the single emitted channel, so the rest of the entry workflow doesn't change.

#### 4.1.3. Run the pipeline

Run the pipeline:

```bash
nextflow run main.nf -resume
```

??? success "Command output"

    ```console
     N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [hopeful_stallman] revision: 94f98ed45b

    WARN: Static typing is a preview feature -- syntax and behavior may change in future releases

    executor >  local (9)
    [c1/83d8aa] GREET:SAY_HELLO (alice) | 4 of 4 ✔
    [0b/0f1a5c] SHOUT (alice)           | 4 of 4 ✔
    [13/4c4553] COLLECT_GREETINGS       | 1 of 1 ✔

    Outputs:

      /workspaces/training/side-quests/types_and_records/results

      greetings:
        - {id: bruno, greeting_file: greetings/bruno.txt}
        - {id: chiara, greeting_file: greetings/chiara.txt}
        - {id: alice, greeting_file: greetings/alice.txt}
        - {id: dieter, greeting_file: greetings/dieter.txt}

      shouted:
        - {id: bruno, shouted: shouted/bruno-shouted.txt}
        - {id: dieter, shouted: shouted/dieter-shouted.txt}
        - {id: chiara, shouted: shouted/chiara-shouted.txt}
        - {id: alice, shouted: shouted/alice-shouted.txt}

      collected: COLLECTED-batch.txt
    ```

`SAY_HELLO` is now listed as `GREET:SAY_HELLO`, since it runs inside the `GREET` workflow.
The process name is part of each task's hash, so every `SAY_HELLO` task ran again under its new name, and the tasks downstream of it followed.

!!! tip "Record types across files"

    A record type is only visible in the file that declares it.
    If you moved `GREET` and `Person` into their own file, say `workflows/greet.nf`, then any other script that uses the `Person` type would need to include it by name, alongside the workflow:

    ```groovy
    include { GREET ; Person } from './workflows/greet.nf'
    ```

    `nextflow lint` reports `` `Person` is not defined `` if you forget.

### 4.2. Emit more than one output

`GREET` only covers the first step.
Move `SHOUT` inside it too.
The workflow then has two results that the caller needs (the greetings and the shouted greetings), so each output gets a name.

#### 4.2.1. Add named outputs

Update the `GREET` workflow:

=== "After"

    ```groovy title="main.nf" linenums="21" hl_lines="7 10 11"
    workflow GREET {
        take:
        people: Channel<Person>

        main:
        greetings = SAY_HELLO(people)
        shouted = SHOUT(greetings)

        emit:
        greetings: Channel<Record> = greetings
        shouted: Channel<Record> = shouted
    }
    ```

=== "Before"

    ```groovy title="main.nf" linenums="21" hl_lines="9"
    workflow GREET {
        take:
        people: Channel<Person>

        main:
        greetings = SAY_HELLO(people)

        emit:
        greetings
    }
    ```

Each named output has the form `name: Type = value`.
Both outputs are channels of records; `Record` is the general type for any record when you don't need a specific record type.

#### 4.2.2. Use the named outputs

When a workflow has more than one output, the call returns a record with one field per output.
Update the entry workflow to keep that result in `res` and read its fields:

=== "After"

    ```groovy title="main.nf" linenums="46" hl_lines="1 2 5 6"
        res = GREET(people)
        collected = COLLECT_GREETINGS(res.shouted.map { r -> r.shouted }.collect(), params.batch)

        publish:
        greetings = res.greetings
        shouted = res.shouted
        collected = collected
    ```

=== "Before"

    ```groovy title="main.nf" linenums="46" hl_lines="1 2 3 6 7"
        greetings = GREET(people)
        shouted = SHOUT(greetings)
        collected = COLLECT_GREETINGS(shouted.map { r -> r.shouted }.collect(), params.batch)

        publish:
        greetings = greetings
        shouted = shouted
        collected = collected
    ```

`res.greetings` and `res.shouted` are the two channels emitted by `GREET`.
In typed code you always access a workflow's outputs this way, through the value the call returns.

#### 4.2.3. Run the pipeline

Run the pipeline:

```bash
nextflow run main.nf -resume
```

??? success "Command output"

    ```console
     N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [amazing_morse] revision: 0574ac6ec1

    WARN: Static typing is a preview feature -- syntax and behavior may change in future releases

    executor >  local (5)
    [9d/97dab2] GREET:SAY_HELLO (chiara) | 4 of 4, cached: 4 ✔
    [e0/7a4e9b] GREET:SHOUT (chiara)     | 4 of 4 ✔
    [87/dc4c70] COLLECT_GREETINGS        | 1 of 1 ✔

    Outputs:

      /workspaces/training/side-quests/types_and_records/results

      greetings:
        - {id: alice, greeting_file: greetings/alice.txt}
        - {id: bruno, greeting_file: greetings/bruno.txt}
        - {id: chiara, greeting_file: greetings/chiara.txt}
        - {id: dieter, greeting_file: greetings/dieter.txt}

      shouted:
        - {id: alice, shouted: shouted/alice-shouted.txt}
        - {id: bruno, shouted: shouted/bruno-shouted.txt}
        - {id: dieter, shouted: shouted/dieter-shouted.txt}
        - {id: chiara, shouted: shouted/chiara-shouted.txt}

      collected: COLLECTED-batch.txt
    ```

Both processes now run inside `GREET`, and all three outputs are still published.
`GREET:SAY_HELLO` came from the cache, while `SHOUT` ran again under its new name, `GREET:SHOUT`, and `COLLECT_GREETINGS` followed.

### 4.3. Pass and return single values

The last step is to move `COLLECT_GREETINGS` into `GREET`, so that the whole greeting logic lives behind one interface.
That brings up a question the workflow hasn't had to answer yet: what type is the batch name, and what type is the collected file?

#### 4.3.1. Channels and values

Nextflow has two dataflow types:

- A **`Channel`** carries any number of items, which arrive one at a time. `people` is a channel: one record per person.
- A **`Value`** carries exactly one item, which may not be available yet. A `Value` can be used any number of times, and every process task that needs it gets the same item.

You've used both already without naming them.
`collect()` turns a channel into a value: one item containing everything the channel emitted.
And a process that is called only with values runs once and returns a value, which is why `COLLECT_GREETINGS` produces a single file.

#### 4.3.2. Move `COLLECT_GREETINGS` into the workflow

Add a `batch` input of type `Value<String>` and a `collected` output of type `Value<Path>`:

=== "After"

    ```groovy title="main.nf" linenums="21" hl_lines="4 9 14"
    workflow GREET {
        take:
        people: Channel<Person>
        batch: Value<String>

        main:
        greetings = SAY_HELLO(people)
        shouted = SHOUT(greetings)
        collected = COLLECT_GREETINGS(shouted.map { r -> r.shouted }.collect(), batch)

        emit:
        greetings: Channel<Record> = greetings
        shouted: Channel<Record> = shouted
        collected: Value<Path> = collected
    }
    ```

=== "Before"

    ```groovy title="main.nf" linenums="21"
    workflow GREET {
        take:
        people: Channel<Person>

        main:
        greetings = SAY_HELLO(people)
        shouted = SHOUT(greetings)

        emit:
        greetings: Channel<Record> = greetings
        shouted: Channel<Record> = shouted
    }
    ```

Then update the entry workflow:

=== "After"

    ```groovy title="main.nf" linenums="49" hl_lines="1 6"
        res = GREET(people, channel.value(params.batch))

        publish:
        greetings = res.greetings
        shouted = res.shouted
        collected = res.collected
    ```

=== "Before"

    ```groovy title="main.nf" linenums="49" hl_lines="1 2 7"
        res = GREET(people)
        collected = COLLECT_GREETINGS(res.shouted.map { r -> r.shouted }.collect(), params.batch)

        publish:
        greetings = res.greetings
        shouted = res.shouted
        collected = collected
    ```

`params.batch` is a plain `String`, not a dataflow value, so `channel.value(params.batch)` wraps it in a `Value<String>` to match the declared input.
Reading the `GREET` interface now tells you everything a caller needs to know: it takes a channel of people and one batch name, and it returns two channels of records and one file.

#### 4.3.3. Run the pipeline

Run the pipeline:

```bash
nextflow run main.nf -resume
```

??? success "Command output"

    ```console
     N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [sad_brattain] revision: 390a881f6d

    WARN: Static typing is a preview feature -- syntax and behavior may change in future releases

    executor >  local (1)
    [74/f8a757] GREET:SAY_HELLO (dieter) | 4 of 4, cached: 4 ✔
    [42/6079d8] GREET:SHOUT (bruno)      | 4 of 4, cached: 4 ✔
    [46/824db1] GREET:COLLECT_GREETINGS  | 1 of 1 ✔

    Outputs:

      /workspaces/training/side-quests/types_and_records/results

      greetings:
        - {id: dieter, greeting_file: greetings/dieter.txt}
        - {id: chiara, greeting_file: greetings/chiara.txt}
        - {id: alice, greeting_file: greetings/alice.txt}
        - {id: bruno, greeting_file: greetings/bruno.txt}

      shouted:
        - {id: dieter, shouted: shouted/dieter-shouted.txt}
        - {id: bruno, shouted: shouted/bruno-shouted.txt}
        - {id: chiara, shouted: shouted/chiara-shouted.txt}
        - {id: alice, shouted: shouted/alice-shouted.txt}

      collected: COLLECTED-batch.txt
    ```

All three processes now run inside `GREET`, and only `GREET:COLLECT_GREETINGS` had to run again.
The collected file is published as before:

```bash
cat results/COLLECTED-batch.txt
```

??? abstract "File contents"

    ```console
    CIAO, CHIARA!
    BONJOUR, BRUNO!
    HELLO, ALICE!
    HALLO, DIETER!
    ```

The order of the lines may differ in your run, since the greetings are collected in whatever order the `SHOUT` tasks finish.

### Takeaway

In this section, you've learned:

- **How to declare a record type**: `#!groovy record Person { ... }` names a shape, so you can refer to it in interfaces
- **How to type a workflow's inputs**: `take:` declares `Channel<T>` and `Value<T>` inputs
- **How to emit outputs**: a single output is unnamed and returned directly, while multiple outputs are named and returned as a record (`res.greetings`)
- **The difference between `Channel` and `Value`**: many items over time versus exactly one item

The pipeline is fully typed now: every process and workflow declares what it takes and what it returns.
In the last section, you'll see how that information lets `nextflow lint` catch mistakes before they ever reach a run.

---

## 5. Let the type checker help

In section 1, a renamed column went unnoticed until someone read the output.
With types in place, mistakes like that are caught much earlier, and some of them before the pipeline runs at all.
In this section you'll make two mistakes on purpose and see where each one is caught.

### 5.1. A missing argument

Suppose you forget to pass the batch name when calling `GREET`:

=== "After"

    ```groovy title="main.nf" linenums="49" hl_lines="1"
        res = GREET(people)
    ```

=== "Before"

    ```groovy title="main.nf" linenums="49" hl_lines="1"
        res = GREET(people, channel.value(params.batch))
    ```

Check the script with `nextflow lint`:

```bash
nextflow lint main.nf
```

??? failure "Command output"

    ```console
    Linting Nextflow code..
    Linting: main.nf
    Error main.nf:49:11: Incorrect number of call arguments, expected 2 but received 1
    │  49 |     res = GREET(people)
    ╰     |           ^^^^^^^^^^^^^


    Nextflow linting complete!
     ❌ 1 file had 1 error
     ✅ 3 files had no errors
    ```

The linter knows that `GREET` declares two inputs, so it reports the problem without running anything.
Put the argument back:

=== "After"

    ```groovy title="main.nf" linenums="49" hl_lines="1"
        res = GREET(people, channel.value(params.batch))
    ```

=== "Before"

    ```groovy title="main.nf" linenums="49" hl_lines="1"
        res = GREET(people)
    ```

### 5.2. A misspelled field

Now misspell a field name where a record is used in the workflow, in the `map` that feeds `COLLECT_GREETINGS`:

=== "After"

    ```groovy title="main.nf" linenums="29" hl_lines="1"
        collected = COLLECT_GREETINGS(shouted.map { r -> r.shout }.collect(), batch)
    ```

=== "Before"

    ```groovy title="main.nf" linenums="29" hl_lines="1"
        collected = COLLECT_GREETINGS(shouted.map { r -> r.shouted }.collect(), batch)
    ```

This is the same kind of mistake as the renamed column in section 1.2: code asks for a field that the data doesn't have.
Check the script:

```bash
nextflow lint main.nf
```

??? success "Command output"

    ```console
    Linting Nextflow code..
    Linting: main.nf
    Nextflow linting complete!
     ✅ 4 files had no errors
    ```

The linter in Nextflow 26.04 doesn't catch this one.
Run the pipeline:

```bash
nextflow run main.nf -resume
```

??? failure "Command output"

    ```console
     N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [curious_avogadro] revision: 32cae9123e

    WARN: Static typing is a preview feature -- syntax and behavior may change in future releases

    [c1/83d8aa] GREET:SAY_HELLO (alice) | 4 of 4, cached: 4 ✔
    [e0/7a4e9b] GREET:SHOUT (chiara)    | 4 of 4, cached: 4 ✔
    [-        ] GREET:COLLECT_GREETINGS -
    ERROR ~ Error executing process > 'GREET:COLLECT_GREETINGS'

    Caused by:
      Path value cannot be null

    Tip: you can replicate the issue by changing to the process work dir and entering the command `bash .command.run`

     -- Check '.nextflow.log' file for details
    ```

`r.shout` evaluates to `null`, and `COLLECT_GREETINGS` rejects the `null` path because its input is typed as `Bag<Path>`.
The mistake is still caught, at a process boundary instead of in the published results, but only after the upstream tasks have run.

<!-- TODO(26.10): recapture this lint output on the 26.10 release -->

Newer versions of Nextflow have a more complete type checker.
They know that `SHOUT` emits records with the fields `id` and `shouted`, so `nextflow lint` reports this mistake before anything runs:

```console title="nextflow lint main.nf (newer Nextflow versions)"
Linting Nextflow code..
Linting: main.nf
Error main.nf:29:54: Unrecognized property `shout` for type Record {
    id: String
    shouted: Path
}
│  29 |     collected = COLLECT_GREETINGS(shouted.map { r -> r.shout }.collect
╰     |                                                      ^^^^^^^


Nextflow linting complete!
 ❌ 1 file had 1 error
 ✅ 3 files had no errors
```

!!! warning "Missing fields and `null` on Nextflow 26.04"

    On Nextflow 26.04, reading a field that a record doesn't have, such as a misspelled field name, gives `null`.
    Destructured process inputs and plain typed inputs reject that `null` when a task starts, as you saw in sections 3.6 and 5.2.
    A process input declared with a record type (for example `person: Person`) does not check its fields on this version, so a missing field can still reach the script as `null`.
    Newer versions of Nextflow fail fast in that case too, and their linter catches most misspelled fields before you run anything.

Fix the typo, then run the linter one last time over the whole project:

=== "After"

    ```groovy title="main.nf" linenums="29" hl_lines="1"
        collected = COLLECT_GREETINGS(shouted.map { r -> r.shouted }.collect(), batch)
    ```

=== "Before"

    ```groovy title="main.nf" linenums="29" hl_lines="1"
        collected = COLLECT_GREETINGS(shouted.map { r -> r.shout }.collect(), batch)
    ```

```bash
nextflow lint .
```

??? success "Command output"

    ```console
    Linting Nextflow code..
    Linting: modules/shout.nf
    Linting: modules/say_hello.nf
    Linting: modules/collect_greetings.nf
    Linting: nextflow.config
    Linting: main.nf
    Nextflow linting complete!
     ✅ 5 files had no errors
    ```

### Takeaway

In this section, you've seen where different mistakes are caught:

- **A call that doesn't match a workflow's interface**: caught by `nextflow lint`, before running
- **A misspelled record field in workflow code**: caught when a typed process input receives `null`, and caught by `nextflow lint` on newer versions
- **Nothing reaches the published results silently**, which was the problem you started with

Run `nextflow lint` after every change: the more of your pipeline is typed, the more it can check.

---

## Summary

In this side quest, you've migrated an untyped pipeline to static types and records.
You started with a pipeline that happily wrote `null, alice!` when a samplesheet column was renamed, and ended with one that refuses bad input at the first process, with interfaces that say exactly what goes in and what comes out.

Compared with meta maps and tuples, typed code with records:

- Names every field, so producers and consumers never depend on element order
- Lets each process declare exactly which fields it needs, and ignore the rest
- Turns missing values into errors instead of `null` in your results
- Gives `nextflow lint` enough information to catch mistakes before you run anything

### Key patterns

1.  **Enabling static typing**: Add the feature flag to each script that uses typed processes or workflows.

    ```groovy
    nextflow.enable.types = true
    ```

2.  **Typed process inputs and outputs**: Inputs are `name: Type`, and outputs are values.

    ```groovy
    input:
    input_files: Bag<Path>
    batch_name: String

    output:
    file("COLLECTED-${batch_name}.txt")
    ```

3.  **Reading a samplesheet into records**: Split the file inside `flatMap`, then build one record per row.

    ```groovy
    channel.of(params.input)
        .flatMap { csv -> csv.splitCsv(header: true) }
        .map { row -> record(id: row.id, name: row.name) }
    ```

4.  **Merging records**: The right-hand side wins when both records have the same field.

    ```groovy
    person + record(name: person.name.capitalize())
    ```

5.  **Record inputs and outputs**: Destructure only the fields a process needs, and name every output field.

    ```groovy
    input:
    record(id: String, greeting_file: Path)

    output:
    record(id: id, shouted: file("${id}-shouted.txt"))
    ```

6.  **Record types and typed workflows**: Name a shape with `record`, and use `Channel<T>` and `Value<T>` in `take:` and `emit:`.

    ```groovy
    record Person {
        id: String
        name: String
    }

    workflow GREET {
        take:
        people: Channel<Person>
        batch: Value<String>

        main:
        // ...

        emit:
        greetings: Channel<Record> = greetings
        collected: Value<Path> = collected
    }
    ```

    The caller reads multiple outputs from the returned record: `res = GREET(people, channel.value(params.batch))`, then `res.greetings`.

### Additional resources

- [Static typing](https://www.nextflow.io/docs/latest/static-typing.html)
- [Typed processes](https://www.nextflow.io/docs/latest/process-typed.html)
- [Typed workflows](https://www.nextflow.io/docs/latest/workflow-typed.html)
- [Record type reference](https://www.nextflow.io/docs/latest/reference/stdlib-types.html#record)
- [Migrating to static typing](https://www.nextflow.io/docs/latest/tutorials/static-types.html)
- [`nextflow lint` command reference](https://www.nextflow.io/docs/latest/reference/cli.html#lint)

---

## What's next?

Return to the [menu of Side Quests](../index.md) or click the button in the bottom right of the page to move on to the next topic in the list.
