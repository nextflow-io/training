<!-- TODO(26.10): recapture console output on 26.10 release -->

# Types and Records

When you're developing a pipeline, most of the data flowing through your channels has no declared shape.
A meta map is a bag of keys, a tuple is a list of values in some order, and nothing tells Nextflow which keys or positions a process expects.
That works until something changes: a collaborator renames a column in a samplesheet, or a module reads `meta.sample` when the map contains `meta.id`.
The pipeline usually keeps running, quietly substitutes `null`, and you only find out when someone reads the results.

Nextflow's **static typing** lets you declare what your data looks like, and **records** give your metadata named, typed fields in place of maps and positional tuples.
Together they turn many of these silent mistakes into errors that name the problem, often before the pipeline runs at all.

### Learning goals

In this side quest, you'll take a small, working, untyped pipeline, watch it produce wrong results without complaint, and migrate it step by step to static types and records.

By the end of this side quest, you'll be able to:

- Enable static typing and convert processes to typed inputs and outputs
- Declare record types and load a samplesheet directly into records
- Carry metadata through a process with `+`, and select a subset of fields with a destructured input
- Give a named workflow typed `take:` and `emit:` blocks, and explain the difference between a `Channel` and a `Value`
- Use `nextflow lint` to catch typing mistakes before you run the pipeline

### Prerequisites

Before taking on this side quest, you should:

- Have completed the [Hello Nextflow](../../hello_nextflow/index.md) tutorial or equivalent beginner's course.
- Be comfortable using basic Nextflow concepts and mechanisms (processes, channels, operators, modules)
- Be familiar with meta maps, as covered in the [Metadata and Meta Maps](../metadata/index.md) side quest

!!! note "Nextflow version"

    This side quest requires Nextflow 26.10 or later.

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

The samplesheet lists a few people and how to greet each of them:

```csv title="data/people.csv"
id,name,greeting
alice,Alice,Hello
bruno,Bruno,Bonjour
chiara,Chiara,Ciao
dieter,Dieter,Hallo
```

The pipeline reads each row into a meta map, writes a greeting file for each person, makes an uppercase copy of each greeting, and collects all the uppercase greetings into a single file.
The processes are deliberately tiny (`echo`, `tr` and `cat`), so you can focus on how the data flows between them.

#### Review the assignment

Your challenge is to migrate this pipeline to static types and records without changing the files it produces, and to make it fail loudly when its input doesn't have the shape it expects.

Along the way, you will:

1. **Watch the untyped pipeline fail silently** on a samplesheet with a renamed column
2. **Enable static typing** and convert a process to typed inputs and outputs
3. **Replace meta maps and tuples with records**, so the same bad samplesheet is rejected immediately
4. **Wrap the processing in a typed workflow** with a declared interface
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
            [id: row.id, name: row.name, greeting: row.greeting]
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
     N E X T F L O W   ~  version 26.09.1-edge

    Launching `main.nf` [high_becquerel] revision: 6d1746eb92

    executor >  local (9)
    [20/8b3387] SAY_HELLO (dieter) | 4 of 4 ✔
    [0d/e1cc86] SHOUT (bruno)      | 4 of 4 ✔
    [46/69d139] COLLECT_GREETINGS  | 1 of 1 ✔

    Outputs:

      /workspaces/training/side-quests/types_and_records/results

      greetings:
        - [{id: bruno, name: Bruno, greeting: Bonjour}, greetings/bruno.txt]
        - [{id: alice, name: Alice, greeting: Hello}, greetings/alice.txt]
        - [{id: dieter, name: Dieter, greeting: Hallo}, greetings/dieter.txt]
        - [{id: chiara, name: Chiara, greeting: Ciao}, greetings/chiara.txt]

      shouted:
        - [{id: dieter, name: Dieter, greeting: Hallo}, shouted/dieter-shouted.txt]
        - [{id: alice, name: Alice, greeting: Hello}, shouted/alice-shouted.txt]
        - [{id: chiara, name: Chiara, greeting: Ciao}, shouted/chiara-shouted.txt]
        - [{id: bruno, name: Bruno, greeting: Bonjour}, shouted/bruno-shouted.txt]

      collected: COLLECTED-batch.txt
    ```

Check one of the greetings it wrote:

```bash
cat results/greetings/alice.txt
```

??? abstract "File contents"

    ```console
    Hello, Alice!
    ```

### 1.2. Run it on a collaborator's samplesheet

A collaborator has exported the same people from a different system and sent you `data/people_v2.csv`.
Take a look at it:

```csv title="data/people_v2.csv"
id,name,salutation
alice,Alice,Hello
bruno,Bruno,Bonjour
chiara,Chiara,Ciao
dieter,Dieter,Hallo
```

The other system calls the greeting column `salutation`.
It's an easy difference to miss, so run the pipeline on the new file and see what happens:

```bash
nextflow run main.nf --input data/people_v2.csv -resume
```

??? success "Command output"

    ```console
     N E X T F L O W   ~  version 26.09.1-edge

    Launching `main.nf` [reverent_salas] revision: 6d1746eb92

    executor >  local (9)
    [c1/30afc6] SAY_HELLO (chiara) | 4 of 4 ✔
    [f3/eeeb9c] SHOUT (chiara)     | 4 of 4 ✔
    [e9/375412] COLLECT_GREETINGS  | 1 of 1 ✔

    Outputs:

      /workspaces/training/side-quests/types_and_records/results

      greetings:
        - [{id: chiara, name: Chiara, greeting: null}, greetings/chiara.txt]
        - [{id: dieter, name: Dieter, greeting: null}, greetings/dieter.txt]
        - [{id: bruno, name: Bruno, greeting: null}, greetings/bruno.txt]
        - [{id: alice, name: Alice, greeting: null}, greetings/alice.txt]

      shouted:
        - [{id: bruno, name: Bruno, greeting: null}, shouted/bruno-shouted.txt]
        - [{id: dieter, name: Dieter, greeting: null}, shouted/dieter-shouted.txt]
        - [{id: alice, name: Alice, greeting: null}, shouted/alice-shouted.txt]
        - [{id: chiara, name: Chiara, greeting: null}, shouted/chiara-shouted.txt]

      collected: COLLECTED-batch.txt
    ```

Every task succeeded and every output was published.
Now look at what was actually written:

```bash
cat results/greetings/alice.txt
```

??? abstract "File contents"

    ```console
    null, Alice!
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
In the next section, you'll switch it on for the first file.

---

## 2. Enable static typing

Static typing in Nextflow is opt-in through a feature flag, set in each script.
That means you can migrate a pipeline gradually, one file at a time.
In this section you'll start with the simplest module, `COLLECT_GREETINGS`, which takes a list of files and a batch name and produces one file.

### 2.1. Turn on static typing in a module

Add the feature flag at the top of `modules/collect_greetings.nf`:

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
     N E X T F L O W   ~  version 26.09.1-edge

    Launching `main.nf` [elated_almeida] revision: 6d1746eb92

    Error modules/collect_greetings.nf:9:5: Invalid input declaration in typed process
    │   9 |     path input_files
    ╰     |     ^^^^^^^^^^^^^^^^

    ERROR ~ Script compilation failed

     -- Check '.nextflow.log' file for details
    ```

The script doesn't compile.
Once a script enables static typing, every process in it must use the typed syntax, and `path input_files` is a legacy input declaration.
Typed and legacy processes can live in the same pipeline, but not in the same file.

### 2.2. Rewrite the inputs and outputs

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
- `path input_files` becomes `input_files: Bag<Path>`.
  `Path` inputs are staged into the task directory like `path` inputs, and `Bag<Path>` says the input is an unordered collection of files, which is what `collect()` produces.
- The output is now a plain value: `#!groovy file("COLLECTED-${batch_name}.txt")` returns the file of that name from the task directory.
  A process with a single output doesn't need to name it.
- `input_files` is a collection, and `#!groovy ${input_files}` would render it in list form, as `[alice-shouted.txt, bruno-shouted.txt, ...]`, brackets and commas included.
  `.join(' ')` turns it into the space-separated file names that `cat` expects.

Run the pipeline to check that the converted process still produces the same result:

```bash
nextflow run main.nf -resume
```

??? success "Command output"

    ```console
     N E X T F L O W   ~  version 26.09.1-edge

    Launching `main.nf` [compassionate_hoover] revision: 6d1746eb92

    executor >  local (1)
    [dc/0c7e7d] SAY_HELLO (alice) | 4 of 4, cached: 4 ✔
    [5b/ed0d3a] SHOUT (alice)     | 4 of 4, cached: 4 ✔
    [72/6ff8b9] COLLECT_GREETINGS | 1 of 1 ✔

    Outputs:

      /workspaces/training/side-quests/types_and_records/results

      greetings:
        - [{id: alice, name: Alice, greeting: Hello}, greetings/alice.txt]
        - [{id: chiara, name: Chiara, greeting: Ciao}, greetings/chiara.txt]
        - [{id: dieter, name: Dieter, greeting: Hallo}, greetings/dieter.txt]
        - [{id: bruno, name: Bruno, greeting: Bonjour}, greetings/bruno.txt]

      shouted:
        - [{id: dieter, name: Dieter, greeting: Hallo}, shouted/dieter-shouted.txt]
        - [{id: alice, name: Alice, greeting: Hello}, shouted/alice-shouted.txt]
        - [{id: chiara, name: Chiara, greeting: Ciao}, shouted/chiara-shouted.txt]
        - [{id: bruno, name: Bruno, greeting: Bonjour}, shouted/bruno-shouted.txt]

      collected: COLLECTED-batch.txt
    ```

`COLLECT_GREETINGS` ran again because its definition changed, while `SAY_HELLO` and `SHOUT` came from the cache.
The run also put the correct greetings back in `results/`, since `--input` is no longer set to the collaborator's file.

`main.nf` and the other two modules are still untyped, and they call the typed `COLLECT_GREETINGS` without any change.
That's what makes a gradual migration possible.

### Takeaway

In this section, you've learned:

- **How to enable static typing**: `nextflow.enable.types = true`, in each script that uses typed code
- **How to migrate incrementally**: typed and untyped files work together, but each file is either typed or untyped
- **How to write a typed process**: inputs are `name: Type`, and outputs are values such as `file(...)`

You now have one typed process.
But the process that wrote the `null` greeting, `SAY_HELLO`, still takes an untyped meta map.
Typing it as `meta: Map` wouldn't help much, because a `Map` can still hold any keys or none.
What you need is a data structure whose fields Nextflow knows about, and that's what records are.

---

## 3. Replace meta maps and tuples with records

A **record** is a set of named fields, like a meta map, but designed for static typing.
Records can hold files and values side by side, so they replace both the meta map and the `[meta, file]` tuple.

Section 1 showed two problems: a missing column silently became `null`, and tuples depend on the position of each element.
In this section you'll fix the first one by loading the samplesheet straight into records of a declared type, then fix the second by switching the processes over to records.

### 3.1. Declare a record type

A **record type** gives a name to a set of fields and their types.
Create a new file, `types.nf`, next to `main.nf`:

```groovy title="types.nf" linenums="1"
nextflow.enable.types = true

record Person {
    id: String
    name: String
    greeting: String
}
```

`Person` says that a person has an `id`, a `name` and a `greeting`, all strings.
The record types live in their own file so that `main.nf` and the modules can all use them.
A record type is only visible in the file that declares it, so every script that uses `Person` includes it by name, as you'll see next.

### 3.2. Load the samplesheet into records

Now update the top of `main.nf`.
Turn on static typing, include `Person`, declare `input` as a channel of people, and pass it straight to `SAY_HELLO`:

=== "After"

    ```groovy title="main.nf" linenums="1" hl_lines="3 8 11 17"
    #!/usr/bin/env nextflow

    nextflow.enable.types = true

    include { SAY_HELLO } from './modules/say_hello.nf'
    include { SHOUT } from './modules/shout.nf'
    include { COLLECT_GREETINGS } from './modules/collect_greetings.nf'
    include { Person } from './types.nf'

    params {
        input: Channel<Person>
        batch: String = 'batch'
    }

    workflow {
        main:
        greetings = SAY_HELLO(params.input)
    ```

=== "Before"

    ```groovy title="main.nf" linenums="1" hl_lines="8 14-21"
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
                [id: row.id, name: row.name, greeting: row.greeting]
            }

        greetings = SAY_HELLO(people)
    ```

The whole block that read the CSV is gone.
When a parameter is declared as `Channel<Person>`, Nextflow loads the samplesheet you pass it (CSV, JSON or YAML) and emits one `Person` record per row.
It checks each row against the record type as it goes, and converts each value to the declared field type.

A default in the `params` block must already have the parameter's declared type, and a file name is a `String`, not a `Channel<Person>`.
Nextflow only loads a samplesheet into records when the value comes from the command line or a config file.
So the default samplesheet goes in `nextflow.config`:

=== "After"

    ```groovy title="nextflow.config" linenums="1" hl_lines="1"
    params.input = 'data/people.csv'

    outputDir = 'results'
    workflow.output.mode = 'copy'
    ```

=== "Before"

    ```groovy title="nextflow.config" linenums="1"
    outputDir = 'results'
    workflow.output.mode = 'copy'
    ```

Run the pipeline on the original samplesheet:

```bash
nextflow run main.nf -resume
```

??? success "Command output"

    ```console
     N E X T F L O W   ~  version 26.09.1-edge

    Launching `main.nf` [dreamy_shockley] revision: 7564fcec8f

    WARN: Type checking found 1 error(s) -- run `nextflow lint` to inspect them

    [40/395cbd] SAY_HELLO (chiara) | 4 of 4, cached: 4 ✔
    [5b/ed0d3a] SHOUT (alice)      | 4 of 4, cached: 4 ✔
    [72/6ff8b9] COLLECT_GREETINGS  | 1 of 1, cached: 1 ✔

    Outputs:

      /workspaces/training/side-quests/types_and_records/results

      greetings:
        - [{id: chiara, name: Chiara, greeting: Ciao}, greetings/chiara.txt]
        - [{id: bruno, name: Bruno, greeting: Bonjour}, greetings/bruno.txt]
        - [{id: alice, name: Alice, greeting: Hello}, greetings/alice.txt]
        - [{id: dieter, name: Dieter, greeting: Hallo}, greetings/dieter.txt]

      shouted:
        - [{id: dieter, name: Dieter, greeting: Hallo}, shouted/dieter-shouted.txt]
        - [{id: alice, name: Alice, greeting: Hello}, shouted/alice-shouted.txt]
        - [{id: bruno, name: Bruno, greeting: Bonjour}, shouted/bruno-shouted.txt]
        - [{id: chiara, name: Chiara, greeting: Ciao}, shouted/chiara-shouted.txt]

      collected: COLLECTED-batch.txt
    ```

The untyped `SAY_HELLO` and `SHOUT` accept the `Person` records without any change: `meta.greeting` and `meta.name` read a record's fields the same way they read a map's keys.
Every task came from the cache, because each `Person` record holds the same fields and values as the meta map it replaces.
The published outputs still show `[meta, file]` tuples, now with a record in place of the map.

The warning comes from the type checker, which Nextflow runs before every pipeline.
`main.nf` is now typed, but `SHOUT` isn't, so the checker doesn't know that `SHOUT` returns a channel of `[meta, file]` tuples.
It treats the result as a generic `Record`, and reports the `#!groovy shouted.map { ... }` call as a method a `Record` doesn't have.
The pipeline still runs correctly.
This is normal partway through a migration, and the warning goes away once `SHOUT` is typed in section 3.5.

??? info "The error behind the warning"

    `nextflow lint` shows the error that the warning counts:

    ```console
    Error main.nf:19:43: Unrecognized method `map` for type Record
    │  19 |     collected = COLLECT_GREETINGS(shouted.map { _meta, file -> file }.
    ╰     |                                           ^^^
    ```

### 3.3. Rerun the collaborator's samplesheet

Now for the real test.
Run the pipeline on `data/people_v2.csv` again, the samplesheet that silently produced `null, Alice!` in section 1.2:

```bash
nextflow run main.nf --input data/people_v2.csv -resume
```

??? failure "Command output"

    ```console
     N E X T F L O W   ~  version 26.09.1-edge

    Launching `main.nf` [sharp_booth] revision: 7564fcec8f

    WARN: Type checking found 1 error(s) -- run `nextflow lint` to inspect them

    Invalid record in samplesheet 'data/people_v2.csv' for parameter `input` -- Input record [id:alice, name:Alice, salutation:Hello] is missing field 'greeting' required by record type 'Person'
    ```

This time no task runs at all.
Nextflow checks every row against `Person` while it loads the samplesheet, and the first row has no `greeting`.
The error names the samplesheet, the row, the missing field and the record type that requires it.

The record type is what catches this.
`greeting: String` is not nullable, so a `Person` without a greeting is not a valid `Person`.
A field that may legitimately be missing would be declared with a `?`, as `greeting: String?`.

### 3.4. Fix the collaborator's samplesheet

The error tells you exactly what to change: the samplesheet needs a `greeting` column.
Rename the `salutation` column in `data/people_v2.csv`:

=== "After"

    ```csv title="data/people_v2.csv" linenums="1" hl_lines="1"
    id,name,greeting
    ```

=== "Before"

    ```csv title="data/people_v2.csv" linenums="1" hl_lines="1"
    id,name,salutation
    ```

Run it again:

```bash
nextflow run main.nf --input data/people_v2.csv -resume
```

??? success "Command output"

    ```console
     N E X T F L O W   ~  version 26.09.1-edge

    Launching `main.nf` [peaceful_venter] revision: 7564fcec8f

    WARN: Type checking found 1 error(s) -- run `nextflow lint` to inspect them

    [40/395cbd] SAY_HELLO (chiara) | 4 of 4, cached: 4 ✔
    [5b/ed0d3a] SHOUT (alice)      | 4 of 4, cached: 4 ✔
    [72/6ff8b9] COLLECT_GREETINGS  | 1 of 1, cached: 1 ✔

    Outputs:

      /workspaces/training/side-quests/types_and_records/results

      greetings:
        - [{id: bruno, name: Bruno, greeting: Bonjour}, greetings/bruno.txt]
        - [{id: alice, name: Alice, greeting: Hello}, greetings/alice.txt]
        - [{id: dieter, name: Dieter, greeting: Hallo}, greetings/dieter.txt]
        - [{id: chiara, name: Chiara, greeting: Ciao}, greetings/chiara.txt]

      shouted:
        - [{id: bruno, name: Bruno, greeting: Bonjour}, shouted/bruno-shouted.txt]
        - [{id: alice, name: Alice, greeting: Hello}, shouted/alice-shouted.txt]
        - [{id: chiara, name: Chiara, greeting: Ciao}, shouted/chiara-shouted.txt]
        - [{id: dieter, name: Dieter, greeting: Hallo}, shouted/dieter-shouted.txt]

      collected: COLLECTED-batch.txt
    ```

Every task came from the cache.
The fixed samplesheet produces exactly the same `Person` records as `data/people.csv`, so every task has the same inputs as in the previous run.

### 3.5. Replace the tuples

The `null` problem is solved, but the processes still pass `[meta, file]` tuples, so every producer and consumer must agree on the order of the elements.
Records fix that too, by naming every field.

#### 3.5.1. Carry metadata through `SAY_HELLO` with `+`

In the metadata side quest, a process kept its metadata by emitting `[meta, file]`, or by adding to the map with `meta + [...]`.
Records support the same idea: `+` merges two records into a new one.

Replace the contents of `modules/say_hello.nf`:

=== "After"

    ```groovy title="modules/say_hello.nf" linenums="1" hl_lines="1 3 9 12 15 19"
    nextflow.enable.types = true

    include { Person } from '../types.nf'

    /*
     * Write a personalised greeting to a file
     */
    process SAY_HELLO {
        tag "${person.id}"

        input:
        person: Person

        output:
        person + record(greeting_file: file("${person.id}.txt"))

        script:
        """
        echo '${person.greeting}, ${person.name}!' > ${person.id}.txt
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

The input is now `person: Person`, a record of the type you declared, and the script reads its fields by name (`person.greeting`, `person.name`).
The module includes `Person` from `types.nf`, because the type isn't visible here otherwise; without the `include`, `nextflow lint` reports `` `Person` is not defined ``.

The output is `#!groovy person + record(greeting_file: file("${person.id}.txt"))`: the incoming person with one new field added.
This is the record version of `meta + [...]`.
The `+` operator builds a new record:

- Fields that appear on only one side are copied into the result.
- Fields that appear on both sides take the value from the **right-hand** record.

So each output record has the person's `id`, `name` and `greeting`, plus the `greeting_file`.
Downstream code asks for `greeting_file` by name instead of remembering that the file is the second element of a tuple.

#### 3.5.2. Select only the fields `SHOUT` needs

`SHOUT` needs much less: the person's `id` and the greeting file.
Instead of taking a whole record, it can declare exactly those two fields.

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

The input is a _destructured_ record: `record(id: String, greeting_file: Path)` says "this process takes a record, and from it I need `id` and `greeting_file`".
Each field becomes a variable in the process, so the script refers to `#!groovy ${id}` and `#!groovy ${greeting_file}` directly.
Because `greeting_file` is declared as a `Path`, Nextflow stages it into the task directory, exactly as `path(greeting_file)` did in the tuple.

A destructured input only sees the fields it names.
The records coming from `SAY_HELLO` also carry `name` and `greeting`, but `SHOUT` can't pass them on, so its output records hold only `id` and `shouted`.
When a process's output should keep all the metadata, take the whole record, as `SAY_HELLO` does, and add to it with `+`.
Destructure when a process needs only a few fields and its output doesn't need the rest.

#### 3.5.3. Update the workflow

Back in `main.nf`, update the one place that unpacked a tuple by position:

=== "After"

    ```groovy title="main.nf" linenums="19" hl_lines="1"
        collected = COLLECT_GREETINGS(shouted.map { s -> s.shouted }.collect(), params.batch)
    ```

=== "Before"

    ```groovy title="main.nf" linenums="19" hl_lines="1"
        collected = COLLECT_GREETINGS(shouted.map { _meta, file -> file }.collect(), params.batch)
    ```

`#!groovy shouted.map { s -> s.shouted }` picks the `shouted` file out of each record by name.
The `output {}` block doesn't need to change: when you publish a channel of records, Nextflow publishes every file field they contain.

Run the pipeline:

```bash
nextflow run main.nf -resume
```

??? success "Command output"

    ```console
     N E X T F L O W   ~  version 26.09.1-edge

    Launching `main.nf` [intergalactic_ride] revision: 303038e3ad

    executor >  local (9)
    [75/49abcf] SAY_HELLO (chiara) | 4 of 4 ✔
    [ce/cd2082] SHOUT (alice)      | 4 of 4 ✔
    [30/f5d6d9] COLLECT_GREETINGS  | 1 of 1 ✔

    Outputs:

      /workspaces/training/side-quests/types_and_records/results

      greetings:
        - {name: Chiara, id: chiara, greeting: Ciao, greeting_file: greetings/chiara.txt}
        - {name: Dieter, id: dieter, greeting: Hallo, greeting_file: greetings/dieter.txt}
        - {name: Bruno, id: bruno, greeting: Bonjour, greeting_file: greetings/bruno.txt}
        - {name: Alice, id: alice, greeting: Hello, greeting_file: greetings/alice.txt}

      shouted:
        - {id: chiara, shouted: shouted/chiara-shouted.txt}
        - {id: bruno, shouted: shouted/bruno-shouted.txt}
        - {id: dieter, shouted: shouted/dieter-shouted.txt}
        - {id: alice, shouted: shouted/alice-shouted.txt}

      collected: COLLECTED-batch.txt
    ```

`SAY_HELLO` and `SHOUT` ran again because their definitions changed, and `COLLECT_GREETINGS` ran again because the files staged into it are new.
The type-checking warning is gone, since every process is now typed.

The published outputs are now listed as records.
Each greeting record still carries the person's `name` and `greeting` next to the file, thanks to `+`, while the slimmer `shouted` records hold only the two fields `SHOUT` produces.

### Takeaway

In this section, you've learned:

- **How to declare a record type**: `record Person { ... }` names a set of typed fields, and other files include it by name
- **How to load a samplesheet into records**: a `Channel<Person>` parameter loads the file and checks every row against the type
- **How to carry metadata through a process**: take a whole record (`person: Person`) and emit `person + record(...)`, the record version of `meta + [...]`
- **How to select a subset of fields**: a destructured input such as `record(id: String, greeting_file: Path)` takes only the fields a process needs, and stages `Path` fields

The data now carries its own field names and types, and the processes declare what they need.
But the workflow itself is still a loose sequence of channels in the entry workflow.
In the next section, you'll give that logic a name and a typed interface.

---

## 4. Give the workflow a typed interface

Right now, all the processing happens in the unnamed entry workflow.
To find out what that logic needs and what it produces, you have to read it line by line: a channel of what, with which fields, and what comes out the other end?

A named workflow answers those questions in one place.
It declares its inputs in a `take:` block and its outputs in an `emit:` block, and the caller passes inputs as arguments and receives the outputs as the call's result.
The [Workflows of Workflows](../workflows_of_workflows/index.md) side quest covers named workflows in depth.
Here, the inputs and outputs will have types, so the interface documents what the workflow accepts and produces, and `nextflow lint` can check every call against it, as you'll see in section 5.

### 4.1. Move the processing into a typed workflow

#### 4.1.1. Declare the output types

The workflow will emit the greeting records and the shouted records, so give each of those shapes a name too.
Add two record types to `types.nf`:

=== "After"

    ```groovy title="types.nf" linenums="1" hl_lines="9-19"
    nextflow.enable.types = true

    record Person {
        id: String
        name: String
        greeting: String
    }

    record Greeting {
        id: String
        name: String
        greeting: String
        greeting_file: Path
    }

    record Shouted {
        id: String
        shouted: Path
    }
    ```

=== "Before"

    ```groovy title="types.nf" linenums="1"
    nextflow.enable.types = true

    record Person {
        id: String
        name: String
        greeting: String
    }
    ```

`Greeting` describes what `SAY_HELLO` emits (a person plus their `greeting_file`), and `Shouted` describes what `SHOUT` emits.

`SAY_HELLO` never mentions `Greeting`, and it doesn't need to.
Records match by their fields: a record type names a shape, and any record with those fields, of those types, is accepted where that type is expected, even if it carries extra fields.
The `person + record(greeting_file: ...)` records have exactly the fields of a `Greeting`, so they count as `Greeting` records.

#### 4.1.2. Add the `GREET` workflow

In `main.nf`, include the new types, and add a `GREET` workflow that runs `SAY_HELLO` and `SHOUT`:

=== "After"

    ```groovy title="main.nf" linenums="8" hl_lines="1 8-19"
    include { Person ; Greeting ; Shouted } from './types.nf'

    params {
        input: Channel<Person>
        batch: String = 'batch'
    }

    workflow GREET {
        take:
        people: Channel<Person>

        main:
        greetings = SAY_HELLO(people)
        shouted = SHOUT(greetings)

        emit:
        greetings: Channel<Greeting> = greetings
        shouted: Channel<Shouted> = shouted
    }

    workflow {
    ```

=== "Before"

    ```groovy title="main.nf" linenums="8" hl_lines="1"
    include { Person } from './types.nf'

    params {
        input: Channel<Person>
        batch: String = 'batch'
    }

    workflow {
    ```

The `take:` block declares one input, `people: Channel<Person>`: a channel whose elements are `Person` records.
The `emit:` block declares two outputs, each in the form `name: Type = value`.
`greetings` is a channel of `Greeting` records and `shouted` is a channel of `Shouted` records, so anyone reading the interface can look up exactly which fields each output carries.

When a workflow has more than one output, the call returns a record with one field per output.
Update the entry workflow to call `GREET`, keep the result in `res`, and read its fields:

=== "After"

    ```groovy title="main.nf" linenums="30" hl_lines="1 2 5 6"
        res = GREET(params.input)
        collected = COLLECT_GREETINGS(res.shouted.map { s -> s.shouted }.collect(), params.batch)

        publish:
        greetings = res.greetings
        shouted = res.shouted
        collected = collected
    ```

=== "Before"

    ```groovy title="main.nf" linenums="30" hl_lines="1 2 3 6 7"
        greetings = SAY_HELLO(params.input)
        shouted = SHOUT(greetings)
        collected = COLLECT_GREETINGS(shouted.map { s -> s.shouted }.collect(), params.batch)

        publish:
        greetings = greetings
        shouted = shouted
        collected = collected
    ```

`res.greetings` and `res.shouted` are the two channels emitted by `GREET`.
A workflow with a single output would leave it unnamed (`emit:` followed by the channel alone), and the call would return that channel directly.

Run the pipeline:

```bash
nextflow run main.nf -resume
```

??? success "Command output"

    ```console
     N E X T F L O W   ~  version 26.09.1-edge

    Launching `main.nf` [cranky_cuvier] revision: f00616cc96

    executor >  local (9)
    [dd/c5939f] GREET:SAY_HELLO (dieter) | 4 of 4 ✔
    [e4/0295eb] GREET:SHOUT (alice)      | 4 of 4 ✔
    [60/860f5e] COLLECT_GREETINGS        | 1 of 1 ✔

    Outputs:

      /workspaces/training/side-quests/types_and_records/results

      greetings:
        - {name: Alice, id: alice, greeting: Hello, greeting_file: greetings/alice.txt}
        - {name: Chiara, id: chiara, greeting: Ciao, greeting_file: greetings/chiara.txt}
        - {name: Bruno, id: bruno, greeting: Bonjour, greeting_file: greetings/bruno.txt}
        - {name: Dieter, id: dieter, greeting: Hallo, greeting_file: greetings/dieter.txt}

      shouted:
        - {id: bruno, shouted: shouted/bruno-shouted.txt}
        - {id: chiara, shouted: shouted/chiara-shouted.txt}
        - {id: dieter, shouted: shouted/dieter-shouted.txt}
        - {id: alice, shouted: shouted/alice-shouted.txt}

      collected: COLLECTED-batch.txt
    ```

`SAY_HELLO` and `SHOUT` are now listed as `GREET:SAY_HELLO` and `GREET:SHOUT`, since they run inside the `GREET` workflow.
The process name is part of each task's hash, so every task ran again under its new name, and `COLLECT_GREETINGS` followed.

### 4.2. Pass and return single values

The last step is to move `COLLECT_GREETINGS` into `GREET`, so that the whole greeting logic sits behind one interface.
That needs two new entries in the interface: the batch name going in, and the collected file coming out.
Neither is a stream of items like `people`: there is exactly one of each, so they need a different type.

#### 4.2.1. Channels and values

Nextflow has two dataflow types:

- A **`Channel`** carries any number of items, which arrive one at a time.
  `people` is a channel: one record per person.
- A **`Value`** carries exactly one item, which may not be available yet.
  A `Value` can be used any number of times, and every process task that needs it gets the same item.

You've used both already without naming them.
`collect()` turns a channel into a value: one item containing everything the channel emitted.
And a process that is called only with values runs once and returns a value, which is why `COLLECT_GREETINGS` produces a single file.

#### 4.2.2. Move `COLLECT_GREETINGS` into the workflow

Add a `batch` input of type `Value<String>` and a `collected` output of type `Value<Path>`:

=== "After"

    ```groovy title="main.nf" linenums="15" hl_lines="4 9 14"
    workflow GREET {
        take:
        people: Channel<Person>
        batch: Value<String>

        main:
        greetings = SAY_HELLO(people)
        shouted = SHOUT(greetings)
        collected = COLLECT_GREETINGS(shouted.map { s -> s.shouted }.collect(), batch)

        emit:
        greetings: Channel<Greeting> = greetings
        shouted: Channel<Shouted> = shouted
        collected: Value<Path> = collected
    }
    ```

=== "Before"

    ```groovy title="main.nf" linenums="15"
    workflow GREET {
        take:
        people: Channel<Person>

        main:
        greetings = SAY_HELLO(people)
        shouted = SHOUT(greetings)

        emit:
        greetings: Channel<Greeting> = greetings
        shouted: Channel<Shouted> = shouted
    }
    ```

Then update the entry workflow:

=== "After"

    ```groovy title="main.nf" linenums="33" hl_lines="1 6"
        res = GREET(params.input, channel.value(params.batch))

        publish:
        greetings = res.greetings
        shouted = res.shouted
        collected = res.collected
    ```

=== "Before"

    ```groovy title="main.nf" linenums="33" hl_lines="1 2 7"
        res = GREET(params.input)
        collected = COLLECT_GREETINGS(res.shouted.map { s -> s.shouted }.collect(), params.batch)

        publish:
        greetings = res.greetings
        shouted = res.shouted
        collected = collected
    ```

Workflow inputs are dataflow types, because a workflow's inputs usually come from upstream processes and arrive while the pipeline runs.
`params.batch` is a plain `String`, known before the run starts, so `#!groovy channel.value(params.batch)` wraps it in a `Value<String>` to match the declared input.
The `GREET` interface now reads as a summary of the workflow: it takes a channel of people and one batch name, and it returns a channel of greetings, a channel of shouted greetings and one collected file.

Run the pipeline:

```bash
nextflow run main.nf -resume
```

??? success "Command output"

    ```console
     N E X T F L O W   ~  version 26.09.1-edge

    Launching `main.nf` [desperate_cantor] revision: 9f57d677f2

    executor >  local (1)
    [59/919495] GREET:SAY_HELLO (alice) | 4 of 4, cached: 4 ✔
    [8c/19b40a] GREET:SHOUT (dieter)    | 4 of 4, cached: 4 ✔
    [d4/0d1ca4] GREET:COLLECT_GREETINGS | 1 of 1 ✔

    Outputs:

      /workspaces/training/side-quests/types_and_records/results

      greetings:
        - {name: Bruno, id: bruno, greeting: Bonjour, greeting_file: greetings/bruno.txt}
        - {name: Chiara, id: chiara, greeting: Ciao, greeting_file: greetings/chiara.txt}
        - {name: Dieter, id: dieter, greeting: Hallo, greeting_file: greetings/dieter.txt}
        - {name: Alice, id: alice, greeting: Hello, greeting_file: greetings/alice.txt}

      shouted:
        - {id: chiara, shouted: shouted/chiara-shouted.txt}
        - {id: dieter, shouted: shouted/dieter-shouted.txt}
        - {id: alice, shouted: shouted/alice-shouted.txt}
        - {id: bruno, shouted: shouted/bruno-shouted.txt}

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
    HALLO, DIETER!
    HELLO, ALICE!
    ```

The order of the lines may differ in your run, since the greetings are collected in whatever order the `SHOUT` tasks finish.

### Takeaway

In this section, you've learned:

- **How to type a workflow's inputs**: `take:` declares `Channel<T>` and `Value<T>` inputs
- **How to type its outputs**: named outputs such as `greetings: Channel<Greeting>` tell callers which fields each output carries, and the call returns them as a record (`res.greetings`)
- **The difference between `Channel` and `Value`**: many items over time versus exactly one item

Every process and workflow in the pipeline now declares what it takes and what it returns.
In the last section, you'll see how that information lets `nextflow lint` catch mistakes before they ever reach a run.

---

## 5. Let the type checker help

In section 1, a mistake went unnoticed until someone read the output.
`nextflow lint` parses and type-checks every script in a directory without running anything, and with types in place it can check how your data is used.
In this section you'll make two mistakes on purpose and let the linter find them.

### 5.1. A misspelled field

Misspell a field in the `SAY_HELLO` script:

=== "After"

    ```groovy title="modules/say_hello.nf" linenums="19" hl_lines="1"
        echo '${person.greting}, ${person.name}!' > ${person.id}.txt
    ```

=== "Before"

    ```groovy title="modules/say_hello.nf" linenums="19" hl_lines="1"
        echo '${person.greeting}, ${person.name}!' > ${person.id}.txt
    ```

This is the same kind of mistake as the missing column in section 1.2: `person.greting` evaluates to `null` when the pipeline runs, so the greeting file would say `null, Alice!`.
The difference is that the type checker can now find it before you run anything.
Check the project:

```bash
nextflow lint .
```

??? failure "Command output"

    ```console
    Linting Nextflow code..
    Linting: types.nf
    Linting: modules/shout.nf
    Linting: modules/say_hello.nf
    Linting: modules/collect_greetings.nf
    Linting: nextflow.config
    Linting: main.nf
    Error modules/say_hello.nf:19:13: Unrecognized property `greting` for type Person
    │  19 |     echo '${person.greting}, ${person.name}!' > ${person.id}.txt
    ╰     |             ^^^^^^^^^^^^^^

    Nextflow linting complete!
     ❌ 1 file had 1 error
     ✅ 5 files had no errors
    ```

`person` is declared as a `Person`, so the linter knows which fields it has, and `greting` isn't one of them.
The mistake is reported before a single task runs, instead of showing up in the results.
Fix the typo:

=== "After"

    ```groovy title="modules/say_hello.nf" linenums="19" hl_lines="1"
        echo '${person.greeting}, ${person.name}!' > ${person.id}.txt
    ```

=== "Before"

    ```groovy title="modules/say_hello.nf" linenums="19" hl_lines="1"
        echo '${person.greting}, ${person.name}!' > ${person.id}.txt
    ```

### 5.2. The wrong kind of data

Now wire up `GREET` with the wrong kind of data: a channel of plain names in place of `Person` records, as you might when trying the workflow out with a few hard-coded inputs:

=== "After"

    ```groovy title="main.nf" linenums="33" hl_lines="1"
        res = GREET(channel.of('Alice', 'Bruno'), channel.value(params.batch))
    ```

=== "Before"

    ```groovy title="main.nf" linenums="33" hl_lines="1"
        res = GREET(params.input, channel.value(params.batch))
    ```

Check the project:

```bash
nextflow lint .
```

??? failure "Command output"

    ```console
    Linting Nextflow code..
    Linting: types.nf
    Linting: modules/shout.nf
    Linting: modules/say_hello.nf
    Linting: modules/collect_greetings.nf
    Linting: nextflow.config
    Linting: main.nf
    Error main.nf:33:17: Argument with type Channel<String> is not compatible with parameter of type Channel<Person>
    │  33 |     res = GREET(channel.of('Alice', 'Bruno'), channel.value(params.bat
    ╰     |                 ^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    Nextflow linting complete!
     ❌ 1 file had 1 error
     ✅ 5 files had no errors
    ```

The linter compares each argument with the type declared in `GREET`'s `take:` block.
A channel of strings is not a channel of `Person` records, so the call is reported without running anything, and the message says which type the workflow expects.
Put the original argument back:

=== "After"

    ```groovy title="main.nf" linenums="33" hl_lines="1"
        res = GREET(params.input, channel.value(params.batch))
    ```

=== "Before"

    ```groovy title="main.nf" linenums="33" hl_lines="1"
        res = GREET(channel.of('Alice', 'Bruno'), channel.value(params.batch))
    ```

Run the linter one last time:

```bash
nextflow lint .
```

??? success "Command output"

    ```console
    Linting Nextflow code..
    Linting: types.nf
    Linting: modules/shout.nf
    Linting: modules/say_hello.nf
    Linting: modules/collect_greetings.nf
    Linting: nextflow.config
    Linting: main.nf
    Nextflow linting complete!
     ✅ 6 files had no errors
    ```

!!! tip

    As you saw in section 3, when you run a pipeline Nextflow also type-checks it, but it only prints a warning that tells you how many errors it found and to run `nextflow lint`, then carries on.
    With the `person.greting` typo from section 5.1, `nextflow run` prints that warning, runs every task, and still writes `null, Alice!`.
    Run `nextflow lint` yourself to see the errors and fix them before you run.

### Takeaway

In this section, you've seen the type checker catch mistakes before a run:

- **A misspelled record field**: reported because the variable's record type lists its fields
- **An argument of the wrong type**: reported because the workflow's `take:` block declares a type for each input

Combined with the samplesheet check from section 3.3, these catch the kinds of mistakes that slipped through silently in section 1.
Run `nextflow lint` before you run a pipeline: the more of it is typed, the more it can check.

---

## Summary

In this side quest, you've migrated an untyped pipeline to static types and records.
You started with a pipeline that wrote `null, Alice!` when a samplesheet column was renamed, and ended with one that rejects that samplesheet before running a task, naming the missing field, with interfaces that say what goes in and what comes out.

Compared with meta maps and tuples, typed code with records:

- Names every field, so producers and consumers never depend on element order
- Lets each process take a whole record or select only the fields it needs
- Checks samplesheet rows against a record type, so a missing value is an error instead of `null` in your results
- Gives `nextflow lint` enough information to catch mistakes before you run anything

### Key patterns

1.  **Enabling static typing**: Add the feature flag to each script that uses typed processes, workflows or record types.

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

3.  **Record types and samplesheet parameters**: Declare the shape once, include it where it's used, and let a `Channel` parameter load and check the samplesheet.

    ```groovy
    record Person {
        id: String
        name: String
        greeting: String
    }

    params {
        input: Channel<Person>
    }
    ```

4.  **Carrying metadata through a process**: Take the whole record and add fields with `+`; the right-hand side wins when both records have the same field.

    ```groovy
    input:
    person: Person

    output:
    person + record(greeting_file: file("${person.id}.txt"))
    ```

5.  **Selecting a subset of fields**: Destructure only the fields a process needs.

    ```groovy
    input:
    record(id: String, greeting_file: Path)
    ```

6.  **Typed workflows**: Use `Channel<T>` and `Value<T>` in `take:` and `emit:`, and read multiple outputs from the returned record.

    ```groovy
    workflow GREET {
        take:
        people: Channel<Person>
        batch: Value<String>

        main:
        // ...

        emit:
        greetings: Channel<Greeting> = greetings
        collected: Value<Path> = collected
    }
    ```

    The caller reads multiple outputs from the returned record: `res = GREET(params.input, channel.value(params.batch))`, then `res.greetings`.

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
