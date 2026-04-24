# Workflows of Workflows

When you're developing a pipeline, you often find yourself creating similar sequences of processes for different data types or analysis steps. You might end up copying and pasting these process sequences, leading to duplicated code that's hard to maintain; or you might create one massive workflow that's difficult to understand and modify.

One of the most powerful features of Nextflow is its ability to compose complex pipelines from smaller, reusable workflow modules. This modular approach makes pipelines easier to develop, test, and maintain.

### Learning goals

In this side quest, we'll explore how to develop workflow modules that can be tested and used separately, compose those modules into a larger pipeline, and manage data flow between modules.

By the end of this side quest, you'll be able to:

- Break down complex pipelines into logical, reusable units
- Test each workflow module independently
- Mix and match workflows to create new pipelines
- Share common workflow modules across different pipelines
- Make your code more maintainable and easier to understand

These skills will help you build complex pipelines while maintaining clean, maintainable code structure.

### Prerequisites

Before taking on this side quest you should:

- Have completed the [Hello Nextflow](../hello_nextflow/README.md) tutorial or equivalent beginner's course.
- Be comfortable using basic Nextflow concepts and mechanisms (processes, channels, operators, modules)

---

## 0. Get started

#### Open the training codespace

If you haven't yet done so, make sure to open the training environment as described in the [Environment Setup](../envsetup/index.md).

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### Move into the project directory

Let's move into the directory where the files for this tutorial are located.

```bash
cd side-quests/workflows_of_workflows
```

You can set VSCode to focus on this directory:

```bash
code .
```

The editor opens with the project directory in focus.

#### Review the materials

You'll find a `modules` directory with process definitions, a `workflows` directory with two pre-written workflow scripts, and a `main.nf` file that you will progressively update:

```console title="Directory contents"
├── main.nf
├── workflows/
│   ├── greeting.nf              # Standalone greeting workflow (to be made composable)
│   └── transform.nf             # Standalone transform workflow (to be made composable)
└── modules/
    ├── say_hello.nf             # Creates a greeting (from Hello Nextflow)
    ├── say_hello_upper.nf       # Converts to uppercase (from Hello Nextflow)
    ├── timestamp_greeting.nf    # Adds timestamps to greetings
    ├── validate_name.nf         # Validates input names
    └── reverse_text.nf          # Reverses text content
```

The `modules/` directory contains the individual process definitions, and the `workflows/` directory contains the two pre-written workflow scripts you will work with in this side quest.

#### Review the assignment

Your challenge is to assemble these modules into two separate workflows that we will then compose into a main workflow:

- A `GREETING_WORKFLOW` that validates names, creates greetings, and adds timestamps
- A `TRANSFORM_WORKFLOW` that converts text to uppercase and reverses it

<!-- TODO: give a bit more details, similar to how it's done in the Metadata side quest -->

#### Readiness checklist

Think you're ready to dive in?

- [ ] I understand the goal of this course and its prerequisites
- [ ] My codespace is up and running
- [ ] I've set my working directory appropriately
- [ ] I understand the assignment

If you can check all the boxes, you're good to go.

---

## 1. The Greeting Workflow

The greeting workflow validates names and generates timestamped greetings.

### 1.1. Review the greeting workflow

Open `workflows/greeting.nf` and take a look at the code:

```groovy title="workflows/greeting.nf" linenums="1"
include { VALIDATE_NAME } from '../modules/validate_name'
include { SAY_HELLO } from '../modules/say_hello'
include { TIMESTAMP_GREETING } from '../modules/timestamp_greeting'

workflow {
    main:
    names_ch = channel.of('Alice', 'Bob', 'Charlie')

    // Chain processes: validate -> create greeting -> add timestamp
    validated_ch = VALIDATE_NAME(names_ch)
    greetings_ch = SAY_HELLO(validated_ch)
    timestamped_ch = TIMESTAMP_GREETING(greetings_ch)

    publish:
    greetings = greetings_ch
    timestamped = timestamped_ch
}

output {
    greetings {
    }
    timestamped {
    }
}
```

This is a complete, self-contained workflow with the same structure you saw in the 'Hello Nextflow' tutorial.
It hardcodes the input names, chains three processes, and publishes two outputs.

### 1.2. Run the standalone greeting workflow

Run it to verify everything works:

```bash
nextflow run workflows/greeting.nf
```

??? success "Command output"

    ```console
    N E X T F L O W  ~  version 24.10.0
    Launching `workflows/greeting.nf` [peaceful_montalcini] DSL2 - revision: 90f61b7093
    executor >  local (9)
    [51/4f980f] process > VALIDATE_NAME (validating Bob)                    [100%] 3 of 3 ✔
    [2b/dd8dc2] process > SAY_HELLO (greeting Bob)                          [100%] 3 of 3 ✔
    [8e/882565] process > TIMESTAMP_GREETING (adding timestamp to greeting) [100%] 3 of 3 ✔
    ```

This works as expected.
To make it composable with other workflows, a few things need to change.

### 1.3. Make the workflow composable

To make a workflow composable, four things need to change:
the workflow gets a name, inputs move to a `take:` block, outputs move to an `emit:` block,
and the standalone `publish:`/`output {}` blocks are removed (they belong in the entry workflow).

Let's walk through these changes one by one.

#### 1.3.1. Name the workflow

Give the workflow a name so it can be imported from a parent workflow.

=== "After"

    ```groovy title="workflows/greeting.nf" linenums="5" hl_lines="1"
    workflow GREETING_WORKFLOW {
    ```

=== "Before"

    ```groovy title="workflows/greeting.nf" linenums="5" hl_lines="1"
    workflow {
    ```

With a name, the workflow can be imported into other scripts.

#### 1.3.2. Declare inputs with `take:`

Replace the hardcoded channel declaration with a `take:` block that declares what inputs the workflow expects.
The `take:` block goes before `main:`, and the `names_ch = channel.of(...)` line is removed.

=== "After"

    ```groovy title="workflows/greeting.nf" linenums="5" hl_lines="2 3 5"
    workflow GREETING_WORKFLOW {
        take:
        names_ch        // Input channel with names

        main:
        // Chain processes: validate -> create greeting -> add timestamp
        validated_ch = VALIDATE_NAME(names_ch)
        greetings_ch = SAY_HELLO(validated_ch)
        timestamped_ch = TIMESTAMP_GREETING(greetings_ch)
    ```

=== "Before"

    ```groovy title="workflows/greeting.nf" linenums="5"
    workflow GREETING_WORKFLOW {
        main:
        names_ch = channel.of('Alice', 'Bob', 'Charlie')

        // Chain processes: validate -> create greeting -> add timestamp
        validated_ch = VALIDATE_NAME(names_ch)
        greetings_ch = SAY_HELLO(validated_ch)
        timestamped_ch = TIMESTAMP_GREETING(greetings_ch)
    ```

The `take:` block declares the channel by name only — the details of what goes into it will be defined by the parent workflow.

#### 1.3.3. Declare outputs with `emit:`

Replace the `publish:` section and remove the `output {}` block, replacing them with an `emit:` block that names the outputs.

=== "After"

    ```groovy title="workflows/greeting.nf" linenums="14" hl_lines="2 3 4"

        emit:
        greetings = greetings_ch      // Original greetings
        timestamped = timestamped_ch  // Timestamped greetings
    }
    ```

=== "Before"

    ```groovy title="workflows/greeting.nf" linenums="14"

        publish:
        greetings = greetings_ch
        timestamped = timestamped_ch
    }

    output {
        greetings {
        }
        timestamped {
        }
    }
    ```

The `emit:` block exposes named outputs that parent workflows can access via `GREETING_WORKFLOW.out.greetings` and `GREETING_WORKFLOW.out.timestamped`.

#### 1.3.4. Verify the result and test it

After all three changes, the complete file should look like this:

```groovy title="workflows/greeting.nf" linenums="1" hl_lines="5 6 7 9 15 16 17"
include { VALIDATE_NAME } from '../modules/validate_name'
include { SAY_HELLO } from '../modules/say_hello'
include { TIMESTAMP_GREETING } from '../modules/timestamp_greeting'

workflow GREETING_WORKFLOW {
    take:
    names_ch        // Input channel with names

    main:
    // Chain processes: validate -> create greeting -> add timestamp
    validated_ch = VALIDATE_NAME(names_ch)
    greetings_ch = SAY_HELLO(validated_ch)
    timestamped_ch = TIMESTAMP_GREETING(greetings_ch)

    emit:
    greetings = greetings_ch      // Original greetings
    timestamped = timestamped_ch  // Timestamped greetings
}
```

Now try running it directly:

```bash
nextflow run workflows/greeting.nf
```

??? failure "Command output"

    ```console
    N E X T F L O W  ~  version 24.10.0
    Launching `workflows/greeting.nf` [high_brahmagupta] DSL2 - revision: 8f5857af25
    No entry workflow specified
    ```

This introduces a key concept: the **entry workflow**.
Nextflow uses an unnamed `workflow {}` block as the entry point when you run a script directly.
`GREETING_WORKFLOW` is named, so Nextflow doesn't know how to run it on its own.

That's intentional — composable workflows are designed to be called from an entry workflow, not run directly.
The solution is an entry workflow in `main.nf` that imports and calls `GREETING_WORKFLOW`.

### 1.4. Update and test the main workflow

Now we will update the main workflow to import and use the `greeting` workflow.

#### 1.4.1. Edit `main.nf`

Open `main.nf` and make the following changes:

=== "After"

    ```groovy title="main.nf" linenums="1"
    include { GREETING_WORKFLOW } from './workflows/greeting'

    workflow {
        main:
        names = channel.of('Alice', 'Bob', 'Charlie')
        GREETING_WORKFLOW(names)

        publish:
        greetings = GREETING_WORKFLOW.out.greetings
    }

    output {
        greetings {
        }
    }
    ```

=== "Before"

    ```groovy title="main.nf" linenums="1"
    workflow {
        main:
        names = channel.of('Alice', 'Bob', 'Charlie')

        publish:
        greetings = channel.empty()
    }

    output {
        greetings {
        }
    }
    ```

The entry workflow is un-named because it's used as the pipeline entry point.
The `publish:` section wires the `greetings` output to the GREETING_WORKFLOW result, replacing the `channel.empty()` placeholder.

#### 1.4.2. Run the workflow

Run the workflow to test that it works:

```bash
nextflow run main.nf
```

??? success "Command output"

    ```console
    N E X T F L O W  ~  version 24.10.0
    Launching `main.nf` [goofy_mayer] DSL2 - revision: 543f8742fe
    executor >  local (9)
    [05/3cc752] process > GREETING_WORKFLOW:VALIDATE_NAME (validating Char... [100%] 3 of 3 ✔
    [b1/b56ecf] process > GREETING_WORKFLOW:SAY_HELLO (greeting Charlie)      [100%] 3 of 3 ✔
    [ea/342168] process > GREETING_WORKFLOW:TIMESTAMP_GREETING (adding tim... [100%] 3 of 3 ✔
    ```

??? abstract "Directory contents"

    ```console
    results/
    └── greetings
        ├── Alice-output.txt
        ├── Bob-output.txt
        └── Charlie-output.txt
    ```

??? abstract "File contents"

    ```console title="results/greetings/Alice-output.txt"
    Hello, Alice!
    ```

The greeting files are published to `results/greetings/`.
The main workflow calls `GREETING_WORKFLOW` and wires its output directly to the `publish:` section.

### Takeaway

In this section, you've learned several important concepts:

- **Named Workflows**: Creating a named workflow (`GREETING_WORKFLOW`) that can be imported and reused
- **Workflow Interfaces**: Defining clear inputs with `take:` and outputs with `emit:` to create a composable workflow
- **Entry Points**: Understanding that Nextflow needs an unnamed entry workflow to run a script
- **Workflow Composition**: Importing and using a named workflow within another workflow
- **Workflow Namespaces**: Accessing workflow outputs using the `.out` namespace (`GREETING_WORKFLOW.out.greetings`)

You now have a working greeting workflow that:

- Takes a channel of names as input
- Validates each name
- Creates a greeting for each valid name
- Adds timestamps to the greetings
- Exposes both original and timestamped greetings as outputs

This modular approach allows you to test the greeting workflow independently or use it as a component in larger pipelines.

---

## 2. Adapt the TRANSFORM workflow

The transform workflow applies text transformations to the timestamped greetings.

### 2.1. Review the workflow

Open `workflows/transform.nf` and take a look at the code:

```groovy title="workflows/transform.nf" linenums="1"
include { SAY_HELLO_UPPER } from '../modules/say_hello_upper'
include { REVERSE_TEXT } from '../modules/reverse_text'

workflow {
    main:
    input_ch = channel.fromPath('results/timestamped/*.txt')

    // Apply transformations in sequence
    upper_ch = SAY_HELLO_UPPER(input_ch)
    reversed_ch = REVERSE_TEXT(upper_ch)

    publish:
    upper = upper_ch
    reversed = reversed_ch
}

output {
    upper {
    }
    reversed {
    }
}
```

This standalone workflow reads timestamped greeting files from the `results/` directory produced by `greeting.nf`, converts them to uppercase, then reverses the text.

### 2.2. Make it composable

Apply the same three changes as in section 1.3: name the workflow, replace the hardcoded input with `take:`, and replace `publish:`/`output {}` with `emit:`.

The finished file should look like this:

```groovy title="workflows/transform.nf" linenums="1" hl_lines="4 5 6 8 13 14 15"
include { SAY_HELLO_UPPER } from '../modules/say_hello_upper'
include { REVERSE_TEXT } from '../modules/reverse_text'

workflow TRANSFORM_WORKFLOW {
    take:
    input_ch         // Input channel with messages

    main:
    // Apply transformations in sequence
    upper_ch = SAY_HELLO_UPPER(input_ch)
    reversed_ch = REVERSE_TEXT(upper_ch)

    emit:
    upper = upper_ch        // Uppercase greetings
    reversed = reversed_ch  // Reversed uppercase greetings
}
```

The transform workflow is now composable and ready to be imported into the main workflow.

### 2.3. Update the main workflow

#### 2.3.1. Edit `main.nf`

Update `main.nf` to use both workflows:

=== "After"

    ```groovy title="main.nf" linenums="1"
    include { GREETING_WORKFLOW } from './workflows/greeting'
    include { TRANSFORM_WORKFLOW } from './workflows/transform'

    workflow {
        main:
        names = channel.of('Alice', 'Bob', 'Charlie')

        // Run the greeting workflow
        GREETING_WORKFLOW(names)

        // Run the transform workflow
        TRANSFORM_WORKFLOW(GREETING_WORKFLOW.out.timestamped)

        publish:
        greetings = GREETING_WORKFLOW.out.greetings
        upper = TRANSFORM_WORKFLOW.out.upper
        reversed = TRANSFORM_WORKFLOW.out.reversed
    }

    output {
        greetings {
        }
        upper {
        }
        reversed {
        }
    }
    ```

=== "Before"

    ```groovy title="main.nf" linenums="1"
    include { GREETING_WORKFLOW } from './workflows/greeting'

    workflow {
        main:
        names = channel.of('Alice', 'Bob', 'Charlie')
        GREETING_WORKFLOW(names)

        publish:
        greetings = GREETING_WORKFLOW.out.greetings
    }

    output {
        greetings {
        }
    }
    ```

#### 2.3.2. Run the complete pipeline

Run the pipeline to test that it all works:

```bash
nextflow run main.nf
```

??? success "Command output"

    ```console
    N E X T F L O W  ~  version 24.10.0
    Launching `main.nf` [sick_kimura] DSL2 - revision: 8dc45fc6a8
    executor >  local (15)
    [83/1b51f4] process > GREETING_WORKFLOW:VALIDATE_NAME (validating Alice)  [100%] 3 of 3 ✔
    [68/556150] process > GREETING_WORKFLOW:SAY_HELLO (greeting Alice)        [100%] 3 of 3 ✔
    [de/511abd] process > GREETING_WORKFLOW:TIMESTAMP_GREETING (adding tim... [100%] 3 of 3 ✔
    [cd/e6a7e0] process > TRANSFORM_WORKFLOW:SAY_HELLO_UPPER (converting t... [100%] 3 of 3 ✔
    [f0/74ba4a] process > TRANSFORM_WORKFLOW:REVERSE_TEXT (reversing UPPER... [100%] 3 of 3 ✔
    ```

??? abstract "Directory contents"

    ```console
    results/
    ├── greetings
    │   ├── Alice-output.txt
    │   ├── Bob-output.txt
    │   └── Charlie-output.txt
    ├── reversed
    │   ├── REVERSED-UPPER-timestamped_Alice-output.txt
    │   ├── REVERSED-UPPER-timestamped_Bob-output.txt
    │   └── REVERSED-UPPER-timestamped_Charlie-output.txt
    └── upper
        ├── UPPER-timestamped_Alice-output.txt
        ├── UPPER-timestamped_Bob-output.txt
        └── UPPER-timestamped_Charlie-output.txt
    ```

??? abstract "File contents"

    ```console title="results/reversed/REVERSED-UPPER-timestamped_Alice-output.txt"
    !ECILA ,OLLEH ]04:50:71 60-30-5202[
    ```

The pipeline is working end-to-end: the greeting has been uppercased and reversed.

### Takeaway

You should now have a complete pipeline that:

- Processes names through the greeting workflow
- Feeds the timestamped greetings into the transform workflow
- Produces both uppercase and reversed versions of the greetings

---

## Summary

In this side quest, we've explored the powerful concept of workflow composition in Nextflow, which allows us to build complex pipelines from smaller, reusable components.

This modular approach offers several advantages over monolithic pipelines:

- Each workflow can be developed, tested, and debugged independently
- Workflows can be reused across different pipelines
- The overall pipeline structure becomes more readable and maintainable
- Changes to one workflow don't necessarily affect others if the interfaces remain consistent
- Entry points can be configured to run different parts of your pipeline as needed

It's important to note that while calling workflows is a bit like calling processes, it's not actually the same thing. You can't, for example, run a workflow N times by calling it with a channel of size N - you would need to pass a channel of size N to the workflow and iterate internally.

Applying these techniques in your own work will enable you to build more sophisticated Nextflow pipelines that can handle complex bioinformatics tasks while remaining maintainable and scalable.

### Key patterns

1.  **Workflow structure**: We defined clear inputs and outputs for each workflow using the `take:` and `emit:` syntax, creating well-defined interfaces between components, and wrapped workflow logic within the `main:` block.

    ```groovy
    workflow EXAMPLE_WORKFLOW {
        take:
            // Input channels are declared here
            input_ch

        main:
            // Workflow logic goes here
            // This is where processes are called and channels are manipulated
            result_ch = SOME_PROCESS(input_ch)

        emit:
            // Output channels are declared here
            output_ch = result_ch
    }
    ```

2.  **Workflow imports:** We built two independent workflow modules and imported them into a main pipeline with include statements.

    - Include a single workflow

    ```groovy
    include { WORKFLOW_NAME } from './path/to/workflow'
    ```

    - Include multiple workflows

    ```groovy
    include { WORKFLOW_A; WORKFLOW_B } from './path/to/workflows'
    ```

    - Include with alias to avoid name conflicts

    ```groovy
    include { WORKFLOW_A as WORKFLOW_A_ALIAS } from './path/to/workflow'
    ```

3.  **Entry points**: Nextflow requires an unnamed entry workflow to know where to start execution. This entry workflow calls your named workflows.

    - Unnamed workflow (entry point)

    ```groovy
    workflow {
        // This is the entry point when the script is run
        NAMED_WORKFLOW(input_ch)
    }
    ```

    - Named workflow (called from entry workflow)

    ```groovy
    workflow NAMED_WORKFLOW {
        // Must be called from the entry workflow
    }
    ```

4.  **Managing data flow:** We learned how to access workflow outputs using the namespace notation (`WORKFLOW_NAME.out.channel_name`) and pass them to other workflows.

    ```nextflow
    WORKFLOW_A(input_ch)
    WORKFLOW_B(WORKFLOW_A.out.some_channel)
    ```

### Additional resources

- [Nextflow Workflow Documentation](https://www.nextflow.io/docs/latest/workflow.html)
- [Channel Operators Reference](https://www.nextflow.io/docs/latest/operator.html)
- [Error Strategy Documentation](https://www.nextflow.io/docs/latest/process.html#errorstrategy)

---

## What's next?

Return to the [menu of Side Quests](../) or click the button in the bottom right of the page to move on to the next topic in the list.
