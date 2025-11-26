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

#### Review the materials

You'll find a `modules` directory containing several process definitions that build upon what you learned in 'Hello Nextflow':

```console title="Directory contents"
modules/
├── say_hello.nf             # Creates a greeting (from Hello Nextflow)
├── say_hello_upper.nf       # Converts to uppercase (from Hello Nextflow)
├── timestamp_greeting.nf    # Adds timestamps to greetings
├── validate_name.nf         # Validates input names
└── reverse_text.nf          # Reverses text content
```

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

## 1. Create the Greeting Workflow

Let's start by creating a workflow that validates names and generates timestamped greetings.

### 1.1. Create the workflow structure

```bash title="Create workflow directory and file"
mkdir -p workflows
touch workflows/greeting.nf
```

### 1.2. Add the first (sub)workflow code

Add this code to `workflows/greeting.nf`:

```groovy title="workflows/greeting.nf" linenums="1"
include { VALIDATE_NAME } from '../modules/validate_name'
include { SAY_HELLO } from '../modules/say_hello'
include { TIMESTAMP_GREETING } from '../modules/timestamp_greeting'

workflow {

    names_ch = channel.of('Alice', 'Bob', 'Charlie')

    // Chain processes: validate -> create greeting -> add timestamp
    validated_ch = VALIDATE_NAME(names_ch)
    greetings_ch = SAY_HELLO(validated_ch)
    timestamped_ch = TIMESTAMP_GREETING(greetings_ch)
}
```

This is a complete workflow, with a structure similar to the ones you saw in the 'Hello Nextflow' tutorial, that we can test independently. Let's try that now:

```bash title="Run the greeting workflow"
nextflow run workflows/greeting.nf
```

```console title="Expected output"
N E X T F L O W  ~  version 24.10.0
Launching `workflows/greeting.nf` [peaceful_montalcini] DSL2 - revision: 90f61b7093
executor >  local (9)
[51/4f980f] process > VALIDATE_NAME (validating Bob)                    [100%] 3 of 3 ✔
[2b/dd8dc2] process > SAY_HELLO (greeting Bob)                          [100%] 3 of 3 ✔
[8e/882565] process > TIMESTAMP_GREETING (adding timestamp to greeting) [100%] 3 of 3 ✔
```

This works as expected, but to make it composable there's a few things we need to change.

### 1.3. Make the workflow composable

Composable workflows have some differences from the ones you saw in the 'Hello Nextflow' tutorial:

- The workflow block needs to be named
- Inputs are declared using the `take:` keyword
- Workflow content is placed inside the `main:` block
- Outputs are declared using the `emit:` keyword

Let's update the greeting workflow to match this structure. Change the code to the following:

```groovy title="workflows/greeting.nf" linenums="1" hl_lines="6 7 9 15 16 17"
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

You can see that the workflow is now named and has a `take:` and `emit:` block, and these are the connections we will use to compose a higher level workflow.
The workflow content is also placed inside the `main:` block. Note also that we have removed the `names_ch` input channel declaration, as it's now passed as an argument to the workflow.

Let's test the workflow again to see if it works as expected:

```bash title="Run the greeting workflow"
nextflow run workflows/greeting.nf
```

What you'll actually see in response is:

```console title="Expected output"
N E X T F L O W  ~  version 24.10.0
Launching `workflows/greeting.nf` [high_brahmagupta] DSL2 - revision: 8f5857af25
WARN: No entry workflow specified
```

This tells you about another new concept, an 'entry workflow'. The entry workflow is the workflow that gets called when you run a Nextflow script. By default, Nextflow will use an un-named workflow as the entry workflow, when present, and that's what you've been doing so far, with workflow blocks starting like this:

```groovy title="hello.nf" linenums="1"
workflow {
```

But our greeting workflow doesn't have an un-named workflow, rather we have a named workflow:

```groovy title="workflows/greeting.nf" linenums="1"
workflow GREETING_WORKFLOW {
```

... so Nextflow will throw an error. We can actually tell Nextflow to use our named workflow as the entry workflow by adding this line to Nextflow's command line:

```bash title="Run the greeting workflow"
nextflow run workflows/greeting.nf -entry GREETING_WORKFLOW
```

This will also throw an error, because the workflow is expecting an input channel:

```console title="Expected output"
N E X T F L O W  ~  version 24.10.0
Launching `workflows/greeting.nf` [compassionate_fermi] DSL2 - revision: 8f5857af25
ERROR ~ Workflow `GREETING_WORKFLOW` declares 1 input channels but 0 were given

 -- Check '.nextflow.log' file for details
```

... but if you wanted to call a named workflow that didn't require inputs, you could call it this way.

But we didn't add that syntax so we could call the workflow directly, we did it so we could compose it with other workflows. Let's start by creating a main workflow that imports and uses the `greeting` workflow.

### 1.4. Create and test the main workflow

Now we will create a main workflow that imports and uses the `greeting` workflow.

Create `main.nf`:

```groovy title="main.nf" linenums="1"
include { GREETING_WORKFLOW } from './workflows/greeting'

workflow {
    names = channel.of('Alice', 'Bob', 'Charlie')
    GREETING_WORKFLOW(names)

    GREETING_WORKFLOW.out.greetings.view { "Original: $it" }
    GREETING_WORKFLOW.out.timestamped.view { "Timestamped: $it" }
}

```

Note that our workflow entry in this file is un-named, and that's because we're going to use it as an entry workflow.

Run this and see the output:

```bash title="Run the workflow"
nextflow run main.nf
```

```console title="Expected output"
N E X T F L O W  ~  version 24.10.0
Launching `main.nf` [goofy_mayer] DSL2 - revision: 543f8742fe
executor >  local (9)
[05/3cc752] process > GREETING_WORKFLOW:VALIDATE_NAME (validating Char... [100%] 3 of 3 ✔
[b1/b56ecf] process > GREETING_WORKFLOW:SAY_HELLO (greeting Charlie)      [100%] 3 of 3 ✔
[ea/342168] process > GREETING_WORKFLOW:TIMESTAMP_GREETING (adding tim... [100%] 3 of 3 ✔
Original: /workspaces/training/side_quests/workflows_of_workflows/work/bb/c8aff3df0ebc15a4d7d35f736db44c/Alice-output.txt
Original: /workspaces/training/side_quests/workflows_of_workflows/work/fb/fa877776e8a5d90b537b1bcd3b6f5b/Bob-output.txt
Original: /workspaces/training/side_quests/workflows_of_workflows/work/b1/b56ecf938fda8bcbec211847c8f0be/Charlie-output.txt
Timestamped: /workspaces/training/side_quests/workflows_of_workflows/work/06/877bc909f140bbf8223343450cea36/timestamped_Alice-output.txt
Timestamped: /workspaces/training/side_quests/workflows_of_workflows/work/aa/bd31b71cdb745b7c155ca7f8837b8a/timestamped_Bob-output.txt
Timestamped: /workspaces/training/side_quests/workflows_of_workflows/work/ea/342168d4ba04cc899a89c56cbfd9b0/timestamped_Charlie-output.txt
```

It works! We've wrapped the named greeting workflow in a main workflow with an un-named entry `workflow` block. The main workflow is using the `GREETING_WORKFLOW` workflow almost (not quite) like a process, and passing the `names` channel as an argument.

### Takeaway

In this section, you've learned several important concepts:

- **Named Workflows**: Creating a named workflow (`GREETING_WORKFLOW`) that can be imported and reused
- **Workflow Interfaces**: Defining clear inputs with `take:` and outputs with `emit:` to create a composable workflow
- **Entry Points**: Understanding that Nextflow needs an entry workflow (either unnamed or specified with `-entry`)
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

## 2. Add the Transform Workflow

Now let's create a workflow that applies text transformations to the greetings.

### 2.1. Create the workflow file

```bash title="Create transform workflow"
touch workflows/transform.nf
```

### 2.2. Add the workflow code

Add this code to `workflows/transform.nf`:

```groovy title="workflows/transform.nf" linenums="1"
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

We won't repeat the explanation of the composable syntax here, but note the named workflow is again declared with a `take:` and `emit:` block, and the workflow content is placed inside the `main:` block.

### 2.3. Update the main workflow

Update `main.nf` to use both workflows:

```groovy title="main.nf" linenums="1"
include { GREETING_WORKFLOW } from './workflows/greeting'
include { TRANSFORM_WORKFLOW } from './workflows/transform'

workflow {
    names = channel.of('Alice', 'Bob', 'Charlie')

    // Run the greeting workflow
    GREETING_WORKFLOW(names)

    // Run the transform workflow
    TRANSFORM_WORKFLOW(GREETING_WORKFLOW.out.timestamped)

    // View results
    TRANSFORM_WORKFLOW.out.upper.view { "Uppercase: $it" }
    TRANSFORM_WORKFLOW.out.reversed.view { "Reversed: $it" }
}
```

Run the complete pipeline:

```bash title="Run the workflow"
nextflow run main.nf
```

```console title="Expected output"
N E X T F L O W  ~  version 24.10.0
Launching `main.nf` [sick_kimura] DSL2 - revision: 8dc45fc6a8
executor >  local (13)
executor >  local (15)
[83/1b51f4] process > GREETING_WORKFLOW:VALIDATE_NAME (validating Alice)  [100%] 3 of 3 ✔
[68/556150] process > GREETING_WORKFLOW:SAY_HELLO (greeting Alice)        [100%] 3 of 3 ✔
[de/511abd] process > GREETING_WORKFLOW:TIMESTAMP_GREETING (adding tim... [100%] 3 of 3 ✔
[cd/e6a7e0] process > TRANSFORM_WORKFLOW:SAY_HELLO_UPPER (converting t... [100%] 3 of 3 ✔
[f0/74ba4a] process > TRANSFORM_WORKFLOW:REVERSE_TEXT (reversing UPPER... [100%] 3 of 3 ✔
Uppercase: /workspaces/training/side_quests/workflows_of_workflows/work/a0/d4f5df4d6344604498fa47a6084a11/UPPER-timestamped_Bob-output.txt
Uppercase: /workspaces/training/side_quests/workflows_of_workflows/work/69/b5e37f6c79c2fd38adb75d0eca8f87/UPPER-timestamped_Charlie-output.txt
Uppercase: /workspaces/training/side_quests/workflows_of_workflows/work/cd/e6a7e0b17e7d5a2f71bb8123cd53a7/UPPER-timestamped_Alice-output.txt
Reversed: /workspaces/training/side_quests/workflows_of_workflows/work/7a/7a222f7957b35d1d121338566a24ac/REVERSED-UPPER-timestamped_Bob-output.txt
Reversed: /workspaces/training/side_quests/workflows_of_workflows/work/46/8d19af62e33a5a6417c773496e0f90/REVERSED-UPPER-timestamped_Charlie-output.txt
Reversed: /workspaces/training/side_quests/workflows_of_workflows/work/f0/74ba4a10d9ef5c82f829d1c154d0f6/REVERSED-UPPER-timestamped_Alice-output.txt
```

If you take a look at one of those reversed files, you'll see that it's the uppercase version of the greeting reversed:

```bash title="Check a final output file"
cat /workspaces/training/side_quests/workflows_of_workflows/work/f0/74ba4a10d9ef5c82f829d1c154d0f6/REVERSED-UPPER-timestamped_Alice-output.txt
```

```console title="Reversed file content"
!ECILA ,OLLEH ]04:50:71 60-30-5202[
```

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

_It's important to note however that while calling workflows is a bit like calling processes, it's not actually the same thing. You can't, for example, run a workflow N times by calling it with a channel of size N - you would need to pass a channel of size N to the workflow and iterate internally._

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

3.  **Entry points**: We learned that Nextflow requires an entry workflow (either an unnamed workflow or a named workflow specified with `-entry`) to know where to start execution.

    - Unnamed workflow (default entry point)

    ```groovy
    workflow {
        // This is automatically the entry point when the script is run
    }
    ```

    - Named workflow (not an entry point by default)

    ```groovy
    workflow NAMED_WORKFLOW {
        // This is not automatically run
    }
    ```

    - Running a named workflow as entry point (bash command)

    ```bash
    nextflow run script.nf -entry NAMED_WORKFLOW
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

Return to the [menu of Side Quests](./index.md) or click the button in the bottom right of the page to move on to the next topic in the list.
