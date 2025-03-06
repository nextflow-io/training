# Side Quest: Workflows of Workflows

One of the most powerful features of Nextflow is its ability to compose complex pipelines from smaller, reusable workflow modules. This modular approach makes pipelines easier to develop, test, and maintain.

Let's explore why workflow composition is so important. When you're developing a pipeline, you often find yourself creating similar sequences of processes for different data types or analysis steps. Without workflow composition, you might end up copying and pasting these process sequences, leading to duplicated code that's hard to maintain. Or you might create one massive workflow that's difficult to understand and modify.

With workflow composition, you can:
- Break down complex pipelines into logical, reusable units
- Test each workflow module independently
- Mix and match workflows to create new pipelines
- Share common workflow modules across different pipelines
- Make your code more maintainable and easier to understand

In this side quest, we'll create a pipeline that demonstrates these benefits by:
1. Creating independent workflow modules that can be tested and used separately
2. Composing these modules into a larger pipeline
3. Using Nextflow's workflow composition features to manage data flow between modules

---
## 0. Warmup

### 0.1 Prerequisites

Before taking on this side quest you should:

- Complete the [Hello Nextflow](../hello_nextflow/README.md) tutorial
- Understand basic Nextflow concepts (processes, channels, operators)

### 0.1. Starting Point

Let's move into the project directory.

```bash
cd side-quests/workflows_of_workflows
```

You'll find a `modules` directory containing several process definitions that build upon what you learned in 'Hello Nextflow':

```console title="Directory contents"
modules/
├── say_hello.nf             # Creates a greeting (from Hello Nextflow)
├── say_hello_upper.nf       # Converts to uppercase (from Hello Nextflow)
├── timestamp_greeting.nf    # Adds timestamps to greetings
├── validate_name.nf         # Validates input names
└── reverse_text.nf          # Reverses text content
```

We going to compose these modules into two separate workflows that we will then compose into a main workflow.

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

This is already a complete workflow that we can test independently. You might notice that it's very similar to the `hello` workflow from the 'Hello Nextflow' tutorial, but with some additional syntax to allow it to receive input ('take:') and produce output ('emit:'). These are the connections we will use to compose a higher level workflow.

### 1.3. Create and test the main workflow

Now we will create a main workflow that imports and uses the `greeting` workflow.

Create `main.nf`:

```groovy title="main.nf" linenums="1"
include { GREETING_WORKFLOW } from './workflows/greeting'

workflow {
    names = Channel.from('Alice', 'Bob', 'Charlie')
    GREETING_WORKFLOW(names)

    GREETING_WORKFLOW.out.greetings.view { "Original: $it" }
    GREETING_WORKFLOW.out.timestamped.view { "Timestamped: $it" }
}

```

Run this and see the output:
```bash title="Run the workflow"
nextflow run main.nf
```

```console title="Expected output"
N E X T F L O W  ~  version 23.10.1
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

### Takeaway

You should now have a working greeting workflow that:
- Takes a channel of names as input
- Validates each name
- Creates a greeting for each valid name
- Adds timestamps to the greetings

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

### 2.3. Update the main workflow

Update `main.nf` to use both workflows:

```groovy title="main.nf" linenums="1"
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

Run the complete pipeline:
```bash title="Run the workflow"
nextflow run main.nf
```

```console title="Expected output"
N E X T F L O W  ~  version 23.10.1
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
cat /workspaces/training/side_quests/workflows_of_workflows/work/f0/74ba4a10d9ef5c82f829d1c154d0f6/REVERSED-UPPER-timestamped_Alice-output.txt```

```console title="Reversed file content"
!ECILA ,OLLEH ]04:50:71 60-30-5202[
```

### Takeaway

You should now have a complete pipeline that:
- Processes names through the greeting workflow
- Feeds the timestamped greetings into the transform workflow
- Produces both uppercase and reversed versions of the greetings

## Summary

In this side quest, we've explored the powerful concept of workflow composition in Nextflow, which allows us to build complex pipelines from smaller, reusable components. Here's what we've accomplished:

1. **Created Modular Workflows**: We built two independent workflow modules:
   - A `GREETING_WORKFLOW` that validates names, creates greetings, and adds timestamps
   - A `TRANSFORM_WORKFLOW` that converts text to uppercase and reverses it

2. **Composed Workflows Together**: We connected these workflows in a main pipeline, demonstrating how data can flow from one workflow to another.

3. **Used Workflow Interfaces**: We defined clear inputs and outputs for each workflow using the `take:` and `emit:` syntax, creating well-defined interfaces between components.

4. **Managed Data Flow**: We learned how to access workflow outputs using the namespace notation (`WORKFLOW_NAME.out.channel_name`) and pass them to other workflows.

5. **Practiced Modular Design**: We experienced firsthand how breaking a pipeline into logical components makes the code more maintainable and easier to understand.

This modular approach offers several advantages over monolithic pipelines:
- Each workflow can be developed, tested, and debugged independently
- Workflows can be reused across different pipelines
- The overall pipeline structure becomes more readable and maintainable
- Changes to one workflow don't necessarily affect others if the interfaces remain consistent

It's important to note that while calling workflows is a bit like calling processes, it's not the same. You can't, for example, run a workflow n times by calling it with a channel of size n- you would need to pass a channel of size n to the workflow and iterate internally.

By mastering workflow composition, you're now equipped to build more sophisticated Nextflow pipelines that can handle complex bioinformatics tasks while remaining maintainable and scalable.

### Key Concepts

1. **Workflow Inclusion**
   ```nextflow
   // Include a single workflow
   include { WORKFLOW_NAME } from './path/to/workflow'

   // Include multiple workflows
   include { WORKFLOW_A; WORKFLOW_B } from './path/to/workflows'

   // Include with alias to avoid name conflicts
   include { WORKFLOW_A as WORKFLOW_A_ALIAS } from './path/to/workflow'
   ```

2. **Workflow Inputs and Outputs**
   ```nextflow
   workflow EXAMPLE {
       take:
           input_ch    // Declare inputs
       emit:
           output_ch   // Declare outputs
   }
   ```

3. **Workflow Composition**
   ```nextflow
   // Using explicit connections
   WORKFLOW_A(input_ch)
   WORKFLOW_B(WORKFLOW_A.out.some_channel)
   ```

## Resources

- [Nextflow Workflow Documentation](https://www.nextflow.io/docs/latest/workflow.html)
- [Channel Operators Reference](https://www.nextflow.io/docs/latest/operator.html)
- [Error Strategy Documentation](https://www.nextflow.io/docs/latest/process.html#errorstrategy)
