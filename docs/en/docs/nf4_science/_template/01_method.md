# Part 1: Method overview and manual testing

{BRIEF_METHOD_DESCRIPTION}

![Pipeline overview](img/{PIPELINE_DIAGRAM}.png)

{OPTIONAL_EXTENDED_METHOD_DESCRIPTION}

### Methods

{DESCRIBE_THE_TWO_APPROACHES: single-sample and multi-sample/aggregation}

Before we dive into writing any workflow code, we are going to try out the commands manually on some test data.

### Dataset

We provide the following data and related resources:

- **{PRIMARY_INPUT_DESCRIPTION}**
- **{SAMPLE_DESCRIPTION}**
- **{ADDITIONAL_RESOURCES_DESCRIPTION}**

### Software

The main tools involved are [{TOOL_A}]({TOOL_A_URL}) and [{TOOL_B}]({TOOL_B_URL}).

These tools are not installed in the GitHub Codespaces environment, so we'll use them via containers (see [Hello Containers](../../hello_nextflow/05_hello_containers.md)).

!!! note

    Make sure you're in the `nf4-science/{DOMAIN_DIR}` directory so that the last part of the path shown when you type `pwd` is `{DOMAIN_DIR}`.

---

## 1. {SINGLE_SAMPLE_SECTION_TITLE}

{BRIEF_DESCRIPTION_OF_SINGLE_SAMPLE_APPROACH}

In this section we test the commands that make up the single-sample processing approach.
These are the commands we'll wrap into a Nextflow workflow in Part 2 of this course.

1. {STEP_1_SUMMARY}
2. {STEP_2_SUMMARY}

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/{DOMAIN_DIR}/img/{DIAGRAM_1}.svg"
</figure>

We start by testing the commands on just one sample.

### 1.1. {FIRST_TOOL_STEP_TITLE}

{BRIEF_DESCRIPTION_OF_STEP}

#### 1.1.1. Pull the container

Run the `docker pull` command to download the container image:

```bash
docker pull {TOOL_A_CONTAINER_URI}
```

??? success "Command output"

    ```console
    {EXPECTED_PULL_OUTPUT}
    ```

#### 1.1.2. Spin up the container interactively

Spin up the container and mount the `data` directory so the tools can access the input files:

```bash
docker run -it -v ./data:/data {TOOL_A_CONTAINER_URI}
```

Your prompt changes to indicate you are inside the container.

#### 1.1.3. Run the command

```bash
{TOOL_A_COMMAND}
```

??? success "Command output"

    ```console
    {EXPECTED_COMMAND_OUTPUT}
    ```

{DESCRIBE_EXPECTED_OUTPUT_FILES}

#### 1.1.4. Exit the container

To exit the container, type `exit`.

```bash
exit
```

Your prompt should be back to normal.

### 1.2. {SECOND_TOOL_STEP_TITLE}

{BRIEF_DESCRIPTION_OF_STEP}

#### 1.2.1. Pull the container

```bash
docker pull {TOOL_B_CONTAINER_URI}
```

??? success "Command output"

    ```console
    {EXPECTED_PULL_OUTPUT}
    ```

#### 1.2.2. Spin up the container interactively

```bash
docker run -it -v ./data:/data {TOOL_B_CONTAINER_URI}
```

#### 1.2.3. Run the command

```bash
{TOOL_B_COMMAND}
```

??? success "Command output"

    ```console
    {EXPECTED_COMMAND_OUTPUT}
    ```

{DESCRIBE_EXPECTED_OUTPUT_FILES}

#### 1.2.4. Exit the container

To exit the container, type `exit`.

```bash
exit
```

Your prompt should be back to normal.
That concludes the single-sample processing test.

---

## 2. {MULTI_SAMPLE_SECTION_TITLE}

{BRIEF_DESCRIPTION_OF_MULTI_SAMPLE_APPROACH}

{EXPLAIN_WHY_MULTI_SAMPLE_IS_NEEDED}

In this section we test the additional commands needed for multi-sample processing.
These are the commands we'll wrap into a Nextflow workflow in Part 3 of this course.

1. {AGGREGATION_STEP_1_SUMMARY}
2. {AGGREGATION_STEP_2_SUMMARY}

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/{DOMAIN_DIR}/img/{DIAGRAM_2}.svg"
</figure>

### 2.1. {AGGREGATION_STEP_TITLE}

{STEP_INSTRUCTIONS_FOLLOWING_SAME_PATTERN_AS_ABOVE}

---

### Takeaway

You know how to test the {TOOL_A} and {TOOL_B} commands in their respective containers, including how to {MULTI_SAMPLE_SUMMARY}.

### What's next?

Learn how to wrap those same commands into workflows that use containers to execute the work.
