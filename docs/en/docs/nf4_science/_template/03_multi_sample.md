# Part 3: Multi-sample aggregation

In Part 2, you built a per-sample processing pipeline that handled each sample independently.
Now we're going to extend it to implement multi-sample {AGGREGATION_METHOD}, as covered in [Part 1](01_method.md).

## Assignment

In this part of the course, we're going to extend the workflow to do the following:

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/{DOMAIN_DIR}/img/{DIAGRAM_2}.svg"
</figure>

1. {PER_SAMPLE_STEP_1}
2. {PER_SAMPLE_STEP_2}
3. {AGGREGATION_STEP_1}
4. {AGGREGATION_STEP_2}

This part builds directly on the workflow produced by Part 2.

??? info "How to begin from this section"

    This section of the course assumes you have completed [Part 2: Single-sample processing](./02_single_sample.md) and have a working `{DOMAIN_DIR}.nf` pipeline.

    If you did not complete Part 2 or want to start fresh for this part, you can use the Part 2 solution as your starting point.
    Run these commands from inside the `nf4-science/{DOMAIN_DIR}/` directory:

    ```bash
    cp solutions/part2/{DOMAIN_DIR}-2.nf {DOMAIN_DIR}.nf
    cp solutions/part2/nextflow.config .
    cp solutions/part2/modules/* modules/
    ```

    This gives you a complete single-sample processing workflow.
    You can test that it runs successfully by running the following command:

    ```bash
    nextflow run {DOMAIN_DIR}.nf -profile test
    ```

## Lesson plan

We've broken this down into two steps:

1. **{MODIFICATION_STEP_SUMMARY}.**
   This covers updating process commands and outputs.
2. **{AGGREGATION_STEP_SUMMARY}.**
   This introduces the `collect()` operator {AND_OTHER_CONCEPTS}.

!!! note

     Make sure you're in the correct working directory:
     `cd /workspaces/training/nf4-science/{DOMAIN_DIR}`

---

## 1. {MODIFICATION_STEP_TITLE}

{DESCRIPTION_OF_MODIFICATION_TO_EXISTING_PROCESSES}

Recall the modified command from [Part 1](01_method.md):

```bash
{MODIFIED_COMMAND}
```

{EXPLAIN_DIFFERENCES_FROM_PART_2}

### 1.1. {MODIFICATION_SUBSTEP}

{INSTRUCTIONS_WITH_BEFORE_AFTER_TABS}

=== "After"

    ```groovy title="modules/{MODULE_FILE}.nf" linenums="{N}" hl_lines="{LINES}"
    {UPDATED_CODE}
    ```

=== "Before"

    ```groovy title="modules/{MODULE_FILE}.nf" linenums="{N}" hl_lines="{LINES}"
    {ORIGINAL_CODE}
    ```

### 1.2. Run the workflow to verify the modification

```bash
nextflow run {DOMAIN_DIR}.nf -profile test
```

??? success "Command output"

    ```console
    {EXPECTED_OUTPUT}
    ```

{VERIFY_OUTPUT_INSTRUCTIONS}

### Takeaway

You know how to modify process commands and outputs to adapt the workflow behavior.

### What's next?

Add the multi-sample aggregation step.

---

## 2. {AGGREGATION_STEP_TITLE}

{DESCRIPTION_OF_AGGREGATION}

### 2.1. Write the aggregation module

{MODULE_INSTRUCTIONS}

### 2.2. Collect per-sample outputs and feed them to the aggregation process

{INSTRUCTIONS_USING_COLLECT_OPERATOR}

### 2.3. Run the completed workflow

```bash
nextflow run {DOMAIN_DIR}.nf -profile test
```

??? success "Command output"

    ```console
    {EXPECTED_OUTPUT}
    ```

{VERIFY_OUTPUT_INSTRUCTIONS}

---

### Takeaway

You have a complete pipeline that processes samples individually and aggregates results across all samples.
You know how to use channel operators like `collect()` to aggregate per-sample outputs for multi-sample analysis.

### What's next?

Congratulations on completing this course! Head to the [course summary](next_steps.md) to review what you've learned and explore next steps.
