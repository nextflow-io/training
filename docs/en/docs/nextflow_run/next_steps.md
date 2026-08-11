# Course summary

Congratulations on completing the Nextflow Run training course! 🎉

<!-- placeholder for video -->

## Your journey

You started with a very basic workflow, and learned to run it, find the outputs, and manage its execution.
Then, you worked your way through increasingly more complex versions of that workflow and learned to recognize the essential concepts and mechanisms that power Nextflow pipelines, including channels, containers, and `-resume`.
Then, you learned how to customize the configuration of a pipeline using `nextflow.config` and profiles.
Finally, you learned how to run pipelines directly from remote repositories such as GitHub.

### What you learned

You are now able to manage the execution of the Hello pipeline, describe how it is structured, and identify the main pieces of code involved.

- The final form of the Hello workflow takes as input a CSV file containing text greetings.
- The four steps are implemented as Nextflow processes (`sayHello`, `convertToUpper`, `collectGreetings`, and `cowpy`) stored in separate module files.
- The results are published to a directory called `results/`.
- The final output of the pipeline is a plain text file containing ASCII art of a character saying the transformed text.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello_pipeline_complete.svg"
</figure>

1. **`sayHello`:** Writes each greeting to its own output file (_e.g._ "Hello-output.txt")
2. **`convertToUpper`:** Converts each greeting to uppercase (_e.g._ "HELLO")
3. **`collectGreetings`:** Collects all uppercase greetings into a single batch file
4. **`cowpy`:** Generates ASCII art using the `cowpy` tool

The workflow configuration supports providing inputs and parameters in a flexible, reproducible way.

### Skills acquired

Through this hands-on course, you've learned how to:

- Launch a Nextflow workflow locally and find its outputs
- Recognize the core Nextflow components that constitute a simple multi-step workflow
- Process multiple inputs in parallel using channels
- Use `-resume` to avoid repeating completed work
- Use containers to manage software dependencies
- Configure pipelines using `nextflow.config` and profiles
- Run a pipeline directly from a remote repository and pin it to a specific revision

You're now equipped with the foundational knowledge to start integrating existing Nextflow pipelines into your own work.

## Next steps to build your skills

Here are our top suggestions for what to do next:

- Dive deeper into pipeline configuration with [Nextflow Dot Config](../nextflow_config/index.md)
- Learn to run nf-core community pipelines with [Run nf-core](../nfcore_run/index.md)
- Launch and monitor pipelines at scale with [Run with Seqera](../seqera_run/index.md)
- Don't just run Nextflow, write it! Become a Nextflow developer with [Hello Nextflow](../hello_nextflow/index.md)
- Apply Nextflow to a scientific analysis use case with [Nextflow for Science](../nf4_science/index.md)
- Learn troubleshooting techniques with the [Debugging Side Quest](../side_quests/debugging/index.md)

## Getting help

For help resources and community support, see the [Help page](../help.md).

## Feedback survey

Before you move on, please take a minute to complete the course survey! Your feedback helps us improve our training materials for everyone.

[Take the survey :material-arrow-right:](survey.md){ .md-button .md-button--primary }
