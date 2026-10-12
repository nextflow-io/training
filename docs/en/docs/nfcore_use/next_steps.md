# Course summary

Congratulations on completing the Use nf-core training course! 🎉

<!-- placeholder for video -->

## Your journey

You started by finding and retrieving the `nf-core/demo` pipeline, then learned to run it using its test profile and examine its outputs.
Next, you configured its execution through pipeline parameters and configuration files, and saw how nf-core pipelines validate parameters and input data.
Finally, you applied those same skills to `nf-core/rnaseq`, a production-scale pipeline, and learned how to override its default resource allocations to fit the hardware available to you.
You then chose a pipeline for a scientific question by evaluating its documentation and a completed run, prepared a samplesheet for your own data, diagnosed and fixed samplesheet errors, and used execution reports to tune resource allocations.

### What you learned

You are now able to find, retrieve, run, and configure nf-core pipelines.

- nf-core pipelines are retrieved with `nextflow pull` and follow a standard code organization.
- Every nf-core pipeline ships with a `test` profile for quick validation on a small dataset.
- Pipeline parameters (set via `--param_name` or `-params-file`) and configuration (set via `-c`) serve different purposes: inputs and analysis options versus execution logistics like resource allocation.
- nf-core pipelines validate parameters and input files automatically, catching errors before any work is done.
- Resource defaults are assigned through labels (`process_low`, `process_medium`, `process_high`) defined in `conf/base.config`, which you can override with a custom configuration file.
- Samplesheet requirements are documented in the usage docs and enforced by the input schema (`assets/schema_input.json`); some pipelines add their own checks across rows.
- `-preview` and `-stub` let you check inputs without running the full analysis.
- Execution reports show how much of its allocation each process actually used, which tells you how to right-size resources.

### Skills acquired

Through this hands-on course, you've learned how to:

- Find an nf-core pipeline on the nf-co.re website and retrieve its source code
- Run a pipeline using its built-in test profile and examine its outputs
- Get help, set parameters, and understand parameter and input validation
- Customize resource allocation and tool arguments through configuration files
- Pull and run a production-scale pipeline, and override its default resource labels
- Find and evaluate pipelines for a specific scientific question
- Write a samplesheet for your own data, including samples sequenced on multiple lanes
- Diagnose and fix input validation errors
- Tune resource allocations based on the execution report

You're now equipped with the foundational knowledge to start running nf-core pipelines for your own analyses.

## Next steps to build your skills

Here are our top suggestions for what to do next:

- Launch and monitor these pipelines at scale with [Scale with Seqera](../seqera_scale/index.md)
- Don't just run nf-core pipelines, develop them! Learn nf-core best practices with [Build with nf-core](../nfcore_build/index.md)
- New to Nextflow itself? Start with [Nextflow Run](../nextflow_run/index.md)
- Apply Nextflow to a scientific analysis use case with [Nextflow for Science](../nf4_science/index.md)
- Explore more advanced Nextflow features with the [Side Quests](../side_quests/index.md)

## Getting help

For help resources and community support, see the [Help page](../help.md).

## Feedback survey

Before you move on, please take a minute to complete the course survey! Your feedback helps us improve our training materials for everyone.

[Take the survey :material-arrow-right:](survey.md){ .md-button .md-button--primary }
