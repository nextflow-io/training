---
name: write-nf-test
description: Generate an nf-test that asserts on stable output and excludes unstable content.
---

# Write nf-test

When asked to write an nf-test for a process or a pipeline, follow the rules below.

## Steps

1. Scaffold the test under `tests/`.
   For a process use a `nextflow_process` block; for a pipeline use a `nextflow_pipeline` block.
   Provide representative inputs or test data so the test runs end to end.

2. Assert on deterministic output only:

   - Numeric or tabular result files, such as a Salmon `quant.sf` table
   - The existence of the expected output files and directories

3. Do NOT snapshot unstable content:

   - Reports that embed timestamps or versions, such as MultiQC HTML
   - Log files, and files that record the command line or tool version (for example Salmon `cmd_info.json` and `meta_info.json`)
   - Any path or value containing the work-directory hash
   - Version strings

4. Run `nf-test test` on the generated file.
   When it fails, narrow the assertion rather than snapshotting unstable content, and re-run until it passes.
