---
name: write-nf-test
description: Generate an nf-test for a Nextflow process that asserts on stable output and excludes unstable content.
---

# Write nf-test

When this skill is invoked with a process name, generate an `nf-test` for that process following the rules below.

## Steps

1. Scaffold an `nf-test` for the named process under `tests/`.
   Use a `nextflow_process` block with the process `name`, `script`, and `process` fields, and provide a representative `input` in the `when` block.

2. Assert on deterministic output only:

   - The per-transcript count columns in `quant.sf`
   - The existence of the expected output files and directories

3. Do NOT snapshot unstable content:

   - MultiQC HTML reports (embedded timestamps and versions)
   - Salmon `cmd_info.json`, `meta_info.json`, and log files (timestamps, command line, version)
   - Any path or value containing the work-directory hash
   - Version strings

4. Run `nf-test test` on the generated file.
   When it fails, narrow the assertion rather than snapshotting unstable content, and re-run until it passes.
