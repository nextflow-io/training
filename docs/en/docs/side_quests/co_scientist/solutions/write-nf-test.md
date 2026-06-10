---
name: write-nf-test
description: Generate an nf-test for a Nextflow process that asserts on stable output and excludes unstable content.
---

<!-- NOTE: illustrative skill file; verify CoScientist's actual skill format/location before use -->

# Write nf-test

When this skill is invoked with a process name, generate an `nf-test` for that process following the rules below.

## Steps

1. Scaffold an `nf-test` for the named process in the `tests/` directory.
   Use the `nextflow_process` block with `name`, `script`, and `process` fields.
   Provide a representative `input` in the `when` block.

2. Assert on deterministic output only:

   - `Salmon` per-transcript count columns (stable given fixed inputs)
   - `FASTQC` pass/fail status (stable for fixed data)
   - File existence and line counts (structural, stable)

3. Do NOT snapshot unstable content:

   - `MultiQC` HTML reports (embed timestamps and software versions)
   - `Salmon` `cmd_info.json` and log files (contain timestamps and absolute paths)
   - Any path or value containing the work-directory hash (changes every run)
   - Version strings (change with tool and container updates)

4. Run `nf-test test` on the generated file.
   If the test fails, read the failure message and narrow the assertion rather than updating the snapshot when the cause is unstable content.
   Repeat until the test reports a pass.
