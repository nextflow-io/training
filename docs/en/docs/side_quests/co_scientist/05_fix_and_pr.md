# Fixing and proposing a change

Apply CoScientist's fix to your fork, confirm the pipeline runs to completion, and open a pull request.

---

## 1. Apply the fix.

Ask CoScientist to correct the configuration it identified as the cause of the failure:

```text
Fix the cause of the failure in my fork and commit the change.
```

CoScientist updates the memory directive for the `QUANT` process and commits the change to your fork.

=== "After"

    ```groovy title="nextflow.config" hl_lines="2"
    process {
        withName: QUANT { memory = 2.GB }
    }
    ```

=== "Before"

    ```groovy title="nextflow.config" hl_lines="2"
    process {
        withName: QUANT { memory = 1.MB }
    }
    ```

<!-- TODO: confirm the exact config directive matches the break recorded in 03_development.md -->

!!! note "Checkpoint"

    A new commit with the fix is on your fork.

## 2. Re-run to green.

Ask CoScientist to submit a new run with the corrected configuration:

```text
Re-launch the pipeline and tell me when it completes.
```

CoScientist launches the pipeline on the Launchpad and monitors its progress, reporting back when the run finishes.

!!! note "Checkpoint"

    A new run reaches status `succeeded` in the Platform **Runs** list.

<!-- TODO: screenshot: Runs list with a succeeded run -->

## 3. Open a pull request.

Ask CoScientist to propose the fix upstream:

```text
Open a pull request from my fork back to nextflow-io/rnaseq-nf describing the fix.
```

!!! note

    In a real contribution the pull request targets the upstream repository.
    In a training setting, if upstream PRs are restricted, it is fine to open the PR against your own fork's default branch.
    The trainer can adjust the target as appropriate for the environment.

!!! note "Checkpoint"

    A pull request is open and visible on GitHub.

<!-- TODO: screenshot: the open PR on GitHub -->

### Takeaway

CoScientist diagnosed the failure, committed a fix, verified the pipeline completed successfully, and proposed the change as a pull request — the full loop from failure to contribution.

### What's next?

[Add a test to the pipeline](06_feature_nf_test.md).
