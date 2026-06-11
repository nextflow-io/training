# Debugging a failed run

The run you launched in the previous lesson has finished, and it failed.
Use CoScientist to find out why.

---

## 1. Find the failed run.

In the Platform **Runs** list, locate the run for your fork.
It now shows a `failed` status.

!!! note "Checkpoint"

    A run with status `failed` is visible in the Platform **Runs** list.

<!-- TODO: screenshot: Runs list with a failed rnaseq-nf run -->

## 2. Link from the Platform back to CoScientist.

From the failed run's detail page, open the action that navigates back to a CoScientist conversation about that run.
This takes you directly into a CoScientist session with context about the run already loaded.

<!-- TODO: verify the exact UI affordance that links a run to CoScientist -->
<!-- TODO: screenshot: the link/button from a run to CoScientist -->

## 3. Ask CoScientist what went wrong.

With the run in context, ask CoScientist to explain the failure:

```text
This run failed. What caused the failure, and which process and exit code are involved?
```

??? example "What CoScientist typically does"

    It reads the run logs, identifies the `QUANT` process as the point of failure, reports that the task was killed for exceeding its memory limit (exit code 137), and points to the `withName: QUANT` memory setting as the cause.
    The exact wording will differ from run to run.

## 4. What just happened.

CoScientist read the failing task's `.command.log` and exit code and reasoned over them.
This is the same information you would read by hand from the work directory.
This is the same workflow described in [Troubleshooting Workflows](../debugging/index.md), performed through a conversational interface instead of manual inspection.
The failure you set up intentionally in the previous lesson, constraining the `QUANT` process to an unrealistically small memory allocation, produced exactly the kind of exit code and log output that CoScientist can interpret.

!!! note "Checkpoint"

    CoScientist has identified the `QUANT` process and a memory-related cause for the failure.

### Takeaway

CoScientist reads a failed run's logs and diagnoses the cause, turning a red `failed` status into an actionable explanation without requiring manual navigation of the work directory.

### What's next?

[Apply the fix and open a pull request](05_fix_and_pr.md).
