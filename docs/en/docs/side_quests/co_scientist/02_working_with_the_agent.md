# Working with the agent

CoScientist takes real actions on your workspace, so how you prompt and check it matters.
The next lessons let it change real code, where a careless prompt can do work you then have to undo.
This lesson covers the habits that keep you in control.

---

## 1. Prompt for intent.

Describe the goal and any constraints, then let the agent work out the steps.
An underspecified prompt forces it to guess, and the guess is often wrong.

=== "Better"

    ```text
    The QUANT process fails with a non-zero exit.
    Read the run log, tell me the cause, and propose a fix before changing anything.
    ```

=== "Vague"

    ```text
    Fix the pipeline.
    ```

The better prompt gives the agent the symptom, points it at the run log, and asks it to explain before acting.
That keeps the change small and lets you see the reasoning before anything happens.

## 2. Verify what the agent did.

The agent often reports success even when it is wrong, so check the real state yourself instead of trusting its summary.

Things to check after an action:

- The diff: read the actual change, not the agent's description of it.
- The Platform state: confirm a run reached the expected status, or a pipeline appears on the Launchpad.
- The artifact: confirm the pull request exists, targets the right branch, and contains what you expect.

<!-- TODO: verify how CoScientist surfaces diffs/changes in the UI -->

!!! warning

    Treat the agent's "Done" as a claim to check, not a fact.
    The Checkpoints in this side quest exist for this reason: each one is a real Platform or GitHub state you confirm yourself, whatever the agent says it did.

## 3. Redirect the agent when it goes wrong.

When the agent proposes something wrong, too broad, or based on a false assumption, correct it straight away instead of accepting the change and cleaning up later.
Telling it exactly what to touch, and asking it to pause before it commits, keeps you in control of an agent that acts quickly.

```text
That changes more than I asked for. Only modify the QUANT process and leave the rest untouched.
```

```text
Before you commit, show me the diff and wait for my confirmation.
```

```text
You assumed the compute environment is AWS Batch. Confirm that before proceeding.
```

## 4. Give the agent the real source of truth.

Hand the agent the failing log, the actual error message, or the exact file path rather than your paraphrase of them.
It reasons far better from the real artifact than from a description of it.
Because CoScientist can reach your Platform and repository assets directly, point it at the artifact instead of summarizing what you saw.

<!-- TODO: verify CoScientist can read raw run logs directly, not just Platform metadata -->

!!! note "Checkpoint"

    You can describe how to stay in control of the agent: ask for the outcome you want, give it the real inputs, and check what it actually did before moving on.

### Takeaway

You ask for an outcome, the agent acts, and you check the result.
You stay responsible for whether the work is correct; the agent makes it faster.

### What's next?

[Start developing against rnaseq-nf](03_development.md).
