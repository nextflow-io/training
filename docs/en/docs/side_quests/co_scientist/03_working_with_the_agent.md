# Working with the agent

CoScientist is an agent that takes real actions, not a search box.
Working well with it is a skill: state your intent clearly, verify what it did, and steer it when it drifts.
The next lessons let it change real code, so these habits matter now.

---

## 1. Prompt for intent, not for keystrokes.

Describe the goal and the constraints, then let the agent choose the steps.
You are not writing commands; you are delegating an outcome.
An underspecified prompt forces the agent to guess, and the guess is often wrong.

=== "Better"

    ```text
    The QUANT process fails with a non-zero exit.
    Read the run log, tell me the cause, and propose a fix before changing anything.
    ```

=== "Vague"

    ```text
    Fix the pipeline.
    ```

The better prompt works because it gives the agent the symptom, names the source of truth (the run log), and asks it to explain before acting.
That limits scope and keeps you in the loop before any change is made.

## 2. Verify what the agent did.

The agent reports success confidently even when it is wrong, so confirm the real-world state yourself rather than trusting the summary.
An agent's natural tendency is to sound certain; your job is to check that certainty against actual artifacts.

Things to verify after an action:

- The diff: read the actual change, not the agent's description of it.
- The Platform state: confirm a run reached the expected status, or a pipeline appears on the Launchpad.
- The artifact: confirm the PR exists and targets the right branch, and that the committed file contains what you expect.

<!-- TODO: verify how CoScientist surfaces diffs/changes in the UI -->

!!! warning

    Treat the agent's "Done" as a claim to check, not a fact.
    The Checkpoints in this side quest exist for exactly this reason: each one is a real Platform or GitHub state you confirm yourself, independent of what the agent says it did.

## 3. Course-correct when it drifts.

When the agent proposes something wrong, over-broad, or based on a false assumption, redirect it immediately rather than accepting the change and cleaning up afterward.
Corrections are cheap before a commit; they are expensive after one.
Giving the agent explicit scope and an explicit pause point before it acts keeps you in control of something that moves quickly.

```text
That changes more than I asked for. Only modify the QUANT process and leave the rest untouched.
```

```text
Before you commit, show me the diff and wait for my confirmation.
```

```text
You assumed the compute environment is AWS Batch. Confirm that before proceeding.
```

## 4. Keep the agent grounded.

Give it the real source of truth (the failing log, the actual error message, the exact file path) rather than letting it reason from your paraphrase.
An agent working from the actual run log is far more reliable than one working from your description of it.
Because CoScientist can read logs and configs directly over MCP (as covered in the previous lesson), prefer pointing it at the real artifact over summarizing what you saw.

!!! note "Checkpoint"

    You can describe the three habits you will use in the rest of this side quest: prompt for intent, verify the agent's actions against real state, and course-correct when it drifts.

### Takeaway

Working with an agent is a loop of intent, action, and verification.
You stay responsible for correctness, and the agent accelerates the work.

### What's next?

[Start developing against rnaseq-nf](04_development.md).
