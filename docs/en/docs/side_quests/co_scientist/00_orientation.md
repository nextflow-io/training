# Orientation

Before starting the exercises, open the CoScientist chat interface and confirm it can see the training workspace.

---

## 1. Open the chat interface.

Navigate to [https://ai.seqera.io/chat](https://ai.seqera.io/chat) and sign in with your Seqera credentials.

<!-- TODO: screenshot: the CoScientist chat landing view after login -->

## 2. Connect to the training workspace.

Use the workspace selector to choose the provided training workspace.

<!-- TODO: verify exact label of the workspace selector control in the UI -->

<!-- TODO: screenshot: workspace selector showing the training workspace -->

## 3. Confirm access with a first prompt.

Copy and send the following prompt:

```text
Which Seqera workspace am I connected to, and what compute environments are available here?
```

??? example "What CoScientist typically does"

    It names the connected workspace and lists the available compute environment(s).
    The exact wording will differ from run to run.

!!! note "Checkpoint"

    CoScientist correctly names the **training workspace** and lists at least one compute environment.
    If it cannot, the workspace connection is not set up; revisit step 2.

### Takeaway

You have a working CoScientist chat session bound to the training workspace and are ready to interact with its Platform assets.

### What's next?

In the next lesson, you will meet the agent and learn how it reaches Platform assets.
