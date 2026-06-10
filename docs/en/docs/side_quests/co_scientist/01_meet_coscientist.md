# Meet CoScientist

This is your first real conversation with the agent.
The one concept that explains everything else is how CoScientist reaches your tools. That mechanism is covered in section 2.

---

## 1. Ask CoScientist what it can do.

Send the following prompt to get a baseline picture of the agent's capabilities in this workspace:

```text
What can you help me with in this workspace? What can you see and what actions can you take?
```

??? example "What CoScientist typically does"

    It describes that it can help develop pipelines, launch and monitor runs, browse data, and act on GitHub.
    The exact wording will differ from run to run.

## 2. How CoScientist reaches your assets.

CoScientist does not generate text and hand it back for you to act on.
It calls real tools on your behalf.
The mechanism that makes this possible is MCP.

!!! info

    CoScientist talks to the Seqera Platform and GitHub through **MCP** (Model Context Protocol).
    MCP is the standard interface that lets the agent call real tools: list compute environments, launch a run, read run logs, open a pull request, rather than only generating text.
    When you ask it to register a pipeline on the Launchpad, it is calling a Platform tool over MCP on your behalf.

## 3. Running on your own foundational model (enterprise).

On enterprise deployments, the foundational model can be swapped for one inside the customer's own cloud account.

!!! info

    On Seqera Enterprise (AWS only), CoScientist can run against a customer's own foundational model via Amazon Bedrock, so prompts and data stay within the customer's AWS account.
    The training environment uses the default hosted model.

<!-- TODO: verify the Bedrock / AWS-only enterprise foundational-model detail against current Seqera product docs before publishing -->

!!! note "Checkpoint"

    You can state, in one sentence, that CoScientist acts on the Platform through MCP tool calls rather than only producing text.

### Takeaway

CoScientist is an agent that takes real actions on your Platform and GitHub assets via MCP.
Responses that modify state (launching a run, opening a PR, registering a pipeline) are the result of tool calls, not generated text alone.

### What's next?

In the next lesson, [put MCP to work on real Platform assets](02_explore_platform.md).
