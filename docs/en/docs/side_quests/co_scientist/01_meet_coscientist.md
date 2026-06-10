# Meet CoScientist

This is your first real conversation with the agent.
CoScientist does not just answer questions; it takes real actions on your Platform and GitHub on your behalf.
This lesson shows what it can do, and introduces the Seqera MCP Server that exposes the Platform to AI tools.

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
It takes real actions on your behalf, using your authenticated access to the workspace and your connected GitHub account.
It can list compute environments, launch a run, read run logs, and open a pull request.

!!! info "MCP and the Seqera Platform"

    Seqera exposes Platform operations to AI assistants through the **Seqera MCP Server** (Model Context Protocol).
    MCP is an open standard that lets an assistant call real tools rather than only describing them.
    It is what lets assistants such as Claude Code reach your Platform, and you can connect it to your own tools as well.

<!-- TODO: verify how CoScientist itself accesses the Platform internally. The Seqera MCP Server is documented as the interface for external AI assistants; the CoScientist docs do not state CoScientist uses it as its own transport. Confirm with engineering before asserting a mechanism. -->

## 3. Running on your own foundational model (enterprise).

On enterprise deployments, the foundational model can be swapped for one inside the customer's own cloud account.

!!! info

    On Seqera Enterprise (AWS only), CoScientist can run against a customer's own foundational model via Amazon Bedrock, so prompts and data stay within the customer's AWS account.
    The training environment uses the default hosted model.

<!-- TODO: verify the Bedrock / AWS-only enterprise foundational-model detail against current Seqera product docs before publishing -->

!!! note "Checkpoint"

    You can state, in one sentence, that CoScientist takes real actions on your Platform and GitHub assets rather than only producing text.

### Takeaway

CoScientist is an agent that takes real actions on your Platform and GitHub assets.
Responses that modify state (launching a run, opening a PR, registering a pipeline) are the result of real tool calls, not generated text alone.

### What's next?

In the next lesson, [put CoScientist to work on real Platform assets](02_explore_platform.md).
