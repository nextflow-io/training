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
It takes real actions on your behalf, using your authenticated Platform access and your connected GitHub account.
It can list compute environments, launch a run, read run logs, and open a pull request.

!!! info "MCP and the Seqera Platform"

    CoScientist reaches the Platform through the **Seqera MCP server** (Model Context Protocol).
    MCP is the interface that exposes Platform operations as tools the agent can call: list compute environments, launch a run, read run logs, query datasets.
    The agent calls these tools with your own Platform token, so it acts with your permissions rather than a shared account.
    The same MCP server is available to other assistants, and in the CLI you can reach it directly through the `/seqera-mcp` skill.

## 3. Running on your own foundational model (enterprise).

On Seqera Enterprise deployments, the model that powers CoScientist can run inside the customer's own cloud account.

!!! info

    Seqera Enterprise CoScientist is currently available on AWS only.
    It can run Claude inference inside the customer's own AWS account through Amazon Bedrock, so prompts and data stay within that account.
    Seqera Cloud, including this training environment, uses Seqera's hosted model and meters usage through credits.

!!! note "Checkpoint"

    You can state, in one sentence, that CoScientist takes real actions on your Platform and GitHub assets rather than only producing text.

### Takeaway

CoScientist is an agent that takes real actions on your Platform and GitHub assets.
Responses that modify state (launching a run, opening a PR, registering a pipeline) are the result of real tool calls, not generated text alone.

### What's next?

In the next lesson, [put CoScientist to work on real Platform assets](02_explore_platform.md).
