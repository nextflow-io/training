# terraform/pipeline

The Maintain tier's pipeline asset. This adds your own pipeline,
`rnaseq-nf-$USERNAME`, to the Launchpad with Terraform. It adds it. It does not
launch it. See the Platform Automation side quest, §2 (Maintain).

## Before you start

You need the **Maintain** role and two values from the Admin tier (or your
presenter):

- `workspace_id`: the shared workspace.
- `compute_env_id`: the shared Azure Batch compute environment.

And create a Platform token (User menu, Your tokens, Add token), then:

```bash
export TOWER_ACCESS_TOKEN=<your-platform-token>
export SEQERA_ACCESS_TOKEN=$TOWER_ACCESS_TOKEN
```

## The flow

```bash
cd side-quests/platform_automation/terraform/pipeline

# 1. Set your inputs. Either copy the example file:
cp terraform.tfvars.example terraform.tfvars
# and edit it, or pass -var flags on each command below.

# 2. Initialise (downloads the seqera 0.30.5 provider).
terraform init

# 3. Preview. One resource to add. Nothing changes yet.
terraform plan -var="username=$USER"

# 4. Add the pipeline to the Launchpad.
terraform apply -var="username=$USER"
```

Open the workspace Launchpad. `rnaseq-nf-$USERNAME` is there. No run has
started.

## Idempotency

Run apply again:

```bash
terraform apply -var="username=$USER"
```

"No changes. Your infrastructure matches the configuration." Terraform manages
existence, not actions. Applying ten times does not launch ten runs. It
launches zero. To actually run the pipeline, use seqerakit/tw (see
`../../seqerakit`). That is the imperative half: do the thing, now.

## Teardown

Remove your pipeline when you are done:

```bash
terraform apply -destroy -var="username=$USER"
```

It deletes only the pipeline you added. The shared compute environment and
workspace are untouched.
