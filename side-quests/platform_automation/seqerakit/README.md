# seqerakit

The imperative asset, used in two tiers. See the Platform Automation side quest, §2 (Maintain) and §3 (Launch). seqerakit wraps `tw` and reads a YAML file:

- `add-rnaseq.yml` (Maintain) adds `rnaseq-nf-$USERNAME` to the Launchpad. Re-running errors unless `--overwrite`, the imperative contrast to Terraform's idempotent `apply`.
- `launch-rnaseq.yml` (Launch) launches that pipeline. Run it twice and you get two runs.

Terraform owns the pipeline's existence; seqerakit and `tw` act on it. This is the imperative contrast to the API Action and to Terraform.

## Declarative vs imperative

|              | Terraform                            | seqerakit / tw       |
| ------------ | ------------------------------------ | -------------------- |
| Manages      | existence (does the pipeline exist?) | actions (run it now) |
| Run twice    | no-op, still one pipeline            | two runs             |
| Mental model | desired state                        | do the thing, now    |

You added `rnaseq-nf-$USERNAME` with Terraform. Now you launch it.

## Setup

```bash
export TOWER_ACCESS_TOKEN=<your-platform-token>
export SEQERA_ACCESS_TOKEN=$TOWER_ACCESS_TOKEN
export USERNAME=$USER
```

Edit `workspace` and `compute-env` in the YAML to match the workshop workspace if they differ.

## Dry run first

Always dry run first. `--dryrun` prints the `tw` command seqerakit would run and changes nothing:

```bash
seqerakit launch-rnaseq.yml --dryrun
```

Then launch for real:

```bash
seqerakit launch-rnaseq.yml
```

## The tw equivalent

seqerakit is a wrapper over `tw`. The dry run shows you the underlying command.
It is the same as:

```bash
tw launch rnaseq-nf-$USERNAME \
  --workspace=my-org/platform-automation \
  --compute-env=azure-batch-manual
```

This is the callback to the Azure mechanics: launching here submits jobs to the Platform, which hands them to Azure Batch, which stacks tasks onto the pool VMs. Same handoff the Admin tier showed in the portal.
