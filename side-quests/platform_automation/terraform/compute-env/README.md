# terraform/compute-env

The Admin tier's compute-environment asset. See the Platform Automation side quest, §1 (Admin). You run this
only if you have the **Admin** (or Owner) role and permissions in the cloud
environment. In a workshop the presenter builds the shared compute environment
before the session; Maintain-role attendees do not run this and instead get a
`compute_env_id` to build pipelines against.

## What it does

The realistic, complex path. One config, two providers, end to end:

1. **azurerm** creates a head pool and a worker pool on the customer's existing
   Azure Batch account, on their existing subnet, running as their existing
   managed identity.
2. **seqera** stores the Azure credentials and creates a manual Azure Batch
   compute environment that uses the head pool and routes tasks to the worker
   pool.

This is what "manage your cloud as code" looks like for a lab that already
owns its Azure footprint. Terraform does not create the Batch account, the
network, or the identity. Those exist under change control and are referenced
with data sources (in `main.tf`). Terraform adds the pools and the Platform
wiring.

## Three tiers, for context

The session shows compute environments at three levels of control:

1. **Batch Forge (UI)**: the Platform creates and scales the pool. Easiest.
2. **Manual CE pointing at an existing pool**: you already have a pool, the
   Platform just submits to it.
3. **This config**: Terraform provisions the pools on your existing account and
   network, then wires the Platform. Full infrastructure as code.

## Manual vs Forge

The compute environment is manual: `head_pool` is set and there is no `forge`
block. The Platform submits to the pools this config created. It never creates
or scales a pool itself. A `forge` block would flip it to Batch Forge.

## Run it

```bash
az login
export TOWER_ACCESS_TOKEN=<your-platform-token>
export SEQERA_ACCESS_TOKEN=$TOWER_ACCESS_TOKEN

cp terraform.tfvars.example terraform.tfvars
# edit terraform.tfvars: subscription, workspace, existing resource names, keys

terraform init
terraform plan
terraform apply
```

After apply, note the `compute_env_id` output. The Maintain tier needs it for
`terraform/pipeline`.

## The cloud side

With Admin and cloud access you can also see the Azure side: the Batch account,
the pool VMs, jobs and tasks stacking onto them, a task's logs and exit code.
That is the whole point of the tier: the Platform submits to Azure Batch, which
runs the work. Maintain- and Launch-role users see only the compute environment
in the workspace, not the cloud behind it.
