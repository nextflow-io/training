# Seqera Platform Automation

The Seqera Platform does not run your work.
It is an API and a control plane over your cloud.
It hands jobs to a compute environment, the compute environment runs them on cloud VMs, and the Platform reads back state, logs, and exit codes.

Everything the web UI does, it does by calling the Platform API. So everything is automatable: with the API, Terraform, and the CLI you can manage compute environments, pipelines, and runs as code, with no clicking. This side quest walks that programmatic surface from the most privileged role to the least.

### Learning goals

In this side quest, we'll drive the Platform across three workspace roles of decreasing permission: **Admin**, **Maintain**, and **Launch**. Each role does one job and hands an artifact to the next. You'll learn how to:

- Create a compute environment two ways, with increasing control over the cloud: the UI (Batch Forge) and Terraform
- Add a pipeline to the Launchpad declaratively with Terraform, and see idempotency
- Launch a pipeline using the GUI, CLI and `seqerakit`.
- Tell declarative existence (Terraform) apart from imperative actions (the UI, `seqerakit`, `tw`) and learn when to use each

### Prerequisites

Before taking on this side quest, you should:

- Have a Seqera Platform account on Seqera Cloud [https://cloud.seqera.io](`https://cloud.seqera.io`), or an Enterprise install
- Be a member of a workspace with a role: Admin, Maintain, or Launch
- Have credentials for your cloud provider already added to that workspace (not covered here, see [credentials overview](https://docs.seqera.io/platform-cloud/credentials/overview))
- Ideally, GitHub credentials are added to your workspace. This is optional, but it would help prevent incidents caused by GitHub rate limits.
- Be comfortable with the command line and basic Nextflow concepts

### Know your role

The work is split by role, from most to least privileged. Start at the highest tier your access allows:

| Tier                                                    | Role          | Can change                               | Produces              |
| ------------------------------------------------------- | ------------- | ---------------------------------------- | --------------------- |
| [Admin](#1-admin-compute-environments-and-cloud)        | Owner / Admin | Compute environments and cloud resources | a Compute environment |
| [Maintain](#2-maintain-add-a-pipeline-to-the-launchpad) | Maintain      | Pipelines on the Launchpad               | a Launchpad pipeline  |
| [Launch](#3-launch-run-a-pipeline)                      | Launch        | Nothing; can only run                    | pipeline runs         |

An Admin user will learn how to create compute environments. If you only have Maintain, your workspace must already include one valid compute environment; start at section 2. If you only have Launch, your workspace must have a pre-configured pipeline (and optionally an Action) to run; start at section 3. Each section opens with what it requires.

---

## 0. Get started

### Open the training codespace

The Codespace contains all the tools you need (Terraform, `tw`, `seqerakit`); you
install nothing yourself. Open it now and read on while it builds; it is ready
when the terminal returns to a prompt.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

### Move into the project directory

The Codespace terminal opens at the repository root (`/workspaces/training`). All the assets for this side quest live under `side-quests/platform_automation/`, so move there now:

```bash
cd /workspaces/training/side-quests/platform_automation
```

Focus VSCode on this directory so the file explorer shows the assets you'll edit:

```bash
code .
```

### Review the materials

The directory holds the Terraform and `seqerakit` configurations for each role. Each section below tells you which subdirectory to `cd` into:

??? abstract "Directory contents"

    ```console
    .
    ├── terraform
    │   ├── compute-env        # section 1 (Admin): provision the cloud + compute environment
    │   │   ├── main.tf
    │   │   ├── variables.tf
    │   │   └── terraform.tfvars.example
    │   └── pipeline           # section 2 (Maintain): add a pipeline to the Launchpad
    │       ├── main.tf
    │       ├── variables.tf
    │       └── terraform.tfvars.example
    └── seqerakit              # sections 2 & 3: add and launch a pipeline
        ├── add-rnaseq.yml
        ├── launch-rnaseq.yml
        └── launch-rnaseq-multiple.yml
    ```

### Check workspace access

This training requires you to be in a workspace. Navigate to [Seqera Platform](https://cloud.seqera.io/) or your Enterprise install and check you can see the relevant workspace which should be under an organization. There should already be credentials pre-configured and if you have Maintain or lower permissions, a compute environment and/or pipeline should already exist. If you cannot see the workspace, please ask your workspace owner to add you to one.

### Create an access token

We authenticate to the Platform to tell it who we are. That is what grants the permissions of our role. Seqera uses a **personal access token**: a bearer token you create once and send with every API call, in an `Authorization: Bearer <token>` header. The token carries your identity and role, so the Platform knows who you are and what you can do. Anyone with the token can act as you, so keep it secret. In this training we save it as an environment variable.

1. In the Platform, open the user menu (top right) and choose **Your tokens**.
2. Click **Add token**, name it (e.g. `platform-automation`), and click **Add**.
3. Copy the token now. The Platform shows it only once.
4. In the Codespace terminal, export it under both names:

```bash
export TOWER_ACCESS_TOKEN=<paste-token>
```

Terraform, `tw`, and `seqerakit` all read the token from these variables. Check it worked:

```bash
tw info
```

This prints the API endpoint and the authenticated user. If it errors, the token is wrong or not exported.

!!! warning

    This only exists in a single terminal session, if you open a new terminal you will need to export them again.

For Seqera Enterprise, also set the environment variable `TOWER_API_ENDPOINT` to your install's API URL. Everything else is identical.

```bash
export TOWER_API_ENDPOINT=<your-api-url>
```

In Terraform, you must set the `server_url` variable to the same value.

```hcl
provider "seqera" {
  server_url = "<your-api-url>"
}
```

### Find your organization and workspace IDs

Later steps need the **numeric workspace ID**. The web UI navigates by name, so the number is not in the URL; find it like this:

1. In the Platform, click your organization name to open the organization page.
2. Open the **Workspaces** tab. Each workspace is listed with its numeric ID. That
   number is the `workspace_id` the Terraform and API steps ask for.
3. You only need the numeric organization ID if you list workspaces over the API; if you read the workspace ID from the Workspaces tab, you can skip it.

Or list everything your token can reach from the terminal:

```bash
tw workspaces list
```

This prints a table with the workspace ID, the workspace name, and the organization name.

Two forms of the same workspace turn up in this module. Terraform and the API want the **numeric** ID (`workspace_id`). `tw` and `seqerakit` accept either the numeric ID or the `Organization/Workspace` **name** (e.g. `my-org/platform-automation`). Export the numeric ID once so `tw` targets the shared workspace without a `--workspace` flag on every command:

```bash
export TOWER_WORKSPACE_ID=<numeric-workspace-id>
```

#### Pick a workshop handle

In a shared workspace, every pipeline you add needs a name no one else is using. Pick a short handle and export it once; everything from section 2 onward reuses this one variable. Use only letters, numbers, and hyphens. Launchpad names cannot contain dots or spaces:

```bash
export WORKSHOP_USER=<your-handle>   # e.g. ada, jdoe, team-rocket
```

We set this by hand rather than reading `$USER`: in the Codespace everyone shares the same Linux user, so `$USER` is identical for all learners and would collide.

---

## 1. Admin: compute environments and cloud

**Requires:** Owner or Admin role, and permissions in the cloud environment. If you are using a Cloud provider, authenticate with them first via their CLI (e.g. `aws sso login`). There must be valid credentials in the workspace to perform the operations in the Cloud and build compute environments, please check the Seqera documentation for further details.

The same compute environment can be created with increasing control over the cloud. We build it two ways.

### 1.1. Click-ops a compute environment with Batch Forge

Batch Forge is the hands-off option: the Platform reaches into your cloud and creates the resources for you. Most convenient, least control.

In the side bar, open **Compute Environments** and click **Create compute environment**. Give it a name, pick your cloud platform, and select the credentials for that cloud. The rest of the form differs by provider; recommended settings:

**AWS**

- Region
- Work directory (S3 bucket)
- Wave, Fusion, Fast Instance should all be enabled
- Config mode should be Batch Forge
- Provisioning model should be Spot, using Spot instances for nodes which run Nextflow tasks.

**Azure**

- Location
- Work directory (Azure Storage container)
- Wave and Fusion should be enabled.
- Select "Separate head and worker pools" to create one pool for the Nextflow head job and a separate pool for the tasks.
- For the head pool VM type, select a VM you have sufficient quota for. Standard_D2s_v3 (2 CPUs) is a good starting choice.
- Set the head pool VM count to 1. The head pool runs a single Nextflow head job; it does not need more.
- For the worker pool, select a VM type you have sufficient quota for. Standard_D4s_v3 (4 CPUs) is a good starting choice.
- Set the worker pool VM count to 4
- Leave autoscaling enabled for both pools

Leave the other settings at their defaults. Keep **Dispose resources** enabled so the pool is torn down when you delete the compute environment. Click **Add** in the top right.

Watch your cloud console and you will see Forge create the resources: an identity and roles for Nextflow (access to blob storage and the Batch service), a more limited role for the worker tasks (storage only), the pools, and their networking. Many moving parts, all handled for you.

That convenience is the trade-off: Forge owns those resources and gives you little control. A lab that already runs on managed infrastructure wants the opposite, to describe the resources itself and wire the Platform to them. That is the Terraform path.

### 1.2. Read it back over the API with `tw`

The form you just filled in did nothing special: it sent a configuration to the Platform API.
Ask the API for what you made, from the terminal, with nothing retyped.

`tw` reads the same workspace the UI does. List the compute environments and the one you created in 1.1 is there:

```bash
tw compute-envs list
```

Now export its full configuration to a file.
Use the name you gave it in 1.1:

```bash
tw compute-envs export forge-ce.json --name="<the name you gave it in 1.1>"
```

`tw` reads the workspace from `TOWER_WORKSPACE_ID` (set in section 0), so no `--workspace` flag is needed.
The file `forge-ce.json` is the compute environment as the API stores it:

```json title="forge-ce.json"
{
  "computeEnv": {
    "name": "azure-batch",
    "platform": "azure-batch",
    "config": {
      "workDir": "az://work",
      "region": "eastus",
      "waveEnabled": true,
      "fusion2Enabled": true,
      "forge": {
        "headPool": {
          "vmType": "Standard_D2s_v3",
          "vmCount": 1,
          "autoScale": true
        },
        "workerPool": {
          "vmType": "Standard_D4s_v3",
          "vmCount": 4,
          "autoScale": true
        },
        "disposeOnDeletion": true
      }
    }
  }
}
```

The `forge` block is exactly what the checkboxes set in 1.1. Its presence is what tells the Platform to create and scale the pools for you. In the next section, the absence of this block is what makes the Terraform compute environment manual.

Close the loop: import the file to recreate the compute environment, with nothing retyped.

```bash
tw compute-envs import forge-ce.json --name="azure-batch-copy"
```

The UI sent this JSON to create the compute environment; `tw export` reads it back, and `tw import` sends it again. The GUI and the CLI are two clients of the same API.

The exported JSON does not include a credentials ID. If your workspace has more than one credential for this cloud, `tw` cannot tell which to use and fails with `Multiple credentials match this compute environment`. Name one explicitly with `-c` (see `tw credentials list`):

```bash
tw compute-envs import forge-ce.json --name="azure-batch-copy" --credentials="<your-azure-credentials-name>"
```

### 1.3. Provision the cloud with Terraform

Terraform manages cloud resources declaratively: you describe what should exist and Terraform makes it so. It does not run your work; it only creates the compute environment. `side-quests/platform_automation/terraform/compute-env/` does it in a single `main.tf`, two providers in one apply:

- data sources: the existing Batch account and the Azure credentials already
  stored in the Platform, referenced, not created.
- `azurerm`: creates a head pool and a worker pool on that account.
- `seqera`: creates the compute environment pointing at those pools.

The manual marker: `head_pool` is set to the pool Terraform just made and there is **no `forge` block**. That one difference is what makes the compute environment manual instead of Forge. `nextflow_config` routes tasks to the worker pool.

Let's walk the key blocks of `terraform/compute-env/main.tf`.

The first block declares the providers and pins their versions. The `seqera` provider is pinned to exactly `0.40.1`:

```terraform
terraform {
  required_version = ">= 1.9"

  required_providers {
    seqera  = { source = "seqeralabs/seqera", version = "0.40.1" }
    azurerm = { source = "hashicorp/azurerm", version = ">= 3.0" }
    random  = { source = "hashicorp/random", version = ">= 3.0" }
  }
}
```

Next we configure the providers. The `seqera` provider reads `TOWER_ACCESS_TOKEN` from the environment, so there is no token argument here:

```terraform
provider "azurerm" {
  features {}
  subscription_id = var.subscription_id
}

provider "seqera" {
  server_url = var.server_url
}
```

Next we can define some variables we can use throughout our deployment. In this case, we will set the node sizes and details once and reuse them throughout the configuration:

```terraform
locals {
  head_vm_size      = "Standard_E4ds_v5"
  worker_vm_size    = "Standard_E16ds_v5"
  worker_max_nodes  = 8
  node_agent_sku_id = "batch.node.ubuntu 24.04"

  # The Platform requires a pool's task slots to equal the node's vCPU count.
  # Azure VM size names embed that count (Standard_E16ds_v5 -> 16), so derive it
  # from the size rather than hardcoding a number that can drift out of sync.
  head_vm_cores   = tonumber(regex("^Standard_[A-Z]+([0-9]+)", local.head_vm_size)[0])
  worker_vm_cores = tonumber(regex("^Standard_[A-Z]+([0-9]+)", local.worker_vm_size)[0])
}
```

Next we can refer to resources that already exist. These are called `data` in Terraform. Here we refer to the Azure Batch account:

```terraform
# Existing Azure resources the customer owns. Referenced, not managed.
data "azurerm_batch_account" "existing" {
  name                = var.batch_account_name
  resource_group_name = var.batch_account_rg
}
```

Now the resources Terraform does create. The pool names must stay unique across re-creates, so we first generate a short random suffix to append to them:

```terraform
# Keeps pool names unique across re-creates.
resource "random_string" "suffix" {
  length  = 6
  special = false
  upper   = false
}
```

Then the node pools themselves. Terraform creates two pools in Azure Batch with the `azurerm_batch_pool` resource. You can add as many as you like, but the details must match the values the Azure provider expects. Here we start a fixed-size head pool and a dynamically sized worker pool:

```terraform
# Head pool: runs the Nextflow head job. One node is enough.
resource "azurerm_batch_pool" "head" {
  name                = "rnaseq-head-${random_string.suffix.result}"
  display_name        = "rnaseq head pool"
  resource_group_name = var.batch_account_rg
  account_name        = var.batch_account_name
  vm_size             = local.head_vm_size
  node_agent_sku_id   = local.node_agent_sku_id
  max_tasks_per_node  = local.head_vm_cores

  fixed_scale {
    target_dedicated_nodes = 1
  }

  storage_image_reference {
    publisher = "microsoft-dsvm"
    offer     = "ubuntu-hpc"
    sku       = "2404"
    version   = "latest"
  }

  container_configuration {
    type = "DockerCompatible"
  }
}

# Worker pool: runs the pipeline tasks. Autoscales 0..worker_max_nodes.
resource "azurerm_batch_pool" "worker" {
  name                = "rnaseq-worker-${random_string.suffix.result}"
  display_name        = "rnaseq worker pool"
  resource_group_name = var.batch_account_rg
  account_name        = var.batch_account_name
  vm_size             = local.worker_vm_size
  node_agent_sku_id   = local.node_agent_sku_id
  max_tasks_per_node  = local.worker_vm_cores

  auto_scale {
    evaluation_interval = "PT5M"
    formula             = <<-FORMULA
      pending = avg($PendingTasks.GetSample(180 * TimeInterval_Second));
      $TargetDedicatedNodes = min(${local.worker_max_nodes}, pending);
      $NodeDeallocationOption = taskcompletion;
    FORMULA
  }

  storage_image_reference {
    publisher = "microsoft-dsvm"
    offer     = "ubuntu-hpc"
    sku       = "2404"
    version   = "latest"
  }

  container_configuration {
    type = "DockerCompatible"
  }

  # Terraform manages the pool, not its live scale.
  lifecycle {
    ignore_changes = [auto_scale, fixed_scale]
  }
}
```

Each pool in `main.tf` also carries a `start_task` (elided above) that runs once per node at startup. Forge adds an equivalent task to the pools it creates; manual pools get nothing, so we replicate the essential parts: install `azcopy` (Nextflow's Azure Batch executor uses it to stage data) and register the AppArmor profile Fusion needs on Ubuntu 24.04 nodes. See the `node_start_task_script` local in the file for the full script.

Once we've built the Azure Batch node pools, we add the Seqera Platform compute environment that points to them. The compute environment needs Azure credentials to talk to your cloud. We added those to the workspace earlier (see Prerequisites), so we look them up by name rather than storing keys here. The `seqera_credentials` data source lists the credentials in the workspace, and we pick out the one whose name matches:

```terraform
# Azure credentials already stored in the Platform, looked up by name.
data "seqera_credentials" "all" {
  workspace_id = var.workspace_id
}

locals {
  azure_credentials_id = one([
    for c in data.seqera_credentials.all.credentials : c.id
    if c.name == var.azure_credential_name
  ])
}
```

Finally, the compute environment itself. This is the manual marker in code: `head_pool` points at the pool we just created and `nextflow_config` routes tasks to the worker pool's queue. `enable_wave` and `enable_fusion` turn on Wave (container provisioning) and the Fusion file system, the same checkboxes you would tick in the UI; with Fusion enabled, Wave must be too. Because it references `azurerm_batch_pool.head.name` and `azurerm_batch_pool.worker.name`, Terraform knows to create the pools first:

```terraform
resource "seqera_compute_env" "main" {
  workspace_id = var.workspace_id

  compute_env = {
    name           = "azure-batch-manual"
    description    = "Azure Batch CE on Terraform-managed pools"
    platform       = "azure-batch"
    credentials_id = local.azure_credentials_id

    config = {
      azure_batch = {
        region          = var.region
        work_dir        = var.work_dir
        head_pool       = azurerm_batch_pool.head.name
        enable_wave     = true
        enable_fusion   = true
        nextflow_config = "process.queue = '${azurerm_batch_pool.worker.name}'\n"
      }
    }
  }
}
```

The last block is the output the next tier needs:

```terraform
output "compute_env_id" {
  value = seqera_compute_env.main.compute_env_id
}
```

Terraform reads those references and works out the order itself: pools first, then the compute environment that depends on them.

Now, we can run Terraform to apply these changes.

```bash
# you may need to log in to the cloud provider with `az login`
terraform init
terraform apply
```

After you run `terraform apply`, it will ask you to fill in the details defined in `variables.tf`

!!! tip "Tearing it down"

    You can save them as a .tfvars file and Terraform will read them automatically, so you don't have to type them each time. See `terraform.tfvars.example` for the full variable reference.
    Alternatively, you can save them as an environment variable preceded by `TF_VAR_`, for example `TF_VAR_subscription_id=00000000-0000-0000-0000-000000000000`.

Terraform will show you a preview, if it looks accurate, type `yes` to apply the changes. It will create the pools and the compute environment, and print the `compute_env_id` you hand to the Maintain tier.

### 1.4. See both sides

The compute environment now exists in two places, and as an Admin with cloud access you can see both.

On the Platform side, read the ID Terraform produced:

```bash
terraform output compute_env_id
```

On the Azure side, open the Batch account in the portal. You will see the head and worker pools Terraform created, sitting idle with no jobs yet. Once you launch a pipeline (sections 2 and 3), jobs and tasks stack onto these pools, and you can drill into a task's logs and exit code. That is the whole point of the Admin tier: the Platform submits to Azure Batch, and Azure Batch runs the work. Maintain- and Launch-role users see only the compute environment in the workspace, not the cloud behind it.

### Takeaway

One compute environment, three levels of control: Forge in the UI (easiest, the Platform owns the pool), `tw` (export and import the same config over the API), and Terraform (you own the pools and the wiring, end to end). The artifact you hand to the Maintain tier is a `compute_env_id`.

!!! note "Tearing it down"

    Terraform manages the cloud resources it created, so it can remove them too. From `terraform/compute-env`, run `terraform apply -destroy` to delete the compute environment and the Batch pools in one step. The existing Batch account and credentials are referenced, not managed, so they are left untouched.

---

## 2. Maintain: add a pipeline to the Launchpad

**Requires:** Maintain role, and a `compute_env_id` (from section 1 or an Admin on your team) plus the numeric `workspace_id`. Maintain manages pipelines, not compute environments. That split is deliberate: Admin owns the compute environment, you own your pipelines.

You add the same pipeline, `rnaseq-nf-$WORKSHOP_USER`, four ways. The first three are imperative, you run a command and it acts. The last, Terraform, is declarative. Watch what happens when you run each one twice.

### 2.1. Add a pipeline via the UI

In the workspace, open the **Launchpad** and click **Add pipeline**. Fill in the form:

- **Name**: `rnaseq-nf-` followed by your handle, e.g. `rnaseq-nf-ada` (the handle you exported as `WORKSHOP_USER` in section 0).
- **Compute environment**: the one from section 1, e.g. `azure-batch-manual`.
- **Pipeline to launch**: `https://github.com/nextflow-io/rnaseq-nf`.
- **Revision**: `master`.
- **Work directory**: your Azure Blob work dir, e.g. `az://work`.

Click **Add**. The pipeline appears on the Launchpad with no run started. Adding a pipeline only saves a launch configuration; it does not run anything.

To remove from the workspace, you can click the hamburger menu on the top right and click **Delete**.

### 2.2. Add a pipeline via `tw`

The CLI does the same thing in one command:

```bash
tw pipelines add \
  --name="rnaseq-nf-$WORKSHOP_USER" \
  --compute-env="azure-batch-manual" \
  --work-dir="az://work" \
  --revision="master" \
  https://github.com/nextflow-io/rnaseq-nf
```

Run it again and `tw` errors: a pipeline with that name already exists. The command is imperative, so each invocation tries to add a pipeline; it has no notion of "already in the desired state".

To remove the pipeline, you can use the `tw` command line again:

```bash
tw pipelines delete --name="rnaseq-nf-$WORKSHOP_USER"
```

### 2.3. Add a pipeline with `seqerakit`

`seqerakit` is a wrapper over `tw` that reads a YAML file and runs the underlying `tw` commands. It keeps the configuration as code, so a teammate can reproduce the exact same pipeline. `seqerakit/add-rnaseq.yml` describes the pipeline:

```yaml
pipelines:
  - name: "rnaseq-nf-${WORKSHOP_USER}"
    url: "https://github.com/nextflow-io/rnaseq-nf"
    workspace: "${TOWER_WORKSPACE_ID}"
    description: "Added with seqerakit"
    compute-env: "azure-batch-manual"
    work-dir: "az://work"
    revision: "master"
```

`seqerakit` expands variables like `${TOWER_WORKSPACE_ID}` and `${WORKSHOP_USER}` from the environment. Both were set in section 0, so just add the pipeline:

```bash
cd /workspaces/training/side-quests/platform_automation/seqerakit
seqerakit add-rnaseq.yml
```

If you run it again, it errors, the same way `tw` did, because the pipeline already exists:

```console
The pipelines resource already exists and will not be created. Please set 'on_exists: overwrite' to replace the resource or set 'on_exists: ignore' to ignore this error.
```

You can force it through with `seqerakit add-rnaseq.yml --overwrite`, which deletes and recreates the pipeline. But that is you telling it to repeat the action. By default, adding twice is an error. To make "add once, and leave it alone after that" the _default_ behaviour, we need a tool that manages existence rather than actions: Terraform.

### 2.4. Add a pipeline with Terraform

`side-quests/platform_automation/terraform/pipeline` adds the same pipeline declaratively. You describe the pipeline that should exist; Terraform makes the workspace match. Like the imperative methods above, the config tracks the `master` branch; you will pin it to a release tag below to see Terraform update the pipeline in place.

```bash
cd /workspaces/training/side-quests/platform_automation/terraform/pipeline
terraform init
terraform apply -var="username=$WORKSHOP_USER" -var="workspace_id=<ID>" -var="compute_env_id=<CE_ID>"
```

!!! tip

    Copy `terraform.tfvars.example` to `terraform.tfvars` and fill it in so you stop passing `-var` flags.

Open the Launchpad: `rnaseq-nf-$WORKSHOP_USER` is there, with no run. Now apply again:

```bash
terraform apply -var="username=$WORKSHOP_USER" -var="workspace_id=<ID>" -var="compute_env_id=<CE_ID>"
```

```console
No changes. Your infrastructure matches the configuration.
```

That is the difference. `tw` and `seqerakit` errored on the second run because they perform an action every time. Terraform manages whether the pipeline **exists**: it is already there, so there is nothing to do. Apply ten times, still one pipeline. This is what keeps you from accumulating competing, half-duplicated resources in a shared workspace. To remove the pipeline, `terraform apply -destroy` with the same vars.

**Update in place.** Terraform does more than create-or-skip; it reconciles the pipeline to whatever the config says. The config tracks the `master` branch, which moves as the pipeline gets new commits. Pin it to a fixed release tag instead, so every launch runs the same code. Edit one line: in `main.tf`, line 35 of the `launch` block, change the revision from the `master` branch to the `v2.4` release tag.

=== "After"

    ```hcl title="main.tf" hl_lines="4" linenums="32"
      launch = {
        pipeline = "https://github.com/nextflow-io/rnaseq-nf"
        # The master branch moves; pin a release tag (e.g. v2.4) for reproducible launches.
        revision       = "v2.4"
        compute_env_id = var.compute_env_id
        work_dir       = var.work_dir
      }
    ```

=== "Before"

    ```hcl title="main.tf" hl_lines="4" linenums="32"
      launch = {
        pipeline = "https://github.com/nextflow-io/rnaseq-nf"
        # The master branch moves; pin a release tag (e.g. v2.4) for reproducible launches.
        revision       = "master"
        compute_env_id = var.compute_env_id
        work_dir       = var.work_dir
      }
    ```

Then apply:

```console
Terraform will perform the following actions:

  # seqera_pipeline.rnaseq_nf will be updated in-place
  ~ resource "seqera_pipeline" "rnaseq_nf" {
      ~ launch = {
          ~ revision = "master" -> "v2.4"
            # (4 unchanged attributes hidden)
        }
    }

Plan: 0 to add, 1 to change, 0 to destroy.
```

One pipeline, modified; not a second pipeline, and not an error. That is the third declarative behaviour, alongside create and no-op: update to match. Moving from a branch to a release tag is exactly the kind of change you want under version control: the pinned revision is now recorded in `main.tf`.

!!! note "Pinning a revision vs. Launchpad versions"

    The Platform also has a separate concept of **Launchpad pipeline versions**: several saved launch configurations on one pipeline, with one marked as the default (most recent). You can use this to maintain and control pipeline versions within one concept of a "pipeline".

    The Terraform provider can also publish and promote those versions, and flag drift if someone changes the default in the UI. See the provider's [pipeline versioning guide](https://github.com/seqeralabs/terraform-provider-seqera/blob/master/docs/guides/pipeline-versioning.md).

### Takeaway

Four ways to add one pipeline, two mental models. `tw` and `seqerakit` are imperative: each run is an action, and adding twice is an error. Terraform is declarative: it manages existence, so a second apply is a no-op or an update-in-place. The artifact you hand to the Launch tier is the Launchpad pipeline `rnaseq-nf-$WORKSHOP_USER`.

---

## 3. Launch: run a pipeline

**Requires:** Launch role and a token with that role, plus a pipeline on the Launchpad (added by a Maintainer in section 2). Launch can only run things; it cannot create or modify compute environments, pipelines, or Actions.

Everything in this section runs the pipeline the Maintainer already configured. None of it changes the pipeline; launching is an action, not a state.

### 3.1. Use the GUI

Open the **Launchpad**, select `rnaseq-nf-$WORKSHOP_USER`, and click **Launch**. The compute environment, revision, and work directory are already filled in by the Maintainer, so a Launch user just clicks the button. Submit it and the run appears under **Runs**.

### 3.2. Use the `tw` CLI

```bash
tw launch rnaseq-nf-$WORKSHOP_USER
```

It should show something like:

```console
  Workflow 2ZXaU1AzEn7Onk submitted at [Organization / Workspace] workspace.

    https://cloud.seqera.io/orgs/Organization/workspaces/Workspace/watch/2ZXaU1AzEn7Onk
```

You can monitor the run by clicking the provided URL.

!!! note "Using a params file"

    If you wish to provide parameters to the pipeline, you can do so with the `--params-file` flag.

### 3.3. Use `seqerakit`

`seqerakit` launches a pipeline that already exists on the Launchpad, filling in the `tw launch` command from `seqerakit/launch-rnaseq.yml`. `WORKSHOP_USER` is set from section 0; set the compute environment name, dry run first, then launch:

```bash
cd /workspaces/training/side-quests/platform_automation/seqerakit
export COMPUTE_ENVIRONMENT=<compute environment name>

seqerakit launch-rnaseq.yml --dryrun   # prints the tw command, changes nothing
seqerakit launch-rnaseq.yml            # launches for real
```

The dry run shows the underlying `tw` command. Run it twice and you get two runs. That is the imperative model: do the thing, now.

One advantage of using Seqerakit is the YAML forms a template of actions to perform. Because of this we can do two more things: firstly, we can save it and re-use it later. Secondly, we can launch the same pipeline several times by adding more launch blocks to the YAML file. This is useful if you want to launch the same pipeline with different parameters or on different compute environments. Let's try and launch the pipeline now with

```bash
seqerakit launch-rnaseq-multiple.yml
```

### 3.4. Side note: Examine the cloud resources

It's worth taking a second to examine the cloud resources here.

A run does not create one neatly named cloud job. The platform submits a single job to the compute environment who's only job is to configure and run Nextflow. The Nextflow running within this job submits **many** Batch jobs and tasks, one per process invocation. Overall, this means you will not find a single resource named after the Platform run. Open the run on the Platform (**Runs** → your run), which mirrors the underlying Batch service: the task table here is the same work you can see in the cloud console.

If you have cloud access, the prefixes Nextflow uses in each Batch service are:

- **Azure Batch**: jobs named `job-<hex>-<process>`, tasks named `nf-<taskhash>`.
- **AWS Batch**: job names of the form `<run-name>_<process>`, with unsupported characters stripped (max 128 chars).
- **GCP Batch**: job IDs of the form `nf-<taskhash>-<timestamp>`.

The relationship is the point: the Platform hands work to Batch, Batch runs it on the pool VMs, and the Platform reads back state, logs, and exit codes.

### 3.5. Side note: launch Actions

There is one more way to launch, built for automation rather than people. An **Action** is a saved launch configuration behind a single URL: call that URL and the pipeline runs, with no Launchpad and no `tw`. A Maintainer creates one in the UI under **Actions** → **Add action**, picks the pipeline, compute environment and other settings, and saves it as a configuration. The Platform then shows a ready-made `curl` command for the Action's endpoint.

The split mirrors the roles: _creating_ an Action needs the Maintain role, but _triggering_ it needs only a Launch token. That is the point of the stratified launch: a Maintainer hands launchers (or an automation system) a URL they can call to run the pipeline, and nothing more.

### Takeaway

Several ways to launch, all imperative: the GUI button, `tw launch`, `seqerakit`, and a launch Action's URL. Each does the thing, now; run it twice and you get two runs. That is the opposite of Terraform, which manages whether something exists.

---

## Summary

You drove the Seqera Platform programmatically across three roles, each handing an artifact to the next: Admin builds the compute environment and owns the cloud, Maintain adds pipelines and Launch triggers a run.

### Key patterns

- **Everything is one API.** The GUI, Terraform, `tw`, and `seqerakit` all call the same Platform API. Anything you can click, you can automate.
- **Use the right tool.** Terraform manages resources as state: what should exist, where. `tw` and `seqerakit` act on the Platform imperatively: they do the thing, now.

|              | Terraform     | `seqerakit` / `tw` / Action |
| ------------ | ------------- | --------------------------- |
| Manages      | existence     | actions                     |
| Run twice    | no-op         | two runs                    |
| Mental model | desired state | do the thing, now           |

- **Roles stratify what each token can do.** A Maintain token defines pipelines and Actions; a Launch token can only trigger an Action and nothing else. The Action is the safe handoff from maintainers to launchers (and to automation).

## What's next?

- The AI half of the workshop is the [Building with Seqera AI (CoScientist)](../co_scientist/index.md) side quest. It uses the same `rnaseq-nf` pipeline and API endpoints, but drives them via AI agents.
