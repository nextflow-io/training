# Seqera Platform Automation

The Seqera Platform does not run your work.
It is an API and a control plane over your cloud.
It hands jobs to a compute environment, the compute environment runs them on cloud VMs, and the Platform reads back state, logs, and exit codes.

Everything the web UI does, it does by calling the Platform API. So everything is automatable: with the API, Terraform, and the CLI you can manage compute environments, pipelines, and runs as code, with no clicking. This side quest walks that programmatic surface from the most privileged role to the least.

<figure class="excalidraw">
--8<-- "docs/en/docs/side_quests/platform_automation/img/one-api.excalidraw.svg"
</figure>

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

For this tutorial you will need the following tools:

- [`tw`](https://github.com/seqeralabs/tower-cli/)
- [`seqerakit`](https://github.com/seqeralabs/seqerakit)
- [Terraform](https://www.terraform.io/)

The Codespace should already contain all of these tools, but you can verify them by trying their respective CLI commands.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

### Move into the project directory

The Codespace terminal opens at the repository root (`/workspaces/training`). All the assets for this side quest live under `side-quests/platform_automation/`, so move there now:

```bash
cd side-quests/platform_automation
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
    └── seqerakit              # section 3 (Launch): launch a pipeline
        ├── launch-rnaseq.yml
        └── launch-rnaseq-multiple.yml
    ```

We'll cover the content of these files in a moment.

### Check workspace access

This training requires you to be in a workspace. Navigate to [Seqera Platform](https://cloud.seqera.io/) or your Enterprise install and check you can see the relevant workspace which should be under an organization. There should already be credentials pre-configured and if you have Maintain or lower permissions, a compute environment and/or pipeline should already exist. If you cannot see the workspace, please ask your workspace owner to add you to one.

### Create an access token

We authenticate to the Platform to tell it who we are. That is what grants the permissions of our role. Seqera uses a **personal access token**: a bearer token you create once and send with every API call, in an `Authorization: Bearer <token>` header. The token carries your identity and role, so the Platform knows who you are and what you can do. Anyone with the token can act as you, so keep it secret. In this training we save it as an environment variable.

1. In the Platform, open the user menu (top right) and choose **Your tokens**.
2. Click **Add token**, name it (e.g. `platform-automation`), and click **Add**.
3. Copy the token now. The Platform shows it only once.
4. In the Codespace terminal, export it:

```bash
export TOWER_ACCESS_TOKEN=<paste-token>
```

Terraform, `tw`, and `seqerakit` all read the token from this variable. Check it worked:

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

The same compute environment can be created with increasing control over the cloud.
You'll build it two ways.

!!! note "This quest uses Azure Batch"

    Section 1.1 covers AWS, Azure, GCP, on-prem HPC, and Kubernetes, but everything after it (the `tw` export, the Terraform, and the `az://work` work directory used throughout) assumes an Azure Batch compute environment.
    Users on other providers can follow the same pattern, but the provider-specific resources and work-directory paths will differ.

### 1.1. Create a compute environment in the UI

On the cloud Batch providers (AWS, Azure, GCP), the UI offers **Batch Forge**: the hands-off option where the Platform reaches into your cloud and creates the resources for you. Most convenient, least control. On-prem HPC clusters and Kubernetes use a different, non-Forge flow: the cluster already exists and you point the Platform at it.

In the side bar, open **Compute Environments** and click **Create compute environment**. Give it a name, pick your platform, and select the credentials for it. The rest of the form differs by provider. Select yours below for recommended settings:

=== "AWS"

    - Region
    - Work directory (S3 bucket)
    - Config mode should be Batch Forge
    - Provisioning model should be Spot, using Spot instances for nodes which run Nextflow tasks.

=== "Azure"

    - Location
    - Work directory (Azure Storage container)
    - Select "Separate head and worker pools" to create one pool for the Nextflow head job and a separate pool for the tasks.
    - For the head pool VM type, select a VM you have sufficient quota for. Standard_D2s_v3 (2 CPUs) is a good starting choice.
    - Set the head pool VM count to 1. The head pool runs a single Nextflow head job; it does not need more.
    - For the worker pool, select a VM type you have sufficient quota for. Standard_D4s_v3 (4 CPUs) is a good starting choice.
    - Set the worker pool VM count to 4
    - Leave autoscaling enabled for both pools

=== "GCP"

    - Location (region)
    - Work directory (Cloud Storage bucket, e.g. `gs://my-bucket`)
    - Enable Spot instances under GCP Resources to reduce cost.
    - Config mode should be Batch Forge

=== "On-prem / HPC"

    This is **not** Batch Forge: no pools are created for you. The cluster already exists and the Platform submits to it over SSH or the Tower Agent.

    - Compute platform: your scheduler (Slurm, LSF, PBS Pro, Grid Engine, Moab)
    - Credentials (SSH key or Tower Agent)
    - Login hostname: the cluster's login node
    - Work directory: an absolute path on the shared filesystem
    - Launch directory (defaults to the work directory if omitted)
    - Head queue name: the queue the Nextflow head job submits to
    - Compute queue name: the queue Nextflow submits tasks to

=== "Kubernetes"

    Also **not** Batch Forge: a manual flow against a cluster you prepare first (apply a YAML manifest creating the service accounts and role bindings).

    - Credentials (service-account token or X509 client certificate)
    - Control-plane URL (from `kubectl cluster-info`)
    - SSL certificate (from `~/.kube/config`)
    - Namespace (default `tower-nf`)
    - Head service account (default `tower-launcher-sa`)
    - Storage claim name (default `tower-scratch`)
    - Storage mount path (default `/scratch`)
    - Work directory: the storage mount path or a subdirectory of it
    - Head job CPUs and memory

On the cloud Batch providers (AWS, Azure, GCP), leave the other settings at their defaults, keep **Dispose resources** enabled so the pool is torn down when you delete the compute environment, and click **Add** in the top right.

Watch your cloud console and you will see Forge create the resources: an identity and roles for Nextflow (access to blob storage and the Batch service), a more limited role for the worker tasks (storage only), the pools, and their networking. Many moving parts, all handled for you. The on-prem and Kubernetes paths have nothing to forge; they submit to infrastructure you already run.

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
  "discriminator": "azure-batch",
  "region": "eastus",
  "workDir": "az://work",
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
    "disposeOnDeletion": true,
    "dualPoolConfig": true
  }
}
```

The `forge` block is exactly what the checkboxes set in 1.1. Its presence is what tells the Platform to create and scale the pools for you. In the next section, the absence of this block is what makes the Terraform compute environment manual.

This file can be used to recreate the compute environment, with nothing retyped.

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

The whole file is in the repo to read at your own pace, and most of it is standard Terraform: the providers, the locals, and the two Azure Batch pools.
Only one block encodes the lesson that makes this compute environment _manual_, so we'll focus there and treat the rest as reference.

The first block declares the providers, pins their versions and sets them up. The `seqera` provider is pinned to exactly `0.40.1`:

```terraform title="terraform/compute-env/main.tf" hl_lines="5"
terraform {
  required_version = ">= 1.9"

  required_providers {
    seqera  = { source = "seqeralabs/seqera", version = "0.40.1" }
    azurerm = { source = "hashicorp/azurerm", version = ">= 3.0" }
    random  = { source = "hashicorp/random", version = ">= 3.0" }
  }
}

provider "azurerm" {
  features {}
  subscription_id = var.subscription_id
}

provider "seqera" {
  server_url = var.server_url
}
```

Let's move into the compute env directory and initialise the Terraform provider:

```bash
cd terraform/compute-env
terraform init
```

You should see Terraform install the relevant providers and create a lockfile. This makes sure anyone using this Terraform is using the same provider versions. Terraform can and should be shared, so it is very important to make sure everyone uses the same versions!

Next we can define some variables we can use throughout our deployment. In this case, we will set the node sizes and details once and reuse them throughout the configuration:

```terraform title="terraform/compute-env/main.tf" hl_lines="2 3"
locals {
  head_vm_size      = "Standard_E4ds_v5"
  worker_vm_size    = "Standard_E8ds_v5"
  worker_max_nodes  = 8
  node_agent_sku_id = "batch.node.ubuntu 24.04"

  # The Platform requires a pool's task slots to equal the node's vCPU count.
  # Azure VM size names embed that count (Standard_E8ds_v5 -> 8), so derive it
  # from the size rather than hardcoding a number that can drift out of sync.
  head_vm_cores   = tonumber(regex("^Standard_[A-Z]+([0-9]+)", local.head_vm_size)[0])
  worker_vm_cores = tonumber(regex("^Standard_[A-Z]+([0-9]+)", local.worker_vm_size)[0])
}
```

Now the resources Terraform does create. In this code, Terraform creates two pools in Azure Batch with the `azurerm_batch_pool` resource. You can add as many as you like, but the details must match the values the Azure provider expects. We won't go into the details of the resource here, you can check them yourself, but you can see how some details refer to variables and outputs of other steps:

```terraform title="terraform/compute-env/main.tf"
# Head pool: runs the Nextflow head job. One node is enough.
resource "azurerm_batch_pool" "head" {
  name                = "rnaseq-head-${random_string.suffix.result}"
  display_name        = "rnaseq head pool"
  resource_group_name = var.batch_account_rg
  account_name        = var.batch_account_name
  vm_size             = local.head_vm_size
  node_agent_sku_id   = local.node_agent_sku_id
  max_tasks_per_node  = local.head_vm_cores
  # Etc...
```

Once we've built the Azure Batch node pools, we add the Seqera Platform compute environment that points to them. This is the manual marker in code: `head_pool` points at the pool we just created and `nextflow_config` routes tasks to the worker pool's queue. Because they reference `azurerm_batch_pool.head.name` and `azurerm_batch_pool.worker.name` respectively, Terraform knows to create the pools on Azure Batch first:

```terraform title="terraform/compute-env/main.tf" hl_lines="14 17"
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
        worker_pool     = azurerm_batch_pool.worker.name
        enable_wave     = true
        enable_fusion   = true
      }
    }
  }
}
```

Now, we can run Terraform to apply these changes, where it will create our Azure Batch pools and then the Seqera Platform compute environment.

```bash
terraform apply
```

After you run `terraform apply`, it will ask you to fill in the details defined in `variables.tf`

- `subscription_id` is the subscription of your Azure account, which you can find in the Azure Portal under **Subscriptions**.
- `workspace_id` is the numeric workspace ID you exported in section 0.
- `region` is the Azure region you want to deploy to, e.g. `eastus`.
- `work_dir` is the Azure Blob work directory you want to use, e.g. `az://work`.
- `batch_account_name` is the name of the Azure Batch account you want to use, which you can find in the Azure Portal under **Batch accounts**.
- `batch_account_rg` is the resource group of the Azure Batch account you want to use, which you can find in the Azure Portal under **Batch accounts**.
- `azure_credentials_name` is the ID of the Azure credentials you want to use, which you can find in the Platform under **Credentials**.
- `server_url` is the URL of your Seqera Platform instance, which you can find in the Platform under **Settings**. This will default to the cloud instance of Seqera Platform by default.

!!! tip "Re-using variables"

    You can save them as a .tfvars file and Terraform will read them automatically, so you don't have to type them each time. See `terraform.tfvars.example` for the full variable reference.
    Alternatively, you can save them as an environment variable preceded by `TF_VAR_`, for example `TF_VAR_subscription_id=00000000-0000-0000-0000-000000000000`.

Terraform will show you a preview, if it looks accurate, type `yes` to apply the changes. It will create the pools and the compute environment, and print the `compute_env_id` you hand to the Maintain tier.

### 1.4. See both sides

The compute environment now exists in two places, and as an Admin with cloud access you can see both.

On the Platform side, read the ID Terraform produced:

```bash
terraform output compute_env_id
```

If you have access to the Azure side, open the Batch account in the portal. You will see the head and worker pools Terraform created, sitting idle with no jobs yet. Once you launch a pipeline (sections 2 and 3), jobs and tasks stack onto these pools, and you can drill into a task's logs and exit code. That is the whole point of the Admin tier: to create cloud resources for other users to work on. Maintain- and Launch-role users see only the compute environment in the workspace, not the cloud behind it.

### Takeaway

One compute environment, three levels of control: Forge in the UI (easiest, the Platform owns the pool), `tw` (export and import the same config over the API), and Terraform (you own the pools and the wiring, end to end). The artifact you hand to the Maintain tier is a `compute_env_id`.

!!! note "Tearing it down"

    Terraform manages the cloud resources it created, so it can remove them too. From `terraform/compute-env`, run `terraform apply -destroy` to delete the compute environment and the Batch pools in one step. The existing Batch account and credentials are referenced, not managed, so they are left untouched.

---

## 2. Maintain: add a pipeline to the Launchpad

**Requires:** Maintain role, and a `compute_env_id` (from section 1 or an Admin on your team) plus the numeric `workspace_id`. Maintain manages pipelines, not compute environments. That split is deliberate: Admin owns the compute environment, you own your pipelines.

You add the same pipeline, `rnaseq-nf-$WORKSHOP_USER`, two ways. The first, `tw`, is imperative: you run a command and it acts. The second, Terraform, is declarative. Watch what happens when you run each one twice.

### 2.1. Add a pipeline via `tw`

The CLI adds a pipeline in one command.

!!! note

    Set the compute-env name to the one you created in section 1.4, or the one your Admin provided. We've left it as `azure-batch-manual` to match the Terraform example above, assuming you won't delete the compute environment in between.

```bash
tw pipelines add \
  --name="rnaseq-nf-$WORKSHOP_USER" \
  --compute-env="azure-batch-manual" \
  --work-dir="az://work" \
  --revision="master" \
  https://github.com/nextflow-io/rnaseq-nf
```

Run it again and `tw` errors: a pipeline with that name already exists. The command is imperative, so each invocation tries to add a pipeline; it has no notion of "already in the desired state". To make "add once, and leave it alone after that" the _default_ behaviour, we need a tool that manages existence rather than actions: Terraform.

The next section adds a pipeline with the same name using Terraform, so delete this one first to start from a clean slate:

```bash
tw pipelines delete --name="rnaseq-nf-$WORKSHOP_USER"
```

### 2.2. Add a pipeline with Terraform

`side-quests/platform_automation/terraform/pipeline` adds the same pipeline declaratively. You describe the pipeline that should exist; Terraform makes the workspace match. Like the imperative `tw` command above, the config tracks the `master` branch; you will pin it to a release tag below to see Terraform update the pipeline in place.

```bash
cd /workspaces/training/side-quests/platform_automation/terraform/pipeline
terraform init
terraform apply -var="username=$WORKSHOP_USER" -var="workspace_id=<ID>" -var="compute_env_id=<CE_ID>"
```

!!! tip

    Copy `terraform.tfvars.example` to `terraform.tfvars` and fill it in so you stop passing `-var` flags.

Confirm it landed, in the UI (open the Launchpad) or from the terminal, the same way you read the compute environment back in section 1.2:

```bash
tw pipelines list
```

`rnaseq-nf-$WORKSHOP_USER` is there, with no run. Now apply again:

```bash
terraform apply -var="username=$WORKSHOP_USER" -var="workspace_id=<ID>" -var="compute_env_id=<CE_ID>"
```

??? success "Command output"

    ```console
    No changes. Your infrastructure matches the configuration.
    ```

That is the difference. `tw` errored on the second run because it performs an action every time. Terraform manages whether the pipeline **exists**: it is already there, so there is nothing to do. Apply ten times, still one pipeline. This is what keeps you from accumulating competing, half-duplicated resources in a shared workspace. To remove the pipeline, `terraform apply -destroy` with the same vars.

**Update in place.** Terraform does more than create-or-skip; it reconciles the pipeline to whatever the config says. The config tracks the `master` branch, which moves as the pipeline gets new commits. Pin it to a fixed release tag instead, so every launch runs the same code. Edit one line: in `main.tf`, line 35 of the `launch` block, change the revision from the `master` branch to the `v2.4` release tag.

=== "After"

    ```hcl title="terraform/pipeline/main.tf" hl_lines="4" linenums="32"
      launch = {
        pipeline = "https://github.com/nextflow-io/rnaseq-nf"
        # The master branch moves; pin a release tag (e.g. v2.4) for reproducible launches.
        revision       = "v2.4"
        compute_env_id = var.compute_env_id
        work_dir       = var.work_dir
      }
    ```

=== "Before"

    ```hcl title="terraform/pipeline/main.tf" hl_lines="4" linenums="32"
      launch = {
        pipeline = "https://github.com/nextflow-io/rnaseq-nf"
        # The master branch moves; pin a release tag (e.g. v2.4) for reproducible launches.
        revision       = "master"
        compute_env_id = var.compute_env_id
        work_dir       = var.work_dir
      }
    ```

Then apply:

??? success "Command output"

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

Two ways to add one pipeline, two mental models. `tw` is imperative: each run is an action, and adding twice is an error. Terraform is declarative: it manages existence, so a second apply is a no-op or an update-in-place. The artifact you hand to the Launch tier is the Launchpad pipeline `rnaseq-nf-$WORKSHOP_USER`.

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

??? success "Command output"

    ```console
      Workflow 2ZXaU1AzEn7Onk submitted at [Organization / Workspace] workspace.

        https://cloud.seqera.io/orgs/Organization/workspaces/Workspace/watch/2ZXaU1AzEn7Onk
    ```

You can monitor the run by clicking the provided URL.

!!! note "Using a params file"

    If you wish to provide parameters to the pipeline, you can do so with the `--params-file` flag.

### 3.3. Use `seqerakit`

`seqerakit` is a command-line tool that wraps `tw`, driving it from YAML instead of flags.
Each YAML file declares the resources or actions to perform: here a `launch:` block lists the pipeline, workspace, and compute environment to run, with `${...}` values pulled from the environment.

`seqerakit` launches a pipeline that already exists on the Launchpad, filling in the `tw launch` command from `seqerakit/launch-rnaseq.yml`. `WORKSHOP_USER` is set from section 0; set the compute environment name, dry run first, then launch:

```bash
cd /workspaces/training/side-quests/platform_automation/seqerakit
export COMPUTE_ENVIRONMENT=<compute environment name>

seqerakit launch-rnaseq.yml --dryrun   # prints the tw command, changes nothing
seqerakit launch-rnaseq.yml            # launches for real
```

The dry run shows the underlying `tw` command. Run it twice and you get two runs. That is the imperative model: do the thing, now.

The file you just ran is a single `launch` block:

```yaml title="seqerakit/launch-rnaseq.yml"
launch:
  - name: "rnaseq-nf-${WORKSHOP_USER}-run"
    workspace: "${TOWER_WORKSPACE_ID}"
    pipeline: "rnaseq-nf-${WORKSHOP_USER}"
    compute-env: "${COMPUTE_ENVIRONMENT}"
```

One advantage of `seqerakit` is that the YAML forms a template of actions to perform.
This buys you two things.
First, you can save it and re-use it later.
Second, you can launch the same pipeline several times by adding more launch blocks to the YAML file, which is useful for launching with different parameters or on different compute environments:

```yaml title="seqerakit/launch-rnaseq-multiple.yml"
launch:
  - name: "rnaseq-nf-${WORKSHOP_USER}-run-1"
    workspace: "${TOWER_WORKSPACE_ID}"
    pipeline: "rnaseq-nf-${WORKSHOP_USER}"
    compute-env: "${COMPUTE_ENVIRONMENT}"
  - name: "rnaseq-nf-${WORKSHOP_USER}-run-2"
    workspace: "${TOWER_WORKSPACE_ID}"
    pipeline: "rnaseq-nf-${WORKSHOP_USER}"
    compute-env: "${COMPUTE_ENVIRONMENT}"
```

Launch both runs now with:

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

- [Building with Seqera AI (CoScientist)](../co_scientist/index.md) side quest uses the same `rnaseq-nf` pipeline and API endpoints, but drives them via AI agents instead of manually.
