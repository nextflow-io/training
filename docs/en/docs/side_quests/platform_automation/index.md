# Seqera Platform Automation

The Seqera Platform does not run your work. It is an API and a control plane over your cloud. It hands jobs to a compute environment, the compute environment runs them on cloud VMs, and the Platform reads back state, logs, and exit codes.

Everything the web UI does, it does by calling the Platform API. So everything is automatable: with the API, Terraform, and the CLI you can manage compute environments, pipelines, and runs as code, with no clicking. This side quest walks that programmatic surface from the most privileged role to the least.

### Learning goals

In this side quest, we'll drive the Platform across three workspace roles of decreasing permission: **Admin**, **Maintain**, and **Launch**. Each role does one job and hands an artifact to the next. You'll learn how to:

- Create a compute environment two ways, with increasing control over the cloud: the UI (Batch Forge) and Terraform
- Add a pipeline to the Launchpad declaratively with Terraform, and see idempotency
- Launch a pipeline using the GUI, CLI and `seqerakit`.
- Tell declarative existence (Terraform) apart from imperative actions (the UI, `seqerakit`, `tw`) and learn when to use each

### Prerequisites

Before taking on this side quest, you should:

- Have a Seqera Platform account on Seqera Cloud (`https://cloud.seqera.io`), or an Enterprise install
- Be a member of a workspace with a role: Admin, Maintain, or Launch
- Have credentials for your cloud provider already added to that workspace (not covered here)
- Ideally, a GitHub access token added too, to avoid GitHub rate limits
- Be comfortable with the command line and basic Nextflow concepts

New to running pipelines on Seqera at all? Start with the gentler "Run pipelines on Seqera" module in the [Nextflow Triathlon](https://training.nextflow.io/) (sign up, launch in the UI, the `tw` CLI). This side quest is the automation layer on top of it: roles, Terraform, `seqerakit`, and Actions.

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
    │   │   ├── terraform.tfvars.example
    │   │   └── README.md
    │   └── pipeline           # section 2 (Maintain): add a pipeline to the Launchpad
    │       ├── main.tf
    │       ├── variables.tf
    │       ├── terraform.tfvars.example
    │       └── README.md
    └── seqerakit              # sections 2 & 3: add and launch a pipeline
        ├── add-rnaseq.yml
        ├── launch-rnaseq.yml
        └── README.md
    ```

### Create an access token

We authenticate to the Platform to tell it who we are. That is what grants the permissions of our role. Seqera uses a **personal access token**: a bearer token you create once and send with every API call, in an `Authorization: Bearer <token>` header. The token carries your identity and role, so the Platform knows who you are and what you can do. Anyone with the token can act as you, so keep it secret. In this training we save it as an environment variable.

1. In the Platform, open the user menu (top right) and choose **Your tokens**.
2. Click **Add token**, name it (e.g. `platform-automation`), and click **Add**.
3. Copy the token now. The Platform shows it only once.
4. In the Codespace terminal, export it under both names:

```bash
export TOWER_ACCESS_TOKEN=<paste-token>
export SEQERA_ACCESS_TOKEN=$TOWER_ACCESS_TOKEN
```

Terraform, `tw`, and `seqerakit` all read the token from these variables. Check it worked:

```bash
tw info
```

This prints the API endpoint and the authenticated user. If it errors, the token is wrong or not exported.

!!! warning

    This only exists in a single terminal session, if you open a new terminal you will need to export them again.

For Seqera Enterprise, also set `SEQERA_API_URL` (and the `server_url` Terraform variable) to your install's API URL. Everything else is identical.

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

#### Know your role

The work is split by role, from most to least privileged. Start at the highest tier your access allows:

| Tier                                                    | Role          | Can change                               | Produces             |
| ------------------------------------------------------- | ------------- | ---------------------------------------- | -------------------- |
| [Admin](#1-admin-compute-environments-and-cloud)        | Owner / Admin | Compute environments and cloud resources | a `compute_env_id`   |
| [Maintain](#2-maintain-add-a-pipeline-to-the-launchpad) | Maintain      | Pipelines on the Launchpad               | a Launchpad pipeline |
| [Launch](#3-launch-run-a-pipeline)                      | Launch        | Nothing; can only run                    | pipeline runs        |

Admin users create compute environments and assign roles to everyone else. If you only have Maintain, an Admin on your team hands you a `compute_env_id`; start at section 2. If you only have Launch, they hand you a Launchpad pipeline (and optionally an Action) to run; start at section 3. Each section opens with what it requires.

---

## 1. Admin: compute environments and cloud

**Requires:** Owner or Admin role, and permissions in the cloud environment. If you are using a Cloud provider, authenticate with them first via their CLI.

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

### 1.2. Provision the cloud with Terraform

Terraform manages cloud resources declaratively: you describe what should exist and Terraform makes it so. It does not run your work; it only creates the compute environment. `side-quests/platform_automation/terraform/compute-env/` does it in a single `main.tf`, two providers in one apply:

- data sources: the existing Batch account, managed identity, and vnet/subnet,
  referenced, not created.
- `azurerm`: creates a head pool and a worker pool on that account and subnet,
  as that identity.
- `seqera`: stores the Azure credentials and creates the compute environment.

The manual marker: `head_pool` is set to the pool Terraform just made and there is **no `forge` block**. That one difference is what makes the compute environment manual instead of Forge. `nextflow_config` routes tasks to the worker pool.

Let's walk `terraform/compute-env/main.tf` from top to bottom.

The first block declares the providers and pins their versions. The `seqera` provider is pinned to exactly `0.30.5`:

```terraform
terraform {
  required_version = ">= 1.9"

  required_providers {
    seqera  = { source = "seqeralabs/seqera", version = "0.30.5" }
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

Pool sizing and the VM image live in a `locals` block, so you edit them in one place rather than wiring every value through as a variable:

```terraform
locals {
  head_vm_size      = "Standard_D4ds_v5"
  worker_vm_size    = "Standard_E16ds_v5"
  worker_max_nodes  = 8
  worker_max_tasks  = 16
  node_agent_sku_id = "batch.node.ubuntu 22.04"
}
```

The cloud resources the lab already owns, the Batch account, the managed identity, and the subnet, are referenced with data sources. Terraform reads them; it does not create or manage them:

```terraform
data "azurerm_batch_account" "existing" {
  name                = var.batch_account_name
  resource_group_name = var.batch_account_rg
}

data "azurerm_user_assigned_identity" "existing" {
  name                = var.managed_identity_name
  resource_group_name = var.managed_identity_rg
}

data "azurerm_subnet" "existing" {
  name                 = var.subnet_name
  virtual_network_name = var.vnet_name
  resource_group_name  = var.vnet_rg
}
```

Now the resources Terraform does create. The head pool runs the Nextflow head job, so one fixed node is enough:

```terraform
resource "azurerm_batch_pool" "head" {
  name                = "rnaseq-head-${random_string.suffix.result}"
  resource_group_name = var.batch_account_rg
  account_name        = var.batch_account_name
  vm_size             = local.head_vm_size
  node_agent_sku_id   = local.node_agent_sku_id

  fixed_scale {
    target_dedicated_nodes = 1
  }

  storage_image_reference {
    publisher = "microsoft-dsvm"
    offer     = "ubuntu-hpc"
    sku       = "2204"
    version   = "latest"
  }
  # identity, container, and network_configuration omitted for brevity
}
```

The worker pool runs the pipeline tasks, so it autoscales from zero up to `worker_max_nodes` based on the number of pending tasks:

```terraform
resource "azurerm_batch_pool" "worker" {
  name               = "rnaseq-worker-${random_string.suffix.result}"
  vm_size            = local.worker_vm_size
  node_agent_sku_id  = local.node_agent_sku_id
  max_tasks_per_node = local.worker_max_tasks
  # ...same account, image, identity, and network as the head pool

  auto_scale {
    evaluation_interval = "PT5M"
    formula             = <<-FORMULA
      pending = avg($PendingTasks.GetSample(180 * TimeInterval_Second));
      $TargetDedicatedNodes = min(${local.worker_max_nodes}, pending);
      $NodeDeallocationOption = taskcompletion;
    FORMULA
  }
}
```

Both pools run the same `ubuntu-hpc` `2204` image, which is why `node_agent_sku_id` is `batch.node.ubuntu 22.04` for both.

The `seqera` provider stores the Azure credentials in the Platform:

```terraform
resource "seqera_azure_credential" "main" {
  name         = "azure-batch"
  workspace_id = var.workspace_id
  batch_name   = var.batch_account_name
  batch_key    = var.azure_batch_key
  storage_name = var.azure_storage_name
  storage_key  = var.azure_storage_key
}
```

Finally, the compute environment itself. This is the manual marker in code: `head_pool` points at the pool we just created, there is no `forge` block, and `nextflow_config` routes tasks to the worker pool's queue. Because it references `azurerm_batch_pool.head.name` and `azurerm_batch_pool.worker.name`, Terraform knows to create the pools first:

```terraform
resource "seqera_compute_env" "main" {
  workspace_id = var.workspace_id

  compute_env = {
    name           = "azure-batch-manual"
    platform       = "azure-batch"
    credentials_id = seqera_azure_credential.main.credentials_id

    config = {
      azure_batch = {
        region                     = var.region
        work_dir                   = var.work_dir
        head_pool                  = azurerm_batch_pool.head.name
        managed_identity_client_id = data.azurerm_user_assigned_identity.existing.client_id
        nextflow_config            = "process.queue = '${azurerm_batch_pool.worker.name}'\n"
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

Terraform reads those references and works out the order itself: credentials and pools first, then the compute environment that depends on them. Run it:

```bash
az login
cd /workspaces/training/side-quests/platform_automation/terraform/compute-env
terraform init
terraform plan
terraform apply
```

See `side-quests/platform_automation/terraform/compute-env/README.md` for the full variable list and the `terraform.tfvars` setup.

### 1.3. See both sides

The compute environment now exists in two places, and as an Admin with cloud access you can see both.

On the Platform side, read the ID Terraform produced:

```bash
terraform output compute_env_id
```

On the Azure side, open the Batch account in the portal. You will see the head and worker pools Terraform created, sitting idle with no jobs yet. Once you launch a pipeline (sections 2 and 3), jobs and tasks stack onto these pools, and you can drill into a task's logs and exit code. That is the whole point of the Admin tier: the Platform submits to Azure Batch, and Azure Batch runs the work. Maintain- and Launch-role users see only the compute environment in the workspace, not the cloud behind it.

### Takeaway

One compute environment, two levels of control: Forge in the UI (easiest, the Platform owns the pool) and Terraform (you own the pools and the wiring, end to end). The artifact you hand to the Maintain tier is a `compute_env_id`.

---

## 2. Maintain: add a pipeline to the Launchpad

**Requires:** Maintain role, and a `compute_env_id` (from section 1 or an Admin on your team) plus the numeric `workspace_id`. Maintain manages pipelines, not compute environments. That split is deliberate: Admin owns the compute environment, you own your pipelines.

You add the same pipeline, `rnaseq-nf-$USER`, four ways. The first three are imperative, you run a command and it acts. The last, Terraform, is declarative. Watch what happens when you run each one twice.

### 2.1. Add a pipeline via the UI

In the workspace, open the **Launchpad** and click **Add pipeline**. Fill in the form:

- **Name**: `rnaseq-nf-<your-username>` (unique in a shared workspace).
- **Compute environment**: the one from section 1, e.g. `azure-batch-manual`.
- **Pipeline to launch**: `https://github.com/nextflow-io/rnaseq-nf`.
- **Revision**: `master`.
- **Work directory**: your Azure Blob work dir, e.g. `az://nf-work/work`.

Click **Add**. The pipeline appears on the Launchpad with no run started. Adding a pipeline only saves a launch configuration; it does not run anything.

### 2.2. Add a pipeline via `tw`

The CLI does the same thing in one command:

```bash
tw pipelines add \
  --name="rnaseq-nf-$USER" \
  --compute-env="azure-batch-manual" \
  --work-dir="az://nf-work/work" \
  --revision="master" \
  https://github.com/nextflow-io/rnaseq-nf
```

Run it again and `tw` errors: a pipeline with that name already exists. The command is imperative, so each invocation tries to add a pipeline; it has no notion of "already in the desired state".

### 2.3. Add a pipeline with `seqerakit`

`seqerakit` is a wrapper over `tw` that reads a YAML file and runs the underlying `tw` commands. It keeps the configuration as code, so a teammate can reproduce the exact same pipeline. `seqerakit/add-rnaseq.yml` describes the pipeline:

```yaml
pipelines:
  - name: "rnaseq-nf-${USERNAME}"
    url: "https://github.com/nextflow-io/rnaseq-nf"
    workspace: "my-org/platform-automation"
    description: "rnaseq-nf for ${USERNAME}, added with seqerakit"
    compute-env: "azure-batch-manual"
    work-dir: "az://nf-work/work"
    revision: "master"
```

`seqerakit` expands `${USERNAME}` from the environment, so set it first, then add the pipeline:

```bash
cd /workspaces/training/side-quests/platform_automation/seqerakit
export USERNAME=$USER
seqerakit add-rnaseq.yml
```

If you run it again, it errors, the same way `tw` did, because the pipeline already exists:

```console
ERROR: A pipeline with name 'rnaseq-nf-<username>' already exists.
```

You can force it through with `seqerakit add-rnaseq.yml --overwrite`, which deletes and recreates the pipeline. But that is you telling it to repeat the action. By default, adding twice is an error. To make "add once, and leave it alone after that" the _default_ behaviour, we need a tool that manages existence rather than actions: Terraform.

### 2.4. Add a pipeline with Terraform

`side-quests/platform_automation/terraform/pipeline` adds the same pipeline declaratively. You describe the pipeline that should exist; Terraform makes the workspace match:

```bash
cd /workspaces/training/side-quests/platform_automation/terraform/pipeline
terraform init
terraform plan  -var="username=$USER" -var="workspace_id=<ID>" -var="compute_env_id=<CE_ID>"
terraform apply -var="username=$USER" -var="workspace_id=<ID>" -var="compute_env_id=<CE_ID>"
```

Tip: copy `terraform.tfvars.example` to `terraform.tfvars` so you stop passing `-var` flags.

Open the Launchpad: `rnaseq-nf-$USER` is there, with no run. Now apply again:

```bash
terraform apply -var="username=$USER" -var="workspace_id=<ID>" -var="compute_env_id=<CE_ID>"
```

```console
No changes. Your infrastructure matches the configuration.
```

That is the difference. `tw` and `seqerakit` errored on the second run because they perform an action every time. Terraform manages whether the pipeline **exists**: it is already there, so there is nothing to do. Apply ten times, still one pipeline. This is what keeps you from accumulating competing, half-duplicated resources in a shared workspace. To remove the pipeline, `terraform apply -destroy` with the same vars.

### Takeaway

Four ways to add one pipeline, two mental models. `tw` and `seqerakit` are imperative: each run is an action, and adding twice is an error. Terraform is declarative: it manages existence, so a second apply is a no-op. The artifact you hand to the Launch tier is the Launchpad pipeline `rnaseq-nf-$USER`.

---

## 3. Launch: run a pipeline

**Requires:** Launch role and a token with that role, plus a pipeline on the Launchpad (added by a Maintainer in section 2). Launch can only run things; it cannot create or modify compute environments, pipelines, or Actions.

Everything in this section runs the pipeline the Maintainer already configured. None of it changes the pipeline; launching is an action, not a state.

### 3.1. Use the GUI

Open the **Launchpad**, select `rnaseq-nf-$USER`, and click **Launch**. The compute environment, revision, and work directory are already filled in by the Maintainer, so a Launch user just clicks the button. Submit it and the run appears under **Runs**.

### 3.2. Use the `tw` CLI

```bash
tw launch rnaseq-nf-$USER
```

One command, no flags: everything is pre-configured on the Launchpad entry. The CLI submits the run and prints its URL.

For an automated launch, pin the parameters in a version-controlled file and pass it with `--params-file params.yaml`. That keeps each run reproducible, which is the whole reason to drive launches from code rather than the form.

### 3.3. Use `seqerakit`

`seqerakit` launches a pipeline that already exists on the Launchpad, filling in the `tw launch` command from `seqerakit/launch-rnaseq.yml`. Set `USERNAME`, dry run first, then launch:

```bash
cd /workspaces/training/side-quests/platform_automation/seqerakit
export USERNAME=$USER

seqerakit launch-rnaseq.yml --dryrun   # prints the tw command, changes nothing
seqerakit launch-rnaseq.yml            # launches for real
```

The dry run shows the underlying `tw` command. Run it twice and you get two runs. That is the imperative model: do the thing, now. It is the opposite of Terraform (section 2.4), which manages existence. See `side-quests/platform_automation/seqerakit/README.md`.

### 3.4. Examine the cloud resources

A run does not create one neatly named cloud job. Nextflow submits **many** Batch jobs and tasks, one per process invocation, so you will not find a single resource named after the Platform run. Open the run on the Platform (**Runs** → your run), which mirrors the underlying Batch service: the task table here is the same work you can see in the cloud console.

If you have cloud access, the prefixes Nextflow uses in each Batch service are:

- **Azure Batch**: jobs named `job-<hex>-<process>`, tasks named `nf-<taskhash>`.
- **AWS Batch**: job names of the form `<run-name>_<process>`, with unsupported characters stripped (max 128 chars).
- **GCP Batch**: job IDs of the form `nf-<taskhash>-<timestamp>`.

The relationship is the point: the Platform hands work to Batch, Batch runs it on the pool VMs, and the Platform reads back state, logs, and exit codes.

### 3.5. Side note: launch Actions

There is one more way to launch, built for automation rather than people. An **Action** is a saved launch configuration behind a single URL: call that URL and the pipeline runs, with no Launchpad and no `tw`. A Maintainer creates one in the UI under **Actions** → **Add action**, picks the pipeline and compute environment, and saves it. The Platform then shows a ready-made `curl` command for the Action's endpoint.

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
