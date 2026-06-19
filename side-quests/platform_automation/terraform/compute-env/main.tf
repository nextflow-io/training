# Admin tier: Terraform creates the Azure Batch pools, then wires a manual
# Seqera compute environment to them. Manual = head_pool is set and there is no
# forge block, so the Platform submits to these pools but never creates or
# scales them. One config, two providers, end to end.
#
# Auth: run `az login`, then export TOWER_ACCESS_TOKEN (the seqera provider
# reads it from the environment). Sensitive keys come from terraform.tfvars or
# TF_VAR_ env vars; never commit them.

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

# Pool sizing and image. Edit here rather than exposing every value as a knob.
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

# Keeps pool names unique across re-creates.
resource "random_string" "suffix" {
  length  = 6
  special = false
  upper   = false
}

# Manual pools (head_pool set, no forge block) don't receive the node start task
# that Forge adds to auto-created pools, so we replicate its essential parts:
#   1. Install azcopy on the shared node path. Nextflow's Azure Batch executor
#      uses it to stage data; manual pools must provide it.
#   2. Register the AppArmor profile Fusion needs on Ubuntu 24.04+ nodes. Without
#      it, container init fails with "unable to apply apparmor profile"
#      (COMP-1248, moby/moby#50013, launchpad #2111105).
locals {
  node_start_task_script = <<-SCRIPT
    set -euo pipefail

    # azcopy: Nextflow's Azure Batch executor uses it to transfer files.
    curl -sL https://aka.ms/downloadazcopy-v10-linux -o /tmp/azcopy.tgz
    tar -xzf /tmp/azcopy.tgz --strip-components=1 -C /tmp
    mkdir -p "$AZ_BATCH_NODE_SHARED_DIR/bin/"
    cp /tmp/azcopy "$AZ_BATCH_NODE_SHARED_DIR/bin/"

    # AppArmor profile for Fusion containers on Ubuntu 24.04+.
    mkdir -p /etc/apparmor.d/containers
    printf '%s\n' \
      'abi <abi/4.0>,' \
      'include <tunables/global>' \
      '' \
      'profile seqera-fusionfs-container flags=(default_allow) {' \
      '  userns,' \
      '  mount fstype=fuse.fusion -> /fusion/,' \
      '  mount fstype=fuse.fusion -> /fusion/**,' \
      '  umount,' \
      '  include <abstractions/base>' \
      '  include <abstractions/nameservice>' \
      '  include if exists <local/seqera-fusionfs-container>' \
      '}' > /etc/apparmor.d/containers/seqera-fusionfs-container
    apparmor_parser -r /etc/apparmor.d/containers/seqera-fusionfs-container
  SCRIPT

  # Batch parses the command line itself (no shell), so multi-line scripts with
  # quoting don't survive. Base64-encode the script and decode it on the node.
  pool_start_command = "/bin/bash -c 'echo ${base64encode(local.node_start_task_script)} | base64 -d | bash'"
}

# Existing Azure resources the customer owns. Referenced, not managed.
data "azurerm_batch_account" "existing" {
  name                = var.batch_account_name
  resource_group_name = var.batch_account_rg
}

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

  # Install azcopy and load the Fusion AppArmor profile on each node at startup.
  start_task {
    command_line       = local.pool_start_command
    task_retry_maximum = 1
    wait_for_success   = true

    user_identity {
      auto_user {
        elevation_level = "Admin"
        scope           = "Pool"
      }
    }
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

  # Install azcopy and load the Fusion AppArmor profile on each node at startup.
  start_task {
    command_line       = local.pool_start_command
    task_retry_maximum = 1
    wait_for_success   = true

    user_identity {
      auto_user {
        elevation_level = "Admin"
        scope           = "Pool"
      }
    }
  }

  # Terraform manages the pool, not its live scale.
  lifecycle {
    ignore_changes = [auto_scale, fixed_scale]
  }
}

# Azure credentials already stored in the Platform. We look them up by name
# rather than creating them, so no keys appear in this config.
data "seqera_credentials" "all" {
  workspace_id = var.workspace_id
}

locals {
  azure_credentials_id = one([
    for c in data.seqera_credentials.all.credentials : c.id
    if c.name == var.azure_credential_name
  ])
}

# Manual Azure Batch compute environment: uses the head pool, routes tasks to
# the worker pool. References to the pool names make Terraform create them first.
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

# The Maintain tier (terraform/pipeline) needs this compute_env_id.
output "compute_env_id" {
  description = "ID of the Azure Batch compute environment."
  value       = seqera_compute_env.main.compute_env_id
}
