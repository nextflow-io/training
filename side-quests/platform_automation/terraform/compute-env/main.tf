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
    seqera  = { source = "seqeralabs/seqera", version = "0.30.5" }
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
  head_vm_size      = "Standard_D4ds_v5"
  worker_vm_size    = "Standard_E16ds_v5"
  worker_max_nodes  = 8
  worker_max_tasks  = 16
  node_agent_sku_id = "batch.node.ubuntu 22.04"
}

# Keeps pool names unique across re-creates.
resource "random_string" "suffix" {
  length  = 6
  special = false
  upper   = false
}

# Existing Azure resources the customer owns. Referenced, not managed.
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

# Let the managed identity read and write the work directory's storage.
resource "azurerm_role_assignment" "batch_data_contributor" {
  scope                = data.azurerm_batch_account.existing.id
  role_definition_name = "Azure Batch Data Contributor"
  principal_id         = data.azurerm_user_assigned_identity.existing.principal_id
}

# Head pool: runs the Nextflow head job. One node is enough.
resource "azurerm_batch_pool" "head" {
  name                = "rnaseq-head-${random_string.suffix.result}"
  display_name        = "rnaseq head pool"
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

  container_configuration {
    type = "DockerCompatible"
  }

  identity {
    type         = "UserAssigned"
    identity_ids = [data.azurerm_user_assigned_identity.existing.id]
  }

  network_configuration {
    subnet_id                        = data.azurerm_subnet.existing.id
    public_address_provisioning_type = "NoPublicIPAddresses"
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
  max_tasks_per_node  = local.worker_max_tasks

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
    sku       = "2204"
    version   = "latest"
  }

  container_configuration {
    type = "DockerCompatible"
  }

  identity {
    type         = "UserAssigned"
    identity_ids = [data.azurerm_user_assigned_identity.existing.id]
  }

  network_configuration {
    subnet_id                        = data.azurerm_subnet.existing.id
    public_address_provisioning_type = "NoPublicIPAddresses"
  }

  # Terraform manages the pool, not its live node count.
  lifecycle {
    ignore_changes = [auto_scale]
  }
}

# Azure credentials stored in the Platform. Keys are write-only.
resource "seqera_azure_credential" "main" {
  name         = "azure-batch"
  workspace_id = var.workspace_id
  batch_name   = var.batch_account_name
  batch_key    = var.azure_batch_key
  storage_name = var.azure_storage_name
  storage_key  = var.azure_storage_key
}

# Manual Azure Batch compute environment: uses the head pool, routes tasks to
# the worker pool. References to the pool names make Terraform create them first.
resource "seqera_compute_env" "main" {
  workspace_id = var.workspace_id

  compute_env = {
    name           = "azure-batch-manual"
    description    = "Azure Batch CE on Terraform-managed pools"
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

# The Maintain tier (terraform/pipeline) needs this compute_env_id.
output "compute_env_id" {
  description = "ID of the Azure Batch compute environment."
  value       = seqera_compute_env.main.compute_env_id
}
