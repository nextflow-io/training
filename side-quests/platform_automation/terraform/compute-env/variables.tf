# Inputs for the Admin tier compute-environment config. The Platform token is
# NOT here; it comes from TOWER_ACCESS_TOKEN. Sensitive Azure keys come from a
# gitignored terraform.tfvars or TF_VAR_ environment variables.

variable "server_url" {
  description = "Seqera Platform API endpoint. Cloud by default."
  type        = string
  default     = "https://api.cloud.seqera.io"
}

variable "subscription_id" {
  description = "Azure subscription ID."
  type        = string
}

variable "workspace_id" {
  description = "Numeric ID of the Platform workspace that owns the compute environment."
  type        = number
}

# Existing Azure resources (referenced via data sources, not created).
variable "batch_account_name" {
  description = "Name of the existing Azure Batch account."
  type        = string
}

variable "batch_account_rg" {
  description = "Resource group of the existing Azure Batch account."
  type        = string
}

variable "managed_identity_name" {
  description = "Name of the existing user-assigned managed identity for Batch nodes."
  type        = string
}

variable "managed_identity_rg" {
  description = "Resource group of the existing managed identity."
  type        = string
}

variable "vnet_name" {
  description = "Name of the existing virtual network the pools attach to."
  type        = string
}

variable "vnet_rg" {
  description = "Resource group of the existing virtual network."
  type        = string
}

variable "subnet_name" {
  description = "Name of the existing subnet the pool nodes join."
  type        = string
}

variable "region" {
  description = "Azure region code of the Batch account, e.g. uaenorth."
  type        = string
}

variable "work_dir" {
  description = "Azure Blob Storage work directory, e.g. az://nf-work/work."
  type        = string
}

# Seqera credential keys (sensitive).
variable "azure_batch_key" {
  description = "Azure Batch account access key."
  type        = string
  sensitive   = true
}

variable "azure_storage_name" {
  description = "Azure Storage account name for the work directory."
  type        = string
}

variable "azure_storage_key" {
  description = "Azure Storage account access key."
  type        = string
  sensitive   = true
}
