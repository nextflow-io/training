# Inputs for the Admin tier compute-environment config. The Platform token is
# NOT here; it comes from TOWER_ACCESS_TOKEN. The Azure credentials are not here
# either: they are already stored in the Platform and referenced by name.

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

# Existing Azure Batch account (referenced via a data source, not created).
variable "batch_account_name" {
  description = "Name of the existing Azure Batch account."
  type        = string
}

variable "batch_account_rg" {
  description = "Resource group of the existing Azure Batch account."
  type        = string
}

# Azure credentials already stored in the Platform (added to the workspace
# beforehand). Referenced by name; the keys never appear in this config.
variable "azure_credential_name" {
  description = "Name of the existing Azure credentials in the Platform workspace."
  type        = string
}

variable "region" {
  description = "Azure region code of the Batch account, e.g. uaenorth."
  type        = string
}

variable "work_dir" {
  description = "Azure Blob Storage work directory, e.g. az://work."
  type        = string
}
