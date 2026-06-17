# Inputs for the Maintain tier pipeline config.

variable "server_url" {
  description = "Seqera Platform API endpoint. Cloud by default."
  type        = string
  default     = "https://api.cloud.seqera.io"
}

# Your workshop handle makes the pipeline name unique in a shared workspace. No
# default on purpose, so you must set it. Set it with -var or in terraform.tfvars:
#
#   terraform apply -var="username=$WORKSHOP_USER"
variable "username" {
  description = "Your workshop handle. Makes the Launchpad pipeline name unique."
  type        = string

  # No dots: Launchpad pipeline names reject them.
  validation {
    condition     = can(regex("^[a-zA-Z0-9_-]{2,80}$", var.username))
    error_message = "username must be 2-80 chars: letters, numbers, dash, underscore (no dots)."
  }
}

variable "workspace_id" {
  description = "Numeric ID of the shared workshop workspace."
  type        = number
}

# The Admin tier built the shared compute environment with terraform/compute-env
# and shared its ID. This is a string ID, not the numeric workspace ID.
variable "compute_env_id" {
  description = "ID of the shared Azure Batch compute environment to launch against."
  type        = string
}

variable "work_dir" {
  description = "Azure Blob Storage work directory for runs of this pipeline."
  type        = string
  default     = "az://work"
}
