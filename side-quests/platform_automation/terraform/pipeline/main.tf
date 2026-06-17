# Maintain tier: add one rnaseq-nf pipeline to the Launchpad. This ADDS the
# pipeline; it does not launch it. Terraform manages existence, not runs:
#
#   terraform apply           creates rnaseq-nf-<username>
#   terraform apply (again)   no-op; idempotent, never launches a run
#   terraform apply -destroy  removes it
#
# Launching is imperative, done with seqerakit/tw (see ../../seqerakit).
#
# Auth: export TOWER_ACCESS_TOKEN; the seqera provider reads it from the env.

terraform {
  required_version = ">= 1.6"

  required_providers {
    seqera = {
      source  = "seqeralabs/seqera"
      version = "0.40.1"
    }
  }
}

provider "seqera" {
  server_url = var.server_url
}

resource "seqera_pipeline" "rnaseq_nf" {
  name         = "rnaseq-nf-${var.username}"
  description  = "rnaseq-nf for ${var.username}, added with Terraform"
  workspace_id = var.workspace_id

  launch = {
    pipeline = "https://github.com/nextflow-io/rnaseq-nf"
    # The master branch moves; pin a release tag (e.g. v2.4) for reproducible launches.
    revision       = "master"
    compute_env_id = var.compute_env_id
    work_dir       = var.work_dir
  }
}

output "pipeline_name" {
  description = "Launchpad pipeline name. Use this when you launch it."
  value       = seqera_pipeline.rnaseq_nf.name
}
