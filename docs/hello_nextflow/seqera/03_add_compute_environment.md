## 3. Add a compute environment

In this section, we will simulate setting up a new compute environment to run our pipeline in Seqera Platform Launchpad.
Seqera Platform uses the concept of **Compute Environments** to define the execution platform where a workflow will run.
It supports the launching of workflows into a growing number of **cloud** and **on-premise** infrastructures.

![Compute environments](seqera/img/compute_env_platforms.png)

Each compute environment must be pre-configured to enable Seqera Platform to submit tasks. You can read more on how to set up each environment using the links below.

!!! tip "The following guides describe how to configure each of these compute environments."

    * [AWS Batch](https://docs.seqera.io/platform/24.1/compute-envs/aws-batch)
    * [Azure Batch](https://docs.seqera.io/platform/24.1/compute-envs/azure-batch)
    * [Google Batch](https://docs.seqera.io/platform/24.1/compute-envs/google-cloud-batch)
    * [Google Life Sciences](https://docs.seqera.io/platform/24.1/compute-envs/google-cloud-lifesciences)
    * [HPC (LSF, Slurm, Grid Engine, Altair PBS Pro)](https://docs.seqera.io/platform/24.1/compute-envs/hpc)
    * [Amazon Kubernetes (EKS)](https://docs.seqera.io/platform/24.1/compute-envs/eks)
    * [Google Kubernetes (GKE)](https://docs.seqera.io/platform/24.1/compute-envs/gke)
    * [Hosted Kubernetes](https://docs.seqera.io/platform/24.1/compute-envs/k8s)

To practice this process, we will simulate setting up a new slurm compute environment on our gitpod using Tower Agent.

### 3.1. Setup Tower Agent

### 3.2. Add a new Slurm compute environment
