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

Most Seqera compute environments require provisioning a [credential](https://docs.seqera.io/platform/latest/credentials/overview) that grants access to those compute resources.
In this case, we will use [Tower Agent ](https://docs.seqera.io/platform/24.1/supported_software/agent), a lightweight program that can be installed on any machine to enable Seqera Platform to run Nextflow workflows on your behalf.

Follow these instructions to configure your Tower Agent.

    1. Create a new Token named "GitpodAgentToken" on the tokens page, following previous instructions.
    1. Export your token into your current terminal:

        ```bash
         export TOWER_ACCESS_TOKEN=<your-token>
        ```

    1. Within Seqera Platform, click on your workspace name on the top left, and change back to your personal user workspace.
    1. Click on the "Credentials" tab
    1. Add the name "gitpodTowerAgent" to the name field
    1. Select "Tower Agent" from the list of providers.
    1. Copy the "Agent Connection ID" shown in the dropdown, it should look similar to `75d74f5f-9454-48b6-8967-cf20b74f6c78`
    1. In your terminal execute the command below replacing with your connection ID:

        ```bash
        tw-agent 75d74f5f-9454-48b6-8967-cf20b74f6c78 --work-dir=./work
        ```

    1. Return to Seqera Platform and click on the "Add" button.

If you completed this successfully, you'll see "gitpodTowerAgent" in the list of credentials in Seqera Platform, and console output in your terminal similar to what's below:

```console title="tw-agent logs"
21:47:33.531 INFO - Established active environments: [cli]
21:47:33.662 INFO - TOWER AGENT v0.5.0
21:47:33.662 INFO - Compatible with TOWER API v1.8
21:47:33.662 INFO - Connecting as user 'gitpod' with default work directory '/workspaces/training/work'
21:47:34.565 INFO - Connecting to Tower
21:47:34.801 INFO - Connection to Tower established
21:48:18.674 INFO - Sending heartbeat
21:48:18.755 INFO - Received heartbeat
21:49:03.680 INFO - Sending heartbeat
21:49:03.741 INFO - Received heartbeat
```

### 3.2. Add a _simulated_ Grid Engine compute environment to Seqera Plaform.

Now that we have set up a credential granting us secure access to our computational resources, we will set up a Compute Environment in Seqera which is a set of configuration that allows us to launch nextflow pipelines.
Follow the steps below:

    1. Navigate to the Compute Environments tab of Seqera Platform.
    1. Click "Add compute environment"
    1. Name your compute environment `gitpodGridEngine`.
    1. Select "Grid Engine" as the target execution platform.
    1. Make sure "gitpodTowerAgent" is selected from the list of credentials.
    1. Click "Add"

### 3.3. Add a pipeline to the launchpad for execution on our "Compute Environment"

Now we need to set up a pipeline in Seqera Platform to run on our simulated Grid Engine compute environment.

Follow the steps below:

    1. Click on the "Launchpad" tab of Seqera Platform.
    1. On the top right click on "Add Pipeline".
    1. For name enter `hello-grid-engine`
    1. For the Compute Environment select `gitpodGridEngine`
    1. For the "Pipeline to launch" enter `nextflow-io/hello`.
    1. Finally click "Add"

You now should see "hello-grid-engine` in the list of pipelines.

### 3.4. Launch your Nextflow pipeline

To create that:
