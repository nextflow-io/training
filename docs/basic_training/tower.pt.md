---
title: Nextflow Tower
description: Comece a usar o Nextflow Tower
---

# Comece a usar o Nextflow Tower

## Conceitos Básicos

Nextflow Tower é o posto de comando centralizado para gerenciamento de dados e pipelines. Ele traz monitoramento, registro e observabilidade para fluxos de trabalho distribuídos e simplifica a implantação de pipelines em qualquer nuvem, cluster ou laptop.

Os principais recursos do Nextflow Tower incluem:

-   O lançamento de pipelines pré-configurados com facilidade.
-   Integração programática para atender às necessidades de uma organização.
-   Publicando pipelines em espaços de trabalho compartilhados.
-   Gerenciamento da infraestrutura necessária para executar análise de dados em escala.

!!! tip

    [Registre-se](https://cloud.tower.nf/) para experimentar o Tower gratuitamente ou solicitar uma [demonstração](https://seqera.io/demo/) para implantações em seu próprio ambiente local ou na nuvem.

## Como usar

Você pode usar o Tower por meio da opção `-with-tower` ao usar o comando Nextflow **run**, por meio da **interface gráfica online** ou da **API**.

### Com o comando Nextflow `run`

Crie uma conta e faça login no Tower.

**1. Crie um novo token**

Você pode acessar seus tokens no menu suspenso **Settings**:

![Criar um token](img/usage_create_token.png)

**2. Dê um nome para seu token**

![Dê um nome para seu token](img/usage_name_token.png)

**3. Salve seu token com segurança**

Copie e guarde seu novo token em um local seguro.

![Salve o seu token](img/usage_token.png)

**4. Exporte seu token**

Uma vez que seu token foi criado, abra um terminal e digite:

```bash
export TOWER_ACCESS_TOKEN=eyxxxxxxxxxxxxxxxQ1ZTE=
export NXF_VER=20.10.0
```

Onde `eyxxxxxxxxxxxxxxxQ1ZTE=` é o token que você acabou de criar.

!!! note

    Verifique seu `nextflow -version`. Os tokens de portador requerem o Nextflow versão 20.10.0 ou posterior e podem ser configurados com o segundo comando mostrado acima. Você pode alterar a versão, se necessário.

Para enviar um pipeline para uma [área de trabalho](https://help.tower.nf/getting-started/workspace/) (workspace) usando a ferramenta de linha de comando Nextflow, inclua o ID da área de trabalho em seu ambiente. Por exemplo:

```bash
export TOWER_WORKSPACE_ID=000000000000000
```

O ID da área de trabalho pode ser encontrado na página de visão geral das áreas de trabalho (workspaces) da organização.

**5. Execute o Nextflow com o Tower**

Execute seus fluxos de trabalho do Nextflow normalmente com a adição do comando `-with-tower`:

```bash
nextflow run hello.nf -with-tower
```

Você verá e poderá monitorar seus **trabalhos do Nextflow** no Tower.

Para configurar e executar tarefas do Nextflow em **ambientes de nuvem**, visite a [seção de ambientes de computação](https://help.tower.nf/compute-envs/overview/) (compute environment).

!!! exercise

    Execute o script de RNA-Seq `script7.nf` usando o sinalizador `-with-tower`, depois de concluir corretamente as configurações de token descritas acima.

    ??? tip

        Vá para <https://tower.nf/>, faça o login, em seguida, clique na guia de execução (Run) e selecione a execução que você acabou de enviar. Se você não conseguir encontrá-lo, verifique novamente se seu token foi inserido corretamente.

### Com uma interface gráfica online

Para executar usando a interface gráfica (GUI), existem três etapas principais:

1. Crie uma conta e faça login no Tower, disponível gratuitamente, em [tower.nf](https://tower.nf).
2. Crie e configure um novo [ambiente de computação](https://help.tower.nf/compute-envs/overview/) (compute environment).
3. Comece a [lançar pipelines](https://help.tower.nf/launch/launchpad/).

#### Configurando seu ambiente de computação

O Tower usa o conceito de **Ambientes de Computação** (compute environment) para definir a plataforma de execução onde um pipeline será executado.

Ele suporta o lançamento de pipelines em um número crescente de infraestruturas de **nuvem** e **on-prem** (infraestrutura dedicada).

![Ambientes de computação](img/compute_env_platforms.png)

Cada ambiente de computação deve ser pré-configurado para permitir que o Tower envie tarefas. Você pode ler mais sobre como configurar cada ambiente usando os links abaixo.

!!! tip "Os guias a seguir descrevem como configurar cada um desses ambientes de computação."

    * [AWS Batch](https://help.tower.nf/compute-envs/aws-batch/)
    * [Azure Batch](https://help.tower.nf/compute-envs/azure-batch/)
    * [Google Cloud](https://help.tower.nf/compute-envs/google-cloud/)
    * [IBM LSF](https://help.tower.nf/compute-envs/lsf/)
    * [Slurm](https://help.tower.nf/compute-envs/slurm/)
    * [Grid Engine](https://help.tower.nf/compute-envs/grid-engine/)
    * [Altair PBS Pro](https://help.tower.nf/compute-envs/altair-pbs-pro/)
    * [Amazon Kubernetes (EKS)](https://help.tower.nf/compute-envs/eks/)
    * [Google Kubernetes (GKE)](https://help.tower.nf/compute-envs/gke/)
    * [Hosted Kubernetes](https://help.tower.nf/compute-envs/k8s/)

#### Selecionando um ambiente de computação padrão

If you have more than one **Compute Environment**, you can select which one will be used by default when launching a pipeline.

1. Navigate to your [compute environments](https://help.tower.nf/compute-envs/overview/).
2. Choose your default environment by selecting the **Make primary** button.

**Congratulations!**

You are now ready to launch pipelines with your primary compute environment.

#### Launchpad

Launchpad makes it easy for any workspace user to launch a pre-configured pipeline.

![Launchpad](img/overview_launch.png)

A pipeline is a repository containing a Nextflow workflow, a compute environment and pipeline parameters.

#### Formulário de Parâmetros de Pipeline

Launchpad automatically detects the presence of a `nextflow_schema.json` in the root of the repository and dynamically creates a form where users can easily update the parameters.

!!! info

    The parameter forms view will appear if the workflow has a Nextflow schema file for the parameters. Please refer to the [Nextflow Schema guide](https://help.tower.nf/pipeline-schema/overview) to learn more about the schema file use-cases and how to create them.

This makes it trivial for users without any expertise in Nextflow to enter their pipeline parameters and launch.

![Pipeline parameters](img/launch_rnaseq_nextflow_schema.png)

#### Adicionando um novo pipeline

Adding a pipeline to the pre-saved workspace launchpad is detailed in full on the [tower webpage docs](https://help.tower.nf/launch/launch/).

In brief, these are the steps you need to follow to set up a pipeline.

1. Select the Launchpad button in the navigation bar. This will open the **Launch Form**.
2. Select a [compute environment](https://help.tower.nf/compute-envs/overview).
3. Enter the repository of the pipeline you want to launch. e.g. <https://github.com/nf-core/rnaseq.git>
4. Select a pipeline **Revision number**. The Git default branch (main/master) or `manifest.defaultBranch` in the Nextflow configuration will be used by default.
5. Set the **Work directory** location of the Nextflow work directory. The location associated with the compute environment will be selected by default.
6. Enter the name(s) of each of the Nextflow **Config profiles** followed by the `Enter` key. See the Nextflow [Config profiles](https://www.nextflow.io/docs/latest/config.html#config-profiles) documentation for more details.
7. Enter any Pipeline parameters in YAML or JSON format. YAML example:

    ```yaml
    reads: "s3://nf-bucket/exome-data/ERR013140_{1,2}.fastq.bz2"
    paired_end: true
    ```

8. Select Launchpad to begin the pipeline execution.

!!! info

    Nextflow pipelines are simply Git repositories and can be changed to any public or private Git-hosting platform. See Git Integration in the Tower docs and Pipeline Sharing in the Nextflow docs for more details.

!!! note

    The credentials associated with the compute environment must be able to access the work directory.

!!! info

    In the configuration, the full path to a bucket must be specified with single quotes around strings and no quotes around booleans or numbers.

!!! tip

    To create your own customized Nextflow Schema for your pipeline, see the examples from the `nf-core` workflows that have adopted this approach. For example, [eager](https://github.com/nf-core/eager/blob/2.3.3/nextflow_schema.json) and [rnaseq](https://github.com/nf-core/rnaseq/blob/3.0/nextflow_schema.json).

For advanced settings options check out this [page](https://help.tower.nf/launch/advanced/).

There is also community support available if you get into trouble, join the Nextflow Slack by following this [link](https://www.nextflow.io/slack-invite.html).

### API

To learn more about using the Tower API, visit the [API section](https://help.tower.nf/api/overview/) in this documentation.

## Áreas de trabalho e Organizações

Nextflow Tower simplifies the development and execution of workflows by providing a centralized interface for users and organizations.

Each user has a unique **workspace** where they can interact and manage all resources such as workflows, compute environments and credentials. Details of this can be found [here](https://help.tower.nf/getting-started/workspace/).

By default, each user has their own private workspace, while organizations have the ability to run and manage users through role-based access as **members** and **collaborators**.

### Recursos de organização

You can create your own organization and participant workspace by following the docs at [tower](https://help.tower.nf/orgs-and-teams/workspace-management/).

Tower allows the creation of multiple organizations, each of which can contain multiple workspaces with shared users and resources. This allows any organization to customize and organize the usage of resources while maintaining an access control layer for users associated with a workspace.

### Usuários da organização

Any user can be added or removed from a particular organization or a workspace and can be allocated a specific access role within that workspace.

The Teams feature provides a way for organizations to group various users and participants together into teams. For example, `workflow-developers` or `analysts`, and apply access control to all the users within this team collectively.

For further information, please refer to the [User Management](https://help.tower.nf/orgs-and-teams/organizations/) section.

#### Configurando uma nova organização

Organizations are the top-level structure and contain Workspaces, Members, Teams and Collaborators.

To create a new Organization:

1.  Click on the dropdown next to your name and select New organization to open the creation dialog.
2.  On the dialog, fill in the fields as per your organization. The Name and Full name fields are compulsory.

    !!! note

        A valid name for the organization must follow a specific pattern. Please refer to the UI for further instructions.

3.  The rest of the fields such as Description, Location, Website URL and Logo Url are optional.
4.  Once the details are filled in, you can access the newly created organization using the organization’s page, which lists all of your organizations.

    !!! note

        It is possible to change the values of the optional fields either using the Edit option on the organization’s page or by using the Settings tab within the organization page, provided that you are the Owner of the organization.

    !!! tip

        A list of all the included Members, Teams and Collaborators can be found on the organization page.
