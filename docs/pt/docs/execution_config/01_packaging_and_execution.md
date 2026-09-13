# Parte 1: Adapte ao seu ambiente de computação

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tradução assistida por IA - [saiba mais e sugira melhorias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Em [Nextflow Run](../nextflow_run/index.md), você configurou as entradas, os parâmetros e as saídas de um pipeline.
Este curso cobre a outra metade do cenário: adaptar a execução de um pipeline a qualquer ambiente de computação em que ele precise rodar, sem alterar o código do fluxo de trabalho.

!!! example "Cenário"

    Você desenvolveu e testou seu pipeline no seu laptop usando Docker.
    Agora você precisa repassá-lo: um colaborador só tem Conda configurado, e o cluster HPC da sua instituição espera que os jobs passem pelo seu próprio agendador com seus próprios limites de recursos.
    Nada disso deve exigir a reescrita do próprio pipeline.

O mesmo código de pipeline pode rodar em todos esses lugares, porque nada disso está embutido no fluxo de trabalho.
Empacotamento de software, plataforma de execução e alocação de recursos são todos controlados por configuração, em camadas sobre o código, e é isso que este curso cobre: como adaptar o mesmo pipeline a um novo ambiente alterando a configuração, não o código.

---

## 1. Selecione uma tecnologia de empacotamento de software

Em [Nextflow Run](../nextflow_run/index.md), você viu um profile `conda` já configurado no `nextflow.config` como alternativa ao Docker.
Aqui você vai construir essa mesma troca você mesmo, e ver o que é necessário para tornar um processo realmente utilizável com Conda.

### 1.1. Desabilite o Docker e habilite o Conda

Mude `docker.enabled` para `false` e adicione uma diretiva habilitando o Conda.

=== "Depois"

    ```groovy title="nextflow.config" linenums="1" hl_lines="1-2"
    docker.enabled = false
    conda.enabled = true
    ```

=== "Antes"

    ```groovy title="nextflow.config" linenums="1"
    docker.enabled = true
    ```

Isso permite que o Nextflow crie e use ambientes Conda para qualquer processo que tenha um pacote Conda especificado.
O processo `cowpy` ainda não tem um, então vamos adicionar um, inteiramente via configuração.

### 1.2. Adicione um pacote Conda via configuração

Uma diretiva `conda` pode ser definida na própria definição do processo, da mesma forma que `container` já está em `modules/cowpy.nf`, mas não precisa ser: `withName` permite defini-la a partir da configuração, com escopo apenas para o processo `cowpy`.

=== "Depois"

    ```groovy title="nextflow.config" linenums="6" hl_lines="3-5"
    process {
        memory = 1.GB
        withName: 'cowpy' {
            conda = 'conda-forge::cowpy==1.1.5'
        }
    }
    ```

=== "Antes"

    ```groovy title="nextflow.config" linenums="6"
    process {
        memory = 1.GB
    }
    ```

Isso não substitui a diretiva `container` já presente no código do pipeline, mas adiciona uma alternativa ao lado dela, sem tocar nesse código.

!!! tip "Dica"

    A busca no [Seqera Containers](https://seqera.io/containers/) é uma forma conveniente de encontrar o URI do pacote Conda para uma determinada ferramenta, mesmo que você não esteja planejando construir um contêiner a partir dele.

### 1.3. Execute o fluxo de trabalho para verificar que ele pode usar Conda

```bash
nextflow run main.nf --batch conda
```

??? success "Saída do comando"

    ```console
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [romantic_pike] revision: c3c85dec78

    executor >  local (8)
    [6d/d48030] sayHello (2)       | 3 of 3 ✔
    [f5/7a9d76] convertToUpper (1) | 3 of 3 ✔
    [1c/79b693] collectGreetings   | 1 of 1 ✔
    Creating env using conda: conda-forge::cowpy==1.1.5 [cache /workspaces/training/execution-config/work/conda/env-898314d566668b6587ad714ae06b8520]
    [bb/64b67c] cowpy              | 1 of 1 ✔

    Outputs:

      /workspaces/training/execution-config/results

      first_output:
        - full_pipeline/intermediates/Hola-output.txt
        - full_pipeline/intermediates/Bonjour-output.txt
        - full_pipeline/intermediates/Hello-output.txt

      uppercased:
        - full_pipeline/intermediates/UPPER-Hello-output.txt
        - full_pipeline/intermediates/UPPER-Bonjour-output.txt
        - full_pipeline/intermediates/UPPER-Hola-output.txt

      collected: full_pipeline/intermediates/COLLECTED-conda-output.txt

      batch_report: full_pipeline/conda-report.txt

      cowpy_art: full_pipeline/cowpy-COLLECTED-conda-output.txt
    ```

Isso produz a mesma saída que rodar com Docker, mesmo que os mecanismos sejam diferentes nos bastidores: o Nextflow recupera o pacote Conda e constrói um ambiente a partir dele, em vez de baixar uma imagem de contêiner.

!!! info "Info"

    Construir um novo ambiente Conda pode demorar um pouco mais do que baixar um contêiner pela primeira vez, mas o pacote usado aqui é pequeno, então deve ser rápido.

Agora volte para Docker para o restante deste curso.

```groovy title="nextflow.config" linenums="1"
docker.enabled = true
```

??? tip "Misturando Docker e Conda"

    Como essas configurações têm escopo por processo, você pode misturá-las: alguns processos usam Docker, outros usam Conda, dependendo do que está disponível para cada ferramenta.
    Se tanto uma diretiva `container` (no código do pipeline) quanto uma diretiva `conda` (aqui, via configuração) estiverem definidas para o mesmo processo e ambos os sistemas de empacotamento estiverem habilitados, o Nextflow prioriza os contêineres.

### Conclusão

Você sabe como configurar qual tecnologia de empacotamento de software um processo deve usar, e como alternar entre Docker e Conda.

### O que vem a seguir?

Aprenda como mudar a plataforma de execução que o Nextflow usa para realmente rodar suas tarefas.

---

## 2. Selecione uma plataforma de execução

Todo pipeline que você executou até agora usou o executor local: cada tarefa roda na mesma máquina que o próprio Nextflow.
O Nextflow verifica os CPUs e a memória disponíveis, e retém as tarefas até que recursos suficientes sejam liberados.

O executor local é conveniente, mas não escala além de uma única máquina.
O Nextflow suporta [muitos outros backends de execução](https://nextflow.io/docs/latest/executor.html), incluindo agendadores HPC (Slurm, LSF, SGE, PBS e outros) e plataformas de nuvem (AWS Batch, Google Cloud Batch, Azure Batch, Kubernetes e mais).

### 2.1. Direcione para um backend diferente

O executor é definido por uma diretiva de processo chamada `executor`.
Por padrão é `local`, então o seguinte está implícito:

```groovy title="Built-in configuration"
process {
    executor = 'local'
}
```

Para direcionar para um backend diferente, defina a diretiva para o executor desejado.

```groovy title="nextflow.config"
process {
    executor = 'slurm'
}
```

!!! warning "Aviso"

    O ambiente de treinamento não está conectado a um cluster HPC, então isso não é algo que você pode executar aqui.

### 2.2. A sintaxe específica de cada backend é abstraída

A maioria das plataformas HPC exige que os envios de jobs especifiquem requisições de recursos, como CPUs, memória e nome da fila, usando sua própria sintaxe.
A mesma requisição de 8 CPUs e 4 GB de RAM em uma fila chamada `my-science-work` parece completamente diferente dependendo do agendador.

??? abstract "Exemplos"

    ```bash title="Config for SLURM / submit using sbatch"
    #SBATCH -o /path/to/my/task/directory/my-task-1.log
    #SBATCH --no-requeue
    #SBATCH -c 8
    #SBATCH --mem 4096M
    #SBATCH -p my-science-work
    ```

    ```bash title="Config for PBS / submit using qsub"
    #PBS -o /path/to/my/task/directory/my-task-1.log
    #PBS -j oe
    #PBS -q my-science-work
    #PBS -l nodes=1:ppn=8
    #PBS -l mem=4gb
    ```

    ```bash title="Config for SGE / submit using qsub"
    #$ -o /path/to/my/task/directory/my-task-1.log
    #$ -j y
    #$ -terse
    #$ -notify
    #$ -q my-science-work
    #$ -l slots=8
    #$ -l h_rss=4096M,mem_free=4096M
    ```

O Nextflow abstrai tudo isso: você especifica propriedades padronizadas como `cpus`, `memory` e `queue` uma vez (veja as [diretivas de processo](https://nextflow.io/docs/latest/reference/process.html#process-directives) para a lista completa), e o Nextflow as traduz para os scripts específicos de cada backend em tempo de execução.

### 2.3. Veja o que o Nextflow realmente executa

Essa tradução não é apenas uma conveniência de arquivo de configuração: ela é respaldada por algo concreto que você pode inspecionar agora mesmo, mesmo com o executor local.
Em [Nextflow Run, seção 1.3](../nextflow_run/01_run_nextflow.md#13-explore-the-work-directory), você olhou dentro de um diretório de tarefa em `work/` e encontrou `.command.sh`, o comando exato que o Nextflow executou.
Esse mesmo diretório também contém um arquivo que você ainda não viu: `.command.run`.

```bash
cat work/0a/0df4a1*/.command.run
```

??? success "Saída do comando (trecho)"

    ```console
    #!/bin/bash
    ### ---
    ### name: 'convertToUpper (3)'
    ### container: 'null'
    ### outputs:
    ### - 'UPPER-Bonjour-output.txt'
    ### ...
    set -e
    set -u
    ...
    nxf_launch() {
        /bin/bash -ue /workspaces/training/nextflow-run/work/0a/0df4a1028c2001758b1841cff92fc7/.command.sh
    }
    ...
    ```

`.command.run` é o script real que o Nextflow entrega para execução.
Ele envolve `.command.sh` com tudo o que é necessário para realmente executá-lo: configuração do ambiente, staging de entrada/saída e reporte do resultado de volta ao Nextflow.
Com o executor `local`, o Nextflow simplesmente executa esse script na mesma máquina.

É exatamente isso que muda quando você define um `executor` diferente.
Para um agendador HPC como Slurm ou PBS, o Nextflow gera esse mesmo tipo de script wrapper, adiciona o cabeçalho específico do agendador que você viu em [2.2](#22-backend-specific-syntax-is-abstracted-away) (traduzido a partir das suas configurações de `cpus`, `memory` e `queue`), e entrega o resultado ao próprio comando de envio desse agendador, por exemplo `sbatch` para Slurm.
A partir daí, o Nextflow consulta o agendador para verificar o status do job em vez de monitorar um processo local diretamente.
Os backends de nuvem funcionam de forma um pouco diferente, pois são acionados por chamadas de API em vez de um comando de envio, mas a mesma ideia subjacente se aplica: o mesmo script de tarefa é executado, apenas como ele é iniciado e rastreado muda.

### Conclusão

Você sabe como mudar o executor para direcionar para diferentes infraestruturas de computação, que o Nextflow abstrai a sintaxe de envio específica de cada backend, e o que realmente acontece nos bastidores quando uma tarefa roda em um backend diferente.

### O que vem a seguir?

Siga para a [Parte 2](./02_resources_and_retries.md), onde você aprenderá como perfilar e alocar recursos de computação, e lidar com falhas de tarefas com novas tentativas.

---

## Resumo

Nesta parte você aprendeu a:

- Alternar a tecnologia de empacotamento de software entre Docker e Conda
- Adicionar uma diretiva `conda` à definição de um processo
- Mudar a plataforma de execução com a diretiva `executor`
- Inspecionar o que o Nextflow realmente gera e executa para uma tarefa, e como isso muda entre os executores
