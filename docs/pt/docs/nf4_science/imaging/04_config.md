# Parte 4: Configuração

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tradução assistida por IA - [saiba mais e sugira melhorias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Nas Partes 1-3, aprendemos como executar Nextflow, executar um pipeline nf-core e gerenciar entradas com arquivos de parâmetros e samplesheets.
Agora vamos explorar como configurar pipelines para diferentes ambientes de computação usando **arquivos de configuração** e **profiles**.

## Objetivos de aprendizagem

Ao final desta parte, você será capaz de:

- Entender como o Nextflow resolve a configuração de múltiplas fontes
- Usar profiles integrados do nf-core para contêineres e testes
- Criar profiles personalizados para diferentes ambientes de computação
- Personalizar solicitações de recursos usando process labels
- Gerenciar limites de recursos em ambientes restritos
- Inspecionar a configuração resolvida com `nextflow config`

---

## 1. Entendendo a configuração do Nextflow

### 1.1. O que é um arquivo de configuração?

O Nextflow usa arquivos de configuração para separar a **lógica do fluxo de trabalho** (o que fazer) das **configurações de execução** (como e onde fazer).

Arquivos de configuração controlam:

- Engines de contêineres (Docker, Singularity, Conda)
- Recursos computacionais (CPUs, memória, tempo)
- Plataformas de execução (local, HPC, nuvem)
- Parâmetros do pipeline

### 1.2. Precedência de configuração

O Nextflow carrega a configuração de múltiplas fontes, com fontes posteriores sobrescrevendo as anteriores:

1. **Configuração do pipeline**: `nextflow.config` no repositório do pipeline
2. **Configuração do diretório**: `nextflow.config` no seu diretório de trabalho atual
3. **Configuração do usuário**: `~/.nextflow/config`
4. **Linha de comando**: Parâmetros e opções passados diretamente

Esta abordagem em camadas permite manter valores padrão no pipeline, sobrescrever com configurações específicas do usuário e fazer ajustes rápidos na linha de comando.

### 1.3. Nossa configuração atual

Vamos olhar a configuração que estivemos usando:

```groovy title="nextflow.config"
docker.enabled = true
process {
    resourceLimits = [
        cpus: 2,
        memory: '7.GB',
    ]
}

```

Vamos comentar ou alterar de volta a linha `docker.enabled = true` da Parte 2, e descobrir como podemos alcançar o mesmo resultado usando um profile no molkart.

---

## 2. Usando profiles

### 2.1. O que são profiles?

Profiles são conjuntos nomeados de configuração que podem ser ativados com a flag `-profile` via o comando `nextflow run`.
Eles facilitam a alternância entre diferentes cenários de computação sem editar arquivos de configuração.

Todos os pipelines nf-core vêm com diversos profiles padrão que podemos usar.

### 2.2. Inspecionando profiles integrados

Vamos inspecioná-los no arquivo `molkart/nextflow.config` associado à base de código do pipeline:

```bash
code molkart/nextflow.config
```

Procure pelo bloco `profiles`:

```groovy title="molkart/nextflow.config (trecho)"
profiles {
    docker {
        docker.enabled          = true
        singularity.enabled     = false
        conda.enabled           = false
    }
    singularity {
        singularity.enabled     = true
        singularity.autoMounts  = true
        docker.enabled          = false
    }
    conda {
        conda.enabled           = true
        docker.enabled          = false
        conda.channels          = ['conda-forge', 'bioconda']
    }
}
```

Profiles de contêineres comuns:

- `docker`: Usa contêineres Docker (mais comum para desenvolvimento local)
- `singularity`: Usa Singularity/Apptainer (comum em HPC)
- `conda`: Usa ambientes Conda
- `apptainer`: Usa contêineres Apptainer

### 2.3. Re-executando com profiles ao invés de nextflow.config

Agora que desabilitamos a configuração do docker no nosso arquivo `nextflow.config` local e entendemos profiles, vamos re-executar o pipeline usando a flag `-profile`.

Anteriormente na Parte 3, criamos um arquivo `params.yaml` com nossos parâmetros personalizados.
Agora podemos combinar isso com o profile Docker integrado:

```bash
nextflow run ./molkart \
  -profile docker \
  -params-file params.yaml \
  -resume
```

Vamos detalhar o que cada flag faz:

- `-profile docker`: Ativa o profile Docker do `nextflow.config` do molkart, que define `docker.enabled = true`
- `-params-file params.yaml`: Carrega todos os parâmetros do pipeline do nosso arquivo YAML
- `-resume`: Reutiliza resultados em cache de execuções anteriores

Como estamos usando `-resume`, o Nextflow vai verificar se algo mudou desde a última execução.
Se os parâmetros, entradas e código forem os mesmos, todas as tarefas serão recuperadas do cache e o pipeline será concluído quase instantaneamente.

```console title="Saída (trecho)"
executor >  local (12)
...
[1a/2b3c4d] NFCORE_MOLKART:MOLKART:MINDAGAP_MINDAGAP (mem_only)   [100%] 2 of 2, cached: 2 ✔
[5e/6f7g8h] NFCORE_MOLKART:MOLKART:CLAHE (mem_only)               [100%] 2 of 2, cached: 2 ✔
...
-[nf-core/molkart] Pipeline completed successfully-
```

Note que todos os processos mostram `cached: 2` ou `cached: 1` - nada foi re-executado!

### 2.4. Profiles de teste

Profiles de teste fornecem maneiras rápidas de especificar parâmetros de entrada padrão e arquivos de dados para permitir que você verifique se o pipeline funciona.
Pipelines nf-core sempre incluirão pelo menos dois profiles de teste:

- `test`: Dataset pequeno com parâmetros rápidos para testes rápidos
- `test_full`: Teste mais abrangente com dados maiores

Vamos dar uma olhada mais de perto no profile `test` no molkart que é incluído usando a diretiva `includeConfig`:

```groovy title="molkart/nextflow.config (trecho)"
profiles {
  ...
    test      { includeConfig 'conf/test.config'      }
}
```

Isso significa que sempre que executamos o pipeline com `-profile test`, o Nextflow carregará a configuração de `conf/test.config`.

```groovy title="molkart/conf/test.config (trecho)"
params {
    config_profile_name        = 'Test profile'
    config_profile_description = 'Minimal test dataset to check pipeline function'

    input = 'https://raw.githubusercontent.com/nf-core/test-datasets/molkart/test_data/samplesheets/samplesheet_membrane.csv'
    mindagap_tilesize = 90
    mindagap_boxsize = 7
    mindagap_loopnum = 100
    clahe_pyramid_tile = 368
    segmentation_method = "mesmer,cellpose,stardist"
}

process {
    resourceLimits = [
        cpus: 4,
        memory: '15.GB',
        time: '1.h'
    ]
}
```

Note que este profile contém os mesmos parâmetros que usamos no nosso arquivo `params.yaml` anteriormente.

Você pode ativar múltiplos profiles separando-os com vírgulas.
Vamos usar isso para testar nosso pipeline sem precisar do nosso arquivo params:

```bash
nextflow run ./molkart -profile docker,test --outdir results -resume
```

Isso combina:

- `docker`: Habilita contêineres Docker
- `test`: Usa dataset e parâmetros de teste

Profiles são aplicados da esquerda para a direita, então profiles posteriores sobrescrevem os anteriores se definirem os mesmos valores.

### Conclusão

Pipelines nf-core vêm com profiles integrados para contêineres, testes e ambientes especiais.
Você pode combinar múltiplos profiles para construir a configuração que precisa.

### Próximos passos

Aprenda como criar seus próprios profiles personalizados para diferentes ambientes de computação.

---

## 3. Criando profiles personalizados

### 3.1. Criar profiles para alternar entre desenvolvimento local e execução em HPC

Vamos criar profiles personalizados para dois cenários:

1. Desenvolvimento local com Docker
2. HPC universitário com agendador Slurm e Singularity

Adicione o seguinte ao seu `nextflow.config`:

```groovy title="nextflow.config"
profiles {
    local_dev {
        docker.enabled = true
        process.executor = 'local'
    }

    hpc_cluster {
        singularity.enabled = true
        process.executor = 'slurm'
        process.queue = 'standard_queue'
        singularity.cacheDir = '/shared/containers'
    }
}
```

Agora você pode alternar entre ambientes facilmente:

```bash
# Para desenvolvimento local
nextflow run ./molkart -profile local_dev --input data/samplesheet.csv --outdir results

# Para HPC (quando disponível)
nextflow run ./molkart -profile hpc_cluster --input data/samplesheet.csv --outdir results
```

!!! note "Nota"

    Não podemos testar o profile HPC neste ambiente de treinamento, já que não temos acesso a um agendador Slurm.
    Mas isso mostra como você configuraria para uso no mundo real.

### 3.2. Use `nextflow config` para inspecionar a configuração

O comando `nextflow config` mostra a configuração totalmente resolvida sem executar o pipeline.

Visualize a configuração padrão:

```bash
nextflow config ./molkart
```

Visualize a configuração com um profile específico:

```bash
nextflow config -profile local_dev ./molkart
```

Isso é extremamente útil para:

- Depurar problemas de configuração
- Entender quais valores serão realmente usados
- Verificar como múltiplos profiles interagem

### Conclusão

Profiles personalizados permitem que você alterne entre diferentes ambientes de computação com uma única flag de linha de comando.
Use `nextflow config` para inspecionar a configuração resolvida antes de executar.

### Próximos passos

Aprenda como personalizar solicitações de recursos para processos individuais usando o sistema de process label do nf-core.

---

## 4. Personalizando solicitações de recursos

### 4.1. Entendendo process labels em pipelines nf-core

Por simplicidade, pipelines nf-core usam [**process labels**](https://www.nextflow.io/docs/latest/reference/process.html#process-label) para padronizar a alocação de recursos em todos os pipelines.
Cada processo é marcado com um label como `process_low`, `process_medium` ou `process_high` para descrever requisitos de recursos computacionais baixos, médios ou altos, respectivamente.
Esses labels são convertidos em solicitações de recursos específicas em um dos arquivos de configuração localizados no diretório `conf/` do pipeline.

```groovy title="molkart/conf/base.config (trecho)"
process {
    cpus   = { 1      * task.attempt }
    memory = { 6.GB   * task.attempt }
    time   = { 4.h    * task.attempt }

    errorStrategy = { task.exitStatus in ((130..145) + 104) ? 'retry' : 'finish' }
    maxRetries    = 1
    maxErrors     = '-1'

    withLabel:process_single {
        cpus   = { 1                   }
        memory = { 6.GB * task.attempt }
        time   = { 4.h  * task.attempt }
    }
    withLabel:process_low {
        cpus   = { 2     * task.attempt }
        memory = { 12.GB * task.attempt }
        time   = { 4.h   * task.attempt }
    }
    withLabel:process_medium {
        cpus   = { 6     * task.attempt }
        memory = { 36.GB * task.attempt }
        time   = { 8.h   * task.attempt }
    }
    withLabel:process_high {
        cpus   = { 12    * task.attempt }
        memory = { 72.GB * task.attempt }
        time   = { 16.h  * task.attempt }
    }
}
```

Note o multiplicador `task.attempt` - isso permite que tentativas subsequentes de tarefas solicitem mais recursos, se o pipeline estiver definido com `process.maxRetries > 1`.

### 4.2. Sobrescrevendo recursos para processos específicos

Para controle refinado, direcione processos individuais pelo nome:

```groovy title="nextflow.config"
process {
    withName: 'NFCORE_MOLKART:MOLKART:CELLPOSE' {
        cpus   = 16
        memory = 32.GB
    }
}
```

Se tentarmos executar este pipeline com a sobrescrita acima, o processo `CELLPOSE` solicitará 16 CPUs e 32 GB de memória ao invés do padrão definido por seu label.
Isso fará com que o pipeline falhe no nosso ambiente atual, já que não temos essa quantidade de RAM disponível.
Aprenderemos como prevenir esses tipos de falhas na próxima seção.

!!! tip "Dica"

    Para encontrar nomes de processos, olhe a saída de execução do pipeline ou verifique `.nextflow.log`.
    Nomes de processos seguem o padrão `WORKFLOW:SUBWORKFLOW:PROCESS`.

### Conclusão

Pipelines nf-core usam process labels para padronizar a alocação de recursos.
Você pode sobrescrever recursos por label (afeta múltiplos processos) ou por nome (afeta um processo específico).

### Próximos passos

Aprenda como gerenciar limites de recursos em ambientes restritos como GitHub Codespaces.

---

## 5. Gerenciando recursos em ambientes restritos

### 5.1. O problema dos limites de recursos

Se tentássemos executar o molkart com um processo solicitando 16 CPUs e 32 GB de memória (como mostrado na seção 4.2), falharia no nosso ambiente atual porque não temos tantos recursos disponíveis.
Em um ambiente de cluster com nós maiores, tais solicitações seriam submetidas ao agendador.

Em ambientes restritos como GitHub Codespaces, sem limites, o Nextflow se recusaria a executar processos que excedessem os recursos disponíveis.

### 5.2. Definindo limites de recursos

A diretiva `resourceLimits` limita as solicitações de recursos a valores especificados:

```groovy title="nextflow.config"
process {
    resourceLimits = [ cpus: 2, memory: 7.GB ]
}
```

Isso diz ao Nextflow: "Se qualquer processo solicitar mais de 2 CPUs ou 7 GB de memória, limite a esses valores."

### 5.3. Adicionando limites de recursos a profiles personalizados

Atualize seus profiles personalizados para incluir limites apropriados:

```groovy title="nextflow.config"
profiles {
    local_dev {
        docker.enabled = true
        process.executor = 'local'
        process.resourceLimits = [
            cpus: 2,
            memory: 7.GB
        ]
    }

    hpc_cluster {
        singularity.enabled = true
        process.executor = 'slurm'
        process.queue = 'batch'
        process.resourceLimits = [
            cpus: 32,
            memory: 128.GB,
            time: 24.h
        ]
    }
}
```

!!! warning "Aviso"

    Definir limites de recursos muito baixos pode causar falhas nos processos ou execução lenta.
    O pipeline pode precisar usar algoritmos menos intensivos em memória ou processar dados em pedaços menores.

### Conclusão

Use `resourceLimits` para executar pipelines em ambientes com recursos restritos limitando as solicitações de recursos dos processos.
Diferentes profiles podem ter diferentes limites apropriados para seu ambiente.

### Próximos passos

Você completou o treinamento principal de Nextflow para Bioimagem!

---

## Conclusão

Agora você entende como configurar pipelines Nextflow para diferentes ambientes de computação.

Habilidades principais que você aprendeu:

- **Precedência de configuração**: Como o Nextflow resolve configurações de múltiplas fontes
- **Profiles nf-core**: Usando profiles integrados para contêineres, testes e utilitários
- **Profiles personalizados**: Criando seus próprios profiles para diferentes ambientes
- **Process labels**: Entendendo e sobrescrevendo solicitações de recursos por label
- **Limites de recursos**: Gerenciando ambientes restritos com `resourceLimits`
- **Inspeção de configuração**: Usando `nextflow config` para depurar e verificar configurações

Essas habilidades de configuração são transferíveis para qualquer pipeline Nextflow e ajudarão você a executar fluxos de trabalho de forma eficiente em máquinas locais, clusters HPC e plataformas de nuvem.

### Próximos passos

Parabéns por completar o curso Nextflow para Bioimagem!

Próximos passos:

- Preencha a pesquisa do curso para fornecer feedback
- Confira [Hello Nextflow](../hello_nextflow/index.md) para aprender mais sobre desenvolvimento de fluxos de trabalho
- Explore [Hello nf-core](../hello_nf-core/index.md) para mergulhar mais fundo nas ferramentas nf-core
- Navegue por outros cursos nas [coleções de treinamento](../training_collections/index.md)
