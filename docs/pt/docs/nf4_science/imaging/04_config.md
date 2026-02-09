# Parte 4: Configuração

Nas Partes 1-3, aprendemos como executar o Nextflow, executar um pipeline nf-core e gerenciar entradas com arquivos de parâmetros e planilhas de amostras.
Agora vamos explorar como configurar pipelines para diferentes ambientes computacionais usando **arquivos de configuração** e **perfis**.

## Objetivos de aprendizagem

Ao final desta parte, você será capaz de:

- Entender como o Nextflow resolve a configuração de múltiplas fontes
- Usar perfis integrados do nf-core para contêineres e testes
- Criar perfis personalizados para diferentes ambientes computacionais
- Personalizar solicitações de recursos usando rótulos de processo
- Gerenciar limites de recursos em ambientes restritos
- Inspecionar a configuração resolvida com `nextflow config`

---

## 1. Entendendo a configuração do Nextflow

### 1.1. O que é um arquivo de configuração?

O Nextflow usa arquivos de configuração para separar a **lógica do fluxo de trabalho** (o que fazer) das **configurações de execução** (como e onde fazer).

Os arquivos de configuração controlam:

- Motores de contêiner (Docker, Singularity, Conda)
- Recursos computacionais (CPUs, memória, tempo)
- Plataformas de execução (local, HPC, nuvem)
- Parâmetros do pipeline

### 1.2. Precedência de configuração

O Nextflow carrega a configuração de múltiplas fontes, com as fontes posteriores sobrescrevendo as anteriores:

1. **Configuração do pipeline**: `nextflow.config` no repositório do pipeline
2. **Configuração do diretório**: `nextflow.config` no seu diretório de trabalho atual
3. **Configuração do usuário**: `~/.nextflow/config`
4. **Linha de comando**: Parâmetros e opções passados diretamente

Essa abordagem em camadas permite que você mantenha padrões no pipeline, sobrescreva com configurações específicas do usuário e faça ajustes rápidos na linha de comando.

### 1.3. Nossa configuração atual

Vamos ver a configuração que estivemos usando:

```groovy title="nextflow.config"
docker.enabled = true
process {
    resourceLimits = [
        cpus: 2,
        memory: '7.GB',
    ]
}

```

Vamos comentar ou alterar de volta a linha `docker.enabled = true` da Parte 2, e descobrir como podemos alcançar o mesmo resultado usando um perfil no molkart.

---

## 2. Usando perfis

### 2.1. O que são perfis?

Perfis são conjuntos nomeados de configuração que podem ser ativados com a flag `-profile` através do comando `nextflow run`.
Eles facilitam a alternância entre diferentes cenários computacionais sem editar arquivos de configuração.

Todos os pipelines nf-core vêm com vários perfis padrão que podemos usar.

### 2.2. Inspecionando perfis integrados

Vamos inspecioná-los no arquivo `molkart/nextflow.config` associado à base de código do pipeline:

```bash
code molkart/nextflow.config
```

Procure pelo bloco `profiles`:

```groovy title="molkart/nextflow.config (excerpt)"
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

Perfis de contêiner comuns:

- `docker`: Usa contêineres Docker (mais comum para desenvolvimento local)
- `singularity`: Usa Singularity/Apptainer (comum em HPC)
- `conda`: Usa ambientes Conda
- `apptainer`: Usa contêineres Apptainer

### 2.3. Re-executando com perfis em vez de nextflow.config

Agora que desabilitamos a configuração do docker em nosso arquivo `nextflow.config` local e entendemos os perfis, vamos re-executar o pipeline usando a flag `-profile`.

Anteriormente na Parte 3, criamos um arquivo `params.yaml` com nossos parâmetros personalizados.
Agora podemos combinar isso com o perfil Docker integrado:

```bash
nextflow run ./molkart \
  -profile docker \
  -params-file params.yaml \
  -resume
```

Vamos detalhar o que cada flag faz:

- `-profile docker`: Ativa o perfil Docker do `nextflow.config` do molkart, que define `docker.enabled = true`
- `-params-file params.yaml`: Carrega todos os parâmetros do pipeline do nosso arquivo YAML
- `-resume`: Reutiliza resultados em cache de execuções anteriores

Como estamos usando `-resume`, o Nextflow verificará se algo mudou desde a última execução.
Se os parâmetros, entradas e código forem os mesmos, todas as tarefas serão recuperadas do cache e o pipeline será concluído quase instantaneamente.

```console title="Output (excerpt)"
executor >  local (12)
...
[1a/2b3c4d] NFCORE_MOLKART:MOLKART:MINDAGAP_MINDAGAP (mem_only)   [100%] 2 of 2, cached: 2 ✔
[5e/6f7g8h] NFCORE_MOLKART:MOLKART:CLAHE (mem_only)               [100%] 2 of 2, cached: 2 ✔
...
-[nf-core/molkart] Pipeline completed successfully-
```

Observe que todos os processos mostram `cached: 2` ou `cached: 1` - nada foi re-executado!

### 2.4. Perfis de teste

Perfis de teste fornecem maneiras rápidas de especificar parâmetros de entrada padrão e arquivos de dados para permitir que você verifique se o pipeline funciona.
Os pipelines nf-core sempre incluirão pelo menos dois perfis de teste:

- `test`: Conjunto de dados pequeno com parâmetros rápidos para testes rápidos
- `test_full`: Teste mais abrangente com dados maiores

Vamos dar uma olhada mais de perto no perfil `test` no molkart, que é incluído usando a diretiva `includeConfig`:

```groovy title="molkart/nextflow.config (excerpt)"
profiles {
  ...
    test      { includeConfig 'conf/test.config'      }
}
```

Isso significa que sempre que executarmos o pipeline com `-profile test`, o Nextflow carregará a configuração de `conf/test.config`.

```groovy title="molkart/conf/test.config (excerpt)"
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

Observe que este perfil contém os mesmos parâmetros que usamos em nosso arquivo `params.yaml` anteriormente.

Você pode ativar múltiplos perfis separando-os com vírgulas.
Vamos usar isso para testar nosso pipeline sem precisar do nosso arquivo de parâmetros:

```bash
nextflow run ./molkart -profile docker,test --outdir results -resume
```

Isso combina:

- `docker`: Habilita contêineres Docker
- `test`: Usa conjunto de dados e parâmetros de teste

Os perfis são aplicados da esquerda para a direita, então perfis posteriores sobrescrevem os anteriores se definirem os mesmos valores.

### Conclusão

Os pipelines nf-core vêm com perfis integrados para contêineres, testes e ambientes especiais.
Você pode combinar múltiplos perfis para construir a configuração que você precisa.

### O que vem a seguir?

Aprenda como criar seus próprios perfis personalizados para diferentes ambientes computacionais.

---

## 3. Criando perfis personalizados

### 3.1. Criar perfis para alternar entre desenvolvimento local e execução em HPC

Vamos criar perfis personalizados para dois cenários:

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

    Não podemos testar o perfil HPC neste ambiente de treinamento, pois não temos acesso a um agendador Slurm.
    Mas isso mostra como você configuraria para uso no mundo real.

### 3.2. Use `nextflow config` para inspecionar a configuração

O comando `nextflow config` mostra a configuração totalmente resolvida sem executar o pipeline.

Visualize a configuração padrão:

```bash
nextflow config ./molkart
```

Visualize a configuração com um perfil específico:

```bash
nextflow config -profile local_dev ./molkart
```

Isso é extremamente útil para:

- Depurar problemas de configuração
- Entender quais valores serão realmente usados
- Verificar como múltiplos perfis interagem

### Conclusão

Perfis personalizados permitem que você alterne entre diferentes ambientes computacionais com uma única flag de linha de comando.
Use `nextflow config` para inspecionar a configuração resolvida antes de executar.

### O que vem a seguir?

Aprenda como personalizar solicitações de recursos para processos individuais usando o sistema de rótulos de processo do nf-core.

---

## 4. Personalizando solicitações de recursos

### 4.1. Entendendo rótulos de processo em pipelines nf-core

Para simplificar, os pipelines nf-core usam [**rótulos de processo**](https://www.nextflow.io/docs/latest/reference/process.html#process-label) para padronizar a alocação de recursos em todos os pipelines.
Cada processo é marcado com um rótulo como `process_low`, `process_medium` ou `process_high` para descrever requisitos de recursos computacionais baixos, médios ou altos, respectivamente.
Esses rótulos são convertidos em solicitações de recursos específicas em um dos arquivos de configuração localizados no diretório `conf/` do pipeline.

```groovy title="molkart/conf/base.config (excerpt)"
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

Observe o multiplicador `task.attempt` - isso permite que tentativas subsequentes de tarefas solicitem mais recursos, se o pipeline estiver configurado com `process.maxRetries > 1`.

### 4.2. Sobrescrevendo recursos para processos específicos

Para controle refinado, direcione processos individuais por nome:

```groovy title="nextflow.config"
process {
    withName: 'NFCORE_MOLKART:MOLKART:CELLPOSE' {
        cpus   = 16
        memory = 32.GB
    }
}
```

Se tentarmos executar este pipeline com a sobrescrita acima, o processo `CELLPOSE` solicitará 16 CPUs e 32 GB de memória em vez do padrão definido por seu rótulo.
Isso fará com que o pipeline falhe em nosso ambiente atual, pois não temos tanta RAM disponível.
Aprenderemos como prevenir esses tipos de falhas na próxima seção.

!!! tip "Dica"

    Para encontrar nomes de processos, veja a saída de execução do pipeline ou verifique `.nextflow.log`.
    Os nomes de processos seguem o padrão `WORKFLOW:SUBWORKFLOW:PROCESS`.

### Conclusão

Os pipelines nf-core usam rótulos de processo para padronizar a alocação de recursos.
Você pode sobrescrever recursos por rótulo (afeta múltiplos processos) ou por nome (afeta um processo específico).

### O que vem a seguir?

Aprenda como gerenciar limites de recursos em ambientes restritos como GitHub Codespaces.

---

## 5. Gerenciando recursos em ambientes restritos

### 5.1. O problema dos limites de recursos

Se tentássemos executar o molkart com um processo solicitando 16 CPUs e 32 GB de memória (como mostrado na seção 4.2), ele falharia em nosso ambiente atual porque não temos tantos recursos disponíveis.
Em um ambiente de cluster com nós maiores, tais solicitações seriam enviadas ao agendador.

Em ambientes restritos como GitHub Codespaces, sem limites, o Nextflow se recusaria a executar processos que excedem os recursos disponíveis.

### 5.2. Definindo limites de recursos

A diretiva `resourceLimits` limita as solicitações de recursos em valores especificados:

```groovy title="nextflow.config"
process {
    resourceLimits = [ cpus: 2, memory: 7.GB ]
}
```

Isso diz ao Nextflow: "Se algum processo solicitar mais de 2 CPUs ou 7 GB de memória, limite-o a esses valores."

### 5.3. Adicionando limites de recursos a perfis personalizados

Atualize seus perfis personalizados para incluir limites apropriados:

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

    Definir limites de recursos muito baixos pode fazer com que os processos falhem ou sejam executados lentamente.
    O pipeline pode precisar usar algoritmos menos intensivos em memória ou processar dados em partes menores.

### Conclusão

Use `resourceLimits` para executar pipelines em ambientes com recursos restritos, limitando as solicitações de recursos dos processos.
Diferentes perfis podem ter diferentes limites apropriados para seu ambiente.

### O que vem a seguir?

Você completou o treinamento principal de Nextflow para Bioimagem!

---

## Conclusão

Agora você entende como configurar pipelines Nextflow para diferentes ambientes computacionais.

Habilidades principais que você aprendeu:

- **Precedência de configuração**: Como o Nextflow resolve configurações de múltiplas fontes
- **Perfis nf-core**: Usando perfis integrados para contêineres, testes e utilitários
- **Perfis personalizados**: Criando seus próprios perfis para diferentes ambientes
- **Rótulos de processo**: Entendendo e sobrescrevendo solicitações de recursos por rótulo
- **Limites de recursos**: Gerenciando ambientes restritos com `resourceLimits`
- **Inspeção de configuração**: Usando `nextflow config` para depurar e verificar configurações

Essas habilidades de configuração são transferíveis para qualquer pipeline Nextflow e ajudarão você a executar fluxos de trabalho de forma eficiente em máquinas locais, clusters HPC e plataformas de nuvem.

### O que vem a seguir?

Parabéns por completar o curso Nextflow para Bioimagem!

Próximos passos:

- Preencha a pesquisa do curso para fornecer feedback
- Confira [Hello Nextflow](../hello_nextflow/index.md) para aprender mais sobre desenvolvimento de fluxos de trabalho
- Explore [Hello nf-core](../hello_nf-core/index.md) para se aprofundar nas ferramentas nf-core
- Navegue por outros cursos nas [coleções de treinamento](../training_collections/index.md)
