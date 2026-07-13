# Parte 2: Executar nf-core/molkart

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tradução assistida por IA - [saiba mais e sugira melhorias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Na Parte 1, executamos um fluxo de trabalho simples Hello World para entender os conceitos básicos da execução do Nextflow.
Agora vamos executar um pipeline de bioimagem do mundo real: **nf-core/molkart**.

Este pipeline processa dados de transcriptômica espacial de Cartografia Molecular da Resolve Bioscience.
No entanto, os padrões do Nextflow que você aprenderá aqui se aplicam a qualquer pipeline nf-core ou fluxo de trabalho de produção.

## 1. Entendendo pipelines nf-core

Antes de executarmos o pipeline, vamos entender o que é o nf-core e por que ele é importante para executar fluxos de trabalho.

### 1.1. O que é nf-core?

[nf-core](https://nf-co.re/) é uma coleção orientada pela comunidade de pipelines Nextflow de alta qualidade.
Todos os pipelines nf-core seguem a mesma estrutura e convenções, o que significa que depois de aprender a executar um, você pode executar qualquer um deles.

Principais características dos pipelines nf-core:

- **Estrutura padronizada**: Todos os pipelines têm nomes de parâmetros e padrões de uso consistentes
- **Dados de teste integrados**: Cada pipeline inclui perfis de teste para validação rápida
- **Documentação abrangente**: Instruções de uso detalhadas e descrições de parâmetros
- **Controle de qualidade**: Relatórios de QC automatizados usando MultiQC
- **Suporte a contêineres**: Contêineres pré-construídos para reprodutibilidade

!!! tip "Quer aprender mais sobre nf-core?"

    Para uma introdução aprofundada ao desenvolvimento de pipelines nf-core, confira o curso de treinamento [Hello nf-core](../../hello_nf-core/index.md).
    Ele abrange como criar e personalizar pipelines nf-core do zero.

### 1.2. O pipeline molkart

![Pipeline nf-core/molkart](img/molkart.png)

O pipeline [nf-core/molkart](https://nf-co.re/molkart) processa dados de imagem de transcriptômica espacial através de vários estágios:

1. **Pré-processamento de imagem**: Preenchimento de padrão de grade e melhoria opcional de contraste
2. **Segmentação celular**: Múltiplas opções de algoritmo (Cellpose, Mesmer, ilastik, Stardist)
3. **Atribuição de pontos**: Atribuir pontos de transcritos a células segmentadas
4. **Controle de qualidade**: Gerar relatórios abrangentes de QC

As principais saídas são:

- Tabelas de contagem célula por transcrito
- Máscaras de segmentação
- Relatório de controle de qualidade MultiQC

---

## 2. Executar molkart com dados de teste

Antes de começarmos, vamos clonar o repositório molkart localmente para que possamos inspecionar seu código:

```bash
cd /workspaces/training/nf4-science/imaging
git clone --branch 1.2.0 --depth 1 https://github.com/nf-core/molkart
```

Isso cria um diretório `molkart/` contendo o código-fonte completo do pipeline.

!!! note "Por que estamos clonando localmente?"

    Normalmente, você executaria pipelines nf-core diretamente do GitHub usando `nextflow run nf-core/molkart -r 1.2.0`.
    O Nextflow baixa automaticamente a versão solicitada do pipeline para você em `$HOME/.nextflow/assets/nf-core/molkart` e o executa de lá.
    No entanto, para este treinamento, estamos clonando o pipeline para um diretório local diferente para que possamos inspecionar o código mais facilmente.

### 2.1. Entendendo requisitos de contêineres

Antes de executar o pipeline completo, vamos aprender por que os contêineres são essenciais para pipelines nf-core.

Vamos fornecer os parâmetros do pipeline usando um arquivo de parâmetros.
Um arquivo de parâmetros é um arquivo YAML que lista cada parâmetro e seu valor, o que mantém os valores tipados (como inteiros) intactos e mantém a linha de comando curta.

Um arquivo `params.yaml` já está disponível no diretório de trabalho:

```yaml title="params.yaml"
input: "data/samplesheet.csv"
outdir: "results"
mindagap_tilesize: 90
mindagap_boxsize: 7
mindagap_loopnum: 100
clahe_pyramid_tile: 368
segmentation_method: "cellpose"
```

Esses parâmetros são:

- `input`: Caminho para a planilha contendo metadados da amostra
- `mindagap_tilesize`, `mindagap_boxsize`, `mindagap_loopnum`: Parâmetros para preenchimento de padrão de grade
- `clahe_pyramid_tile`: Tamanho do kernel para melhoria de contraste
- `segmentation_method`: Qual(is) algoritmo(s) usar para segmentação celular
- `outdir`: Onde salvar os resultados

Vamos tentar executar o pipeline usando esses parâmetros:

```bash
nextflow run ./molkart -params-file params.yaml
```

!!! Warning "Este comando falhará - isso é intencional!"

    Estamos executando deliberadamente sem contêineres para demonstrar por que eles são necessários.

Após alguns momentos, você verá um erro como este:

??? failure "Saída do comando"

    ```console
    ERROR ~ Error executing process > 'NFCORE_MOLKART:MOLKART:MINDAGAP_DUPLICATEFINDER (mem_only)'

    Caused by:
      Process `NFCORE_MOLKART:MOLKART:MINDAGAP_DUPLICATEFINDER (mem_only)` terminated with an error exit status (127)

    Command executed:

      duplicate_finder.py \
          spots.txt \
          90

    Command exit status:
      127

    Command error:
      .command.sh: line 3: duplicate_finder.py: command not found
    ```

**O que está acontecendo aqui?**

O erro `command not found` (status de saída 127) significa que o Nextflow tentou executar `duplicate_finder.py` mas não conseguiu encontrá-lo no seu sistema.
Isso ocorre porque:

1. O pipeline espera que software especializado de bioinformática esteja instalado
2. Essas ferramentas (como `duplicate_finder.py`, `apply_clahe.dask.py`, etc.) não fazem parte de distribuições Linux padrão
3. Sem contêineres, o Nextflow tenta executar comandos diretamente na sua máquina local

**De onde essas ferramentas deveriam vir?**

Vamos inspecionar um dos módulos de processo para ver como ele declara seus requisitos de software.

Abra o módulo de pré-processamento CLAHE:

```bash
code molkart/modules/local/clahe/main.nf
```

Observe a linha 5 - você verá:

```groovy
container 'ghcr.io/schapirolabor/molkart-local:v0.0.4'
```

Esta linha diz ao Nextflow: "Para executar este processo, use a imagem Docker `ghcr.io/schapirolabor/molkart-local:v0.0.4`, que contém todo o software necessário."

Cada processo declara qual imagem de contêiner fornece suas ferramentas necessárias.
No entanto, o Nextflow só usa esses contêineres se você disser para ele fazer isso!

**A solução: Habilitar Docker na configuração**

### 2.2. Configurar Docker e iniciar o pipeline

Para habilitar o Docker, precisamos mudar `docker.enabled` de `false` para `true` no arquivo `nextflow.config`.

Abra o arquivo de configuração:

```bash
code nextflow.config
```

Altere `docker.enabled = false` para `docker.enabled = true`:

```groovy
docker.enabled = true
process {
    resourceLimits = [
        cpus: 2,
        memory: '7.GB',
    ]
}
```

Agora execute o pipeline novamente, desta vez executando todos os três métodos de segmentação para que possamos compará-los mais tarde:

```bash
nextflow run ./molkart -params-file params.yaml --segmentation_method "mesmer,cellpose,stardist"
```

Desta vez, o Nextflow irá:

1. Ler a configuração `docker.enabled = true` do arquivo de configuração
2. Baixar as imagens Docker necessárias (apenas na primeira vez)
3. Executar cada processo dentro de seu contêiner especificado
4. Executar com sucesso porque todas as ferramentas estão disponíveis dentro dos contêineres

!!! Tip "Por que os contêineres são importantes"

    A maioria dos pipelines nf-core **requer** containerização (Docker, Singularity, Podman, etc.) porque:

    - Eles usam software especializado de bioinformática não disponível em ambientes padrão
    - Os contêineres garantem reprodutibilidade - as mesmas versões exatas de software são executadas em todos os lugares
    - Você não precisa instalar manualmente dezenas de ferramentas e suas dependências

    Para mais detalhes sobre contêineres no Nextflow, consulte [Hello Containers](../../hello_nextflow/05_hello_containers.md) do treinamento Hello Nextflow.

### 2.3. Monitorar a execução

Conforme o pipeline é executado, você verá uma saída semelhante a esta:

??? success "Saída do comando"

    ```console
     N E X T F L O W   ~  version 26.04.4

    Launching `./molkart/main.nf` [exotic_banach] revision: e4152308ec

    WARN: Unrecognized config option 'validation.help.enabled'
    WARN: Unrecognized config option 'validation.defaultIgnoreParams'
    WARN: Unrecognized config option 'validation.monochromeLogs'

    ------------------------------------------------------
                                            ,--./,-.
            ___     __   __   __   ___     /,-._.--~'
      |\ | |__  __ /  ` /  \ |__) |__         }  {
      | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                            `._,._,'
      nf-core/molkart 1.2.0
    ------------------------------------------------------
    Segmentation methods and options
      segmentation_method: mesmer,cellpose,stardist

    Image preprocessing
      mindagap_boxsize   : 7
      mindagap_loopnum   : 100
      clahe_kernel       : 25
      mindagap_tilesize  : 90
      clahe_pyramid_tile : 368

    Input/output options
      input              : data/samplesheet.csv
      outdir             : results

    Generic options
      trace_report_suffix: 2026-06-23_16-25-30

    Core Nextflow options
      runName            : exotic_banach
      containerEngine    : docker
      launchDir          : /workspaces/training/nf4-science/imaging
      workDir            : /workspaces/training/nf4-science/imaging/work
      projectDir         : /workspaces/training/nf4-science/imaging/molkart
      userName           : root
      profile            : standard
      configFiles        : /workspaces/training/nf4-science/imaging/molkart/nextflow.config, /workspaces/training/nf4-science/imaging/nextflow.config

    !! Only displaying parameters that differ from the pipeline defaults !!
    ------------------------------------------------------
    * The pipeline
        https://doi.org/10.5281/zenodo.10650748

    * The nf-core framework
        https://doi.org/10.1038/s41587-020-0439-x

    * Software dependencies
        https://github.com/nf-core/molkart/blob/master/CITATIONS.md

    executor >  local (22)
    [b4/e57ff1] NFC…NDAGAP_MINDAGAP (mem_only) | 2 of 2 ✔
    [2e/cf8910] NFC…T:MOLKART:CLAHE (mem_only) | 2 of 2 ✔
    [94/b1ef70] NFC…RT:CREATE_STACK (mem_only) | 1 of 1 ✔
    [58/4b6426] NFC…DUPLICATEFINDER (mem_only) | 1 of 1 ✔
    [c2/ef7002] NFC…DEEPCELL_MESMER (mem_only) | 1 of 1 ✔
    [4f/ddd328] NFC…OLKART:STARDIST (mem_only) | 1 of 1 ✔
    [8c/dd0362] NFC…OLKART:CELLPOSE (mem_only) | 1 of 1 ✔
    [08/2ddcb5] NFC…KART:MASKFILTER (mem_only) | 3 of 3 ✔
    [56/e30de0] NFC…LKART:SPOT2CELL (mem_only) | 3 of 3 ✔
    [b0/5aa635] NFC…:CREATE_ANNDATA (mem_only) | 3 of 3 ✔
    [5e/39d5d0] NFC…LKART:MOLKARTQC (mem_only) | 3 of 3 ✔
    [8e/8bd365] NFCORE_MOLKART:MOLKART:MULTIQC | 1 of 1 ✔
    -[nf-core/molkart] Pipeline completed successfully-
    Completed at: 23-Jun-2026 16:31:40
    Duration    : 6m 9s
    CPU hours   : (a few seconds)
    Succeeded   : 22
    ```

Observe como esta saída é mais detalhada do que nosso exemplo Hello World por causa das convenções nf-core que o pipeline segue:

- O pipeline mostra sua versão e logotipo
- Os parâmetros de configuração são exibidos
- Múltiplos processos são executados em paralelo (indicado por várias linhas de processo)
- Os nomes dos processos incluem o caminho completo do módulo (ex.: `NFCORE_MOLKART:MOLKART:MINDAGAP_MINDAGAP`)

### 2.4. Entendendo a execução de processos

A linha do executor `executor > local (22)` informa:

- **executor**: Qual ambiente de computação está sendo usado (`local` = sua máquina)
- **(22)**: Número total de tarefas iniciadas

Cada linha de processo mostra:

- **Hash** (`[b4/e57ff1]`): Identificador do diretório de trabalho (como antes)
- **Nome do processo**: Caminho completo do módulo e nome do processo
- **Identificador de entrada**: Nome da amostra entre parênteses
- **Progresso**: Contagem de tarefas e status de conclusão (ex.: `1 of 1 ✔`)

### Conclusão

Você sabe como iniciar um pipeline nf-core com dados de teste e interpretar sua saída de execução.

### O que vem a seguir?

Aprenda onde encontrar os resultados e como interpretá-los.

---

## 3. Encontrar e examinar as saídas

Quando o pipeline é concluído com sucesso, você verá uma mensagem de conclusão e resumo de execução.

### 3.1. Localizar o diretório de resultados

Por padrão, os pipelines nf-core gravam saídas em um diretório especificado pelo parâmetro `outdir`, que definimos como `results/`.

Liste o conteúdo:

```bash
tree results/
```

Você deve ver vários subdiretórios:

```console title="results/"
results/
├── anndata/
├── clahe/
├── mindagap/
├── molkartqc/
├── multiqc/
├── pipeline_info/
├── segmentation/
├── spot2cell/
└── stack/
```

Cada subdiretório contém saídas de um estágio específico do pipeline:

- **mindagap/**: Imagens preenchidas com grade da etapa de pré-processamento MindaGap
- **clahe/**: Imagens com contraste aprimorado do pré-processamento CLAHE
- **stack/**: Pilhas de imagens multicanal criadas para segmentação
- **segmentation/**: Resultados de segmentação de diferentes algoritmos (cellpose/, mesmer/, stardist/, filtered_masks/)
- **spot2cell/**: Tabelas de contagem célula por transcrito
- **anndata/**: Objetos AnnData contendo matrizes célula por transcrito e coordenadas espaciais
- **molkartqc/**: Métricas de controle de qualidade para atribuição de pontos
- **multiqc/**: Relatório abrangente de controle de qualidade
- **pipeline_info/**: Relatórios de execução e logs

### 3.2. Examinar o relatório MultiQC

O relatório MultiQC é um arquivo HTML abrangente que agrega métricas de qualidade de todas as etapas do pipeline.

Abra o relatório no navegador de arquivos e clique no botão "Show Preview" para vê-lo renderizado diretamente no VS Code.

O relatório inclui:

- Estatísticas gerais para todas as amostras
- Métricas de pré-processamento
- Métricas de qualidade de segmentação
- Número de células e pontos detectados

!!! Tip "Dica"

    Os relatórios MultiQC são normalmente incluídos em todos os pipelines nf-core.
    Eles sempre fornecem uma visão geral de alto nível da execução do pipeline e qualidade dos dados.

### 3.3. Examinar as tabelas célula por transcrito

A saída científica mais importante é a tabela de contagem célula por transcrito.
Isso informa quantos de cada transcrito foram detectados em cada célula.

Navegue até o diretório spot2cell:

```bash
ls results/spot2cell/
```

Você encontrará arquivos como:

- `cellxgene_mem_only_cellpose.csv`: Tabela célula por transcrito usando segmentação Cellpose
- `cellxgene_mem_only_mesmer.csv`: Tabela célula por transcrito usando segmentação Mesmer
- `cellxgene_mem_only_stardist.csv`: Tabela célula por transcrito usando segmentação Stardist

Executamos apenas 1 amostra neste conjunto de dados de teste, mas em um experimento real teríamos essas tabelas para cada amostra.
Observe como o Nextflow é capaz de processar múltiplos métodos de segmentação em paralelo, tornando fácil comparar resultados.

### 3.4. Visualizar relatórios de execução

O Nextflow gera vários relatórios de execução automaticamente.

Verifique o diretório pipeline_info:

```bash
ls results/pipeline_info/
```

Arquivos principais:

- **execution_report.html**: Visualização de linha do tempo e uso de recursos
- **execution_timeline.html**: Gráfico de Gantt da execução de processos
- **execution_trace.txt**: Métricas detalhadas de execução de tarefas
- **pipeline_dag.html**: Grafo acíclico dirigido mostrando a estrutura do fluxo de trabalho

Abra o relatório de execução para ver o uso de recursos:

```bash
code results/pipeline_info/execution_report.html
```

Isso mostra:

- Quanto tempo cada processo levou
- Uso de CPU e memória
- Quais tarefas foram armazenadas em cache vs. executadas

!!! Tip "Dica"

    Esses relatórios são incrivelmente úteis para otimizar alocação de recursos e solucionar problemas de desempenho.

### Conclusão

Você sabe como localizar saídas do pipeline, examinar relatórios de controle de qualidade e acessar métricas de execução.

### O que vem a seguir?

Aprenda sobre o diretório de trabalho e como o Nextflow gerencia arquivos intermediários.

---

## 4. Explorar o diretório de trabalho

Assim como no nosso exemplo Hello World, todo o trabalho real acontece no diretório `work/`.

### 4.1. Entendendo a estrutura do diretório de trabalho

O diretório de trabalho contém um subdiretório para cada tarefa que foi executada.
Para este pipeline com 22 tarefas, haverá 22 subdiretórios de trabalho.

Liste o diretório de trabalho:

```bash
ls -d work/*/*/ | head -5
```

Isso mostra os primeiros 5 diretórios de tarefas.

### 4.2. Inspecionar um diretório de tarefa

Escolha um dos hashes de processo de segmentação da saída do console (ex.: `[3m/4n5o6p]`) e olhe dentro:

```bash
ls -la work/3m/4n5o6p*/
```

Você verá:

- **Arquivos .command.\***: Scripts e logs de execução do Nextflow (como antes)
- **Arquivos de entrada preparados**: Links simbólicos para os arquivos de entrada reais
- **Arquivos de saída**: Máscaras de segmentação, resultados intermediários, etc.

A diferença principal do Hello World:

- Pipelines reais preparam arquivos de entrada grandes (imagens, dados de referência)
- Arquivos de saída podem ser bastante grandes (máscaras de segmentação, imagens processadas)
- Múltiplos arquivos de entrada e saída por tarefa

!!! Tip "Dica"

    Se um processo falhar, você pode navegar até seu diretório de trabalho, examinar `.command.err` para mensagens de erro e até mesmo executar `.command.sh` manualmente para depurar o problema.

### 4.3. Limpeza do diretório de trabalho

O diretório de trabalho pode se tornar bastante grande ao longo de múltiplas execuções do pipeline.
Como aprendemos na Parte 1, você pode usar `nextflow clean` para remover diretórios de trabalho de execuções antigas.

No entanto, para pipelines nf-core com arquivos intermediários grandes, é especialmente importante limpar regularmente.

### Conclusão

Você entende como os pipelines nf-core organizam seus diretórios de trabalho e como inspecionar tarefas individuais para depuração.

### O que vem a seguir?

Aprenda sobre o cache do Nextflow e como retomar execuções de pipeline que falharam.

---

## 5. Retomar uma execução de pipeline

Um dos recursos mais poderosos do Nextflow é a capacidade de retomar um pipeline do ponto de falha.

### 5.1. O mecanismo de cache

Quando você executa um pipeline com `-resume`, o Nextflow:

1. Verifica o cache para cada tarefa
2. Se as entradas, código e parâmetros são idênticos, reutiliza o resultado em cache
3. Apenas executa novamente as tarefas que mudaram ou falharam

Isso é essencial para pipelines de longa duração onde falhas podem ocorrer tarde na execução.

### 5.2. Tentar resume com molkart

Execute o mesmo comando novamente, mas adicione `-resume`:

```bash
nextflow run ./molkart -params-file params.yaml --segmentation_method "mesmer,cellpose,stardist" -resume
```

Você deve ver uma saída como:

```console
executor >  local (1)
[43/e03702] NFC…NDAGAP_MINDAGAP (mem_only) | 2 of 2, cached: 2 ✔
[2e/cf8910] NFC…T:MOLKART:CLAHE (mem_only) | 2 of 2, cached: 2 ✔
[94/b1ef70] NFC…RT:CREATE_STACK (mem_only) | 1 of 1, cached: 1 ✔
[58/4b6426] NFC…DUPLICATEFINDER (mem_only) | 1 of 1, cached: 1 ✔
[c2/ef7002] NFC…DEEPCELL_MESMER (mem_only) | 1 of 1, cached: 1 ✔
[4f/ddd328] NFC…OLKART:STARDIST (mem_only) | 1 of 1, cached: 1 ✔
[8c/dd0362] NFC…OLKART:CELLPOSE (mem_only) | 1 of 1, cached: 1 ✔
[08/2ddcb5] NFC…KART:MASKFILTER (mem_only) | 3 of 3, cached: 3 ✔
[56/e30de0] NFC…LKART:SPOT2CELL (mem_only) | 3 of 3, cached: 3 ✔
[b0/5aa635] NFC…:CREATE_ANNDATA (mem_only) | 3 of 3, cached: 3 ✔
[5e/39d5d0] NFC…LKART:MOLKARTQC (mem_only) | 3 of 3, cached: 3 ✔
[73/239f45] NFCORE_MOLKART:MOLKART:MULTIQC | 1 of 1 ✔
-[nf-core/molkart] Pipeline completed successfully-
```

Observe a anotação `cached: N` em cada processo de pré-processamento e segmentação — essas tarefas foram reutilizadas em vez de executadas novamente.

### 5.3. Quando resume é útil

Resume é particularmente valioso quando:

- Um pipeline falha devido a limites de recursos (memória insuficiente, limite de tempo excedido)
- Você precisa modificar processos posteriores sem executar novamente etapas anteriores
- Sua conexão de rede cai durante o download de dados
- Você quer adicionar saídas adicionais sem refazer a computação

!!! Warning "Aviso"

    Resume só funciona se você não alterou os dados de entrada, código do pipeline ou parâmetros.
    Se você alterar qualquer um destes, o Nextflow irá corretamente executar novamente as tarefas afetadas.

### Conclusão

Você sabe como usar `-resume` para executar pipelines de forma eficiente sem repetir tarefas bem-sucedidas.

### O que vem a seguir?

Agora que você pode executar nf-core/molkart com dados de teste, você está pronto para aprender como configurá-lo para seus próprios conjuntos de dados.
