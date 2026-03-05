# Parte 1: Executar um pipeline de demonstração

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tradução assistida por IA - [saiba mais e sugira melhorias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Nesta primeira parte do curso de treinamento Hello nf-core, mostramos como encontrar e experimentar um pipeline do nf-core, entender como o código é organizado e reconhecer como ele difere do código Nextflow simples, conforme mostrado em [Hello Nextflow](../hello_nextflow/index.md).

Vamos usar um pipeline chamado nf-core/demo que é mantido pelo projeto nf-core como parte de seu inventário de pipelines para demonstrar estrutura de código e operações de ferramentas.

Certifique-se de que seu diretório de trabalho esteja definido como `hello-nf-core/` conforme instruído na página [Primeiros passos](./00_orientation.md).

---

## 1. Encontrar e recuperar o pipeline nf-core/demo

Vamos começar localizando o pipeline nf-core/demo no site do projeto em [nf-co.re](https://nf-co.re), que centraliza todas as informações, como: documentação geral e artigos de ajuda, documentação para cada um dos pipelines, posts de blog, anúncios de eventos e assim por diante.

### 1.1. Encontrar o pipeline no site

No seu navegador web, vá para [https://nf-co.re/pipelines/](https://nf-co.re/pipelines/) e digite `demo` na barra de pesquisa.

![search results](./img/search-results.png)

Clique no nome do pipeline, `demo`, para acessar a página de documentação do pipeline.

Cada pipeline lançado tem uma página dedicada que inclui as seguintes seções de documentação:

- **Introduction:** Uma introdução e visão geral do pipeline
- **Usage:** Descrições de como executar o pipeline
- **Parameters:** Parâmetros do pipeline agrupados com descrições
- **Output:** Descrições e exemplos dos arquivos de saída esperados
- **Results:** Exemplos de arquivos de saída gerados a partir do conjunto de dados de teste completo
- **Releases & Statistics:** Histórico de versões do pipeline e estatísticas

Sempre que você estiver considerando adotar um novo pipeline, você deve ler a documentação do pipeline cuidadosamente primeiro para entender o que ele faz e como deve ser configurado antes de tentar executá-lo.

Dê uma olhada agora e veja se você consegue descobrir:

- Quais ferramentas o pipeline executará (Verifique a aba: `Introduction`)
- Quais entradas e parâmetros o pipeline aceita ou requer (Verifique a aba: `Parameters`)
- Quais são as saídas produzidas pelo pipeline (Verifique a aba: `Output`)

#### 1.1.1. Visão geral do pipeline

A aba `Introduction` fornece uma visão geral do pipeline, incluindo uma representação visual (chamada de mapa de metrô) e uma lista de ferramentas que são executadas como parte do pipeline.

![pipeline subway map](./img/nf-core-demo-subway-cropped.png)

1. Read QC (FASTQC)
2. Adapter and quality trimming (SEQTK_TRIM)
3. Present QC for raw reads (MULTIQC)

#### 1.1.2. Exemplo de linha de comando

A documentação também fornece um arquivo de entrada de exemplo (discutido mais adiante) e um exemplo de linha de comando.

```bash
nextflow run nf-core/demo \
  -profile <docker/singularity/.../institute> \
  --input samplesheet.csv \
  --outdir <OUTDIR>
```

Você notará que o comando de exemplo NÃO especifica um arquivo de fluxo de trabalho, apenas a referência ao repositório do pipeline, `nf-core/demo`.

Quando invocado dessa forma, o Nextflow assumirá que o código está organizado de uma certa maneira.
Vamos recuperar o código para que possamos examinar essa estrutura.

### 1.2. Recuperar o código do pipeline

Depois de determinarmos que o pipeline parece ser adequado para nossos propósitos, vamos experimentá-lo.
Felizmente, o Nextflow facilita a recuperação de pipelines de repositórios formatados corretamente sem precisar baixar nada manualmente.

Vamos retornar ao terminal e executar o seguinte:

```bash
nextflow pull nf-core/demo
```

??? success "Saída do comando"

    ```console
    Checking nf-core/demo ...
    downloaded from https://github.com/nf-core/demo.git - revision: 04060b4644 [master]
    ```

O Nextflow faz um `pull` do código do pipeline, o que significa que ele baixa o repositório completo para sua unidade local.

Para ser claro, você pode fazer isso com qualquer pipeline Nextflow que esteja configurado adequadamente no GitHub, não apenas pipelines do nf-core.
No entanto, o nf-core é a maior coleção de código aberto de pipelines Nextflow.

Você pode fazer o Nextflow fornecer uma lista de quais pipelines você recuperou dessa maneira:

```bash
nextflow list
```

??? success "Saída do comando"

    ```console
    nf-core/demo
    ```

Você notará que os arquivos não estão no seu diretório de trabalho atual.
Por padrão, o Nextflow os salva em `$NXF_HOME/assets`.

```bash
tree -L 2 $NXF_HOME/assets/
```

```console title="Conteúdo do diretório"
/workspaces/.nextflow/assets/
└── nf-core
    └── demo

2 directories, 0 files
```

!!! note "Nota"

    O caminho completo pode ser diferente no seu sistema se você não estiver usando nosso ambiente de treinamento.

O Nextflow mantém o código-fonte baixado intencionalmente 'fora do caminho' com base no princípio de que esses pipelines devem ser usados mais como bibliotecas do que código com o qual você interagiria diretamente.

No entanto, para os propósitos deste treinamento, queremos poder explorar e ver o que há lá dentro.
Então, para facilitar isso, vamos criar um link simbólico para esse local a partir do nosso diretório de trabalho atual.

```bash
ln -s $NXF_HOME/assets pipelines
```

Isso cria um atalho que facilita a exploração do código que acabamos de baixar.

```bash
tree -L 2 pipelines
```

```console title="Conteúdo do diretório"
pipelines
└── nf-core
    └── demo

2 directories, 0 files
```

Agora podemos espiar o código-fonte mais facilmente conforme necessário.

Mas primeiro, vamos tentar executar nosso primeiro pipeline do nf-core!

### Conclusão

Agora você sabe como encontrar um pipeline através do site do nf-core e recuperar uma cópia local do código-fonte.

### Qual é o próximo passo?

Aprenda como experimentar um pipeline do nf-core com o mínimo de esforço.

---

## 2. Experimentar o pipeline com seu perfil de teste

Convenientemente, todo pipeline do nf-core vem com um perfil de teste.
Este é um conjunto mínimo de configurações para o pipeline executar usando um pequeno conjunto de dados de teste hospedado no repositório [nf-core/test-datasets](https://github.com/nf-core/test-datasets).
É uma ótima maneira de experimentar rapidamente um pipeline em pequena escala.

!!! note "Nota"

    O sistema de perfil de configuração do Nextflow permite que você alterne facilmente entre diferentes motores de contêiner ou ambientes de execução.
    Para mais detalhes, consulte [Hello Nextflow Parte 6: Configuração](../hello_nextflow/06_hello_config.md).

### 2.1. Examinar o perfil de teste

É uma boa prática verificar o que o perfil de teste de um pipeline especifica antes de executá-lo.
O perfil `test` para `nf-core/demo` está no arquivo de configuração `conf/test.config` e é mostrado abaixo.

```groovy title="conf/test.config" linenums="1" hl_lines="8 26"
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Nextflow config file for running minimal tests
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Defines input files and everything required to run a fast and simple pipeline test.

    Use as follows:
        nextflow run nf-core/demo -profile test,<docker/singularity> --outdir <OUTDIR>

----------------------------------------------------------------------------------------
*/

process {
    resourceLimits = [
        cpus: 4,
        memory: '4.GB',
        time: '1.h'
    ]
}

params {
    config_profile_name        = 'Test profile'
    config_profile_description = 'Minimal test dataset to check pipeline function'

    // Dados de entrada
    input  = 'https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/samplesheet/samplesheet_test_illumina_amplicon.csv'

}
```

Você notará imediatamente que o bloco de comentário no topo inclui um exemplo de uso mostrando como executar o pipeline com este perfil de teste.

```groovy title="conf/test.config" linenums="7"
Use as follows:
        nextflow run nf-core/demo -profile test,<docker/singularity> --outdir <OUTDIR>
```

As únicas coisas que precisamos fornecer são o que é mostrado entre colchetes angulares no comando de exemplo: `<docker/singularity>` e `<OUTDIR>`.

Como lembrete, `<docker/singularity>` refere-se à escolha do sistema de contêiner. Todos os pipelines do nf-core são projetados para serem usáveis com contêineres (Docker, Singularity, etc.) para garantir reprodutibilidade e eliminar problemas de instalação de software.
Então precisaremos especificar se queremos usar Docker ou Singularity para testar o pipeline.

A parte `--outdir <OUTDIR>` refere-se ao diretório onde o Nextflow escreverá as saídas do pipeline.
Precisamos fornecer um nome para ele, que podemos simplesmente inventar.
Se ainda não existir, o Nextflow o criará para nós em tempo de execução.

Seguindo para a seção após o bloco de comentário, o perfil de teste nos mostra o que foi pré-configurado para teste: mais notavelmente, o parâmetro `input` já está configurado para apontar para um conjunto de dados de teste, então não precisamos fornecer nossos próprios dados.
Se você seguir o link para a entrada pré-configurada, verá que é um arquivo csv contendo identificadores de amostra e caminhos de arquivo para várias amostras experimentais.

```csv title="samplesheet_test_illumina_amplicon.csv"
sample,fastq_1,fastq_2
SAMPLE1_PE,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample1_R1.fastq.gz,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample1_R2.fastq.gz
SAMPLE2_PE,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample2_R1.fastq.gz,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample2_R2.fastq.gz
SAMPLE3_SE,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample1_R1.fastq.gz,
SAMPLE3_SE,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample2_R1.fastq.gz,
```

Isso é chamado de planilha de amostras e é a forma mais comum de entrada para pipelines do nf-core.

!!! note "Nota"

    Não se preocupe se você não estiver familiarizado com os formatos e tipos de dados, isso não é importante para o que se segue.

Então isso confirma que temos tudo o que precisamos para experimentar o pipeline.

### 2.2. Executar o pipeline

Vamos decidir usar Docker para o sistema de contêiner e `demo-results` como o diretório de saída, e estamos prontos para executar o comando de teste:

```bash
nextflow run nf-core/demo -profile docker,test --outdir demo-results
```

??? success "Saída do comando"

    ```console
     N E X T F L O W   ~  version 25.04.3

    Launching `https://github.com/nf-core/demo` [magical_pauling] DSL2 - revision: db7f526ce1 [master]


    ------------------------------------------------------
                                            ,--./,-.
            ___     __   __   __   ___     /,-._.--~'
      |\ | |__  __ /  ` /  \ |__) |__         }  {
      | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                            `._,._,'
      nf-core/demo 1.0.2
    ------------------------------------------------------
    Input/output options
      input                     : https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/samplesheet/samplesheet_test_illumina_amplicon.csv
      outdir                    : demo-results

    Institutional config options
      config_profile_name       : Test profile
      config_profile_description: Minimal test dataset to check pipeline function

    Generic options
      trace_report_suffix       : 2025-11-21_04-57-41

    Core Nextflow options
      revision                  : master
      runName                   : magical_pauling
      containerEngine           : docker
      launchDir                 : /workspaces/training/hello-nf-core
      workDir                   : /workspaces/training/hello-nf-core/work
      projectDir                : /workspaces/.nextflow/assets/nf-core/demo
      userName                  : root
      profile                   : docker,test
      configFiles               : /workspaces/.nextflow/assets/nf-core/demo/nextflow.config

    !! Only displaying parameters that differ from the pipeline defaults !!
    ------------------------------------------------------
    * The pipeline
        https://doi.org/10.5281/zenodo.12192442

    * The nf-core framework
        https://doi.org/10.1038/s41587-020-0439-x

    * Software dependencies
        https://github.com/nf-core/demo/blob/master/CITATIONS.md


    executor >  local (7)
    [ff/a6976b] NFCORE_DEMO:DEMO:FASTQC (SAMPLE3_SE)     | 3 of 3 ✔
    [39/731ab7] NFCORE_DEMO:DEMO:SEQTK_TRIM (SAMPLE3_SE) | 3 of 3 ✔
    [7c/78d96e] NFCORE_DEMO:DEMO:MULTIQC                 | 1 of 1 ✔
    -[nf-core/demo] Pipeline completed successfully-
    ```

Se sua saída corresponder a isso, parabéns! Você acabou de executar seu primeiro pipeline do nf-core.

Você notará que há muito mais saída no console do que quando você executa um pipeline Nextflow básico.
Há um cabeçalho que inclui um resumo da versão do pipeline, entradas e saídas, e alguns elementos de configuração.

!!! note "Nota"

    Sua saída mostrará carimbos de data/hora, nomes de execução e caminhos de arquivo diferentes, mas a estrutura geral e a execução do processo devem ser semelhantes.

Seguindo para a saída de execução, vamos dar uma olhada nas linhas que nos dizem quais processos foram executados:

```console
    [ff/a6976b] NFCORE_DEMO:DEMO:FASTQC (SAMPLE3_SE)     | 3 of 3 ✔
    [39/731ab7] NFCORE_DEMO:DEMO:SEQTK_TRIM (SAMPLE3_SE) | 3 of 3 ✔
    [7c/78d96e] NFCORE_DEMO:DEMO:MULTIQC                 | 1 of 1 ✔
```

Isso nos diz que três processos foram executados, correspondendo às três ferramentas mostradas na página de documentação do pipeline no site do nf-core: FASTQC, SEQTK_TRIM e MULTIQC.

Os nomes completos dos processos como mostrado aqui, como `NFCORE_DEMO:DEMO:MULTIQC`, são mais longos do que o que você pode ter visto no material introdutório do Hello Nextflow.
Estes incluem os nomes de seus fluxos de trabalho pai e refletem a modularidade do código do pipeline.
Vamos entrar em mais detalhes sobre isso daqui a pouco.

### 2.3. Examinar as saídas do pipeline

Finalmente, vamos dar uma olhada no diretório `demo-results` produzido pelo pipeline.

```bash
tree -L 2 demo-results
```

??? abstract "Conteúdo do diretório"

    ```console
    demo-results
    ├── fastqc
    │   ├── SAMPLE1_PE
    │   ├── SAMPLE2_PE
    │   └── SAMPLE3_SE
    ├── fq
    │   ├── SAMPLE1_PE
    │   ├── SAMPLE2_PE
    │   └── SAMPLE3_SE
    ├── multiqc
    │   ├── multiqc_data
    │   ├── multiqc_plots
    │   └── multiqc_report.html
    └── pipeline_info
        ├── execution_report_2025-11-21_04-57-41.html
        ├── execution_timeline_2025-11-21_04-57-41.html
        ├── execution_trace_2025-11-21_04-57-41.txt
        ├── nf_core_demo_software_mqc_versions.yml
        ├── params_2025-11-21_04-57-46.json
        └── pipeline_dag_2025-11-21_04-57-41.html
    ```

Isso pode parecer muito.
Para saber mais sobre as saídas do pipeline `nf-core/demo`, consulte sua [página de documentação](https://nf-co.re/demo/1.0.2/docs/output/).

Nesta etapa, o que é importante observar é que os resultados são organizados por módulo, e há adicionalmente um diretório chamado `pipeline_info` contendo vários relatórios com carimbos de data/hora sobre a execução do pipeline.

Por exemplo, o arquivo `execution_timeline_*` mostra quais processos foram executados, em que ordem e quanto tempo levaram para executar:

![execution timeline report](./img/execution_timeline.png)

!!! note "Nota"

    Aqui as tarefas não foram executadas em paralelo porque estamos executando em uma máquina minimalista no Github Codespaces.
    Para ver essas execuções em paralelo, tente aumentar a alocação de CPU do seu codespace e os limites de recursos na configuração de teste.

Esses relatórios são gerados automaticamente para todos os pipelines do nf-core.

### Conclusão

Você sabe como executar um pipeline do nf-core usando seu perfil de teste integrado e onde encontrar suas saídas.

### Qual é o próximo passo?

Aprenda como o código do pipeline é organizado.

---

## 3. Examinar a estrutura do código do pipeline

Agora que executamos o pipeline com sucesso como usuários, vamos mudar nossa perspectiva para ver como os pipelines do nf-core são estruturados internamente.

O projeto nf-core impõe diretrizes fortes sobre como os pipelines são estruturados e como o código é organizado, configurado e documentado.
Entender como tudo isso é organizado é o primeiro passo para desenvolver seus próprios pipelines compatíveis com o nf-core, que abordaremos na Parte 2 deste curso.

Vamos dar uma olhada em como o código do pipeline está organizado no repositório `nf-core/demo`, usando o link simbólico `pipelines` que criamos anteriormente.

Você pode usar `tree` ou usar o explorador de arquivos para encontrar e abrir o diretório `nf-core/demo`.

```bash
tree -L 1 pipelines/nf-core/demo
```

??? abstract "Conteúdo do diretório"

    ```console
    pipelines/nf-core/demo
    ├── assets
    ├── CHANGELOG.md
    ├── CITATIONS.md
    ├── CODE_OF_CONDUCT.md
    ├── conf
    ├── docs
    ├── LICENSE
    ├── main.nf
    ├── modules
    ├── modules.json
    ├── nextflow.config
    ├── nextflow_schema.json
    ├── nf-test.config
    ├── README.md
    ├── ro-crate-metadata.json
    ├── subworkflows
    ├── tests
    ├── tower.yml
    └── workflows
    ```

Há muita coisa acontecendo lá, então vamos abordar isso passo a passo.

Primeiro, vamos observar que no nível superior, você pode encontrar um arquivo README com informações resumidas, bem como arquivos acessórios que resumem informações do projeto, como licenciamento, diretrizes de contribuição, citação e código de conduta.
A documentação detalhada do pipeline está localizada no diretório `docs`.
Todo esse conteúdo é usado para gerar as páginas da web no site do nf-core programaticamente, então elas estão sempre atualizadas com o código.

Agora, para o resto, vamos dividir nossa exploração em três etapas:

1. Componentes de código do pipeline (`main.nf`, `workflows`, `subworkflows`, `modules`)
2. Configuração do pipeline
3. Entradas e validação

Vamos começar com os componentes de código do pipeline.
Vamos nos concentrar na hierarquia de arquivos e na organização estrutural, em vez de mergulhar no código dentro de arquivos individuais.

### 3.1. Componentes de código do pipeline

A organização padrão de código de pipeline do nf-core segue uma estrutura modular que é projetada para maximizar a reutilização de código, conforme introduzido em [Hello Modules](../hello_nextflow/04_hello_modules.md), Parte 4 do curso [Hello Nextflow](../hello_nextflow/index.md), embora no verdadeiro estilo nf-core, isso seja implementado com um pouco de complexidade adicional.
Especificamente, os pipelines do nf-core fazem uso abundante de subworkflows, ou seja, scripts de fluxo de trabalho que são importados por um fluxo de trabalho pai.

Isso pode parecer um pouco abstrato, então vamos dar uma olhada em como isso é usado na prática no pipeline `nf-core/demo`.

!!! note "Nota"

    Não vamos passar pelo código real de _como_ esses componentes modulares são conectados, porque há alguma complexidade adicional associada ao uso de subworkflows que pode ser confusa, e entender isso não é necessário nesta etapa do treinamento.
    Por enquanto, vamos nos concentrar na organização geral e na lógica.

#### 3.1.1. Visão geral geral

Aqui está como são as relações entre os componentes de código relevantes para o pipeline `nf-core/demo`:

<figure class="excalidraw">
    --8<-- "docs/en/docs/hello_nf-core/img/nf-core_demo_code_organization.svg"
</figure>

Há um script chamado _ponto de entrada_ chamado `main.nf`, que atua como um wrapper para dois tipos de fluxos de trabalho aninhados: o fluxo de trabalho contendo a lógica de análise real, localizado em `workflows/` e chamado `demo.nf`, e um conjunto de fluxos de trabalho de manutenção localizados em `subworkflows/`.
O fluxo de trabalho `demo.nf` chama **módulos** localizados em `modules/`; estes contêm os **processos** que realizarão as etapas de análise reais.

!!! note "Nota"

    Subworkflows não estão limitados a funções de manutenção, e eles podem fazer uso de módulos de processo.

    O pipeline `nf-core/demo` mostrado aqui acontece de estar no lado mais simples do espectro, mas outros pipelines do nf-core (como `nf-core/rnaseq`) utilizam subworkflows que estão envolvidos na análise real.

Agora, vamos revisar esses componentes por vez.

#### 3.1.2. O script de ponto de entrada: `main.nf`

O script `main.nf` é o ponto de entrada de onde o Nextflow começa quando executamos `nextflow run nf-core/demo`.
Isso significa que quando você executa `nextflow run nf-core/demo` para executar o pipeline, o Nextflow automaticamente encontra e executa o script `main.nf`.
Isso funciona para qualquer pipeline Nextflow que segue essa nomeação e estrutura convencional, não apenas pipelines do nf-core.

Usar um script de ponto de entrada facilita a execução de subworkflows de 'manutenção' padronizados antes e depois da execução do script de análise real.
Vamos passar por eles depois de revisarmos o fluxo de trabalho de análise real e seus módulos.

#### 3.1.3. O script de análise: `workflows/demo.nf`

O fluxo de trabalho `workflows/demo.nf` é onde a lógica central do pipeline é armazenada.
Ele é estruturado muito como um fluxo de trabalho Nextflow normal, exceto que é projetado para ser chamado de um fluxo de trabalho pai, o que requer alguns recursos extras.
Vamos cobrir as diferenças relevantes na próxima parte deste curso, quando abordaremos a conversão do pipeline Hello simples do Hello Nextflow em uma forma compatível com o nf-core.

O fluxo de trabalho `demo.nf` chama **módulos** localizados em `modules/`, que revisaremos a seguir.

!!! note "Nota"

    Alguns fluxos de trabalho de análise do nf-core exibem níveis adicionais de aninhamento ao chamar subworkflows de nível inferior.
    Isso é usado principalmente para envolver dois ou mais módulos que são comumente usados juntos em segmentos de pipeline facilmente reutilizáveis.
    Você pode ver alguns exemplos navegando pelos [subworkflows do nf-core](https://nf-co.re/subworkflows/) disponíveis no site do nf-core.

    Quando o script de análise usa subworkflows, eles são armazenados no diretório `subworkflows/`.

#### 3.1.4. Os módulos

Os módulos são onde o código do processo reside, conforme descrito na [Parte 4 do curso de treinamento Hello Nextflow](../hello_nextflow/04_hello_modules.md).

No projeto nf-core, os módulos são organizados usando uma estrutura aninhada de vários níveis que reflete tanto sua origem quanto seu conteúdo.
No nível superior, os módulos são diferenciados como `nf-core` ou `local` (não parte do projeto nf-core), e depois colocados em um diretório nomeado com base na(s) ferramenta(s) que eles envolvem.
Se a ferramenta pertence a um kit de ferramentas (ou seja, um pacote contendo várias ferramentas), então há um nível de diretório intermediário nomeado com base no kit de ferramentas.

Você pode ver isso aplicado na prática aos módulos do pipeline `nf-core/demo`:

```bash
tree -L 3 pipelines/nf-core/demo/modules
```

??? abstract "Conteúdo do diretório"

    ```console
    pipelines/nf-core/demo/modules
    └── nf-core
        ├── fastqc
        │   ├── environment.yml
        │   ├── main.nf
        │   ├── meta.yml
        │   └── tests
        ├── multiqc
        │   ├── environment.yml
        │   ├── main.nf
        │   ├── meta.yml
        │   └── tests
        └── seqtk
            └── trim

    7 directories, 6 files
    ```

Aqui você vê que os módulos `fastqc` e `multiqc` estão no nível superior dentro dos módulos `nf-core`, enquanto o módulo `trim` está sob o kit de ferramentas ao qual pertence, `seqtk`.
Neste caso, não há módulos `local`.

O arquivo de código do módulo que descreve o processo sempre se chama `main.nf`, e é acompanhado por testes e arquivos `.yml` que vamos ignorar por enquanto.

Considerados em conjunto, o fluxo de trabalho de ponto de entrada, fluxo de trabalho de análise e módulos são suficientes para executar as partes 'interessantes' do pipeline.
No entanto, sabemos que também há subworkflows de manutenção lá, então vamos olhar para eles agora.

#### 3.1.5. Os subworkflows de manutenção

Como módulos, subworkflows são diferenciados em diretórios `local` e `nf-core`, e cada subworkflow tem sua própria estrutura de diretório aninhada com seu próprio script `main.nf`, testes e arquivo `.yml`.

```bash
tree -L 3 pipelines/nf-core/demo/subworkflows
```

??? abstract "Conteúdo do diretório"

    ```console
    pipelines/nf-core/demo/subworkflows
    ├── local
    │   └── utils_nfcore_demo_pipeline
    │       └── main.nf
    └── nf-core
        ├── utils_nextflow_pipeline
        │   ├── main.nf
        │   ├── meta.yml
        │   └── tests
        ├── utils_nfcore_pipeline
        │   ├── main.nf
        │   ├── meta.yml
        │   └── tests
        └── utils_nfschema_plugin
            ├── main.nf
            ├── meta.yml
            └── tests

    9 directories, 7 files
    ```

Como observado acima, o pipeline `nf-core/demo` não inclui nenhum subworkflow específico de análise, então todos os subworkflows que vemos aqui são chamados fluxos de trabalho de 'manutenção' ou 'utilitário', como denotado pelo prefixo `utils_` em seus nomes.
Esses subworkflows são o que produz o cabeçalho sofisticado do nf-core na saída do console, entre outras funções acessórias.

!!! tip "Dica"

    Além de seu padrão de nomenclatura, outra indicação de que esses subworkflows não executam nenhuma função verdadeiramente relacionada à análise é que eles não chamam nenhum processo.

Isso completa o resumo dos componentes de código principais que constituem o pipeline `nf-core/demo`.
Agora vamos dar uma olhada nos elementos restantes que você deve saber um pouco antes de mergulhar no desenvolvimento: configuração do pipeline e validação de entrada.

### 3.2. Configuração do pipeline

Você aprendeu anteriormente que o Nextflow oferece muitas opções para configurar a execução do pipeline, seja em termos de entradas e parâmetros, recursos de computação e outros aspectos de orquestração.
O projeto nf-core aplica diretrizes altamente padronizadas para configuração de pipeline que visam construir sobre as opções de personalização flexíveis do Nextflow de uma maneira que forneça maior consistência e manutenibilidade entre os pipelines.

O arquivo de configuração central `nextflow.config` é usado para definir valores padrão para parâmetros e outras opções de configuração.
A maioria dessas opções de configuração são aplicadas por padrão, enquanto outras (por exemplo, perfis de dependência de software) são incluídas como perfis opcionais.

Existem vários arquivos de configuração adicionais que são armazenados na pasta `conf` e que podem ser adicionados à configuração por padrão ou opcionalmente como perfis:

- `base.config`: Um arquivo de configuração 'em branco', apropriado para uso geral na maioria dos ambientes de computação de alto desempenho. Isso define bins amplos de uso de recursos, por exemplo, que são convenientes para aplicar aos módulos.
- `modules.config`: Diretivas de módulo adicionais e argumentos.
- `test.config`: Um perfil para executar o pipeline com dados de teste mínimos, que usamos quando executamos o pipeline de demonstração.
- `test_full.config`: Um perfil para executar o pipeline com um conjunto de dados de teste de tamanho completo.

Vamos tocar em alguns desses arquivos mais tarde no curso.

### 3.3. Entradas e validação

Como observamos anteriormente, quando examinamos o perfil de teste do pipeline `nf-core/demo`, ele foi projetado para receber como entrada uma planilha de amostras contendo caminhos de arquivo e identificadores de amostra.
Os caminhos de arquivo vinculados a dados reais localizados no repositório `nf-core/test-datasets`.

Um exemplo de planilha de amostras também é fornecido no diretório `assets`, embora os caminhos neste não sejam reais.

```csv title="assets/samplesheet.csv" linenums="1"
sample,fastq_1,fastq_2
SAMPLE_PAIRED_END,/path/to/fastq/files/AEG588A1_S1_L002_R1_001.fastq.gz,/path/to/fastq/files/AEG588A1_S1_L002_R2_001.fastq.gz
SAMPLE_SINGLE_END,/path/to/fastq/files/AEG588A4_S4_L003_R1_001.fastq.gz,

```

Esta planilha de amostras específica é bastante simples, mas alguns pipelines são executados em planilhas de amostras que são mais complexas, com muito mais metadados associados às entradas primárias.

Infelizmente, como esses arquivos podem ser difíceis de verificar visualmente, a formatação inadequada de dados de entrada é uma fonte muito comum de falhas de pipeline.
Um problema relacionado é quando os parâmetros são fornecidos incorretamente.

A solução para esses problemas é executar verificações de validação automatizadas em todos os arquivos de entrada para garantir que eles contenham os tipos esperados de informação, formatados corretamente, e em parâmetros para garantir que sejam do tipo esperado.
Isso é chamado de validação de entrada e, idealmente, deve ser feito _antes_ de tentar executar um pipeline, em vez de esperar que o pipeline falhe para descobrir que havia um problema com as entradas.

Assim como para configuração, o projeto nf-core é muito opinativo sobre validação de entrada e recomenda o uso do [plugin nf-schema](https://nextflow-io.github.io/nf-schema/latest/), um plugin Nextflow que fornece recursos abrangentes de validação para pipelines Nextflow.

Vamos cobrir este tópico com mais detalhes na Parte 5 deste curso.
Por enquanto, apenas esteja ciente de que existem dois arquivos JSON fornecidos para esse propósito, `nextflow_schema.json` e `assets/schema_input.json`.

O `nextflow_schema.json` é um arquivo usado para armazenar informações sobre os parâmetros do pipeline, incluindo tipo, descrição e texto de ajuda em um formato legível por máquina.
Isso é usado para vários propósitos, incluindo validação automatizada de parâmetros, geração de texto de ajuda e renderização interativa de formulário de parâmetros em interfaces de UI.

O `schema_input.json` é um arquivo usado para definir a estrutura da planilha de amostras de entrada.
Cada coluna pode ter um tipo, padrão, descrição e texto de ajuda em um formato legível por máquina.
O schema é usado para vários propósitos, incluindo validação automatizada e fornecimento de mensagens de erro úteis.

### Conclusão

Você sabe quais são os principais componentes de um pipeline do nf-core e como o código é organizado; onde os elementos principais de configuração estão localizados; e está ciente de para que serve a validação de entrada.

### Qual é o próximo passo?

Faça uma pausa! Foi muita coisa. Quando você estiver se sentindo renovado e pronto, passe para a próxima seção para aplicar o que você aprendeu para escrever um pipeline compatível com o nf-core.

!!! tip "Dica"

    Se você gostaria de aprender como compor fluxos de trabalho com subworkflows antes de passar para a próxima parte, confira a [Missão Secundária Workflows de Workflows](../side_quests/workflows_of_workflows.md).
