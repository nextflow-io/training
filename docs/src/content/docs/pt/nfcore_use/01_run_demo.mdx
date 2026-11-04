# Parte 1: Executar um pipeline de demonstração

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tradução assistida por IA - [saiba mais e sugira melhorias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Nesta primeira parte do curso Use nf-core, mostramos como encontrar um pipeline nf-core e testá-lo usando seu perfil de teste integrado.

Vamos usar um pipeline chamado nf-core/demo, mantido pelo projeto nf-core como parte de seu inventário de pipelines para fins de demonstração e treinamento.

Certifique-se de que seu diretório de trabalho está definido como `nfcore-use/`, conforme instruído na página [Primeiros passos](./00_orientation.md).

---

## 1. Encontrar e recuperar o pipeline nf-core/demo

Vamos começar localizando o pipeline nf-core/demo no site do projeto em [nf-co.re](https://nf-co.re), que centraliza todas as informações, como: documentação geral e artigos de ajuda, documentação de cada pipeline, posts de blog, anúncios de eventos e muito mais.

### 1.1. Encontrar o pipeline no site

No seu navegador, acesse [https://nf-co.re/pipelines/](https://nf-co.re/pipelines/) e digite `demo` na barra de pesquisa.

![resultados da pesquisa](./img/search-results.png)

Clique no nome do pipeline, `demo`, para acessar a página de documentação do pipeline.

Cada pipeline lançado tem uma página dedicada que inclui as seguintes seções de documentação:

- **Introduction:** Uma introdução e visão geral do pipeline
- **Usage:** Descrições de como executar o pipeline
- **Parameters:** Parâmetros do pipeline agrupados com descrições
- **Output:** Descrições e exemplos dos arquivos de saída esperados
- **Results:** Exemplos de arquivos de saída gerados a partir do conjunto de dados de teste completo
- **Releases & Statistics:** Histórico de versões e estatísticas do pipeline

Sempre que estiver considerando adotar um novo pipeline, você deve ler a documentação do pipeline com atenção primeiro para entender o que ele faz e como deve ser configurado antes de tentar executá-lo.

Dê uma olhada agora e veja se consegue descobrir:

- Quais ferramentas o pipeline vai executar (Verifique a aba: `Introduction`)
- Quais entradas e parâmetros o pipeline aceita ou requer (Verifique a aba: `Parameters`)
- Quais são as saídas produzidas pelo pipeline (Verifique a aba: `Output`)

#### 1.1.1. Visão geral do pipeline

A aba `Introduction` fornece uma visão geral do pipeline, incluindo uma representação visual (chamada de mapa de metrô) e uma lista de ferramentas que são executadas como parte do pipeline.

![mapa de metrô do pipeline](./img/nf-core-demo-subway-cropped.png)

1. Controle de qualidade das leituras ([FASTQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
2. Corte de adaptadores e qualidade ([SEQTK_TRIM](https://github.com/lh3/seqtk))
3. Apresentação do controle de qualidade das leituras brutas ([MULTIQC](http://multiqc.info/))
4. Geração de uma mensagem de texto divertida de uma vaca ([COWPY](https://github.com/jeffbuttars/cowpy))

#### 1.1.2. Exemplo de linha de comando

A documentação também fornece um arquivo de entrada de exemplo (discutido mais adiante) e um exemplo de linha de comando.

```bash
nextflow run nf-core/demo \
  -profile <docker/singularity/.../institute> \
  --input samplesheet.csv \
  --outdir <OUTDIR>
```

Você vai notar que o comando de exemplo NÃO especifica um arquivo de fluxo de trabalho, apenas a referência ao repositório do pipeline, `nf-core/demo`.

Quando invocado dessa forma, o Nextflow assume que o código está organizado de uma determinada maneira.
Vamos recuperar o código para que possamos examinar essa estrutura.

### 1.2. Recuperar o código do pipeline

Depois de determinarmos que o pipeline parece ser adequado para nossos propósitos, vamos testá-lo.
Felizmente, o Nextflow facilita a recuperação de pipelines de repositórios corretamente formatados sem precisar baixar nada manualmente.

#### 1.2.1. Usar `nextflow pull`

Vamos voltar ao terminal e executar o seguinte:

```bash
nextflow pull nf-core/demo
```

??? success "Saída do comando"

    ```console
    Checking nf-core/demo ...
     downloaded from https://github.com/nf-core/demo.git - revision: 32893afef8 [master]
    ```

O Nextflow faz um `pull` do código do pipeline, ou seja, baixa o repositório completo para o seu disco local.

Para ficar claro, você pode fazer isso com qualquer pipeline Nextflow que esteja configurado adequadamente no GitHub, não apenas pipelines nf-core.
No entanto, o nf-core é a maior coleção open-source de pipelines Nextflow.

#### 1.2.2. Usar `nextflow list`

Você pode pedir ao Nextflow que liste os pipelines que você recuperou dessa forma:

```bash
nextflow list
```

??? success "Saída do comando"

    ```console
    nf-core/demo
    ```

Você pode tentar fazer o pull de alguns outros pipelines para ver como eles aparecem na lista quando você tem mais de um.

#### 1.2.3. Encontrar onde o pipeline foi baixado

Você vai notar que os arquivos não estão no seu diretório de trabalho atual.
Por padrão, o Nextflow salva os pipelines baixados em `$NXF_HOME/assets`.

Para encontrar onde um pipeline específico está localizado, pergunte diretamente ao Nextflow:

```bash
nextflow info nf-core/demo
```

??? success "Saída do comando"

    ```console
     project name: nf-core/demo
     repository  : https://github.com/nf-core/demo
     local path  : /workspaces/.nextflow/assets/.repos/nf-core/demo
     main script : main.nf
     description : An nf-core demo pipeline
     revisions   :
       TEMPLATE
       bumper
       dev
       fix-nxfversion
       manually-merge-3_0_2
     > master (default)
       nf-core-template-merge-2.13.2.dev0
       nf-core-template-merge-2.14.0
       nf-core-template-merge-2.14.1
       nf-core-template-merge-3.0.0
       nf-core-template-merge-3.0.1
       nf-core-template-merge-3.0.2
       nf-core-template-merge-3.1.0
       nf-core-template-merge-3.1.2
       nf-core-template-merge-3.2.0
       nf-core-template-merge-3.2.1
       nf-core-template-merge-3.3.1
       nf-core-template-merge-3.3.2
       nf-core-template-merge-4.0.0
       nf-core-template-merge-4.0.3
       nf-core-template-merge-4.1.0
       nf-core-template-merge-4.1.0-2
       patch
       1.0.0 [t]
       1.0.1 [t]
       1.0.2 [t]
       1.1.0 [t]
     > 1.2.0 [t]
    ```

!!! info "Info"

    O caminho completo pode ser diferente no seu sistema se você não estiver usando nosso ambiente de treinamento.

O Nextflow mantém o código-fonte baixado intencionalmente "fora do caminho", com base no princípio de que esses pipelines devem ser usados mais como bibliotecas do que como código com o qual você interagiria diretamente.

Por baixo dos panos, o Nextflow armazena cada pipeline baixado como um repositório git em `$NXF_HOME/assets/.repos/`, e faz o checkout do código de cada revisão em um subdiretório `clones/<commit>/`.
Como `.repos` é um diretório oculto, um simples `tree -L 2 $NXF_HOME/assets/` parecerá vazio.

#### 1.2.4. Criar um link simbólico para acessar o código-fonte facilmente

Não vamos examinar o código em detalhes, mas vamos dar uma rápida olhada apenas para ter uma ideia de como é a organização geral.

Para facilitar a navegação pelo código-fonte do pipeline, crie um link simbólico apontando para a cópia do pipeline que foi feito o checkout:

```bash
mkdir -p pipelines/nf-core
ln -s "$(echo $NXF_HOME/assets/.repos/nf-core/demo/clones/*/)" pipelines/nf-core/demo
```

Isso cria um atalho para que você possa explorar o código com `tree -L 2 pipelines/nf-core/demo` ou abrir arquivos diretamente.

#### 1.2.5. Visão geral da organização do código

Você pode usar `tree` ou o explorador de arquivos para encontrar e abrir o diretório `nf-core/demo`.

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

    7 directories, 12 files
    ```

Como você pode ver, há muita coisa lá dentro, mas a maior parte não precisa ser sua preocupação.

Resumidamente, vamos observar que no nível superior você pode encontrar um arquivo README com informações resumidas, além de arquivos auxiliares que resumem informações do projeto, como licenciamento, diretrizes de contribuição, citação e código de conduta.
A documentação detalhada do pipeline está localizada no diretório `docs`.
Todo esse conteúdo é usado para gerar as páginas web no site do nf-core de forma programática, portanto, elas estão sempre atualizadas com o código.

Para o restante, podemos distinguir três grupos funcionais de arquivos de código:

1. Componentes de código do pipeline (`main.nf`, `workflows`, `subworkflows`, `modules`)
2. Configuração do pipeline
3. Parâmetros / entradas e validação do pipeline

Não vamos analisar os componentes de código do pipeline nesta parte do curso, mas vamos abordar elementos de configuração e validação que provavelmente serão relevantes para você como usuário final de pipelines nf-core.

!!! tip "Dica"

    Você também pode navegar pelo código-fonte de qualquer pipeline nf-core no GitHub, por exemplo, [github.com/nf-core/demo](https://github.com/nf-core/demo).
    Todo pipeline nf-core segue o mesmo layout de diretórios, então, uma vez que você conhece a estrutura, pode encontrar arquivos de configuração, módulos e fluxos de trabalho de qualquer pipeline da mesma forma.

Por ora, vamos executar o pipeline!

### Conclusão

Agora você sabe como encontrar um pipeline no site do nf-core e recuperar uma cópia local do código-fonte.

### O que vem a seguir?

Aprenda como testar um pipeline nf-core com o mínimo de esforço.

---

## 2. Testar o pipeline com seu perfil de teste

Convenientemente, todo pipeline nf-core vem com um perfil de teste.
Esse é um conjunto mínimo de configurações para o pipeline ser executado usando um pequeno conjunto de dados de teste hospedado no repositório [nf-core/test-datasets](https://github.com/nf-core/test-datasets).
É uma ótima maneira de testar rapidamente um pipeline em pequena escala.

!!! tip "Dica"

    O sistema de perfis de configuração do Nextflow permite que você alterne facilmente entre diferentes motores de contêiner ou ambientes de execução.
    Para mais detalhes, consulte [Hello Nextflow Parte 6: Configuração](../hello_nextflow/06_hello_config.md).

### 2.1. Examinar o perfil de teste

É uma boa prática verificar o que o perfil de teste de um pipeline especifica antes de executá-lo.
O perfil `test` do `nf-core/demo` está no arquivo de configuração `conf/test.config`.
Você pode encontrá-lo localmente dentro do código-fonte do pipeline que o `nextflow pull` baixou, por meio do link simbólico `pipelines` criado na seção 1.2.4:

```bash
code pipelines/nf-core/demo/conf/test.config
```

Aqui está o conteúdo desse arquivo:

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
        cpus: 2,
        memory: '4.GB',
        time: '1.h',
    ]
}

params {
    config_profile_name        = 'Test profile'
    config_profile_description = 'Minimal test dataset to check pipeline function'

    // Input data
    input                      = 'https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/samplesheet/samplesheet_test_illumina_amplicon.csv'
}
```

Você vai notar imediatamente que o bloco de comentários no topo inclui um exemplo de uso mostrando como executar o pipeline com esse perfil de teste.

```groovy title="conf/test.config" linenums="7"
    Use as follows:
        nextflow run nf-core/demo -profile test,<docker/singularity> --outdir <OUTDIR>
```

As únicas coisas que precisamos fornecer são o que está mostrado entre colchetes angulares no comando de exemplo: `<docker/singularity>` e `<OUTDIR>`.

Como lembrete, `<docker/singularity>` refere-se à escolha do sistema de contêiner. Todos os pipelines nf-core são projetados para serem usados com contêineres (Docker, Singularity, etc.) para garantir reprodutibilidade e eliminar problemas de instalação de software.
Portanto, precisamos especificar se queremos usar Docker ou Singularity para testar o pipeline.

A parte `--outdir <OUTDIR>` refere-se ao diretório onde o Nextflow vai escrever as saídas do pipeline.
Precisamos fornecer um nome para ele, que podemos simplesmente inventar.
Se ele ainda não existir, o Nextflow vai criá-lo para nós em tempo de execução.

Passando para a seção após o bloco de comentários, o perfil de teste nos mostra o que foi pré-configurado para os testes: mais notavelmente, o parâmetro `input` já está definido para apontar para um conjunto de dados de teste, então não precisamos fornecer nossos próprios dados.
Se você seguir o link para a entrada pré-configurada, verá que é um arquivo CSV contendo identificadores de amostras e caminhos de arquivos para várias amostras experimentais.

```csv title="samplesheet_test_illumina_amplicon.csv"
sample,fastq_1,fastq_2
SAMPLE1_PE,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample1_R1.fastq.gz,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample1_R2.fastq.gz
SAMPLE2_PE,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample2_R1.fastq.gz,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample2_R2.fastq.gz
SAMPLE3_SE,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample1_R1.fastq.gz,
SAMPLE3_SE,https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/illumina/amplicon/sample2_R1.fastq.gz,
```

Isso é chamado de samplesheet, e é a forma mais comum de entrada para pipelines nf-core.
Não se preocupe se você não estiver familiarizado com os formatos e tipos de dados, isso não é importante para o que vem a seguir.

Agora temos tudo o que precisamos para testar o pipeline.

### 2.2. Executar o pipeline

Conforme observado acima, podemos usar o comando de teste de exemplo quase como está; só precisamos especificar qual empacotamento de software usar e como nomear o diretório de saída.
Aqui vamos usar Docker como sistema de contêiner e `demo-results`, respectivamente.

Com isso, podemos executar o comando de teste:

```bash
nextflow run nf-core/demo -profile test,docker --outdir demo-results
```

??? success "Saída do comando"

    ```console
     N E X T F L O W   ~  version 26.04.4

    Downloading plugin nf-schema@2.7.2
    Launching `https://github.com/nf-core/demo` [cranky_curry] revision: 32893afef8 [master]


    ------------------------------------------------------
                                            ,--./,-.
            ___     __   __   __   ___     /,-._.--~'
      |\ | |__  __ /  ` /  \ |__) |__         }  {
      | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                            `._,._,'
      nf-core/demo 1.2.0
    ------------------------------------------------------

    Input/output options
      input                     : https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/samplesheet/samplesheet_test_illumina_amplicon.csv
      outdir                    : demo-results

    Institutional config options
      config_profile_name       : Test profile
      config_profile_description: Minimal test dataset to check pipeline function

    Generic options
      trace_report_suffix       : 2026-07-03_21-31-35

    Core Nextflow options
      revision                  : master
      runName                   : cranky_curry
      containerEngine           : docker
      launchDir                 : /workspaces/training/nfcore-use
      workDir                   : /workspaces/training/nfcore-use/work
      projectDir                : /workspaces/.nextflow/assets/.repos/nf-core/demo/clones/32893afef8076a03a2767a020b3f0cab2e0b40b2
      userName                  : root
      profile                   : test,docker
      configFiles               : /workspaces/.nextflow/assets/.repos/nf-core/demo/clones/32893afef8076a03a2767a020b3f0cab2e0b40b2/nextflow.config

    !! Only displaying parameters that differ from the pipeline defaults !!
    ------------------------------------------------------

    * The pipeline
        https://doi.org/10.5281/zenodo.12192442

    * The nf-core framework
        https://doi.org/10.1038/s41587-020-0439-x

    * Software dependencies
        https://github.com/nf-core/demo/blob/master/CITATIONS.md

    executor >  local (8)
    [ca/5b0f3e] NFCORE_DEMO:DEMO:FASTQC (SAMPLE3_SE)     | 3 of 3 ✔
    [b7/cb6812] NFCORE_DEMO:DEMO:SEQTK_TRIM (SAMPLE3_SE) | 3 of 3 ✔
    [ff/6ebd98] NFCORE_DEMO:DEMO:COWPY                   | 1 of 1 ✔
    [09/bbd1b4] NFCORE_DEMO:DEMO:MULTIQC (demo)          | 1 of 1 ✔
    -[nf-core/demo] Pipeline completed successfully-
    ```

Se a sua saída corresponder a essa, parabéns! Você acabou de executar seu primeiro pipeline nf-core.

Você vai notar que há muito mais saída no console do que quando você executa um pipeline Nextflow básico.
Há um cabeçalho que inclui um resumo da versão do pipeline, entradas e saídas, e alguns elementos de configuração.

!!! info "Info"

    Sua saída vai mostrar timestamps, nomes de execução e caminhos de arquivos diferentes, mas a estrutura geral e a execução dos processos devem ser semelhantes.

Observe a linha próxima ao topo da saída:

```console
Launching `https://github.com/nf-core/demo` [cranky_curry] revision: 32893afef8 [master]
```

Isso informa qual revisão do pipeline foi usada.
Como não especificamos uma versão, o Nextflow usou o commit mais recente em `master`.
Para execuções reproduzíveis, você deve fixar uma versão específica usando a flag `-r`:

```bash
nextflow run nf-core/demo -r 1.2.0 -profile test,docker --outdir demo-results
```

Isso garante que o mesmo código do pipeline seja usado sempre, independentemente de novos commits ou lançamentos.
Para este treinamento, omitimos `-r` por simplicidade, mas em produção você deve sempre especificá-lo.

Passando para a saída de execução, vamos dar uma olhada nas linhas que nos dizem quais processos foram executados:

```console
executor >  local (8)
[ca/5b0f3e] NFCORE_DEMO:DEMO:FASTQC (SAMPLE3_SE)     | 3 of 3 ✔
[b7/cb6812] NFCORE_DEMO:DEMO:SEQTK_TRIM (SAMPLE3_SE) | 3 of 3 ✔
[ff/6ebd98] NFCORE_DEMO:DEMO:COWPY                   | 1 of 1 ✔
[09/bbd1b4] NFCORE_DEMO:DEMO:MULTIQC (demo)          | 1 of 1 ✔
-[nf-core/demo] Pipeline completed successfully-
```

Isso nos diz que quatro processos foram executados, correspondendo às quatro ferramentas mostradas na página de documentação do pipeline no site do nf-core: `FASTQC`, `SEQTK_TRIM`, `MULTIQC` e `COWPY`.

Os nomes completos dos processos como mostrados aqui, como `NFCORE_DEMO:DEMO:MULTIQC`, são mais longos do que o que você pode ter visto no material introdutório Hello Nextflow.
Eles incluem os nomes de seus fluxos de trabalho pai e refletem a modularidade do código do pipeline.
Se você quiser aprender a desenvolver pipelines no estilo nf-core, consulte o curso [Build with nf-core](../nfcore_build/index.md).

### 2.3. Examinar as saídas do pipeline

Por fim, vamos dar uma olhada no diretório `demo-results` produzido pelo pipeline.

```bash
tree -L 2 demo-results
```

??? abstract "Conteúdo do diretório"

    ```console
    demo-results
    ├── cowpy
    │   └── cowpy.txt
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
    │   └── multiqc_report.html
    └── pipeline_info
        ├── execution_report_2026-07-03_21-31-35.html
        ├── execution_timeline_2026-07-03_21-31-35.html
        ├── execution_trace_2026-07-03_21-31-35.txt
        ├── nf_core_demo_software_mqc_versions.yml
        ├── params_2026-07-03_21-31-43.json
        └── pipeline_dag_2026-07-03_21-31-35.html

    12 directories, 8 files
    ```

Pode parecer bastante coisa.
Para saber mais sobre as saídas do pipeline `nf-core/demo`, consulte sua [página de documentação](https://nf-co.re/demo/1.2.0/docs/output/).

Neste momento, o que é importante observar é que os resultados estão organizados por módulo, e há adicionalmente um diretório chamado `pipeline_info` contendo vários relatórios com timestamps sobre a execução do pipeline.

Por exemplo, o arquivo `execution_timeline_*` mostra quais processos foram executados, em que ordem e quanto tempo levaram para ser executados:

![relatório de linha do tempo de execução](./img/execution_timeline.png)

!!! info "Info"

    Aqui as tarefas não foram executadas em paralelo porque estamos rodando em uma máquina minimalista no Github Codespaces.
    Para ver essas tarefas sendo executadas em paralelo, tente aumentar a alocação de CPU do seu codespace e os limites de recursos na configuração de teste.

Esses relatórios são gerados automaticamente para todos os pipelines nf-core.

### Conclusão

Você sabe como executar um pipeline nf-core usando seu perfil de teste integrado e onde encontrar suas saídas.

### O que vem a seguir?

Siga para a [Parte 2](./02_configure_execution.md), onde você vai aprender como configurar a execução do pipeline.

---

## Resumo

Nesta parte você aprendeu a:

- Encontrar e recuperar um pipeline nf-core e examinar sua estrutura de código
- Executar um pipeline usando seu perfil de teste integrado
