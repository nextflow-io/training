# Parte 2: Implementação de amostra única

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tradução assistida por IA - [saiba mais e sugira melhorias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Nesta parte do curso, vamos escrever o fluxo de trabalho mais simples possível que envolve todos os comandos que executamos na Parte 1 para automatizar sua execução, e vamos processar apenas uma amostra por vez.

!!! warning "Pré-requisito"

    Você deve trabalhar na [Parte 1: Visão geral do método](./01_method.md) antes de iniciar esta lição.
    Especificamente, trabalhar na seção 1.2.3 cria o arquivo de índice do genoma (`data/genome_index.tar.gz`) necessário para a etapa de alinhamento nesta lição.

## Tarefa

Nesta parte do curso, vamos desenvolver um fluxo de trabalho que faz o seguinte:

1. Executar controle de qualidade (FastQC) nos reads de entrada
2. Cortar adaptadores e executar QC pós-corte (Trim Galore)
3. Alinhar reads cortados a um genoma de referência (HISAT2)

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/rnaseq/img/rnaseq-wf-02.svg"
</figure>

Isso automatiza as etapas da primeira seção da [Parte 1: Visão geral do método](./01_method.md#1-single-sample-processing), onde você executou esses comandos manualmente em seus contêineres.

Como ponto de partida, fornecemos um arquivo de fluxo de trabalho, `rnaseq.nf`, que descreve as principais partes do fluxo de trabalho, bem como quatro arquivos de módulo no diretório `modules/` (`fastqc.nf`, `trim_galore.nf`, `hisat2_align.nf` e `multiqc.nf`) que descrevem a estrutura de cada processo.

??? full-code "Arquivos de esqueleto"

    ```groovy title="rnaseq.nf"
    #!/usr/bin/env nextflow

    // Declarações de INCLUDE de módulo

    /*
     * Pipeline parameters
     */

    // Entrada primária

    workflow {

        main:
        // Cria canal de entrada

        // Chama processos

        publish:
        // Declara saídas para publicar
    }

    output {
        // Configura destinos de publicação
    }
    ```

    ```groovy title="modules/fastqc.nf"
    #!/usr/bin/env nextflow

    /*
     * Run FastQC on input reads
     */
    process FASTQC {

        container

        input:

        output:

        script:
        """

        """
    }
    ```

    ```groovy title="modules/trim_galore.nf"
    #!/usr/bin/env nextflow

    /*
     * Trim adapters and run post-trimming QC
     */
    process TRIM_GALORE {

        container

        input:

        output:

        script:
        """

        """
    }
    ```

    ```groovy title="modules/hisat2_align.nf"
    #!/usr/bin/env nextflow

    /*
     * Align reads to a reference genome
     */
    process HISAT2_ALIGN {

        container

        input:

        output:

        script:
        """

        """
    }
    ```

    ```groovy title="modules/multiqc.nf"
    #!/usr/bin/env nextflow

    /*
     * Aggregate QC reports with MultiQC
     */
    process MULTIQC {

        container

        input:

        output:

        script:
        """

        """
    }
    ```

Esses arquivos não são funcionais; seu propósito é apenas servir como esqueletos para você preencher com as partes interessantes do código.

## Plano da lição

Para tornar o processo de desenvolvimento mais educacional, dividimos isso em três etapas:

1. **Escrever um fluxo de trabalho de estágio único que executa a etapa de QC inicial.**
   Isso cobre a configuração de um parâmetro CLI, criação de um canal de entrada, escrita de um módulo de processo e configuração de publicação de saída.
2. **Adicionar corte de adaptadores e QC pós-corte.**
   Isso introduz o encadeamento de processos conectando a saída de um processo à entrada de outro.
3. **Adicionar alinhamento ao genoma de referência.**
   Isso cobre o tratamento de entradas de referência adicionais e trabalho com arquivos compactados.

Cada etapa se concentra em um aspecto específico do desenvolvimento de fluxo de trabalho.

!!! tip "Dica"

     Certifique-se de estar no diretório de trabalho correto:
     `cd /workspaces/training/nf4-science/rnaseq`

---

## 1. Escrever um fluxo de trabalho de estágio único que executa o QC inicial

Esta primeira etapa se concentra no básico: carregar um arquivo FASTQ e executar controle de qualidade nele.

Lembre-se do comando `fastqc` da [Parte 1](01_method.md):

```bash
fastqc <reads>
```

O comando recebe um arquivo FASTQ como entrada e produz um relatório de controle de qualidade como um arquivo `.zip` e um resumo `.html`.
O URI do contêiner era `community.wave.seqera.io/library/trim-galore:0.6.10--1bf8ca4e1967cd18`.

Vamos pegar essas informações e envolvê-las no Nextflow em três etapas:

1. Configurar a entrada
2. Escrever o processo de QC e chamá-lo no fluxo de trabalho
3. Configurar o tratamento de saída

### 1.1. Configurar a entrada

Precisamos declarar um parâmetro de entrada, criar um perfil de teste para fornecer um valor padrão conveniente e criar um canal de entrada.

#### 1.1.1. Adicionar uma declaração de parâmetro de entrada

Em `rnaseq.nf`, na seção `Pipeline parameters`, declare um parâmetro chamado `input` com o tipo `Path`.

=== "Depois"

    ```groovy title="rnaseq.nf" linenums="5" hl_lines="4-7"
    /*
     * Pipeline parameters
     */
    params {
        // Entrada primária
        input: Path
    }
    ```

=== "Antes"

    ```groovy title="rnaseq.nf" linenums="5"
    /*
     * Pipeline parameters
     */

    // Entrada primária
    ```

Isso configura o parâmetro CLI, mas não queremos digitar o caminho do arquivo toda vez que executamos o fluxo de trabalho durante o desenvolvimento.
Existem várias opções para fornecer um valor padrão; aqui usamos um perfil de teste.

#### 1.1.2. Criar um perfil de teste com um valor padrão em `nextflow.config`

Um perfil de teste fornece valores padrão convenientes para experimentar um fluxo de trabalho sem especificar entradas na linha de comando.
Esta é uma convenção comum no ecossistema Nextflow (veja [Hello Config](../../hello_nextflow/06_hello_config.md) para mais detalhes).

Adicione um bloco `profiles` ao `nextflow.config` com um perfil `test` que define o parâmetro `input` para um dos arquivos FASTQ de teste.

=== "Depois"

    ```groovy title="nextflow.config" linenums="1" hl_lines="3-7"
    docker.enabled = true

    profiles {
        test {
            params.input = "${projectDir}/data/reads/ENCSR000COQ1_1.fastq.gz"
        }
    }
    ```

=== "Antes"

    ```groovy title="nextflow.config" linenums="1"
    docker.enabled = true
    ```

Aqui, estamos usando `#!groovy ${projectDir}`, uma variável integrada do Nextflow que aponta para o diretório onde o script do fluxo de trabalho está localizado.
Isso facilita a referência a arquivos de dados e outros recursos sem codificar caminhos absolutos.

O parâmetro agora tem um padrão conveniente. Em seguida, precisamos criar um canal a partir dele.

#### 1.1.3. Configurar o canal de entrada

No bloco workflow, crie um canal de entrada a partir do valor do parâmetro usando a factory de canal `.fromPath` (como usado em [Hello Channels](../../hello_nextflow/02_hello_channels.md)).

=== "Depois"

    ```groovy title="rnaseq.nf" linenums="13" hl_lines="4-5"
    workflow {

        main:
        // Cria canal de entrada a partir de um caminho de arquivo
        read_ch = channel.fromPath(params.input)

        // Chama processos

        publish:
        // Declara saídas para publicar
    }
    ```

=== "Antes"

    ```groovy title="rnaseq.nf" linenums="13"
    workflow {

        main:
        // Cria canal de entrada

        // Chama processos

        publish:
        // Declara saídas para publicar
    }
    ```

Em seguida, precisaremos criar o processo para executar QC nesta entrada.

### 1.2. Escrever o processo de QC e chamá-lo no fluxo de trabalho

Precisamos preencher a definição do processo no arquivo de módulo, importá-lo para o fluxo de trabalho usando uma instrução include e chamá-lo na entrada.

#### 1.2.1. Preencher o módulo para o processo de QC

Abra `modules/fastqc.nf` e examine o esboço da definição do processo.
Você deve reconhecer os principais elementos estruturais; caso contrário, considere ler [Hello Nextflow](../../hello_nextflow/01_hello_world.md) para uma revisão.

Vá em frente e preencha a definição do processo por conta própria usando as informações fornecidas acima, depois verifique seu trabalho contra a solução na aba "Depois" abaixo.

=== "Antes"

    ```groovy title="modules/fastqc.nf" linenums="1"
    #!/usr/bin/env nextflow

    /*
     * Run FastQC on input reads
     */
    process FASTQC {

        container

        input:

        output:

        script:
        """

        """
    }
    ```

=== "Depois"

    ```groovy title="modules/fastqc.nf" linenums="1" hl_lines="8 11 14 15 19"
    #!/usr/bin/env nextflow

    /*
     * Run FastQC on input reads
     */
    process FASTQC {

        container "community.wave.seqera.io/library/trim-galore:0.6.10--1bf8ca4e1967cd18"

        input:
        path reads

        output:
        path "${reads.simpleName}_fastqc.zip", emit: zip
        path "${reads.simpleName}_fastqc.html", emit: html

        script:
        """
        fastqc ${reads}
        """
    }
    ```

O acessor `simpleName` remove todas as extensões do nome do arquivo, então `ENCSR000COQ1_1.fastq.gz` se torna `ENCSR000COQ1_1`.
Usamos a sintaxe `emit:` para atribuir nomes a cada canal de saída, o que será útil para conectar saídas ao bloco publish.

Uma vez que você tenha completado isso, o processo está completo.
Para usá-lo no fluxo de trabalho, você precisará importar o módulo e adicionar uma chamada de processo.

#### 1.2.2. Incluir o módulo

Em `rnaseq.nf`, adicione uma instrução `include` para tornar o processo disponível para o fluxo de trabalho:

=== "Depois"

    ```groovy title="rnaseq.nf" linenums="3" hl_lines="2"
    // Declarações de INCLUDE de módulo
    include { FASTQC } from './modules/fastqc.nf'
    ```

=== "Antes"

    ```groovy title="rnaseq.nf" linenums="3"
    // Declarações de INCLUDE de módulo
    ```

O processo agora está disponível no escopo do fluxo de trabalho.

#### 1.2.3. Chamar o processo de QC na entrada

Adicione uma chamada para `FASTQC` no bloco workflow, passando o canal de entrada como argumento.

=== "Depois"

    ```groovy title="rnaseq.nf" linenums="14" hl_lines="7-8"
    workflow {

        main:
        // Cria canal de entrada a partir de um caminho de arquivo
        read_ch = channel.fromPath(params.input)

        // Controle de qualidade inicial
        FASTQC(read_ch)

        publish:
        // Declara saídas para publicar
    }
    ```

=== "Antes"

    ```groovy title="rnaseq.nf" linenums="14"
    workflow {

        main:
        // Cria canal de entrada a partir de um caminho de arquivo
        read_ch = channel.fromPath(params.input)

        // Chama processos

        publish:
        // Declara saídas para publicar
    }
    ```

O fluxo de trabalho agora carrega a entrada e executa o processo de QC nela.
Em seguida, precisamos configurar como a saída é publicada.

### 1.3. Configurar o tratamento de saída

Precisamos declarar quais saídas de processo publicar e especificar onde elas devem ir.

#### 1.3.1. Declarar saídas na seção `publish:`

A seção `publish:` dentro do bloco workflow declara quais saídas de processo devem ser publicadas.
Atribua as saídas de `FASTQC` a destinos nomeados.

=== "Depois"

    ```groovy title="rnaseq.nf" linenums="23" hl_lines="2-3"
        publish:
        fastqc_zip = FASTQC.out.zip
        fastqc_html = FASTQC.out.html
    }
    ```

=== "Antes"

    ```groovy title="rnaseq.nf" linenums="23"
        publish:
        // Declara saídas para publicar
    }
    ```

Em seguida, precisaremos dizer ao Nextflow onde colocar as saídas publicadas.

#### 1.3.2. Configurar os destinos de saída no bloco `output {}`

O bloco `output {}` fica fora do fluxo de trabalho e especifica onde cada destino nomeado é publicado.
Configure ambos os destinos para publicar em um subdiretório `fastqc/`.

=== "Depois"

    ```groovy title="rnaseq.nf" linenums="28" hl_lines="2-7"
    output {
        fastqc_zip {
            path 'fastqc'
        }
        fastqc_html {
            path 'fastqc'
        }
    }
    ```

=== "Antes"

    ```groovy title="rnaseq.nf" linenums="28"
    output {
        // Configura destinos de publicação
    }
    ```

!!! note "Nota"

    Por padrão, o Nextflow publica arquivos de saída como links simbólicos, o que evita duplicação desnecessária.
    Embora os arquivos de dados que estamos usando aqui sejam muito pequenos, em genômica eles podem ficar muito grandes.
    Links simbólicos quebrarão quando você limpar seu diretório `work`, então para fluxos de trabalho de produção você pode querer substituir o modo de publicação padrão para `'copy'`.

### 1.4. Executar o fluxo de trabalho

Neste ponto, temos um fluxo de trabalho de QC de uma etapa que deve ser totalmente funcional.

Executamos com `-profile test` para usar o valor padrão configurado no perfil de teste, evitando a necessidade de escrever o caminho na linha de comando.

```bash
nextflow run rnaseq.nf -profile test
```

??? success "Saída do comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `rnaseq.nf` [mad_lorenz] DSL2 - revision: 5846a164d2

    executor >  local (1)
    [7b/8ee79e] FASTQC (1) | 1 of 1 ✔
    ```

Isso deve executar muito rapidamente se você trabalhou na Parte 1 e já baixou o contêiner.
Se você pulou essa parte, o Nextflow baixará o contêiner para você; você não precisa fazer nada para que isso aconteça, mas pode precisar esperar até um minuto.

Você pode verificar as saídas no diretório results.

```bash
ls results/fastqc
```

```console title="Saída"
ENCSR000COQ1_1_fastqc.html  ENCSR000COQ1_1_fastqc.zip
```

Os relatórios de QC para a amostra agora estão publicados no subdiretório `fastqc/`.

### Conclusão

Você sabe como criar um módulo contendo um processo, importá-lo para um fluxo de trabalho, chamá-lo com um canal de entrada e publicar os resultados usando o bloco de saída no nível do fluxo de trabalho.

### O que vem a seguir?

Adicionar corte de adaptadores com QC pós-corte como uma segunda etapa no fluxo de trabalho.

---

## 2. Adicionar corte de adaptadores e QC pós-corte

Agora que temos o QC inicial em vigor, podemos adicionar a etapa de corte de adaptadores com seu QC pós-corte integrado.

Lembre-se do comando `trim_galore` da [Parte 1](01_method.md):

```bash
trim_galore --fastqc <reads>
```

O comando corta adaptadores de um arquivo FASTQ e executa FastQC na saída cortada.
Ele produz reads cortados, um relatório de corte e relatórios FastQC para os reads cortados.
O URI do contêiner era `community.wave.seqera.io/library/trim-galore:0.6.10--1bf8ca4e1967cd18`.

Só precisamos escrever a definição do processo, importá-lo, chamá-lo no fluxo de trabalho e atualizar o tratamento de saída.

### 2.1. Escrever o processo de corte e chamá-lo no fluxo de trabalho

Como antes, precisamos preencher a definição do processo, importar o módulo e adicionar a chamada do processo.

#### 2.1.1. Preencher o módulo para o processo de corte

Abra `modules/trim_galore.nf` e examine o esboço da definição do processo.

Vá em frente e preencha a definição do processo por conta própria usando as informações fornecidas acima, depois verifique seu trabalho contra a solução na aba "Depois" abaixo.

=== "Antes"

    ```groovy title="modules/trim_galore.nf" linenums="1"
    #!/usr/bin/env nextflow

    /*
     * Trim adapters and run post-trimming QC
     */
    process TRIM_GALORE {

        container

        input:

        output:

        script:
        """

        """
    }
    ```

=== "Depois"

    ```groovy title="modules/trim_galore.nf" linenums="1" hl_lines="8 11 14 15 16 20"
    #!/usr/bin/env nextflow

    /*
     * Trim adapters and run post-trimming QC
     */
    process TRIM_GALORE {

        container "community.wave.seqera.io/library/trim-galore:0.6.10--1bf8ca4e1967cd18"

        input:
        path reads

        output:
        path "${reads.simpleName}_trimmed.fq.gz", emit: trimmed_reads
        path "${reads}_trimming_report.txt", emit: trimming_reports
        path "${reads.simpleName}_trimmed_fastqc.{zip,html}", emit: fastqc_reports

        script:
        """
        trim_galore --fastqc ${reads}
        """
    }
    ```

Este processo tem três saídas nomeadas: os reads cortados que alimentam a etapa de alinhamento, o relatório de corte e os relatórios FastQC pós-corte.
A flag `--fastqc` diz ao Trim Galore para executar automaticamente o FastQC na saída cortada.

#### 2.1.2. Incluir o módulo

Atualize `rnaseq.nf` para importar o novo módulo:

=== "Depois"

    ```groovy title="rnaseq.nf" linenums="3" hl_lines="3"
    // Declarações de INCLUDE de módulo
    include { FASTQC } from './modules/fastqc.nf'
    include { TRIM_GALORE } from './modules/trim_galore.nf'
    ```

=== "Antes"

    ```groovy title="rnaseq.nf" linenums="3"
    // Declarações de INCLUDE de módulo
    include { FASTQC } from './modules/fastqc.nf'
    ```

Em seguida, adicionaremos a chamada do processo ao fluxo de trabalho.

#### 2.1.3. Chamar o processo de corte na entrada

Adicione a chamada do processo no bloco workflow:

=== "Depois"

    ```groovy title="rnaseq.nf" linenums="15" hl_lines="10-11"
    workflow {

        main:
        // Cria canal de entrada a partir de um caminho de arquivo
        read_ch = channel.fromPath(params.input)

        // Controle de qualidade inicial
        FASTQC(read_ch)

        // Corte de adaptador e QC pós-corte
        TRIM_GALORE(read_ch)

        publish:
        fastqc_zip = FASTQC.out.zip
        fastqc_html = FASTQC.out.html
    }
    ```

=== "Antes"

    ```groovy title="rnaseq.nf" linenums="15"
    workflow {

        main:
        // Cria canal de entrada a partir de um caminho de arquivo
        read_ch = channel.fromPath(params.input)

        // Controle de qualidade inicial
        FASTQC(read_ch)

        publish:
        fastqc_zip = FASTQC.out.zip
        fastqc_html = FASTQC.out.html
    }
    ```

O processo de corte agora está conectado ao fluxo de trabalho.

### 2.2. Atualizar o tratamento de saída

Precisamos adicionar as saídas de corte à declaração de publicação e configurar para onde elas vão.

#### 2.2.1. Adicionar destinos de publicação para as saídas de corte

Adicione as saídas de corte à seção `publish:`:

=== "Depois"

    ```groovy title="rnaseq.nf" linenums="27" hl_lines="4-6"
        publish:
        fastqc_zip = FASTQC.out.zip
        fastqc_html = FASTQC.out.html
        trimmed_reads = TRIM_GALORE.out.trimmed_reads
        trimming_reports = TRIM_GALORE.out.trimming_reports
        trimming_fastqc = TRIM_GALORE.out.fastqc_reports
    }
    ```

=== "Antes"

    ```groovy title="rnaseq.nf" linenums="27"
        publish:
        fastqc_zip = FASTQC.out.zip
        fastqc_html = FASTQC.out.html
    }
    ```

Em seguida, precisaremos dizer ao Nextflow onde colocar essas saídas.

#### 2.2.2. Configurar os novos destinos de saída

Adicione entradas para os destinos de corte no bloco `output {}`, publicando-os em um subdiretório `trimming/`:

=== "Depois"

    ```groovy title="rnaseq.nf" linenums="35" hl_lines="8-16"
    output {
        fastqc_zip {
            path 'fastqc'
        }
        fastqc_html {
            path 'fastqc'
        }
        trimmed_reads {
            path 'trimming'
        }
        trimming_reports {
            path 'trimming'
        }
        trimming_fastqc {
            path 'trimming'
        }
    }
    ```

=== "Antes"

    ```groovy title="rnaseq.nf" linenums="35"
    output {
        fastqc_zip {
            path 'fastqc'
        }
        fastqc_html {
            path 'fastqc'
        }
    }
    ```

A configuração de saída está completa.

### 2.3. Executar o fluxo de trabalho

O fluxo de trabalho agora inclui tanto QC inicial quanto corte de adaptadores.

```bash
nextflow run rnaseq.nf -profile test
```

??? success "Saída do comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `rnaseq.nf` [gloomy_becquerel] DSL2 - revision: bb11055736

    executor >  local (2)
    [f6/c8ef2e] FASTQC (1)      | 1 of 1 ✔
    [58/c58d8a] TRIM_GALORE (1) | 1 of 1 ✔
    ```

Isso também deve executar muito rapidamente, já que estamos executando em um arquivo de entrada tão pequeno.

Você pode encontrar as saídas de corte no diretório results.

```bash
ls results/trimming
```

```console title="Saída"
ENCSR000COQ1_1.fastq.gz_trimming_report.txt  ENCSR000COQ1_1_trimmed_fastqc.zip
ENCSR000COQ1_1_trimmed_fastqc.html           ENCSR000COQ1_1_trimmed.fq.gz
```

As saídas de corte e relatórios de QC pós-corte agora estão no subdiretório `trimming/`.

### Conclusão

Você sabe como adicionar uma segunda etapa de processamento que executa independentemente na mesma entrada, produzindo múltiplas saídas nomeadas.

### O que vem a seguir?

Adicionar a etapa de alinhamento que se encadeia a partir da saída de reads cortados.

---

## 3. Adicionar alinhamento ao genoma de referência

Finalmente podemos adicionar a etapa de alinhamento do genoma usando HISAT2.

Lembre-se do comando de alinhamento da [Parte 1](01_method.md):

```bash
hisat2 -x <genome_index> -U <reads> \
    --new-summary --summary-file <reads>.hisat2.log | \
    samtools view -bS -o <reads>.bam
```

O comando alinha reads a um genoma de referência e converte a saída para formato BAM.
Ele requer um arquivo de índice do genoma pré-construído e produz um arquivo BAM e um log de resumo de alinhamento.
O URI do contêiner era `community.wave.seqera.io/library/hisat2_samtools:5e49f68a37dc010e`.

Este processo requer uma entrada adicional (o arquivo de índice do genoma), então precisamos configurar isso primeiro, depois escrever e conectar o processo.

### 3.1. Configurar as entradas

Precisamos declarar um parâmetro para o arquivo de índice do genoma.

#### 3.1.1. Adicionar um parâmetro para o índice do genoma

Adicione uma declaração de parâmetro para o arquivo de índice do genoma em `rnaseq.nf`:

=== "Depois"

    ```groovy title="rnaseq.nf" linenums="11" hl_lines="5-6"
    params {
        // Entrada primária
        input: Path

        // Arquivo do genoma de referência
        hisat2_index_zip: Path
    }
    ```

=== "Antes"

    ```groovy title="rnaseq.nf" linenums="11"
    params {
        // Entrada primária
        input: Path
    }
    ```

#### 3.1.2. Adicionar o padrão do índice do genoma ao perfil de teste

Assim como fizemos para `input` na seção 1.1.2, adicione um valor padrão para o índice do genoma ao perfil de teste em `nextflow.config`:

=== "Depois"

    ```groovy title="nextflow.config" linenums="1" hl_lines="6"
    docker.enabled = true

    profiles {
        test {
            params.input = "${projectDir}/data/reads/ENCSR000COQ1_1.fastq.gz"
            params.hisat2_index_zip = "${projectDir}/data/genome_index.tar.gz"
        }
    }
    ```

=== "Antes"

    ```groovy title="nextflow.config" linenums="1"
    docker.enabled = true

    profiles {
        test {
            params.input = "${projectDir}/data/reads/ENCSR000COQ1_1.fastq.gz"
        }
    }
    ```

O parâmetro está pronto; agora podemos criar o processo de alinhamento.

### 3.2. Escrever o processo de alinhamento e chamá-lo no fluxo de trabalho

Como antes, precisamos preencher a definição do processo, importar o módulo e adicionar a chamada do processo.

#### 3.2.1. Preencher o módulo para o processo de alinhamento

Abra `modules/hisat2_align.nf` e examine o esboço da definição do processo.

Vá em frente e preencha a definição do processo por conta própria usando as informações fornecidas acima, depois verifique seu trabalho contra a solução na aba "Depois" abaixo.

=== "Antes"

    ```groovy title="modules/hisat2_align.nf" linenums="1"
    #!/usr/bin/env nextflow

    /*
     * Align reads to a reference genome
     */
    process HISAT2_ALIGN {

        container

        input:

        output:

        script:
        """

        """
    }
    ```

=== "Depois"

    ```groovy title="modules/hisat2_align.nf" linenums="1" hl_lines="8 11 12 15 16 20 21 22 23"
    #!/usr/bin/env nextflow

    /*
     * Align reads to a reference genome
     */
    process HISAT2_ALIGN {

        container "community.wave.seqera.io/library/hisat2_samtools:5e49f68a37dc010e"

        input:
        path reads
        path index_zip

        output:
        path "${reads.simpleName}.bam", emit: bam
        path "${reads.simpleName}.hisat2.log", emit: log

        script:
        """
        tar -xzvf ${index_zip}
        hisat2 -x ${index_zip.simpleName} -U ${reads} \
            --new-summary --summary-file ${reads.simpleName}.hisat2.log | \
            samtools view -bS -o ${reads.simpleName}.bam
        """
    }
    ```

Este processo recebe duas entradas: os reads e o arquivo de índice do genoma.
O bloco script primeiro extrai o índice do arquivo, depois executa o alinhamento HISAT2 canalizado para `samtools view` para converter a saída para formato BAM.
O acessor `simpleName` em `index_zip` extrai o nome base do arquivo (`genome_index`) para usar como prefixo do índice.

#### 3.2.2. Incluir o módulo

Atualize `rnaseq.nf` para importar o novo módulo:

=== "Depois"

    ```groovy title="rnaseq.nf" linenums="3" hl_lines="4"
    // Declarações de INCLUDE de módulo
    include { FASTQC } from './modules/fastqc.nf'
    include { TRIM_GALORE } from './modules/trim_galore.nf'
    include { HISAT2_ALIGN } from './modules/hisat2_align.nf'
    ```

=== "Antes"

    ```groovy title="rnaseq.nf" linenums="3"
    // Declarações de INCLUDE de módulo
    include { FASTQC } from './modules/fastqc.nf'
    include { TRIM_GALORE } from './modules/trim_galore.nf'
    ```

Em seguida, adicionaremos a chamada do processo ao fluxo de trabalho.

#### 3.2.3. Chamar o processo de alinhamento

Os reads cortados estão no canal de saída `TRIM_GALORE.out.trimmed_reads` produzido pela etapa anterior.
Usamos `#!groovy file(params.hisat2_index_zip)` para fornecer o arquivo de índice do genoma.

=== "Depois"

    ```groovy title="rnaseq.nf" linenums="19" hl_lines="14-15"
    workflow {

        main:
        // Cria canal de entrada a partir de um caminho de arquivo
        read_ch = channel.fromPath(params.input)

        // Controle de qualidade inicial
        FASTQC(read_ch)

        // Corte de adaptador e QC pós-corte
        TRIM_GALORE(read_ch)

        // Alinhamento a um genoma de referência
        HISAT2_ALIGN(TRIM_GALORE.out.trimmed_reads, file(params.hisat2_index_zip))
    ```

=== "Antes"

    ```groovy title="rnaseq.nf" linenums="19"
    workflow {

        main:
        // Cria canal de entrada a partir de um caminho de arquivo
        read_ch = channel.fromPath(params.input)

        // Controle de qualidade inicial
        FASTQC(read_ch)

        // Corte de adaptador e QC pós-corte
        TRIM_GALORE(read_ch)
    ```

O processo de alinhamento agora está conectado ao fluxo de trabalho.

### 3.3. Atualizar o tratamento de saída

Precisamos adicionar as saídas de alinhamento à declaração de publicação e configurar para onde elas vão.

#### 3.3.1. Adicionar destinos de publicação para as saídas de alinhamento

Adicione as saídas de alinhamento à seção `publish:`:

=== "Depois"

    ```groovy title="rnaseq.nf" linenums="34" hl_lines="7-8"
        publish:
        fastqc_zip = FASTQC.out.zip
        fastqc_html = FASTQC.out.html
        trimmed_reads = TRIM_GALORE.out.trimmed_reads
        trimming_reports = TRIM_GALORE.out.trimming_reports
        trimming_fastqc = TRIM_GALORE.out.fastqc_reports
        bam = HISAT2_ALIGN.out.bam
        align_log = HISAT2_ALIGN.out.log
    }
    ```

=== "Antes"

    ```groovy title="rnaseq.nf" linenums="34"
        publish:
        fastqc_zip = FASTQC.out.zip
        fastqc_html = FASTQC.out.html
        trimmed_reads = TRIM_GALORE.out.trimmed_reads
        trimming_reports = TRIM_GALORE.out.trimming_reports
        trimming_fastqc = TRIM_GALORE.out.fastqc_reports
    }
    ```

Em seguida, precisaremos dizer ao Nextflow onde colocar essas saídas.

#### 3.3.2. Configurar os novos destinos de saída

Adicione entradas para os destinos de alinhamento no bloco `output {}`, publicando-os em um subdiretório `align/`:

=== "Depois"

    ```groovy title="rnaseq.nf" linenums="44" hl_lines="17-22"
    output {
        fastqc_zip {
            path 'fastqc'
        }
        fastqc_html {
            path 'fastqc'
        }
        trimmed_reads {
            path 'trimming'
        }
        trimming_reports {
            path 'trimming'
        }
        trimming_fastqc {
            path 'trimming'
        }
        bam {
            path 'align'
        }
        align_log {
            path 'align'
        }
    }
    ```

=== "Antes"

    ```groovy title="rnaseq.nf" linenums="44"
    output {
        fastqc_zip {
            path 'fastqc'
        }
        fastqc_html {
            path 'fastqc'
        }
        trimmed_reads {
            path 'trimming'
        }
        trimming_reports {
            path 'trimming'
        }
        trimming_fastqc {
            path 'trimming'
        }
    }
    ```

A configuração de saída está completa.

### 3.4. Executar o fluxo de trabalho

O fluxo de trabalho agora inclui todas as três etapas de processamento: QC, corte e alinhamento.

```bash
nextflow run rnaseq.nf -profile test
```

??? success "Saída do comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `rnaseq.nf` [elated_stonebraker] DSL2 - revision: e8e57d0cdd

    executor >  local (3)
    [e8/fa29d6] FASTQC (1)       | 1 of 1 ✔
    [ca/ffdde2] TRIM_GALORE (1)  | 1 of 1 ✔
    [b6/1c6ca3] HISAT2_ALIGN (1) | 1 of 1 ✔
    ```

Você pode encontrar as saídas de alinhamento no diretório results.

```bash
ls results/align
```

```console title="Saída"
ENCSR000COQ1_1_trimmed.bam  ENCSR000COQ1_1_trimmed.hisat2.log
```

Isso completa o processamento básico que precisamos aplicar a cada amostra.

_Vamos adicionar a agregação de relatórios MultiQC na Parte 3, depois de fazer o fluxo de trabalho aceitar várias amostras de uma vez._

---

### Conclusão

Você sabe como envolver todas as etapas principais para processar amostras de RNAseq single-end individualmente.

### O que vem a seguir?

Faça uma pausa! Isso foi muito.

Quando estiver se sentindo revigorado, vá para a [Parte 3](./03_multi-sample.md), onde você aprenderá como modificar o fluxo de trabalho para processar várias amostras em paralelo, agregar relatórios de QC em todas as etapas para todas as amostras e permitir a execução do fluxo de trabalho em dados de RNAseq paired-end.
