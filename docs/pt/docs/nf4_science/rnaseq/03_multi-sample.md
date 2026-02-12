# Parte 3: Implementação de múltiplas amostras com dados paired-end

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tradução assistida por IA - [saiba mais e sugira melhorias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Anteriormente, você construiu um pipeline de chamada de variantes por amostra que processava os dados de cada amostra independentemente.
Nesta parte do curso, vamos levar nosso fluxo de trabalho simples ao próximo nível, transformando-o em uma poderosa ferramenta de automação em lote para lidar com números arbitrários de amostras.
E enquanto fazemos isso, também vamos atualizá-lo para esperar dados paired-end, que são mais comuns em estudos mais recentes.

??? info "Como começar a partir desta seção"

    Esta seção do curso assume que você completou a [Parte 1: Visão Geral do Método](./01_method.md), [Parte 2: Implementação de amostra única](./02_single-sample.md) e tem um pipeline `rnaseq.nf` funcional com arquivos de módulo preenchidos.

    Se você não completou a Parte 2 ou quer começar do zero para esta parte, você pode usar a solução da Parte 2 como seu ponto de partida.
    Execute estes comandos de dentro do diretório `nf4-science/rnaseq/`:

    ```bash
    cp solutions/part2/rnaseq-2.nf rnaseq.nf
    cp solutions/part2/modules/fastqc.nf modules/
    cp solutions/part2/modules/trim_galore.nf modules/
    cp solutions/part2/modules/hisat2_align.nf modules/
    cp solutions/part2/nextflow.config .
    ```

    Isso lhe dá um fluxo de trabalho completo de processamento de amostra única.
    Você pode testar se ele executa com sucesso:

    ```bash
    nextflow run rnaseq.nf -profile test
    ```

## Tarefa

Nesta parte do curso, vamos estender o fluxo de trabalho para fazer o seguinte:

1. Ler informações de amostra de uma planilha CSV
2. Executar QC por amostra, corte e alinhamento em todas as amostras em paralelo
3. Agregar todos os relatórios de QC em um relatório MultiQC abrangente

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/rnaseq/img/rnaseq-wf-03.svg"
</figure>

Isso automatiza as etapas da segunda seção da [Parte 1: Visão Geral do Método](./01_method.md#2-multi-sample-qc-aggregation), onde você executou esses comandos manualmente em seus contêineres.

## Plano de aula

Dividimos isso em três estágios:

1. **Fazer o fluxo de trabalho aceitar múltiplas amostras de entrada.**
   Isso cobre a mudança de um único caminho de arquivo para uma planilha CSV, analisando-a com `splitCsv()`, e executando todos os processos existentes em múltiplas amostras.
2. **Adicionar geração abrangente de relatório de QC.**
   Isso introduz o operador `collect()` para agregar saídas entre amostras, e adiciona um processo MultiQC para produzir um relatório combinado.
3. **Mudar para dados de RNAseq paired-end.**
   Isso cobre a adaptação de processos para entradas paired-end (usando tuplas), criação de módulos paired-end, e configuração de um perfil de teste separado.

Isso implementa o método descrito na [Parte 1: Visão Geral do Método](./01_method.md) (segunda seção cobrindo o caso de uso de múltiplas amostras) e se baseia diretamente no fluxo de trabalho produzido pela Parte 2.

!!! tip "Dica"

     Certifique-se de estar no diretório de trabalho correto:
     `cd /workspaces/training/nf4-science/rnaseq`

---

## 1. Fazer o fluxo de trabalho aceitar múltiplas amostras de entrada

Para executar em múltiplas amostras, precisamos mudar como gerenciamos a entrada: em vez de fornecer um único caminho de arquivo, vamos ler informações de amostra de um arquivo CSV.

Fornecemos um arquivo CSV contendo IDs de amostra e caminhos de arquivos FASTQ no diretório `data/`.

```csv title="data/single-end.csv" linenums="1"
sample_id,fastq_path
ENCSR000COQ1,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COQ1_1.fastq.gz
ENCSR000COQ2,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COQ2_1.fastq.gz
ENCSR000COR1,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COR1_1.fastq.gz
ENCSR000COR2,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COR2_1.fastq.gz
ENCSR000CPO1,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000CPO1_1.fastq.gz
ENCSR000CPO2,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000CPO2_1.fastq.gz
```

Este arquivo CSV inclui uma linha de cabeçalho que nomeia as colunas.

Note que estes ainda são dados de leitura single-end.

!!! warning "Aviso"

    Os caminhos de arquivo no CSV são caminhos absolutos que devem corresponder ao seu ambiente.
    Se você não está executando isso no ambiente de treinamento que fornecemos, você precisará atualizar os caminhos para corresponder ao seu sistema.

### 1.1. Mudar a entrada primária para um CSV de caminhos de arquivos no perfil de teste

Primeiro, precisamos atualizar o perfil de teste em `nextflow.config` para fornecer o caminho do arquivo CSV em vez do caminho FASTQ único.

=== "Depois"

    ```groovy title="nextflow.config" linenums="1" hl_lines="5"
    docker.enabled = true

    profiles {
        test {
            params.input = "${projectDir}/data/single-end.csv"
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
            params.hisat2_index_zip = "${projectDir}/data/genome_index.tar.gz"
        }
    }
    ```

Em seguida, precisaremos atualizar a criação do canal para ler deste CSV.

### 1.2. Atualizar a fábrica de canal para analisar entrada CSV

Precisamos carregar o conteúdo do arquivo no canal em vez de apenas o caminho do arquivo em si.

Podemos fazer isso usando o mesmo padrão que usamos na [Parte 2 do Hello Nextflow](../../hello_nextflow/02_hello_channels.md#42-use-the-splitcsv-operator-to-parse-the-file): aplicando o operador [`splitCsv()`](https://nextflow.io/docs/latest/reference/operator.html#splitcsv) para analisar o arquivo, depois uma operação `map` para extrair o caminho do arquivo FASTQ de cada linha.

=== "Depois"

    ```groovy title="rnaseq.nf" linenums="19" hl_lines="4-7"
    workflow {

        main:
        // Cria canal de entrada a partir do conteúdo de um arquivo CSV
        read_ch = channel.fromPath(params.input)
            .splitCsv(header: true)
            .map { row -> file(row.fastq_path) }
    ```

=== "Antes"

    ```groovy title="rnaseq.nf" linenums="19"
    workflow {

        main:
        // Cria canal de entrada a partir de um caminho de arquivo
        read_ch = channel.fromPath(params.input)
    ```

Uma coisa que é nova comparada ao que você encontrou no curso Hello Nextflow é que este CSV tem uma linha de cabeçalho, então adicionamos `#!groovy header: true` à chamada `splitCsv()`.
Isso nos permite referenciar colunas por nome na operação `map`: `#!groovy row.fastq_path` extrai o caminho do arquivo da coluna `fastq_path` de cada linha.

O tratamento de entrada está atualizado e o fluxo de trabalho está pronto para testar.

### 1.3. Executar o fluxo de trabalho

O fluxo de trabalho agora lê informações de amostra de um arquivo CSV e processa todas as amostras em paralelo.

```bash
nextflow run rnaseq.nf -profile test
```

??? success "Saída do comando"

    ```console
    N E X T F L O W   ~  version 24.10.0

    Launching `rnaseq.nf` [golden_curry] DSL2 - revision: 2a5ba5be1e

    executor >  local (18)
    [07/3ff9c5] FASTQC (6)       [100%] 6 of 6 ✔
    [cc/16859f] TRIM_GALORE (6)  [100%] 6 of 6 ✔
    [68/4c27b5] HISAT2_ALIGN (6) [100%] 6 of 6 ✔
    ```

Desta vez cada etapa é executada 6 vezes, uma vez para cada amostra no arquivo CSV.

Foi só isso que precisamos para fazer o fluxo de trabalho executar em múltiplos arquivos.
O Nextflow lida com todo o paralelismo para nós.

### Conclusão

Você sabe como mudar de uma entrada de arquivo único para entrada de múltiplas amostras baseada em CSV que o Nextflow processa em paralelo.

### O que vem a seguir?

Adicionar uma etapa de agregação de relatório de QC que combina métricas de todas as amostras.

---

## 2. Agregar métricas de QC de pré-processamento em um único relatório MultiQC

Tudo isso produz muitos relatórios de QC, e não queremos ter que vasculhar relatórios individuais.
Este é o ponto perfeito para colocar uma etapa de agregação de relatório MultiQC.

Lembre-se do comando `multiqc` da [Parte 1](01_method.md):

```bash
multiqc . -n <output_name>.html
```

O comando escaneia o diretório atual em busca de arquivos de saída de QC reconhecidos e os agrega em um único relatório HTML.
O URI do contêiner era `community.wave.seqera.io/library/pip_multiqc:a3c26f6199d64b7c`.

Precisamos configurar um parâmetro adicional, preparar as entradas, escrever o processo, conectá-lo, e atualizar o tratamento de saída.

### 2.1. Configurar as entradas

O processo MultiQC precisa de um parâmetro de nome de relatório e as saídas de QC coletadas de todas as etapas anteriores agrupadas juntas.

#### 2.1.1. Adicionar um parâmetro `report_id`

Adicione um parâmetro para nomear o relatório de saída.

=== "Depois"

    ```groovy title="rnaseq.nf" linenums="11" hl_lines="8-9"
    params {
        // Entrada primária
        input: Path

        // Arquivo do genoma de referência
        hisat2_index_zip: Path

        // ID do relatório
        report_id: String
    }
    ```

=== "Antes"

    ```groovy title="rnaseq.nf" linenums="11"
    params {
        // Entrada primária
        input: Path

        // Arquivo do genoma de referência
        hisat2_index_zip: Path
    }
    ```

Adicione o padrão do ID do relatório ao perfil de teste:

=== "Depois"

    ```groovy title="nextflow.config" linenums="1" hl_lines="7"
    docker.enabled = true

    profiles {
        test {
            params.input = "${projectDir}/data/single-end.csv"
            params.hisat2_index_zip = "${projectDir}/data/genome_index.tar.gz"
            params.report_id = "all_single-end"
        }
    }
    ```

=== "Antes"

    ```groovy title="nextflow.config" linenums="1"
    docker.enabled = true

    profiles {
        test {
            params.input = "${projectDir}/data/single-end.csv"
            params.hisat2_index_zip = "${projectDir}/data/genome_index.tar.gz"
        }
    }
    ```

Em seguida, precisaremos preparar as entradas para o processo MultiQC.

#### 2.1.2. Coletar e combinar saídas de QC das etapas anteriores

Precisamos dar ao processo `MULTIQC` todas as saídas relacionadas a QC das etapas anteriores agrupadas juntas.

Para isso, usamos o operador `.mix()`, que agrega múltiplos canais em um único.
Começamos de `channel.empty()` e misturamos todos os canais de saída que queremos combinar.
Isso é mais limpo do que encadear `.mix()` diretamente em um dos canais de saída, porque trata todas as entradas simetricamente.

No nosso fluxo de trabalho, as saídas relacionadas a QC para agregar são:

- `FASTQC.out.zip`
- `FASTQC.out.html`
- `TRIM_GALORE.out.trimming_reports`
- `TRIM_GALORE.out.fastqc_reports`
- `HISAT2_ALIGN.out.log`

Misturamos elas em um único canal, depois usamos `.collect()` para agregar os relatórios de todas as amostras em uma única lista.

Adicione estas linhas ao corpo do fluxo de trabalho após a chamada `HISAT2_ALIGN`:

=== "Depois"

    ```groovy title="rnaseq.nf" linenums="38" hl_lines="4-13"
        // Alinhamento a um genoma de referência
        HISAT2_ALIGN(TRIM_GALORE.out.trimmed_reads, file(params.hisat2_index_zip))

        // Geração de relatório de QC abrangente
        multiqc_files_ch = channel.empty().mix(
            FASTQC.out.zip,
            FASTQC.out.html,
            TRIM_GALORE.out.trimming_reports,
            TRIM_GALORE.out.fastqc_reports,
            HISAT2_ALIGN.out.log,
        )
        multiqc_files_list = multiqc_files_ch.collect()
    ```

=== "Antes"

    ```groovy title="rnaseq.nf" linenums="38"
        // Alinhamento a um genoma de referência
        HISAT2_ALIGN(TRIM_GALORE.out.trimmed_reads, file(params.hisat2_index_zip))
    ```

Usar variáveis intermediárias torna cada etapa clara: `multiqc_files_ch` contém todos os arquivos de QC individuais misturados em um canal, e `multiqc_files_list` é o pacote coletado pronto para passar ao MultiQC.

### 2.2. Escrever o processo de agregação de QC e chamá-lo no fluxo de trabalho

Como antes, precisamos preencher a definição do processo, importar o módulo e adicionar a chamada do processo.

#### 2.2.1. Preencher o módulo para o processo de agregação de QC

Abra `modules/multiqc.nf` e examine o esboço da definição do processo.

Vá em frente e preencha a definição do processo por si mesmo usando as informações fornecidas acima, depois verifique seu trabalho contra a solução na aba "Depois" abaixo.

=== "Antes"

    ```groovy title="modules/multiqc.nf" linenums="1"
    #!/usr/bin/env nextflow

    /*
     * Agregar relatórios de QC com MultiQC
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

=== "Depois"

    ```groovy title="modules/multiqc.nf" linenums="1" hl_lines="8 11 12 15 16 20"
    #!/usr/bin/env nextflow

    /*
     * Agregar relatórios de QC com MultiQC
     */
    process MULTIQC {

        container "community.wave.seqera.io/library/pip_multiqc:a3c26f6199d64b7c"

        input:
        path '*'
        val output_name

        output:
        path "${output_name}.html", emit: report
        path "${output_name}_data", emit: data

        script:
        """
        multiqc . -n ${output_name}.html
        """
    }
    ```

Este processo usa `#!groovy path '*'` como o qualificador de entrada para os arquivos de QC.
O curinga `'*'` diz ao Nextflow para preparar todos os arquivos coletados no diretório de trabalho sem exigir nomes específicos.
A entrada `val output_name` é uma string que controla o nome do arquivo do relatório.

O comando `multiqc .` escaneia o diretório atual (onde todos os arquivos de QC preparados estão) e gera o relatório.

Uma vez que você completou isso, o processo está pronto para usar.

#### 2.2.2. Incluir o módulo

Adicione a instrução de importação a `rnaseq.nf`:

=== "Depois"

    ```groovy title="rnaseq.nf" linenums="3" hl_lines="5"
    // Declarações de INCLUDE de módulo
    include { FASTQC } from './modules/fastqc.nf'
    include { TRIM_GALORE } from './modules/trim_galore.nf'
    include { HISAT2_ALIGN } from './modules/hisat2_align.nf'
    include { MULTIQC } from './modules/multiqc.nf'
    ```

=== "Antes"

    ```groovy title="rnaseq.nf" linenums="3"
    // Declarações de INCLUDE de módulo
    include { FASTQC } from './modules/fastqc.nf'
    include { TRIM_GALORE } from './modules/trim_galore.nf'
    include { HISAT2_ALIGN } from './modules/hisat2_align.nf'
    ```

Agora adicione a chamada do processo ao fluxo de trabalho.

#### 2.2.3. Adicionar a chamada do processo

Passe os arquivos de QC coletados e o ID do relatório ao processo `MULTIQC`:

=== "Depois"

    ```groovy title="rnaseq.nf" linenums="49" hl_lines="2"
        multiqc_files_list = multiqc_files_ch.collect()
        MULTIQC(multiqc_files_list, params.report_id)
    ```

=== "Antes"

    ```groovy title="rnaseq.nf" linenums="49"
        multiqc_files_list = multiqc_files_ch.collect()
    ```

O processo MultiQC agora está conectado ao fluxo de trabalho.

### 2.3. Atualizar o tratamento de saída

Precisamos adicionar as saídas do MultiQC à declaração de publicação e configurar para onde elas vão.

#### 2.3.1. Adicionar alvos de publicação para as saídas do MultiQC

Adicione as saídas do MultiQC à seção `publish:`:

=== "Depois"

    ```groovy title="rnaseq.nf" linenums="52" hl_lines="9-10"
        publish:
        fastqc_zip = FASTQC.out.zip
        fastqc_html = FASTQC.out.html
        trimmed_reads = TRIM_GALORE.out.trimmed_reads
        trimming_reports = TRIM_GALORE.out.trimming_reports
        trimming_fastqc = TRIM_GALORE.out.fastqc_reports
        bam = HISAT2_ALIGN.out.bam
        align_log = HISAT2_ALIGN.out.log
        multiqc_report = MULTIQC.out.report
        multiqc_data = MULTIQC.out.data
    }
    ```

=== "Antes"

    ```groovy title="rnaseq.nf" linenums="52"
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

Em seguida, precisaremos dizer ao Nextflow onde colocar essas saídas.

#### 2.3.2. Configurar os novos alvos de saída

Adicione entradas para os alvos do MultiQC no bloco `output {}`, publicando-os em um subdiretório `multiqc/`:

=== "Depois"

    ```groovy title="rnaseq.nf" linenums="82" hl_lines="7-12"
        align_log {
            path 'align'
        }
        multiqc_report {
            path 'multiqc'
        }
        multiqc_data {
            path 'multiqc'
        }
    }
    ```

=== "Antes"

    ```groovy title="rnaseq.nf" linenums="82"
        align_log {
            path 'align'
        }
    }
    ```

A configuração de saída está completa.

### 2.4. Executar o fluxo de trabalho

Usamos `-resume` para que as etapas de processamento anteriores sejam armazenadas em cache e apenas a nova etapa MultiQC seja executada.

```bash
nextflow run rnaseq.nf -profile test -resume
```

??? success "Saída do comando"

    ```console
    N E X T F L O W   ~  version 24.10.0

    Launching `rnaseq.nf` [modest_pare] DSL2 - revision: fc724d3b49

    executor >  local (1)
    [07/3ff9c5] FASTQC (6)       [100%] 6 of 6, cached: 6 ✔
    [2c/8d8e1e] TRIM_GALORE (5)  [100%] 6 of 6, cached: 6 ✔
    [a4/7f9c44] HISAT2_ALIGN (6) [100%] 6 of 6, cached: 6 ✔
    [56/e1f102] MULTIQC          [100%] 1 of 1 ✔
    ```

Uma única chamada ao MULTIQC foi adicionada após as chamadas de processo em cache.

Você pode encontrar as saídas do MultiQC no diretório de resultados.

```bash
tree -L 2 results/multiqc
```

```console title="Saída"
results/multiqc
├── all_single-end_data
│   ├── cutadapt_filtered_reads_plot.txt
│   ├── cutadapt_trimmed_sequences_plot_3_Counts.txt
│   ├── cutadapt_trimmed_sequences_plot_3_Obs_Exp.txt
│   ├── fastqc_adapter_content_plot.txt
│   ├── fastqc_overrepresented_sequences_plot.txt
│   ├── fastqc_per_base_n_content_plot.txt
│   ├── fastqc_per_base_sequence_quality_plot.txt
│   ├── fastqc_per_sequence_gc_content_plot_Counts.txt
│   ├── fastqc_per_sequence_gc_content_plot_Percentages.txt
│   ├── fastqc_per_sequence_quality_scores_plot.txt
│   ├── fastqc_sequence_counts_plot.txt
│   ├── fastqc_sequence_duplication_levels_plot.txt
│   ├── fastqc_sequence_length_distribution_plot.txt
│   ├── fastqc-status-check-heatmap.txt
│   ├── fastqc_top_overrepresented_sequences_table.txt
│   ├── hisat2_se_plot.txt
│   ├── multiqc_citations.txt
│   ├── multiqc_cutadapt.txt
│   ├── multiqc_data.json
│   ├── multiqc_fastqc.txt
│   ├── multiqc_general_stats.txt
│   ├── multiqc_hisat2.txt
│   ├── multiqc.log
│   ├── multiqc_software_versions.txt
│   └── multiqc_sources.txt
└── all_single-end.html
```

Esse último arquivo `all_single-end.html` é o relatório agregado completo, convenientemente empacotado em um arquivo HTML fácil de navegar.

### Conclusão

Você sabe como coletar saídas de múltiplos canais, agrupá-las com `.mix()` e `.collect()`, e passá-las para um processo de agregação.

### O que vem a seguir?

Adaptar o fluxo de trabalho para lidar com dados de RNAseq paired-end.

---

## 3. Habilitar o processamento de dados de RNAseq paired-end

Agora nosso fluxo de trabalho só pode lidar com dados de RNAseq single-end.
É cada vez mais comum ver dados de RNAseq paired-end, então queremos ser capazes de lidar com isso.

Fazer o fluxo de trabalho completamente agnóstico do tipo de dado exigiria o uso de recursos de linguagem Nextflow um pouco mais avançados, então não vamos fazer isso aqui, mas podemos fazer uma versão de processamento paired-end para demonstrar o que precisa ser adaptado.

### 3.1. Copiar o fluxo de trabalho e atualizar as entradas

Começamos copiando o arquivo do fluxo de trabalho single-end e atualizando-o para dados paired-end.

#### 3.1.1. Copiar o arquivo do fluxo de trabalho

Crie uma cópia do arquivo do fluxo de trabalho para usar como ponto de partida para a versão paired-end.

```bash
cp rnaseq.nf rnaseq_pe.nf
```

Agora atualize os parâmetros e o tratamento de entrada no novo arquivo.

#### 3.1.2. Adicionar um perfil de teste paired-end

Fornecemos um segundo arquivo CSV contendo IDs de amostra e caminhos de arquivos FASTQ pareados no diretório `data/`.

```csv title="data/paired-end.csv" linenums="1"
sample_id,fastq_1,fastq_2
ENCSR000COQ1,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COQ1_1.fastq.gz,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COQ1_2.fastq.gz
ENCSR000COQ2,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COQ2_1.fastq.gz,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COQ2_2.fastq.gz
ENCSR000COR1,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COR1_1.fastq.gz,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COR1_2.fastq.gz
ENCSR000COR2,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COR2_1.fastq.gz,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COR2_2.fastq.gz
ENCSR000CPO1,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000CPO1_1.fastq.gz,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000CPO1_2.fastq.gz
ENCSR000CPO2,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000CPO2_1.fastq.gz,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000CPO2_2.fastq.gz
```

Adicione um perfil `test_pe` a `nextflow.config` que aponta para este arquivo e usa um ID de relatório paired-end.

=== "Depois"

    ```groovy title="nextflow.config" linenums="1" hl_lines="9-13"
    docker.enabled = true

    profiles {
        test {
            params.input = "${projectDir}/data/single-end.csv"
            params.hisat2_index_zip = "${projectDir}/data/genome_index.tar.gz"
            params.report_id = "all_single-end"
        }
        test_pe {
            params.input = "${projectDir}/data/paired-end.csv"
            params.hisat2_index_zip = "${projectDir}/data/genome_index.tar.gz"
            params.report_id = "all_paired-end"
        }
    }
    ```

=== "Antes"

    ```groovy title="nextflow.config" linenums="1"
    docker.enabled = true

    profiles {
        test {
            params.input = "${projectDir}/data/single-end.csv"
            params.hisat2_index_zip = "${projectDir}/data/genome_index.tar.gz"
            params.report_id = "all_single-end"
        }
    }
    ```

O perfil de teste para dados paired-end está pronto.

#### 3.1.3. Atualizar a fábrica de canal

O operador `.map()` precisa pegar ambos os caminhos de arquivo FASTQ e retorná-los como uma lista.

=== "Depois"

    ```groovy title="rnaseq_pe.nf" linenums="25" hl_lines="4"
        // Cria canal de entrada a partir do conteúdo de um arquivo CSV
        read_ch = channel.fromPath(params.input)
            .splitCsv(header: true)
            .map { row -> [file(row.fastq_1), file(row.fastq_2)] }
    ```

=== "Antes"

    ```groovy title="rnaseq_pe.nf" linenums="25"
        // Cria canal de entrada a partir do conteúdo de um arquivo CSV
        read_ch = channel.fromPath(params.input)
            .splitCsv(header: true)
            .map { row -> file(row.fastq_path) }
    ```

O tratamento de entrada está configurado para dados paired-end.

### 3.2. Adaptar o módulo FASTQC para dados paired-end

Copie o módulo para criar uma versão paired-end:

```bash
cp modules/fastqc.nf modules/fastqc_pe.nf
```

A entrada do processo FASTQC não precisa mudar — quando o Nextflow recebe uma lista de dois arquivos, ele prepara ambos e `reads` se expande para ambos os nomes de arquivo.
A única mudança necessária é no bloco de saída: já que agora obtemos dois relatórios FastQC por amostra, mudamos de padrões baseados em `simpleName` para curingas.

=== "Depois"

    ```groovy title="modules/fastqc_pe.nf" linenums="10" hl_lines="2 3"
        output:
        path "*_fastqc.zip", emit: zip
        path "*_fastqc.html", emit: html
    ```

=== "Antes"

    ```groovy title="modules/fastqc_pe.nf" linenums="10"
        output:
        path "${reads.simpleName}_fastqc.zip", emit: zip
        path "${reads.simpleName}_fastqc.html", emit: html
    ```

Isso generaliza o processo de uma forma que o torna capaz de lidar com dados de RNAseq single-end ou paired-end.

Atualize a importação em `rnaseq_pe.nf` para usar a versão paired-end:

=== "Depois"

    ```groovy title="rnaseq_pe.nf" linenums="4" hl_lines="1"
    include { FASTQC } from './modules/fastqc_pe.nf'
    ```

=== "Antes"

    ```groovy title="rnaseq_pe.nf" linenums="4"
    include { FASTQC } from './modules/fastqc.nf'
    ```

O módulo FASTQC e sua importação estão atualizados para dados paired-end.

### 3.3. Adaptar o módulo TRIM_GALORE para dados paired-end

Copie o módulo para criar uma versão paired-end:

```bash
cp modules/trim_galore.nf modules/trim_galore_pe.nf
```

Este módulo precisa de mudanças mais substanciais:

- A entrada muda de um único caminho para uma tupla de dois caminhos
- O comando adiciona a flag `--paired` e recebe ambos os arquivos de leitura
- A saída muda para refletir os arquivos adicionados e diferentes convenções de nomenclatura do Trim Galore, produzindo relatórios FastQC separados para cada arquivo de leitura

=== "Depois"

    ```groovy title="modules/trim_galore_pe.nf" linenums="8" hl_lines="2 5 7 8 12"
        input:
        tuple path(read1), path(read2)

        output:
        tuple path("*_val_1.fq.gz"), path("*_val_2.fq.gz"), emit: trimmed_reads
        path "*_trimming_report.txt", emit: trimming_reports
        path "*_val_1_fastqc.{zip,html}", emit: fastqc_reports_1
        path "*_val_2_fastqc.{zip,html}", emit: fastqc_reports_2

        script:
        """
        trim_galore --fastqc --paired ${read1} ${read2}
        """
    ```

=== "Antes"

    ```groovy title="modules/trim_galore_pe.nf" linenums="8"
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
    ```

Atualize a importação em `rnaseq_pe.nf`:

=== "Depois"

    ```groovy title="rnaseq_pe.nf" linenums="5" hl_lines="1"
    include { TRIM_GALORE } from './modules/trim_galore_pe.nf'
    ```

=== "Antes"

    ```groovy title="rnaseq_pe.nf" linenums="5"
    include { TRIM_GALORE } from './modules/trim_galore.nf'
    ```

O módulo TRIM_GALORE e sua importação estão atualizados para dados paired-end.

### 3.4. Adaptar o módulo HISAT2_ALIGN para dados paired-end

Copie o módulo para criar uma versão paired-end:

```bash
cp modules/hisat2_align.nf modules/hisat2_align_pe.nf
```

Este módulo precisa de mudanças similares:

- A entrada muda de um único caminho para uma tupla de dois caminhos
- O comando HISAT2 muda de `-U` (não pareado) para argumentos de leitura `-1` e `-2` (pareados)
- Todos os usos de `reads.simpleName` mudam para `read1.simpleName` já que agora referenciamos um membro específico do par

=== "Depois"

    ```groovy title="modules/hisat2_align_pe.nf" linenums="8" hl_lines="2 6 7 12 13 14"
        input:
        tuple path(read1), path(read2)
        path index_zip

        output:
        path "${read1.simpleName}.bam", emit: bam
        path "${read1.simpleName}.hisat2.log", emit: log

        script:
        """
        tar -xzvf ${index_zip}
        hisat2 -x ${index_zip.simpleName} -1 ${read1} -2 ${read2} \
            --new-summary --summary-file ${read1.simpleName}.hisat2.log | \
            samtools view -bS -o ${read1.simpleName}.bam
        """
    ```

=== "Antes"

    ```groovy title="modules/hisat2_align_pe.nf" linenums="8"
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
    ```

Atualize a importação em `rnaseq_pe.nf`:

=== "Depois"

    ```groovy title="rnaseq_pe.nf" linenums="6" hl_lines="1"
    include { HISAT2_ALIGN } from './modules/hisat2_align_pe.nf'
    ```

=== "Antes"

    ```groovy title="rnaseq_pe.nf" linenums="6"
    include { HISAT2_ALIGN } from './modules/hisat2_align.nf'
    ```

O módulo HISAT2_ALIGN e sua importação estão atualizados para dados paired-end.

### 3.5. Atualizar a agregação MultiQC para saídas paired-end

O processo `TRIM_GALORE` paired-end agora produz dois canais de relatório FastQC separados (`fastqc_reports_1` e `fastqc_reports_2`) em vez de um.
Atualize o bloco `.mix()` em `rnaseq_pe.nf` para incluir ambos:

=== "Depois"

    ```groovy title="rnaseq_pe.nf" linenums="40" hl_lines="5 6"
        multiqc_files_ch = channel.empty().mix(
            FASTQC.out.zip,
            FASTQC.out.html,
            TRIM_GALORE.out.trimming_reports,
            TRIM_GALORE.out.fastqc_reports_1,
            TRIM_GALORE.out.fastqc_reports_2,
            HISAT2_ALIGN.out.log,
        )
    ```

=== "Antes"

    ```groovy title="rnaseq_pe.nf" linenums="40"
        multiqc_files_ch = channel.empty().mix(
            FASTQC.out.zip,
            FASTQC.out.html,
            TRIM_GALORE.out.trimming_reports,
            TRIM_GALORE.out.fastqc_reports,
            HISAT2_ALIGN.out.log,
        )
    ```

A agregação MultiQC agora inclui ambos os conjuntos de relatórios FastQC paired-end.

### 3.6. Atualizar o tratamento de saída para saídas paired-end

A seção `publish:` e o bloco `output {}` também precisam refletir os dois canais de relatório FastQC separados do processo `TRIM_GALORE` paired-end.

Atualize a seção `publish:` em `rnaseq_pe.nf`:

=== "Depois"

    ```groovy title="rnaseq_pe.nf" linenums="52" hl_lines="6-7"
        publish:
        fastqc_zip = FASTQC.out.zip
        fastqc_html = FASTQC.out.html
        trimmed_reads = TRIM_GALORE.out.trimmed_reads
        trimming_reports = TRIM_GALORE.out.trimming_reports
        trimming_fastqc_1 = TRIM_GALORE.out.fastqc_reports_1
        trimming_fastqc_2 = TRIM_GALORE.out.fastqc_reports_2
        bam = HISAT2_ALIGN.out.bam
        align_log = HISAT2_ALIGN.out.log
        multiqc_report = MULTIQC.out.report
        multiqc_data = MULTIQC.out.data
    }
    ```

=== "Antes"

    ```groovy title="rnaseq_pe.nf" linenums="52"
        publish:
        fastqc_zip = FASTQC.out.zip
        fastqc_html = FASTQC.out.html
        trimmed_reads = TRIM_GALORE.out.trimmed_reads
        trimming_reports = TRIM_GALORE.out.trimming_reports
        trimming_fastqc = TRIM_GALORE.out.fastqc_reports
        bam = HISAT2_ALIGN.out.bam
        align_log = HISAT2_ALIGN.out.log
        multiqc_report = MULTIQC.out.report
        multiqc_data = MULTIQC.out.data
    }
    ```

Atualize as entradas correspondentes no bloco `output {}`:

=== "Depois"

    ```groovy title="rnaseq_pe.nf" linenums="77" hl_lines="4-9"
        trimming_reports {
            path 'trimming'
        }
        trimming_fastqc_1 {
            path 'trimming'
        }
        trimming_fastqc_2 {
            path 'trimming'
        }
        bam {
            path 'align'
        }
    ```

=== "Antes"

    ```groovy title="rnaseq_pe.nf" linenums="77"
        trimming_reports {
            path 'trimming'
        }
        trimming_fastqc {
            path 'trimming'
        }
        bam {
            path 'align'
        }
    ```

O fluxo de trabalho paired-end agora está totalmente atualizado e pronto para executar.

### 3.7. Executar o fluxo de trabalho

Não usamos `-resume` já que isso não usaria cache, e há duas vezes mais dados para processar do que antes, mas ainda deve ser completado em menos de um minuto.

```bash
nextflow run rnaseq_pe.nf -profile test_pe
```

??? success "Saída do comando"

    ```console
    N E X T F L O W   ~  version 24.10.0

    Launching `rnaseq_pe.nf` [reverent_kare] DSL2 - revision: 9c376cc219

    executor >  local (19)
    [c5/cbde15] FASTQC (5)       [100%] 6 of 6 ✔
    [e4/fa2784] TRIM_GALORE (5)  [100%] 6 of 6 ✔
    [3a/e23049] HISAT2_ALIGN (5) [100%] 6 of 6 ✔
    [e6/a3ccd9] MULTIQC          [100%] 1 of 1 ✔
    ```

Agora temos duas versões ligeiramente divergentes do nosso fluxo de trabalho, uma para dados de leitura single-end e uma para dados paired-end.
O próximo passo lógico seria fazer o fluxo de trabalho aceitar qualquer tipo de dado dinamicamente, o que está fora do escopo deste curso, mas podemos abordar isso em um acompanhamento.

---

### Conclusão

Você sabe como adaptar um fluxo de trabalho de amostra única para paralelizar o processamento de múltiplas amostras, gerar um relatório de QC abrangente e adaptar o fluxo de trabalho para usar dados de leitura paired-end.

### O que vem a seguir?

Dê um grande tapinha nas costas! Você completou o curso Nextflow para RNAseq.

Vá para o [resumo final do curso](./next_steps.md) para revisar o que você aprendeu e descobrir o que vem a seguir.
