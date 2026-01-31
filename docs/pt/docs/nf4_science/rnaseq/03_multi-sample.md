# Parte 3: Implementação de múltiplas amostras com dados paired-end

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tradução assistida por IA - [saiba mais e sugira melhorias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Nesta parte final do curso, vamos levar nosso fluxo de trabalho simples ao próximo nível, transformando-o em uma poderosa ferramenta de automação em lote para lidar com números arbitrários de amostras.
E enquanto fazemos isso, também vamos alterá-lo para esperar dados paired-end, que são mais comuns em estudos mais recentes.

Faremos isso em três etapas:

1. Fazer o fluxo de trabalho aceitar múltiplas amostras de entrada e paralelizar a execução
2. Adicionar geração abrangente de relatório de QC
3. Mudar para dados de RNAseq paired-end

---

## 1. Fazer o fluxo de trabalho aceitar múltiplas amostras de entrada e paralelizar a execução

Vamos precisar mudar como gerenciamos a entrada.

### 1.1. Mudar a entrada primária para ser um CSV de caminhos de arquivos em vez de um único arquivo

Fornecemos um arquivo CSV contendo IDs de amostra e caminhos de arquivos FASTQ no diretório `data/`.
Este arquivo CSV inclui uma linha de cabeçalho.
Note que os caminhos dos arquivos FASTQ são caminhos absolutos.

```csv title="data/single-end.csv" linenums="1"
sample_id,fastq_path
ENCSR000COQ1,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COQ1_1.fastq.gz
ENCSR000COQ2,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COQ2_1.fastq.gz
ENCSR000COR1,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COR1_1.fastq.gz
ENCSR000COR2,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COR2_1.fastq.gz
ENCSR000CPO1,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000CPO1_1.fastq.gz
ENCSR000CPO2,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000CPO2_1.fastq.gz
```

Vamos renomear o parâmetro de entrada primária para `input_csv` e mudar o valor padrão para ser o caminho para o arquivo `single-end.csv`.

```groovy title="rnaseq.nf" linenums="13"
params {
    // Entrada primária
    input_csv: Path = "data/single-end.csv"

    // Arquivo do genoma de referência
    hisat2_index_zip: Path = "data/genome_index.tar.gz"
}
```

### 1.2. Atualizar a fábrica de canal de entrada para lidar com um CSV como entrada

Vamos querer carregar o conteúdo do arquivo no canal em vez de apenas o caminho do arquivo em si, então usamos o operador `.splitCsv()` para analisar o formato CSV, depois o operador `.map()` para pegar a informação específica que queremos (o caminho do arquivo FASTQ).

```groovy title="rnaseq.nf" linenums="16"
    // Cria canal de entrada a partir do conteúdo de um arquivo CSV
    read_ch = channel.fromPath(params.input_csv)
        .splitCsv(header:true)
        .map { row -> file(row.fastq_path) }
```

### 1.3. Executar o fluxo de trabalho para testar se funciona

```bash
nextflow run rnaseq.nf
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

Desta vez vemos cada etapa sendo executada 6 vezes, em cada um dos 6 arquivos de dados que fornecemos.

Foi só isso que precisamos para fazer o fluxo de trabalho executar em múltiplos arquivos!
O Nextflow lida com todo o paralelismo para nós.

---

## 2. Agregar métricas de QC de pré-processamento em um único relatório MultiQC

Tudo isso produz muitos relatórios de QC, e não queremos ter que vasculhar relatórios individuais.
Este é o ponto perfeito para colocar uma etapa de agregação de relatório MultiQC!

### 2.1. Criar um módulo para o processo de agregação de QC

Vamos criar um arquivo de módulo chamado `modules/multiqc.nf` para abrigar o processo `MULTIQC`:

```bash
touch modules/multiqc.nf
```

Abra o arquivo no editor de código e copie o seguinte código nele:

```groovy title="modules/multiqc.nf" linenums="1"
#!/usr/bin/env nextflow

process MULTIQC {

    container "community.wave.seqera.io/library/pip_multiqc:a3c26f6199d64b7c"
    publishDir "results/multiqc", mode: 'symlink'

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

### 2.2. Importar o módulo no arquivo do fluxo de trabalho

Adicione a instrução `include { MULTIQC } from './modules/multiqc.nf'` ao arquivo `rnaseq.nf`:

```groovy title="rnaseq.nf" linenums="3"
// Declarações de INCLUDE de módulo
include { FASTQC } from './modules/fastqc.nf'
include { TRIM_GALORE } from './modules/trim_galore.nf'
include { HISAT2_ALIGN } from './modules/hisat2_align.nf'
include { MULTIQC } from './modules/multiqc.nf'
```

### 2.3. Adicionar um parâmetro `report_id` e dar a ele um padrão sensato

```groovy title="rnaseq.nf" linenums="9"
params {
    // Entrada primária
    input_csv: Path = "data/single-end.csv"

    // Arquivo do genoma de referência
    hisat2_index_zip: Path = "data/genome_index.tar.gz"

    // ID do relatório
    report_id: String = "all_single-end"
}
```

### 2.4. Chamar o processo nas saídas das etapas anteriores

Precisamos dar ao processo `MULTIQC` todas as saídas relacionadas a QC das etapas anteriores.

Para isso, vamos usar o operador `.mix()`, que agrega múltiplos canais em um único.

Se tivéssemos quatro processos chamados A, B, C e D com um canal simples `.out` cada, a sintaxe seria assim: `A.out.mix( B.out, C.out, D.out )`. Como você pode ver, você aplica ao primeiro dos canais que deseja combinar (não importa qual) e apenas adiciona todos os outros, separados por vírgulas, nos parênteses que seguem.

No caso do nosso fluxo de trabalho, temos as seguintes saídas para agregar:

- `FASTQC.out.zip`
- `FASTQC.out.html`
- `TRIM_GALORE.out.trimming_reports`
- `TRIM_GALORE.out.fastqc_reports`
- `HISAT2_ALIGN.out.log`

Então o exemplo de sintaxe se torna:

```groovy title="Aplicando .mix() na chamada MULTIQC"
        FASTQC.out.zip.mix(
        FASTQC.out.html,
        TRIM_GALORE.out.trimming_reports,
        TRIM_GALORE.out.fastqc_reports,
        HISAT2_ALIGN.out.log
        )
```

Isso coletará relatórios de QC por amostra.
Mas como queremos agregá-los em todas as amostras, precisamos adicionar o operador `collect()` para reunir os relatórios de todas as amostras em uma única chamada ao `MULTIQC`.
E também precisamos dar a ele o parâmetro `report_id`.

Isso nos dá o seguinte:

```groovy title="A chamada MULTIQC completa" linenums="33"
    // Geração de relatório de QC abrangente
    MULTIQC(
        FASTQC.out.zip.mix(
        FASTQC.out.html,
        TRIM_GALORE.out.trimming_reports,
        TRIM_GALORE.out.fastqc_reports,
        HISAT2_ALIGN.out.log
        ).collect(),
        params.report_id
    )
```

No contexto do bloco completo do fluxo de trabalho, acaba ficando assim:

```groovy title="rnaseq.nf" linenums="18"
workflow {
    // Cria canal de entrada a partir do conteúdo de um arquivo CSV
    read_ch = channel.fromPath(params.input_csv)
        .splitCsv(header:true)
        .map { row -> file(row.fastq_path) }

    // Controle de qualidade inicial
    FASTQC(read_ch)

    // Corte de adaptador e QC pós-corte
    TRIM_GALORE(read_ch)

    // Alinhamento a um genoma de referência
    HISAT2_ALIGN(TRIM_GALORE.out.trimmed_reads, file (params.hisat2_index_zip))

    // Geração de relatório de QC abrangente
    MULTIQC(
        FASTQC.out.zip.mix(
        FASTQC.out.html,
        TRIM_GALORE.out.trimming_reports,
        TRIM_GALORE.out.fastqc_reports,
        HISAT2_ALIGN.out.log
        ).collect(),
        params.report_id
    )
}
```

### 2.5. Executar o fluxo de trabalho para testar se funciona

```bash
nextflow run rnaseq.nf -resume
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

Desta vez vemos uma única chamada ao MULTIQC adicionada após as chamadas de processo em cache:

Você pode encontrar as saídas em `results/trimming` conforme especificado no processo `TRIM_GALORE` pela diretiva `publishDir`.

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

---

## 3. Habilitar o processamento de dados de RNAseq paired-end

Agora nosso fluxo de trabalho só pode lidar com dados de RNAseq single-end.
É cada vez mais comum ver dados de RNAseq paired-end, então queremos ser capazes de lidar com isso.

Fazer o fluxo de trabalho completamente agnóstico do tipo de dado exigiria o uso de recursos de linguagem Nextflow um pouco mais avançados, então não vamos fazer isso aqui, mas podemos fazer uma versão de processamento paired-end para demonstrar o que precisa ser adaptado.

### 3.1. Fazer uma cópia do fluxo de trabalho chamada `rnaseq_pe.nf`

```bash
cp rnaseq.nf rnaseq_pe.nf
```

### 3.2. Modificar o `input_csv` padrão para apontar para os dados paired-end

Fornecemos um segundo arquivo CSV contendo IDs de amostra e caminhos de arquivos FASTQ pareados no diretório `data/`

```csv title="data/paired-end.csv" linenums="1"
sample_id,fastq_1,fastq_2
ENCSR000COQ1,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COQ1_1.fastq.gz,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COQ1_2.fastq.gz
ENCSR000COQ2,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COQ2_1.fastq.gz,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COQ2_2.fastq.gz
ENCSR000COR1,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COR1_1.fastq.gz,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COR1_2.fastq.gz
ENCSR000COR2,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COR2_1.fastq.gz,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000COR2_2.fastq.gz
ENCSR000CPO1,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000CPO1_1.fastq.gz,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000CPO1_2.fastq.gz
ENCSR000CPO2,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000CPO2_1.fastq.gz,/workspaces/training/nf4-science/rnaseq/data/reads/ENCSR000CPO2_2.fastq.gz
```

Vamos mudar o padrão de `input_csv` para ser o caminho para o arquivo `paired-end.csv`.

```groovy title="rnaseq_pe.nf" linenums="15"
params {
    // Entrada primária
    input_csv: Path = "data/paired-end.csv"

    // Arquivo do genoma de referência
    hisat2_index_zip: Path = "data/genome_index.tar.gz"

    // ID do relatório
    report_id: String = "all_single-end"
}
```

### 3.3. Atualizar a fábrica de canal

Precisamos dizer ao operador `.map()` para pegar ambos os caminhos de arquivo FASTQ agora.

Então `row -> file(row.fastq_path)` se torna `row -> [file(row.fastq_1), file(row.fastq_2)]`

```groovy title="rnaseq_pe.nf" linenums="19"
    // Cria canal de entrada a partir do conteúdo de um arquivo CSV
    read_ch = channel.fromPath(params.input_csv)
        .splitCsv(header:true)
        .map { row -> [file(row.fastq_1), file(row.fastq_2)] }
```

### 3.4. Fazer uma versão paired-end do processo FASTQC

Vamos fazer uma cópia do módulo para que possamos ter ambas as versões à mão.

```bash
cp modules/fastqc.nf modules/fastqc_pe.nf
```

Abra o novo arquivo de módulo `fastqc_pe.nf` no editor de código e faça as seguintes mudanças de código:

- Mude `fastqc $reads` para `fastqc ${reads}` no bloco `script` (linha 17) para que a entrada `reads` seja desempacotada, já que agora é uma tupla de dois caminhos em vez de um único caminho.
- Substitua `${reads.simpleName}` por um curinga (`*`) para evitar ter que lidar com os arquivos de saída individualmente.

```groovy title="modules/fastqc_pe.nf" linenums="8"
    input:
    path reads

    output:
    path "*_fastqc.zip", emit: zip
    path "*_fastqc.html", emit: html

    script:
    """
    fastqc ${reads}
    """
```

Tecnicamente isso generaliza o processo `FASTQC` de uma forma que o torna capaz de lidar com dados de RNAseq single-end ou paired-end.

Finalmente, atualize a instrução de importação do módulo para usar a versão paired-end do módulo.

```groovy title="rnaseq_pe.nf" linenums="4"
include { FASTQC } from './modules/fastqc_pe.nf'
```

### 3.5. Fazer uma versão paired-end do processo TRIM_GALORE

Faça uma cópia do módulo para que possamos ter ambas as versões à mão.

```bash
cp modules/trim_galore.nf modules/trim_galore_pe.nf
```

Abra o novo arquivo de módulo `trim_galore_pe.nf` no editor de código e faça as seguintes mudanças de código:

- Mude a declaração de entrada de `path reads` para `tuple path(read1), path(read2)`
- Atualize o comando no bloco `script`, substituindo `$reads` por `--paired ${read1} ${read2}`
- Atualize as declarações de saída para refletir os arquivos adicionados e diferentes convenções de nomenclatura, usando curingas para evitar ter que listar tudo.

```groovy title="modules/trim_galore_pe.nf" linenums="8"
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

Finalmente, atualize a instrução de importação do módulo para usar a versão paired-end do módulo.

```groovy title="rnaseq_pe.nf" linenums="5"
include { TRIM_GALORE } from './modules/trim_galore_pe.nf'
```

### 3.6. Atualizar a chamada ao processo MULTIQC para esperar dois relatórios de TRIM_GALORE

O processo `TRIM_GALORE` agora produz um canal de saída adicional, então precisamos alimentar isso ao MultiQC.

Substitua `TRIM_GALORE.out.fastqc_reports,` por `TRIM_GALORE.out.fastqc_reports_1,` mais `TRIM_GALORE.out.fastqc_reports_2,`:

```groovy title="rnaseq_pe.nf" linenums="33"
    // Geração de relatório de QC abrangente
    MULTIQC(
        FASTQC.out.zip.mix(
        FASTQC.out.html,
        TRIM_GALORE.out.trimming_reports,
        TRIM_GALORE.out.fastqc_reports_1,
        TRIM_GALORE.out.fastqc_reports_2,
        HISAT2_ALIGN.out.log
        ).collect(),
        params.report_id
    )
```

Já que estamos no MultiQC, vamos também atualizar o padrão do parâmetro `report_id` de `"all_single-end"` para `"all_paired-end"`.

```groovy title="rnaseq_pe.nf" linenums="9"
params {
    // Entrada primária
    input_csv: Path = "data/paired-end.csv"

    // Arquivo do genoma de referência
    hisat2_index_zip: Path = "data/genome_index.tar.gz"

    // ID do relatório
    report_id: String = "all_paired-end"
}
```

### 3.7. Fazer uma versão paired-end do processo HISAT2_ALIGN

Faça uma cópia do módulo para que possamos ter ambas as versões à mão.

```bash
cp modules/hisat2_align.nf modules/hisat2_align_pe.nf
```

Abra o novo arquivo de módulo `hisat2_align_pe.nf` no editor de código e faça as seguintes mudanças de código:

- Mude a declaração de entrada de `path reads` para `tuple path(read1), path(read2)`
- Atualize o comando no bloco `script`, substituindo `-U $reads` por `-1 ${read1} -2 ${read2}`
- Substitua todas as instâncias de `${reads.simpleName}` por `${read1.simpleName}` no comando no bloco `script` bem como nas declarações de saída.

```groovy title="modules/hisat2_align_pe.nf" linenums="8"
    input:
    tuple path(read1), path(read2)
    path index_zip

    output:
    path "${read1.simpleName}.bam", emit: bam
    path "${read1.simpleName}.hisat2.log", emit: log

    script:
    """
    tar -xzvf $index_zip
    hisat2 -x ${index_zip.simpleName} -1 ${read1} -2 ${read2} \
        --new-summary --summary-file ${read1.simpleName}.hisat2.log | \
        samtools view -bS -o ${read1.simpleName}.bam
    """
```

Finalmente, atualize a instrução de importação do módulo para usar a versão paired-end do módulo.

```groovy title="rnaseq_pe.nf" linenums="5"
include { HISAT2_ALIGN } from './modules/hisat2_align_pe.nf'
```

### 3.8. Executar o fluxo de trabalho para testar se funciona

Não usamos `-resume` já que isso não usaria cache, e há duas vezes mais dados para processar do que antes, mas ainda deve ser completado em menos de um minuto.

```bash
nextflow run rnaseq_pe.nf
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

E é isso! Agora temos duas versões ligeiramente divergentes do nosso fluxo de trabalho, uma para dados de leitura single-end e uma para dados paired-end.
O próximo passo lógico seria fazer o fluxo de trabalho aceitar qualquer tipo de dado dinamicamente, o que está fora do escopo deste curso, mas podemos abordar isso em um acompanhamento.

---

### Conclusão

Você sabe como adaptar um fluxo de trabalho de amostra única para paralelizar o processamento de múltiplas amostras, gerar um relatório de QC abrangente e adaptar o fluxo de trabalho para usar dados de leitura paired-end se necessário.

### O que vem a seguir?

Parabéns, você completou o mini-curso Nextflow Para RNAseq! Celebre seu sucesso e faça uma pausa bem merecida!

Em seguida, pedimos que você complete uma pesquisa muito breve sobre sua experiência com este curso de treinamento, então vamos levá-lo a uma página com links para recursos de treinamento adicionais e links úteis.
