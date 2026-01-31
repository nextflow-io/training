# Parte 2: Implementação de amostra única

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tradução assistida por IA - [saiba mais e sugira melhorias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Nesta parte do curso, vamos escrever o fluxo de trabalho mais simples possível que envolve todos os comandos que executamos na Parte 1 para automatizar sua execução, e vamos processar apenas uma amostra por vez.

Faremos isso em três etapas:

1. Escrever um fluxo de trabalho de estágio único que executa a etapa de QC inicial
2. Adicionar a remoção de adaptadores e QC pós-remoção
3. Adicionar o alinhamento ao genoma de referência

!!! warning "Pré-requisito"

    Você deve trabalhar na Parte 1 do curso antes de iniciar esta lição.
    Especificamente, trabalhar nas seções 2.1-3 cria o arquivo de índice do genoma (`data/genome_index.tar.gz`) necessário para a etapa de alinhamento nesta lição.

---

## 1. Escrever um fluxo de trabalho de estágio único que executa o QC inicial

Vamos começar escrevendo um fluxo de trabalho simples que executa a ferramenta FastQC em um arquivo FASTQ contendo reads de RNAseq single-end.

Fornecemos um arquivo de fluxo de trabalho, `rnaseq.nf`, que descreve as principais partes do fluxo de trabalho.

```groovy title="rnaseq.nf" linenums="1"
#!/usr/bin/env nextflow

// Declarações de INCLUDE de módulo

/*
 * Pipeline parameters
 */

// Entrada primária

workflow {

    // Cria canal de entrada

    // Chama processos

}
```

Tenha em mente que este código de fluxo de trabalho está correto, mas não é funcional; seu propósito é apenas servir como um esqueleto que você usará para escrever o fluxo de trabalho real.

### 1.1. Criar um diretório para armazenar módulos

Vamos criar módulos independentes para cada processo para facilitar o gerenciamento e reutilização, então vamos criar um diretório para armazená-los.

```bash
mkdir modules
```

### 1.2. Criar um módulo para o processo de coleta de métricas de QC

Vamos criar um arquivo de módulo chamado `modules/fastqc.nf` para abrigar o processo `FASTQC`:

```bash
touch modules/fastqc.nf
```

Abra o arquivo no editor de código e copie o seguinte código nele:

```groovy title="modules/fastqc.nf" linenums="1"
#!/usr/bin/env nextflow

process FASTQC {

    container "community.wave.seqera.io/library/trim-galore:0.6.10--1bf8ca4e1967cd18"
    publishDir "results/fastqc", mode: 'symlink'

    input:
    path reads

    output:
    path "${reads.simpleName}_fastqc.zip", emit: zip
    path "${reads.simpleName}_fastqc.html", emit: html

    script:
    """
    fastqc $reads
    """
}
```

Você deve reconhecer todas as peças do que aprendeu na Parte 1 e Parte 2 desta série de treinamento; a única mudança notável é que desta vez estamos usando `mode: symlink` para a diretiva `publishDir`, e estamos usando um parâmetro para definir o `publishDir`.

!!! note "Nota"

    Embora os arquivos de dados que estamos usando aqui sejam muito pequenos, em genômica eles podem ficar muito grandes. Para fins de demonstração no ambiente de ensino, estamos usando o modo de publicação 'symlink' para evitar cópias desnecessárias de arquivos. Você não deve fazer isso em seus fluxos de trabalho finais, pois perderá resultados quando limpar seu diretório `work`.

### 1.3. Importar o módulo para o arquivo de fluxo de trabalho

Adicione a instrução `include { FASTQC } from './modules/fastqc.nf'` ao arquivo `rnaseq.nf`:

```groovy title="rnaseq.nf" linenums="3"
// Declarações de INCLUDE de módulo
include { FASTQC } from './modules/fastqc.nf'
```

### 1.4. Adicionar uma declaração de entrada

Declare um parâmetro de entrada com um valor padrão:

```groovy title="rnaseq.nf" linenums="10"
params {
    // Entrada primária
    reads: Path = "data/reads/ENCSR000COQ1_1.fastq.gz"
}
```

### 1.5. Criar um canal de entrada no bloco workflow

Use uma factory de canal básica `.fromPath()` para criar o canal de entrada:

```groovy title="rnaseq.nf" linenums="13"
workflow {

    // Cria canal de entrada a partir de um caminho de arquivo
    read_ch = channel.fromPath(params.reads)

    // Chama processos

}
```

### 1.6. Chamar o processo `FASTQC` no canal de entrada

```groovy title="rnaseq.nf" linenums="13"
workflow {

    // Cria canal de entrada a partir de um caminho de arquivo
    read_ch = channel.fromPath(params.reads)

    // Controle de qualidade inicial
    FASTQC(read_ch)

}
```

### 1.7. Executar o fluxo de trabalho para testar se funciona

Poderíamos usar o parâmetro `--reads` para especificar uma entrada da linha de comando, mas durante o desenvolvimento podemos ser preguiçosos e apenas usar o padrão de teste que configuramos.

```bash
nextflow run rnaseq.nf
```

??? success "Saída do comando"

    ```console
    N E X T F L O W   ~  version 24.10.0

    Launching `rnaseq.nf` [fabulous_snyder] DSL2 - revision: 3394c725ee

    executor >  local (1)
    [d6/d94c3a] FASTQC (1) [100%] 1 of 1 ✔
    ```

Isso deve executar muito rapidamente se você trabalhou na Parte 1 e já baixou o contêiner.
Se você pulou essa parte, o Nextflow baixará o contêiner para você; você não precisa fazer nada para que isso aconteça, mas pode precisar esperar até um minuto.

Você pode encontrar as saídas em `results/fastqc` conforme especificado no processo `FASTQC` pela diretiva `publishDir`.

```bash
ls results/fastqc
```

```console title="Saída"
ENCSR000COQ1_1_fastqc.html  ENCSR000COQ1_1_fastqc.zip
```

---

## 2. Adicionar remoção de adaptadores e controle de qualidade pós-remoção

Vamos usar o wrapper Trim_Galore, que agrupa o Cutadapt para a remoção em si e o FastQC para o controle de qualidade pós-remoção.

### 2.1. Criar um módulo para o processo de remoção e QC

Vamos criar um arquivo de módulo chamado `modules/trim_galore.nf` para abrigar o processo `TRIM_GALORE`:

```bash
touch modules/trim_galore.nf
```

Abra o arquivo no editor de código e copie o seguinte código nele:

```groovy title="modules/trim_galore.nf" linenums="1"
#!/usr/bin/env nextflow

process TRIM_GALORE {

    container "community.wave.seqera.io/library/trim-galore:0.6.10--1bf8ca4e1967cd18"
    publishDir "results/trimming", mode: 'symlink'

    input:
    path reads

    output:
    path "${reads.simpleName}_trimmed.fq.gz", emit: trimmed_reads
    path "${reads}_trimming_report.txt", emit: trimming_reports
    path "${reads.simpleName}_trimmed_fastqc.{zip,html}", emit: fastqc_reports

    script:
    """
    trim_galore --fastqc $reads
    """
}
```

### 2.2. Importar o módulo para o arquivo de fluxo de trabalho

Adicione a instrução `include { TRIM_GALORE } from './modules/trim_galore.nf'` ao arquivo `rnaseq.nf`:

```groovy title="rnaseq.nf" linenums="3"
// Declarações de INCLUDE de módulo
include { FASTQC } from './modules/fastqc.nf'
include { TRIM_GALORE } from './modules/trim_galore.nf'
```

### 2.3. Chamar o processo no canal de entrada

```groovy title="rnaseq.nf" linenums="14"
workflow {

    // Cria canal de entrada a partir de um caminho de arquivo
    read_ch = channel.fromPath(params.reads)

    // Controle de qualidade inicial
    FASTQC(read_ch)

    // Corte de adaptador e QC pós-corte
    TRIM_GALORE(read_ch)
}
```

### 2.4. Executar o fluxo de trabalho para testar se funciona

```bash
nextflow run rnaseq.nf
```

??? success "Saída do comando"

    ```console
    N E X T F L O W   ~  version 24.10.0

    Launching `rnaseq.nf` [fabulous_snyder] DSL2 - revision: 3394c725ee

    executor >  local (1)
    [d6/d94c3a] FASTQC (1) [100%] 1 of 1 ✔
    [c2/e4a9bb] TRIM_GALORE (1)  [100%] 1 of 1 ✔
    ```

Isso também deve executar muito rapidamente, já que estamos executando em um arquivo de entrada tão pequeno.

Você pode encontrar as saídas em `results/trimming` conforme especificado no processo `TRIM_GALORE` pela diretiva `publishDir`.

```bash
ls results/trimming
```

```console title="Saída"
ENCSR000COQ1_1.fastq.gz_trimming_report.txt  ENCSR000COQ1_1_trimmed_fastqc.zip
ENCSR000COQ1_1_trimmed_fastqc.html           ENCSR000COQ1_1_trimmed.fq.gz
```

---

## 3. Alinhar os reads ao genoma de referência

Finalmente podemos executar a etapa de alinhamento do genoma usando o Hisat2, que também emitirá métricas de controle de qualidade no estilo FastQC.

### 3.1. Criar um módulo para o processo HiSat2

Vamos criar um arquivo de módulo chamado `modules/hisat2_align.nf` para abrigar o processo `HISAT2_ALIGN`:

```bash
touch modules/hisat2_align.nf
```

Abra o arquivo no editor de código e copie o seguinte código nele:

```groovy title="modules/hisat2_align.nf" linenums="1"
#!/usr/bin/env nextflow

process HISAT2_ALIGN {

    container "community.wave.seqera.io/library/hisat2_samtools:5e49f68a37dc010e"
    publishDir "results/align", mode: 'symlink'

    input:
    path reads
    path index_zip

    output:
    path "${reads.simpleName}.bam", emit: bam
    path "${reads.simpleName}.hisat2.log", emit: log

    script:
    """
    tar -xzvf $index_zip
    hisat2 -x ${index_zip.simpleName} -U $reads \
        --new-summary --summary-file ${reads.simpleName}.hisat2.log | \
        samtools view -bS -o ${reads.simpleName}.bam
    """
}
```

### 3.2. Importar o módulo para o arquivo de fluxo de trabalho

Adicione a instrução `include { HISAT2_ALIGN } from './modules/hisat2_align.nf'` ao arquivo `rnaseq.nf`:

```groovy title="rnaseq.nf" linenums="3"
// Declarações de INCLUDE de módulo
include { FASTQC } from './modules/fastqc.nf'
include { TRIM_GALORE } from './modules/trim_galore.nf'
include { HISAT2_ALIGN } from './modules/hisat2_align.nf'
```

### 3.3. Adicionar uma declaração de parâmetro para fornecer o índice do genoma

Declare um parâmetro de entrada com um valor padrão:

```groovy title="rnaseq.nf" linenums="8"
params {
    // Entrada primária
    reads: Path = "data/reads/ENCSR000COQ1_1.fastq.gz"

    // Arquivo do genoma de referência
    hisat2_index_zip: Path = "data/genome_index.tar.gz"
}
```

### 3.4. Chamar o processo `HISAT2_ALIGN` nos reads processados pela saída de `TRIM_GALORE`

Os reads processados estão no canal de saída `TRIM_GALORE.out.trimmed_reads` produzido pela etapa anterior.

Além disso, usamos `file (params.hisat2_index_zip)` para fornecer à ferramenta Hisat2 o arquivo tarball compactado do índice do genoma.

```groovy title="rnaseq.nf" linenums="16"
workflow {

    // Cria canal de entrada a partir de um caminho de arquivo
    read_ch = channel.fromPath(params.reads)

    // Controle de qualidade inicial
    FASTQC(read_ch)

    // Corte de adaptador e QC pós-corte
    TRIM_GALORE(read_ch)

    // Alinhamento a um genoma de referência
    HISAT2_ALIGN(TRIM_GALORE.out.trimmed_reads, file (params.hisat2_index_zip))
}
```

### 3.5. Executar o fluxo de trabalho para testar se funciona

```bash
nextflow run rnaseq.nf
```

??? success "Saída do comando"

    ```console
    N E X T F L O W   ~  version 24.10.0

    Launching `rnaseq.nf` [extravagant_khorana] DSL2 - revision: 701b41bd16

    executor >  local (3)
    [e4/d15ad4] FASTQC (1)       [100%] 1 of 1 ✔
    [c6/12b2be] TRIM_GALORE (1)  [100%] 1 of 1 ✔
    [c6/7a9f13] HISAT2_ALIGN (1) [100%] 1 of 1 ✔
    ```

Você pode encontrar as saídas em `results/align` conforme especificado no processo `HISAT2_ALIGN` pela diretiva `publishDir`.

```bash
ls results/align
```

```console title="Saída"
ENCSR000COQ1_1_trimmed.bam  ENCSR000COQ1_1_trimmed.hisat2.log
```

Isso completa o processamento básico que precisamos aplicar a cada amostra.

_Vamos adicionar a agregação de relatórios MultiQC na Parte 2, depois de fazer o fluxo de trabalho aceitar várias amostras de uma vez._

---

### Conclusão

Você sabe como envolver todas as etapas principais para processar amostras de RNAseq single-end individualmente.

### O que vem a seguir?

Aprenda como modificar o fluxo de trabalho para processar várias amostras em paralelo, agregar relatórios de QC em todas as etapas para todas as amostras e permitir a execução do fluxo de trabalho em dados de RNAseq paired-end.
