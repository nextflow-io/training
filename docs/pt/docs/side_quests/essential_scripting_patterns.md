---
title: Padrões Essenciais de Scripting
description: Aprenda técnicas avançadas de programação em Nextflow.
weight: 1200
---

# Padrões Essenciais de Scripting

Esta lição irá guiá-lo através dos padrões essenciais de programação que são fundamentais para criar fluxos de trabalho Nextflow eficazes. Vamos cobrir padrões para lidar com entrada de dados, transformar valores, controlar a lógica do fluxo de trabalho, alocar recursos dinamicamente e mais.

## Objetivos de aprendizado

- Entender as diferenças entre paradigmas de fluxo de dados e scripting
- Aplicar closures, operadores ternários e outras técnicas de Groovy
- Dominar técnicas para manipular metadados e extrair informações de arquivos
- Usar expressões regulares e processamento de strings para analisar nomes de arquivos
- Criar funções reutilizáveis para lógica complexa
- Implementar alocação dinâmica de recursos e estratégias de repetição
- Adicionar lógica condicional para controlar a execução do fluxo de trabalho
- Escrever código robusto usando operadores de navegação segura e Elvis
- Validar entradas com mensagens de erro claras
- Usar manipuladores de eventos para gerenciar a conclusão do fluxo de trabalho

## Pré-requisitos

- Compreensão básica dos fluxos de trabalho Nextflow
- Familiaridade com a sintaxe DSL2
- Conhecimento básico de canais e processos
- Ambiente de desenvolvimento configurado com Nextflow v23.04.0 ou posterior
- Acesso ao Docker (ou Conda) para contêineres de software

## Como começar

Este tutorial assume que você tem conhecimento básico da sintaxe Nextflow e das operações de canal.

Para começar, navegue até a pasta `side-quests/essential_scripting_patterns`:

```bash
cd side-quests/essential_scripting_patterns
```

O repositório contém vários arquivos:

- `main.nf`: O fluxo de trabalho principal
- `modules/fastp.nf`: Módulo para processamento de qualidade com FASTP
- `modules/generate_report.nf`: Módulo para gerar relatórios
- `modules/trimgalore.nf`: Módulo para um programa alternativo de trimming
- `data/samples.csv`: Um CSV com informações de amostras
- `data/sequences/*.fastq`: Arquivos fastq de exemplo

Vamos examinar os arquivos:

```bash
cat main.nf
```

Você deve ver um fluxo de trabalho mínimo que usa inclusões de módulos:

```groovy title="main.nf" linenums="1"
include { FASTP } from './modules/fastp.nf'
include { GENERATE_REPORT } from './modules/generate_report.nf'

workflow {
    ch_samples = channel.fromPath("./data/samples.csv")
        .splitCsv(header: true)
        .map { row ->
            tuple(
                [id: row.sample_id],
                file(row.file_path)
            )
        }

    ch_fastp = FASTP(ch_samples)
    GENERATE_REPORT(ch_samples)
}
```

Os módulos também são simples:

```bash
cat modules/fastp.nf
```

```groovy title="modules/fastp.nf" linenums="1"
process FASTP {
    container 'community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}_trimmed.fastq.gz"), emit: fastq
    tuple val(meta), path("${meta.id}.fastp.json"), emit: json

    script:
    """
    fastp --in1 ${reads} --out1 ${meta.id}_trimmed.fastq.gz --json ${meta.id}.fastp.json --html ${meta.id}.fastp.html
    """
}
```

```bash
cat modules/generate_report.nf
```

```groovy title="modules/generate_report.nf" linenums="1"
process GENERATE_REPORT {
    container 'community.wave.seqera.io/library/bash:5.0'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}_report.txt"), emit: report

    script:
    """
    echo "Sample ID: ${meta.id}" > ${meta.id}_report.txt
    echo "File processed: ${reads}" >> ${meta.id}_report.txt
    echo "Process command: fastp --in1 ${reads} --out1 ${meta.id}_trimmed.fastq.gz" >> ${meta.id}_report.txt
    """
}
```

Finalmente, vamos dar uma olhada nos dados de exemplo:

```bash
cat data/samples.csv
```

```csv title="data/samples.csv"
sample_id,organism,tissue_type,sequencing_depth,quality_score,file_path
SAMPLE_001,human,liver,30000000,38.5,./data/sequences/SAMPLE_001_S1_L001_R1_001.fastq
SAMPLE_002,mouse,brain,25000000,35.2,./data/sequences/SAMPLE_002_S2_L001_R1_001.fastq
SAMPLE_003,human,kidney,45000000,42.1,./data/sequences/SAMPLE_003_S3_L001_R1_001.fastq
```

Vamos executar primeiro o fluxo de trabalho básico para garantir que funcione:

```bash
nextflow run main.nf
```

Você deve ver o fluxo de trabalho executar e completar com sucesso.

Agora vamos começar a melhorar o fluxo de trabalho com padrões de scripting avançados!

---

## 1. Entendendo o Fluxo de Dados vs Scripting

Nextflow utiliza **dois paradigmas de programação distintos**:

1. **Fluxo de dados**: A orquestração de canais através de operadores (`.map`, `.filter`, `.branch`)
2. **Scripting**: Código Groovy executado dentro de closures ou blocos de script em processos

Ambos são cruciais, mas funcionam de maneira diferente, então esta primeira seção esclarecerá as diferenças.

### 1.1. Closures e operadores ternários

Um conceito fundamental no Nextflow é a **closure** - um bloco de código que pode ser passado como um objeto e executado posteriormente. Closures são essenciais para operadores de canal (`map`, `filter`, etc.).

Vamos pegar nosso operador `.map` atual e melhorá-lo para extrair mais metadados:

```groovy
.map { row ->
    tuple(
        [id: row.sample_id],
        file(row.file_path)
    )
}
```

O bloco `{ row -> ... }` é uma closure que recebe um argumento `row`.

Agora, vamos melhorá-lo para extrair mais metadados do CSV:

=== "After"

    ```groovy title="main.nf" linenums="7" hl_lines="2-8"
        .map { row ->
            def sample_meta = [
                id: row.sample_id.toLowerCase(),
                organism: row.organism,
                tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                depth: row.sequencing_depth.toInteger(),
                quality: row.quality_score.toDouble()
            ]
            tuple(sample_meta, file(row.file_path))
        }
    ```

=== "Before"

    ```groovy title="main.nf" linenums="7" hl_lines="2-5"
        .map { row ->
            tuple(
                [id: row.sample_id],
                file(row.file_path)
            )
        }
    ```

Agora, vamos modificar o `main.nf` para usar estes metadados adicionais. Primeiro, vamos mudar a forma como geramos o relatório. Vamos rotular as amostras de alta qualidade como tendo "alta prioridade" usando um **operador ternário**.

=== "After"

    ```groovy title="main.nf" linenums="7" hl_lines="9-10"
        .map { row ->
            def sample_meta = [
                id: row.sample_id.toLowerCase(),
                organism: row.organism,
                tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                depth: row.sequencing_depth.toInteger(),
                quality: row.quality_score.toDouble()
            ]
            def priority = sample_meta.quality > 40 ? 'high' : 'normal'
            tuple(sample_meta + [priority: priority], file(row.file_path))
        }
    ```

=== "Before"

    ```groovy title="main.nf" linenums="7" hl_lines="9"
        .map { row ->
            def sample_meta = [
                id: row.sample_id.toLowerCase(),
                organism: row.organism,
                tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                depth: row.sequencing_depth.toInteger(),
                quality: row.quality_score.toDouble()
            ]
            tuple(sample_meta, file(row.file_path))
        }
    ```

O operador ternário `condition ? valor_se_verdadeiro : valor_se_falso` é uma forma concisa de escrever uma instrução if/else. Se `sample_meta.quality > 40`, então `priority` será `'high'`; caso contrário, será `'normal'`.

Vamos também modificar o processo `GENERATE_REPORT` para incluir os metadados adicionais:

=== "After"

    ```groovy title="modules/generate_report.nf" linenums="11" hl_lines="2-7"
    script:
    """
    echo "Sample ID: ${meta.id}" > ${meta.id}_report.txt
    echo "Organism: ${meta.organism}" >> ${meta.id}_report.txt
    echo "Tissue: ${meta.tissue}" >> ${meta.id}_report.txt
    echo "Sequencing depth: ${meta.depth}" >> ${meta.id}_report.txt
    echo "Quality score: ${meta.quality}" >> ${meta.id}_report.txt
    echo "Priority: ${meta.priority}" >> ${meta.id}_report.txt
    echo "File processed: ${reads}" >> ${meta.id}_report.txt
    echo "Process command: fastp --in1 ${reads} --out1 ${meta.id}_trimmed.fastq.gz" >> ${meta.id}_report.txt
    """
    ```

=== "Before"

    ```groovy title="modules/generate_report.nf" linenums="11" hl_lines="2-4"
    script:
    """
    echo "Sample ID: ${meta.id}" > ${meta.id}_report.txt
    echo "File processed: ${reads}" >> ${meta.id}_report.txt
    echo "Process command: fastp --in1 ${reads} --out1 ${meta.id}_trimmed.fastq.gz" >> ${meta.id}_report.txt
    """
    ```

Agora, vamos executar o fluxo de trabalho atualizado:

```bash
nextflow run main.nf
```

??? success "Saída do comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [jolly_fourier] DSL2 - revision: d3e76a7fce

    executor >  local (6)
    [a3/6b9e80] process > FASTP (sample_003)           [100%] 3 of 3 ✔
    [29/54f4b6] process > GENERATE_REPORT (sample_001) [100%] 3 of 3 ✔
    ```

Podemos verificar o resultado examinando um dos arquivos de relatório:

```console
cat work/29/54f4b6b0eb90fed9e3a673b8e47629/sample_001_report.txt
```

Agora você deve ver todos os metadados incluídos:

```
Sample ID: sample_001
Organism: human
Tissue: liver
Sequencing depth: 30000000
Quality score: 38.5
Priority: normal
File processed: SAMPLE_001_S1_L001_R1_001.fastq
Process command: fastp --in1 SAMPLE_001_S1_L001_R1_001.fastq --out1 sample_001_trimmed.fastq.gz
```

### 1.2. A coleção 'meta' vs 'reads'

No Nextflow, canais, tuplas e coleções (como mapas e listas) são estruturas de dados fundamentais. Existe uma diferença importante entre **coleções no fluxo de dados** (canais/operadores) e **coleções em blocos de script**.

Para demonstrar isso, vamos modificar nosso fluxo de trabalho para extrair metadados dos nomes dos arquivos FASTQ usando expressões regulares. É comum que os nomes de arquivos FASTQ sigam uma convenção como: `SAMPLE_001_S1_L001_R1_001.fastq`.

Neste formato:

- `SAMPLE_001`: ID da amostra
- `S1`: Número da amostra
- `L001`: Número da lane
- `R1`: Número da read
- `001`: Fragmento ou chunk

Vamos extrair esses metadados do nome do arquivo FASTQ:

=== "After"

    ```groovy title="main.nf" linenums="7" hl_lines="12-20"
        .map { row ->
            def sample_meta = [
                id: row.sample_id.toLowerCase(),
                organism: row.organism,
                tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                depth: row.sequencing_depth.toInteger(),
                quality: row.quality_score.toDouble()
            ]
            def priority = sample_meta.quality > 40 ? 'high' : 'normal'
            def fastq_path = file(row.file_path)

            // Extract metadata from filename using regex
            def m = (fastq_path.name =~ /^(.+)_S(\d+)_L(\d{3})_(R[12])_(\d{3})\.fastq(?:\.gz)?$/)
            def file_meta = m ? [
                sample_num: m[0][2].toInteger(),
                lane: m[0][3],
                read: m[0][4],
                chunk: m[0][5]
            ] : [:]

            tuple(sample_meta + file_meta + [priority: priority], fastq_path)
        }
    ```

=== "Before"

    ```groovy title="main.nf" linenums="7" hl_lines="10"
        .map { row ->
            def sample_meta = [
                id: row.sample_id.toLowerCase(),
                organism: row.organism,
                tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                depth: row.sequencing_depth.toInteger(),
                quality: row.quality_score.toDouble()
            ]
            def priority = sample_meta.quality > 40 ? 'high' : 'normal'
            tuple(sample_meta + [priority: priority], file(row.file_path))
        }
    ```

Há muito a descompactar aqui!

1. `fastq_path.name` obtém o nome do arquivo (sem o caminho)
2. O operador `=~` é para correspondência de padrões regex
3. O padrão `/^(.+)_S(\d+)_L(\d{3})_(R[12])_(\d{3})\.fastq(?:\.gz)?$/` captura os componentes
4. `m ? ... : ...` é um operador ternário que lida com o caso onde o nome do arquivo não corresponde ao padrão
5. `m[0][2]` acessa o segundo grupo de captura (os índices começam em 1 para os grupos)
6. A função `toInteger()` converte a string capturada em um número inteiro
7. `sample_meta + file_meta + [priority: priority]` mescla os três mapas em um

Vamos também modificar o processo `GENERATE_REPORT` para incluir esses novos metadados:

=== "After"

    ```groovy title="modules/generate_report.nf" linenums="11" hl_lines="7-10"
    script:
    """
    echo "Sample ID: ${meta.id}" > ${meta.id}_report.txt
    echo "Organism: ${meta.organism}" >> ${meta.id}_report.txt
    echo "Tissue: ${meta.tissue}" >> ${meta.id}_report.txt
    echo "Sequencing depth: ${meta.depth}" >> ${meta.id}_report.txt
    echo "Quality score: ${meta.quality}" >> ${meta.id}_report.txt
    echo "Sample number: ${meta.sample_num}" >> ${meta.id}_report.txt
    echo "Lane: ${meta.lane}" >> ${meta.id}_report.txt
    echo "Read: ${meta.read}" >> ${meta.id}_report.txt
    echo "Chunk: ${meta.chunk}" >> ${meta.id}_report.txt
    echo "Priority: ${meta.priority}" >> ${meta.id}_report.txt
    echo "File processed: ${reads}" >> ${meta.id}_report.txt
    echo "Process command: fastp --in1 ${reads} --out1 ${meta.id}_trimmed.fastq.gz" >> ${meta.id}_report.txt
    """
    ```

=== "Before"

    ```groovy title="modules/generate_report.nf" linenums="11" hl_lines="7"
    script:
    """
    echo "Sample ID: ${meta.id}" > ${meta.id}_report.txt
    echo "Organism: ${meta.organism}" >> ${meta.id}_report.txt
    echo "Tissue: ${meta.tissue}" >> ${meta.id}_report.txt
    echo "Sequencing depth: ${meta.depth}" >> ${meta.id}_report.txt
    echo "Quality score: ${meta.quality}" >> ${meta.id}_report.txt
    echo "Priority: ${meta.priority}" >> ${meta.id}_report.txt
    echo "File processed: ${reads}" >> ${meta.id}_report.txt
    echo "Process command: fastp --in1 ${reads} --out1 ${meta.id}_trimmed.fastq.gz" >> ${meta.id}_report.txt
    """
    ```

Vamos executar o fluxo de trabalho novamente:

```bash
nextflow run main.nf
```

??? success "Saída do comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [jolly_hawking] DSL2 - revision: cd0a5b0d29

    executor >  local (6)
    [b3/1cb89f] process > FASTP (sample_001)           [100%] 3 of 3 ✔
    [e6/c2f254] process > GENERATE_REPORT (sample_001) [100%] 3 of 3 ✔
    ```

Vamos verificar um dos arquivos de relatório atualizados:

```console
cat work/e6/c2f2542ec23ee6e7f0aa9c66c12e30/sample_001_report.txt
```

Agora você deve ver os metadados adicionais extraídos do nome do arquivo:

```
Sample ID: sample_001
Organism: human
Tissue: liver
Sequencing depth: 30000000
Quality score: 38.5
Sample number: 1
Lane: 001
Read: R1
Chunk: 001
Priority: normal
File processed: SAMPLE_001_S1_L001_R1_001.fastq
Process command: fastp --in1 SAMPLE_001_S1_L001_R1_001.fastq --out1 sample_001_trimmed.fastq.gz
```

### 1.3. Operações de coleção em closures vs manipulação de canais

Um ponto frequentemente confuso no Nextflow são as operações de coleções dentro das closures. Métodos como `collect()` funcionam de maneira diferente dependendo do contexto:

1. Nos **canais Nextflow**: `channel.collect()` é um operador que reúne todos os elementos de um canal em uma única lista
2. Nas **listas e mapas Groovy**: `list.collect {...}` aplica uma função a cada elemento e retorna uma nova lista

Vamos modificar o processo FASTP para ilustrar essa diferença:

```groovy title="modules/fastp.nf" linenums="11" hl_lines="3-5"
script:
"""
# Demonstração de operações de coleção no script
def options = ["--in1", reads, "--out1", "${meta.id}_trimmed.fastq.gz", "--json", "${meta.id}.fastp.json", "--html", "${meta.id}.fastp.html"]
fastp ${options.join(' ')}
"""
```

Este exemplo usa `options.join(' ')` para juntar os elementos da lista em uma única string com espaços entre eles.

No entanto, isso gerará um erro porque estamos tentando executar código Groovy dentro de um bloco de script bash. Vamos alterá-lo para mover a lógica de coleção para fora do bloco de script:

=== "After"

    ```groovy title="modules/fastp.nf" linenums="11" hl_lines="1-3"
    def options = ["--in1", reads, "--out1", "${meta.id}_trimmed.fastq.gz", "--json", "${meta.id}.fastp.json", "--html", "${meta.id}.fastp.html"]
    def cmd = "fastp ${options.join(' ')}"

    script:
    """
    $cmd
    """
    ```

=== "Before"

    ```groovy title="modules/fastp.nf" linenums="11" hl_lines="2-3"
    script:
    """
    fastp --in1 ${reads} --out1 ${meta.id}_trimmed.fastq.gz --json ${meta.id}.fastp.json --html ${meta.id}.fastp.html
    """
    ```

Execute o fluxo de trabalho novamente e deve funcionar! Isso demonstra o uso de scripting Groovy para manipulação de coleções antes de passar o comando para o bloco de script.

```bash
nextflow run main.nf
```

??? success "Saída do comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [compassionate_shaw] DSL2 - revision: 3471dc57d9

    executor >  local (6)
    [32/8e94af] process > FASTP (sample_001)           [100%] 3 of 3 ✔
    [6e/e3e56d] process > GENERATE_REPORT (sample_001) [100%] 3 of 3 ✔
    ```

### Conclusão

Nesta seção, aprendemos sobre as diferenças entre o **fluxo de dados** (operações de canal) e o **scripting** (código dentro de closures e blocos de script). Usamos:

- **Closures** para extrair e transformar metadados dos dados de entrada
- **Operadores ternários** (`condição ? valor_verdadeiro : valor_falso`) para lógica condicional compacta
- **Expressões regulares** com o operador `=~` para extrair componentes dos nomes de arquivos
- **Manipulação de coleções** como `join()` para criar strings de comandos

Essas técnicas são fundamentais para escrever fluxos de trabalho Nextflow que sejam limpos, eficazes e fáceis de manter.

---

## 2. Processamento de strings para lidar com nomes de arquivos e metadados

O processamento de strings é uma tarefa comum em fluxos de trabalho bioinformáticos, especialmente ao extrair metadados de nomes de arquivos ou gerar scripts dinamicamente. Nextflow tem poderosas capacidades de manipulação de strings que são cruciais para fluxos de trabalho robustos.

### 2.1. Expressões regulares para analisar nomes de arquivos

Já usamos expressões regulares para extrair metadados dos nomes dos arquivos FASTQ. Vamos entender melhor como funciona o operador de correspondência de padrões `=~`.

Quando você escreve `x =~ /padrão/`:

1. Cria um objeto `java.util.regex.Matcher`
2. Se avaliado em um contexto booleano, verifica se há **alguma correspondência**
3. Se atribuído a uma variável, você pode acessar os **grupos de captura**

Vamos ver como poderíamos complicar nosso regex de fastq para lidar com mais variações:

=== "After"

    ```groovy title="main.nf" linenums="14" hl_lines="2"
            // Extract metadata from filename using regex
            def m = (fastq_path.name =~ /^(.+)_S(\d+)_L(\d{3})_(R[12])_(\d{3})\.fastq(?:\.gz)?$/)
            def file_meta = m ? [
                sample_num: m[0][2].toInteger(),
                lane: m[0][3],
                read: m[0][4],
                chunk: m[0][5]
            ] : [:]
    ```

=== "Before"

    ```groovy title="main.nf" linenums="14" hl_lines="2"
            // Extract metadata from filename using regex
            def m = (fastq_path.name =~ /^(.+)_S(\d+)_L(\d{3})_(R[12])_(\d{3})\.fastq$/)
            def file_meta = m ? [
                sample_num: m[0][2].toInteger(),
                lane: m[0][3],
                read: m[0][4],
                chunk: m[0][5]
            ] : [:]
    ```

Agora, `/^(.+)_S(\d+)_L(\d{3})_(R[12])_(\d{3})\.fastq(?:\.gz)?$/` corresponderá a `.fastq` e também a `.fastq.gz`.

!!! note "Grupos de captura em regex" - Os parênteses `(...)` definem "grupos de captura" - `m[0]` é a correspondência completa - `m[0][1]`, `m[0][2]`, etc. são os grupos de captura (começando com 1)

Se você tivesse nomes de arquivos com diferentes convenções, poderia usar o operador OR (`|`) em seu regex:

```groovy title="exemplo de regex"
def m = (filename =~ /^(.+)_S(\d+)_L(\d{3})_(R[12])_(\d{3})\.fastq(?:\.gz)?$/ |
                      /^(.+)_(\d{6})_([ACGT]+)_L(\d{3})_(R[12])\.fastq(?:\.gz)?$/)
```

### 2.2. Geração dinâmica de scripts

Outra aplicação poderosa do processamento de strings é a geração dinâmica de scripts bash baseados em metadados ou entradas. Isso é particularmente útil para lógica condicional dentro dos processos.

Vamos modificar o processo `GENERATE_REPORT` para gerar diferentes tipos de relatórios com base na prioridade:

=== "After"

    ```groovy title="modules/generate_report.nf" linenums="10" hl_lines="3-11"
    script:
    // Generate different report content based on priority
    def report_content = meta.priority == 'high' ?
    """
    echo "AMOSTRA DE ALTA PRIORIDADE" > ${meta.id}_report.txt
    echo "===============================================" >> ${meta.id}_report.txt
    """ :
    """
    echo "Amostra padrão" > ${meta.id}_report.txt
    echo "---------------------------------------------" >> ${meta.id}_report.txt
    """

    """
    ${report_content}
    echo "Sample ID: ${meta.id}" >> ${meta.id}_report.txt
    echo "Organism: ${meta.organism}" >> ${meta.id}_report.txt
    echo "Tissue: ${meta.tissue}" >> ${meta.id}_report.txt
    echo "Sequencing depth: ${meta.depth}" >> ${meta.id}_report.txt
    echo "Quality score: ${meta.quality}" >> ${meta.id}_report.txt
    echo "Sample number: ${meta.sample_num}" >> ${meta.id}_report.txt
    echo "Lane: ${meta.lane}" >> ${meta.id}_report.txt
    echo "Read: ${meta.read}" >> ${meta.id}_report.txt
    echo "Chunk: ${meta.chunk}" >> ${meta.id}_report.txt
    echo "Priority: ${meta.priority}" >> ${meta.id}_report.txt
    echo "File processed: ${reads}" >> ${meta.id}_report.txt
    echo "Process command: fastp --in1 ${reads} --out1 ${meta.id}_trimmed.fastq.gz" >> ${meta.id}_report.txt
    """
    ```

=== "Before"

    ```groovy title="modules/generate_report.nf" linenums="11" hl_lines="2-13"
    script:
    """
    echo "Sample ID: ${meta.id}" > ${meta.id}_report.txt
    echo "Organism: ${meta.organism}" >> ${meta.id}_report.txt
    echo "Tissue: ${meta.tissue}" >> ${meta.id}_report.txt
    echo "Sequencing depth: ${meta.depth}" >> ${meta.id}_report.txt
    echo "Quality score: ${meta.quality}" >> ${meta.id}_report.txt
    echo "Sample number: ${meta.sample_num}" >> ${meta.id}_report.txt
    echo "Lane: ${meta.lane}" >> ${meta.id}_report.txt
    echo "Read: ${meta.read}" >> ${meta.id}_report.txt
    echo "Chunk: ${meta.chunk}" >> ${meta.id}_report.txt
    echo "Priority: ${meta.priority}" >> ${meta.id}_report.txt
    echo "File processed: ${reads}" >> ${meta.id}_report.txt
    echo "Process command: fastp --in1 ${reads} --out1 ${meta.id}_trimmed.fastq.gz" >> ${meta.id}_report.txt
    """
    ```

Vamos executar o fluxo de trabalho para ver o resultado:

```bash
nextflow run main.nf
```

??? success "Saída do comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [condescending_northcutt] DSL2 - revision: 1a3c16a96f

    executor >  local (6)
    [95/e633a2] process > FASTP (sample_001)           [100%] 3 of 3 ✔
    [a8/e4c214] process > GENERATE_REPORT (sample_001) [100%] 3 of 3 ✔
    ```

Vamos examinar os arquivos de relatório gerados - procurando especificamente o de alta prioridade:

```console
find work -name "sample_003_report.txt" -exec cat {} \;
```

Você deve ver o cabeçalho especial:

```
AMOSTRA DE ALTA PRIORIDADE
===============================================
Sample ID: sample_003
Organism: human
...
```

### 2.3. Interpolação de variáveis: quando Nextflow interpreta vs quando bash interpreta

Um ponto sutil, mas crucial para entender, é quando ocorre a interpolação de variáveis:

1. `${var}` - Interpolado pelo Nextflow durante a compilação do script
2. `\${var}` - Escapado, passado literalmente para bash como `${var}` (para variáveis de ambiente bash)
3. `\$(cmd)` - Substituição de comandos shell (escapado, executado por bash em tempo de execução)

Podemos ver isso em ação atualizando o processo `GENERATE_REPORT` para usar uma variável de ambiente:

=== "After"

    ```groovy title="modules/generate_report.nf" linenums="11" hl_lines="15"
    script:
    // Generate different report content based on priority
    def report_content = meta.priority == 'high' ?
    """
    echo "AMOSTRA DE ALTA PRIORIDADE" > ${meta.id}_report.txt
    echo "===============================================" >> ${meta.id}_report.txt
    """ :
    """
    echo "Amostra padrão" > ${meta.id}_report.txt
    echo "---------------------------------------------" >> ${meta.id}_report.txt
    """

    """
    ${report_content}
    echo "Sample ID: ${meta.id}" >> ${meta.id}_report.txt
    echo "Processed by: \${USER}" >> ${meta.id}_report.txt
    echo "Organism: ${meta.organism}" >> ${meta.id}_report.txt
    echo "Tissue: ${meta.tissue}" >> ${meta.id}_report.txt
    echo "Sequencing depth: ${meta.depth}" >> ${meta.id}_report.txt
    echo "Quality score: ${meta.quality}" >> ${meta.id}_report.txt
    echo "Sample number: ${meta.sample_num}" >> ${meta.id}_report.txt
    echo "Lane: ${meta.lane}" >> ${meta.id}_report.txt
    echo "Read: ${meta.read}" >> ${meta.id}_report.txt
    echo "Chunk: ${meta.chunk}" >> ${meta.id}_report.txt
    echo "Priority: ${meta.priority}" >> ${meta.id}_report.txt
    echo "File processed: ${reads}" >> ${meta.id}_report.txt
    echo "Process command: fastp --in1 ${reads} --out1 ${meta.id}_trimmed.fastq.gz" >> ${meta.id}_report.txt
    """
    ```

=== "Before"

    ```groovy title="modules/generate_report.nf" linenums="11"
    script:
    // Generate different report content based on priority
    def report_content = meta.priority == 'high' ?
    """
    echo "AMOSTRA DE ALTA PRIORIDADE" > ${meta.id}_report.txt
    echo "===============================================" >> ${meta.id}_report.txt
    """ :
    """
    echo "Amostra padrão" > ${meta.id}_report.txt
    echo "---------------------------------------------" >> ${meta.id}_report.txt
    """

    """
    ${report_content}
    echo "Sample ID: ${meta.id}" >> ${meta.id}_report.txt
    echo "Organism: ${meta.organism}" >> ${meta.id}_report.txt
    echo "Tissue: ${meta.tissue}" >> ${meta.id}_report.txt
    echo "Sequencing depth: ${meta.depth}" >> ${meta.id}_report.txt
    echo "Quality score: ${meta.quality}" >> ${meta.id}_report.txt
    echo "Sample number: ${meta.sample_num}" >> ${meta.id}_report.txt
    echo "Lane: ${meta.lane}" >> ${meta.id}_report.txt
    echo "Read: ${meta.read}" >> ${meta.id}_report.txt
    echo "Chunk: ${meta.chunk}" >> ${meta.id}_report.txt
    echo "Priority: ${meta.priority}" >> ${meta.id}_report.txt
    echo "File processed: ${reads}" >> ${meta.id}_report.txt
    echo "Process command: fastp --in1 ${reads} --out1 ${meta.id}_trimmed.fastq.gz" >> ${meta.id}_report.txt
    """
    ```

Se executarmos o fluxo de trabalho como está, ele falhará:

```bash
nextflow run main.nf
```

??? failure "Saída do comando"

    ```console
    Error modules/generate_report.nf:15:27: `USER` is not defined
    │  15 |     echo "Processed by: ${USER}" >> ${meta.id}_report.txt
    ╰     |                           ^^^^

    ERROR ~ Script compilation failed
    ```

O problema é que Nextflow tenta interpretar `${USER}` como uma variável Nextflow, mas não existe nenhuma variável chamada `USER`. Precisamos escapar o `$` para que ele seja passado para o bash, que tem uma variável de ambiente `USER`:

```groovy title="modules/generate_report.nf" linenums="15" hl_lines="1"
    echo "Processed by: \${USER}" >> ${meta.id}_report.txt
```

Agora funciona! A barra invertida (`\`) diz ao Nextflow "não interprete isso, passe para o Bash."

### Conclusão

Nesta seção, você aprendeu técnicas de **processamento de strings**:

- **Expressões regulares para análise de arquivos**: Usando o operador `=~` e padrões regex (`~/padrão/`) para extrair metadados de convenções complexas de nomenclatura de arquivos
- **Geração dinâmica de scripts**: Usando lógica condicional (if/else, operadores ternários) para gerar diferentes strings de script com base nas características da entrada
- **Interpolação de variáveis**: Entendendo quando o Nextflow interpreta strings vs quando o shell interpreta
  - `${var}` - Variáveis Nextflow (interpoladas por Nextflow em tempo de compilação do workflow)
  - `\${var}` - Variáveis de ambiente shell (escapadas, passadas para bash em tempo de execução)
  - `\$(cmd)` - Substituição de comandos shell (escapada, executada por bash em tempo de execução)

Estes padrões de processamento e geração de strings são essenciais para lidar com os diversos formatos de arquivo e convenções de nomenclatura que você encontrará em fluxos de trabalho bioinformáticos do mundo real.

---

## 3. Criando Funções Reutilizáveis

Lógica complexa de fluxo de trabalho inline em operadores de canal ou definições de processo reduz a legibilidade e a manutenibilidade. **Funções** permitem que você extraia essa lógica em componentes nomeados e reutilizáveis.

Nossa operação map cresceu longa e complexa. Vamos extraí-la em uma função reutilizável usando a palavra-chave `def`.

Para ilustrar como isso fica com nosso fluxo de trabalho existente, faça a modificação abaixo, usando `def` para definir uma função reutilizável chamada `separateMetadata`:

=== "After"

    ```groovy title="main.nf" linenums="1" hl_lines="4-24 29"
    include { FASTP } from './modules/fastp.nf'
    include { GENERATE_REPORT } from './modules/generate_report.nf'

    def separateMetadata(row) {
        def sample_meta = [
            id: row.sample_id.toLowerCase(),
            organism: row.organism,
            tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
            depth: row.sequencing_depth.toInteger(),
            quality: row.quality_score.toDouble()
        ]
        def fastq_path = file(row.file_path)

        def m = (fastq_path.name =~ /^(.+)_S(\d+)_L(\d{3})_(R[12])_(\d{3})\.fastq(?:\.gz)?$/)
        def file_meta = m ? [
            sample_num: m[0][2].toInteger(),
            lane: m[0][3],
            read: m[0][4],
            chunk: m[0][5]
        ] : [:]

        def priority = sample_meta.quality > 40 ? 'high' : 'normal'
        return tuple(sample_meta + file_meta + [priority: priority], fastq_path)
    }

    workflow {
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map{ row -> separateMetadata(row) }

        ch_fastp = FASTP(ch_samples)
        GENERATE_REPORT(ch_samples)
    }
    ```

=== "Before"

    ```groovy title="main.nf" linenums="1" hl_lines="7-27"
    include { FASTP } from './modules/fastp.nf'
    include { GENERATE_REPORT } from './modules/generate_report.nf'

    workflow {
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row ->
                def sample_meta = [
                    id: row.sample_id.toLowerCase(),
                    organism: row.organism,
                    tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                    depth: row.sequencing_depth.toInteger(),
                    quality: row.quality_score.toDouble()
                ]
                def fastq_path = file(row.file_path)

                def m = (fastq_path.name =~ /^(.+)_S(\d+)_L(\d{3})_(R[12])_(\d{3})\.fastq(?:\.gz)?$/)
                def file_meta = m ? [
                    sample_num: m[0][2].toInteger(),
                    lane: m[0][3],
                    read: m[0][4],
                    chunk: m[0][5]
                ] : [:]

                def priority = sample_meta.quality > 40 ? 'high' : 'normal'
                return tuple(sample_meta + file_meta + [priority: priority], fastq_path)
            }

        ch_fastp = FASTP(ch_samples)
        GENERATE_REPORT(ch_samples)
    }
    ```

Ao extrair esta lógica em uma função, reduzimos a lógica real do fluxo de trabalho para algo muito mais limpo:

```groovy title="workflow mínimo"
    ch_samples = channel.fromPath("./data/samples.csv")
        .splitCsv(header: true)
        .map{ row -> separateMetadata(row) }

    ch_fastp = FASTP(ch_samples)
    GENERATE_REPORT(ch_samples)
```

Isso torna a lógica do fluxo de trabalho muito mais fácil de ler e entender de relance. A função `separateMetadata` encapsula toda a lógica complexa para analisar e enriquecer metadados, tornando-a reutilizável e testável.

Execute o fluxo de trabalho para garantir que ainda funcione:

```bash
nextflow run main.nf
```

??? success "Saída do comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [admiring_panini] DSL2 - revision: 8cc832e32f

    executor >  local (6)
    [8c/2e3f91] process > FASTP (3)           [100%] 3 of 3 ✔
    [7a/1b4c92] process > GENERATE_REPORT (3) [100%] 3 of 3 ✔
    ```

A saída deve mostrar que ambos os processos completam com sucesso. O fluxo de trabalho agora está muito mais limpo e fácil de manter, com toda a lógica complexa de processamento de metadados encapsulada na função `separateMetadata`.

### Conclusão

Nesta seção, você aprendeu sobre **criação de funções**:

- **Definindo funções com `def`**: A palavra-chave para criar funções nomeadas (como `def` em Python ou `function` em JavaScript)
- **Escopo de funções**: Funções definidas no nível do script são acessíveis em todo o seu fluxo de trabalho Nextflow
- **Valores de retorno**: Funções automaticamente retornam a última expressão, ou use `return` explicitamente
- **Código mais limpo**: Extrair lógica complexa em funções é uma prática fundamental de engenharia de software em qualquer linguagem

A seguir, exploraremos como usar closures em diretivas de processo para alocação dinâmica de recursos.

---

## 4. Diretivas de Recursos Dinâmicas com Closures

Até agora, usamos scripting no bloco `script` de processos. Mas as **closures** (introduzidas na Seção 1.1) também são incrivelmente úteis nas diretivas de processo, especialmente para alocação dinâmica de recursos. Vamos adicionar diretivas de recursos ao nosso processo FASTP que se adaptam com base nas características da amostra.

### 4.1. Alocação de recursos específica por amostra

Atualmente, nosso processo FASTP usa recursos padrão. Vamos torná-lo mais inteligente, alocando mais CPUs para amostras de alta profundidade. Edite `modules/fastp.nf` para incluir uma diretiva `cpus` dinâmica e uma diretiva `memory` estática:

=== "After"

    ```groovy title="modules/fastp.nf" linenums="1" hl_lines="4-5"
    process FASTP {
        container 'community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690'

        cpus { meta.depth > 40000000 ? 2 : 1 }
        memory 2.GB

        input:
        tuple val(meta), path(reads)
    ```

=== "Before"

    ```groovy title="modules/fastp.nf" linenums="1"
    process FASTP {
        container 'community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690'

        input:
        tuple val(meta), path(reads)
    ```

A closure `{ meta.depth > 40000000 ? 2 : 1 }` usa o **operador ternário** (abordado na Seção 1.1) e é avaliada para cada tarefa, permitindo alocação de recursos por amostra. Amostras de alta profundidade (>40M leituras) recebem 2 CPUs, enquanto outras recebem 1 CPU.

!!! note "Acessando Variáveis de Entrada nas Diretivas"

    A closure pode acessar quaisquer variáveis de entrada (como `meta` aqui) porque Nextflow avalia essas closures no contexto de cada execução de tarefa.

Execute o fluxo de trabalho novamente com a opção `-ansi-log false` para tornar mais fácil ver os hashes de tarefas.

```bash
nextflow run main.nf -ansi-log false
```

??? success "Saída do comando"

    ```console
    N E X T F L O W  ~  version 25.10.2
    Launching `main.nf` [fervent_albattani] DSL2 - revision: fa8f249759
    [bd/ff3d41] Submitted process > FASTP (2)
    [a4/a3aab2] Submitted process > FASTP (1)
    [48/6db0c9] Submitted process > FASTP (3)
    [ec/83439d] Submitted process > GENERATE_REPORT (3)
    [bd/15d7cc] Submitted process > GENERATE_REPORT (2)
    [42/699357] Submitted process > GENERATE_REPORT (1)
    ```

Você pode verificar o comando exato do `docker` que foi executado para ver a alocação de CPU para qualquer tarefa específica:

```console title="Verificar comando docker"
cat work/48/6db0c9e9d8aa65e4bb4936cd3bd59e/.command.run | grep "docker run"
```

Você deve ver algo como:

```bash title="comando docker"
    docker run -i --cpu-shares 4096 --memory 2048m -e "NXF_TASK_WORKDIR" -v /workspaces/training/side-quests/essential_scripting_patterns:/workspaces/training/side-quests/essential_scripting_patterns -w "$NXF_TASK_WORKDIR" --name $NXF_BOXID community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690 /bin/bash -ue /workspaces/training/side-quests/essential_scripting_patterns/work/48/6db0c9e9d8aa65e4bb4936cd3bd59e/.command.sh
```

Neste exemplo, escolhemos um caso que solicitou 2 CPUs (`--cpu-shares 2048`), porque era uma amostra de alta profundidade, mas você deve ver diferentes alocações de CPU dependendo da profundidade da amostra. Tente isso também para as outras tarefas.

### 4.2. Estratégias de repetição

Outro padrão poderoso é usar `task.attempt` para estratégias de repetição. Para mostrar por que isso é útil, vamos começar reduzindo a alocação de memória para FASTP para menos do que ele precisa. Altere a diretiva `memory` em `modules/fastp.nf` para `1.GB`:

=== "After"

    ```groovy title="modules/fastp.nf" linenums="1" hl_lines="5"
    process FASTP {
        container 'community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690'

        cpus { meta.depth > 40000000 ? 4 : 2 }
        memory 1.GB

        input:
        tuple val(meta), path(reads)
    ```

=== "Before"

    ```groovy title="modules/fastp.nf" linenums="1" hl_lines="5"
    process FASTP {
        container 'community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690'

        cpus { meta.depth > 40000000 ? 4 : 2 }
        memory 2.GB

        input:
        tuple val(meta), path(reads)
    ```

... e execute o fluxo de trabalho novamente:

```bash
nextflow run main.nf
```

??? failure "Saída do comando"

    ```console hl_lines="2 11"
    Command exit status:
      137

    Command output:
      (empty)

    Command error:
      Detecting adapter sequence for read1...
      No adapter detected for read1

      .command.sh: line 7:   101 Killed                  fastp --in1 SAMPLE_002_S2_L001_R1_001.fastq --out1 sample_002_trimmed.fastq.gz --json sample_002.fastp.json --html sample_002.fastp.html --thread 2
    ```

Isso indica que o processo foi encerrado por exceder os limites de memória.

Este é um cenário muito comum em fluxos de trabalho do mundo real - às vezes você simplesmente não sabe quanta memória uma tarefa vai precisar até executá-la.

Para tornar nosso fluxo de trabalho mais robusto, podemos implementar uma estratégia de repetição que aumenta a alocação de memória em cada tentativa, novamente usando uma closure Groovy. Modifique a diretiva `memory` para multiplicar a memória base por `task.attempt`, e adicione as diretivas `errorStrategy 'retry'` e `maxRetries 2`:

=== "After"

    ```groovy title="modules/fastp.nf" linenums="1" hl_lines="5-7"
    process FASTP {
        container 'community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690'

        cpus { meta.depth > 40000000 ? 4 : 2 }
        memory { 1.GB * task.attempt }
        errorStrategy 'retry'
        maxRetries 2

        input:
        tuple val(meta), path(reads)
    ```

=== "Before"

    ```groovy title="modules/fastp.nf" linenums="1" hl_lines="5"
    process FASTP {
        container 'community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690'

        cpus { meta.depth > 40000000 ? 4 : 2 }
        memory 2.GB

        input:
        tuple val(meta), path(reads)
    ```

Agora, se o processo falhar devido a memória insuficiente, Nextflow tentará novamente com mais memória:

- Primeira tentativa: 1 GB (task.attempt = 1)
- Segunda tentativa: 2.GB (task.attempt = 2)

... e assim por diante, até o limite `maxRetries`.

### Conclusão

Diretivas dinâmicas com closures permitem que você:

- Aloque recursos baseados nas características de entrada
- Implemente estratégias automáticas de repetição com recursos crescentes
- Combine múltiplos fatores (metadados, número de tentativa, prioridades)
- Use lógica condicional para cálculos complexos de recursos

Isso torna seus fluxos de trabalho tanto mais eficientes (sem superalocação) quanto mais robustos (repetição automática com mais recursos).

---

## 5. Lógica Condicional e Controle de Processo

Anteriormente, usamos `.map()` com scripting para transformar dados de canal. Agora, usaremos lógica condicional para controlar quais processos executar com base nos dados — essencial para fluxos de trabalho flexíveis adaptando-se a diferentes tipos de amostra.

Os [operadores de fluxo de dados](https://www.nextflow.io/docs/latest/reference/operator.html) do Nextflow aceitam closures avaliadas em tempo de execução, permitindo lógica condicional para direcionar decisões de fluxo de trabalho com base no conteúdo do canal.

### 5.1. Roteamento com `.branch()`

Por exemplo, vamos imaginar que nossas amostras de sequenciamento precisam ser recortadas com FASTP apenas se forem amostras humanas com uma cobertura acima de um certo limiar. Amostras de camundongos ou amostras de baixa cobertura devem ser executadas com Trimgalore (este é um exemplo inventado, mas ilustra o ponto).

Fornecemos um processo simples de Trimgalore em `modules/trimgalore.nf`, dê uma olhada se quiser, mas os detalhes não são importantes para este exercício. O ponto-chave é que queremos rotear amostras com base em seus metadados.

Inclua o novo módulo de `modules/trimgalore.nf`:

=== "After"

    ```groovy title="main.nf" linenums="1" hl_lines="2"
    include { FASTP } from './modules/fastp.nf'
    include { TRIMGALORE } from './modules/trimgalore.nf'
    ```

=== "Before"

    ```groovy title="main.nf" linenums="1"
    include { FASTP } from './modules/fastp.nf'
    ```

... e então modifique seu fluxo de trabalho `main.nf` para ramificar amostras com base em seus metadados e roteá-las através do processo de trimming apropriado, assim:

=== "After"

    ```groovy title="main.nf" linenums="28" hl_lines="5-12"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row -> separateMetadata(row) }

        trim_branches = ch_samples
            .branch { meta, reads ->
                fastp: meta.organism == 'human' && meta.depth >= 30000000
                trimgalore: true
            }

        ch_fastp = FASTP(trim_branches.fastp)
        ch_trimgalore = TRIMGALORE(trim_branches.trimgalore)
        GENERATE_REPORT(ch_samples)
    ```

=== "Before"

    ```groovy title="main.nf" linenums="28" hl_lines="5"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row -> separateMetadata(row) }

        ch_fastp = FASTP(ch_samples)
        GENERATE_REPORT(ch_samples)
    ```

Execute este fluxo de trabalho modificado:

```bash
nextflow run main.nf
```

??? success "Saída do comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [adoring_galileo] DSL2 - revision: c9e83aaef1

    executor >  local (6)
    [1d/0747ac] process > FASTP (2)           [100%] 2 of 2 ✔
    [cc/c44caf] process > TRIMGALORE (1)      [100%] 1 of 1 ✔
    [34/bd5a9f] process > GENERATE_REPORT (1) [100%] 3 of 3 ✔
    ```

Aqui, usamos expressões condicionais pequenas, mas poderosas, dentro do operador `.branch{}` para rotear amostras com base em seus metadados. Amostras humanas com alta cobertura passam por `FASTP`, enquanto todas as outras amostras passam por `TRIMGALORE`.

### 5.2. Usando `.filter()` com Truthiness

Outro padrão poderoso para controlar a execução do fluxo de trabalho é o operador `.filter()`, que usa uma closure para determinar quais itens devem continuar pelo pipeline. Dentro da closure de filtro, você escreverá **expressões booleanas** que decidem quais itens passam.

Nextflow (como muitas linguagens dinâmicas) tem um conceito de **"truthiness"** que determina quais valores são avaliados como `true` ou `false` em contextos booleanos:

- **Truthy**: Valores não nulos, strings não vazias, números não zero, coleções não vazias
- **Falsy**: `null`, strings vazias `""`, zero `0`, coleções vazias `[]` ou `[:]`, `false`

Isso significa que `meta.id` sozinho (sem um explícito `!= null`) verifica se o ID existe e não está vazio. Vamos usar isso para filtrar amostras que não atendem aos nossos requisitos de qualidade.

Adicione o seguinte antes da operação branch:

=== "After"

    ```groovy title="main.nf" linenums="28" hl_lines="5-11"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row -> separateMetadata(row) }

        // Filtrar amostras inválidas ou de baixa qualidade
        ch_valid_samples = ch_samples
            .filter { meta, reads ->
                meta.id && meta.organism && meta.depth >= 25000000
            }

        trim_branches = ch_valid_samples
            .branch { meta, reads ->
                fastp: meta.organism == 'human' && meta.depth >= 30000000
                trimgalore: true
            }
    ```

=== "Before"

    ```groovy title="main.nf" linenums="28" hl_lines="5"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row -> separateMetadata(row) }

        trim_branches = ch_samples
            .branch { meta, reads ->
                fastp: meta.organism == 'human' && meta.depth >= 30000000
                trimgalore: true
            }
    ```

Execute o fluxo de trabalho novamente:

```bash
nextflow run main.nf
```

??? success "Saída do comando"

    ```console
    N E X T F L O W  ~  version 25.10.2
    Launching `main.nf` [lonely_williams] DSL2 - revision: d0b3f121ec
    [94/b48eac] Submitted process > FASTP (2)
    [2c/d2b28f] Submitted process > GENERATE_REPORT (2)
    [65/2e3be4] Submitted process > GENERATE_REPORT (1)
    [94/b48eac] NOTE: Process `FASTP (2)` terminated with an error exit status (137) -- Execution is retried (1)
    [3e/0d8664] Submitted process > TRIMGALORE (1)
    [6a/9137b0] Submitted process > FASTP (1)
    [6a/9137b0] NOTE: Process `FASTP (1)` terminated with an error exit status (137) -- Execution is retried (1)
    [83/577ac0] Submitted process > GENERATE_REPORT (3)
    [a2/5117de] Re-submitted process > FASTP (1)
    [1f/a1a4ca] Re-submitted process > FASTP (2)
    ```

Como escolhemos um filtro que exclui algumas amostras, menos tarefas foram executadas.

A expressão de filtro `meta.id && meta.organism && meta.depth >= 25000000` combina truthiness com comparações explícitas:

- `meta.id && meta.organism` verifica se ambos os campos existem e não estão vazios (usando truthiness)
- `meta.depth >= 25000000` garante profundidade de sequenciamento suficiente com uma comparação explícita

!!! note "Truthiness na Prática"

    A expressão `meta.id && meta.organism` é mais concisa do que escrever:
    ```groovy
    meta.id != null && meta.id != '' && meta.organism != null && meta.organism != ''
    ```

    Isso torna a lógica de filtragem muito mais limpa e fácil de ler.

### Conclusão

Nesta seção, você aprendeu a usar lógica condicional para controlar a execução do fluxo de trabalho usando as interfaces de closure de operadores Nextflow como `.branch{}` e `.filter{}`, aproveitando truthiness para escrever expressões condicionais concisas.

Nosso pipeline agora roteia inteligentemente amostras através dos processos apropriados, mas fluxos de trabalho de produção precisam lidar graciosamente com dados inválidos. Vamos tornar nosso fluxo de trabalho robusto contra valores ausentes ou nulos.

---

## 6. Operadores de Navegação Segura e Elvis

Nossa função `separateMetadata` atualmente assume que todos os campos CSV estão presentes e válidos. Mas o que acontece com dados incompletos? Vamos descobrir.

### 6.1. O Problema: Acessando Propriedades que Não Existem

Digamos que queremos adicionar suporte para informações opcionais de execução de sequenciamento. Em alguns laboratórios, as amostras podem ter um campo adicional para o ID de execução de sequenciamento ou número de lote, mas nosso CSV atual não tem essa coluna. Vamos tentar acessá-la mesmo assim.

Modifique a função `separateMetadata` para incluir um campo run_id:

=== "After"

    ```groovy title="main.nf" linenums="5" hl_lines="9"
    def separateMetadata(row) {
        def sample_meta = [
            id: row.sample_id.toLowerCase(),
            organism: row.organism,
            tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
            depth: row.sequencing_depth.toInteger(),
            quality: row.quality_score.toDouble()
        ]
        def run_id = row.run_id.toUpperCase()
    ```

=== "Before"

    ```groovy title="main.nf" linenums="5"
    def separateMetadata(row) {
        def sample_meta = [
            id: row.sample_id.toLowerCase(),
            organism: row.organism,
            tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
            depth: row.sequencing_depth.toInteger(),
            quality: row.quality_score.toDouble()
        ]
    ```

Agora execute o fluxo de trabalho:

```bash
nextflow run main.nf
```

??? failure "Saída do comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [trusting_torvalds] DSL2 - revision: b56fbfbce2

    ERROR ~ Cannot invoke method toUpperCase() on null object

    -- Check script 'main.nf' at line: 13 or see '.nextflow.log' file for more details
    ```

Isto falha com uma NullPointerException.

O problema é que `row.run_id` retorna `null` porque a coluna `run_id` não existe em nosso CSV. Quando tentamos chamar `.toUpperCase()` em `null`, ocorre um erro. É aqui que o operador de navegação segura salva o dia.

### 6.2. Operador de Navegação Segura (`?.`)

O operador de navegação segura (`?.`) retorna `null` em vez de lançar uma exceção quando chamado em um valor `null`. Se o objeto antes de `?.` for `null`, toda a expressão é avaliada como `null` sem executar o método.

Atualize a função para usar a navegação segura:

=== "After"

    ```groovy title="main.nf" linenums="4" hl_lines="9"
    def separateMetadata(row) {
        def sample_meta = [
            id: row.sample_id.toLowerCase(),
            organism: row.organism,
            tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
            depth: row.sequencing_depth.toInteger(),
            quality: row.quality_score.toDouble()
        ]
        def run_id = row.run_id?.toUpperCase()
    ```

=== "Before"

    ```groovy title="main.nf" linenums="4" hl_lines="9"
    def separateMetadata(row) {
        def sample_meta = [
            id: row.sample_id.toLowerCase(),
            organism: row.organism,
            tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
            depth: row.sequencing_depth.toInteger(),
            quality: row.quality_score.toDouble()
        ]
        def run_id = row.run_id.toUpperCase()
    ```

Execute novamente:

```bash
nextflow run main.nf
```

??? success "Saída do comando"

    ```console
    <!-- TODO: output -->
    ```

Sem erro! O fluxo de trabalho agora lida com o campo ausente de forma elegante. Quando `row.run_id` é `null`, o operador `?.` impede a chamada `.toUpperCase()`, e `run_id` se torna `null` em vez de causar uma exceção.

### 6.3. Operador Elvis (`?:`) para Valores Padrão

O operador Elvis (`?:`) fornece valores padrão quando o lado esquerdo é "falsy" (como explicado anteriormente). Ele recebeu esse nome devido a Elvis Presley porque `?:` se parece com seu famoso cabelo e olhos quando visto de lado!

Agora que estamos usando navegação segura, `run_id` será `null` para amostras sem esse campo. Vamos usar o operador Elvis para fornecer um valor padrão e adicioná-lo ao nosso mapa `sample_meta`:

=== "After"

    ```groovy title="main.nf" linenums="5" hl_lines="9-10"
    def separateMetadata(row) {
        def sample_meta = [
            id: row.sample_id.toLowerCase(),
            organism: row.organism,
            tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
            depth: row.sequencing_depth.toInteger(),
            quality: row.quality_score.toDouble()
        ]
        def run_id = row.run_id?.toUpperCase() ?: 'UNSPECIFIED'
        sample_meta.run = run_id
    ```

=== "Before"

    ```groovy title="main.nf" linenums="5" hl_lines="9"
    def separateMetadata(row) {
        def sample_meta = [
            id: row.sample_id.toLowerCase(),
            organism: row.organism,
            tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
            depth: row.sequencing_depth.toInteger(),
            quality: row.quality_score.toDouble()
        ]
        def run_id = row.run_id?.toUpperCase()
    ```

Adicione também um operador `view()` no fluxo de trabalho para ver os resultados:

=== "After"

    ```groovy title="main.nf" linenums="30" hl_lines="4"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map{ row -> separateMetadata(row) }
            .view()
    ```

=== "Before"

    ```groovy title="main.nf" linenums="30"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map{ row -> separateMetadata(row) }
    ```

E execute o fluxo de trabalho:

```bash
nextflow run main.nf
```

??? success "Saída do comando"

    ```console
    [[id:sample_001, organism:human, tissue:liver, depth:30000000, quality:38.5, run:UNSPECIFIED, sample_num:1, lane:001, read:R1, chunk:001, priority:normal], /workspaces/training/side-quests/essential_scripting_patterns/data/sequences/SAMPLE_001_S1_L001_R1_001.fastq]
    [[id:sample_002, organism:mouse, tissue:brain, depth:25000000, quality:35.2, run:UNSPECIFIED, sample_num:2, lane:001, read:R1, chunk:001, priority:normal], /workspaces/training/side-quests/essential_scripting_patterns/data/sequences/SAMPLE_002_S2_L001_R1_001.fastq]
    [[id:sample_003, organism:human, tissue:kidney, depth:45000000, quality:42.1, run:UNSPECIFIED, sample_num:3, lane:001, read:R1, chunk:001, priority:high], /workspaces/training/side-quests/essential_scripting_patterns/data/sequences/SAMPLE_003_S3_L001_R1_001.fastq]
    ```

Perfeito! Agora todas as amostras têm um campo `run` com o ID de execução real (em maiúsculas) ou o valor padrão 'UNSPECIFIED'. A combinação de `?.` e `?:` fornece tanto segurança (sem falhas) quanto valores padrão sensatos.

Remova o operador `.view()` agora que confirmamos que funciona.

!!! tip "Combinando Navegação Segura e Elvis"

    O padrão `value?.method() ?: 'default'` é comum em fluxos de trabalho de produção:

    - `value?.method()` - Chama o método com segurança, retorna `null` se `value` for `null`
    - `?: 'default'` - Fornece fallback se o resultado for `null`

    Este padrão lida com dados ausentes/incompletos com elegância.

Use esses operadores consistentemente em funções, closures de operadores (`.map{}`, `.filter{}`), scripts de processo e arquivos de configuração. Eles evitam falhas ao lidar com dados do mundo real.

### Conclusão

- **Navegação segura (`?.`)**: Evita falhas em valores nulos - retorna null em vez de lançar exceção
- **Operador Elvis (`?:`)**: Fornece valores padrão - `value ?: 'default'`
- **Combinando**: `value?.method() ?: 'default'` é o padrão comum

Esses operadores tornam os fluxos de trabalho resilientes a dados incompletos - essenciais para trabalho no mundo real.

---

## 7. Validação com `error()` e `log.warn`

Às vezes você precisa parar o fluxo de trabalho imediatamente se os parâmetros de entrada forem inválidos. No Nextflow, você pode usar funções integradas como `error()` e `log.warn`, bem como construções de programação padrão como declarações `if` e lógica booleana, para implementar lógica de validação. Vamos adicionar validação ao nosso fluxo de trabalho.

Crie uma função de validação antes do bloco do fluxo de trabalho, chame-a a partir do fluxo de trabalho e altere a criação do canal para usar um parâmetro para o caminho do arquivo CSV. Se o parâmetro estiver ausente ou o arquivo não existir, chame `error()` para interromper a execução com uma mensagem clara.

=== "After"

    ```groovy title="main.nf" linenums="1" hl_lines="5-20 23-24"
    include { FASTP } from './modules/fastp.nf'
    include { TRIMGALORE } from './modules/trimgalore.nf'
    include { GENERATE_REPORT } from './modules/generate_report.nf'

    def validateInputs() {
        // Check input parameter is provided
        if (!params.input) {
            error("Input CSV file path not provided. Please specify --input <file.csv>")
        }

        // Check CSV file exists
        if (!file(params.input).exists()) {
            error("Input CSV file not found: ${params.input}")
        }
    }
    ...
    workflow {
        validateInputs()
        ch_samples = channel.fromPath(params.input)
    ```

=== "Before"

    ```groovy title="main.nf" linenums="1"
    include { FASTP } from './modules/fastp.nf'
    include { TRIMGALORE } from './modules/trimgalore.nf'
    include { GENERATE_REPORT } from './modules/generate_report.nf'

    ...
    workflow {
        ch_samples = channel.fromPath("./data/samples.csv")
    ```

Agora tente executar sem o arquivo CSV:

```bash
nextflow run main.nf
```

??? failure "Saída do comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [confident_coulomb] DSL2 - revision: 07059399ed

    WARN: Access to undefined parameter `input` -- Initialise it to a default value eg. `params.input = some_value`
    Input CSV file path not provided. Please specify --input <file.csv>
    ```

O fluxo de trabalho para imediatamente com uma mensagem de erro clara em vez de falhar misteriosamente depois.

Agora execute-o com um arquivo inexistente:

```bash
nextflow run main.nf --input ./data/nonexistent.csv
```

??? failure "Saída do comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [cranky_gates] DSL2 - revision: 26839ae3eb

    Input CSV file not found: ./data/nonexistent.csv
    ```

Finalmente, execute-o com o arquivo correto:

```bash
nextflow run main.nf --input ./data/samples.csv
```

??? success "Saída do comando"

    ```console
    <!-- TODO: output -->
    ```

Desta vez executa com sucesso.

Você também pode adicionar validação dentro da função `separateMetadata`. Vamos usar o `log.warn` não fatal para emitir avisos para amostras com baixa profundidade de sequenciamento, mas ainda permitir que o fluxo de trabalho continue:

=== "After"

    ```groovy title="main.nf" linenums="1" hl_lines="3-6"
        def priority = sample_meta.quality > 40 ? 'high' : 'normal'

        // Validate data makes sense
        if (sample_meta.depth < 30000000) {
            log.warn "Low sequencing depth for ${sample_meta.id}: ${sample_meta.depth}"
        }

        return tuple(sample_meta + file_meta + [priority: priority], fastq_path)
    }
    ```

=== "Before"

    ```groovy title="main.nf" linenums="1"
        def priority = sample_meta.quality > 40 ? 'high' : 'normal'

        return tuple(sample_meta + file_meta + [priority: priority], fastq_path)
    }
    ```

Execute o fluxo de trabalho novamente com o CSV original:

```bash
nextflow run main.nf --input ./data/samples.csv
```

??? warning "Saída do comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [awesome_goldwasser] DSL2 - revision: a31662a7c1

    executor >  local (5)
    [ce/df5eeb] process > FASTP (2)           [100%] 2 of 2 ✔
    [-        ] process > TRIMGALORE          -
    [d1/7d2b4b] process > GENERATE_REPORT (3) [100%] 3 of 3 ✔
    WARN: Low sequencing depth for sample_002: 25000000
    ```

Vemos um aviso sobre baixa profundidade de sequenciamento para uma das amostras.

### Conclusão

- **`error()`**: Para o fluxo de trabalho imediatamente com mensagem clara
- **`log.warn`**: Emite avisos sem parar o fluxo de trabalho
- **Validação antecipada**: Verifica entradas antes do processamento para falhar rapidamente com erros úteis
- **Funções de validação**: Crie lógica de validação reutilizável que pode ser chamada no início do fluxo de trabalho

A validação adequada torna os fluxos de trabalho mais robustos e amigáveis ao usuário, detectando problemas precocemente com mensagens de erro claras.

---

## 8. Manipuladores de Eventos do Fluxo de Trabalho

Até agora, escrevemos código em nossos scripts de fluxo de trabalho e definições de processos. Mas há mais um recurso importante que você deve conhecer: manipuladores de eventos do fluxo de trabalho.

Os manipuladores de eventos são closures que são executadas em pontos específicos do ciclo de vida do seu fluxo de trabalho. Eles são perfeitos para adicionar registro, notificações ou operações de limpeza. Esses manipuladores devem ser definidos no seu script de fluxo de trabalho junto com a definição do fluxo de trabalho.

### 8.1. O Manipulador `onComplete`

O manipulador de eventos mais comumente usado é `onComplete`, que é executado quando seu fluxo de trabalho termina (seja com sucesso ou falha). Vamos adicionar um para resumir os resultados do nosso pipeline.

Adicione o manipulador de eventos ao seu arquivo `main.nf`, dentro da definição do fluxo de trabalho:

=== "After"

    ```groovy title="main.nf" linenums="66" hl_lines="5-16"
        ch_fastp = FASTP(trim_branches.fastp)
        ch_trimgalore = TRIMGALORE(trim_branches.trimgalore)
        GENERATE_REPORT(ch_samples)

        workflow.onComplete = {
            println ""
            println "Pipeline execution summary:"
            println "=========================="
            println "Completed at: ${workflow.complete}"
            println "Duration    : ${workflow.duration}"
            println "Success     : ${workflow.success}"
            println "workDir     : ${workflow.workDir}"
            println "exit status : ${workflow.exitStatus}"
            println ""
        }
    }
    ```

=== "Before"

    ```groovy title="main.nf" linenums="66" hl_lines="4"
        ch_fastp = FASTP(trim_branches.fastp)
        ch_trimgalore = TRIMGALORE(trim_branches.trimgalore)
        GENERATE_REPORT(ch_samples)
    }
    ```

Esta closure é executada quando o fluxo de trabalho é concluído. Dentro dela, você tem acesso ao objeto `workflow`, que fornece propriedades úteis sobre a execução.

Execute seu fluxo de trabalho e você verá este resumo aparecer no final!

```bash
nextflow run main.nf --input ./data/samples.csv -ansi-log false
```

??? success "Saída do comando"

    ```console
    N E X T F L O W  ~  version 25.10.2
    Launching `main.nf` [marvelous_boltzmann] DSL2 - revision: a31662a7c1
    WARN: Low sequencing depth for sample_002: 25000000
    [9b/d48e40] Submitted process > FASTP (2)
    [6a/73867a] Submitted process > GENERATE_REPORT (2)
    [79/ad0ac5] Submitted process > GENERATE_REPORT (1)
    [f3/bda6cb] Submitted process > FASTP (1)
    [34/d5b52f] Submitted process > GENERATE_REPORT (3)

    Pipeline execution summary:
    ==========================
    Completed at: 2025-10-10T12:14:24.885384+01:00
    Duration    : 2.9s
    Success     : true
    workDir     : /workspaces/training/side-quests/essential_scripting_patterns/work
    exit status : 0
    ```

Vamos torná-lo mais útil adicionando lógica condicional:

=== "After"

    ```groovy title="main.nf" linenums="66" hl_lines="5-22"
        ch_fastp = FASTP(trim_branches.fastp)
        ch_trimgalore = TRIMGALORE(trim_branches.trimgalore)
        GENERATE_REPORT(ch_samples)

        workflow.onComplete = {
            println ""
            println "Pipeline execution summary:"
            println "=========================="
            println "Completed at: ${workflow.complete}"
            println "Duration    : ${workflow.duration}"
            println "Success     : ${workflow.success}"
            println "workDir     : ${workflow.workDir}"
            println "exit status : ${workflow.exitStatus}"
            println ""

            if (workflow.success) {
                println "✅ Pipeline completed successfully!"
            } else {
                println "❌ Pipeline failed!"
                println "Error: ${workflow.errorMessage}"
            }
        }
    }
    ```

=== "Before"

    ```groovy title="main.nf" linenums="66" hl_lines="5-16"
        ch_fastp = FASTP(trim_branches.fastp)
        ch_trimgalore = TRIMGALORE(trim_branches.trimgalore)
        GENERATE_REPORT(ch_samples)

        workflow.onComplete = {
            println ""
            println "Pipeline execution summary:"
            println "=========================="
            println "Completed at: ${workflow.complete}"
            println "Duration    : ${workflow.duration}"
            println "Success     : ${workflow.success}"
            println "workDir     : ${workflow.workDir}"
            println "exit status : ${workflow.exitStatus}"
            println ""
        }
    }
    ```

Agora obtemos um resumo ainda mais informativo, incluindo uma mensagem de sucesso/falha e o diretório de saída, se especificado:

<!-- TODO: add run command -->

??? success "Saída do comando"

    ```console
    N E X T F L O W  ~  version 25.10.2
    Launching `main.nf` [boring_linnaeus] DSL2 - revision: a31662a7c1
    WARN: Low sequencing depth for sample_002: 25000000
    [e5/242efc] Submitted process > FASTP (2)
    [3b/74047c] Submitted process > GENERATE_REPORT (3)
    [8a/7a57e6] Submitted process > GENERATE_REPORT (1)
    [a8/b1a31f] Submitted process > GENERATE_REPORT (2)
    [40/648429] Submitted process > FASTP (1)

    Pipeline execution summary:
    ==========================
    Completed at: 2025-10-10T12:16:00.522569+01:00
    Duration    : 3.6s
    Success     : true
    workDir     : /workspaces/training/side-quests/essential_scripting_patterns/work
    exit status : 0

    ✅ Pipeline completed successfully!
    ```

Você também pode escrever o resumo em um arquivo usando operações de arquivo:

```groovy title="main.nf - Writing summary to file"
workflow {
    // ... your workflow code ...

    workflow.onComplete = {
        def summary = """
        Pipeline Execution Summary
        ===========================
        Completed: ${workflow.complete}
        Duration : ${workflow.duration}
        Success  : ${workflow.success}
        Command  : ${workflow.commandLine}
        """

        println summary

        // Write to a log file
        def log_file = file("${workflow.launchDir}/pipeline_summary.txt")
        log_file.text = summary
    }
}
```

### 8.2. O Manipulador `onError`

Além de `onComplete`, há outro manipulador de eventos que você pode usar: `onError`, que é executado apenas se o fluxo de trabalho falhar:

```groovy title="main.nf - onError handler"
workflow {
    // ... your workflow code ...

    workflow.onError = {
        println "="* 50
        println "Pipeline execution failed!"
        println "Error message: ${workflow.errorMessage}"
        println "="* 50

        // Write detailed error log
        def error_file = file("${workflow.launchDir}/error.log")
        error_file.text = """
        Workflow Error Report
        =====================
        Time: ${new Date()}
        Error: ${workflow.errorMessage}
        Error report: ${workflow.errorReport ?: 'No detailed report available'}
        """

        println "Error details written to: ${error_file}"
    }
}
```

Você pode usar vários manipuladores juntos no seu script de fluxo de trabalho:

```groovy title="main.nf - Combined handlers"
workflow {
    // ... your workflow code ...

    workflow.onError = {
        println "Workflow failed: ${workflow.errorMessage}"
    }

    workflow.onComplete = {
        def duration_mins = workflow.duration.toMinutes().round(2)
        def status = workflow.success ? "SUCCESS ✅" : "FAILED ❌"

        println """
        Pipeline finished: ${status}
        Duration: ${duration_mins} minutes
        """
    }
}
```

### Conclusão

Nesta seção, você aprendeu:

- **Closures de manipulador de eventos**: Closures no seu script de fluxo de trabalho que são executadas em diferentes pontos do ciclo de vida
- **Manipulador `onComplete`**: Para resumos de execução e relatórios de resultados
- **Manipulador `onError`**: Para tratamento de erros e registro de falhas
- **Propriedades do objeto workflow**: Acessando `workflow.success`, `workflow.duration`, `workflow.errorMessage`, etc.

Os manipuladores de eventos mostram como você pode usar todo o poder da linguagem Nextflow dentro de seus scripts de fluxo de trabalho para adicionar capacidades sofisticadas de registro e notificação.

---

## Resumo

Parabéns, você conseguiu!

Ao longo desta busca lateral, você construiu um pipeline abrangente de processamento de amostras que evoluiu do manuseio básico de metadados para um fluxo de trabalho sofisticado e pronto para produção.
Cada seção construiu sobre a anterior, demonstrando como construções de programação transformam fluxos de trabalho simples em poderosos sistemas de processamento de dados, com os seguintes benefícios:

- **Código mais claro**: Entender o fluxo de dados vs scripting ajuda você a escrever fluxos de trabalho mais organizados
- **Manipulação robusta**: Operadores de navegação segura e Elvis tornam os fluxos de trabalho resilientes a dados ausentes
- **Processamento flexível**: Lógica condicional permite que seus fluxos de trabalho processem diferentes tipos de amostras adequadamente
- **Recursos adaptativos**: Diretivas dinâmicas otimizam o uso de recursos com base nas características de entrada

Esta progressão reflete a evolução do mundo real dos pipelines bioinformáticos, de protótipos de pesquisa que lidam com algumas amostras a sistemas de produção que processam milhares de amostras em laboratórios e instituições.
Cada desafio que você resolveu e padrão que aprendeu reflete problemas reais que os desenvolvedores enfrentam ao escalar fluxos de trabalho Nextflow.

Aplicar esses padrões em seu próprio trabalho permitirá que você construa fluxos de trabalho robustos e prontos para produção.

### Padrões-chave

1.  **Fluxo de dados vs Scripting:** Você aprendeu a distinguir entre operações de fluxo de dados (orquestração de canais) e scripting (código que manipula dados), incluindo as diferenças cruciais entre operações em diferentes tipos como `collect` em Canal vs Lista.

    - Fluxo de dados: orquestração de canais

    ```groovy
    channel.fromPath('*.fastq').splitCsv(header: true)
    ```

    - Scripting: processamento de dados em coleções

    ```groovy
    sample_data.collect { it.toUpperCase() }
    ```

2.  **Processamento Avançado de Strings**: Você dominou expressões regulares para análise de nomes de arquivos, geração dinâmica de scripts em processos e interpolação de variáveis (Nextflow vs Bash vs Shell).

    - Correspondência de padrões

    ```groovy
    filename =~ ~/^(\w+)_(\w+)_(\d+)\.fastq$/
    ```

    - Função com retorno condicional

    ```groovy
    def parseSample(filename) {
        def matcher = filename =~ pattern
        return matcher ? [valid: true, data: matcher[0]] : [valid: false]
    }
    ```

    - Coleção de arquivos para argumentos de comando (em bloco de script de processo)

    ```groovy
    script:
    def file_args = input_files.collect { file -> "--input ${file}" }.join(' ')
    """
    analysis_tool ${file_args} --output results.txt
    """
    ```

3.  **Criando Funções Reutilizáveis**: Você aprendeu a extrair lógica complexa em funções nomeadas que podem ser chamadas de operadores de canal, tornando os fluxos de trabalho mais legíveis e de fácil manutenção.

    - Definir uma função nomeada

    ```groovy
    def separateMetadata(row) {
        def sample_meta = [ /* code hidden for brevity */ ]
        def fastq_path = file(row.file_path)
        def m = (fastq_path.name =~ /^(.+)_S(\d+)_L(\d{3})_(R[12])_(\d{3})\.fastq(?:\.gz)?$/)
        def file_meta = m ? [ /* code hidden for brevity */ ] : [:]
        def priority = sample_meta.quality > 40 ? 'high' : 'normal'

        return tuple(sample_meta + file_meta + [priority: priority], fastq_path)
    }
    ```

    - Chamar a função nomeada em um fluxo de trabalho

    ```groovy
    workflow {
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map{ row -> separateMetadata(row) }

        ch_fastp = FASTP(ch_samples)
    }
    ```

4.  **Diretivas de Recursos Dinâmicos com Closures**: Você explorou o uso de closures em diretivas de processo para alocação adaptativa de recursos com base nas características de entrada.

    - Closures nomeadas e composição

    ```groovy
    def enrichData = normalizeId >> addQualityCategory >> addFlags
    def processor = generalFunction.curry(fixedParam)
    ```

    - Closures com acesso a escopo

    ```groovy
    def collectStats = { data -> stats.count++; return data }
    ```

5.  **Lógica Condicional e Controle de Processo**: Você adicionou roteamento inteligente usando os operadores `.branch()` e `.filter()`, aproveitando a truthiness para expressões condicionais concisas.

    - Use `.branch()` para rotear dados através de diferentes ramificações do fluxo de trabalho

    ```groovy
    trim_branches = ch_samples
    .branch { meta, reads ->
        fastp: meta.organism == 'human' && meta.depth >= 30000000
        trimgalore: true
    }

    ch_fastp = FASTP(trim_branches.fastp)
    ch_trimgalore = TRIMGALORE(trim_branches.trimgalore)
    ```

    - Avaliação booleana com Groovy Truth

    ```groovy
    if (sample.files) println "Has files"
    ```

    - Use `filter()` para filtrar dados com 'truthiness'

    ```groovy
    ch_valid_samples = ch_samples
        .filter { meta, reads ->
            meta.id && meta.organism && meta.depth >= 25000000
        }
    ```

6.  **Operadores de Navegação Segura e Elvis**: Você tornou o pipeline robusto contra dados ausentes usando `?.` para acesso seguro a propriedades e `?:` para fornecer valores padrão.

    ```groovy
    def id = data?.sample?.id ?: 'unknown'
    ```

7.  **Validação com error() e log.warn**: Você aprendeu a validar entradas cedo e falhar rapidamente com mensagens de erro claras.

    ```groovy
    try {
        def errors = validateSample(sample)
        if (errors) throw new RuntimeException("Invalid: ${errors.join(', ')}")
    } catch (Exception e) {
        println "Error: ${e.message}"
    }
    ```

8.  **Manipuladores de Eventos de Configuração**: Você aprendeu a usar manipuladores de eventos de fluxo de trabalho (`onComplete` e `onError`) para registro, notificações e gerenciamento do ciclo de vida.

    - Usando `onComplete` para registrar e notificar

    ```groovy
    workflow.onComplete = {
        println "Success     : ${workflow.success}"
        println "exit status : ${workflow.exitStatus}"

        if (workflow.success) {
            println "✅ Pipeline completed successfully!"
        } else {
            println "❌ Pipeline failed!"
            println "Error: ${workflow.errorMessage}"
        }
    }
    ```

    - Usando `onError` para tomar ações especificamente em caso de falha

    ```groovy
    workflow.onError = {
        // Write detailed error log
        def error_file = file("${workflow.launchDir}/error.log")
        error_file.text = """
        Time: ${new Date()}
        Error: ${workflow.errorMessage}
        Error report: ${workflow.errorReport ?: 'No detailed report available'}
        """

        println "Error details written to: ${error_file}"
    }
    ```

### Recursos adicionais

- [Referência da Linguagem Nextflow](https://nextflow.io/docs/latest/reference/syntax.html)
- [Operadores Nextflow](https://www.nextflow.io/docs/latest/operator.html)
- [Sintaxe de Script Nextflow](https://www.nextflow.io/docs/latest/script.html)
- [Biblioteca Padrão Nextflow](https://nextflow.io/docs/latest/reference/stdlib.html)

Certifique-se de verificar esses recursos quando precisar explorar recursos mais avançados.

Você se beneficiará praticando e expandindo suas habilidades para:

- Escrever fluxos de trabalho mais limpos com separação adequada entre fluxo de dados e scripting
- Dominar a interpolação de variáveis para evitar armadilhas comuns com variáveis Nextflow, Bash e shell
- Usar diretivas de recursos dinâmicos para fluxos de trabalho eficientes e adaptativos
- Transformar coleções de arquivos em argumentos de linha de comando formatados corretamente
- Lidar com diferentes convenções de nomenclatura de arquivos e formatos de entrada graciosamente usando regex e processamento de strings
- Criar código reutilizável e de fácil manutenção usando padrões avançados de closure e programação funcional
- Processar e organizar conjuntos de dados complexos usando operações de coleção
- Adicionar validação, tratamento de erros e registro para tornar seus fluxos de trabalho prontos para produção
- Implementar gerenciamento do ciclo de vida do fluxo de trabalho com manipuladores de eventos

---

## O que vem a seguir?

Retorne ao [menu de Buscas Laterais](./index.md) ou clique no botão no canto inferior direito da página para avançar para o próximo tópico da lista.
