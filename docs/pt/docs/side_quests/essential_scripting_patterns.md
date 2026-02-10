# Padrões Essenciais de Script Nextflow

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tradução assistida por IA - [saiba mais e sugira melhorias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Nextflow é uma linguagem de programação que roda na Máquina Virtual Java. Embora Nextflow seja construído em [Groovy](http://groovy-lang.org/) e compartilhe muito de sua sintaxe, Nextflow é mais do que apenas "Groovy com extensões" -- é uma linguagem independente com uma [sintaxe](https://nextflow.io/docs/latest/reference/syntax.html) e [biblioteca padrão](https://nextflow.io/docs/latest/reference/stdlib.html) totalmente especificadas.

Você pode escrever muito Nextflow sem se aventurar além da sintaxe básica para variáveis, mapas e listas. A maioria dos tutoriais de Nextflow foca na orquestração de fluxo de trabalho (canais, processos e fluxo de dados), e você pode ir surpreendentemente longe apenas com isso.

No entanto, quando você precisa manipular dados, analisar nomes de arquivos complexos, implementar lógica condicional ou construir fluxos de trabalho robustos para produção, ajuda pensar em dois aspectos distintos do seu código: **fluxo de dados** (canais, operadores, processos e fluxos de trabalho) e **script** (o código dentro de closures, funções e scripts de processo). Embora essa distinção seja um tanto arbitrária—é tudo código Nextflow—ela fornece um modelo mental útil para entender quando você está orquestrando seu pipeline versus quando está manipulando dados. Dominar ambos melhora dramaticamente sua capacidade de escrever fluxos de trabalho claros e mantíveis.

### Objetivos de aprendizado

Esta missão secundária leva você em uma jornada prática desde conceitos básicos até padrões prontos para produção.
Vamos transformar um fluxo de trabalho simples de leitura de CSV em um pipeline sofisticado de bioinformática, evoluindo-o passo a passo através de desafios realistas:

- **Entendendo limites:** Distinguir entre operações de fluxo de dados e script, e entender como eles trabalham juntos
- **Manipulação de dados:** Extrair, transformar e subconjuntar mapas e coleções usando operadores poderosos
- **Processamento de strings:** Analisar esquemas complexos de nomenclatura de arquivos com padrões regex e dominar interpolação de variáveis
- **Funções reutilizáveis:** Extrair lógica complexa em funções nomeadas para fluxos de trabalho mais limpos e mantíveis
- **Lógica dinâmica:** Construir processos que se adaptam a diferentes tipos de entrada e usar closures para alocação dinâmica de recursos
- **Roteamento condicional:** Rotear amostras inteligentemente através de diferentes processos com base em suas características de metadados
- **Operações seguras:** Lidar graciosamente com dados ausentes usando operadores null-safe e validar entradas com mensagens de erro claras
- **Manipuladores baseados em configuração:** Usar manipuladores de eventos de fluxo de trabalho para logging, notificações e gerenciamento de ciclo de vida

### Pré-requisitos

Antes de assumir esta missão secundária, você deve:

- Ter completado o tutorial [Hello Nextflow](../hello_nextflow/README.md) ou curso equivalente para iniciantes.
- Estar confortável usando conceitos e mecanismos básicos do Nextflow (processos, canais, operadores, trabalhando com arquivos, metadados)
- Ter familiaridade básica com construções comuns de programação (variáveis, mapas, listas)

Este tutorial explicará conceitos de programação conforme os encontramos, então você não precisa de experiência extensa em programação.
Começaremos com conceitos fundamentais e construiremos até padrões avançados.

---

## 0. Primeiros passos

#### Abra o codespace de treinamento

Se você ainda não o fez, certifique-se de abrir o ambiente de treinamento conforme descrito na [Configuração do Ambiente](../envsetup/index.md).

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### Mova-se para o diretório do projeto

Vamos mover para o diretório onde os arquivos para este tutorial estão localizados.

```bash
cd side-quests/essential_scripting_patterns
```

#### Revise os materiais

Você encontrará um arquivo de fluxo de trabalho principal e um diretório `data` contendo arquivos de dados de exemplo.

```console title="Directory contents"
.
├── collect.nf
├── data
│   ├── samples.csv
│   └── sequences
│       ├── SAMPLE_001_S1_L001_R1_001.fastq
│       ├── SAMPLE_002_S2_L001_R1_001.fastq
│       └── SAMPLE_003_S3_L001_R1_001.fastq
├── main.nf
├── modules
│   ├── fastp.nf
│   ├── generate_report.nf
│   └── trimgalore.nf
└── nextflow.config
```

Nosso CSV de amostra contém informações sobre amostras biológicas que precisam de processamento diferente com base em suas características:

```console title="samples.csv"
sample_id,organism,tissue_type,sequencing_depth,file_path,quality_score
SAMPLE_001,human,liver,30000000,data/sequences/SAMPLE_001_S1_L001_R1_001.fastq,38.5
SAMPLE_002,mouse,brain,25000000,data/sequences/SAMPLE_002_S2_L001_R1_001.fastq,35.2
SAMPLE_003,human,kidney,45000000,data/sequences/SAMPLE_003_S3_L001_R1_001.fastq,42.1
```

Usaremos este conjunto de dados realista para explorar técnicas práticas de programação que você encontrará em fluxos de trabalho reais de bioinformática.

<!-- TODO: Can we make this more domain-agnostic? -->

<!-- TODO: add an assignment statement? #### Review the assignment -->

#### Lista de verificação de prontidão

Acha que está pronto para mergulhar?

- [ ] Entendo o objetivo deste curso e seus pré-requisitos
- [ ] Meu codespace está funcionando
- [ ] Defini meu diretório de trabalho apropriadamente
<!-- - [ ] I understand the assignment -->

Se você pode marcar todas as caixas, está pronto para começar.

---

## 1. Fluxo de Dados vs Script: Entendendo os Limites

### 1.1. Identificando o Que é o Quê

Ao escrever fluxos de trabalho Nextflow, é importante distinguir entre **fluxo de dados** (como os dados se movem através de canais e processos) e **script** (o código que manipula dados e toma decisões). Vamos construir um fluxo de trabalho demonstrando como eles trabalham juntos.

#### 1.1.1. Fluxo de Trabalho Básico Nextflow

Comece com um fluxo de trabalho simples que apenas lê o arquivo CSV (já fizemos isso para você em `main.nf`):

```groovy title="main.nf" linenums="1"
workflow {
    ch_samples = channel.fromPath("./data/samples.csv")
        .splitCsv(header: true)
        .view()
}
```

O bloco `workflow` define nossa estrutura de pipeline, enquanto `channel.fromPath()` cria um canal a partir de um caminho de arquivo. O operador `.splitCsv()` processa o arquivo CSV e converte cada linha em uma estrutura de dados de mapa.

Execute este fluxo de trabalho para ver os dados brutos do CSV:

```bash
nextflow run main.nf
```

??? success "Saída do comando"

    ```console
    Launching `main.nf` [marvelous_tuckerman] DSL2 - revision: 6113e05c17

    [sample_id:SAMPLE_001, organism:human, tissue_type:liver, sequencing_depth:30000000, file_path:data/sequences/SAMPLE_001_S1_L001_R1_001.fastq, quality_score:38.5]
    [sample_id:SAMPLE_002, organism:mouse, tissue_type:brain, sequencing_depth:25000000, file_path:data/sequences/SAMPLE_002_S2_L001_R1_001.fastq, quality_score:35.2]
    [sample_id:SAMPLE_003, organism:human, tissue_type:kidney, sequencing_depth:45000000, file_path:data/sequences/SAMPLE_003_S3_L001_R1_001.fastq, quality_score:42.1]
    ```

#### 1.1.2. Adicionando o Operador Map

Agora vamos adicionar script para transformar os dados, usando o operador `.map()` que você provavelmente já conhece. Este operador recebe uma 'closure' onde podemos escrever código para transformar cada item.

!!! note "Nota"

    Uma **closure** é um bloco de código que pode ser passado e executado depois. Pense nela como uma função que você define inline. Closures são escritas com chaves `{ }` e podem receber parâmetros. Elas são fundamentais para como os operadores Nextflow funcionam e se você tem escrito Nextflow há algum tempo, pode já ter estado usando-as sem perceber!

Aqui está como essa operação map se parece:

=== "Depois"

    ```groovy title="main.nf" linenums="2" hl_lines="3-6"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row ->
                return row
            }
            .view()
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="2" hl_lines="3"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .view()
    ```

Esta é nossa primeira **closure** - uma função anônima que você pode passar como argumento (similar a lambdas em Python ou arrow functions em JavaScript). Closures são essenciais para trabalhar com operadores Nextflow.

A closure `{ row -> return row }` recebe um parâmetro `row` (poderia ser qualquer nome: `item`, `sample`, etc.).

Quando o operador `.map()` processa cada item do canal, ele passa esse item para sua closure. Aqui, `row` contém uma linha CSV por vez.

Aplique esta mudança e execute o fluxo de trabalho:

```bash
nextflow run main.nf
```

Você verá a mesma saída de antes, porque estamos simplesmente retornando a entrada inalterada. Isso confirma que o operador map está funcionando corretamente. Agora vamos começar a transformar os dados.

#### 1.1.3. Criando uma Estrutura de Dados Map

Agora vamos escrever lógica de **script** dentro de nossa closure para transformar cada linha de dados. É aqui que processamos itens de dados individuais em vez de orquestrar o fluxo de dados.

=== "Depois"

    ```groovy title="main.nf" linenums="2" hl_lines="4-12"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row ->
                // Script para transformação de dados
                def sample_meta = [
                    id: row.sample_id.toLowerCase(),
                    organism: row.organism,
                    tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                    depth: row.sequencing_depth.toInteger(),
                    quality: row.quality_score.toDouble()
                ]
                return sample_meta
            }
            .view()
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="2" hl_lines="4"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row ->
                return row
            }
            .view()
    ```

O mapa `sample_meta` é uma estrutura de dados chave-valor (como dicionários em Python, objetos em JavaScript, ou hashes em Ruby) armazenando informações relacionadas: ID da amostra, organismo, tipo de tecido, profundidade de sequenciamento e pontuação de qualidade.

Usamos métodos de manipulação de strings como `.toLowerCase()` e `.replaceAll()` para limpar nossos dados, e métodos de conversão de tipo como `.toInteger()` e `.toDouble()` para converter dados de string do CSV nos tipos numéricos apropriados.

Aplique esta mudança e execute o fluxo de trabalho:

```bash
nextflow run main.nf
```

??? success "Saída do comando"

    ```console
    [id:sample_001, organism:human, tissue:liver, depth:30000000, quality:38.5]
    [id:sample_002, organism:mouse, tissue:brain, depth:25000000, quality:35.2]
    [id:sample_003, organism:human, tissue:kidney, depth:45000000, quality:42.1]
    ```

#### 1.1.4. Adicionando Lógica Condicional

Agora vamos adicionar mais script - desta vez usando um operador ternário para tomar decisões com base em valores de dados.

Faça a seguinte mudança:

=== "Depois"

    ```groovy title="main.nf" linenums="2" hl_lines="11-12"
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
                def priority = sample_meta.quality > 40 ? 'high' : 'normal'
                return sample_meta + [priority: priority]
            }
            .view()
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="2" hl_lines="11"
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
                return sample_meta
            }
            .view()
    ```

O operador ternário é uma forma abreviada de uma instrução if/else que segue o padrão `condição ? valor_se_verdadeiro : valor_se_falso`. Esta linha significa: "Se a qualidade for maior que 40, use 'high', caso contrário use 'normal'". Seu primo, o **operador Elvis** (`?:`), fornece valores padrão quando algo é null ou vazio - exploraremos esse padrão mais tarde neste tutorial.

O operador de adição de mapa `+` cria um **novo mapa** em vez de modificar o existente. Esta linha cria um novo mapa que contém todos os pares chave-valor de `sample_meta` mais a nova chave `priority`.

!!! Note "Nota"

    Nunca modifique mapas passados para closures - sempre crie novos usando `+` (por exemplo). No Nextflow, os mesmos dados frequentemente fluem através de múltiplas operações simultaneamente. Modificar um mapa in-place pode causar efeitos colaterais imprevisíveis quando outras operações referenciam esse mesmo objeto. Criar novos mapas garante que cada operação tenha sua própria cópia limpa.

Execute o fluxo de trabalho modificado:

```bash
nextflow run main.nf
```

??? success "Saída do comando"

    ```console
    [id:sample_001, organism:human, tissue:liver, depth:30000000, quality:38.5, priority:normal]
    [id:sample_002, organism:mouse, tissue:brain, depth:25000000, quality:35.2, priority:normal]
    [id:sample_003, organism:human, tissue:kidney, depth:45000000, quality:42.1, priority:high]
    ```

Adicionamos com sucesso lógica condicional para enriquecer nossos metadados com um nível de prioridade baseado em pontuações de qualidade.

#### 1.1.5. Subconjuntando Mapas com `.subMap()`

Enquanto o operador `+` adiciona chaves a um mapa, às vezes você precisa fazer o oposto - extrair apenas chaves específicas. O método `.subMap()` é perfeito para isso.

Vamos adicionar uma linha para criar uma versão simplificada de nossos metadados que contém apenas campos de identificação:

=== "Depois"

    ```groovy title="main.nf" linenums="2" hl_lines="12-15"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row ->
                // Script para transformação de dados
                def sample_meta = [
                    id: row.sample_id.toLowerCase(),
                    organism: row.organism,
                    tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                    depth: row.sequencing_depth.toInteger(),
                    quality: row.quality_score.toDouble()
                ]
                def id_only = sample_meta.subMap(['id', 'organism', 'tissue'])
                println "ID fields only: ${id_only}"

                def priority = sample_meta.quality > 40 ? 'high' : 'normal'
                return sample_meta + [priority: priority]
            }
            .view()
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="2" hl_lines="12"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row ->
                // Script para transformação de dados
                def sample_meta = [
                    id: row.sample_id.toLowerCase(),
                    organism: row.organism,
                    tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                    depth: row.sequencing_depth.toInteger(),
                    quality: row.quality_score.toDouble()
                ]
                def priority = sample_meta.quality > 40 ? 'high' : 'normal'
                return sample_meta + [priority: priority]
            }
            .view()
    ```

Execute o fluxo de trabalho modificado:

```bash
nextflow run main.nf
```

??? success "Saída do comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [peaceful_cori] DSL2 - revision: 4cc4a8340f

    ID fields only: [id:sample_001, organism:human, tissue:liver]
    ID fields only: [id:sample_002, organism:mouse, tissue:brain]
    ID fields only: [id:sample_003, organism:human, tissue:kidney]
    [id:sample_001, organism:human, tissue:liver, depth:30000000, quality:38.5, priority:normal]
    [id:sample_002, organism:mouse, tissue:brain, depth:25000000, quality:35.2, priority:normal]
    [id:sample_003, organism:human, tissue:kidney, depth:45000000, quality:42.1, priority:high]
    ```

Isso mostra tanto os metadados completos exibidos pela operação `view()` quanto o subconjunto extraído que imprimimos com `println`.

O método `.subMap()` recebe uma lista de chaves e retorna um novo mapa contendo apenas essas chaves. Se uma chave não existe no mapa original, ela simplesmente não é incluída no resultado.

Isso é particularmente útil quando você precisa criar diferentes versões de metadados para diferentes processos - alguns podem precisar de metadados completos enquanto outros precisam apenas de campos mínimos de identificação.

Agora remova essas instruções println para restaurar seu fluxo de trabalho ao estado anterior, pois não precisamos delas daqui para frente.

!!! tip "Resumo de Operações com Mapas"

    - **Adicionar chaves**: `map1 + [new_key: value]` - Cria novo mapa com chaves adicionais
    - **Extrair chaves**: `map1.subMap(['key1', 'key2'])` - Cria novo mapa com apenas chaves especificadas
    - **Ambas operações criam novos mapas** - Mapas originais permanecem inalterados

#### 1.1.6. Combinando Mapas e Retornando Resultados

Até agora, só temos retornado o que a comunidade Nextflow chama de 'meta map', e temos ignorado os arquivos aos quais esses metadados se relacionam. Mas se você está escrevendo fluxos de trabalho Nextflow, provavelmente quer fazer algo com esses arquivos.

Vamos produzir uma estrutura de canal compreendendo uma tupla de 2 elementos: o mapa de metadados enriquecido e o caminho de arquivo correspondente. Este é um padrão comum no Nextflow para passar dados para processos.

=== "Depois"

    ```groovy title="main.nf" linenums="2" hl_lines="12"
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
                def priority = sample_meta.quality > 40 ? 'high' : 'normal'
                return tuple( sample_meta + [priority: priority], file(row.file_path) )
            }
            .view()
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="2" hl_lines="12"
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
                def priority = sample_meta.quality > 40 ? 'high' : 'normal'
                return sample_meta + [priority: priority]
            }
            .view()
    ```

Aplique esta mudança e execute o fluxo de trabalho:

```bash
nextflow run main.nf
```

??? success "Saída do comando"

    ```console
    [[id:sample_001, organism:human, tissue:liver, depth:30000000, quality:38.5, priority:normal], /workspaces/training/side-quests/essential_scripting_patterns/data/sequences/SAMPLE_001_S1_L001_R1_001.fastq]
    [[id:sample_002, organism:mouse, tissue:brain, depth:25000000, quality:35.2, priority:normal], /workspaces/training/side-quests/essential_scripting_patterns/data/sequences/SAMPLE_002_S2_L001_R1_001.fastq]
    [[id:sample_003, organism:human, tissue:kidney, depth:45000000, quality:42.1, priority:high], /workspaces/training/side-quests/essential_scripting_patterns/data/sequences/SAMPLE_003_S3_L001_R1_001.fastq]
    ```

Esta estrutura de tupla `[meta, file]` é um padrão comum no Nextflow para passar tanto metadados quanto arquivos associados para processos.

!!! note "Nota"

    **Mapas e Metadados**: Mapas são fundamentais para trabalhar com metadados no Nextflow. Para uma explicação mais detalhada sobre trabalhar com mapas de metadados, veja a missão secundária [Trabalhando com metadados](./metadata.md).

Nosso fluxo de trabalho demonstra o padrão central: **operações de fluxo de dados** (`workflow`, `channel.fromPath()`, `.splitCsv()`, `.map()`, `.view()`) orquestram como os dados se movem através do pipeline, enquanto **script** (mapas `[key: value]`, métodos de string, conversões de tipo, operadores ternários) dentro da closure `.map()` lida com a transformação de itens de dados individuais.

### 1.2. Entendendo Diferentes Tipos: Canal vs Lista

Até agora, tudo bem, podemos distinguir entre operações de fluxo de dados e script. Mas e quando o mesmo nome de método existe em ambos os contextos?

Um exemplo perfeito é o método `collect`, que existe tanto para tipos de canal quanto para tipos de Lista na biblioteca padrão do Nextflow. O método `collect()` em uma Lista transforma cada elemento, enquanto o operador `collect()` em um canal reúne todas as emissões do canal em um canal de item único.

Vamos demonstrar isso com alguns dados de amostra, começando por refrescar o que o operador `collect()` de canal faz. Confira `collect.nf`:

```groovy title="collect.nf" linenums="1"
def sample_ids = ['sample_001', 'sample_002', 'sample_003']

// channel.collect() - agrupa múltiplas emissões de canal em uma
ch_input = channel.fromList(sample_ids)
ch_input.view { sample -> "Individual channel item: ${sample}" }
ch_collected = ch_input.collect()
ch_collected.view { list -> "channel.collect() result: ${list} (${list.size()} items grouped into 1)" }
```

Passos:

- Define uma Lista de IDs de amostra
- Cria um canal com `fromList()` que emite cada ID de amostra separadamente
- Imprime cada item com `view()` conforme flui
- Reúne todos os itens em uma única lista com o operador `collect()` do canal
- Imprime o resultado coletado (item único contendo todos os IDs de amostra) com um segundo `view()`

Mudamos a estrutura do canal, mas não mudamos os dados em si.

Execute o fluxo de trabalho para confirmar isso:

```bash
nextflow run collect.nf
```

??? success "Saída do comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `collect.nf` [loving_mendel] DSL2 - revision: e8d054a46e

    Individual channel item: sample_001
    Individual channel item: sample_002
    Individual channel item: sample_003
    channel.collect() result: [sample_001, sample_002, sample_003] (3 items grouped into 1)
    ```

`view()` retorna uma saída para cada emissão de canal, então sabemos que esta saída única contém todos os 3 itens originais agrupados em uma lista.

Agora vamos ver o método `collect` em uma Lista em ação. Modifique `collect.nf` para aplicar o método `collect` da Lista à lista original de IDs de amostra:

=== "Depois"

    ```groovy title="main.nf" linenums="1" hl_lines="9-13"
    def sample_ids = ['sample_001', 'sample_002', 'sample_003']

    // channel.collect() - agrupa múltiplas emissões de canal em uma
    ch_input = channel.fromList(sample_ids)
    ch_input.view { sample -> "Individual channel item: ${sample}" }
    ch_collected = ch_input.collect()
    ch_collected.view { list -> "channel.collect() result: ${list} (${list.size()} items grouped into 1)" }

    // List.collect() - transforma cada elemento, preserva estrutura
    def formatted_ids = sample_ids.collect { id ->
        id.toUpperCase().replace('SAMPLE_', 'SPECIMEN_')
    }
    println "List.collect() result: ${formatted_ids} (${sample_ids.size()} items transformed into ${formatted_ids.size()})"
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="1"
    def sample_ids = ['sample_001', 'sample_002', 'sample_003']

    // channel.collect() - agrupa múltiplas emissões de canal em uma
    ch_input = channel.fromList(sample_ids)
    ch_input.view { sample -> "Individual channel item: ${sample}" }
    ch_collected = ch_input.collect()
    ch_collected.view { list -> "channel.collect() result: ${list} (${list.size()} items grouped into 1)" }
    ```

Neste novo trecho nós:

- Definimos uma nova variável `formatted_ids` que usa o método `collect` da Lista para transformar cada ID de amostra na lista original
- Imprimimos o resultado usando `println`

Execute o fluxo de trabalho modificado:

```bash
nextflow run collect.nf
```

??? success "Saída do comando"

    ```console hl_lines="5"
    N E X T F L O W   ~  version 25.10.2

    Launching `collect.nf` [cheeky_stonebraker] DSL2 - revision: 2d5039fb47

    List.collect() result: [SPECIMEN_001, SPECIMEN_002, SPECIMEN_003] (3 items transformed into 3)
    Individual channel item: sample_001
    Individual channel item: sample_002
    Individual channel item: sample_003
    channel.collect() result: [sample_001, sample_002, sample_003] (3 items grouped into 1)
    ```

Desta vez, NÃO mudamos a estrutura dos dados, ainda temos 3 itens na lista, mas TRANSFORMAMOS cada item usando o método `collect` da Lista para produzir uma nova lista com valores modificados. Isso é similar a usar o operador `map` em um canal, mas está operando em uma estrutura de dados Lista em vez de um canal.

`collect` é um caso extremo que estamos usando aqui para fazer um ponto. A lição chave é que quando você está escrevendo fluxos de trabalho, sempre distinga entre **estruturas de dados** (Listas, Mapas, etc.) e **canais** (construções de fluxo de dados). Operações podem compartilhar nomes mas se comportar completamente diferente dependendo do tipo em que são chamadas.

### 1.3. O Operador Spread (`*.`) - Atalho para Extração de Propriedades

Relacionado ao método `collect` da Lista está o operador spread (`*.`), que fornece uma maneira concisa de extrair propriedades de coleções. É essencialmente açúcar sintático para um padrão comum de `collect`.

Vamos adicionar uma demonstração ao nosso arquivo `collect.nf`:

=== "Depois"

    ```groovy title="collect.nf" linenums="1" hl_lines="15-18"
    def sample_ids = ['sample_001', 'sample_002', 'sample_003']

    // channel.collect() - agrupa múltiplas emissões de canal em uma
    ch_input = channel.fromList(sample_ids)
    ch_input.view { sample -> "Individual channel item: ${sample}" }
    ch_collected = ch_input.collect()
    ch_collected.view { list -> "channel.collect() result: ${list} (${list.size()} items grouped into 1)" }

    // List.collect() - transforma cada elemento, preserva estrutura
    def formatted_ids = sample_ids.collect { id ->
        id.toUpperCase().replace('SAMPLE_', 'SPECIMEN_')
    }
    println "List.collect() result: ${formatted_ids} (${sample_ids.size()} items transformed into ${formatted_ids.size()})"

    // Operador spread - acesso conciso a propriedades
    def sample_data = [[id: 's1', quality: 38.5], [id: 's2', quality: 42.1], [id: 's3', quality: 35.2]]
    def all_ids = sample_data*.id
    println "Spread operator result: ${all_ids}"
    ```

=== "Antes"

    ```groovy title="collect.nf" linenums="1"
    def sample_ids = ['sample_001', 'sample_002', 'sample_003']

    // channel.collect() - agrupa múltiplas emissões de canal em uma
    ch_input = channel.fromList(sample_ids)
    ch_input.view { sample -> "Individual channel item: ${sample}" }
    ch_collected = ch_input.collect()
    ch_collected.view { list -> "channel.collect() result: ${list} (${list.size()} items grouped into 1)" }

    // List.collect() - transforma cada elemento, preserva estrutura
    def formatted_ids = sample_ids.collect { id ->
        id.toUpperCase().replace('SAMPLE_', 'SPECIMEN_')
    }
    println "List.collect() result: ${formatted_ids} (${sample_ids.size()} items transformed into ${formatted_ids.size()})"
    ```

Execute o fluxo de trabalho atualizado:

```bash title="Test spread operator"
nextflow run collect.nf
```

??? success "Saída do comando"

    ```console hl_lines="6"
    N E X T F L O W   ~  version 25.10.2

    Launching `collect.nf` [cranky_galileo] DSL2 - revision: 5f3c8b2a91

    List.collect() result: [SPECIMEN_001, SPECIMEN_002, SPECIMEN_003] (3 items transformed into 3)
    Spread operator result: [s1, s2, s3]
    Individual channel item: sample_001
    Individual channel item: sample_002
    Individual channel item: sample_003
    channel.collect() result: [sample_001, sample_002, sample_003] (3 items grouped into 1)
    ```

O operador spread `*.` é um atalho para um padrão comum de collect:

```groovy
// Estes são equivalentes:
def ids = samples*.id
def ids = samples.collect { it.id }

// Também funciona com chamadas de método:
def names = files*.getName()
def names = files.collect { it.getName() }
```

O operador spread é particularmente útil quando você precisa extrair uma única propriedade de uma lista de objetos - é mais legível do que escrever a closure `collect` completa.

!!! tip "Quando Usar Spread vs Collect"

    - **Use spread (`*.`)** para acesso simples a propriedades: `samples*.id`, `files*.name`
    - **Use collect** para transformações ou lógica complexa: `samples.collect { it.id.toUpperCase() }`, `samples.collect { [it.id, it.quality > 40] }`

### Conclusão

Nesta seção, você aprendeu:

- **Fluxo de dados vs script**: Operadores de canal orquestram como os dados fluem através do seu pipeline, enquanto script transforma itens de dados individuais
- **Entendendo tipos**: O mesmo nome de método (como `collect`) pode se comportar diferentemente dependendo do tipo em que é chamado (Canal vs Lista)
- **Contexto importa**: Sempre esteja ciente se você está trabalhando com canais (fluxo de dados) ou estruturas de dados (script)

Entender esses limites é essencial para depuração, documentação e escrever fluxos de trabalho mantíveis.

A seguir, mergulharemos mais profundamente nas capacidades de processamento de strings, que são essenciais para lidar com dados do mundo real.

---

## 2. Processamento de Strings e Geração Dinâmica de Scripts

Dominar o processamento de strings separa fluxos de trabalho frágeis de pipelines robustos. Esta seção cobre análise de nomes de arquivos complexos, geração dinâmica de scripts e interpolação de variáveis.

### 2.1. Correspondência de Padrões e Expressões Regulares

Arquivos de bioinformática frequentemente têm convenções de nomenclatura complexas codificando metadados. Vamos extrair isso automaticamente usando correspondência de padrões com expressões regulares.

Vamos retornar ao nosso fluxo de trabalho `main.nf` e adicionar alguma lógica de correspondência de padrões para extrair informações adicionais de amostra dos nomes de arquivos. Os arquivos FASTQ em nosso conjunto de dados seguem convenções de nomenclatura estilo Illumina com nomes como `SAMPLE_001_S1_L001_R1_001.fastq.gz`. Estes podem parecer crípticos, mas na verdade codificam metadados úteis como ID da amostra, número da lane e direção de leitura. Vamos usar capacidades regex para analisar esses nomes.

Faça a seguinte mudança no seu fluxo de trabalho `main.nf` existente:

=== "Depois"

    ```groovy title="main.nf" linenums="4" hl_lines="10-21"
            .map { row ->
                // Script para transformação de dados
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
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="4" hl_lines="10-11"
            .map { row ->
                // Script para transformação de dados
                def sample_meta = [
                    id: row.sample_id.toLowerCase(),
                    organism: row.organism,
                    tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                    depth: row.sequencing_depth.toInteger(),
                    quality: row.quality_score.toDouble()
                ]
                def priority = sample_meta.quality > 40 ? 'high' : 'normal'
                return tuple(sample_meta + [priority: priority], file(row.file_path))
            }
    ```

Isso demonstra **conceitos chave de processamento de strings**:

1. **Literais de expressão regular** usando sintaxe `~/padrão/` - isso cria um padrão regex sem precisar escapar barras invertidas
2. **Correspondência de padrões** com o operador `=~` - isso tenta corresponder uma string contra um padrão regex
3. **Objetos Matcher** que capturam grupos com `[0][1]`, `[0][2]`, etc. - `[0]` refere-se à correspondência inteira, `[1]`, `[2]`, etc. referem-se a grupos capturados em parênteses

Vamos decompor o padrão regex `^(.+)_S(\d+)_L(\d{3})_(R[12])_(\d{3})\.fastq(?:\.gz)?$`:

| Padrão              | Corresponde                                 | Captura                             |
| ------------------- | ------------------------------------------- | ----------------------------------- |
| `^(.+)`             | Nome da amostra desde o início              | Grupo 1: nome da amostra            |
| `_S(\d+)`           | Número da amostra `_S1`, `_S2`, etc.        | Grupo 2: número da amostra          |
| `_L(\d{3})`         | Número da lane `_L001`                      | Grupo 3: lane (3 dígitos)           |
| `_(R[12])`          | Direção de leitura `_R1` ou `_R2`           | Grupo 4: direção de leitura         |
| `_(\d{3})`          | Número do chunk `_001`                      | Grupo 5: chunk (3 dígitos)          |
| `\.fastq(?:\.gz)?$` | Extensão de arquivo `.fastq` ou `.fastq.gz` | Não capturado (?: é não-capturante) |

Isso analisa convenções de nomenclatura estilo Illumina para extrair metadados automaticamente.

Execute o fluxo de trabalho modificado:

```bash title="Test pattern matching"
nextflow run main.nf
```

??? success "Saída do comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [clever_pauling] DSL2 - revision: 605d2058b4

    [[id:sample_001, organism:human, tissue:liver, depth:30000000, quality:38.5, sample_num:1, lane:001, read:R1, chunk:001, priority:normal], /workspaces/training/side-quests/essential_scripting_patterns/data/sequences/SAMPLE_001_S1_L001_R1_001.fastq]
    [[id:sample_002, organism:mouse, tissue:brain, depth:25000000, quality:35.2, sample_num:2, lane:001, read:R1, chunk:001, priority:normal], /workspaces/training/side-quests/essential_scripting_patterns/data/sequences/SAMPLE_002_S2_L001_R1_001.fastq]
    [[id:sample_003, organism:human, tissue:kidney, depth:45000000, quality:42.1, sample_num:3, lane:001, read:R1, chunk:001, priority:high], /workspaces/training/side-quests/essential_scripting_patterns/data/sequences/SAMPLE_003_S3_L001_R1_001.fastq]
    ```

Isso mostra os metadados enriquecidos a partir dos nomes de arquivos.

### 2.2. Geração Dinâmica de Scripts em Processos

Blocos de script de processo são essencialmente strings multi-linha que são passadas para o shell. Você pode usar **lógica condicional** (if/else, operadores ternários) para gerar dinamicamente diferentes strings de script baseadas em características de entrada. Isso é essencial para lidar com tipos de entrada diversos—como leituras single-end vs paired-end—sem duplicar definições de processo.

Vamos adicionar um processo ao nosso fluxo de trabalho que demonstra esse padrão. Abra `modules/fastp.nf` e dê uma olhada:

```groovy title="modules/fastp.nf" linenums="1"
process FASTP {
    container 'community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*_trimmed*.fastq.gz"), emit: reads

    script:
    """
    fastp \\
        --in1 ${reads[0]} \\
        --in2 ${reads[1]} \\
        --out1 ${meta.id}_trimmed_R1.fastq.gz \\
        --out2 ${meta.id}_trimmed_R2.fastq.gz \\
        --json ${meta.id}.fastp.json \\
        --html ${meta.id}.fastp.html \\
        --thread $task.cpus
    """
}
```

O processo recebe arquivos FASTQ como entrada e executa a ferramenta `fastp` para cortar adaptadores e filtrar leituras de baixa qualidade. Infelizmente, a pessoa que escreveu este processo não permitiu as leituras single-end que temos em nosso conjunto de dados de exemplo. Vamos adicioná-lo ao nosso fluxo de trabalho e ver o que acontece:

Primeiro, inclua o módulo na primeira linha do seu fluxo de trabalho `main.nf`:

```groovy title="main.nf" linenums="1"
include { FASTP } from './modules/fastp.nf'
```

Então modifique o bloco `workflow` para conectar o canal `ch_samples` ao processo `FASTP`:

=== "Depois"

    ```groovy title="main.nf" linenums="25" hl_lines="27"
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
    }
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="25" hl_lines="26"
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
                return [sample_meta + file_meta + [priority: priority], file(row.file_path)]
            }
            .view()
    }
    ```

Execute este fluxo de trabalho modificado:

```bash
nextflow run main.nf
```

??? failure "Saída do comando"

    ```console
    ERROR ~ Error executing process > 'FASTP (3)'

    Caused by:
      Process `FASTP (3)` terminated with an error exit status (255)


    Command executed:

      fastp \
          --in1 SAMPLE_003_S3_L001_R1_001.fastq \
          --in2 null \
          --out1 sample_003_trimmed_R1.fastq.gz \
          --out2 sample_003_trimmed_R2.fastq.gz \
          --json sample_003.fastp.json \
          --html sample_003.fastp.html \
          --thread 2

    Command exit status:
      255

    Command output:
      (empty)
    ```

Você pode ver que o processo está tentando executar `fastp` com um valor `null` para o segundo arquivo de entrada, o que está causando falha. Isso ocorre porque nosso conjunto de dados contém leituras single-end, mas o processo está codificado para esperar leituras paired-end (dois arquivos de entrada por vez).

Corrija isso adicionando lógica condicional ao bloco `script:` do processo `FASTP`. Uma instrução if/else verifica a contagem de arquivos de leitura e ajusta o comando de acordo.

=== "Depois"

    ```groovy title="main.nf" linenums="10" hl_lines="3-27"
        script:
        // Detecção simples de single-end vs paired-end
        def is_single = reads instanceof List ? reads.size() == 1 : true

        if (is_single) {
            def input_file = reads instanceof List ? reads[0] : reads
            """
            fastp \\
                --in1 ${input_file} \\
                --out1 ${meta.id}_trimmed.fastq.gz \\
                --json ${meta.id}.fastp.json \\
                --html ${meta.id}.fastp.html \\
                --thread $task.cpus
            """
        } else {
            """
            fastp \\
                --in1 ${reads[0]} \\
                --in2 ${reads[1]} \\
                --out1 ${meta.id}_trimmed_R1.fastq.gz \\
                --out2 ${meta.id}_trimmed_R2.fastq.gz \\
                --json ${meta.id}.fastp.json \\
                --html ${meta.id}.fastp.html \\
                --thread $task.cpus
            """
        }
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="10" hl_lines="2-11"
            script:
            """
            fastp \\
                --in1 ${reads[0]} \\
                --in2 ${reads[1]} \\
                --out1 ${meta.id}_trimmed_R1.fastq.gz \\
                --out2 ${meta.id}_trimmed_R2.fastq.gz \\
                --json ${meta.id}.fastp.json \\
                --html ${meta.id}.fastp.html \\
                --thread $task.cpus
            """
        }
    ```

Agora o fluxo de trabalho pode lidar com leituras single-end e paired-end graciosamente. A lógica condicional verifica o número de arquivos de entrada e constrói o comando apropriado para `fastp`. Vamos ver se funciona:

```bash
nextflow run main.nf
```

??? success "Saída do comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [adoring_rosalind] DSL2 - revision: 04b1cd93e9

    executor >  local (3)
    [31/a8ad4d] process > FASTP (3) [100%] 3 of 3 ✔
    ```

Parece bom! Se verificarmos os comandos reais que foram executados (personalize para seu hash de tarefa):

```console title="Check commands executed"
cat work/31/a8ad4d95749e685a6d842d3007957f/.command.sh
```

Podemos ver que o Nextflow escolheu corretamente o comando certo para leituras single-end:

```bash title=".command.sh"
#!/bin/bash -ue
fastp \
    --in1 SAMPLE_003_S3_L001_R1_001.fastq \
    --out1 sample_003_trimmed.fastq.gz \
    --json sample_003.fastp.json \
    --html sample_003.fastp.html \
    --thread 2
```

Outro uso comum de lógica de script dinâmica pode ser visto em [o módulo Nextflow for Science Genomics](../../nf4science/genomics/02_joint_calling). Nesse módulo, o processo GATK sendo chamado pode receber múltiplos arquivos de entrada, mas cada um deve ser prefixado com `-V` para formar uma linha de comando correta. O processo usa script para transformar uma coleção de arquivos de entrada (`all_gvcfs`) nos argumentos de comando corretos:

```groovy title="command line manipulation for GATK" linenums="1" hl_lines="2 5"
    script:
    def gvcfs_line = all_gvcfs.collect { gvcf -> "-V ${gvcf}" }.join(' ')
    """
    gatk GenomicsDBImport \
        ${gvcfs_line} \
        -L ${interval_list} \
        --genomicsdb-workspace-path ${cohort_name}_gdb
    """
```

Esses padrões de usar script em blocos de script de processo são extremamente poderosos e podem ser aplicados em muitos cenários - desde lidar com tipos de entrada variáveis até construir argumentos complexos de linha de comando a partir de coleções de arquivos, tornando seus processos verdadeiramente adaptáveis aos requisitos diversos de dados do mundo real.

### 2.3. Interpolação de Variáveis: Variáveis Nextflow e Shell

Scripts de processo misturam variáveis Nextflow, variáveis shell e substituições de comando, cada uma com sintaxe de interpolação diferente. Usar a sintaxe errada causa erros. Vamos explorar isso com um processo que cria um relatório de processamento.

Dê uma olhada no arquivo de módulo `modules/generate_report.nf`:

```groovy title="modules/generate_report.nf" linenums="1"
process GENERATE_REPORT {

    publishDir 'results/reports', mode: 'copy'

    input:
    tuple val(meta), path(reads)

    output:
    path "${meta.id}_report.txt"

    script:
    """
    echo "Processing ${reads}" > ${meta.id}_report.txt
    echo "Sample: ${meta.id}" >> ${meta.id}_report.txt
    """
}
```

Este processo escreve um relatório simples com o ID da amostra e nome do arquivo. Agora vamos executá-lo para ver o que acontece quando precisamos misturar diferentes tipos de variáveis.

Inclua o processo em seu `main.nf` e adicione-o ao fluxo de trabalho:

=== "Depois"

    ```groovy title="main.nf" linenums="1" hl_lines="2 30"
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

=== "Antes"

    ```groovy title="main.nf" linenums="1" hl_lines="1 10-29"
    include { FASTP } from './modules/fastp.nf'

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
    }
    ```

Agora execute o fluxo de trabalho e verifique os relatórios gerados em `results/reports/`. Eles devem conter informações básicas sobre cada amostra.

<!-- TODO: add the run command -->

??? success "Saída do comando"

    ```console
    <!-- TODO: output -->
    ```

Mas e se quisermos adicionar informações sobre quando e onde o processamento ocorreu? Vamos modificar o processo para usar variáveis **shell** e um pouco de substituição de comando para incluir o usuário atual, hostname e data no relatório:

=== "Depois"

    ```groovy title="modules/generate_report.nf" linenums="10" hl_lines="5-7"
        script:
        """
        echo "Processing ${reads}" > ${meta.id}_report.txt
        echo "Sample: ${meta.id}" >> ${meta.id}_report.txt
        echo "Processed by: ${USER}" >> ${meta.id}_report.txt
        echo "Hostname: $(hostname)" >> ${meta.id}_report.txt
        echo "Date: $(date)" >> ${meta.id}_report.txt
        """
    ```

=== "Antes"

    ```groovy title="modules/generate_report.nf" linenums="10"
        script:
        """
        echo "Processing ${reads}" > ${meta.id}_report.txt
        echo "Sample: ${meta.id}" >> ${meta.id}_report.txt
        """
    ```

Se você executar isso, notará um erro - Nextflow tenta interpretar `${USER}` como uma variável Nextflow que não existe.

??? failure "Saída do comando"

    ```console
    Error modules/generate_report.nf:15:27: `USER` is not defined
    │  15 |     echo "Processed by: ${USER}" >> ${meta.id}_report.txt
    ╰     |                           ^^^^

    ERROR ~ Script compilation failed
    ```

Precisamos escapá-la para que o Bash possa lidar com ela.

Corrija isso escapando as variáveis shell e substituições de comando com uma barra invertida (`\`):

=== "Depois"

    ```groovy title="modules/generate_report.nf" linenums="10" hl_lines="5-7"
        script:
        """
        echo "Processing ${reads}" > ${meta.id}_report.txt
        echo "Sample: ${meta.id}" >> ${meta.id}_report.txt
        echo "Processed by: \${USER}" >> ${meta.id}_report.txt
        echo "Hostname: \$(hostname)" >> ${meta.id}_report.txt
        echo "Date: \$(date)" >> ${meta.id}_report.txt
        """
    ```

=== "Antes"

    ```groovy title="modules/generate_report.nf" linenums="10"
        script:
        """
        echo "Processing ${reads}" > ${meta.id}_report.txt
        echo "Sample: ${meta.id}" >> ${meta.id}_report.txt
        echo "Processed by: ${USER}" >> ${meta.id}_report.txt
        echo "Hostname: $(hostname)" >> ${meta.id}_report.txt
        echo "Date: $(date)" >> ${meta.id}_report.txt
        """
    ```

Agora funciona! A barra invertida (`\`) diz ao Nextflow "não interprete isso, passe para o Bash."

### Conclusão

Nesta seção, você aprendeu técnicas de **processamento de strings**:

- **Expressões regulares para análise de arquivos**: Usando o operador `=~` e padrões regex (`~/padrão/`) para extrair metadados de convenções complexas de nomenclatura de arquivos
- **Geração dinâmica de scripts**: Usando lógica condicional (if/else, operadores ternários) para gerar diferentes strings de script baseadas em características de entrada
- **Interpolação de variáveis**: Entendendo quando o Nextflow interpreta strings vs quando o shell interpreta
  - `${var}` - Variáveis Nextflow (interpoladas pelo Nextflow em tempo de compilação do fluxo de trabalho)
  - `\${var}` - Variáveis de ambiente shell (escapadas, passadas para bash em tempo de execução)
  - `\$(cmd)` - Substituição de comando shell (escapada, executada pelo bash em tempo de execução)

Esses padrões de processamento e geração de strings são essenciais para lidar com os diversos formatos de arquivo e convenções de nomenclatura que você encontrará em fluxos de trabalho reais de bioinformática.

---

## 3. Criando Funções Reutilizáveis

Lógica complexa de fluxo de trabalho inline em operadores de canal ou definições de processo reduz legibilidade e manutenibilidade. **Funções** permitem extrair essa lógica em componentes nomeados e reutilizáveis.

Nossa operação map cresceu longa e complexa. Vamos extraí-la em uma função reutilizável usando a palavra-chave `def`.

Para ilustrar como isso se parece com nosso fluxo de trabalho existente, faça a modificação abaixo, usando `def` para definir uma função reutilizável chamada `separateMetadata`:

=== "Depois"

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

=== "Antes"

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

Ao extrair essa lógica em uma função, reduzimos a lógica real do fluxo de trabalho para algo muito mais limpo:

```groovy title="minimal workflow"
    ch_samples = channel.fromPath("./data/samples.csv")
        .splitCsv(header: true)
        .map{ row -> separateMetadata(row) }

    ch_fastp = FASTP(ch_samples)
    GENERATE_REPORT(ch_samples)
```

Isso torna a lógica do fluxo de trabalho muito mais fácil de ler e entender rapidamente. A função `separateMetadata` encapsula toda a lógica complexa para analisar e enriquecer metadados, tornando-a reutilizável e testável.

Execute o fluxo de trabalho para ter certeza de que ainda funciona:

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

A saída deve mostrar ambos os processos completando com sucesso. O fluxo de trabalho agora está muito mais limpo e fácil de manter, com toda a lógica complexa de processamento de metadados encapsulada na função `separateMetadata`.

### Conclusão

Nesta seção, você aprendeu **criação de funções**:

- **Definindo funções com `def`**: A palavra-chave para criar funções nomeadas (como `def` em Python ou `function` em JavaScript)
- **Escopo de função**: Funções definidas no nível do script são acessíveis em todo o seu fluxo de trabalho Nextflow
- **Valores de retorno**: Funções retornam automaticamente a última expressão, ou usam `return` explícito
- **Código mais limpo**: Extrair lógica complexa em funções é uma prática fundamental de engenharia de software em qualquer linguagem

A seguir, exploraremos como usar closures em diretivas de processo para alocação dinâmica de recursos.

---

## 4. Diretivas de Recursos Dinâmicas com Closures

Até agora usamos script no bloco `script` de processos. Mas **closures** (introduzidas na Seção 1.1) também são incrivelmente úteis em diretivas de processo, especialmente para alocação dinâmica de recursos. Vamos adicionar diretivas de recursos ao nosso processo FASTP que se adaptam com base nas características da amostra.

### 4.1. Alocação de recursos específica por amostra

Atualmente, nosso processo FASTP usa recursos padrão. Vamos torná-lo mais inteligente alocando mais CPUs para amostras de alta profundidade. Edite `modules/fastp.nf` para incluir uma diretiva `cpus` dinâmica e uma diretiva `memory` estática:

=== "Depois"

    ```groovy title="modules/fastp.nf" linenums="1" hl_lines="4-5"
    process FASTP {
        container 'community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690'

        cpus { meta.depth > 40000000 ? 2 : 1 }
        memory 2.GB

        input:
        tuple val(meta), path(reads)
    ```

=== "Antes"

    ```groovy title="modules/fastp.nf" linenums="1"
    process FASTP {
        container 'community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690'

        input:
        tuple val(meta), path(reads)
    ```

A closure `{ meta.depth > 40000000 ? 2 : 1 }` usa o **operador ternário** (coberto na Seção 1.1) e é avaliada para cada tarefa, permitindo alocação de recursos por amostra. Amostras de alta profundidade (>40M leituras) recebem 2 CPUs, enquanto outras recebem 1 CPU.

!!! note "Acessando Variáveis de Entrada em Diretivas"

    A closure pode acessar quaisquer variáveis de entrada (como `meta` aqui) porque o Nextflow avalia essas closures no contexto de cada execução de tarefa.

Execute o fluxo de trabalho novamente com a opção `-ansi-log false` para facilitar a visualização dos hashes de tarefa.

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

Você pode verificar o comando `docker` exato que foi executado para ver a alocação de CPU para qualquer tarefa dada:

```console title="Check docker command"
cat work/48/6db0c9e9d8aa65e4bb4936cd3bd59e/.command.run | grep "docker run"
```

Você deve ver algo como:

```bash title="docker command"
    docker run -i --cpu-shares 4096 --memory 2048m -e "NXF_TASK_WORKDIR" -v /workspaces/training/side-quests/essential_scripting_patterns:/workspaces/training/side-quests/essential_scripting_patterns -w "$NXF_TASK_WORKDIR" --name $NXF_BOXID community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690 /bin/bash -ue /workspaces/training/side-quests/essential_scripting_patterns/work/48/6db0c9e9d8aa65e4bb4936cd3bd59e/.command.sh
```

Neste exemplo escolhemos um exemplo que solicitou 2 CPUs (`--cpu-shares 2048`), porque era uma amostra de alta profundidade, mas você deve ver diferentes alocações de CPU dependendo da profundidade da amostra. Tente isso para as outras tarefas também.

### 4.2. Estratégias de retry

Outro padrão poderoso é usar `task.attempt` para estratégias de retry. Para mostrar por que isso é útil, vamos começar reduzindo a alocação de memória para FASTP para menos do que ele precisa. Mude a diretiva `memory` em `modules/fastp.nf` para `1.GB`:

=== "Depois"

    ```groovy title="modules/fastp.nf" linenums="1" hl_lines="5"
    process FASTP {
        container 'community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690'

        cpus { meta.depth > 40000000 ? 4 : 2 }
        memory 1.GB

        input:
        tuple val(meta), path(reads)
    ```

=== "Antes"

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

Isso indica que o processo foi morto por exceder limites de memória.

Este é um cenário muito comum em fluxos de trabalho do mundo real - às vezes você simplesmente não sabe quanta memória uma tarefa precisará até executá-la.

Para tornar nosso fluxo de trabalho mais robusto, podemos implementar uma estratégia de retry que aumenta a alocação de memória em cada tentativa, mais uma vez usando uma closure Groovy. Modifique a diretiva `memory` para multiplicar a memória base por `task.attempt`, e adicione diretivas `errorStrategy 'retry'` e `maxRetries 2`:

=== "Depois"

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

=== "Antes"

    ```groovy title="modules/fastp.nf" linenums="1" hl_lines="5"
    process FASTP {
        container 'community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690'

        cpus { meta.depth > 40000000 ? 4 : 2 }
        memory 2.GB

        input:
        tuple val(meta), path(reads)
    ```

Agora se o processo falhar devido a memória insuficiente, o Nextflow tentará novamente com mais memória:

- Primeira tentativa: 1 GB (task.attempt = 1)
- Segunda tentativa: 2.GB (task.attempt = 2)

... e assim por diante, até o limite `maxRetries`.

### Conclusão

Diretivas dinâmicas com closures permitem:

- Alocar recursos baseados em características de entrada
- Implementar estratégias automáticas de retry com recursos crescentes
- Combinar múltiplos fatores (metadados, número de tentativa, prioridades)
- Usar lógica condicional para cálculos complexos de recursos

Isso torna seus fluxos de trabalho tanto mais eficientes (não sobre-alocando) quanto mais robustos (retry automático com mais recursos).

---

## 5. Lógica Condicional e Controle de Processo

Anteriormente, usamos `.map()` com script para transformar dados de canal. Agora usaremos lógica condicional para controlar quais processos executam baseados em dados—essencial para fluxos de trabalho flexíveis adaptando-se a diferentes tipos de amostra.

Os [operadores de fluxo de dados](https://www.nextflow.io/docs/latest/reference/operator.html) do Nextflow recebem closures avaliadas em tempo de execução, habilitando lógica condicional para dirigir decisões de fluxo de trabalho baseadas no conteúdo do canal.

### 5.1. Roteamento com `.branch()`

Por exemplo, vamos fingir que nossas amostras de sequenciamento precisam ser cortadas com FASTP apenas se forem amostras humanas com cobertura acima de um certo limiar. Amostras de camundongo ou amostras de baixa cobertura devem ser executadas com Trimgalore (este é um exemplo artificial, mas ilustra o ponto).

Fornecemos um processo Trimgalore simples em `modules/trimgalore.nf`, dê uma olhada se quiser, mas os detalhes não são importantes para este exercício. O ponto chave é que queremos rotear amostras baseadas em seus metadados.

Inclua o novo módulo de `modules/trimgalore.nf`:

=== "Depois"

    ```groovy title="main.nf" linenums="1" hl_lines="2"
    include { FASTP } from './modules/fastp.nf'
    include { TRIMGALORE } from './modules/trimgalore.nf'
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="1"
    include { FASTP } from './modules/fastp.nf'
    ```

... e então modifique seu fluxo de trabalho `main.nf` para ramificar amostras baseadas em seus metadados e roteá-las através do processo de corte apropriado, assim:

=== "Depois"

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

=== "Antes"

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

Aqui, usamos expressões condicionais pequenas mas poderosas dentro do operador `.branch{}` para rotear amostras baseadas em seus metadados. Amostras humanas com alta cobertura passam por `FASTP`, enquanto todas as outras amostras passam por `TRIMGALORE`.

### 5.2. Usando `.filter()` com Truthiness

Outro padrão poderoso para controlar a execução do fluxo de trabalho é o operador `.filter()`, que usa uma closure para determinar quais itens devem continuar pelo pipeline. Dentro da closure de filtro, você escreverá **expressões booleanas** que decidem quais itens passam.

Nextflow (como muitas linguagens dinâmicas) tem um conceito de **"truthiness"** que determina quais valores avaliam para `true` ou `false` em contextos booleanos:

- **Truthy**: Valores não-null, strings não-vazias, números não-zero, coleções não-vazias
- **Falsy**: `null`, strings vazias `""`, zero `0`, coleções vazias `[]` ou `[:]`, `false`

Isso significa que `meta.id` sozinho (sem `!= null` explícito) verifica se o ID existe e não está vazio. Vamos usar isso para filtrar amostras que não atendem nossos requisitos de qualidade.

Adicione o seguinte antes da operação branch:

=== "Depois"

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

=== "Antes"

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

- `meta.id && meta.organism` verifica que ambos os campos existem e não estão vazios (usando truthiness)
- `meta.depth >= 25000000` garante profundidade de sequenciamento suficiente com uma comparação explícita

!!! note "Truthiness na Prática"

    A expressão `meta.id && meta.organism` é mais concisa do que escrever:
    ```groovy
    meta.id != null && meta.id != '' && meta.organism != null && meta.organism != ''
    ```

    Isso torna a lógica de filtragem muito mais limpa e fácil de ler.

### Conclusão

Nesta seção, você aprendeu a usar lógica condicional para controlar a execução do fluxo de trabalho usando as interfaces de closure de operadores Nextflow como `.branch{}` e `.filter{}`, aproveitando truthiness para escrever expressões condicionais concisas.

Nosso pipeline agora roteia amostras inteligentemente através de processos apropriados, mas fluxos de trabalho de produção precisam lidar com dados inválidos graciosamente. Vamos tornar nosso fluxo de trabalho robusto contra valores ausentes ou null.

---

## 6. Operadores de Navegação Segura e Elvis

Nossa função `separateMetadata` atualmente assume que todos os campos CSV estão presentes e válidos. Mas o que acontece com dados incompletos? Vamos descobrir.

### 6.1. O Problema: Acessando Propriedades Que Não Existem

Digamos que queremos adicionar suporte para informações opcionais de execução de sequenciamento. Em alguns laboratórios, amostras podem ter um campo adicional para o ID de execução de sequenciamento ou número de lote, mas nosso CSV atual não tem essa coluna. Vamos tentar acessá-la de qualquer forma.

Modifique a função `separateMetadata` para incluir um campo run_id:

=== "Depois"

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

=== "Antes"

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

Isso falha com uma NullPointerException.

O problema é que `row.run_id` retorna `null` porque a coluna `run_id` não existe em nosso CSV. Quando tentamos chamar `.toUpperCase()` em `null`, ele falha. É aqui que o operador de navegação segura salva o dia.

### 6.2. Operador de Navegação Segura (`?.`)

O operador de navegação segura (`?.`) retorna `null` em vez de lançar uma exceção quando chamado em um valor `null`. Se o objeto antes de `?.` for `null`, a expressão inteira avalia para `null` sem executar o método.

Atualize a função para usar navegação segura:

=== "Depois"

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

=== "Antes"

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

Sem falha! O fluxo de trabalho agora lida com o campo ausente graciosamente. Quando `row.run_id` é `null`, o operador `?.` previne a chamada `.toUpperCase()`, e `run_id` se torna `null` em vez de causar uma exceção.

### 6.3. Operador Elvis (`?:`) para Padrões

O operador Elvis (`?:`) fornece valores padrão quando o lado esquerdo é "falsy" (como explicado anteriormente). É nomeado após Elvis Presley porque `?:` parece com seu famoso cabelo e olhos quando visto de lado!

Agora que estamos usando navegação segura, `run_id` será `null` para amostras sem esse campo. Vamos usar o operador Elvis para fornecer um valor padrão e adicioná-lo ao nosso mapa `sample_meta`:

=== "Depois"

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

=== "Antes"

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

Também adicione um operador `view()` no fluxo de trabalho para ver os resultados:

=== "Depois"

    ```groovy title="main.nf" linenums="30" hl_lines="4"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map{ row -> separateMetadata(row) }
            .view()
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="30"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map{ row -> separateMetadata(row) }
    ```

e execute o fluxo de trabalho:

```bash
nextflow run main.nf
```

??? success "Saída do comando"

    ```console
    [[id:sample_001, organism:human, tissue:liver, depth:30000000, quality:38.5, run:UNSPECIFIED, sample_num:1, lane:001, read:R1, chunk:001, priority:normal], /workspaces/training/side-quests/essential_scripting_patterns/data/sequences/SAMPLE_001_S1_L001_R1_001.fastq]
    [[id:sample_002, organism:mouse, tissue:brain, depth:25000000, quality:35.2, run:UNSPECIFIED, sample_num:2, lane:001, read:R1, chunk:001, priority:normal], /workspaces/training/side-quests/essential_scripting_patterns/data/sequences/SAMPLE_002_S2_L001_R1_001.fastq]
    [[id:sample_003, organism:human, tissue:kidney, depth:45000000, quality:42.1, run:UNSPECIFIED, sample_num:3, lane:001, read:R1, chunk:001, priority:high], /workspaces/training/side-quests/essential_scripting_patterns/data/sequences/SAMPLE_003_S3_L001_R1_001.fastq]
    ```

Perfeito! Agora todas as amostras têm um campo `run` com seu ID de execução real (em maiúsculas) ou o valor padrão 'UNSPECIFIED'. A combinação de `?.` e `?:` fornece tanto segurança (sem falhas) quanto padrões sensatos.

Remova o operador `.view()` agora que confirmamos que funciona.

!!! tip "Combinando Navegação Segura e Elvis"

    O padrão `value?.method() ?: 'default'` é comum em fluxos de trabalho de produção:

    - `value?.method()` - Chama método com segurança, retorna `null` se `value` for `null`
    - `?: 'default'` - Fornece fallback se o resultado for `null`

    Este padrão lida com dados ausentes/incompletos graciosamente.

Use esses operadores consistentemente em funções, closures de operadores (`.map{}`, `.filter{}`), scripts de processo e arquivos de configuração. Eles previnem falhas ao lidar com dados do mundo real.

### Conclusão

- **Navegação segura (`?.`)**: Previne falhas em valores null - retorna null em vez de lançar exceção
- **Operador Elvis (`?:`)**: Fornece padrões - `value ?: 'default'`
- **Combinando**: `value?.method() ?: 'default'` é o padrão comum

Esses operadores tornam fluxos de trabalho resilientes a dados incompletos - essencial para trabalho do mundo real.

---

## 7. Validação com `error()` e `log.warn`

Às vezes você precisa parar o fluxo de trabalho imediatamente se parâmetros de entrada forem inválidos. No Nextflow, você pode usar funções integradas como `error()` e `log.warn`, bem como construções padrão de programação como instruções `if` e lógica booleana, para implementar lógica de validação. Vamos adicionar validação ao nosso fluxo de trabalho.

Crie uma função de validação antes do seu bloco workflow, chame-a do workflow, e mude a criação do canal para usar um parâmetro para o caminho do arquivo CSV. Se o parâmetro estiver ausente ou o arquivo não existir, chame `error()` para parar a execução com uma mensagem clara.

=== "Depois"

    ```groovy title="main.nf" linenums="1" hl_lines="5-20 23-24"
    include { FASTP } from './modules/fastp.nf'
    include { TRIMGALORE } from './modules/trimgalore.nf'
    include { GENERATE_REPORT } from './modules/generate_report.nf'

    def validateInputs() {
        // Verificar se parâmetro de entrada foi fornecido
        if (!params.input) {
            error("Input CSV file path not provided. Please specify --input <file.csv>")
        }

        // Verificar se arquivo CSV existe
        if (!file(params.input).exists()) {
            error("Input CSV file not found: ${params.input}")
        }
    }
    ...
    workflow {
        validateInputs()
        ch_samples = channel.fromPath(params.input)
    ```

=== "Antes"

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

O fluxo de trabalho para imediatamente com uma mensagem de erro clara em vez de falhar misteriosamente depois

Agora execute com um arquivo inexistente:

```bash
nextflow run main.nf --input ./data/nonexistent.csv
```

??? failure "Saída do comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [cranky_gates] DSL2 - revision: 26839ae3eb

    Input CSV file not found: ./data/nonexistent.csv
    ```

Finalmente, execute com o arquivo correto:

```bash
nextflow run main.nf --input ./data/samples.csv
```

??? success "Saída do comando"

    ```console
    <!-- TODO: output -->
    ```

Desta vez executa com sucesso.

Você também pode adicionar validação dentro da função `separateMetadata`. Vamos usar o não-fatal `log.warn` para emitir avisos para amostras com baixa profundidade de sequenciamento, mas ainda permitir que o fluxo de trabalho continue:

=== "Depois"

    ```groovy title="main.nf" linenums="1" hl_lines="3-6"
        def priority = sample_meta.quality > 40 ? 'high' : 'normal'

        // Validar que dados fazem sentido
        if (sample_meta.depth < 30000000) {
            log.warn "Low sequencing depth for ${sample_meta.id}: ${sample_meta.depth}"
        }

        return tuple(sample_meta + file_meta + [priority: priority], fastq_path)
    }
    ```

=== "Antes"

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

- **`error()`**: Para fluxo de trabalho imediatamente com mensagem clara
- **`log.warn`**: Emite avisos sem parar o fluxo de trabalho
- **Validação precoce**: Verificar entradas antes de processar para falhar rápido com erros úteis
- **Funções de validação**: Criar lógica de validação reutilizável que pode ser chamada no início do fluxo de trabalho

Validação adequada torna fluxos de trabalho mais robustos e amigáveis ao usuário capturando problemas cedo com mensagens de erro claras.

---

## 8. Manipuladores de Eventos de Fluxo de Trabalho

Até agora, temos escrito código em nossos scripts de fluxo de trabalho e definições de processo. Mas há mais um recurso importante que você deve conhecer: manipuladores de eventos de fluxo de trabalho.

Manipuladores de eventos são closures que executam em pontos específicos no ciclo de vida do seu fluxo de trabalho. Eles são perfeitos para adicionar logging, notificações ou operações de limpeza. Esses manipuladores devem ser definidos em seu script de fluxo de trabalho junto com sua definição de fluxo de trabalho.

### 8.1. O Manipulador `onComplete`

O manipulador de eventos mais comumente usado é `onComplete`, que executa quando seu fluxo de trabalho termina (seja com sucesso ou falha). Vamos adicionar um para resumir os resultados do nosso pipeline.

Adicione o manipulador de eventos ao seu arquivo `main.nf`, dentro da sua definição de fluxo de trabalho:

=== "Depois"

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

=== "Antes"

    ```groovy title="main.nf" linenums="66" hl_lines="4"
        ch_fastp = FASTP(trim_branches.fastp)
        ch_trimgalore = TRIMGALORE(trim_branches.trimgalore)
        GENERATE_REPORT(ch_samples)
    }
    ```

Esta closure executa quando o fluxo de trabalho completa. Dentro, você tem acesso ao objeto `workflow` que fornece propriedades úteis sobre a execução.

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

=== "Depois"

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

=== "Antes"

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

Agora obtemos um resumo ainda mais informativo, incluindo uma mensagem de sucesso/falha e o diretório de saída se especificado:

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
    // ... seu código de fluxo de trabalho ...

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

        // Escrever em um arquivo de log
        def log_file = file("${workflow.launchDir}/pipeline_summary.txt")
        log_file.text = summary
    }
}
```

### 8.2. O Manipulador `onError`

Além de `onComplete`, há um outro manipulador de eventos que você pode usar: `onError`, que executa apenas se o fluxo de trabalho falhar:

```groovy title="main.nf - onError handler"
workflow {
    // ... seu código de fluxo de trabalho ...

    workflow.onError = {
        println "="* 50
        println "Pipeline execution failed!"
        println "Error message: ${workflow.errorMessage}"
        println "="* 50

        // Escrever log de erro detalhado
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

Você pode usar múltiplos manipuladores juntos em seu script de fluxo de trabalho:

```groovy title="main.nf - Combined handlers"
workflow {
    // ... seu código de fluxo de trabalho ...

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

- **Closures de manipuladores de eventos**: Closures em seu script de fluxo de trabalho que executam em diferentes pontos do ciclo de vida
- **Manipulador `onComplete`**: Para resumos de execução e relatórios de resultados
- **Manipulador `onError`**: Para tratamento de erros e logging de falhas
- **Propriedades do objeto workflow**: Acessando `workflow.success`, `workflow.duration`, `workflow.errorMessage`, etc.

Manipuladores de eventos mostram como você pode usar o poder completo da linguagem Nextflow dentro de seus scripts de fluxo de trabalho para adicionar capacidades sofisticadas de logging e notificação.

---

## Resumo

Parabéns, você conseguiu!

Ao longo desta missão secundária, você construiu um pipeline abrangente de processamento de amostras que evoluiu de manipulação básica de metadados para um fluxo de trabalho sofisticado e pronto para produção.
Cada seção construiu sobre a anterior, demonstrando como construções de programação transformam fluxos de trabalho simples em sistemas poderosos de processamento de dados, com os seguintes benefícios:

- **Código mais claro**: Entender fluxo de dados vs script ajuda você a escrever fluxos de trabalho mais organizados
- **Manipulação robusta**: Navegação segura e operadores Elvis tornam fluxos de trabalho resilientes a dados ausentes
- **Processamento flexível**: Lógica condicional permite que seus fluxos de trabalho processem diferentes tipos de amostra apropriadamente
- **Recursos adaptativos**: Diretivas dinâmicas otimizam uso de recursos baseado em características de entrada

Esta progressão espelha a evolução do mundo real de pipelines de bioinformática, de protótipos de pesquisa lidando com poucas amostras a sistemas de produção processando milhares de amostras através de laboratórios e instituições.
Cada desafio que você resolveu e padrão que aprendeu reflete problemas reais que desenvolvedores enfrentam ao escalar fluxos de trabalho Nextflow.

Aplicar esses padrões em seu próprio trabalho permitirá que você construa fluxos de trabalho robustos e prontos para produção.

### Padrões chave

1.  **Fluxo de Dados vs Script:** Você aprendeu a distinguir entre operações de fluxo de dados (orquestração de canal) e script (código que manipula dados), incluindo as diferenças cruciais entre operações em diferentes tipos como `collect` em Canal vs Lista.

    - Fluxo de dados: orquestração de canal

    ```groovy
    channel.fromPath('*.fastq').splitCsv(header: true)
    ```

    - Script: processamento de dados em coleções

    ```groovy
    sample_data.collect { it.toUpperCase() }
    ```

2.  **Processamento Avançado de Strings**: Você dominou expressões regulares para analisar nomes de arquivos, geração dinâmica de scripts em processos e interpolação de variáveis (Nextflow vs Bash vs Shell).

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

3.  **Criando Funções Reutilizáveis**: Você aprendeu a extrair lógica complexa em funções nomeadas que podem ser chamadas de operadores de canal, tornando fluxos de trabalho mais legíveis e mantíveis.

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

4.  **Diretivas de Recursos Dinâmicas com Closures**: Você explorou usar closures em diretivas de processo para alocação adaptativa de recursos baseada em características de entrada.

    - Closures nomeadas e composição

    ```groovy
    def enrichData = normalizeId >> addQualityCategory >> addFlags
    def processor = generalFunction.curry(fixedParam)
    ```

    - Closures com acesso a escopo

    ```groovy
    def collectStats = { data -> stats.count++; return data }
    ```

5.  **Lógica Condicional e Controle de Processo**: Você adicionou roteamento inteligente usando operadores `.branch()` e `.filter()`, aproveitando truthiness para expressões condicionais concisas.

    - Use `.branch()` para rotear dados através de diferentes ramos de fluxo de trabalho

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

    - Use `filter()` para subconjuntar dados com 'truthiness'

    ```groovy
    ch_valid_samples = ch_samples
        .filter { meta, reads ->
            meta.id && meta.organism && meta.depth >= 25000000
        }
    ```

6.  **Operadores de Navegação Segura e Elvis**: Você tornou o pipeline robusto contra dados ausentes usando `?.` para acesso null-safe a propriedades e `?:` para fornecer valores padrão.

    ```groovy
    def id = data?.sample?.id ?: 'unknown'
    ```

7.  **Validação com error() e log.warn**: Você aprendeu a validar entradas cedo e falhar rápido com mensagens de erro claras.

    ```groovy
    try {
        def errors = validateSample(sample)
        if (errors) throw new RuntimeException("Invalid: ${errors.join(', ')}")
    } catch (Exception e) {
        println "Error: ${e.message}"
    }
    ```

8.  **Manipuladores de Eventos de Configuração**: Você aprendeu a usar manipuladores de eventos de fluxo de trabalho (`onComplete` e `onError`) para logging, notificações e gerenciamento de ciclo de vida.

    - Usando `onComplete` para log e notificação

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

    - Usando `onError` para tomar ação especificamente em caso de falha

    ```groovy
    workflow.onError = {
        // Escrever log de erro detalhado
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

Você se beneficiará de praticar e expandir suas habilidades para:

- Escrever fluxos de trabalho mais limpos com separação adequada entre fluxo de dados e script
- Dominar interpolação de variáveis para evitar armadilhas comuns com variáveis Nextflow, Bash e shell
- Usar diretivas de recursos dinâmicas para fluxos de trabalho eficientes e adaptativos
- Transformar coleções de arquivos em argumentos de linha de comando formatados adequadamente
- Lidar com diferentes convenções de nomenclatura de arquivos e formatos de entrada graciosamente usando regex e processamento de strings
- Construir código reutilizável e mantível usando padrões avançados de closure e programação funcional
- Processar e organizar conjuntos de dados complexos usando operações de coleção
- Adicionar validação, tratamento de erros e logging para tornar seus fluxos de trabalho prontos para produção
- Implementar gerenciamento de ciclo de vida de fluxo de trabalho com manipuladores de eventos

---

## O que vem a seguir?

Retorne ao [menu de Missões Secundárias](./index.md) ou clique no botão no canto inferior direito da página para passar para o próximo tópico da lista.
