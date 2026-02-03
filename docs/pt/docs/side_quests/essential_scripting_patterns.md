# Padrões Essenciais de Script em Nextflow

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tradução assistida por IA - [saiba mais e sugira melhorias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Nextflow é uma linguagem de programação que roda na Máquina Virtual Java. Embora Nextflow seja construído sobre [Groovy](http://groovy-lang.org/) e compartilhe muito de sua sintaxe, Nextflow é mais do que apenas "Groovy com extensões" -- é uma linguagem autônoma com uma [sintaxe](https://nextflow.io/docs/latest/reference/syntax.html) e [biblioteca padrão](https://nextflow.io/docs/latest/reference/stdlib.html) totalmente especificadas.

Você pode escrever muito código Nextflow sem ir além da sintaxe básica para variáveis, mapas e listas. A maioria dos tutoriais de Nextflow foca em orquestração de fluxo de trabalho (canais, processos e fluxo de dados), e você pode ir surpreendentemente longe com apenas isso.

No entanto, quando você precisa manipular dados, analisar nomes de arquivo complexos, implementar lógica condicional ou construir fluxos de trabalho robustos para produção, ajuda pensar sobre dois aspectos distintos do seu código: **fluxo de dados** (canais, operadores, processos e fluxos de trabalho) e **scripting** (o código dentro de closures, funções e scripts de processo). Embora essa distinção seja um tanto arbitrária—é tudo código Nextflow—ela fornece um modelo mental útil para entender quando você está orquestrando seu pipeline versus quando você está manipulando dados. Dominar ambos melhora dramaticamente sua capacidade de escrever fluxos de trabalho claros e de fácil manutenção.

### Objetivos de aprendizagem

Esta missão secundária leva você em uma jornada prática desde conceitos básicos até padrões prontos para produção.
Vamos transformar um fluxo de trabalho simples de leitura de CSV em um pipeline sofisticado de bioinformática, evoluindo-o passo a passo através de desafios realistas:

- **Entendendo limites:** Distinguir entre operações de fluxo de dados e scripting, e entender como eles trabalham juntos
- **Manipulação de dados:** Extrair, transformar e subconjuntos de mapas e coleções usando operadores poderosos
- **Processamento de strings:** Analisar esquemas complexos de nomenclatura de arquivos com padrões regex e dominar interpolação de variáveis
- **Funções reutilizáveis:** Extrair lógica complexa em funções nomeadas para fluxos de trabalho mais limpos e de fácil manutenção
- **Lógica dinâmica:** Construir processos que se adaptam a diferentes tipos de entrada e usar closures para alocação dinâmica de recursos
- **Roteamento condicional:** Rotear inteligentemente amostras através de diferentes processos com base em suas características de metadados
- **Operações seguras:** Lidar graciosamente com dados ausentes usando operadores seguros contra null e validar entradas com mensagens de erro claras
- **Handlers baseados em configuração:** Usar handlers de eventos de fluxo de trabalho para logging, notificações e gerenciamento de ciclo de vida

### Pré-requisitos

Antes de iniciar esta missão secundária, você deve:

- Ter concluído o tutorial [Hello Nextflow](../hello_nextflow/README.md) ou curso equivalente para iniciantes.
- Estar confortável usando conceitos e mecanismos básicos de Nextflow (processos, canais, operadores, trabalhar com arquivos, metadados)
- Ter familiaridade básica com construções comuns de programação (variáveis, mapas, listas)

Este tutorial explicará conceitos de programação conforme os encontramos, então você não precisa de experiência extensa em programação.
Começaremos com conceitos fundamentais e construiremos até padrões avançados.

---

## 0. Começando

#### Abrir o codespace de treinamento

Se você ainda não o fez, certifique-se de abrir o ambiente de treinamento conforme descrito em [Configuração do Ambiente](../envsetup/index.md).

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### Mover para o diretório do projeto

Vamos mover para o diretório onde os arquivos para este tutorial estão localizados.

```bash
cd side-quests/essential_scripting_patterns
```

#### Revisar os materiais

Você encontrará um arquivo de fluxo de trabalho principal e um diretório `data` contendo arquivos de dados de exemplo.

```console title="Conteúdo do diretório"
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

Nosso CSV de exemplo contém informações sobre amostras biológicas que precisam de processamento diferente com base em suas características:

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

## 1. Fluxo de Dados vs Scripting: Entendendo os Limites

### 1.1. Identificando o Que É o Quê

Ao escrever fluxos de trabalho Nextflow, é importante distinguir entre **fluxo de dados** (como os dados se movem através de canais e processos) e **scripting** (o código que manipula dados e toma decisões). Vamos construir um fluxo de trabalho demonstrando como eles trabalham juntos.

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

Agora vamos adicionar scripting para transformar os dados, usando o operador `.map()` que você provavelmente já conhece. Este operador recebe uma 'closure' onde podemos escrever código para transformar cada item.

!!! note

    Uma **closure** é um bloco de código que pode ser passado adiante e executado depois. Pense nela como uma função que você define inline. Closures são escritas com chaves `{ }` e podem receber parâmetros. Elas são fundamentais para como os operadores Nextflow funcionam e, se você tem escrito Nextflow há algum tempo, pode já ter estado usando-as sem perceber!

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

Agora vamos escrever lógica de **scripting** dentro de nossa closure para transformar cada linha de dados. É aqui que processamos itens de dados individuais em vez de orquestrar fluxo de dados.

=== "Depois"

    ```groovy title="main.nf" linenums="2" hl_lines="4-12"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row ->
                // Scripting para transformação de dados
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

Usamos métodos de manipulação de string como `.toLowerCase()` e `.replaceAll()` para limpar nossos dados, e métodos de conversão de tipo como `.toInteger()` e `.toDouble()` para converter dados string do CSV nos tipos numéricos apropriados.

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

Agora vamos adicionar mais scripting - desta vez usando um operador ternário para tomar decisões baseadas em valores de dados.

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

!!! Note

    Nunca modifique mapas passados em closures - sempre crie novos usando `+` (por exemplo). No Nextflow, os mesmos dados frequentemente fluem através de múltiplas operações simultaneamente. Modificar um mapa in-place pode causar efeitos colaterais imprevisíveis quando outras operações referenciam esse mesmo objeto. Criar novos mapas garante que cada operação tenha sua própria cópia limpa.

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

#### 1.1.5. Criando Subconjuntos de Mapas com `.subMap()`

Enquanto o operador `+` adiciona chaves a um mapa, às vezes você precisa fazer o oposto - extrair apenas chaves específicas. O método `.subMap()` é perfeito para isso.

Vamos adicionar uma linha para criar uma versão simplificada de nossos metadados que contém apenas campos de identificação:

=== "Depois"

    ```groovy title="main.nf" linenums="2" hl_lines="12-15"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row ->
                // Scripting para transformação de dados
                def sample_meta = [
                    id: row.sample_id.toLowerCase(),
                    organism: row.organism,
                    tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                    depth: row.sequencing_depth.toInteger(),
                    quality: row.quality_score.toDouble()
                ]
                def id_only = sample_meta.subMap(['id', 'organism', 'tissue'])
                println "Apenas campos de ID: ${id_only}"

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
                // Scripting para transformação de dados
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

    Apenas campos de ID: [id:sample_001, organism:human, tissue:liver]
    Apenas campos de ID: [id:sample_002, organism:mouse, tissue:brain]
    Apenas campos de ID: [id:sample_003, organism:human, tissue:kidney]
    [id:sample_001, organism:human, tissue:liver, depth:30000000, quality:38.5, priority:normal]
    [id:sample_002, organism:mouse, tissue:brain, depth:25000000, quality:35.2, priority:normal]
    [id:sample_003, organism:human, tissue:kidney, depth:45000000, quality:42.1, priority:high]
    ```

Isso mostra tanto os metadados completos exibidos pela operação `view()` quanto o subconjunto extraído que imprimimos com `println`.

O método `.subMap()` recebe uma lista de chaves e retorna um novo mapa contendo apenas essas chaves. Se uma chave não existe no mapa original, ela simplesmente não é incluída no resultado.

Isso é particularmente útil quando você precisa criar diferentes versões de metadados para diferentes processos - alguns podem precisar de metadados completos enquanto outros precisam apenas de campos mínimos de identificação.

Agora remova essas instruções println para restaurar seu fluxo de trabalho ao estado anterior, já que não precisamos delas daqui para frente.

!!! tip "Resumo de Operações com Map"

    - **Adicionar chaves**: `map1 + [new_key: value]` - Cria novo mapa com chaves adicionais
    - **Extrair chaves**: `map1.subMap(['key1', 'key2'])` - Cria novo mapa com apenas as chaves especificadas
    - **Ambas as operações criam novos mapas** - Mapas originais permanecem inalterados

#### 1.1.6. Combinando Mapas e Retornando Resultados

Até agora, apenas retornamos o que a comunidade Nextflow chama de 'meta map', e ignoramos os arquivos aos quais esses metadados se relacionam. Mas se você está escrevendo fluxos de trabalho Nextflow, provavelmente quer fazer algo com esses arquivos.

Vamos gerar uma estrutura de canal compreendendo uma tupla de 2 elementos: o mapa de metadados enriquecido e o caminho de arquivo correspondente. Este é um padrão comum em Nextflow para passar dados para processos.

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

Essa estrutura de tupla `[meta, file]` é um padrão comum em Nextflow para passar tanto metadados quanto arquivos associados para processos.

!!! note

    **Mapas e Metadados**: Mapas são fundamentais para trabalhar com metadados no Nextflow. Para uma explicação mais detalhada sobre trabalhar com mapas de metadados, veja a missão secundária [Trabalhando com metadados](./metadata.md).

Nosso fluxo de trabalho demonstra o padrão central: **operações de fluxo de dados** (`workflow`, `channel.fromPath()`, `.splitCsv()`, `.map()`, `.view()`) orquestram como os dados se movem através do pipeline, enquanto **scripting** (mapas `[key: value]`, métodos de string, conversões de tipo, operadores ternários) dentro da closure `.map()` lida com a transformação de itens de dados individuais.

### 1.2. Entendendo Diferentes Tipos: Channel vs List

Até aqui, tudo bem, podemos distinguir entre operações de fluxo de dados e scripting. Mas e quando o mesmo nome de método existe em ambos os contextos?

Um exemplo perfeito é o método `collect`, que existe tanto para tipos de canal quanto para tipos List na biblioteca padrão Nextflow. O método `collect()` em uma List transforma cada elemento, enquanto o operador `collect()` em um canal reúne todas as emissões do canal em um canal de item único.

Vamos demonstrar isso com alguns dados de exemplo, começando por refrescar o que o operador `collect()` de canal faz. Confira `collect.nf`:

```groovy title="collect.nf" linenums="1"
def sample_ids = ['sample_001', 'sample_002', 'sample_003']

// channel.collect() - agrupa múltiplas emissões de canal em uma
ch_input = channel.fromList(sample_ids)
ch_input.view { sample -> "Item individual do canal: ${sample}" }
ch_collected = ch_input.collect()
ch_collected.view { list -> "Resultado de channel.collect(): ${list} (${list.size()} itens agrupados em 1)" }
```

Passos:

- Definir uma List de IDs de amostra
- Criar um canal com `fromList()` que emite cada ID de amostra separadamente
- Imprimir cada item com `view()` conforme ele flui
- Reunir todos os itens em uma única lista com o operador `collect()` do canal
- Imprimir o resultado coletado (item único contendo todos os IDs de amostra) com um segundo `view()`

Mudamos a estrutura do canal, mas não mudamos os dados em si.

Execute o fluxo de trabalho para confirmar isso:

```bash
nextflow run collect.nf
```

??? success "Saída do comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `collect.nf` [loving_mendel] DSL2 - revision: e8d054a46e

    Item individual do canal: sample_001
    Item individual do canal: sample_002
    Item individual do canal: sample_003
    Resultado de channel.collect(): [sample_001, sample_002, sample_003] (3 itens agrupados em 1)
    ```

`view()` retorna uma saída para cada emissão de canal, então sabemos que esta única saída contém todos os 3 itens originais agrupados em uma lista.

Agora vamos ver o método `collect` em uma List em ação. Modifique `collect.nf` para aplicar o método `collect` da List à lista original de IDs de amostra:

=== "Depois"

    ```groovy title="main.nf" linenums="1" hl_lines="9-13"
    def sample_ids = ['sample_001', 'sample_002', 'sample_003']

    // channel.collect() - agrupa múltiplas emissões de canal em uma
    ch_input = channel.fromList(sample_ids)
    ch_input.view { sample -> "Item individual do canal: ${sample}" }
    ch_collected = ch_input.collect()
    ch_collected.view { list -> "Resultado de channel.collect(): ${list} (${list.size()} itens agrupados em 1)" }

    // List.collect() - transforma cada elemento, preserva estrutura
    def formatted_ids = sample_ids.collect { id ->
        id.toUpperCase().replace('SAMPLE_', 'SPECIMEN_')
    }
    println "Resultado de List.collect(): ${formatted_ids} (${sample_ids.size()} itens transformados em ${formatted_ids.size()})"
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="1"
    def sample_ids = ['sample_001', 'sample_002', 'sample_003']

    // channel.collect() - agrupa múltiplas emissões de canal em uma
    ch_input = channel.fromList(sample_ids)
    ch_input.view { sample -> "Item individual do canal: ${sample}" }
    ch_collected = ch_input.collect()
    ch_collected.view { list -> "Resultado de channel.collect(): ${list} (${list.size()} itens agrupados em 1)" }
    ```

Neste novo trecho nós:

- Definimos uma nova variável `formatted_ids` que usa o método `collect` da List para transformar cada ID de amostra na lista original
- Imprimimos o resultado usando `println`

Execute o fluxo de trabalho modificado:

```bash
nextflow run collect.nf
```

??? success "Saída do comando"

    ```console hl_lines="5"
    N E X T F L O W   ~  version 25.10.2

    Launching `collect.nf` [cheeky_stonebraker] DSL2 - revision: 2d5039fb47

    Resultado de List.collect(): [SPECIMEN_001, SPECIMEN_002, SPECIMEN_003] (3 itens transformados em 3)
    Item individual do canal: sample_001
    Item individual do canal: sample_002
    Item individual do canal: sample_003
    Resultado de channel.collect(): [sample_001, sample_002, sample_003] (3 itens agrupados em 1)
    ```

Desta vez, NÃO mudamos a estrutura dos dados, ainda temos 3 itens na lista, mas transformamos cada item usando o método `collect` da List para produzir uma nova lista com valores modificados. Isso é similar a usar o operador `map` em um canal, mas está operando em uma estrutura de dados List em vez de um canal.

`collect` é um caso extremo que estamos usando aqui para enfatizar um ponto. A lição chave é que quando você está escrevendo fluxos de trabalho, sempre distinga entre **estruturas de dados** (Lists, Maps, etc.) e **canais** (construções de fluxo de dados). Operações podem compartilhar nomes mas se comportar completamente diferente dependendo do tipo em que são chamadas.

### 1.3. O Operador Spread (`*.`) - Atalho para Extração de Propriedades

Relacionado ao método `collect` da List está o operador spread (`*.`), que fornece uma maneira concisa de extrair propriedades de coleções. É essencialmente açúcar sintático para um padrão comum de `collect`.

Vamos adicionar uma demonstração ao nosso arquivo `collect.nf`:

=== "Depois"

    ```groovy title="collect.nf" linenums="1" hl_lines="15-18"
    def sample_ids = ['sample_001', 'sample_002', 'sample_003']

    // channel.collect() - agrupa múltiplas emissões de canal em uma
    ch_input = channel.fromList(sample_ids)
    ch_input.view { sample -> "Item individual do canal: ${sample}" }
    ch_collected = ch_input.collect()
    ch_collected.view { list -> "Resultado de channel.collect(): ${list} (${list.size()} itens agrupados em 1)" }

    // List.collect() - transforma cada elemento, preserva estrutura
    def formatted_ids = sample_ids.collect { id ->
        id.toUpperCase().replace('SAMPLE_', 'SPECIMEN_')
    }
    println "Resultado de List.collect(): ${formatted_ids} (${sample_ids.size()} itens transformados em ${formatted_ids.size()})"

    // Operador spread - acesso conciso a propriedades
    def sample_data = [[id: 's1', quality: 38.5], [id: 's2', quality: 42.1], [id: 's3', quality: 35.2]]
    def all_ids = sample_data*.id
    println "Resultado do operador spread: ${all_ids}"
    ```

=== "Antes"

    ```groovy title="collect.nf" linenums="1"
    def sample_ids = ['sample_001', 'sample_002', 'sample_003']

    // channel.collect() - agrupa múltiplas emissões de canal em uma
    ch_input = channel.fromList(sample_ids)
    ch_input.view { sample -> "Item individual do canal: ${sample}" }
    ch_collected = ch_input.collect()
    ch_collected.view { list -> "Resultado de channel.collect(): ${list} (${list.size()} itens agrupados em 1)" }

    // List.collect() - transforma cada elemento, preserva estrutura
    def formatted_ids = sample_ids.collect { id ->
        id.toUpperCase().replace('SAMPLE_', 'SPECIMEN_')
    }
    println "Resultado de List.collect(): ${formatted_ids} (${sample_ids.size()} itens transformados em ${formatted_ids.size()})"
    ```

Execute o fluxo de trabalho atualizado:

```bash title="Testar operador spread"
nextflow run collect.nf
```

??? success "Saída do comando"

    ```console hl_lines="6"
    N E X T F L O W   ~  version 25.10.2

    Launching `collect.nf` [cranky_galileo] DSL2 - revision: 5f3c8b2a91

    Resultado de List.collect(): [SPECIMEN_001, SPECIMEN_002, SPECIMEN_003] (3 itens transformados em 3)
    Resultado do operador spread: [s1, s2, s3]
    Item individual do canal: sample_001
    Item individual do canal: sample_002
    Item individual do canal: sample_003
    Resultado de channel.collect(): [sample_001, sample_002, sample_003] (3 itens agrupados em 1)
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

- **Fluxo de dados vs scripting**: Operadores de canal orquestram como os dados fluem através do seu pipeline, enquanto scripting transforma itens de dados individuais
- **Entendendo tipos**: O mesmo nome de método (como `collect`) pode se comportar diferentemente dependendo do tipo em que é chamado (Channel vs List)
- **Contexto importa**: Sempre esteja ciente se você está trabalhando com canais (fluxo de dados) ou estruturas de dados (scripting)

Entender esses limites é essencial para debugging, documentação e escrever fluxos de trabalho de fácil manutenção.

A seguir vamos nos aprofundar em capacidades de processamento de strings, que são essenciais para lidar com dados do mundo real.

---

## 2. Processamento de Strings e Geração Dinâmica de Scripts

Dominar o processamento de strings separa fluxos de trabalho frágeis de pipelines robustos. Esta seção cobre análise de nomes de arquivo complexos, geração dinâmica de scripts e interpolação de variáveis.

### 2.1. Correspondência de Padrões e Expressões Regulares

Arquivos de bioinformática frequentemente têm convenções de nomenclatura complexas que codificam metadados. Vamos extrair isso automaticamente usando correspondência de padrões com expressões regulares.

Vamos retornar ao nosso fluxo de trabalho `main.nf` e adicionar alguma lógica de correspondência de padrões para extrair informações adicionais da amostra de nomes de arquivo. Os arquivos FASTQ em nosso conjunto de dados seguem convenções de nomenclatura estilo Illumina com nomes como `SAMPLE_001_S1_L001_R1_001.fastq.gz`. Estes podem parecer crípticos, mas na verdade codificam metadados úteis como ID da amostra, número de lane e direção de leitura. Vamos usar capacidades regex para analisar esses nomes.

Faça a seguinte mudança em seu fluxo de trabalho `main.nf` existente:

=== "Depois"

    ```groovy title="main.nf" linenums="4" hl_lines="10-21"
            .map { row ->
                // Scripting para transformação de dados
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
                // Scripting para transformação de dados
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
2. **Correspondência de padrões** com o operador `=~` - isso tenta corresponder uma string a um padrão regex
3. **Objetos Matcher** que capturam grupos com `[0][1]`, `[0][2]`, etc. - `[0]` refere-se à correspondência inteira, `[1]`, `[2]`, etc. referem-se a grupos capturados entre parênteses

Vamos detalhar o padrão regex `^(.+)_S(\d+)_L(\d{3})_(R[12])_(\d{3})\.fastq(?:\.gz)?$`:

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

```bash title="Testar correspondência de padrões"
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

Isso mostra os metadados enriquecidos a partir dos nomes de arquivo.

### 2.2. Geração Dinâmica de Scripts em Processos

Blocos script de processo são essencialmente strings multi-linha que são passadas para o shell. Você pode usar **lógica condicional** (if/else, operadores ternários) para gerar dinamicamente diferentes strings de script com base nas características da entrada. Isso é essencial para lidar com diversos tipos de entrada—como leituras single-end vs paired-end—sem duplicar definições de processo.

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

O processo recebe arquivos FASTQ como entrada e executa a ferramenta `fastp` para aparar adaptadores e filtrar leituras de baixa qualidade. Infelizmente, a pessoa que escreveu este processo não permitiu as leituras single-end que temos em nosso conjunto de dados de exemplo. Vamos adicioná-lo ao nosso fluxo de trabalho e ver o que acontece:

Primeiro, inclua o módulo na primeira linha do seu fluxo de trabalho `main.nf`:

```groovy title="main.nf" linenums="1"
include { FASTP } from './modules/fastp.nf'
```

Em seguida, modifique o bloco `workflow` para conectar o canal `ch_samples` ao processo `FASTP`:

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

Você pode ver que o processo está tentando executar `fastp` com um valor `null` para o segundo arquivo de entrada, o que está causando a falha. Isso ocorre porque nosso conjunto de dados contém leituras single-end, mas o processo está codificado para esperar leituras paired-end (dois arquivos de entrada por vez).

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

Agora o fluxo de trabalho pode lidar graciosamente com leituras tanto single-end quanto paired-end. A lógica condicional verifica o número de arquivos de entrada e constrói o comando apropriado para `fastp`. Vamos ver se funciona:

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

```console title="Verificar comandos executados"
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

Outro uso comum de lógica dinâmica de script pode ser visto em [o módulo Genomics do Nextflow for Science](../../nf4science/genomics/02_joint_calling). Nesse módulo, o processo GATK sendo chamado pode receber múltiplos arquivos de entrada, mas cada um deve ser prefixado com `-V` para formar uma linha de comando correta. O processo usa scripting para transformar uma coleção de arquivos de entrada (`all_gvcfs`) nos argumentos de comando corretos:

```groovy title="manipulação de linha de comando para GATK" linenums="1" hl_lines="2 5"
    script:
    def gvcfs_line = all_gvcfs.collect { gvcf -> "-V ${gvcf}" }.join(' ')
    """
    gatk GenomicsDBImport \
        ${gvcfs_line} \
        -L ${interval_list} \
        --genomicsdb-workspace-path ${cohort_name}_gdb
    """
```

Esses padrões de usar scripting em blocos script de processo são extremamente poderosos e podem ser aplicados em muitos cenários - desde lidar com tipos de entrada variáveis até construir argumentos complexos de linha de comando a partir de coleções de arquivos, tornando seus processos verdadeiramente adaptáveis aos requisitos diversos de dados do mundo real.

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
- **Geração dinâmica de scripts**: Usando lógica condicional (if/else, operadores ternários) para gerar diferentes strings de script com base nas características da entrada
- **Interpolação de variáveis**: Entendendo quando o Nextflow interpreta strings vs quando o shell interpreta
  - `${var}`
