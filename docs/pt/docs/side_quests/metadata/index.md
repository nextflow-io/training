# Metadados e meta maps

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tradução assistida por IA - [saiba mais e sugira melhorias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Em qualquer análise científica, raramente trabalhamos apenas com os arquivos de dados brutos.
Cada arquivo vem com suas próprias informações adicionais: o que é, de onde veio e o que o torna especial.
Essa informação extra é o que chamamos de metadados.

Metadados são dados que descrevem outros dados.
Os metadados rastreiam detalhes importantes sobre arquivos e condições experimentais, e ajudam a adaptar as análises às características únicas de cada conjunto de dados.

Pense nisso como um catálogo de biblioteca: enquanto os livros contêm o conteúdo real (dados brutos), as fichas catalográficas fornecem informações essenciais sobre cada livro — quando foi publicado, quem o escreveu, onde encontrá-lo (metadados).
Em pipelines Nextflow, os metadados podem ser usados para:

- Rastrear informações específicas de cada arquivo ao longo do fluxo de trabalho
- Configurar processos com base nas características dos arquivos
- Agrupar arquivos relacionados para análise conjunta

### Objetivos de aprendizado

Nesta missão paralela, vamos explorar como lidar com metadados em fluxos de trabalho.
Começando com uma planilha simples (frequentemente chamada de samplesheet em bioinformática) contendo informações básicas sobre os arquivos, você aprenderá como:

- Ler e analisar metadados de arquivos CSV
- Criar e manipular meta maps
- Adicionar novos campos de metadados durante a execução do fluxo de trabalho
- Usar metadados para personalizar o comportamento dos processos

Essas habilidades vão ajudá-lo a construir pipelines mais robustos e flexíveis, capazes de lidar com relações complexas entre arquivos e requisitos de processamento.

### Pré-requisitos

Antes de embarcar nesta missão paralela, você deve:

- Ter concluído o tutorial [Hello Nextflow](../hello_nextflow/README.md) ou um curso equivalente para iniciantes.
- Estar confortável com os conceitos e mecanismos básicos do Nextflow (processos, canais, operadores)

---

## 0. Primeiros passos

#### Abra o codespace de treinamento

Se ainda não tiver feito isso, certifique-se de abrir o ambiente de treinamento conforme descrito em [Configuração do Ambiente](../envsetup/index.md).

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### Acesse o diretório do projeto

Vamos acessar o diretório onde estão os arquivos deste tutorial.

```bash
cd side-quests/metadata
```

Você pode configurar o VSCode para focar neste diretório:

```bash
code .
```

#### Revise os materiais

Você encontrará um arquivo de fluxo de trabalho principal e um diretório `data` contendo uma planilha e alguns arquivos de dados.

??? abstract "Conteúdo do diretório"

    ```console
    .
    ├── data
    │   ├── bonjour.txt
    │   ├── ciao.txt
    │   ├── guten_tag.txt
    │   ├── hallo.txt
    │   ├── hello.txt
    │   ├── hola.txt
    │   ├── salut.txt
    │   └── datasheet.csv
    ├── main.nf
    └── nextflow.config
    ```

O fluxo de trabalho no arquivo `main.nf` é um esboço que você vai expandir gradualmente até se tornar um fluxo de trabalho totalmente funcional.

A planilha lista os caminhos para os arquivos de dados e alguns metadados associados, organizados em 3 colunas:

- `id`: autoexplicativo, um ID atribuído ao arquivo
- `character`: um nome de personagem, que usaremos mais tarde para desenhar diferentes criaturas
- `data`: caminhos para arquivos `.txt` que contêm saudações em diferentes idiomas

```console title="datasheet.csv"
id,character,recording
sampleA,squirrel,/workspaces/training/side-quests/metadata/data/bonjour.txt
sampleB,tux,/workspaces/training/side-quests/metadata/data/guten_tag.txt
sampleC,sheep,/workspaces/training/side-quests/metadata/data/hallo.txt
sampleD,turkey,/workspaces/training/side-quests/metadata/data/hello.txt
sampleE,stegosaurus,/workspaces/training/side-quests/metadata/data/hola.txt
sampleF,moose,/workspaces/training/side-quests/metadata/data/salut.txt
sampleG,turtle,/workspaces/training/side-quests/metadata/data/ciao.txt
```

Cada arquivo de dados contém algum texto de saudação em um dos cinco idiomas (fr: francês, de: alemão, es: espanhol, it: italiano, en: inglês).

Também forneceremos uma ferramenta de análise de idiomas em contêiner chamada `langid`.

#### Revise a tarefa

Seu desafio é escrever um fluxo de trabalho Nextflow que irá:

1. **Identificar** o idioma em cada arquivo automaticamente
2. **Agrupar** os arquivos por família linguística (línguas germânicas vs. românicas)
3. **Personalizar** o processamento de cada arquivo com base em seu idioma e metadados
4. **Organizar** as saídas por grupo linguístico

Isso representa um padrão típico de fluxo de trabalho onde metadados específicos de cada arquivo orientam as decisões de processamento — exatamente o tipo de problema que os meta maps resolvem de forma elegante.

#### Lista de verificação de prontidão

Acha que está pronto para mergulhar de cabeça?

- [ ] Entendo o objetivo deste curso e seus pré-requisitos
- [ ] Meu codespace está funcionando
- [ ] Defini meu diretório de trabalho corretamente
- [ ] Entendo a tarefa

Se você conseguir marcar todas as caixas, pode começar.

---

## 1. Carregar metadados de uma planilha

Abra o arquivo de fluxo de trabalho `main.nf` para examinar o esboço que estamos fornecendo como ponto de partida.

```groovy title="main.nf" linenums="1"
#!/usr/bin/env nextflow

workflow  {

    ch_datasheet = channel.fromPath("./data/datasheet.csv")

}
```

Você pode ver que configuramos um canal básico para carregar a planilha de exemplo como um arquivo, mas isso ainda não vai ler o conteúdo do arquivo.
Vamos começar adicionando isso.

### 1.1. Ler o conteúdo com `splitCsv`

Precisamos escolher um operador que analise o conteúdo do arquivo de forma adequada com o mínimo de esforço da nossa parte.
Como nossa planilha está no formato CSV, este é um trabalho para o operador [`splitCsv`](https://www.nextflow.io/docs/latest/reference/operator.html#splitcsv), que carrega cada linha do arquivo como um elemento no canal.

Faça as seguintes alterações para adicionar uma operação `splitCsv()` ao código de construção do canal, além de uma operação `view()` para verificar se o conteúdo do arquivo está sendo carregado corretamente no canal.

=== "Depois"

    ```groovy title="main.nf" linenums="3" hl_lines="4-5"
    workflow  {

        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .view()

    }
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="3"
    workflow  {

        ch_datasheet = channel.fromPath("./data/datasheet.csv")

    }
    ```

Observe que estamos usando a opção `header: true` para informar ao Nextflow que deve ler a primeira linha do arquivo CSV como a linha de cabeçalho.

Vamos ver o que sai disso, não é?
Execute o fluxo de trabalho:

```bash
nextflow run main.nf
```

??? success "Saída do comando"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [exotic_albattani] DSL2 - revision: c0d03cec83

    [id:sampleA, character:squirrel, recording:/workspaces/training/side-quests/metadata/data/bonjour.txt]
    [id:sampleB, character:tux, recording:/workspaces/training/side-quests/metadata/data/guten_tag.txt]
    [id:sampleC, character:sheep, recording:/workspaces/training/side-quests/metadata/data/hallo.txt]
    [id:sampleD, character:turkey, recording:/workspaces/training/side-quests/metadata/data/hello.txt]
    [id:sampleE, character:stegosaurus, recording:/workspaces/training/side-quests/metadata/data/hola.txt]
    [id:sampleF, character:moose, recording:/workspaces/training/side-quests/metadata/data/salut.txt]
    [id:sampleG, character:turtle, recording:/workspaces/training/side-quests/metadata/data/ciao.txt]
    ```

Podemos ver que o operador construiu um map de pares chave-valor para cada linha do arquivo CSV, com os cabeçalhos das colunas como chaves para os valores correspondentes.

Cada entrada do map corresponde a uma coluna em nossa planilha:

- `id`
- `character`
- `recording`

Ótimo! Isso facilita o acesso a campos específicos de cada arquivo.
Por exemplo, poderíamos acessar o ID do arquivo com `id` ou o caminho do arquivo txt com `recording`.

??? info "(Opcional) Mais sobre maps"

    Em Groovy, a linguagem de programação sobre a qual o Nextflow é construído, um map é uma estrutura de dados de chave-valor semelhante a dicionários em Python, objetos em JavaScript ou hashes em Ruby.

    Aqui está um script executável que mostra como você pode definir um map e acessar seu conteúdo na prática:

    ```groovy title="examples/map_demo.nf"
    #!/usr/bin/env nextflow

    // Cria um map simples
    def my_map = [id:'sampleA', character:'squirrel']

    // Imprime o map completo
    println "map: ${my_map}"

    // Acessa valores individuais usando notação de ponto
    println "id: ${my_map.id}"
    println "character: ${my_map.character}"
    ```

    Mesmo sem um bloco `workflow` adequado, o Nextflow pode executar isso como se fosse um fluxo de trabalho:

    ```bash
    nextflow run examples/map_demo.nf
    ```

    E aqui está o que você pode esperar ver na saída:

    ```console title="Output"
     N E X T F L O W   ~  version 25.10.2

    Launching `map_demo.nf` [cheesy_plateau] DSL2 - revision: fae5b8496e

    map: [id:sampleA, character:squirrel]
    id: sampleA
    character: squirrel
    ```

### 1.2. Selecionar campos específicos com `map`

Digamos que queremos acessar a coluna `character` da planilha e imprimi-la.
Podemos usar o operador `map` do Nextflow para iterar sobre cada item em nosso canal e selecionar especificamente a entrada `character` do objeto map.

Faça as seguintes edições no fluxo de trabalho:

=== "Depois"

    ```groovy title="main.nf" linenums="3" hl_lines="5-7"
    workflow  {

        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map{ row ->
                row.character
            }
            .view()

    }
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="3"
    workflow  {

        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .view()

    }
    ```

Agora execute o fluxo de trabalho novamente:

```bash
nextflow run main.nf
```

??? success "Saída do comando"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [exotic_albattani] DSL2 - revision: c0d03cec83

    squirrel
    tux
    sheep
    turkey
    stegosaurus
    moose
    turtle
    ```

Sucesso! Aproveitamos a estrutura de map derivada de nossa planilha para acessar os valores de colunas individuais para cada linha.

Agora que lemos com sucesso a planilha e temos acesso aos dados em cada linha, podemos começar a implementar a lógica do nosso pipeline.

### 1.3. Organizar os metadados em um 'meta map'

No estado atual do fluxo de trabalho, os arquivos de entrada (sob a chave `recording`) e os metadados associados (`id`, `character`) estão todos no mesmo nível, como se estivessem todos em um grande saco.
A consequência prática é que todo processo que consome este canal precisaria ser configurado com essa estrutura em mente:

```groovy
    input:
    tuple val(id), val(character), file(recording)
```

Isso funciona bem enquanto o número de colunas na planilha não mudar.
No entanto, se você adicionar apenas uma coluna à planilha, a forma do canal não corresponderá mais ao que o processo espera, e o fluxo de trabalho produzirá erros.
Isso também torna o processo difícil de compartilhar com outras pessoas que possam ter dados de entrada ligeiramente diferentes, e você pode acabar tendo que codificar variáveis diretamente no processo que não são necessárias pelo bloco de script.

Para evitar esse problema, precisamos encontrar uma maneira de manter a estrutura do canal consistente independentemente de quantas colunas a planilha contenha.

Podemos fazer isso coletando todos os metadados em um item dentro da tupla, que chamaremos de meta map de metadados, ou mais simplesmente 'meta map'.

Faça as seguintes edições na operação `map`:

=== "Depois"

    ```groovy title="main.nf" linenums="5" hl_lines="4"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                    [ [id: row.id, character: row.character], row.recording ]
            }
            .view()
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="5" hl_lines="4"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map{ row ->
                row.character
            }
            .view()
    ```

Reestruturamos os elementos do canal em uma tupla composta por dois elementos: o meta map e o objeto de arquivo correspondente.

Vamos executar o fluxo de trabalho:

```bash
nextflow run main.nf
```

??? success "Saída do comando"

    ```console title="View meta map"
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [lethal_booth] DSL2 - revision: 0d8f844c07

    [[id:sampleA, character:squirrel], /workspaces/training/side-quests/metadata/data/bonjour.txt]
    [[id:sampleB, character:tux], /workspaces/training/side-quests/metadata/data/guten_tag.txt]
    [[id:sampleC, character:sheep], /workspaces/training/side-quests/metadata/data/hallo.txt]
    [[id:sampleD, character:turkey], /workspaces/training/side-quests/metadata/data/hello.txt]
    [[id:sampleE, character:stegosaurus], /workspaces/training/side-quests/metadata/data/hola.txt]
    [[id:sampleF, character:moose], /workspaces/training/side-quests/metadata/data/salut.txt]
    [[id:sampleG, character:turtle], /workspaces/training/side-quests/metadata/data/ciao.txt]
    ```

Agora, cada elemento no canal contém o meta map primeiro e o objeto de arquivo correspondente em segundo:

```console title="Example output structure"
[
  [id:sampleA, character:squirrel],
  /workspaces/training/side-quests/metadata/data/bonjour.txt
]
```

Como resultado, adicionar mais colunas na planilha disponibilizará mais metadados no map `meta`, mas não mudará a forma do canal.
Isso nos permite escrever processos que consomem o canal sem precisar codificar os itens de metadados diretamente na especificação de entrada:

```groovy title="Syntax example"
    input:
    tuple val(meta), file(recording)
```

Esta é uma convenção amplamente utilizada para organizar metadados em fluxos de trabalho Nextflow.

### Conclusão

Nesta seção, você aprendeu:

- **Por que os metadados são importantes:** Manter os metadados junto com seus dados preserva informações importantes sobre os arquivos ao longo do fluxo de trabalho.
- **Como ler planilhas:** Usando `splitCsv` para ler arquivos CSV com informações de cabeçalho e transformar linhas em dados estruturados
- **Como criar um meta map:** Separando metadados dos dados de arquivo usando a estrutura de tupla `[ [id:valor, ...], arquivo ]`

---

## 2. Manipulando metadados

Agora que temos nossos metadados carregados, vamos fazer algo com eles!

Vamos usar uma ferramenta chamada [`langid`](https://github.com/saffsd/langid.py) para identificar o idioma contido no arquivo de gravação de cada criatura.
A ferramenta vem pré-treinada em um conjunto de idiomas e, dado um trecho de texto, ela produzirá uma previsão de idioma e uma pontuação de probabilidade associada, ambas para `stdout`.

### 2.1. Importar o processo e examinar o código

Fornecemos um módulo de processo pré-escrito chamado `IDENTIFY_LANGUAGE` que encapsula a ferramenta `langid`, então você só precisa adicionar uma instrução include antes do bloco workflow.

Faça a seguinte edição no fluxo de trabalho:

=== "Depois"

    ```groovy title="main.nf" linenums="1" hl_lines="3"
    #!/usr/bin/env nextflow

    include { IDENTIFY_LANGUAGE } from './modules/langid.nf'

    workflow {
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="1"
    #!/usr/bin/env nextflow

    workflow {
    ```

Você pode abrir o arquivo do módulo para examinar seu código:

```groovy title="modules/langid.nf" linenums="1" hl_lines="9 12"
#!/usr/bin/env nextflow

// Usa langid para prever o idioma de cada arquivo de entrada
process IDENTIFY_LANGUAGE {

    container 'community.wave.seqera.io/library/pip_langid:b2269f456a5629ff'

    input:
    tuple val(meta), path(file)

    output:
    tuple val(meta), path(file), stdout

    script:
    """
    langid < ${file} -l en,de,fr,es,it | sed -E "s/.*\\('([a-z]+)'.*/\\1/" | tr -d '\\n'
    """
}
```

Como você pode ver, a definição de entrada usa a mesma estrutura `tuple val(meta), path(file)` que acabamos de aplicar ao nosso canal de entrada.

A definição de saída é estruturada como uma tupla com estrutura semelhante à da entrada, exceto que também contém `stdout` como terceiro elemento.
Esse padrão `tuple val(meta), path(file), <saída>` mantém os metadados associados tanto aos dados de entrada quanto às saídas à medida que fluem pelo pipeline.

Observe que estamos usando o qualificador de saída [`stdout`](https://www.nextflow.io/docs/latest/process.html#outputs) do Nextflow aqui porque a ferramenta imprime sua saída diretamente no console em vez de escrever um arquivo; e usamos `sed` na linha de comando para remover a pontuação de probabilidade, limpar a string removendo caracteres de nova linha e retornar apenas a previsão de idioma.

### 2.2. Adicionar uma chamada ao `IDENTIFY_LANGUAGE`

Agora que o processo está disponível para o fluxo de trabalho, podemos adicionar uma chamada ao processo `IDENTIFY_LANGUAGE` para executá-lo no canal de dados.

Faça as seguintes edições no fluxo de trabalho:

=== "Depois"

    ```groovy title="main.nf" linenums="7" hl_lines="7-9"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                [[id: row.id, character: row.character], row.recording]
            }

        // Executa langid para identificar o idioma de cada saudação
        IDENTIFY_LANGUAGE(ch_datasheet)
        IDENTIFY_LANGUAGE.out.view()
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="7" hl_lines="6"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                [[id: row.id, character: row.character], row.recording]
            }
            .view()
    ```

Observe que removemos a operação `.view()` original na construção do canal.

Agora podemos executar o fluxo de trabalho:

```bash
nextflow run main.nf
```

??? success "Saída do comando"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [voluminous_mcnulty] DSL2 - revision: f9bcfebabb

    executor >  local (7)
    [4e/f722fe] IDENTIFY_LANGUAGE (7) [100%] 7 of 7 ✔
    [[id:sampleA, character:squirrel], /workspaces/training/side-quests/metadata/work/eb/f7148ebdd898fbe1136bec6a714acb/bonjour.txt, fr]
    [[id:sampleB, character:tux], /workspaces/training/side-quests/metadata/work/16/71d72410952c22cd0086d9bca03680/guten_tag.txt, de]
    [[id:sampleD, character:turkey], /workspaces/training/side-quests/metadata/work/c4/b7562adddc1cc0b7d414ec45d436eb/hello.txt, en]
    [[id:sampleC, character:sheep], /workspaces/training/side-quests/metadata/work/ea/04f5d979429e4455e14b9242fb3b45/hallo.txt, de]
    [[id:sampleF, character:moose], /workspaces/training/side-quests/metadata/work/5a/6c2b84bf8fadb98e28e216426be079/salut.txt, fr]
    [[id:sampleE, character:stegosaurus], /workspaces/training/side-quests/metadata/work/af/ee7c69bcab891c40d0529305f6b9e7/hola.txt, es]
    [[id:sampleG, character:turtle], /workspaces/training/side-quests/metadata/work/4e/f722fe47271ba7ebcd69afa42964ca/ciao.txt, it]
    ```

Excelente! Agora temos uma previsão de qual idioma cada personagem fala.

E como observado anteriormente, também incluímos o arquivo de entrada e o meta map na saída, o que significa que ambos permanecem associados às novas informações que acabamos de produzir.
Isso será útil no próximo passo.

!!! note "Nota"

    De forma mais geral, esse padrão de manter o meta map associado aos resultados facilita a associação de resultados relacionados que compartilham os mesmos identificadores.

    Como você já deve ter aprendido, não é possível confiar na ordem dos itens nos canais para associar resultados entre eles.
    Em vez disso, você deve usar chaves para associar os dados corretamente, e os meta maps fornecem uma estrutura ideal para esse propósito.

    Exploramos esse caso de uso em detalhes na missão paralela [Splitting & Grouping](../splitting_and_grouping/).

### 2.3. Enriquecer os metadados com saídas de processos

Dado que os resultados que acabamos de produzir são em si uma forma de metadados sobre o conteúdo dos arquivos, seria útil adicioná-los ao nosso meta map.

No entanto, não queremos modificar o meta map existente diretamente.
Do ponto de vista técnico, é _possível_ fazer isso, mas não é seguro.

Então, em vez disso, criaremos um novo meta map contendo o conteúdo do meta map existente mais um novo par chave-valor `lang: lang_id` contendo as novas informações, usando o operador `+` (um recurso do Groovy).
E combinaremos isso com uma operação [`map`](https://www.nextflow.io/docs/latest/operator.html#map) para substituir o map antigo pelo novo.

Aqui estão as edições que você precisa fazer no fluxo de trabalho:

=== "Depois"

    ```groovy title="main.nf" linenums="13" hl_lines="3-7"
        // Executa langid para identificar o idioma de cada saudação
        IDENTIFY_LANGUAGE(ch_datasheet)
        IDENTIFY_LANGUAGE.out
            .map { meta, file, lang_id ->
                [meta + [lang: lang_id], file]
            }
            .view()
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="13" hl_lines="3"
        // Executa langid para identificar o idioma de cada saudação
        IDENTIFY_LANGUAGE(ch_datasheet)
        IDENTIFY_LANGUAGE.out.view()
    ```

Se você ainda não está familiarizado com o operador `+`, ou se isso parece confuso, reserve alguns minutos para ler a explicação detalhada abaixo.

??? info "Criação do novo meta map usando o operador `+`"

    **Primeiro, você precisa saber que podemos mesclar o conteúdo de dois maps usando o operador Groovy `+`.**

    Digamos que temos os seguintes maps:

    ```groovy
    map1 = [id: 'sampleA', character: 'squirrel']
    map2 = [lang: 'fr']
    ```

    Podemos mesclá-los assim:

    ```groovy
    new_map = map1 + map2
    ```

    O conteúdo de `new_map` será:

    ```groovy
    [id: 'sampleA', character: 'squirrel', lang: 'fr']
    ```

    Ótimo!

    **Mas e se você precisar adicionar um campo que ainda não faz parte de um map?**

    Digamos que você começa novamente a partir de `map1`, mas a previsão de idioma não está em seu próprio map (não existe `map2`).
    Em vez disso, ela está armazenada em uma variável chamada `lang_id`, e você sabe que quer armazenar seu valor (`'fr'`) com a chave `lang`.

    Você pode fazer o seguinte:

    ```groovy
    new_map = [map1 + [lang: lang_id]]
    ```

    Aqui, `[lang: new_info]` cria um novo map sem nome na hora, e `map1 + ` mescla `map1` com o novo map sem nome, produzindo o mesmo conteúdo de `new_map` que antes.

    Legal, não é?

    **Agora vamos transpor isso para o contexto de uma operação `channel.map()` do Nextflow.**

    O código se torna:

    ```groovy
    .map { map1, lang_id ->
        [map1 + [lang: lang_id]]
    }
    ```

    Isso faz o seguinte:

    - `map1, lang_id ->` recebe os dois itens na tupla
    - `[map1 + [lang: lang_id]]` cria o novo map conforme detalhado acima

    A saída é um único map sem nome com o mesmo conteúdo que `new_map` em nosso exemplo acima.
    Então efetivamente transformamos:

    ```groovy
    [id: 'sampleA', character: 'squirrel'], 'fr'
    ```

    em:

    ```groovy
    [id: 'sampleA', character: 'squirrel', lang: 'fr']
    ```

    Esperamos que você possa ver que, se mudarmos `map1` para `meta`, isso é basicamente tudo que precisamos para adicionar a previsão de idioma ao nosso meta map no fluxo de trabalho.

    Exceto por uma coisa!

    No caso do nosso fluxo de trabalho, **também precisamos levar em conta a presença do objeto `file` na tupla**, que é composta por `meta, file, lang_id`.

    Então o código aqui se tornaria:

    ```groovy
    .map { meta, file, lang_id ->
        [meta + [lang: lang_id], file]
    }
    ```

    Se você está tendo dificuldade em entender por que o `file` parece estar se movendo na operação `map`, imagine que em vez de `[meta + [lang: lang_id], file]`, essa linha lê `[new_map, file]`.
    Isso deve deixar mais claro que estamos simplesmente deixando o `file` em seu lugar original na segunda posição na tupla. Apenas pegamos o valor `new_info` e o incorporamos ao map que está na primeira posição.

    **E isso nos traz de volta à estrutura de canal `tuple val(meta), path(file)`!**

Quando você estiver confiante de que entendeu o que esse código está fazendo, execute o fluxo de trabalho para ver se funcionou:

```bash
nextflow run main.nf -resume
```

??? success "Saída do comando"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [cheeky_fermat] DSL2 - revision: d096281ee4

    [4e/f722fe] IDENTIFY_LANGUAGE (7) [100%] 7 of 7, cached: 7 ✔
    [[id:sampleA, character:squirrel, lang:fr], /workspaces/training/side-quests/metadata/work/eb/f7148ebdd898fbe1136bec6a714acb/bonjour.txt]
    [[id:sampleB, character:tux, lang:de], /workspaces/training/side-quests/metadata/work/16/71d72410952c22cd0086d9bca03680/guten_tag.txt]
    [[id:sampleC, character:sheep, lang:de], /workspaces/training/side-quests/metadata/work/ea/04f5d979429e4455e14b9242fb3b45/hallo.txt]
    [[id:sampleD, character:turkey, lang:en], /workspaces/training/side-quests/metadata/work/c4/b7562adddc1cc0b7d414ec45d436eb/hello.txt]
    [[id:sampleF, character:moose, lang:fr], /workspaces/training/side-quests/metadata/work/5a/6c2b84bf8fadb98e28e216426be079/salut.txt]
    [[id:sampleE, character:stegosaurus, lang:es], /workspaces/training/side-quests/metadata/work/af/ee7c69bcab891c40d0529305f6b9e7/hola.txt]
    [[id:sampleG, character:turtle, lang:it], /workspaces/training/side-quests/metadata/work/4e/f722fe47271ba7ebcd69afa42964ca/ciao.txt]
    ```

Sim, está correto!
Reorganizamos cuidadosamente a saída do processo de `meta, file, lang_id` para que `lang_id` seja agora uma das chaves no meta map, e as tuplas do canal se encaixam novamente no modelo `meta, file`.

<!-- TODO (future) Should we also show how to remove a key using subMap?! Or note where to find that. -->

### 2.4. Atribuir um grupo linguístico usando condicionais

Agora que temos nossas previsões de idioma, vamos usar as informações para atribuir novos agrupamentos.

Em nossos dados de exemplo, os idiomas usados pelos nossos personagens podem ser agrupados em línguas germânicas (inglês, alemão) e línguas românicas (francês, espanhol, italiano).
Pode ser útil ter essa classificação prontamente disponível em algum momento posterior no pipeline, então vamos adicionar essa informação no meta map.

E, boas notícias, este é mais um caso que se presta perfeitamente ao uso do operador `map`!

Especificamente, vamos definir uma variável chamada `lang_group`, usar alguma lógica condicional simples para determinar qual valor atribuir ao `lang_group` para cada dado.

A sintaxe geral será assim:

```groovy
.map { meta, file ->

    // lógica condicional definindo lang_group vai aqui

    [meta + [lang_group: lang_group], file]
}
```

Você pode ver que isso é muito semelhante à operação de mesclagem de maps que usamos na etapa anterior.
Só precisamos escrever as instruções condicionais.

Aqui está a lógica condicional que queremos aplicar:

- Definir uma variável chamada `lang_group` com valor padrão `'unknown'`.
- Se `lang` for alemão (`'de'`) ou inglês (`'en'`), alterar `lang_group` para `germanic`.
- Caso contrário, se `lang` estiver incluído em uma lista contendo francês (`'fr'`), espanhol (`'es'`) e italiano (`'it'`), alterar `lang_group` para `romance`.

Tente escrever você mesmo se já souber como escrever instruções condicionais em Nextflow.

!!! tip "Dica"

    Você pode acessar o valor de `lang` dentro da operação map com `meta.lang`.

Você deve acabar fazendo as seguintes alterações no fluxo de trabalho:

=== "Depois"

    ```groovy title="main.nf" linenums="13" hl_lines="7-19"
        // Executa langid para identificar o idioma de cada saudação
        IDENTIFY_LANGUAGE(ch_datasheet)
        IDENTIFY_LANGUAGE.out
            .map { meta, file, lang_id ->
                [meta + [lang: lang_id], file]
            }
            .map { meta, file ->

                def lang_group = "unknown"
                if (meta.lang.equals("de") || meta.lang.equals('en')) {
                    lang_group = "germanic"
                }
                else if (meta.lang in ["fr", "es", "it"]) {
                    lang_group = "romance"
                }

                [meta + [lang_group: lang_group], file]
            }
            .view()
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="13" hl_lines="7"
        // Executa langid para identificar o idioma de cada saudação
        IDENTIFY_LANGUAGE(ch_datasheet)
        IDENTIFY_LANGUAGE.out
            .map { meta, file, lang_id ->
                [meta + [lang: lang_id], file]
            }
            .view()
    ```

Aqui estão os pontos principais:

- Usamos `def lang_group = "unknown"` para criar a variável `lang_group` com valor padrão definido como `unknown`.
- Usamos uma estrutura `if {} else if {}` para a lógica condicional, com testes alternativos `.equals()` para os dois idiomas germânicos, e um teste de existência em uma lista para os três idiomas românicos.
- Usamos a operação de mesclagem `meta + [lang_group:lang_group]` como anteriormente para gerar o meta map atualizado.

<!-- TODO (future) Add note/links to relevant docs in additional resources section -->

Quando tudo isso fizer sentido, execute o fluxo de trabalho novamente para ver o resultado:

```bash
nextflow run main.nf -resume
```

??? success "Saída do comando"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [wise_almeida] DSL2 - revision: 46778c3cd0

    [da/652cc6] IDENTIFY_LANGUAGE (7) [100%] 7 of 7, cached: 7 ✔
    [[id:sampleA, character:squirrel, lang:fr, lang_group:romance], /workspaces/training/side-quests/metadata/data/bonjour.txt]
    [[id:sampleB, character:tux, lang:de, lang_group:germanic], /workspaces/training/side-quests/metadata/data/guten_tag.txt]
    [[id:sampleC, character:sheep, lang:de, lang_group:germanic], /workspaces/training/side-quests/metadata/data/hallo.txt]
    [[id:sampleD, character:turkey, lang:en, lang_group:germanic], /workspaces/training/side-quests/metadata/data/hello.txt]
    [[id:sampleE, character:stegosaurus, lang:es, lang_group:romance], /workspaces/training/side-quests/metadata/data/hola.txt]
    [[id:sampleF, character:moose, lang:fr, lang_group:romance], /workspaces/training/side-quests/metadata/data/salut.txt]
    [[id:sampleG, character:turtle, lang:it, lang_group:romance], /workspaces/training/side-quests/metadata/data/ciao.txt]
    ```

Como você pode ver, os elementos do canal mantêm sua estrutura `[meta, file]`, mas o meta map agora inclui essa nova classificação.

### Conclusão

Nesta seção, você aprendeu como:

- **Aplicar metadados de entrada a canais de saída**: Copiar metadados dessa forma nos permite associar resultados posteriormente com base no conteúdo dos metadados.
- **Criar chaves personalizadas**: Você criou duas novas chaves no seu meta map, mesclando-as com `meta + [nova_chave:valor]` no meta map existente. Uma baseada em um valor calculado por um processo, e outra baseada em uma condição definida no operador `map`.

Isso permite que você associe metadados novos e existentes a arquivos à medida que avança pelo pipeline.
Mesmo que você não esteja usando metadados como parte de um processo, manter o meta map associado aos dados dessa forma facilita manter todas as informações relevantes juntas.

---

## 3. Usando informações do meta map em um processo

Agora que você sabe como criar e atualizar o meta map, podemos chegar à parte realmente divertida: usar os metadados de fato em um processo.

Mais especificamente, vamos adicionar uma segunda etapa ao nosso fluxo de trabalho para desenhar cada animal como arte ASCII e fazê-lo dizer o texto gravado em um balão de fala.
Vamos fazer isso usando uma ferramenta chamada [`cowpy`](https://github.com/jeffbuttars/cowpy).

??? info "O que o `cowpy` faz?"

    `cowpy` é uma ferramenta de linha de comando que gera arte ASCII para exibir entradas de texto arbitrárias de forma divertida.
    É uma implementação em Python da clássica ferramenta [cowsay](https://en.wikipedia.org/wiki/Cowsay) de Tony Monroe.

    ```console
    cowpy "Hello Nextflow"
    ```

    ```console
    ______________________________________________________
    < Hello Nextflow >
    ------------------------------------------------------
        \   ^__^
          \  (oo)\_______
            (__)\       )\/\
              ||----w |
              ||     ||
    ```

    Opcionalmente, você pode selecionar um personagem (ou 'cowacter') para usar em vez da vaca padrão.

    ```console
    cowpy "Hello Nextflow" -c tux
    ```

    ```console
    __________________
    < Hello Nextflow >
    ------------------
      \
        \
            .--.
          |o_o |
          |:_/ |
          //   \ \
        (|     | )
        /'\_   _/`\
        \___)=(___/
    ```

Se você trabalhou no curso Hello Nextflow, já viu essa ferramenta em ação.
Se não, não se preocupe; vamos cobrir tudo que você precisa saber à medida que avançamos.

### 3.1. Importar o processo e examinar o código

Fornecemos um módulo de processo pré-escrito chamado `COWPY` que encapsula a ferramenta `cowpy`, então você só precisa adicionar uma instrução include antes do bloco workflow.

Faça a seguinte edição no fluxo de trabalho:

=== "Depois"

    ```groovy title="main.nf" linenums="1" hl_lines="4"
    #!/usr/bin/env nextflow

    include { IDENTIFY_LANGUAGE } from './modules/langid.nf'
    include { COWPY } from './modules/cowpy.nf'

    workflow {
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="1"
    #!/usr/bin/env nextflow

    include { IDENTIFY_LANGUAGE } from './modules/langid.nf'

    workflow {
    ```

Você pode abrir o arquivo do módulo para examinar seu código:

```groovy title="modules/cowpy.nf" linenums="1"
#!/usr/bin/env nextflow

// Gera arte ASCII com cowpy
process COWPY {

    publishDir "results/", mode: 'copy'

    container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'

    input:
    path input_file
    val character

    output:
    path "cowpy-${input_file}"

    script:
    """
    cat ${input_file} | cowpy -c ${character} > cowpy-${input_file}
    """
}
```

Como você pode ver, este processo está atualmente projetado para receber um arquivo de entrada (contendo o texto a ser exibido) e um valor especificando o personagem que deve ser desenhado em arte ASCII, geralmente fornecido no nível do fluxo de trabalho por um parâmetro de linha de comando.

### 3.2. Passar um campo do meta map como entrada

Quando usamos a ferramenta `cowpy` no curso Hello Nextflow, usamos um parâmetro de linha de comando para determinar qual personagem usar para desenhar a imagem final.
Isso fazia sentido, porque estávamos gerando apenas uma imagem por execução do pipeline.

No entanto, neste tutorial, queremos gerar uma imagem apropriada para cada sujeito que estamos processando, então usar um parâmetro de linha de comando seria muito limitante.

Boas notícias: temos uma coluna `character` em nossa planilha e, portanto, em nosso meta map.
Vamos usar isso para definir o personagem que o processo deve usar para cada entrada.

Para isso, precisaremos fazer três coisas:

1. Dar um nome ao canal de saída do processo anterior para que possamos operar nele de forma mais conveniente.
2. Determinar como acessar as informações de interesse
3. Adicionar uma chamada ao segundo processo e fornecer as informações adequadamente.

Vamos começar.

#### 3.2.1. Nomear o canal de saída anterior

Aplicamos as manipulações anteriores diretamente no canal de saída do primeiro processo, `IDENTIFY_LANGUAGE.out`.
Para alimentar o conteúdo desse canal para o próximo processo (e fazê-lo de uma forma clara e fácil de ler), queremos dar a ele seu próprio nome, `ch_languages`.

Podemos fazer isso usando o operador [`set`](https://www.nextflow.io/docs/latest/reference/operator.html#set).

No fluxo de trabalho principal, substitua o operador `.view()` por `.set { ch_languages }`, e adicione uma linha testando que podemos nos referir ao canal pelo nome.

=== "Depois"

    ```groovy title="main.nf" linenums="14" hl_lines="19 21 22"
        // Executa langid para identificar o idioma de cada saudação
        IDENTIFY_LANGUAGE(ch_datasheet)
        IDENTIFY_LANGUAGE.out
            .map { meta, file, lang_id ->
                [meta + [lang: lang_id], file]
            }
            .map { meta, file ->

                def lang_group = "unknown"
                if (meta.lang.equals("de") || meta.lang.equals('en')) {
                    lang_group = "germanic"
                }
                else if (meta.lang in ["fr", "es", "it"]) {
                    lang_group = "romance"
                }

                [meta + [lang_group: lang_group], file]
            }
            .set { ch_languages }

        // Temporário: inspecionar ch_languages
        ch_languages.view()
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="14" hl_lines="19"
        // Executa langid para identificar o idioma de cada saudação
        IDENTIFY_LANGUAGE(ch_datasheet)
        IDENTIFY_LANGUAGE.out
            .map { meta, file, lang_id ->
                [meta + [lang: lang_id], file]
            }
            .map { meta, file ->

                def lang_group = "unknown"
                if (meta.lang.equals("de") || meta.lang.equals('en')) {
                    lang_group = "germanic"
                }
                else if (meta.lang in ["fr", "es", "it"]) {
                    lang_group = "romance"
                }

                [meta + [lang_group: lang_group], file]
            }
            .view()
    ```

Vamos executar isso:

```bash
nextflow run main.nf
```

??? success "Saída do comando"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `./main.nf` [friendly_austin] DSL2 - revision: 3dbe460fd6

    [36/cca6a7] IDENTIFY_LANGUAGE (7) | 7 of 7 ✔
    [[id:sampleB, character:tux, lang:de, lang_group:germanic], /workspaces/training/side-quests/metadata/work/e2/6db2402d83cf72081bcd2d11784714/guten_tag.txt]
    [[id:sampleA, character:squirrel, lang:fr, lang_group:romance], /workspaces/training/side-quests/metadata/work/6c/114c818317d169457d6e7336d5d55b/bonjour.txt]
    [[id:sampleC, character:sheep, lang:de, lang_group:germanic], /workspaces/training/side-quests/metadata/work/55/68c69c5efb527f3604ddb3daab8057/hallo.txt]
    [[id:sampleD, character:turkey, lang:en, lang_group:germanic], /workspaces/training/side-quests/metadata/work/2a/4752055ccb5d1370b0ef9da41d3993/hello.txt]
    [[id:sampleE, character:stegosaurus, lang:es, lang_group:romance], /workspaces/training/side-quests/metadata/work/f4/fcd3186dc666d5d239ffa6c37d125d/hola.txt]
    [[id:sampleF, character:moose, lang:fr, lang_group:romance], /workspaces/training/side-quests/metadata/work/c3/3b2627f733f278a7088332a5806108/salut.txt]
    [[id:sampleG, character:turtle, lang:it, lang_group:romance], /workspaces/training/side-quests/metadata/work/36/cca6a7dbfa26ac24f9329787a32e9d/ciao.txt]
    ```

Isso confirma que agora podemos nos referir ao canal pelo nome.

#### 3.2.2. Acessar o arquivo e os metadados do personagem

Sabemos, ao olhar o código do módulo, que o processo `COWPY` espera receber um arquivo de texto e um valor `character`.
Para escrever a chamada ao processo `COWPY`, só precisamos saber como extrair o objeto de arquivo correspondente e os metadados de cada elemento no canal.

Como é frequentemente o caso, a maneira mais simples de fazer isso é usar uma operação `map`.

Nosso canal contém tuplas estruturadas como `[meta, file]`, então podemos acessar o objeto `file` diretamente, e podemos acessar o valor `character` armazenado dentro do meta map referindo-se a ele como `meta.character`.

No fluxo de trabalho principal, faça as seguintes alterações no código:

=== "Depois"

    ```groovy title="main.nf" linenums="34"
        // Temporário: acessar o arquivo e o personagem
        ch_languages.map { meta, file -> file }.view { file -> "File: " + file }
        ch_languages.map { meta, file -> meta.character }.view { character -> "Character: " + character }
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="34"
        // Temporário: inspecionar ch_languages
        ch_languages.view()
    ```

Observe que estamos usando closures (como `{ file -> "File: " + file }`) para tornar a saída das operações `.view` mais legível.

Vamos executar isso:

```bash
nextflow run main.nf -resume
```

??? success "Saída do comando"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `./main.nf` [cheesy_cantor] DSL2 - revision: 15af9c1ec7

    [43/05df08] IDENTIFY_LANGUAGE (7) [100%] 7 of 7, cached: 7 ✔
    Character: squirrel
    File: /workspaces/training/side-quests/metadata/work/8d/4b9498bbccb7a74f04e41877cdc3e5/bonjour.txt
    File: /workspaces/training/side-quests/metadata/work/d3/604274985406e40d79021dea658e60/guten_tag.txt
    Character: tux
    Character: turkey
    File: /workspaces/training/side-quests/metadata/work/d4/fafcc9415b61d2b0fea872e6a05e8a/hello.txt
    File: /workspaces/training/side-quests/metadata/work/02/468ac9efb27f636715e8144b37e9a7/hallo.txt
    Character: sheep
    Character: moose
    Character: stegosaurus
    File: /workspaces/training/side-quests/metadata/work/d4/61a7e1188b4f2742bc72004e226eca/salut.txt
    File: /workspaces/training/side-quests/metadata/work/ae/68364be238c11149c588bf6fc858b1/hola.txt
    File: /workspaces/training/side-quests/metadata/work/43/05df081af5d879ab52e5828fa0357e/ciao.txt
    Character: turtle
    ```

_Os caminhos dos arquivos e os valores de personagem podem aparecer em uma ordem diferente na sua saída._

Isso confirma que somos capazes de acessar o arquivo e o personagem para cada elemento no canal.

#### 3.2.3. Chamar o processo `COWPY`

Agora vamos juntar tudo e de fato chamar o processo `COWPY` no canal `ch_languages`.

No fluxo de trabalho principal, faça as seguintes alterações no código:

=== "Depois"

    ```groovy title="main.nf" linenums="34"
        // Executa cowpy para gerar arte ASCII
        COWPY(
            ch_languages.map { meta, file -> file },
            ch_languages.map { meta, file -> meta.character }
        )
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="34"
        // Temporário: acessar o arquivo e o personagem
        ch_languages.map { meta, file -> [file, meta.character] }
            .view()
    ```

Você pode ver que simplesmente copiamos as duas operações map (sem as instruções `.view()`) como entradas para a chamada do processo.
Só certifique-se de não esquecer a vírgula entre elas!

É um pouco desajeitado, mas veremos como melhorar isso na próxima seção.

Vamos executar isso:

```bash
nextflow run main.nf -resume
```

??? success "Saída do comando"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [suspicious_crick] DSL2 - revision: 25541014c5

    executor >  local (7)
    [43/05df08] IDENTIFY_LANGUAGE (7) [100%] 7 of 7, cached: 7 ✔
    [e7/317c18] COWPY (6)             [100%] 7 of 7 ✔
    ```

Se você olhar no diretório de resultados, deverá ver os arquivos individuais contendo a arte ASCII de cada saudação dita pelo personagem correspondente.

??? abstract "Conteúdo do diretório e exemplo de arquivo"

    ```console
    results/
    ├── cowpy-bonjour.txt
    ├── cowpy-ciao.txt
    ├── cowpy-guten_tag.txt
    ├── cowpy-hallo.txt
    ├── cowpy-hello.txt
    ├── cowpy-hola.txt
    └── cowpy-salut.txt
    ```

    ```text title="results/cowpy-bonjour.txt"
     _________________
    / Bonjour         \
    \ Salut, à demain /
    -----------------
      \
        \
                      _ _
          | \__/|  .~    ~.
          /oo `./      .'
          {o__,   \    {
            / .  . )    \
            `-` '-' \    }
          .(   _(   )_.'
          '---.~_ _ _|
    ```

Isso mostra que conseguimos usar as informações no meta map para parametrizar o comando na segunda etapa do pipeline.

No entanto, como observado acima, parte do código envolvido foi um pouco desajeitada, já que tivemos que desempacotar os metadados ainda no contexto do corpo do fluxo de trabalho.
Essa abordagem funciona bem para usar um pequeno número de campos do meta map, mas escalaria mal se quiséssemos usar muito mais.

Existe outro operador chamado `multiMap()` que nos permite simplificar um pouco isso, mas mesmo assim não é ideal.

??? info "(Opcional) Versão alternativa com `multiMap()`"

    Caso você esteja se perguntando, não poderíamos simplesmente escrever uma única operação `map()` que produz tanto o `file` quanto o `character`, porque isso os retornaria como uma tupla.
    Tivemos que escrever duas operações `map()` separadas para alimentar os elementos `file` e `character` ao processo separadamente.

    Tecnicamente existe outra maneira de fazer isso através de uma única operação de mapeamento, usando o operador `multiMap()`, que é capaz de emitir múltiplos canais.
    Por exemplo, você poderia substituir a chamada ao `COWPY` acima pelo seguinte código:

    === "Depois"

        ```groovy title="main.nf" linenums="34"
            // Executa cowpy para gerar arte ASCII
            COWPY(
                ch_languages.multiMap { meta, file ->
                    file: file
                    character: meta.character
                }
            )
        ```

    === "Antes"

        ```groovy title="main.nf" linenums="34"
            // Executa cowpy para gerar arte ASCII
            COWPY(
                ch_languages.map { meta, file -> file },
                ch_languages.map { meta, file -> meta.character }
            )
        ```

    Isso produz exatamente o mesmo resultado.

Em qualquer caso, é desajeitado ter que fazer algum desempacotamento no nível do fluxo de trabalho.

Seria melhor se pudéssemos passar o meta map inteiro para o processo e selecionar o que precisamos lá dentro.

### 3.3. Passar e usar o meta map completo

O objetivo do meta map é, afinal, passar todos os metadados juntos como um pacote.
A única razão pela qual não pudemos fazer isso acima é que o processo não está configurado para aceitar um meta map.
Mas como controlamos o código do processo, podemos mudar isso.

Vamos modificar o processo `COWPY` para aceitar a estrutura de tupla `[meta, file]` que usamos no primeiro processo, para que possamos simplificar o fluxo de trabalho.

Para isso, precisaremos fazer três coisas:

1. Modificar as definições de entrada do módulo do processo `COWPY`
2. Atualizar o comando do processo para usar o meta map
3. Atualizar a chamada do processo no corpo do fluxo de trabalho

Pronto? Vamos lá!

#### 3.3.1. Modificar a entrada do módulo `COWPY`

Faça as seguintes edições no arquivo de módulo `cowpy.nf`:

=== "Depois"

    ```groovy title="cowpy.nf" linenums="10" hl_lines="2"
    input:
    tuple val(meta), path(input_file)
    ```

=== "Antes"

    ```groovy title="cowpy.nf" linenums="10" hl_lines="2-3"
    input:
    path(input_file)
    val character
    ```

Isso nos permite usar a estrutura de tupla `[meta, file]` que abordamos anteriormente no tutorial.

Observe que não atualizamos a definição de saída do processo para produzir o meta map, a fim de manter o tutorial simplificado, mas sinta-se à vontade para fazer isso você mesmo como exercício seguindo o modelo do processo `IDENTIFY_LANGUAGE`.

#### 3.3.2. Atualizar o comando para usar o campo do meta map

O meta map completo agora está disponível dentro do processo, então podemos nos referir às informações que ele contém diretamente de dentro do bloco de comando.

Faça as seguintes edições no arquivo de módulo `cowpy.nf`:

=== "Depois"

    ```groovy title="cowpy.nf" linenums="16" hl_lines="3"
    script:
    """
    cat ${input_file} | cowpy -c ${meta.character} > cowpy-${input_file}
    """
    ```

=== "Antes"

    ```groovy title="cowpy.nf" linenums="16" hl_lines="3"
    script:
    """
    cat ${input_file} | cowpy -c ${character} > cowpy-${input_file}
    """
    ```

Substituímos a referência ao valor `character` anteriormente passado como entrada independente pelo valor armazenado no meta map, ao qual nos referimos usando `meta.character`.

Agora vamos atualizar a chamada do processo de acordo.

#### 3.3.3. Atualizar a chamada do processo e executá-lo

O processo agora espera que sua entrada use a estrutura de tupla `[meta, file]`, que é o que o processo anterior produz, então podemos simplesmente passar o canal `ch_languages` inteiro para o processo `COWPY`.

Faça as seguintes edições no fluxo de trabalho principal:

=== "Depois"

    ```groovy title="main.nf" linenums="34" hl_lines="2"
    // Executa cowpy para gerar arte ASCII
    COWPY(ch_languages)
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="34" hl_lines="3-4"
    // Executa cowpy para gerar arte ASCII
    COWPY(
        ch_languages.map { meta, file -> file },
        ch_languages.map { meta, file -> meta.character }
    )
    ```

Isso simplifica significativamente a chamada!

Vamos apagar os resultados da execução anterior e executar:

```bash
rm -r results
nextflow run main.nf
```

??? success "Saída do comando"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [wise_sammet] DSL2 - revision: 99797b1e92

    executor >  local (14)
    [5d/dffd4e] process > IDENTIFY_LANGUAGE (7) [100%] 7 of 7 ✔
    [25/9243df] process > COWPY (7)             [100%] 7 of 7 ✔
    ```

Se você olhar no diretório de resultados, deverá ver as mesmas saídas de antes, _ou seja_, arquivos individuais contendo a arte ASCII de cada saudação dita pelo personagem correspondente.

??? abstract "Conteúdo do diretório"

    ```console
    ./results/
    ├── cowpy-bonjour.txt
    ├── cowpy-ciao.txt
    ├── cowpy-guten_tag.txt
    ├── cowpy-hallo.txt
    ├── cowpy-hello.txt
    ├── cowpy-hola.txt
    └── cowpy-salut.txt
    ```

Então isso produz os mesmos resultados que antes com um código mais simples.

Claro, isso pressupõe que você seja capaz de modificar o código do processo.
Em alguns casos, você pode ter que depender de processos existentes que não está em posição de modificar, o que limita suas opções.
A boa notícia, se você planeja usar módulos do projeto [nf-core](https://nf-co.re/), é que os módulos nf-core são todos configurados para usar a estrutura de tupla `[meta, file]` como padrão.

### 3.4. Solução de problemas com entradas obrigatórias ausentes

O valor `character` é obrigatório para que o processo `COWPY` seja executado com sucesso.
Se não definirmos um valor padrão para ele em um arquivo de configuração, DEVEMOS fornecer um valor para ele na planilha.

**O que acontece se não fornecermos?**
Depende do que a planilha de entrada contém e qual versão do fluxo de trabalho estamos executando.

#### 3.4.1. A coluna character existe mas está vazia

Digamos que deletamos o valor do personagem para uma das entradas em nossa planilha para simular um erro de coleta de dados:

```csv title="datasheet.csv" linenums="1" hl_lines="2"
id,character,recording
sampleA,,/workspaces/training/side-quests/metadata/data/bonjour.txt
sampleB,tux,/workspaces/training/side-quests/metadata/data/guten_tag.txt
sampleC,sheep,/workspaces/training/side-quests/metadata/data/hallo.txt
sampleD,turkey,/workspaces/training/side-quests/metadata/data/hello.txt
sampleE,stegosaurus,/workspaces/training/side-quests/metadata/data/hola.txt
sampleF,moose,/workspaces/training/side-quests/metadata/data/salut.txt
sampleG,turtle,/workspaces/training/side-quests/metadata/data/ciao.txt
```

Para qualquer versão do fluxo de trabalho que usamos acima, a chave `character` será criada para todas as entradas quando a planilha for lida, mas para `sampleA` o valor será uma string vazia.

Isso causará um erro.

??? failure "Saída do comando"

    ```console hl_lines="8 11 16 28"
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [marvelous_hirsch] DSL2 - revision: 0dfeee3cc1

    executor >  local (9)
    [c1/c5dd4f] process > IDENTIFY_LANGUAGE (7) [ 85%] 6 of 7
    [d3/b7c415] process > COWPY (2)             [  0%] 0 of 6
    ERROR ~ Error executing process > 'COWPY (1)'

    Caused by:
      Process `COWPY (1)` terminated with an error exit status (2)


    Command executed:

      cat bonjour.txt | cowpy -c  > cowpy-bonjour.txt

    Command exit status:
      2

    Command output:
      (empty)

    Command error:
      usage: cowpy [-h] [-l] [-L] [-t] [-u] [-e EYES] [-c COWACTER] [-E] [-r] [-x]
                  [-C]
                  [msg ...]
      cowpy: error: argument -c/--cowacter: expected one argument

    Work dir:
      /workspaces/training/side-quests/metadata/work/ca/9d49796612a54dec5ed466063c809b

    Container:
      community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273

    Tip: you can try to figure out what's wrong by changing to the process work dir and showing the script file named `.command.sh`

    -- Check '.nextflow.log' file for details
    ```

Quando o Nextflow executa a linha de comando `cowpy` para essa amostra, `${meta.character}` é preenchido com uma string vazia na linha de comando do `cowpy`, então a ferramenta `cowpy` lança um erro dizendo que nenhum valor foi fornecido para o argumento `-c`.

#### 3.4.2. A coluna character não existe na planilha

Agora digamos que deletamos a coluna `character` inteiramente da nossa planilha:

```csv title="datasheet.csv" linenums="1"
id,recording
sampleA,/workspaces/training/side-quests/metadata/data/bonjour.txt
sampleB,/workspaces/training/side-quests/metadata/data/guten_tag.txt
sampleC,/workspaces/training/side-quests/metadata/data/hallo.txt
sampleD,/workspaces/training/side-quests/metadata/data/hello.txt
sampleE,/workspaces/training/side-quests/metadata/data/hola.txt
sampleF,/workspaces/training/side-quests/metadata/data/salut.txt
sampleG,/workspaces/training/side-quests/metadata/data/ciao.txt
```

Neste caso, a chave `character` não será criada quando a planilha for lida.

##### 3.4.2.1. Valor acessado no nível do fluxo de trabalho

Se estamos usando a versão do código que escrevemos na seção 3.2, o Nextflow tentará acessar a chave `character` no meta map ANTES de chamar o processo `COWPY`.

Ele não encontrará nenhum elemento que corresponda à instrução, então não executará `COWPY` de forma alguma.

??? success "Saída do comando"

    ```console hl_lines="7"
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [desperate_montalcini] DSL2 - revision: 0dfeee3cc1

    executor >  local (7)
    [1a/df2544] process > IDENTIFY_LANGUAGE (7) [100%] 7 of 7 ✔
    [-        ] process > COWPY                 -
    ```

Do ponto de vista do Nextflow, este fluxo de trabalho foi executado com sucesso!
No entanto, nenhuma das saídas que queremos será produzida.

##### 3.4.2.2. Valor acessado no nível do processo

Se estamos usando a versão da seção 3.3, o Nextflow passará o meta map inteiro para o processo `COWPY` e tentará executar o comando.

Isso causará um erro, mas diferente do primeiro caso.

??? failure "Saída do comando"

    ```console hl_lines="8 11 16"
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [jovial_bohr] DSL2 - revision: eaaf375827

    executor >  local (9)
    [0d/ada9db] process > IDENTIFY_LANGUAGE (5) [ 85%] 6 of 7
    [06/28065f] process > COWPY (2)             [  0%] 0 of 6
    ERROR ~ Error executing process > 'COWPY (2)'

    Caused by:
      Process `COWPY (2)` terminated with an error exit status (1)


    Command executed:

      cat guten_tag.txt | cowpy -c null > cowpy-guten_tag.txt

    Command exit status:
      1

    Command output:
      (empty)

    Command error:
      Traceback (most recent call last):
        File "/opt/conda/bin/cowpy", line 10, in <module>
          sys.exit(main())
                  ~~~~^^
        File "/opt/conda/lib/python3.13/site-packages/cowpy/cow.py", line 1215, in main
          print(cow(eyes=args.eyes,
                ~~~^^^^^^^^^^^^^^^^
                tongue=args.tongue,
                ^^^^^^^^^^^^^^^^^^^
                thoughts=args.thoughts
                ^^^^^^^^^^^^^^^^^^^^^^
                    ).milk(msg)
                    ^
      TypeError: 'str' object is not callable

    Work dir:
      /workspaces/training/side-quests/metadata/work/06/28065f7d9fd7d22bba084aa941b6d6

    Container:
      community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273

    Tip: you can replicate the issue by changing to the process work dir and entering the command `bash .command.run`

    -- Check '.nextflow.log' file for details
    ```

Isso acontece porque `meta.character` não existe, então nossa tentativa de acessá-lo retorna `null`. Como resultado, o Nextflow literalmente insere `null` na linha de comando, que obviamente não é reconhecido pela ferramenta `cowpy`.

#### 3.4.3. Soluções

Além de fornecer um valor padrão como parte da configuração do fluxo de trabalho, há duas coisas que podemos fazer para lidar com isso de forma mais robusta:

1. Implementar validação de entrada no seu fluxo de trabalho para garantir que a planilha contenha todas as informações necessárias. Você pode encontrar uma [introdução à validação de entrada](../hello_nf-core/05_input_validation.md) no curso de treinamento Hello nf-core. <!-- TODO (future) pending a proper Validation side quest -->

2. Se você quiser garantir que qualquer pessoa que use seu módulo de processo possa identificar imediatamente as entradas obrigatórias, você também pode tornar a propriedade de metadados obrigatória uma entrada explícita.

Aqui está um exemplo de como isso funcionaria.

Primeiro, no nível do processo, atualize a definição de entrada da seguinte forma:

=== "Depois"

    ```groovy title="cowpy.nf" linenums="12" hl_lines="2"
        input:
        tuple val(meta), val(character), path(input_file)
    ```

=== "Antes"

    ```groovy title="cowpy.nf" linenums="12" hl_lines="2"
        input:
        tuple val(meta), path(input_file)
    ```

Em seguida, no nível do fluxo de trabalho, use uma operação de mapeamento para extrair a propriedade `character` dos metadados e torná-la um componente explícito da tupla de entrada:

=== "Depois"

    ```groovy title="main.nf" linenums="37" hl_lines="1"
        COWPY(ch_languages.map{meta, file -> [meta, meta.character, file]})
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="37" hl_lines="1"
        COWPY(ch_languages)
    ```

Essa abordagem tem a vantagem de mostrar explicitamente que `character` é obrigatório, e torna o processo mais fácil de reimplantar em outros contextos.

Isso destaca um princípio de design importante:

**Use o meta map para informações opcionais e descritivas, mas extraia os valores obrigatórios como entradas explícitas.**

O meta map é excelente para manter as estruturas de canal limpas e evitar estruturas de canal arbitrárias, mas para elementos obrigatórios que são diretamente referenciados em um processo, extraí-los como entradas explícitas cria um código mais robusto e fácil de manter.

### Conclusão

Nesta seção, você aprendeu como utilizar metadados para personalizar a execução de um processo, acessando-os tanto no nível do fluxo de trabalho quanto no nível do processo.

---

## Exercício suplementar

Se você quiser praticar o uso de informações do meta map dentro de um processo, tente usar outras informações do meta map, como `lang` e `lang_group`, para personalizar como as saídas são nomeadas e/ou organizadas.

Por exemplo, tente modificar o código para produzir este resultado:

```console title="Results directory contents"
results/
├── germanic
│   ├── de-guten_tag.txt
│   ├── de-hallo.txt
│   └── en-hello.txt
└── romance
    ├── es-hola.txt
    ├── fr-bonjour.txt
    ├── fr-salut.txt
    └── it-ciao.txt
```

<!-- TODO (future) Provide worked out solution -->
<!-- the renaming should use the meta inside the process -->
<!-- the output org should use the meta in the workflow outputs -->

---

## Resumo

Nesta missão paralela, você explorou como trabalhar efetivamente com metadados em fluxos de trabalho Nextflow.

Esse padrão de manter os metadados explícitos e associados aos dados é uma prática recomendada fundamental no Nextflow, oferecendo várias vantagens em relação à codificação direta de informações de arquivos:

- Os metadados dos arquivos permanecem associados aos arquivos ao longo do fluxo de trabalho
- O comportamento dos processos pode ser personalizado por arquivo
- A organização das saídas pode refletir os metadados dos arquivos
- As informações dos arquivos podem ser expandidas durante a execução do pipeline

Aplicar esse padrão no seu próprio trabalho permitirá que você construa fluxos de trabalho de bioinformática robustos e fáceis de manter.

### Padrões principais

1.  **Leitura e estruturação de metadados:** Leitura de arquivos CSV e criação de meta maps organizados que permanecem associados aos seus arquivos de dados.

    ```groovy
    channel.fromPath('datasheet.csv')
      .splitCsv(header: true)
      .map { row ->
          [ [id:row.id, character:row.character], row.recording ]
      }
    ```

2.  **Expansão de metadados durante o fluxo de trabalho:** Adição de novas informações aos seus metadados à medida que o pipeline avança, adicionando saídas de processos e derivando valores por meio de lógica condicional.

    - Adicionando novas chaves com base na saída de um processo

    ```groovy
    .map { meta, file, lang ->
      [ meta + [lang:lang], file ]
    }
    ```

    - Adicionando novas chaves usando uma cláusula condicional

    ```groovy
    .map{ meta, file ->
        if ( meta.lang.equals("de") || meta.lang.equals('en') ){
            lang_group = "germanic"
        } else if ( meta.lang in ["fr", "es", "it"] ) {
            lang_group = "romance"
        } else {
            lang_group = "unknown"
        }
    }
    ```

3.  **Personalização do comportamento dos processos:** Usando metadados dentro do processo.

    ```groovy
    cat $input_file | cowpy -c ${meta.character} > cowpy-${input_file}
    ```

### Recursos adicionais

- [map](https://www.nextflow.io/docs/latest/operator.html#map)
- [stdout](https://www.nextflow.io/docs/latest/process.html#outputs)

---

## O que vem a seguir?

Volte ao [menu de Missões Paralelas](../) ou clique no botão no canto inferior direito da página para avançar para o próximo tópico da lista.
