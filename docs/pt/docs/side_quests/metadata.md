# Metadados e mapas meta

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tradução assistida por IA - [saiba mais e sugira melhorias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Em qualquer análise científica, raramente trabalhamos apenas com os arquivos de dados brutos.
Cada arquivo vem com suas próprias informações adicionais: o que é, de onde veio e o que o torna especial.
Essas informações extras são o que chamamos de metadados.

Metadados são dados que descrevem outros dados.
Metadados rastreiam detalhes importantes sobre arquivos e condições experimentais, e ajudam a adaptar análises às características únicas de cada conjunto de dados.

Pense nisso como um catálogo de biblioteca: enquanto os livros contêm o conteúdo real (dados brutos), as fichas do catálogo fornecem informações essenciais sobre cada livro—quando foi publicado, quem o escreveu, onde encontrá-lo (metadados).
Em pipelines Nextflow, os metadados podem ser usados para:

- Rastrear informações específicas de arquivos ao longo do fluxo de trabalho
- Configurar processos com base nas características dos arquivos
- Agrupar arquivos relacionados para análise conjunta

### Objetivos de aprendizado

Nesta missão secundária, exploraremos como lidar com metadados em fluxos de trabalho.
Começando com uma planilha de dados simples (geralmente chamada de samplesheet em bioinformática) contendo informações básicas de arquivos, você aprenderá como:

- Ler e analisar metadados de arquivos a partir de arquivos CSV
- Criar e manipular mapas de metadados
- Adicionar novos campos de metadados durante a execução do fluxo de trabalho
- Usar metadados para personalizar o comportamento de processos

Essas habilidades ajudarão você a construir pipelines mais robustos e flexíveis que podem lidar com relações complexas de arquivos e requisitos de processamento.

### Pré-requisitos

Antes de assumir esta missão secundária, você deve:

- Ter concluído o tutorial [Hello Nextflow](../hello_nextflow/README.md) ou curso equivalente para iniciantes.
- Estar confortável usando conceitos e mecanismos básicos do Nextflow (processos, canais, operadores)

---

## 0. Começando

#### Abra o codespace de treinamento

Se você ainda não o fez, certifique-se de abrir o ambiente de treinamento conforme descrito em [Configuração do Ambiente](../envsetup/index.md).

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### Mova-se para o diretório do projeto

Vamos mover para o diretório onde os arquivos para este tutorial estão localizados.

```bash
cd side-quests/metadata
```

Você pode configurar o VSCode para focar neste diretório:

```bash
code .
```

#### Revise os materiais

Você encontrará um arquivo de fluxo de trabalho principal e um diretório `data` contendo uma planilha de dados e alguns arquivos de dados.

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

O fluxo de trabalho no arquivo `main.nf` é um esboço que você expandirá gradualmente em um fluxo de trabalho totalmente funcional.

A planilha de dados lista os caminhos para os arquivos de dados e alguns metadados associados, organizados em 3 colunas:

- `id`: autoexplicativo, um ID dado ao arquivo
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

Também forneceremos uma ferramenta de análise de idiomas containerizada chamada `langid`.

#### Revise a tarefa

Seu desafio é escrever um fluxo de trabalho Nextflow que irá:

1. **Identificar** o idioma em cada arquivo automaticamente
2. **Agrupar** arquivos por família de idiomas (idiomas germânicos vs. idiomas românicos)
3. **Personalizar** o processamento para cada arquivo com base em seu idioma e metadados
4. **Organizar** as saídas por grupo de idiomas

Isso representa um padrão típico de fluxo de trabalho onde metadados específicos de arquivos orientam decisões de processamento; exatamente o tipo de problema que os mapas de metadados resolvem elegantemente.

#### Lista de verificação de prontidão

Acha que está pronto para mergulhar?

- [ ] Entendo o objetivo deste curso e seus pré-requisitos
- [ ] Meu codespace está funcionando
- [ ] Defini meu diretório de trabalho adequadamente
- [ ] Entendo a tarefa

Se você pode marcar todas as caixas, está pronto para começar.

---

## 1. Carregar metadados de uma planilha de dados

Abra o arquivo de fluxo de trabalho `main.nf` para examinar o esboço de fluxo de trabalho que estamos fornecendo como ponto de partida.

```groovy title="main.nf" linenums="1"
#!/usr/bin/env nextflow

workflow  {

    ch_datasheet = channel.fromPath("./data/datasheet.csv")

}
```

Você pode ver que configuramos uma fábrica de canal básica para carregar a planilha de dados de exemplo como um arquivo, mas isso ainda não lerá o conteúdo do arquivo.
Vamos começar adicionando isso.

### 1.1. Ler o conteúdo com `splitCsv`

Precisamos escolher um operador que analisará o conteúdo do arquivo adequadamente com o mínimo esforço de nossa parte.
Como nossa planilha de dados está no formato CSV, este é um trabalho para o operador [`splitCsv`](https://www.nextflow.io/docs/latest/reference/operator.html#splitcsv), que carrega cada linha do arquivo como um elemento no canal.

Faça as seguintes alterações para adicionar uma operação `splitCsv()` ao código de construção do canal, mais uma operação `view()` para verificar se o conteúdo do arquivo está sendo carregado corretamente no canal.

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

Observe que estamos usando a opção `header: true` para dizer ao Nextflow para ler a primeira linha do arquivo CSV como a linha de cabeçalho.

Vamos ver o que sai disso, certo?
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

Podemos ver que o operador construiu um mapa de pares chave-valor para cada linha no arquivo CSV, com os cabeçalhos das colunas como chaves para os valores correspondentes.

Cada entrada do mapa corresponde a uma coluna em nossa planilha de dados:

- `id`
- `character`
- `recording`

Isso é ótimo! Facilita o acesso a campos específicos de cada arquivo.
Por exemplo, poderíamos acessar o ID do arquivo com `id` ou o caminho do arquivo txt com `recording`.

??? info "(Opcional) Mais sobre mapas"

    Em Groovy, a linguagem de programação sobre a qual o Nextflow é construído, um mapa é uma estrutura de dados chave-valor semelhante a dicionários em Python, objetos em JavaScript ou hashes em Ruby.

    Aqui está um script executável que mostra como você pode definir um mapa e acessar seu conteúdo na prática:

    ```groovy title="examples/map_demo.nf"
    #!/usr/bin/env nextflow

    // Criar um mapa simples
    def my_map = [id:'sampleA', character:'squirrel']

    // Imprimir o mapa inteiro
    println "map: ${my_map}"

    // Acessar valores individuais usando notação de ponto
    println "id: ${my_map.id}"
    println "character: ${my_map.character}"
    ```

    Mesmo que não tenha um bloco `workflow` adequado, o Nextflow pode executar isso como se fosse um fluxo de trabalho:

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

Digamos que queremos acessar a coluna `character` da planilha de dados e imprimi-la.
Podemos usar o operador `map` do Nextflow para iterar sobre cada item em nosso canal e especificamente selecionar a entrada `character` do objeto mapa.

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

Sucesso! Aproveitamos a estrutura de mapa derivada de nossa planilha de dados para acessar os valores de colunas individuais para cada linha.

Agora que lemos com sucesso a planilha de dados e temos acesso aos dados em cada linha, podemos começar a implementar nossa lógica de pipeline.

### 1.3. Organizar os metadados em um 'mapa meta'

No estado atual do fluxo de trabalho, os arquivos de entrada (sob a chave `recording`) e metadados associados (`id`, `character`) estão todos no mesmo nível, como se estivessem todos em uma grande sacola.
A consequência prática é que todo processo que consome este canal precisaria ser configurado com essa estrutura em mente:

```groovy
    input:
    tuple val(id), val(character), file(recording)
```

Isso é bom desde que o número de colunas na planilha de dados não mude.
No entanto, se você adicionar apenas uma coluna à planilha de dados, a forma do canal não corresponderá mais ao que o processo espera, e o fluxo de trabalho produzirá erros.
Também torna o processo difícil de compartilhar com outras pessoas que podem ter dados de entrada ligeiramente diferentes, e você pode acabar tendo que codificar variáveis no processo que não são necessárias pelo bloco de script.

Para evitar esse problema, precisamos encontrar uma maneira de manter a estrutura do canal consistente independentemente de quantas colunas a planilha de dados contém.

Podemos fazer isso coletando todos os metadados em um item dentro da tupla, que chamaremos de mapa de metadados, ou mais simplesmente 'mapa meta'.

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

Reestruturamos nossos elementos de canal em uma tupla consistindo de dois elementos, o mapa meta e o objeto de arquivo correspondente.

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

Agora, cada elemento no canal contém o mapa meta primeiro e o objeto de arquivo correspondente em segundo:

```console title="Estrutura de saída de exemplo"
[
  [id:sampleA, character:squirrel],
  /workspaces/training/side-quests/metadata/data/bonjour.txt
]
```

Como resultado, adicionar mais colunas na planilha de dados tornará mais metadados disponíveis no mapa `meta`, mas não mudará a forma do canal.
Isso nos permite escrever processos que consomem o canal sem ter que codificar os itens de metadados na especificação de entrada:

```groovy title="Exemplo de sintaxe"
    input:
    tuple val(meta), file(recording)
```

Esta é uma convenção amplamente usada para organizar metadados em fluxos de trabalho Nextflow.

### Conclusão

Nesta seção, você aprendeu:

- **Por que os metadados são importantes:** Manter metadados com seus dados preserva informações importantes sobre arquivos ao longo do fluxo de trabalho.
- **Como ler planilhas de dados:** Usando `splitCsv` para ler arquivos CSV com informações de cabeçalho e transformar linhas em dados estruturados
- **Como criar um mapa meta:** Separando metadados de dados de arquivo usando a estrutura de tupla `[ [id:value, ...], file ]`

---

## 2. Manipulando metadados

Agora que temos nossos metadados carregados, vamos fazer algo com eles!

Vamos usar uma ferramenta chamada [`langid`](https://github.com/saffsd/langid.py) para identificar o idioma contido em cada arquivo de gravação da criatura.
A ferramenta vem pré-treinada em um conjunto de idiomas e, dado um trecho de texto, ela produzirá uma previsão de idioma e uma pontuação de probabilidade associada, ambas para `stdout`.

### 2.1. Importar o processo e examinar o código

Fornecemos a você um módulo de processo pré-escrito chamado `IDENTIFY_LANGUAGE` que envolve a ferramenta `langid`, então você só precisa adicionar uma instrução include antes do bloco workflow.

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

// Usar langid para prever o idioma de cada arquivo de entrada
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

A definição de saída é estruturada como uma tupla com estrutura similar à da entrada, exceto que também contém `stdout` como um terceiro elemento.
Este padrão `tuple val(meta), path(file), <output>` mantém os metadados associados tanto aos dados de entrada quanto às saídas conforme fluem pelo pipeline.

Observe que estamos usando o qualificador de saída [`stdout`](https://www.nextflow.io/docs/latest/process.html#outputs) do Nextflow aqui porque a ferramenta imprime sua saída diretamente no console em vez de escrever um arquivo; e usamos `sed` na linha de comando para remover a pontuação de probabilidade, limpar a string removendo caracteres de nova linha e retornar apenas a previsão de idioma.

### 2.2. Adicionar uma chamada para `IDENTIFY_LANGUAGE`

Agora que o processo está disponível para o fluxo de trabalho, podemos adicionar uma chamada ao processo `IDENTIFY_LANGUAGE` para executá-lo no canal de dados.

Faça as seguintes edições no fluxo de trabalho:

=== "Depois"

    ```groovy title="main.nf" linenums="7" hl_lines="7-9"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                [[id: row.id, character: row.character], row.recording]
            }

        // Executar langid para identificar o idioma de cada saudação
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

Excelente! Agora temos uma previsão de que idioma cada personagem fala.

E conforme observado anteriormente, também incluímos o arquivo de entrada e o mapa meta na saída, o que significa que ambos permanecem associados às novas informações que acabamos de produzir.
Isso se mostrará útil no próximo passo.

!!! note

    De forma mais geral, esse padrão de manter o mapa meta associado aos resultados facilita associar resultados relacionados que compartilham os mesmos identificadores.

    Como você já deve ter aprendido, você não pode confiar na ordem dos itens nos canais para combinar resultados entre eles.
    Em vez disso, você deve usar chaves para associar dados corretamente, e os mapas meta fornecem uma estrutura ideal para esse propósito.

    Exploramos esse caso de uso em detalhes na missão secundária [Dividindo e Agrupando](./splitting_and_grouping.md).

### 2.3. Aumentar metadados com saídas de processo

Dado que os resultados que acabamos de produzir são em si mesmos uma forma de metadados sobre o conteúdo dos arquivos, seria útil adicioná-los ao nosso mapa meta.

No entanto, não queremos modificar o mapa meta existente no local.
Do ponto de vista técnico, é _possível_ fazer isso, mas não é seguro.

Então, em vez disso, criaremos um novo mapa meta contendo o conteúdo do mapa meta existente mais um novo par chave-valor `lang: lang_id` contendo a nova informação, usando o operador `+` (um recurso do Groovy).
E combinaremos isso com uma operação [`map`](https://www.nextflow.io/docs/latest/operator.html#map) para substituir o mapa antigo pelo novo.

Aqui estão as edições que você precisa fazer no fluxo de trabalho:

=== "Depois"

    ```groovy title="main.nf" linenums="13" hl_lines="3-7"
        // Executar langid para identificar o idioma de cada saudação
        IDENTIFY_LANGUAGE(ch_datasheet)
        IDENTIFY_LANGUAGE.out
            .map { meta, file, lang_id ->
                [meta + [lang: lang_id], file]
            }
            .view()
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="13" hl_lines="3"
        // Executar langid para identificar o idioma de cada saudação
        IDENTIFY_LANGUAGE(ch_datasheet)
        IDENTIFY_LANGUAGE.out.view()
    ```

Se você ainda não está familiarizado com o operador `+`, ou se isso parece confuso, reserve alguns minutos para examinar a explicação detalhada abaixo.

??? info "Criação do novo mapa meta usando o operador `+`"

    **Primeiro, você precisa saber que podemos mesclar o conteúdo de dois mapas usando o operador Groovy `+`.**

    Digamos que temos os seguintes mapas:

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

    **Mas e se você precisar adicionar um campo que ainda não faz parte de um mapa?**

    Digamos que você comece novamente de `map1`, mas a previsão de idioma não está em seu próprio mapa (não há `map2`).
    Em vez disso, ela é mantida em uma variável chamada `lang_id`, e você sabe que quer armazenar seu valor (`'fr'`) com a chave `lang`.

    Você pode realmente fazer o seguinte:

    ```groovy
    new_map = [map1 + [lang: lang_id]]
    ```

    Aqui, `[lang: new_info]` cria um novo mapa sem nome dinamicamente, e `map1 + ` mescla `map1` com o novo mapa sem nome, produzindo o mesmo conteúdo `new_map` de antes.

    Legal, não é?

    **Agora vamos transpor isso para o contexto de uma operação `channel.map()` do Nextflow.**

    O código se torna:

    ```groovy
    .map { map1, lang_id ->
        [map1 + [lang: lang_id]]
    }
    ```

    Isso faz o seguinte:

    - `map1, lang_id ->` pega os dois itens na tupla
    - `[map1 + [lang: lang_id]]` cria o novo mapa conforme detalhado acima

    A saída é um único mapa sem nome com o mesmo conteúdo de `new_map` em nosso exemplo acima.
    Então, efetivamente transformamos:

    ```groovy
    [id: 'sampleA', character: 'squirrel'], 'it'
    ```

    em:

    ```groovy
    [id: 'sampleA', character: 'squirrel', lang: 'fr']
    ```

    Esperamos que você possa ver que se mudarmos `map1` para `meta`, isso é basicamente tudo o que precisamos para adicionar a previsão de idioma ao nosso mapa meta em nosso fluxo de trabalho.

    Exceto por uma coisa!

    No caso do nosso fluxo de trabalho, **também precisamos levar em conta a presença do objeto `file` na tupla**, que é composta de `meta, file, lang_id`.

    Então o código aqui se tornaria:

    ```groovy
    .map { meta, file, lang_id ->
        [meta + [lang: lang_id], file]
    }
    ```

    Se você está tendo dificuldade em entender por que o `file` parece estar se movendo na operação `map`, imagine que em vez de `[meta + [lang: lang_id], file]`, essa linha diz `[new_map, file]`.
    Isso deve deixar mais claro que estamos simplesmente deixando o `file` em seu lugar original na segunda posição da tupla. Apenas pegamos o valor `new_info` e o incorporamos no mapa que está na primeira posição.

    **E isso nos traz de volta à estrutura de canal `tuple val(meta), path(file)`!**

Uma vez que você esteja confiante de que entende o que este código está fazendo, execute o fluxo de trabalho para ver se funcionou:

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

Sim, isso está correto!
Reorganizamos ordenadamente a saída do processo de `meta, file, lang_id` para que `lang_id` seja agora uma das chaves no mapa meta, e as tuplas do canal se encaixam no modelo `meta, file` mais uma vez.

### 2.4. Atribuir um grupo de idiomas usando condicionais

Agora que temos nossas previsões de idioma, vamos usar as informações para atribuir alguns novos agrupamentos.

Em nossos dados de exemplo, os idiomas usados por nossos personagens podem ser agrupados em idiomas germânicos (inglês, alemão) e idiomas românicos (francês, espanhol, italiano).
Pode ser útil ter essa classificação prontamente disponível em algum lugar mais adiante no pipeline, então vamos adicionar essa informação no mapa meta.

E, boa notícia, este é mais um caso que se presta perfeitamente ao uso do operador `map`!

Especificamente, vamos definir uma variável chamada `lang_group`, usar alguma lógica condicional simples para determinar qual valor atribuir ao `lang_group` para cada pedaço de dados.

A sintaxe geral vai se parecer com isto:

```groovy
.map { meta, file ->

    // lógica condicional definindo lang_group vai aqui

    [meta + [lang_group: lang_group], file]
}
```

Você pode ver que isso é muito semelhante à operação de mesclagem de mapa dinâmica que usamos no passo anterior.
Só precisamos escrever as declarações condicionais.

Aqui está a lógica condicional que queremos aplicar:

- Definir uma variável chamada `lang_group` com valor padrão `'unknown'`.
- Se `lang` for alemão (`'de'`) ou inglês (`'en'`), mudar `lang_group` para `germanic`.
- Senão, se `lang` estiver incluído em uma lista contendo francês (`'fr'`), espanhol (`'es'`) e italiano (`'it'`), mudar `lang_group` para `romance`.

Tente escrever você mesmo se já souber como escrever declarações condicionais em Nextflow.

!!! tip

    Você pode acessar o valor de `lang` dentro da operação map com `meta.lang`.

Você deve acabar fazendo as seguintes alterações no fluxo de trabalho:

=== "Depois"

    ```groovy title="main.nf" linenums="13" hl_lines="7-19"
        // Executar langid para identificar o idioma de cada saudação
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
        // Executar langid para identificar o idioma de cada saudação
        IDENTIFY_LANGUAGE(ch_datasheet)
        IDENTIFY_LANGUAGE.out
            .map { meta, file, lang_id ->
                [meta + [lang: lang_id], file]
            }
            .view()
    ```

Aqui estão os pontos-chave:

- Usamos `def lang_group = "unknown"` para criar a variável `lang_group` com valor padrão definido como `unknown`.
- Usamos uma estrutura `if {} else if {}` para a lógica condicional, com testes `.equals()` alternativos para os dois idiomas germânicos, e um teste de existência em uma lista para os três idiomas românicos.
- Usamos a operação de mesclagem `meta + [lang_group:lang_group]` como anteriormente para gerar o mapa meta atualizado.

Assim que tudo fizer sentido, execute o fluxo de trabalho novamente para ver o resultado:

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

Como você pode ver, os elementos do canal mantêm sua estrutura `[meta, file]`, mas o mapa meta agora inclui esta nova classificação.

### Conclusão

Nesta seção, você aprendeu como:

- **Aplicar metadados de entrada a canais de saída**: Copiar metadados dessa forma nos permite associar resultados posteriormente com base no conteúdo dos metadados.
- **Criar chaves personalizadas**: Você criou duas novas chaves em seu mapa meta, mesclando-as com `meta + [new_key:value]` no mapa meta existente. Uma baseada em um valor computado de um processo, e uma baseada em uma condição que você definiu no operador `map`.

Isso permite que você associe metadados novos e existentes com arquivos conforme progride através do seu pipeline.
Mesmo que você não esteja usando metadados como parte de um processo, manter o mapa meta associado aos dados dessa forma facilita manter todas as informações relevantes juntas.

---

## 3. Usando informações do mapa meta em um processo

Agora que você sabe como criar e atualizar o mapa meta, podemos chegar à parte realmente divertida: realmente usar os metadados em um processo.

Mais especificamente, vamos adicionar um segundo passo ao nosso fluxo de trabalho para desenhar cada animal como arte ASCII e fazê-lo dizer o texto gravado em um balão de fala.
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

Se você trabalhou no curso Hello Nextflow, você já viu esta ferramenta em ação.
Se não, não se preocupe; cobriremos tudo o que você precisa saber enquanto avançamos.

### 3.1. Importar o processo e examinar o código

Fornecemos a você um módulo de processo pré-escrito chamado `COWPY` que envolve a ferramenta `cowpy`, então você só precisa adicionar uma instrução include antes do bloco workflow.

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

// Gerar arte ASCII com cowpy
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

### 3.2. Passar um campo do mapa meta como entrada

Quando usamos a ferramenta `cowpy` no curso Hello Nextflow, usamos um parâmetro de linha de comando para determinar qual personagem usar para desenhar a imagem final.
Isso fazia sentido, porque estávamos gerando apenas uma imagem por execução do pipeline.

No entanto, neste tutorial, queremos gerar uma imagem apropriada para cada sujeito que estamos processando, então usar um parâmetro de linha de comando seria muito limitante.

Boa notícia: temos uma coluna `character` em nossa planilha de dados e, portanto, em nosso mapa meta.
Vamos usar isso para definir o personagem que o processo deve usar para cada entrada.

Para isso, precisaremos fazer três coisas:

1. Dar um nome ao canal de saída vindo do processo anterior para que possamos operá-lo mais convenientemente.
2. Determinar como acessar a informação de interesse
3. Adicionar uma chamada ao segundo processo e alimentar a informação adequadamente.

Vamos começar.

#### 3.2.1. Nomear o canal de saída anterior

Aplicamos as manipulações anteriores diretamente no canal de saída do primeiro processo, `IDENTIFY_LANGUAGE.out`.
Para alimentar o conteúdo desse canal para o próximo processo (e fazer isso de uma forma que seja clara e fácil de ler), queremos dar a ele seu próprio nome, `ch_languages`.

Podemos fazer isso usando o operador [`set`](https://www.nextflow.io/docs/latest/reference/operator.html#set).

No fluxo de trabalho principal, substitua o operador `.view()` por `.set { ch_languages }`, e adicione uma linha testando que podemos nos referir ao canal pelo nome.

=== "Depois"

    ```groovy title="main.nf" linenums="14" hl_lines="19 21 22"
        // Executar langid para identificar o idioma de cada saudação
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

        // Temporário: espiar em ch_languages
        ch_languages.view()
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="14" hl_lines="19"
        // Executar langid para identificar o idioma de cada saudação
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

Vamos executar isto:

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

Sabemos ao olhar o código do módulo que o processo `COWPY` espera receber um arquivo de texto e um valor `character`.
Para escrever a chamada para o processo `COWPY`, só precisamos saber como extrair o objeto de arquivo correspondente e os metadados de cada elemento no canal.

Como é frequentemente o caso, a maneira mais simples de fazer isso é usar uma operação `map`.

Nosso canal contém tuplas estruturadas como `[meta, file]`, então podemos acessar o objeto `file` diretamente, e podemos acessar o valor `character` armazenado dentro do mapa meta referindo-se a ele como `meta.character`.

No fluxo de trabalho principal, faça as seguintes alterações de código:

=== "Depois"

    ```groovy title="main.nf" linenums="34"
        // Temporário: acessar o arquivo e o personagem
        ch_languages.map { meta, file -> file }.view { file -> "File: " + file }
        ch_languages.map { meta, file -> meta.character }.view { character -> "Character: " + character }
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="34"
        // Temporário: espiar em ch_languages
        ch_languages.view()
    ```

Observe que estamos usando closures (como `{ file -> "File: " + file }`) para tornar a saída das operações `.view` mais legível.

Vamos executar isto:

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

_Os caminhos de arquivo e valores de personagem podem aparecer em uma ordem diferente em sua saída._

Isso confirma que somos capazes de acessar o arquivo e o personagem para cada elemento no canal.

#### 3.2.3. Chamar o processo `COWPY`

Agora vamos juntar tudo e realmente chamar o processo `COWPY` no canal `ch_languages`.

No fluxo de trabalho principal, faça as seguintes alterações de código:

=== "Depois"

    ```groovy title="main.nf" linenums="34"
        // Executar cowpy para gerar arte ASCII
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

Você vê que simplesmente copiamos as duas operações map (menos as declarações `.view()`) como entradas para a chamada do processo.
Apenas certifique-se de não esquecer a vírgula entre elas!

É um pouco desajeitado, mas veremos como melhorar isso na próxima seção.

Vamos executar isto:

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

Se você olhar no diretório de resultados, você deve ver os arquivos individuais contendo a arte ASCII de cada saudação falada pelo personagem correspondente.

??? abstract "Diretório e conteúdo de arquivo de exemplo"

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

Isso mostra que conseguimos usar as informações no mapa meta para parametrizar o comando no segundo passo do pipeline.

No entanto, conforme observado acima, parte do código envolvido foi um pouco desajeitado, já que tivemos que desempacotar metadados ainda no contexto do corpo do fluxo de trabalho.
Essa abordagem funciona bem para usar um pequeno número de campos do mapa meta, mas escalaria mal se quiséssemos usar muito mais.

Existe outro operador chamado `multiMap()` que nos permite simplificar isso um pouco, mas mesmo assim não é ideal.

??? info "(Opcional) Versão alternativa com `multiMap()`"

    Caso você esteja se perguntando, não poderíamos apenas escrever uma única operação `map()` que gera tanto o `file` quanto o `character`, porque isso os retornaria como uma tupla.
    Tivemos que escrever duas operações `map()` separadas para alimentar os elementos `file` e `character` ao processo separadamente.

    Tecnicamente existe outra maneira de fazer isso através de uma única operação de mapeamento, usando o operador `multiMap()`, que é capaz de emitir múltiplos canais.
    Por exemplo, você poderia substituir a chamada para `COWPY` acima pelo seguinte código:

    === "Depois"

        ```groovy title="main.nf" linenums="34"
            // Executar cowpy para gerar arte ASCII
            COWPY(
                ch_languages.multiMap { meta, file ->
                    file: file
                    character: meta.character
                }
            )
        ```

    === "Antes"

        ```groovy title="main.nf" linenums="34"
            // Executar cowpy para gerar arte ASCII
            COWPY(
                ch_languages.map { meta, file -> file },
                ch_languages.map { meta, file -> meta.character }
            )
        ```

    Isso produz exatamente o mesmo resultado.

De qualquer forma, é estranho que tenhamos que fazer algum desempacotamento no nível do fluxo de trabalho.

Seria melhor se pudéssemos alimentar o mapa meta inteiro para o processo e escolher o que precisamos uma vez lá.

### 3.3. Passar e usar o mapa meta inteiro

O ponto do mapa meta é, afinal, passar todos os metadados juntos como um pacote.
A única razão pela qual não pudemos fazer isso acima é que o processo não está configurado para aceitar um mapa meta.
Mas como controlamos o código do processo, podemos mudar isso.

Vamos modificar o processo `COWPY` para aceitar a estrutura de tupla `[meta, file]` que usamos no primeiro processo para que possamos simplificar o fluxo de trabalho.

Para isso, precisaremos fazer três coisas:

1. Modificar as definições de entrada do módulo de processo `COWPY`
2. Atualizar o comando do processo para usar o mapa meta
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

Isso nos permite usar a estrutura de tupla `[meta, file]` que cobrimos anteriormente no tutorial.

Observe que não atualizamos a definição de saída do processo para produzir o mapa meta, a fim de manter o tutorial simplificado, mas sinta-se à vontade para fazer isso você mesmo como exercício seguindo o modelo do processo `IDENTIFY_LANGUAGE`.

#### 3.3.2. Atualizar o comando para usar o campo do mapa meta

O mapa meta inteiro agora está disponível dentro do processo, então podemos nos referir às informações que ele contém diretamente de dentro do bloco de comando.

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

Substituímos a referência ao valor `character` previamente passado como uma entrada independente pelo valor mantido no mapa meta, ao qual nos referimos usando `meta.character`.

Agora vamos atualizar a chamada do processo adequadamente.

#### 3.3.3. Atualizar a chamada do processo e executá-lo

O processo agora espera que sua entrada use a estrutura de tupla `[meta, file]`, que é o que o processo anterior produz, então podemos simplesmente passar o canal `ch_languages` inteiro para o processo `COWPY`.

Faça as seguintes edições no fluxo de trabalho principal:

=== "Depois"

    ```groovy title="main.nf" linenums="34" hl_lines="2"
    // Executar cowpy para gerar arte ASCII
    COWPY(ch_languages)
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="34" hl_lines="3-4"
    // Executar cowpy para gerar arte ASCII
    COWPY(
        ch_languages.map { meta, file -> file },
        ch_languages.map { meta, file
