# Metadados e Meta Maps

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
- Entender por que a interface "meta map + arquivo de dados" é uma convenção amplamente utilizada
- Adicionar novos campos de metadados durante a execução do fluxo de trabalho
- Usar metadados para personalizar o comportamento dos processos e organizar as saídas

Essas habilidades vão ajudá-lo a construir pipelines mais robustos e flexíveis, capazes de lidar com relações complexas entre arquivos e requisitos de processamento.

### Pré-requisitos

Antes de embarcar nesta missão paralela, você deve:

- Ter concluído o tutorial [Hello Nextflow](../../hello_nextflow/index.md) ou um curso equivalente para iniciantes.
- Estar confortável com os conceitos e mecanismos básicos do Nextflow (processos, canais, operadores)

---

## 0. Primeiros passos

#### Abra o codespace de treinamento

Se ainda não tiver feito isso, certifique-se de abrir o ambiente de treinamento conforme descrito em [Configuração do Ambiente](../../envsetup/index.md).

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

O editor abre com o diretório do projeto em foco.

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

Vamos usar uma ferramenta chamada [`COWPY`](https://github.com/jeffbuttars/cowpy) para gerar arte ASCII de cada personagem dizendo a saudação gravada.

??? info "O que o `COWPY` faz?"

    `COWPY` é uma ferramenta de linha de comando que gera arte ASCII para exibir entradas de texto arbitrárias de forma divertida.
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

Além disso, vamos usar uma ferramenta de análise de idiomas chamada `langid` para identificar o idioma que cada personagem fala e organizar as saídas do pipeline de acordo.

#### Revise a tarefa

Seu desafio é escrever um fluxo de trabalho Nextflow que irá:

1. **Gerar arte ASCII** de cada personagem
2. **Organizar** as saídas por família linguística (línguas germânicas vs. românicas)

Isso representa um padrão típico de fluxo de trabalho onde metadados específicos de cada arquivo orientam as decisões de processamento — exatamente o tipo de problema que os meta maps resolvem de forma elegante.

#### Lista de verificação de prontidão

Acha que está pronto para mergulhar de cabeça?

- [ ] Entendo o objetivo deste curso e seus pré-requisitos
- [ ] Meu codespace está funcionando
- [ ] Defini meu diretório de trabalho corretamente
- [ ] Entendo a tarefa

Se você conseguir marcar todas as caixas, pode começar.

---

## 1. Opções básicas para carregar e usar metadados

Abra o arquivo de fluxo de trabalho `main.nf` para examinar o esboço que estamos fornecendo como ponto de partida.

```groovy title="main.nf" linenums="1"
#!/usr/bin/env nextflow

workflow  {
    main:
    ch_datasheet = channel.fromPath("./data/datasheet.csv")
        .splitCsv(header: true)
        .view()

    publish:
    cowpy_art = channel.empty()
}

output {
    cowpy_art {
    }
}
```

O operador [`splitCsv`](https://www.nextflow.io/docs/latest/reference/operator.html#splitcsv) lê cada linha do arquivo como um elemento do canal.
Essa é a mesma abordagem que usamos para carregar dados CSV no Hello Nextflow, nosso curso para iniciantes.
Consulte [esta seção](../../hello_nextflow/02_hello_channels.md#4-read-input-values-from-a-csv-file) se precisar de um lembrete de como isso funciona.

Com `header: true`, a primeira linha é tratada como cabeçalhos de coluna, então cada elemento se torna um map de pares chave-valor indexados pelo nome da coluna.

Observe que, como ainda não estamos executando nenhum processo nos dados, os blocos `publish` e `output` são apenas esboços.

### 1.1. Executar o fluxo de trabalho

Execute o fluxo de trabalho para ver como o conteúdo do canal é estruturado após o carregamento:

```bash
nextflow run main.nf
```

??? success "Saída do comando"

    ```console
     N E X T F L O W   ~  version 25.10.4

    Launching `main.nf` [exotic_albattani] DSL2 - revision: c0d03cec83

    [id:sampleA, character:squirrel, recording:/workspaces/training/side-quests/metadata/data/bonjour.txt]
    [id:sampleB, character:tux, recording:/workspaces/training/side-quests/metadata/data/guten_tag.txt]
    [id:sampleC, character:sheep, recording:/workspaces/training/side-quests/metadata/data/hallo.txt]
    [id:sampleD, character:turkey, recording:/workspaces/training/side-quests/metadata/data/hello.txt]
    [id:sampleE, character:stegosaurus, recording:/workspaces/training/side-quests/metadata/data/hola.txt]
    [id:sampleF, character:moose, recording:/workspaces/training/side-quests/metadata/data/salut.txt]
    [id:sampleG, character:turtle, recording:/workspaces/training/side-quests/metadata/data/ciao.txt]
    ```

Como você pode ver, o operador construiu um map de pares chave-valor para cada linha do arquivo CSV, com os cabeçalhos das colunas como chaves para os valores correspondentes.

Cada entrada do map corresponde a uma coluna em nossa planilha:

- `id`
- `character`
- `recording`

Isso facilita o acesso a campos específicos de cada linha.
Por exemplo, poderíamos acessar o ID do arquivo com `id` ou o caminho do arquivo txt com `recording`.

??? info "(Opcional) Mais sobre maps Groovy"

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
     N E X T F L O W   ~  version 25.10.4

    Launching `map_demo.nf` [cheesy_plateau] DSL2 - revision: fae5b8496e

    map: [id:sampleA, character:squirrel]
    id: sampleA
    character: squirrel
    ```

### 1.2. Selecionar um campo específico com `map`

Vamos usar o operador `map` para iterar sobre cada elemento em um canal e selecionar apenas o campo `character`, que podemos acessar pelo nome usando notação de ponto.

#### 1.2.1. Adicionar a operação map

Para acessar a coluna `character`, adicione a operação `map` antes da operação `.view()` da seguinte forma:

=== "Depois"

    ```groovy title="main.nf" linenums="5" hl_lines="3-5"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                row.character
            }
            .view()
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="5"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .view()
    ```

Essa forma de acessar um campo específico é explicada com mais detalhes em [esta seção](../../hello_nextflow/02_hello_channels.md#43-use-the-map-operator-to-extract-the-greetings) do Hello Nextflow, caso precise de um lembrete.

#### 1.2.2. Executar o fluxo de trabalho

Execute o fluxo de trabalho para verificar que você consegue visualizar os nomes dos personagens extraídos.

```bash
nextflow run main.nf
```

??? success "Saída do comando"

    ```console
     N E X T F L O W   ~  version 25.10.4

    Launching `main.nf` [exotic_albattani] DSL2 - revision: c0d03cec83

    squirrel
    tux
    sheep
    turkey
    stegosaurus
    moose
    turtle
    ```

Isso mostra que conseguimos acessar os valores da coluna `character` para cada linha.

Agora vamos fazer algo com esses dados: usar os campos `character` e `recording` juntos para gerar arte ASCII usando `COWPY`.

### 1.3. Emitir sub-canais com `multiMap`

Fornecemos um módulo de processo pré-escrito chamado `COWPY`, então primeiro você precisa examinar os requisitos de entrada do processo.

Você pode abrir o arquivo para ver como o processo é definido:

```groovy title="modules/cowpy.nf" linenums="1" hl_lines="7 8"
// Gera arte ASCII com cowpy
process COWPY {

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

Como você pode ver, o processo recebe duas entradas separadas: um arquivo de gravação e um nome de personagem.
É importante notar que temos valores para ambos, mas eles estão atualmente agrupados dentro de cada elemento no canal.

Uma forma de extrair múltiplos campos em canais separados é o operador [`multiMap`](https://www.nextflow.io/docs/latest/reference/operator.html#multimap), que divide um canal em múltiplos sub-canais nomeados em uma única operação.

#### 1.3.1. Adicionar a operação multiMap

Substitua a operação `map` por `multiMap`:

=== "Depois"

    ```groovy title="main.nf" linenums="5" hl_lines="3-6"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .multiMap { row ->
                file: row.recording
                character: row.character
            }
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="5" hl_lines="3-6"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                row.character
            }
            .view()
    ```

O bloco `multiMap` define dois sub-canais nomeados (`file` e `character`) a partir de cada linha, que podemos acessar como `ch_datasheet.file` e `ch_datasheet.character`.

#### 1.3.2. Chamar COWPY nos sub-canais

Agora, inclua o processo `COWPY` e forneça cada sub-canal como um argumento separado:

=== "Depois"

    ```groovy title="main.nf" linenums="1" hl_lines="3 14"
    #!/usr/bin/env nextflow

    include { COWPY } from './modules/cowpy.nf'

    workflow {
        main:
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .multiMap { row ->
                file: row.recording
                character: row.character
            }

        COWPY(ch_datasheet.file, ch_datasheet.character)

        publish:
        cowpy_art = channel.empty()
    }

    output {
        cowpy_art {
        }
    }
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="1"
    #!/usr/bin/env nextflow

    workflow {
        main:
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .multiMap { row ->
                file: row.recording
                character: row.character
            }

        publish:
        cowpy_art = channel.empty()
    }

    output {
        cowpy_art {
        }
    }
    ```

Isso nos permite passar os dois campos separadamente, como o `COWPY` requer.

#### 1.3.3. Configurar a publicação das saídas

Por fim, adicione a saída do `COWPY` ao bloco `publish:`:

=== "Depois"

    ```groovy title="main.nf" linenums="14" hl_lines="4"
        COWPY(ch_datasheet.file, ch_datasheet.character)

        publish:
        cowpy_art = COWPY.out
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="14" hl_lines="4"
        COWPY(ch_datasheet.file, ch_datasheet.character)

        publish:
        cowpy_art = channel.empty()
    ```

Isso nos permitirá visualizar facilmente as saídas produzidas pelo fluxo de trabalho.

#### 1.3.4. Executar o fluxo de trabalho

Execute o fluxo de trabalho para verificar que o `COWPY` é executado com as entradas fornecidas:

```bash
nextflow run main.nf
```

??? success "Saída do comando"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [clever_dijkstra] DSL2 - revision: a1b2c3d4e5

    executor >  local (7)
    [3a/f1c290] COWPY (7) [100%] 7 of 7 ✔
    ```

Como você pode ver, o `COWPY` foi executado em cada arquivo usando o personagem correto para cada um.

??? abstract "Conteúdo do diretório de resultados"

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

??? example "Conteúdo de results/cowpy-guten_tag.txt"

    ```console
    $ cat results/cowpy-guten_tag.txt
     _____________________________
    / Guten Tag, wie geht es dir? \
    \ Auf Wiedersehen, bis morgen /
     -----------------------------
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

Essa abordagem funciona, mas tem uma limitação: tivemos que dividir o canal em dois sub-canais separados.
Se quiséssemos passar mais campos para o processo, precisaríamos dividi-los em ainda mais sub-canais.
Isso pode ficar chato e confuso.

Boas notícias: existe uma forma mais simples de fazer isso.

### 1.4. Agrupar tudo como uma única entrada para o processo

Em vez de dividir os campos em canais separados, podemos atualizar o processo para receber todas as entradas como uma única tupla, o que simplifica a chamada ao processo.

#### 1.4.1. Atualizar o processo COWPY

Atualize o `COWPY` para aceitar uma tupla correspondente aos três elementos de cada linha:

=== "Depois"

    ```groovy title="modules/cowpy.nf" linenums="1" hl_lines="7 10 14"
    // Gera arte ASCII com cowpy
    process COWPY {

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'

        input:
        tuple val(id), val(character), path(recording)

        output:
        path "cowpy-${recording}"

        script:
        """
        cat ${recording} | cowpy -c ${character} > cowpy-${recording}
        """
    }
    ```

=== "Antes"

    ```groovy title="modules/cowpy.nf" linenums="1" hl_lines="7-8 11 15"
    // Gera arte ASCII com cowpy
    process COWPY {

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

Agora o processo recebe apenas uma entrada contendo tudo que precisamos fornecer.

#### 1.4.2. Usar `map()` para criar a tupla de entrada

Ainda precisamos usar uma operação de mapeamento para enumerar os elementos que queremos passar na tupla para o processo:

=== "Depois"

    ```groovy title="main.nf" linenums="5" hl_lines="3-5"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                [row.id, row.character, row.recording]
            }
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="5" hl_lines="3-6"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .multiMap { row ->
                file: row.recording
                character: row.character
            }
    ```

Você pode se perguntar por que não podemos simplesmente passar o map Groovy inteiro vindo do `splitCsv` como está.
É porque precisamos informar ao Nextflow explicitamente que o arquivo de gravação precisa ser tratado como um caminho (ou seja, precisa ser preparado corretamente).
Isso acontece no nível da interface de entrada do `COWPY`, onde o elemento `recording` é explicitamente designado como um `path`.

#### 1.4.3. Atualizar a chamada ao processo

Por fim, vamos substituir as duas entradas separadas na chamada ao processo pela única tupla que acabamos de criar:

=== "Depois"

    ```groovy title="main.nf" linenums="5" hl_lines="7"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                [row.id, row.character, row.recording]
            }

        COWPY(ch_datasheet)
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="5" hl_lines="7"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                [row.id, row.character, row.recording]
            }

        COWPY(ch_datasheet.file, ch_datasheet.character)
    ```

Isso simplifica um pouco a chamada ao processo.

#### 1.4.4. Executar o fluxo de trabalho

Execute o fluxo de trabalho para verificar que o `COWPY` ainda consegue processar os dados corretamente:

```bash
nextflow run main.nf
```

??? success "Saída do comando"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [pedantic_lovelace] DSL2 - revision: b2c3d4e5f6

    executor >  local (7)
    [5e/2a1b34] COWPY (7) [100%] 7 of 7 ✔
    ```

A saída são os mesmos sete arquivos `cowpy-*.txt` de antes, agora produzidos com uma chamada mais simples ao `COWPY`.

??? abstract "Conteúdo do diretório de resultados"

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

??? example "Conteúdo de results/cowpy-guten_tag.txt"

    ```console
    $ cat results/cowpy-guten_tag.txt
     _____________________________
    / Guten Tag, wie geht es dir? \
    \ Auf Wiedersehen, bis morgen /
     -----------------------------
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

Isso é uma pequena melhoria em relação à abordagem com `multiMap`.
Mas ainda tivemos que desempacotar o map Groovy original para criar a tupla de entrada, e há um acoplamento estreito entre o processo e a planilha: a definição de entrada do `COWPY` agora referencia diretamente os nomes das colunas `id`, `character` e `recording`.

```groovy
input:
tuple val(id), val(character), path(recording)
```

Se um colaborador usar uma planilha com estrutura diferente — com colunas adicionais ou em ordem diferente — esse processo não funcionará sem modificações.
Isso torna o processo frágil, pois sua estrutura de entrada está vinculada à composição exata da planilha.

Para resolver isso, precisamos de uma forma de passar todos os metadados como um pacote sem codificar sua estrutura exata na interface do processo.

### 1.5. Usar a interface meta map + arquivo

A solução é separar duas preocupações distintas no canal: os **metadados sobre uma amostra** e o **arquivo de dados** em si.
Ao agrupar todos os metadados em um único map — o "meta map" — obtemos uma tupla consistente de dois elementos independentemente de quantas colunas de metadados a planilha contenha:

```groovy title="Syntax example"
input:
tuple val(meta), path(file)
```

Adicionar ou remover colunas da planilha muda o que está dentro de `meta`, mas a forma da tupla `[meta, file]` permanece constante.
Processos que aceitam essa estrutura não precisam saber nem se importar com quantos campos de metadados existem.

#### 1.5.1. Reorganizar o conteúdo da tupla em um meta map

Vamos reestruturar a operação `map` para produzir uma tupla `[meta, file]`:

=== "Depois"

    ```groovy title="main.nf" linenums="5" hl_lines="4 6 8 11"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                [[id: row.id, character: row.character], row.recording]
            }
            .view()

        // COWPY(ch_datasheet)  // Será atualizado no próximo passo

        publish:
        cowpy_art = channel.empty() // COWPY.out
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="5" hl_lines="4 7"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                [row.id, row.character, row.recording]
            }

        COWPY(ch_datasheet)

        publish:
        cowpy_art = COWPY.out
    ```

Você vai notar que também adicionamos uma instrução `view()`, comentamos a chamada ao `COWPY` e substituímos `COWPY.out` por `channel.empty()` porque a definição de entrada do processo ainda não corresponde à nova estrutura.

#### 1.5.2. Executar o fluxo de trabalho para inspecionar o conteúdo reorganizado

Execute o fluxo de trabalho para ver a nova forma do canal:

```bash
nextflow run main.nf
```

??? success "Saída do comando"

    ```console title="View meta map"
     N E X T F L O W   ~  version 25.10.4

    Launching `main.nf` [lethal_booth] DSL2 - revision: 0d8f844c07

    [[id:sampleA, character:squirrel], /workspaces/training/side-quests/metadata/data/bonjour.txt]
    [[id:sampleB, character:tux], /workspaces/training/side-quests/metadata/data/guten_tag.txt]
    [[id:sampleC, character:sheep], /workspaces/training/side-quests/metadata/data/hallo.txt]
    [[id:sampleD, character:turkey], /workspaces/training/side-quests/metadata/data/hello.txt]
    [[id:sampleE, character:stegosaurus], /workspaces/training/side-quests/metadata/data/hola.txt]
    [[id:sampleF, character:moose], /workspaces/training/side-quests/metadata/data/salut.txt]
    [[id:sampleG, character:turtle], /workspaces/training/side-quests/metadata/data/ciao.txt]
    ```

Cada elemento no canal agora é uma tupla de dois elementos: o meta map primeiro e o arquivo em segundo.

```console title="Example element structure"
[
  [id:sampleA, character:squirrel],
  /workspaces/training/side-quests/metadata/data/bonjour.txt
]
```

Se mais tarde adicionarmos uma coluna `language` à planilha, ela ficará disponível como `meta.language` sem exigir nenhuma alteração na definição de entrada do processo.

#### 1.5.3. Atualizar o processo `COWPY` para usar o meta map

Atualize o `COWPY` para aceitar a estrutura de tupla `[meta, file]`:

=== "Depois"

    ```groovy title="modules/cowpy.nf" linenums="1" hl_lines="7 10 14"
    // Gera arte ASCII com cowpy
    process COWPY {

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'

        input:
        tuple val(meta), path(input_file)

        output:
        path "cowpy-${input_file}"

        script:
        """
        cat ${input_file} | cowpy -c ${meta.character} > cowpy-${input_file}
        """
    }
    ```

=== "Antes"

    ```groovy title="modules/cowpy.nf" linenums="1" hl_lines="7 10 14"
    // Gera arte ASCII com cowpy
    process COWPY {

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'

        input:
        tuple val(id), val(character), path(recording)

        output:
        path "cowpy-${recording}"

        script:
        """
        cat ${recording} | cowpy -c ${character} > cowpy-${recording}
        """
    }
    ```

Dentro do bloco de script, `meta.character` acessa o campo `character` do meta map.
Qualquer campo no meta map é acessível da mesma forma.

#### 1.5.4. Atualizar a chamada ao processo

Restaure a chamada ao `COWPY` e conecte sua saída para publicação:

=== "Depois"

    ```groovy title="main.nf" linenums="5" hl_lines="7 10"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                [[id: row.id, character: row.character], row.recording]
            }

        COWPY(ch_datasheet)

        publish:
        cowpy_art = COWPY.out
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="5" hl_lines="6 8 11"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                [[id: row.id, character: row.character], row.recording]
            }
            .view()

        // COWPY(ch_datasheet)  // Será atualizado no próximo passo

        publish:
        cowpy_art = channel.empty() // COWPY.out
    ```

Também restauramos a publicação das saídas.

#### 1.5.5. Executar o fluxo de trabalho

Execute o fluxo de trabalho para verificar que tudo funciona:

```bash
nextflow run main.nf
```

??? success "Saída do comando"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [wise_sammet] DSL2 - revision: 99797b1e92

    executor >  local (7)
    [5d/dffd4e] COWPY (7) [100%] 7 of 7 ✔
    ```

O diretório de resultados agora contém os arquivos de arte ASCII.

??? abstract "Conteúdo do diretório"

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

??? example "Conteúdo de results/cowpy-guten_tag.txt"

    ```console
    $ cat results/cowpy-guten_tag.txt
     _____________________________
    / Guten Tag, wie geht es dir? \
    \ Auf Wiedersehen, bis morgen /
     -----------------------------
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

O processo agora recebe todos os metadados como um pacote via `meta`, usa o que precisa (`meta.character`) e ignora o restante.

Essa é a interface padrão usada por todos os módulos do [nf-core](https://nf-co.re/).
O padrão `tuple val(meta), path(file)` aparece de forma consistente em toda a biblioteca de módulos nf-core, razão pela qual fluxos de trabalho que adotam essa convenção podem incorporar módulos nf-core com mínimo atrito.

### Conclusão

Nesta seção, você aprendeu:

- **Como ler planilhas:** Usando `splitCsv` para analisar arquivos CSV com informações de cabeçalho
- **Por que a convenção do meta map existe:** Separar metadados dos arquivos de dados em tuplas `[meta, file]` mantém a estrutura do canal estável à medida que a planilha evolui
- **Como usar campos do meta map dentro de um processo:** Qualquer campo no meta map é acessível via notação de ponto no bloco de script

---

## 2. Manipulações adicionais de metadados

Agora que a interface do meta map está em vigor, podemos enriquecê-la à medida que os dados fluem pelo pipeline.

Vamos usar uma ferramenta chamada [`langid`](https://github.com/saffsd/langid.py) para identificar o idioma em cada arquivo de gravação.
Dado um trecho de texto, ela produz uma previsão de idioma e uma pontuação de probabilidade para `stdout`.

### 2.1. Adicionar uma etapa de identificação de idioma

Fornecemos um módulo de processo pré-escrito chamado `IDENTIFY_LANGUAGE` que encapsula a ferramenta `langid`.

Abra o arquivo do módulo para examinar seu código:

```groovy title="modules/langid.nf" linenums="1" hl_lines="7 10"
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

A definição de entrada usa a mesma estrutura `tuple val(meta), path(file)` que acabamos de construir na Seção 1, então `ch_datasheet` pode alimentar diretamente esse processo sem nenhuma adaptação.

A saída adiciona `stdout` como terceiro elemento: isso captura a previsão de idioma que o `langid` imprime no console.
O comando `sed` remove a pontuação de probabilidade e a nova linha final, deixando apenas o código de idioma de duas letras.

#### 2.1.1. Adicionar uma chamada ao `IDENTIFY_LANGUAGE`

Inclua o módulo do processo `IDENTIFY_LANGUAGE` e chame-o no canal da planilha:

=== "Depois"

    ```groovy title="main.nf" linenums="1" hl_lines="4 14-16"
    #!/usr/bin/env nextflow

    include { COWPY } from './modules/cowpy.nf'
    include { IDENTIFY_LANGUAGE } from './modules/langid.nf'

    workflow {
        main:
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                [[id: row.id, character: row.character], row.recording]
            }

        // Executa langid para identificar o idioma de cada saudação
        IDENTIFY_LANGUAGE(ch_datasheet)
        IDENTIFY_LANGUAGE.out.view()

        COWPY(ch_datasheet)
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="1"
    #!/usr/bin/env nextflow

    include { COWPY } from './modules/cowpy.nf'

    workflow {
        main:
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                [[id: row.id, character: row.character], row.recording]
            }

        COWPY(ch_datasheet)
    ```

A saída principal desse processo é apenas uma string, então não há arquivos de saída para publicar.
Em vez disso, usamos `IDENTIFY_LANGUAGE.out.view()` para visualizar os resultados da operação.

#### 2.1.2. Executar o fluxo de trabalho

Execute o fluxo de trabalho para produzir a identificação de idioma, usando `-resume` para evitar re-executar as tarefas do `COWPY`:

```bash
nextflow run main.nf -resume
```

??? success "Saída do comando"

    ```console
     N E X T F L O W   ~  version 25.10.4

    Launching `main.nf` [voluminous_mcnulty] DSL2 - revision: f9bcfebabb

    executor >  local (14)
    [5d/dffd4e] COWPY (7)             [100%] 7 of 7, cached: 7 ✔
    [4e/f722fe] IDENTIFY_LANGUAGE (7) [100%] 7 of 7 ✔
    [[id:sampleA, character:squirrel], /workspaces/training/side-quests/metadata/work/eb/f7148ebdd898fbe1136bec6a714acb/bonjour.txt, fr]
    [[id:sampleB, character:tux], /workspaces/training/side-quests/metadata/work/16/71d72410952c22cd0086d9bca03680/guten_tag.txt, de]
    [[id:sampleD, character:turkey], /workspaces/training/side-quests/metadata/work/c4/b7562adddc1cc0b7d414ec45d436eb/hello.txt, en]
    [[id:sampleC, character:sheep], /workspaces/training/side-quests/metadata/work/ea/04f5d979429e4455e14b9242fb3b45/hallo.txt, de]
    [[id:sampleF, character:moose], /workspaces/training/side-quests/metadata/work/5a/6c2b84bf8fadb98e28e216426be079/salut.txt, fr]
    [[id:sampleE, character:stegosaurus], /workspaces/training/side-quests/metadata/work/af/ee7c69bcab891c40d0529305f6b9e7/hola.txt, es]
    [[id:sampleG, character:turtle], /workspaces/training/side-quests/metadata/work/4e/f722fe47271ba7ebcd69afa42964ca/ciao.txt, it]
    ```

Agora temos uma previsão de idioma para cada arquivo no conjunto de dados.

Observe que a tupla de saída é composta por `[meta, file, lang_id]`, o que significa que o meta map e o arquivo são carregados junto com o novo resultado.

!!! note "Nota"

    Esse padrão de manter o meta map associado aos resultados facilita a associação de resultados entre canais posteriormente.
    Não é possível confiar na ordem dos itens nos canais para associar dados corretamente.
    Em vez disso, você deve usar chaves.
    Os meta maps fornecem uma estrutura ideal para esse propósito.

    Esse caso de uso é explorado em detalhes na missão paralela [Splitting & Grouping](../splitting_and_grouping/index.md).

### 2.2. Enriquecer os metadados com saídas de processos

A previsão de idioma é em si um metadado sobre os dados no arquivo.
Em vez de mantê-la como um elemento separado, vamos incorporá-la de volta ao meta map.

#### 2.2.1. Criar um meta map novo e expandido

Podemos criar um novo meta map para substituir o original usando o operador Groovy `+`:

=== "Depois"

    ```groovy title="main.nf" linenums="14" hl_lines="3-7"
        // Executa langid para identificar o idioma de cada saudação
        IDENTIFY_LANGUAGE(ch_datasheet)
        IDENTIFY_LANGUAGE.out
            .map { meta, file, lang_id ->
                [meta + [lang: lang_id], file]
            }
            .view()
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="14" hl_lines="3"
        // Executa langid para identificar o idioma de cada saudação
        IDENTIFY_LANGUAGE(ch_datasheet)
        IDENTIFY_LANGUAGE.out.view()
    ```

O coração dessa operação é `#!groovy meta + [lang: lang_id]`.

Esse código essencialmente cria um map temporário com um único par chave-valor contendo o código de idioma (`[lang: lang_id]`), e então usa o operador Groovy `+` para combiná-lo com o map `meta` original contendo os metadados pré-existentes, produzindo um meta map novo e expandido.

Para uma explicação mais detalhada, veja o box abaixo.

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
    new_map = map1 + [lang: lang_id]
    ```

    Aqui, `[lang: lang_id]` cria um novo map sem nome na hora, e `map1 + ` mescla `map1` com o novo map sem nome, produzindo o mesmo conteúdo de `new_map` que antes.

    Legal, não é?

    **Agora vamos transpor isso para o contexto de uma operação `channel.map()` do Nextflow.**

    O código se torna:

    ```groovy
    .map { map1, lang_id ->
        map1 + [lang: lang_id]
    }
    ```

    Isso faz o seguinte:

    - `#!groovy map1, lang_id ->` recebe os dois itens na tupla
    - `#!groovy map1 + [lang: lang_id]` cria o novo map conforme detalhado acima

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

    Se você está tendo dificuldade em entender por que o `file` parece estar se movendo na operação `map`, imagine que em vez de `#!groovy [meta + [lang: lang_id], file]`, essa linha lê `[new_map, file]`.
    Isso deve deixar mais claro que estamos simplesmente deixando o `file` em seu lugar original na segunda posição na tupla. Apenas pegamos o valor `new_info` e o incorporamos ao map que está na primeira posição.

    **E isso nos traz de volta à estrutura de canal `tuple val(meta), path(file)`!**

#### 2.2.2. Executar o fluxo de trabalho

Quando você estiver confiante de que entendeu o que esse código está fazendo, execute o fluxo de trabalho para ver se funcionou:

```bash
nextflow run main.nf -resume
```

??? success "Saída do comando"

    ```console
     N E X T F L O W   ~  version 25.10.4

    Launching `main.nf` [cheeky_fermat] DSL2 - revision: d096281ee4

    [5d/dffd4e] COWPY (7)             [100%] 7 of 7, cached: 7 ✔
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

!!! tip "Removendo chaves de um meta map"

    Você pode remover uma chave de um meta map usando o método Groovy [`subMap`](https://docs.groovy-lang.org/latest/html/groovy-jdk/java/util/Map.html#subMap(java.util.Collection)), que retorna um novo map contendo apenas as chaves que você especificar:

    ```groovy
    meta.subMap(['id', 'character'])  // retorna um map com apenas 'id' e 'character'
    ```

    Isso é útil quando um processo ou módulo downstream não precisa de todos os campos que se acumularam no meta map.

### 2.3. Atribuir um grupo linguístico usando condicionais

Com a previsão de idioma no meta map, podemos derivar mais metadados a partir dela.
Os idiomas em nosso conjunto de dados se enquadram em duas famílias: germânica (inglês, alemão) e românica (francês, espanhol, italiano).
Adicionar um campo `lang_group` tornará essa classificação disponível mais adiante no pipeline.

#### 2.3.1. Adicionar uma operação `map` com a lógica condicional

Vamos usar uma segunda operação `map` com lógica condicional para atribuir a família linguística:

```groovy
.map { meta, file ->

    // lógica condicional definindo lang_group vai aqui

    [meta + [lang_group: lang_group], file]
}
```

Aqui está a lógica a aplicar:

- Começar com `lang_group = 'unknown'` como padrão.
- Se `meta.lang` for `'de'` ou `'en'`, definir `lang_group` como `'germanic'`.
- Caso contrário, se `meta.lang` estiver em `['fr', 'es', 'it']`, definir `lang_group` como `'romance'`.

!!! tip "Dica"

    Você pode acessar o valor de `lang` dentro da operação map com `meta.lang`.

Faça as seguintes alterações no fluxo de trabalho:

=== "Depois"

    ```groovy title="main.nf" linenums="14" hl_lines="7-19 21"
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

        ch_languages.view()
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="14" hl_lines="7"
        // Executa langid para identificar o idioma de cada saudação
        IDENTIFY_LANGUAGE(ch_datasheet)
        IDENTIFY_LANGUAGE.out
            .map { meta, file, lang_id ->
                [meta + [lang: lang_id], file]
            }
            .view()
    ```

Pontos principais:

- `def lang_group = "unknown"` inicializa a variável com um valor padrão seguro.
- A estrutura `if / else if` trata as duas famílias linguísticas; qualquer outro caso permanece como `'unknown'`.
- `#!groovy .set { ch_languages }` dá um nome ao canal resultante para uso na próxima etapa.

<!-- TODO (future) Add note/links to relevant docs in additional resources section -->

#### 2.3.2. Executar o fluxo de trabalho

Execute o fluxo de trabalho para verificar que funciona:

```bash
nextflow run main.nf -resume
```

??? success "Saída do comando"

    ```console
     N E X T F L O W   ~  version 25.10.4

    Launching `main.nf` [wise_almeida] DSL2 - revision: 46778c3cd0

    [5d/dffd4e] COWPY (7)             [100%] 7 of 7, cached: 7 ✔
    [da/652cc6] IDENTIFY_LANGUAGE (7) [100%] 7 of 7, cached: 7 ✔
    [[id:sampleA, character:squirrel, lang:fr, lang_group:romance], /workspaces/training/side-quests/metadata/data/bonjour.txt]
    [[id:sampleB, character:tux, lang:de, lang_group:germanic], /workspaces/training/side-quests/metadata/data/guten_tag.txt]
    [[id:sampleC, character:sheep, lang:de, lang_group:germanic], /workspaces/training/side-quests/metadata/data/hallo.txt]
    [[id:sampleD, character:turkey, lang:en, lang_group:germanic], /workspaces/training/side-quests/metadata/data/hello.txt]
    [[id:sampleE, character:stegosaurus, lang:es, lang_group:romance], /workspaces/training/side-quests/metadata/data/hola.txt]
    [[id:sampleF, character:moose, lang:fr, lang_group:romance], /workspaces/training/side-quests/metadata/data/salut.txt]
    [[id:sampleG, character:turtle, lang:it, lang_group:romance], /workspaces/training/side-quests/metadata/data/ciao.txt]
    ```

O meta map agora carrega quatro campos: `id`, `character`, `lang` e `lang_group`.
A estrutura do canal ainda é `[meta, file]`.

### 2.4. Usar metadados para nomear e organizar as saídas

Com `lang` e `lang_group` agora disponíveis no meta map, podemos usá-los para adicionar um código de idioma aos nomes dos arquivos de saída e organizá-los em subdiretórios por família linguística.

Isso requer três alterações: atualizar o processo `COWPY` para renomear sua saída e incluir `meta` no que ele emite, atualizar a chamada ao `COWPY` para executar em `ch_languages`, e atualizar o bloco de saída para especificar o caminho do subdiretório.

#### 2.4.1. Atualizar o processo `COWPY`

Renomeie o arquivo de saída usando o código de idioma do meta map, e adicione `meta` à saída para que o bloco de saída possa acessar `lang_group` para roteamento de subdiretórios:

=== "Depois"

    ```groovy title="modules/cowpy.nf" linenums="9" hl_lines="2 6"
        output:
        tuple val(meta), path("${meta.lang}-${input_file}")

        script:
        """
        cat ${input_file} | cowpy -c ${meta.character} > ${meta.lang}-${input_file}
        """
    ```

=== "Antes"

    ```groovy title="modules/cowpy.nf" linenums="9" hl_lines="2 6"
        output:
        path "cowpy-${input_file}"

        script:
        """
        cat ${input_file} | cowpy -c ${meta.character} > cowpy-${input_file}
        """
    ```

Isso mostra como podemos aproveitar outros campos de metadados para personalizar o comportamento de um processo, sem precisar modificar a definição de entrada.

#### 2.4.2. Atualizar a chamada ao `COWPY` para executar em `ch_languages`

Substitua `COWPY(ch_datasheet)` por `COWPY(ch_languages)`:

=== "Depois"

    ```groovy title="main.nf" linenums="32" hl_lines="3"
        .set { ch_languages }

        COWPY(ch_languages)

        publish:
        cowpy_art = COWPY.out
    }
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="32" hl_lines="3 5"
        .set { ch_languages }

        ch_languages.view()

        COWPY(ch_datasheet)

        publish:
        cowpy_art = COWPY.out
    }
    ```

Também removemos a linha `ch_languages.view()` já que não precisamos mais inspecionar o conteúdo do canal.

#### 2.4.3. Atualizar o bloco de saída

Adicione uma closure `path` ao bloco `output {}` para rotear cada arquivo para o subdiretório do seu grupo linguístico:

=== "Depois"

    ```groovy title="main.nf" linenums="40" hl_lines="3"
    output {
        cowpy_art {
            path { meta, file -> meta.lang_group }
        }
    }
    ```

=== "Antes"

    ```groovy title="main.nf" linenums="40" hl_lines="2 3"
    output {
        cowpy_art {
        }
    }
    ```

Isso mostra como podemos usar metadados para organizar as saídas com grande flexibilidade.

#### 2.4.4. Executar o pipeline completo

Apague os resultados anteriores e execute o pipeline completo:

```bash
rm -r results
nextflow run main.nf
```

??? success "Saída do comando"

    ```console
     N E X T F L O W   ~  version 25.10.4

    Launching `main.nf` [suspicious_crick] DSL2 - revision: 25541014c5

    executor >  local (14)
    [5d/dffd4e] IDENTIFY_LANGUAGE (7) [100%] 7 of 7 ✔
    [e7/317c18] COWPY (7)             [100%] 7 of 7 ✔
    ```

O diretório de resultados agora está organizado por família linguística, com cada arquivo nomeado de acordo com o idioma detectado:

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

A closure `path` no bloco `output {}` recebe cada tupla `[meta, file]` e retorna `meta.lang_group` como nome do subdiretório.
O nome do arquivo em si vem do que o processo produz (`#!groovy "${meta.lang}-${input_file}"`).
Ambas as informações de metadados (código de idioma e grupo linguístico) vêm do meta map enriquecido construído nesta seção.

### Conclusão

Nesta seção, você aprendeu:

- **Como enriquecer o meta map com saídas de processos:** Adicionar novas chaves com `#!groovy meta + [chave: valor]` mantém a estrutura do canal `[meta, file]` intacta enquanto enriquece os metadados.
- **Como derivar metadados a partir de metadados:** A lógica condicional dentro de uma operação `map` pode calcular novos campos a partir dos existentes.
- **Como usar metadados para organizar as saídas:** A closure `path` no bloco `output {}` pode ler do meta map para rotear arquivos para subdiretórios.

---

## 3. Considerações sobre robustez

Quando valores de metadados orientam o comportamento dos processos, dados ausentes ou incompletos podem causar problemas difíceis de diagnosticar.
Veja o que esperar e como lidar com isso.

### 3.1. O que acontece quando um campo de metadados obrigatório está ausente

O valor `character` é obrigatório para que o processo `COWPY` produza um resultado válido.
O modo de falha depende se a coluna existe na planilha mas está vazia, ou se está completamente ausente.

#### 3.1.1. A coluna existe mas um valor está vazio

Suponha que uma entrada na planilha tenha o campo `character` em branco:

```csv title="datasheet.csv" linenums="1" hl_lines="2"
id,character,recording
sampleA,,/workspaces/training/side-quests/metadata/data/bonjour.txt
sampleB,tux,/workspaces/training/side-quests/metadata/data/guten_tag.txt
...
```

A chave `character` é criada para todas as entradas quando a planilha é analisada, mas `meta.character` para `sampleA` será uma string vazia.
Quando o Nextflow substitui `#!groovy ${meta.character}` no comando, a ferramenta `COWPY` recebe um argumento vazio para `-c` e falha:

??? failure "Saída do comando"

    ```console hl_lines="8 11 16 28"
     N E X T F L O W   ~  version 25.10.4

    Launching `main.nf` [marvelous_hirsch] DSL2 - revision: 0dfeee3cc1

    executor >  local (9)
    [c1/c5dd4f] process > IDENTIFY_LANGUAGE (7) [ 85%] 6 of 7
    [d3/b7c415] process > COWPY (2)             [  0%] 0 of 6
    ERROR ~ Error executing process > 'COWPY (1)'

    Caused by:
      Process `COWPY (1)` terminated with an error exit status (2)


    Command executed:

      cat bonjour.txt | cowpy -c  > fr-bonjour.txt

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

A mensagem de erro (`expected one argument`) aponta para o flag `-c` vazio.
Verificar o arquivo `.command.sh` no diretório de trabalho confirma que o comando foi executado com um valor vazio.

#### 3.1.2. A coluna não existe na planilha

Se a coluna `character` estiver completamente ausente:

```csv title="datasheet.csv" linenums="1"
id,recording
sampleA,/workspaces/training/side-quests/metadata/data/bonjour.txt
sampleB,/workspaces/training/side-quests/metadata/data/guten_tag.txt
...
```

A chave `character` nunca é criada no meta map.
Quando o script do processo avalia `#!groovy ${meta.character}`, a chave ausente retorna `null`, e o Nextflow literalmente substitui a string `null` no comando:

??? failure "Saída do comando"

    ```console hl_lines="8 11 16"
     N E X T F L O W   ~  version 25.10.4

    Launching `main.nf` [jovial_bohr] DSL2 - revision: eaaf375827

    executor >  local (9)
    [0d/ada9db] process > IDENTIFY_LANGUAGE (5) [ 85%] 6 of 7
    [06/28065f] process > COWPY (2)             [  0%] 0 of 6
    ERROR ~ Error executing process > 'COWPY (2)'

    Caused by:
      Process `COWPY (2)` terminated with an error exit status (1)


    Command executed:

      cat guten_tag.txt | cowpy -c null > de-guten_tag.txt

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

O `cowpy -c null` no comando executado é a pista diagnóstica.

### 3.2. Estratégias para lidar com metadados ausentes

Existem duas abordagens complementares para tornar os fluxos de trabalho mais robustos contra metadados ausentes.

**1. Validação de entrada**

A solução mais confiável é validar a planilha antes de qualquer processamento começar, para que os problemas sejam detectados cedo com uma mensagem de erro clara, em vez de aparecerem como uma falha críptica de processo no meio da execução.
O treinamento [Hello nf-core](../../hello_nf-core/05_input_validation.md) aborda como adicionar validação de entrada usando o plugin nf-schema. <!-- TODO (future) pending a proper Validation side quest -->

**2. Entradas explícitas do processo para valores obrigatórios**

Se você quiser que a interface do processo comunique que um determinado valor é obrigatório, considere extraí-lo do meta map como uma entrada explícita:

=== "Definição do processo"

    ```groovy title="modules/cowpy.nf" linenums="6"
    input:
    tuple val(meta), val(character), path(input_file)
    ```

=== "Chamada no fluxo de trabalho"

    ```groovy title="main.nf"
    COWPY(ch_languages.map { meta, file -> [meta, meta.character, file] })
    ```

Essa abordagem torna `character` uma parte visível e obrigatória do contrato do processo.
Qualquer pessoa que leia o módulo pode ver imediatamente que um valor de personagem deve ser fornecido.
Se o campo estiver ausente, o fluxo de trabalho falha claramente no nível do canal antes mesmo de o processo ser executado.

Isso destaca um princípio de design útil:

**Use o meta map para informações opcionais ou descritivas; extraia valores obrigatórios como entradas explícitas.**

O meta map mantém as estruturas de canal limpas e estáveis, mas para valores genuinamente obrigatórios por um processo, torná-los entradas nomeadas melhora a clareza e facilita o uso correto do módulo em outros contextos.

### Conclusão

Nesta seção, você viu:

- **Como metadados ausentes se manifestam:** Um campo vazio produz um argumento vazio; um campo ausente produz `null` substituído literalmente no comando.
- **Duas estratégias complementares:** Validação de entrada para detectar problemas cedo, e entradas explícitas do processo para comunicar os requisitos claramente.

---

## Resumo

Nesta missão paralela, você explorou como trabalhar efetivamente com metadados em fluxos de trabalho Nextflow.

O padrão de tupla "meta map + arquivo de dados" é uma convenção central no Nextflow, oferecendo várias vantagens em relação a passar metadados como valores individuais:

- A estrutura do canal permanece estável à medida que a planilha evolui
- O comportamento dos processos pode ser personalizado por amostra sem codificar nomes de campos diretamente
- Os metadados ficam disponíveis ao longo do pipeline para nomear, agrupar e organizar as saídas
- Módulos escritos para essa interface são intercambiáveis, incluindo módulos nf-core

### Padrões principais

1.  **Leitura e estruturação de metadados:** Analisar uma planilha CSV e criar um meta map.

    ```groovy
    channel.fromPath('datasheet.csv')
        .splitCsv(header: true)
        .map { row ->
            [ [id: row.id, character: row.character], row.recording ]
        }
    ```

2.  **Expansão de metadados durante o fluxo de trabalho:** Adicionar novas chaves a partir de saídas de processos ou lógica derivada.

    ```groovy
    // A partir de uma saída de processo
    .map { meta, file, lang ->
        [ meta + [lang: lang], file ]
    }

    // A partir de lógica condicional
    .map { meta, file ->
        def lang_group = "unknown"
        if (meta.lang in ["de", "en"]) { lang_group = "germanic" }
        else if (meta.lang in ["fr", "es", "it"]) { lang_group = "romance" }
        [ meta + [lang_group: lang_group], file ]
    }
    ```

3.  **Usando metadados dentro de um processo:** Acesse qualquer campo via notação de ponto no bloco de script.

    ```groovy
    cat ${input_file} | cowpy -c ${meta.character} > ${meta.lang}-${input_file}
    ```

4.  **Organizando saídas por valor de metadados:** Use a closure `path` no bloco `output {}`.

    ```groovy
    output {
        cowpy_art {
            path { meta, file -> meta.lang_group }
        }
    }
    ```

### Recursos adicionais

- [operador map](https://www.nextflow.io/docs/latest/operator.html#map)
- [operador multiMap](https://www.nextflow.io/docs/latest/reference/operator.html#multimap)
- [qualificador de saída stdout](https://www.nextflow.io/docs/latest/process.html#outputs)

---

## O que vem a seguir?

Volte ao [menu de Missões Paralelas](../index.md) ou clique no botão no canto inferior direito da página para avançar para o próximo tópico da lista.
