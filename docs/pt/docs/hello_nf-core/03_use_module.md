# Parte 3: Usar um módulo nf-core

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tradução assistida por IA - [saiba mais e sugira melhorias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Nesta terceira parte do curso de treinamento Hello nf-core, mostramos como encontrar, instalar e usar um módulo nf-core existente no seu pipeline.

Um dos grandes benefícios de trabalhar com nf-core é a capacidade de aproveitar módulos pré-construídos e testados do repositório [nf-core/modules](https://github.com/nf-core/modules).
Em vez de escrever cada processo do zero, você pode instalar e usar módulos mantidos pela comunidade que seguem as melhores práticas.

Para demonstrar como isso funciona, vamos substituir o módulo customizado `collectGreetings` pelo módulo `cat/cat` do nf-core/modules no pipeline `core-hello`.

??? info "Como começar a partir desta seção"

    Esta seção do curso assume que você completou a [Parte 2: Reescrever Hello para nf-core](./02_rewrite_hello.md) e tem um pipeline `core-hello` funcionando.

    Se você não completou a Parte 2 ou quer começar do zero para esta parte, você pode usar a solução `core-hello-part2` como ponto de partida.
    Execute este comando dentro do diretório `hello-nf-core/`:

    ```bash
    cp -r solutions/core-hello-part2 core-hello
    cd core-hello
    ```

    Isso fornece um pipeline nf-core totalmente funcional pronto para adicionar módulos.
    Você pode testar que ele executa com sucesso rodando o seguinte comando:

    ```bash
    nextflow run . --outdir core-hello-results -profile test,docker --validate_params false
    ```

---

## 1. Encontrar e instalar um módulo nf-core adequado

Primeiro, vamos aprender como encontrar um módulo nf-core existente e instalá-lo no nosso pipeline.

Nosso objetivo é substituir o processo `collectGreetings`, que usa o comando Unix `cat` para concatenar múltiplos arquivos de saudações em um só.
Concatenar arquivos é uma operação muito comum, então é razoável que já exista um módulo no nf-core projetado para esse propósito.

Vamos começar.

### 1.1. Navegar pelos módulos disponíveis no site nf-core

O projeto nf-core mantém um catálogo centralizado de módulos em [https://nf-co.re/modules](https://nf-co.re/modules).

Navegue até a página de módulos no seu navegador e use a barra de busca para pesquisar 'concatenate'.

![resultados da busca de módulos](./img/module-search-results.png)

Como você pode ver, há muitos resultados, vários deles módulos projetados para concatenar tipos muito específicos de arquivos.
Entre eles, você deve ver um chamado `cat_cat` que é de uso geral.

!!! note "Convenção de nomenclatura de módulos"

    O sublinhado (`_`) é usado como substituto para o caractere barra (`/`) nos nomes dos módulos.

    Os módulos nf-core seguem a convenção de nomenclatura `software/comando` quando uma ferramenta fornece múltiplos comandos, como `samtools/view` (pacote samtools, comando view) ou `gatk/haplotypecaller` (pacote GATK, comando HaplotypeCaller).
    Para ferramentas que fornecem apenas um comando principal, os módulos usam um único nível como `fastqc` ou `multiqc`.

Clique na caixa do módulo `cat_cat` para visualizar a documentação do módulo.

A página do módulo mostra:

- Uma breve descrição: "A module for concatenation of gzipped or uncompressed files"
- Comando de instalação: `nf-core modules install cat/cat`
- Estrutura dos canais de entrada e saída
- Parâmetros disponíveis

### 1.2. Listar módulos disponíveis pela linha de comando

Alternativamente, você também pode pesquisar módulos diretamente da linha de comando usando as ferramentas nf-core.

```bash
nf-core modules list remote
```

Isso exibirá uma lista de todos os módulos disponíveis no repositório nf-core/modules, embora seja um pouco menos conveniente se você ainda não sabe o nome do módulo que está procurando.
No entanto, se você souber, pode canalizar a lista para `grep` para encontrar módulos específicos:

```bash
nf-core modules list remote | grep 'cat/cat'
```

??? success "Saída do comando"

    ```console
    │ cat/cat
    ```

Tenha em mente que a abordagem com `grep` só extrairá resultados com o termo de busca em seu nome, o que não funcionaria para `cat_cat`.

### 1.3. Obter informações detalhadas sobre o módulo

Para ver informações detalhadas sobre um módulo específico da linha de comando, use o comando `info`:

```bash
nf-core modules info cat/cat
```

Isso exibe a documentação sobre o módulo, incluindo suas entradas, saídas e informações básicas de uso.

??? success "Saída do comando"

    ```console

                                              ,--./,-.
              ___     __   __   __   ___     /,-._.--~\
        |\ | |__  __ /  ` /  \ |__) |__         }  {
        | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                              `._,._,'

        nf-core/tools version 3.4.1 - https://nf-co.re


    ╭─ Module: cat/cat  ─────────────────────────────────────────────────╮
    │ 🌐 Repository: https://github.com/nf-core/modules.git              │
    │ 🔧 Tools: cat                                                      │
    │ 📖 Description: A module for concatenation of gzipped or           │
    │ uncompressed files                                                 │
    ╰────────────────────────────────────────────────────────────────────╯
                      ╷                                          ╷
    📥 Inputs        │Description                               │Pattern
    ╺━━━━━━━━━━━━━━━━━┿━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┿━━━━━━━╸
    input[0]         │                                          │
    ╶─────────────────┼──────────────────────────────────────────┼───────╴
      meta  (map)     │Groovy Map containing sample information  │
                      │e.g. [ id:'test', single_end:false ]      │
    ╶─────────────────┼──────────────────────────────────────────┼───────╴
      files_in  (file)│List of compressed / uncompressed files   │      *
                      ╵                                          ╵
                          ╷                                 ╷
    📥 Outputs           │Description                      │     Pattern
    ╺━━━━━━━━━━━━━━━━━━━━━┿━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┿━━━━━━━━━━━━╸
    file_out             │                                 │
    ╶─────────────────────┼─────────────────────────────────┼────────────╴
      meta  (map)         │Groovy Map containing sample     │
                          │information                      │
    ╶─────────────────────┼─────────────────────────────────┼────────────╴
      ${prefix}  (file)   │Concatenated file. Will be       │ ${file_out}
                          │gzipped if file_out ends with    │
                          │".gz"                            │
    ╶─────────────────────┼─────────────────────────────────┼────────────╴
    versions             │                                 │
    ╶─────────────────────┼─────────────────────────────────┼────────────╴
      versions.yml  (file)│File containing software versions│versions.yml
                          ╵                                 ╵

    💻  Installation command: nf-core modules install cat/cat

    ```

Esta é exatamente a mesma informação que você pode encontrar no site.

### 1.4. Instalar o módulo cat/cat

Agora que encontramos o módulo que queremos, precisamos adicioná-lo ao código-fonte do nosso pipeline.

A boa notícia é que o projeto nf-core inclui ferramentas para facilitar esta parte.
Especificamente, o comando `nf-core modules install` torna possível automatizar a recuperação do código e torná-lo disponível para seu projeto em um único passo.

Navegue até o diretório do seu pipeline e execute o comando de instalação:

```bash
cd core-hello
nf-core modules install cat/cat
```

A ferramenta pode primeiro solicitar que você especifique um tipo de repositório.
(Se não, pule para "Finalmente, a ferramenta procederá para instalar o módulo.")

??? success "Saída do comando"

    ```console

                                          ,--./,-.
          ___     __   __   __   ___     /,-._.--~\
    |\ | |__  __ /  ` /  \ |__) |__         }  {
    | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                          `._,._,'

    nf-core/tools version 3.4.1 - https://nf-co.re


    WARNING  'repository_type' not defined in .nf-core.yml
    ? Is this repository a pipeline or a modules repository? (Use arrow keys)
    » Pipeline
      Modules repository
    ```

Se aparecer, pressione enter para aceitar a resposta padrão (`Pipeline`) e continuar.

A ferramenta então oferecerá alterar a configuração do seu projeto para evitar esse prompt no futuro.

??? success "Saída do comando"

    ```console
        INFO     To avoid this prompt in the future, add the 'repository_type' key to your .nf-core.yml file.
        ? Would you like me to add this config now? [y/n] (y):
    ```

Por que não aproveitar essa ferramenta conveniente!
Pressione enter para aceitar a resposta padrão (sim).

Finalmente, a ferramenta procederá para instalar o módulo.

??? success "Saída do comando"

    ```console
    INFO Config added to '.nf-core.yml'
    INFO Reinstalling modules found in 'modules.json' but missing from directory:
    INFO Installing 'cat/cat'
    INFO Use the following statement to include this module:

        include { CAT_CAT } from '../modules/nf-core/cat/cat/main'
    ```

O comando automaticamente:

- Baixa os arquivos do módulo para `modules/nf-core/cat/cat/`
- Atualiza `modules.json` para rastrear o módulo instalado
- Fornece a declaração `include` correta para usar no seu fluxo de trabalho

!!! tip

    Sempre certifique-se de que seu diretório de trabalho atual seja a raiz do projeto do seu pipeline antes de executar o comando de instalação de módulo.

Vamos verificar que o módulo foi instalado corretamente:

```bash
tree -L 4 modules
```

??? abstract "Conteúdo do diretório"

    ```console
    modules
    ├── local
    │   ├── collectGreetings.nf
    │   ├── convertToUpper.nf
    │   ├── cowpy.nf
    │   └── sayHello.nf
    └── nf-core
        └── cat
            └── cat
                ├── environment.yml
                ├── main.nf
                ├── meta.yml
                └── tests

    5 directories, 7 files
    ```

Você também pode verificar a instalação pedindo ao utilitário nf-core para listar os módulos instalados localmente:

```bash
nf-core modules list local
```

??? success "Saída do comando"

    ```console
    INFO     Repository type: pipeline
    INFO     Modules installed in '.':

    ┏━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━┓
    ┃ Module Name ┃ Repository      ┃ Version SHA ┃ Message                                ┃ Date       ┃
    ┡━━━━━━━━━━━━━╇━━━━━━━━━━━━━━━━━╇━━━━━━━━━━━━━╇━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╇━━━━━━━━━━━━┩
    │ cat/cat     │ nf-core/modules │ 41dfa3f     │ update meta.yml of all modules (#8747) │ 2025-07-07 │
    └─────────────┴─────────────────┴─────────────┴────────────────────────────────────────┴────────────┘
    ```

Isso confirma que o módulo `cat/cat` agora faz parte do código-fonte do seu projeto.

No entanto, para realmente usar o novo módulo, precisamos importá-lo no nosso pipeline.

### 1.5. Atualizar as importações de módulos

Vamos substituir a declaração `include` do módulo `collectGreetings` pela do `CAT_CAT` na seção de importações do fluxo de trabalho `workflows/hello.nf`.

Como lembrete, a ferramenta de instalação de módulos nos deu a declaração exata para usar:

```groovy title="Declaração de importação produzida pelo comando de instalação"
include { CAT_CAT } from '../modules/nf-core/cat/cat/main'`
```

Observe que a convenção nf-core é usar letras maiúsculas para nomes de módulos ao importá-los.

Abra [core-hello/workflows/hello.nf](core-hello/workflows/hello.nf) e faça a seguinte substituição:

=== "Depois"

    ```groovy title="core-hello/workflows/hello.nf" linenums="1" hl_lines="10"
    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    include { paramsSummaryMap       } from 'plugin/nf-schema'
    include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
    include { sayHello               } from '../modules/local/sayHello.nf'
    include { convertToUpper         } from '../modules/local/convertToUpper.nf'
    include { CAT_CAT                } from '../modules/nf-core/cat/cat/main'
    include { cowpy                  } from '../modules/local/cowpy.nf'
    ```

=== "Antes"

    ```groovy title="core-hello/workflows/hello.nf" linenums="1" hl_lines="10"
    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    include { paramsSummaryMap       } from 'plugin/nf-schema'
    include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
    include { sayHello               } from '../modules/local/sayHello.nf'
    include { convertToUpper         } from '../modules/local/convertToUpper.nf'
    include { collectGreetings       } from '../modules/local/collectGreetings.nf'
    include { cowpy                  } from '../modules/local/cowpy.nf'
    ```

Observe como o caminho para o módulo nf-core difere dos módulos locais:

- **Módulo nf-core**: `'../modules/nf-core/cat/cat/main'` (referencia `main.nf`)
- **Módulo local**: `'../modules/local/collectGreetings.nf'` (referência de arquivo único)

O módulo agora está disponível para o fluxo de trabalho, então tudo o que precisamos fazer é trocar a chamada de `collectGreetings` para usar `CAT_CAT`. Certo?

Não tão rápido.

Neste ponto, você pode estar tentado a mergulhar e começar a editar o código, mas vale a pena dedicar um momento para examinar cuidadosamente o que o novo módulo espera e o que ele produz.

Vamos abordar isso como uma seção separada porque envolve um novo mecanismo que ainda não cobrimos: mapas de metadados.

!!! note

    Você pode opcionalmente deletar o arquivo `collectGreetings.nf`:

    ```bash
    rm modules/local/collectGreetings.nf
    ```

    No entanto, você pode querer mantê-lo como referência para entender as diferenças entre módulos locais e nf-core.

### Conclusão

Você sabe como encontrar um módulo nf-core e torná-lo disponível para o seu projeto.

### Qual é o próximo passo?

Avaliar o que um novo módulo requer e identificar quaisquer mudanças importantes necessárias para integrá-lo em um pipeline.

---

## 2. Avaliar os requisitos do novo módulo

Especificamente, precisamos examinar a **interface** do módulo, ou seja, suas definições de entrada e saída, e compará-la com a interface do módulo que estamos buscando substituir.
Isso nos permitirá determinar se podemos simplesmente tratar o novo módulo como uma substituição direta ou se precisaremos adaptar parte da conexão.

Idealmente, isso é algo que você deveria fazer _antes_ mesmo de instalar o módulo, mas ei, melhor tarde do que nunca.
(Vale ressaltar que existe um comando `uninstall` para se livrar de módulos que você decidir que não quer mais.)

!!! note

    O processo CAT_CAT inclui um tratamento bastante inteligente de diferentes tipos de compressão, extensões de arquivo e assim por diante, que não são estritamente relevantes para o que estamos tentando mostrar aqui, então vamos ignorar a maior parte e focar apenas nas partes que são importantes.

### 2.1. Comparar as interfaces dos dois módulos

Como lembrete, esta é a aparência da interface do nosso módulo `collectGreetings`:

```groovy title="modules/local/collectGreetings.nf (trecho)" linenums="1" hl_lines="6-7 10"
process collectGreetings {

    publishDir 'results', mode: 'copy'

    input:
        path input_files
        val batch_name

    output:
        path "COLLECTED-${batch_name}-output.txt" , emit: outfile
```

O módulo `collectGreetings` recebe duas entradas:

- `input_files` contém um ou mais arquivos de entrada para processar;
- `batch_name` é um valor que usamos para atribuir um nome específico da execução ao arquivo de saída, que é uma forma de metadados.

Após a conclusão, `collectGreetings` produz um único caminho de arquivo, emitido com a tag `outfile`.

Em comparação, a interface do módulo `cat/cat` é mais complexa:

```groovy title="modules/nf-core/cat/cat/main.nf (trecho)" linenums="1" hl_lines="11 14"
process CAT_CAT {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pigz:2.3.4' :
        'biocontainers/pigz:2.3.4' }"

    input:
    tuple val(meta), path(files_in)

    output:
    tuple val(meta), path("${prefix}"), emit: file_out
    path "versions.yml"               , emit: versions
```

O módulo CAT_CAT recebe uma única entrada, mas essa entrada é uma tupla contendo duas coisas:

- `meta` é uma estrutura contendo metadados, chamada de mapa de metadados (metamap);
- `files_in` contém um ou mais arquivos de entrada para processar, equivalente ao `input_files` de `collectGreetings`.

Após a conclusão, CAT_CAT entrega suas saídas em duas partes:

- Outra tupla contendo o mapa de metadados e o arquivo de saída concatenado, emitido com a tag `file_out`;
- Um arquivo `versions.yml` que captura informações sobre a versão do software que foi usada, emitido com a tag `versions`.

Observe também que, por padrão, o arquivo de saída será nomeado com base em um identificador que faz parte dos metadados (código não mostrado aqui).

Isso pode parecer muita coisa para acompanhar apenas olhando o código, então aqui está um diagrama para ajudá-lo a visualizar como tudo se encaixa.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nf-core/img/module_comparison.svg"
</figure>

Você pode ver que os dois módulos têm requisitos de entrada semelhantes em termos de conteúdo (um conjunto de arquivos de entrada mais alguns metadados), mas expectativas muito diferentes sobre como esse conteúdo é empacotado.
Ignorando o arquivo de versões por enquanto, sua saída principal também é equivalente (um arquivo concatenado), exceto que CAT_CAT também emite o mapa de metadados em conjunto com o arquivo de saída.

As diferenças de empacotamento serão relativamente fáceis de lidar, como você verá daqui a pouco.
No entanto, para entender a parte do mapa de metadados, precisamos apresentar algum contexto adicional.

### 2.2. Entendendo mapas de metadados

Acabamos de dizer que o módulo CAT_CAT espera um mapa de metadados como parte de sua tupla de entrada.
Vamos dedicar alguns minutos para examinar mais de perto o que isso é.

O **mapa de metadados**, frequentemente referido como **metamap** para abreviar, é um mapa no estilo Groovy contendo informações sobre unidades de dados.
No contexto de pipelines Nextflow, unidades de dados podem ser qualquer coisa que você queira: amostras individuais, lotes de amostras ou conjuntos de dados inteiros.

Por convenção, um mapa de metadados nf-core é nomeado `meta` e contém o campo obrigatório `id`, que é usado para nomear saídas e rastrear unidades de dados.

Por exemplo, um mapa de metadados típico pode se parecer com isto:

```groovy title="Exemplo de metamap em nível de amostra"
[id: 'sample1', single_end: false, strandedness: 'forward']
```

Ou em um caso onde os metadados são anexados no nível de lote:

```groovy title="Exemplo de metamap em nível de lote"
[id: 'batch1', date: '25.10.01']
```

Agora vamos colocar isso no contexto do processo `CAT_CAT`, que espera que os arquivos de entrada sejam empacotados em uma tupla com um mapa de metadados, e também produz o mapa de metadados como parte da tupla de saída.

```groovy title="modules/nf-core/cat/cat/main.nf (trecho)" linenums="1" hl_lines="2 5"
input:
tuple val(meta), path(files_in)

output:
tuple val(meta), path("${prefix}"), emit: file_out
```

Como resultado, cada unidade de dados viaja através do pipeline com os metadados relevantes anexados.
Processos subsequentes podem então acessar prontamente esses metadados também.

Lembra como dissemos que o arquivo produzido por `CAT_CAT` será nomeado com base em um identificador que faz parte dos metadados?
Este é o código relevante:

```groovy title="modules/nf-core/cat/cat/main.nf (trecho)" linenums="35"
prefix   = task.ext.prefix ?: "${meta.id}${getFileSuffix(file_list[0])}"
```

Isso se traduz aproximadamente da seguinte forma: se um `prefix` for fornecido através do sistema de parâmetros externos de tarefa (`task.ext`), use-o para nomear o arquivo de saída; caso contrário, crie um usando `${meta.id}`, que corresponde ao campo `id` no mapa de metadados.

Você pode imaginar o canal de entrada chegando a este módulo com conteúdo assim:

```groovy title="Exemplo de conteúdo do canal de entrada"
ch_input = [[[id: 'batch1', date: '25.10.01'], ['file1A.txt', 'file1B.txt']],
            [[id: 'batch2', date: '25.10.26'], ['file2A.txt', 'file2B.txt']],
            [[id: 'batch3', date: '25.11.14'], ['file3A.txt', 'file3B.txt']]]
```

Então o conteúdo do canal de saída saindo assim:

```groovy title="Exemplo de conteúdo do canal de saída"
ch_input = [[[id: 'batch1', date: '25.10.01'], 'batch1.txt'],
            [[id: 'batch2', date: '25.10.26'], 'batch2.txt'],
            [[id: 'batch3', date: '25.11.14'], 'batch3.txt']]
```

Como mencionado anteriormente, a configuração de entrada `tuple val(meta), path(files_in)` é um padrão usado em todos os módulos nf-core.

Esperamos que você possa começar a ver como isso pode ser útil.
Não só permite nomear saídas com base em metadados, mas você também pode fazer coisas como usá-lo para aplicar diferentes valores de parâmetros e, em combinação com operadores específicos, você pode até agrupar, ordenar ou filtrar dados conforme eles fluem pelo pipeline.

!!! note "Saiba mais sobre metadados"

    Para uma introdução abrangente ao trabalho com metadados em fluxos de trabalho Nextflow, incluindo como ler metadados de planilhas de amostras e usá-los para personalizar o processamento, consulte a missão paralela [Metadados em fluxos de trabalho](../side_quests/metadata).

### 2.3. Resumir mudanças a serem feitas

Com base no que revisamos, estas são as principais mudanças que precisamos fazer no nosso pipeline para utilizar o módulo `cat/cat`:

- Criar um mapa de metadados contendo o nome do lote;
- Empacotar o mapa de metadados em uma tupla com o conjunto de arquivos de entrada para concatenar (vindo de `convertToUpper`);
- Trocar a chamada de `collectGreetings()` para `CAT_CAT`;
- Extrair o arquivo de saída da tupla produzida pelo processo `CAT_CAT` antes de passá-lo para `cowpy`.

Isso deve resolver! Agora que temos um plano, estamos prontos para mergulhar.

### Conclusão

Você sabe como avaliar a interface de entrada e saída de um novo módulo para identificar seus requisitos, e aprendeu como os mapas de metadados são usados por pipelines nf-core para manter metadados intimamente associados aos dados conforme eles fluem por um pipeline.

### Qual é o próximo passo?

Integrar o novo módulo em um fluxo de trabalho.

---

## 3. Integrar CAT_CAT no fluxo de trabalho `hello.nf`

Agora que você sabe tudo sobre mapas de metadados (ou o suficiente para os propósitos deste curso, de qualquer forma), é hora de realmente implementar as mudanças que delineamos acima.

Para maior clareza, vamos dividir isso e cobrir cada passo separadamente.

!!! note

    Todas as mudanças mostradas abaixo são feitas na lógica do fluxo de trabalho no bloco `main` no arquivo de fluxo de trabalho `core-hello/workflows/hello.nf`.

### 3.1. Criar um mapa de metadados

Primeiro, precisamos criar um mapa de metadados para `CAT_CAT`, lembrando que os módulos nf-core exigem que o mapa de metadados tenha pelo menos um campo `id`.

Como não precisamos de outros metadados, podemos manter simples e usar algo assim:

```groovy title="Exemplo de sintaxe"
def cat_meta = [id: 'test']
```

Exceto que não queremos codificar o valor `id`; queremos usar o valor do parâmetro `params.batch`.
Então o código fica:

```groovy title="Exemplo de sintaxe"
def cat_meta = [id: params.batch]
```

Sim, é literalmente assim tão simples criar um mapa de metadados básico.

Vamos adicionar essas linhas após a chamada de `convertToUpper`, removendo a chamada de `collectGreetings`:

=== "Depois"

    ```groovy title="core-hello/workflows/hello.nf" linenums="26" hl_lines="7-8"
        // emitir uma saudação
        sayHello(ch_samplesheet)

        // converter a saudação para maiúsculas
        convertToUpper(sayHello.out)

        // criar mapa de metadados com o nome do lote como ID
        def cat_meta = [ id: params.batch ]

        // gerar arte ASCII das saudações com cowpy
        cowpy(collectGreetings.out.outfile, params.character)
    ```

=== "Antes"

    ```groovy title="core-hello/workflows/hello.nf" linenums="26" hl_lines="7-8"
        // emitir uma saudação
        sayHello(ch_samplesheet)

        // converter a saudação para maiúsculas
        convertToUpper(sayHello.out)

        // coletar todas as saudações em um arquivo
        collectGreetings(convertToUpper.out.collect(), params.batch)

        // gerar arte ASCII das saudações com cowpy
        cowpy(collectGreetings.out.outfile, params.character)
    ```

Isso cria um mapa de metadados simples onde o `id` é definido como nosso nome de lote (que será `test` ao usar o perfil de teste).

### 3.2. Criar um canal com tuplas de metadados

Em seguida, transforme o canal de arquivos em um canal de tuplas contendo metadados e arquivos:

=== "Depois"

    ```groovy title="core-hello/workflows/hello.nf" linenums="26" hl_lines="10-11"
        // emitir uma saudação
        sayHello(ch_samplesheet)

        // converter a saudação para maiúsculas
        convertToUpper(sayHello.out)

        // criar mapa de metadados com o nome do lote como ID
        def cat_meta = [ id: params.batch ]

        // criar um canal com metadados e arquivos no formato de tupla
        ch_for_cat = convertToUpper.out.collect().map { files -> tuple(cat_meta, files) }

        // gerar arte ASCII das saudações com cowpy
        cowpy(collectGreetings.out.outfile, params.character)
    ```

=== "Antes"

    ```groovy title="core-hello/workflows/hello.nf" linenums="26"
        // emitir uma saudação
        sayHello(ch_samplesheet)

        // converter a saudação para maiúsculas
        convertToUpper(sayHello.out)

        // criar mapa de metadados com o nome do lote como ID
        def cat_meta = [ id: params.batch ]

        // gerar arte ASCII das saudações com cowpy
        cowpy(collectGreetings.out.outfile, params.character)
    ```

A linha que adicionamos realiza duas coisas:

- `.collect()` reúne todos os arquivos da saída de `convertToUpper` em uma única lista
- `.map { files -> tuple(cat_meta, files) }` cria uma tupla de `[metadados, arquivos]` no formato que `CAT_CAT` espera

Isso é tudo que precisamos fazer para configurar a tupla de entrada para `CAT_CAT`.

### 3.3. Chamar o módulo CAT_CAT

Agora chame `CAT_CAT` no canal recém-criado:

=== "Depois"

    ```groovy title="core-hello/workflows/hello.nf" linenums="26" hl_lines="13-14"
        // emitir uma saudação
        sayHello(ch_samplesheet)

        // converter a saudação para maiúsculas
        convertToUpper(sayHello.out)

        // criar mapa de metadados com o nome do lote como ID
        def cat_meta = [ id: params.batch ]

        // criar um canal com metadados e arquivos no formato de tupla
        ch_for_cat = convertToUpper.out.collect().map { files -> tuple(cat_meta, files) }

        // concatenar arquivos usando o módulo nf-core cat/cat
        CAT_CAT(ch_for_cat)

        // gerar arte ASCII das saudações com cowpy
        cowpy(collectGreetings.out.outfile, params.character)
    ```

=== "Antes"

    ```groovy title="core-hello/workflows/hello.nf" linenums="26"
        // emitir uma saudação
        sayHello(ch_samplesheet)

        // converter a saudação para maiúsculas
        convertToUpper(sayHello.out)

        // criar mapa de metadados com o nome do lote como ID
        def cat_meta = [ id: params.batch ]

        // criar um canal com metadados e arquivos no formato de tupla
        ch_for_cat = convertToUpper.out.collect().map { files -> tuple(cat_meta, files) }

        // gerar arte ASCII das saudações com cowpy
        cowpy(collectGreetings.out.outfile, params.character)
    ```

Isso completa a parte mais complicada desta substituição, mas ainda não terminamos: ainda precisamos atualizar como passamos a saída concatenada para o processo `cowpy`.

### 3.4. Extrair o arquivo de saída da tupla para `cowpy`

Anteriormente, o processo `collectGreetings` simplesmente produzia um arquivo que podíamos passar para `cowpy` diretamente.
No entanto, o processo `CAT_CAT` produz uma tupla que inclui o mapa de metadados além do arquivo de saída.

Como `cowpy` ainda não aceita tuplas de metadados (vamos corrigir isso na próxima parte do curso), precisamos extrair o arquivo de saída da tupla produzida por `CAT_CAT` antes de entregá-lo a `cowpy`:

=== "Depois"

    ```groovy title="core-hello/workflows/hello.nf" linenums="26" hl_lines="16-17 20"
        // emitir uma saudação
        sayHello(ch_samplesheet)

        // converter a saudação para maiúsculas
        convertToUpper(sayHello.out)

        // criar mapa de metadados com o nome do lote como ID
        def cat_meta = [ id: params.batch ]

        // criar um canal com metadados e arquivos no formato de tupla
        ch_for_cat = convertToUpper.out.collect().map { files -> tuple(cat_meta, files) }

        // concatenar as saudações
        CAT_CAT(ch_for_cat)

        // extrair o arquivo da tupla já que cowpy ainda não usa metadados
        ch_for_cowpy = CAT_CAT.out.file_out.map{ meta, file -> file }

        // gerar arte ASCII das saudações com cowpy
        cowpy(ch_for_cowpy, params.character)
    ```

=== "Antes"

    ```groovy title="core-hello/workflows/hello.nf" linenums="26" hl_lines="17"
        // emitir uma saudação
        sayHello(ch_samplesheet)

        // converter a saudação para maiúsculas
        convertToUpper(sayHello.out)

        // criar mapa de metadados com o nome do lote como ID
        def cat_meta = [ id: params.batch ]

        // criar um canal com metadados e arquivos no formato de tupla
        ch_for_cat = convertToUpper.out.collect().map { files -> tuple(cat_meta, files) }

        // concatenar as saudações
        CAT_CAT(ch_for_cat)

        // gerar arte ASCII das saudações com cowpy
        cowpy(collectGreetings.out.outfile, params.character)
    ```

A operação `.map{ meta, file -> file }` extrai o arquivo da tupla `[metadados, arquivo]` produzida por `CAT_CAT` em um novo canal, `ch_for_cowpy`.

Então é só uma questão de passar `ch_for_cowpy` para `cowpy` em vez de `collectGreetings.out.outfile` naquela última linha.

!!! note

    Na próxima parte do curso, atualizaremos `cowpy` para trabalhar diretamente com tuplas de metadados, então este passo de extração não será mais necessário.

### 3.5. Testar o fluxo de trabalho

Vamos testar que o fluxo de trabalho funciona com o módulo `cat/cat` recém-integrado:

```bash
nextflow run . --outdir core-hello-results -profile test,docker --validate_params false
```

Isso deve executar razoavelmente rápido.

??? success "Saída do comando"

    ```console
    N E X T F L O W ~ version 25.04.3

        Launching `./main.nf` [evil_pike] DSL2 - revision: b9e9b3b8de

        Input/output options
          input                     : /workspaces/training/hello-nf-core/core-hello/assets/greetings.csv
          outdir                    : core-hello-results

        Institutional config options
          config_profile_name       : Test profile
          config_profile_description: Minimal test dataset to check pipeline function

        Generic options
          validate_params           : false
          trace_report_suffix       : 2025-10-30_18-50-58

        Core Nextflow options
          runName                   : evil_pike
          containerEngine           : docker
          launchDir                 : /workspaces/training/hello-nf-core/core-hello
          workDir                   : /workspaces/training/hello-nf-core/core-hello/work
          projectDir                : /workspaces/training/hello-nf-core/core-hello
          userName                  : root
          profile                   : test,docker
          configFiles               : /workspaces/training/hello-nf-core/core-hello/nextflow.config

        !! Only displaying parameters that differ from the pipeline defaults !!
        ------------------------------------------------------
        executor >  local (8)
        [b3/f005fd] CORE_HELLO:HELLO:sayHello (3)       [100%] 3 of 3 ✔
        [08/f923d0] CORE_HELLO:HELLO:convertToUpper (3) [100%] 3 of 3 ✔
        [34/3729a9] CORE_HELLO:HELLO:CAT_CAT (test)     [100%] 1 of 1 ✔
        [24/df918a] CORE_HELLO:HELLO:cowpy              [100%] 1 of 1 ✔
        -[core/hello] Pipeline completed successfully-
    ```

Observe que `CAT_CAT` agora aparece na lista de execução de processos em vez de `collectGreetings`.

E é isso! Agora estamos usando um módulo robusto curado pela comunidade em vez de código customizado de nível de protótipo para essa etapa no pipeline.

### Conclusão

Você agora sabe como:

- Encontrar e instalar módulos nf-core
- Avaliar os requisitos de um módulo nf-core
- Criar um mapa de metadados simples para usar com um módulo nf-core
- Integrar um módulo nf-core no seu fluxo de trabalho

### Qual é o próximo passo?

Aprender a adaptar seus módulos locais para seguir as convenções nf-core.
Também mostraremos como criar novos módulos nf-core a partir de um modelo usando as ferramentas nf-core.
