# Parte 4: Criar um módulo nf-core

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tradução assistida por IA - [saiba mais e sugira melhorias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Nesta quarta parte do curso de treinamento Hello nf-core, mostramos como criar um módulo nf-core aplicando as convenções principais que tornam os módulos portáveis e de fácil manutenção.

O projeto nf-core fornece um comando (`nf-core modules create`) que gera templates de módulos estruturados adequadamente de forma automática, semelhante ao que usamos para o fluxo de trabalho na Parte 2.
No entanto, para fins didáticos, vamos começar fazendo isso manualmente: transformando o módulo local `cowpy` em seu pipeline `core-hello` em um módulo no estilo nf-core passo a passo.
Depois disso, mostraremos como usar a criação de módulos baseada em template para trabalhar de forma mais eficiente no futuro.

??? info "Como começar a partir desta seção"

    Esta seção assume que você completou a [Parte 3: Usar um módulo nf-core](./03_use_module.md) e integrou o módulo `CAT_CAT` em seu pipeline.

    Se você não completou a Parte 3 ou quer começar do zero para esta parte, pode usar a solução `core-hello-part3` como ponto de partida.
    Execute estes comandos de dentro do diretório `hello-nf-core/`:

    ```bash
    cp -r solutions/core-hello-part3 core-hello
    cd core-hello
    ```

    Isso fornece um pipeline com o módulo `CAT_CAT` já integrado.
    Você pode testar que ele executa com sucesso executando o seguinte comando:

    ```bash
    nextflow run . --outdir core-hello-results -profile test,docker --validate_params false
    ```

---

## 1. Transformar `cowpy` em um módulo nf-core

Nesta seção, aplicaremos as convenções nf-core ao módulo local `cowpy` em seu pipeline `core-hello`, transformando-o em um módulo que segue os padrões da comunidade nf-core.

Este é o código atual para o módulo de processo `cowpy`:

```groovy title="core-hello/modules/local/cowpy.nf" linenums="1"
#!/usr/bin/env nextflow

// Generate ASCII art with cowpy (https://github.com/jeffbuttars/cowpy)
process cowpy {

    publishDir 'results', mode: 'copy'

    container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
    conda 'conda-forge::cowpy==1.1.5'

    input:
    path input_file
    val character

    output:
    path "cowpy-${input_file}"

    script:
    """
    cat $input_file | cowpy -c "$character" > cowpy-${input_file}
    """
}
```

Aplicaremos as seguintes convenções nf-core de forma incremental:

1. **Colocar o nome do processo em maiúsculas para `COWPY`** para seguir a convenção.
2. **Atualizar `COWPY` para usar tuplas de metadados** para propagar metadados de amostra através do fluxo de trabalho.
3. **Centralizar a configuração de argumentos da ferramenta com `ext.args`** para aumentar a versatilidade do módulo mantendo a interface mínima.
4. **Padronizar a nomenclatura de saída com `ext.prefix`** para promover consistência.
5. **Centralizar a configuração de publicação** para promover consistência.

Após cada etapa, executaremos o pipeline para testar que tudo funciona como esperado.

!!! warning "Diretório de trabalho"

    Certifique-se de estar no diretório `core-hello` (a raiz do seu pipeline) para todas as edições de arquivo e execuções de comando nesta seção.

    ```bash
    cd core-hello
    ```

### 1.1. Colocar o nome do processo em maiúsculas

Esta é puramente uma convenção estilística (não há justificativa técnica), mas como é a norma para módulos nf-core, vamos seguir.

Precisamos fazer três conjuntos de mudanças:

1. Atualizar o nome do processo no módulo
2. Atualizar a declaração de importação do módulo no cabeçalho do fluxo de trabalho
3. Atualizar a chamada do processo e a declaração emit no corpo do fluxo de trabalho

Vamos começar!

#### 1.1.1. Atualizar o nome do processo no módulo

Abra o arquivo do módulo `cowpy.nf` (em `core-hello/modules/local/`) e modifique o nome do processo para maiúsculas:

=== "Depois"

    ```groovy title="core-hello/modules/local/cowpy.nf" linenums="3" hl_lines="2"
    // Generate ASCII art with cowpy (https://github.com/jeffbuttars/cowpy)
    process COWPY {
    ```

=== "Antes"

    ```groovy title="core-hello/modules/local/cowpy.nf" linenums="3" hl_lines="2"
    // Generate ASCII art with cowpy (https://github.com/jeffbuttars/cowpy)
    process cowpy {
    ```

Neste caso, colocar em maiúsculas é completamente direto.

Se o nome do processo fosse composto de várias palavras, por exemplo, se tivéssemos um processo chamado MyCowpyTool originalmente em camel case, a convenção nf-core seria usar underscores para separá-las, resultando em MY_COWPY_TOOL.

#### 1.1.2. Atualizar a declaração de importação do módulo

Nomes de processos são sensíveis a maiúsculas/minúsculas, então agora que mudamos o nome do processo, precisamos atualizar a declaração de importação do módulo de acordo no cabeçalho do fluxo de trabalho de `hello.nf`:

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
    include { COWPY                  } from '../modules/local/cowpy/main.nf'
    include { CAT_CAT                } from '../modules/nf-core/cat/cat/main'
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
    include { cowpy                  } from '../modules/local/cowpy/main.nf'
    include { CAT_CAT                } from '../modules/nf-core/cat/cat/main'
    ```

Poderíamos usar um alias na declaração de importação para evitar ter que atualizar as chamadas ao processo, mas isso de certa forma anularia o propósito de adotar a convenção de maiúsculas.

#### 1.1.3. Atualizar a chamada do processo e a declaração emit

Então agora vamos atualizar as duas referências ao processo no bloco workflow de `hello.nf`:

=== "Depois"

    ```groovy title="core-hello/workflows/hello.nf" linenums="43" hl_lines="2 17"
    // generate ASCII art of the greetings with cowpy
    COWPY(CAT_CAT.out.file_out)

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name:  'hello_software_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    emit:
    cowpy_hellos   = COWPY.out.cowpy_output
    versions       = ch_versions
    ```

=== "Antes"

    ```groovy title="core-hello/workflows/hello.nf" linenums="43" hl_lines="2 17"
    // generate ASCII art of the greetings with cowpy
    cowpy(CAT_CAT.out.file_out)

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name:  'hello_software_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    emit:
    cowpy_hellos   = cowpy.out.cowpy_output
    versions       = ch_versions
    ```

Certifique-se de fazer **ambas** as mudanças, caso contrário você receberá um erro ao executar isso.

#### 1.1.4. Executar o pipeline para testá-lo

Vamos executar o fluxo de trabalho para testar que tudo está funcionando corretamente após essas mudanças.

```bash
nextflow run . --outdir core-hello-results -profile test,docker --validate_params false
```

??? success "Saída do comando"

    ```console
     N E X T F L O W   ~  version 25.04.3

    Launching `./main.nf` [elegant_plateau] DSL2 - revision: b9e9b3b8de

    Input/output options
      input                     : /workspaces/training/hello-nf-core/core-hello/assets/greetings.csv
      outdir                    : core-hello-results

    Institutional config options
      config_profile_name       : Test profile
      config_profile_description: Minimal test dataset to check pipeline function

    Generic options
      validate_params           : false
      trace_report_suffix       : 2026-01-06_04-51-29

    Core Nextflow options
      runName                   : elegant_plateau
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
    [7b/66ceb5] CORE_HELLO:HELLO:sayHello (3)       | 3 of 3 ✔
    [8e/1bafb9] CORE_HELLO:HELLO:convertToUpper (3) | 3 of 3 ✔
    [bb/203575] CORE_HELLO:HELLO:CAT_CAT (test)     | 1 of 1 ✔
    [39/715489] CORE_HELLO:HELLO:COWPY              | 1 of 1 ✔
    -[core/hello] Pipeline completed successfully-
    ```

Certo, isso funciona! Agora vamos passar a fazer mudanças mais substanciais.

### 1.2. Atualizar `COWPY` para usar tuplas de metadados

Na versão atual do pipeline `core-hello`, estamos extraindo o arquivo da tupla de saída de `CAT_CAT` para passar para `COWPY`, como mostrado na metade superior do diagrama abaixo.

<figure class="excalidraw">
    --8<-- "docs/en/docs/hello_nf-core/img/cowpy-inputs.svg"
</figure>

Seria melhor ter `COWPY` aceitando tuplas de metadados diretamente, permitindo que os metadados fluam através do fluxo de trabalho, como mostrado na metade inferior do diagrama.

Para esse fim, precisaremos fazer as seguintes mudanças:

1. Atualizar as definições de entrada e saída
2. Atualizar a chamada do processo no fluxo de trabalho
3. Atualizar o bloco emit no fluxo de trabalho

Uma vez que tivermos feito tudo isso, executaremos o pipeline para testar que tudo ainda funciona como antes.

#### 1.2.1. Atualizar as definições de entrada e saída

Retorne ao arquivo do módulo `cowpy.nf` e modifique-o para aceitar tuplas de metadados como mostrado abaixo.

=== "Depois"

    ```groovy title="core-hello/modules/local/cowpy.nf" linenums="11" hl_lines="2 6"
        input:
            tuple val(meta), path(input_file)
            val character

        output:
            tuple val(meta), path("cowpy-${input_file}"), emit: cowpy_output
    ```

=== "Antes"

    ```groovy title="core-hello/modules/local/cowpy.nf" linenums="11" hl_lines="2 6"
        input:
            path input_file
            val character

        output:
            path "cowpy-${input_file}"
    ```

Como você pode ver, mudamos tanto a **entrada principal** quanto a **saída** para uma tupla que segue o padrão `tuple val(meta), path(input_file)` introduzido na Parte 3 deste treinamento.
Para a saída, também aproveitamos esta oportunidade para adicionar `emit: cowpy_output` a fim de dar ao canal de saída um nome descritivo.

Agora que mudamos o que o processo espera, precisamos atualizar o que fornecemos a ele na chamada do processo.

#### 1.2.2. Atualizar a chamada do processo no fluxo de trabalho

A boa notícia é que essa mudança simplificará a chamada do processo.
Agora que a saída de `CAT_CAT` e a entrada de `COWPY` têm a mesma 'forma', ou seja, ambas consistem em uma estrutura `tuple val(meta), path(input_file)`, podemos simplesmente conectá-las diretamente em vez de ter que extrair o arquivo explicitamente da saída do processo `CAT_CAT`.

Abra o arquivo de fluxo de trabalho `hello.nf` (em `core-hello/workflows/`) e atualize a chamada para `COWPY` como mostrado abaixo.

=== "Depois"

    ```groovy title="core-hello/workflows/hello.nf" linenums="43" hl_lines="2"
        // generate ASCII art of the greetings with cowpy
        COWPY(CAT_CAT.out.file_out, params.character)
    ```

=== "Antes"

    ```groovy title="core-hello/workflows/hello.nf" linenums="43" hl_lines="1-2 5"
        // extract the file from the tuple since cowpy doesn't use metadata yet
        ch_for_cowpy = CAT_CAT.out.file_out.map{ meta, file -> file }

        // generate ASCII art of the greetings with cowpy
        COWPY(ch_for_cowpy, params.character)
    ```

Agora chamamos `COWPY` em `CAT_CAT.out.file_out` diretamente.

Como resultado, não precisamos mais construir o canal `ch_for_cowpy`, então essa linha (e sua linha de comentário) pode ser deletada inteiramente.

#### 1.2.3. Atualizar o bloco emit no fluxo de trabalho

Como `COWPY` agora emite uma saída nomeada, `cowpy_output`, podemos atualizar o bloco `emit:` do fluxo de trabalho `hello.nf` para usar isso.

=== "Depois"

    ```groovy title="core-hello/workflows/hello.nf" linenums="60" hl_lines="2"
        emit:
        cowpy_hellos   = COWPY.out.cowpy_output
        versions       = ch_versions
    ```

=== "Antes"

    ```groovy title="core-hello/workflows/hello.nf" linenums="60" hl_lines="2"
        emit:
        cowpy_hellos   = COWPY.out
        versions       = ch_versions
    ```

Isso tecnicamente não é necessário, mas é uma boa prática referir-se a saídas nomeadas sempre que possível.

#### 1.2.4. Executar o pipeline para testá-lo

Vamos executar o fluxo de trabalho para testar que tudo está funcionando corretamente após essas mudanças.

```bash
nextflow run . --outdir core-hello-results -profile test,docker --validate_params false
```

??? success "Saída do comando"

    ```console
     N E X T F L O W   ~  version 25.04.3

    Launching `./main.nf` [modest_saha] DSL2 - revision: b9e9b3b8de

    Downloading plugin nf-schema@2.5.1
    Input/output options
      input                     : /workspaces/training/hello-nf-core/core-hello/assets/greetings.csv
      outdir                    : core-hello-results

    Institutional config options
      config_profile_name       : Test profile
      config_profile_description: Minimal test dataset to check pipeline function

    Generic options
      validate_params           : false
      trace_report_suffix       : 2025-12-27_06-16-55

    Core Nextflow options
      runName                   : modest_saha
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
    [a8/447993] CORE_HELLO:HELLO:sayHello (3)       [100%] 3 of 3 ✔
    [00/1fc59c] CORE_HELLO:HELLO:convertToUpper (3) [100%] 3 of 3 ✔
    [57/ac800d] CORE_HELLO:HELLO:CAT_CAT (test)     [100%] 1 of 1 ✔
    [b7/092f2b] CORE_HELLO:HELLO:COWPY              [100%] 1 of 1 ✔
    -[core/hello] Pipeline completed successfully-
    ```

O pipeline deve executar com sucesso, com os metadados agora fluindo de `CAT_CAT` através de `COWPY`.

Isso completa o que precisávamos fazer para que `COWPY` manipule tuplas de metadados.
Agora, vamos ver o que mais podemos fazer para aproveitar os padrões de módulos nf-core.

### 1.3. Centralizar a configuração de argumentos da ferramenta com `ext.args`

Em seu estado atual, o processo `COWPY` espera receber um valor para o parâmetro `character`.
Como resultado, temos que fornecer um valor toda vez que chamamos o processo, mesmo que estejamos satisfeitos com os padrões definidos pela ferramenta.
Para `COWPY` isso admitidamente não é um grande problema, mas para ferramentas com muitos parâmetros opcionais, pode se tornar bastante trabalhoso.

O projeto nf-core recomenda usar um recurso do Nextflow chamado [`ext.args`](https://www.nextflow.io/docs/latest/reference/process.html#ext) para gerenciar argumentos de ferramentas de forma mais conveniente através de arquivos de configuração.

Em vez de declarar entradas de processo para cada opção da ferramenta, você escreve o módulo para referenciar `ext.args` na construção de sua linha de comando.
Então é apenas uma questão de configurar a variável `ext.args` para conter os argumentos e valores que você quer usar no arquivo `modules.config`, que consolida detalhes de configuração para todos os módulos.
O Nextflow adicionará esses argumentos com seus valores na linha de comando da ferramenta em tempo de execução.

Vamos aplicar essa abordagem ao módulo `COWPY`.
Vamos precisar fazer as seguintes mudanças:

1. Atualizar o módulo `COWPY`
2. Configurar `ext.args` no arquivo `modules.config`
3. Atualizar o fluxo de trabalho `hello.nf`

Uma vez que tivermos feito tudo isso, executaremos o pipeline para testar que tudo ainda funciona como antes.

#### 1.3.1. Atualizar o módulo `COWPY`

Vamos fazer isso.
Abra o arquivo do módulo `cowpy.nf` (em `core-hello/modules/local/`) e modifique-o para referenciar `ext.args` como mostrado abaixo.

=== "Depois"

    ```groovy title="modules/local/cowpy.nf" linenums="1" hl_lines="18 20"
    #!/usr/bin/env nextflow

    // Generate ASCII art with cowpy (https://github.com/jeffbuttars/cowpy)
    process COWPY {

        publishDir 'results', mode: 'copy'

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
        conda 'conda-forge::cowpy==1.1.5'

        input:
            tuple val(meta), path(input_file)

        output:
            tuple val(meta), path("cowpy-${input_file}"), emit: cowpy_output

        script:
        def args = task.ext.args ?: ''
        """
        cat $input_file | cowpy $args > cowpy-${input_file}
        """
    }
    ```

=== "Antes"

    ```groovy title="core-hello/modules/local/cowpy.nf" linenums="1" hl_lines="13 20"
    #!/usr/bin/env nextflow

    // Generate ASCII art with cowpy (https://github.com/jeffbuttars/cowpy)
    process COWPY {

        publishDir 'results', mode: 'copy'

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
        conda 'conda-forge::cowpy==1.1.5'

        input:
            tuple val(meta), path(input_file)
            val character

        output:
            tuple val(meta), path("cowpy-${input_file}"), emit: cowpy_output

        script:
        """
        cat $input_file | cowpy -c "$character" > cowpy-${input_file}
        """
    }
    ```

Você pode ver que fizemos três mudanças.

1. **No bloco `input:`, removemos a entrada `val character`.**
   Daqui em diante, forneceremos esse argumento através da configuração `ext.args` conforme descrito mais abaixo.

2. **No bloco `script:`, adicionamos a linha `def args = task.ext.args ?: ''`.**
   Essa linha usa o operador `?:` para determinar o valor da variável `args`: o conteúdo de `task.ext.args` se não estiver vazio, ou uma string vazia se estiver.
   Note que embora geralmente nos refiramos a `ext.args`, este código deve referenciar `task.ext.args` para extrair a configuração `ext.args` no nível do módulo.

3. **Na linha de comando, substituímos `-c "$character"` por `$args`.**
   É aqui que o Nextflow injetará quaisquer argumentos de ferramenta definidos em `ext.args` no arquivo `modules.config`.

Como resultado, a interface do módulo agora é mais simples: ela espera apenas as entradas essenciais de metadados e arquivo.

!!! note

    O operador `?:` é frequentemente chamado de 'operador Elvis' porque parece com um rosto de Elvis Presley de lado, com o caractere `?` simbolizando a onda em seu cabelo.

#### 1.3.2. Configurar `ext.args` no arquivo `modules.config`

Agora que tiramos a declaração `character` do módulo, temos que adicioná-la a `ext.args` no arquivo de configuração `modules.config`.

Especificamente, vamos adicionar este pequeno pedaço de código ao bloco `process {}`:

```groovy title="Código a adicionar"
withName: 'COWPY' {
    ext.args = { "-c ${params.character}" }
}
```

A sintaxe `withName:` atribui esta configuração apenas ao processo `COWPY`, e `ext.args = { "-c ${params.character}" }` simplesmente compõe uma string que incluirá o valor do parâmetro `character`.
Note o uso de chaves, que dizem ao Nextflow para avaliar o valor do parâmetro em tempo de execução.

Faz sentido? Vamos adicionar.

Abra `conf/modules.config` e adicione o código de configuração dentro do bloco `process {}` como mostrado abaixo.

=== "Depois"

    ```groovy title="core-hello/conf/modules.config" linenums="13" hl_lines="8-10"
    process {
        publishDir = [
            path: { "${params.outdir}/${task.process.tokenize(':')[-1].tokenize('_')[0].toLowerCase()}" },
            mode: params.publish_dir_mode,
            saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
        ]

        withName: 'COWPY' {
            ext.args = { "-c ${params.character}" }
        }
    }
    ```

=== "Antes"

    ```groovy title="core-hello/conf/modules.config" linenums="13"
    process {
        publishDir = [
            path: { "${params.outdir}/${task.process.tokenize(':')[-1].tokenize('_')[0].toLowerCase()}" },
            mode: params.publish_dir_mode,
            saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
        ]
    }
    ```

Esperamos que você possa imaginar ter todos os módulos em um pipeline com seus `ext.args` especificados neste arquivo, com os seguintes benefícios:

- A **interface do módulo permanece simples** - Ela aceita apenas as entradas essenciais de metadados e arquivo
- O **pipeline ainda expõe `params.character`** - Os usuários finais ainda podem configurá-lo como antes
- O **módulo agora é portável** - Pode ser reutilizado em outros pipelines sem esperar um nome de parâmetro específico
- A configuração é **centralizada** em `modules.config`, mantendo a lógica do fluxo de trabalho limpa

Ao usar o arquivo `modules.config` como o lugar onde todos os pipelines centralizam a configuração por módulo, tornamos nossos módulos mais reutilizáveis em diferentes pipelines.

#### 1.3.3. Atualizar o fluxo de trabalho `hello.nf`

Como o módulo `COWPY` não requer mais o parâmetro `character` como entrada, precisamos atualizar a chamada do fluxo de trabalho de acordo.

Abra o arquivo de fluxo de trabalho `hello.nf` (em `core-hello/workflows/`) e atualize a chamada para `COWPY` como mostrado abaixo.

=== "Depois"

    ```groovy title="core-hello/workflows/hello.nf" linenums="39" hl_lines="2"
        // generate ASCII art of the greetings with cowpy
        COWPY(CAT_CAT.out.file_out)
    ```

=== "Antes"

    ```groovy title="core-hello/workflows/hello.nf" linenums="39" hl_lines="2"
        // generate ASCII art of the greetings with cowpy
        COWPY(CAT_CAT.out.file_out, params.character)
    ```

O código do fluxo de trabalho agora está mais limpo: não precisamos passar `params.character` diretamente para o processo.
A interface do módulo é mantida mínima, tornando-o mais portável, enquanto o pipeline ainda fornece a opção explícita através da configuração.

#### 1.3.4. Executar o pipeline para testá-lo

Vamos testar que o fluxo de trabalho ainda funciona como esperado, especificando um caractere diferente para verificar que a configuração `ext.args` está funcionando.

Execute este comando usando `kosh`, uma das opções mais... enigmáticas:

```bash
nextflow run . --outdir core-hello-results -profile test,docker --validate_params false --character kosh
```

??? success "Saída do comando"

    ```console
     N E X T F L O W   ~  version 25.04.3

    Launching `./main.nf` [exotic_planck] DSL2 - revision: b9e9b3b8de

    Input/output options
      input                     : /workspaces/training/hello-nf-core/core-hello/assets/greetings.csv
      outdir                    : core-hello-results

    Institutional config options
      config_profile_name       : Test profile
      config_profile_description: Minimal test dataset to check pipeline function

    Generic options
      validate_params           : false
      trace_report_suffix       : 2025-12-27_06-23-13

    Core Nextflow options
      runName                   : exotic_planck
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
    [13/9e3c0e] CORE_HELLO:HELLO:sayHello (3)       [100%] 3 of 3 ✔
    [e2/5b0ee5] CORE_HELLO:HELLO:convertToUpper (3) [100%] 3 of 3 ✔
    [b6/4fb569] CORE_HELLO:HELLO:CAT_CAT (test)     [100%] 1 of 1 ✔
    [38/eb29ea] CORE_HELLO:HELLO:COWPY              [100%] 1 of 1 ✔
    -[core/hello] Pipeline completed successfully-
    ```

Isso deve executar com sucesso como anteriormente.

Vamos verificar que a configuração `ext.args` funcionou verificando a saída.
Encontre a saída no navegador de arquivos ou use o hash da tarefa (a parte `38/eb29ea` no exemplo acima) para olhar o arquivo de saída:

```bash
cat work/38/eb29ea*/cowpy-test.txt
```

??? success "Saída do comando"

    ```console
    _________
    / HELLO   \
    | HOLà    |
    \ BONJOUR /
    ---------
        \
        \
          \
      ___       _____     ___
    /   \     /    /|   /   \
    |     |   /    / |  |     |
    |     |  /____/  |  |     |
    |     |  |    |  |  |     |
    |     |  | {} | /   |     |
    |     |  |____|/    |     |
    |     |    |==|     |     |
    |      \___________/      |
    |                         |
    |                         |
    ```

Você deve ver a arte ASCII exibida com o caractere `kosh`, confirmando que a configuração `ext.args` funcionou!

??? info "(Opcional) Inspecionar o arquivo de comando"

    Se você quiser ver exatamente como a configuração foi aplicada, pode inspecionar o arquivo `.command.sh`:

    ```bash
    cat work/38/eb29ea*/.command.sh
    ```

    Você verá o comando `cowpy` com o argumento `-c kosh`:

    ```console
    #!/usr/bin/env bash
    ...
    cat test.txt | cowpy -c kosh > cowpy-test.txt
    ```

    Isso mostra que o arquivo `.command.sh` foi gerado corretamente com base na configuração `ext.args`.

Reserve um momento para pensar sobre o que alcançamos aqui.
Essa abordagem mantém a interface do módulo focada em dados essenciais (arquivos, metadados e quaisquer parâmetros obrigatórios por amostra), enquanto opções que controlam o comportamento da ferramenta são tratadas separadamente através da configuração.

Isso pode parecer desnecessário para uma ferramenta simples como `cowpy`, mas pode fazer uma grande diferença para ferramentas de análise de dados que têm muitos argumentos opcionais.

Para resumir os benefícios desta abordagem:

- **Interface limpa**: O módulo foca em entradas de dados essenciais (metadados e arquivos)
- **Flexibilidade**: Usuários podem especificar argumentos de ferramenta via configuração, incluindo valores específicos por amostra
- **Consistência**: Todos os módulos nf-core seguem este padrão
- **Portabilidade**: Módulos podem ser reutilizados sem opções de ferramenta codificadas
- **Sem mudanças no fluxo de trabalho**: Adicionar ou mudar opções de ferramenta não requer atualizar o código do fluxo de trabalho

!!! note

    O sistema `ext.args` tem capacidades adicionais poderosas não cobertas aqui, incluindo alternar valores de argumentos dinamicamente com base em metadados. Veja as [especificações de módulos nf-core](https://nf-co.re/docs/guidelines/components/modules) para mais detalhes.

### 1.4. Padronizar a nomenclatura de saída com `ext.prefix`

Agora que demos ao processo `COWPY` acesso ao metamap, podemos começar a aproveitar outro padrão útil do nf-core: nomear arquivos de saída com base em metadados.

Aqui vamos usar um recurso do Nextflow chamado `ext.prefix` que nos permitirá padronizar a nomenclatura de arquivos de saída entre módulos usando `meta.id` (o identificador incluído no metamap), enquanto ainda podemos configurar módulos individualmente se desejado.

Isso será semelhante ao que fizemos com `ext.args`, com algumas diferenças que detalharemos conforme avançamos.

Vamos aplicar essa abordagem ao módulo `COWPY`.
Vamos precisar fazer as seguintes mudanças:

1. Atualizar o módulo `COWPY`
2. Configurar `ext.prefix` no arquivo `modules.config`

(Nenhuma mudança necessária no fluxo de trabalho.)

Uma vez que tivermos feito isso, executaremos o pipeline para testar que tudo ainda funciona como antes.

#### 1.4.1. Atualizar o módulo `COWPY`

Abra o arquivo do módulo `cowpy.nf` (em `core-hello/modules/local/`) e modifique-o para referenciar `ext.prefix` como mostrado abaixo.

=== "Depois"

    ```groovy title="core-hello/modules/local/cowpy.nf" linenums="1" hl_lines="2 6 8"
        output:
            tuple val(meta), path("${prefix}.txt"), emit: cowpy_output

        script:
        def args = task.ext.args ?: ''
        prefix = task.ext.prefix ?: "${meta.id}"
        """
        cat $input_file | cowpy $args > ${prefix}.txt
        """
    }
    ```

=== "Antes"

    ```groovy title="core-hello/modules/local/cowpy.nf" linenums="1" hl_lines="2 7"
        output:
            tuple val(meta), path("cowpy-${input_file}"), emit: cowpy_output

        script:
        def args = task.ext.args ?: ''
        """
        cat $input_file | cowpy $args > cowpy-${input_file}
        """
    }
    ```

Você pode ver que fizemos três mudanças.

1. **No bloco `script:`, adicionamos a linha `prefix = task.ext.prefix ?: "${meta.id}"`.**
   Essa linha usa o operador `?:` para determinar o valor da variável `prefix`: o conteúdo de `task.ext.prefix` se não estiver vazio, ou o identificador do metamap (`meta.id`) se estiver.
   Note que embora geralmente nos refiramos a `ext.prefix`, este código deve referenciar `task.ext.prefix` para extrair a configuração `ext.prefix` no nível do módulo.

2. **Na linha de comando, substituímos `cowpy-${input_file}` por `${prefix}.txt`.**
   É aqui que o Nextflow injetará o valor de `prefix` determinado pela linha acima.

3. **No bloco `output:`, substituímos `path("cowpy-${input_file}")` por `path("${prefix}.txt")`.**
   Isso simplesmente reitera qual será o caminho do arquivo de acordo com o que está escrito na linha de comando.

Como resultado, o nome do arquivo de saída agora é construído usando um padrão sensato (o identificador do metamap) combinado com a extensão de formato de arquivo apropriada.

#### 1.4.2. Configurar `ext.prefix` no arquivo `modules.config`

Neste caso, o padrão sensato não é suficientemente expressivo para nosso gosto; queremos usar um padrão de nomenclatura personalizado que inclua o nome da ferramenta, `cowpy-<id>.txt`, como tínhamos antes.

Faremos isso configurando `ext.prefix` em `modules.config`, assim como fizemos para o parâmetro `character` com `ext.args`, exceto que desta vez o bloco `withName: 'COWPY' {}` já existe, e só precisamos adicionar a seguinte linha:

```groovy title="Código a adicionar"
ext.prefix = { "cowpy-${meta.id}" }
```

Isso comporá a string que queremos.
Note que mais uma vez usamos chaves, desta vez para dizer ao Nextflow para avaliar o valor de `meta.id` em tempo de execução.

Vamos adicionar.

Abra `conf/modules.config` e adicione o código de configuração dentro do bloco `process {}` como mostrado abaixo.

=== "Depois"

    ```groovy title="core-hello/conf/modules.config" linenums="21" hl_lines="3"
        withName: 'COWPY' {
            ext.args = { "-c ${params.character}" }
            ext.prefix = { "cowpy-${meta.id}" }
        }
    ```

=== "Antes"

    ```groovy title="core-hello/conf/modules.config" linenums="21"
        withName: 'COWPY' {
            ext.args = { "-c ${params.character}" }
        }
    ```

Caso você esteja se perguntando, o closure `ext.prefix` tem acesso à peça correta de metadados porque a configuração é avaliada no contexto da execução do processo, onde os metadados estão disponíveis.

#### 1.4.3. Executar o pipeline para testá-lo

Vamos testar que o fluxo de trabalho ainda funciona como esperado.

```bash
nextflow run . --outdir core-hello-results -profile test,docker --validate_params false
```

??? success "Saída do comando"

    ```console
     N E X T F L O W   ~  version 25.04.3

    Launching `./main.nf` [admiring_turing] DSL2 - revision: b9e9b3b8de

    Input/output options
      input                     : /workspaces/training/hello-nf-core/core-hello/assets/greetings.csv
      outdir                    : core-hello-results

    Institutional config options
      config_profile_name       : Test profile
      config_profile_description: Minimal test dataset to check pipeline function

    Generic options
      validate_params           : false
      trace_report_suffix       : 2025-12-27_06-29-02

    Core Nextflow options
      runName                   : admiring_turing
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
    [b2/e08524] CORE_HELLO:HELLO:sayHello (3)       [100%] 3 of 3 ✔
    [13/88939f] CORE_HELLO:HELLO:convertToUpper (3) [100%] 3 of 3 ✔
    [23/4554e1] CORE_HELLO:HELLO:CAT_CAT (test)     [100%] 1 of 1 ✔
    [a3/c6cbe9] CORE_HELLO:HELLO:COWPY              [100%] 1 of 1 ✔
    -[core/hello] Pipeline completed successfully-
    ```

Dê uma olhada na saída no diretório de resultados.
Você deve ver o arquivo de saída do cowpy com a mesma nomenclatura de antes: `cowpy-test.txt`, baseado no nome de lote padrão.

??? abstract "Conteúdo do diretório"

    ```console hl_lines="3"
    results
    ├── Bonjour-output.txt
    ├── cowpy-test.txt
    ├── Hello-output.txt
    ├── Holà-output.txt
    ├── UPPER-Bonjour-output.txt
    ├── UPPER-Hello-output.txt
    └── UPPER-Holà-output.txt
    ```

Sinta-se à vontade para mudar a configuração `ext.prefix` em `conf/modules.config` para se convencer de que você pode mudar o padrão de nomenclatura sem ter que fazer quaisquer mudanças no código do módulo ou do fluxo de trabalho.

Alternativamente, você também pode tentar executar isso novamente com um parâmetro `--batch` diferente especificado na linha de comando para se convencer de que essa parte ainda é personalizável em tempo real.

Isso demonstra como `ext.prefix` permite que você mantenha sua convenção de nomenclatura preferida mantendo a interface do módulo flexível.

Para resumir os benefícios desta abordagem:

- **Nomenclatura padronizada**: Arquivos de saída são tipicamente nomeados usando IDs de amostra dos metadados
- **Configurável**: Usuários podem sobrescrever a nomenclatura padrão se necessário
- **Consistente**: Todos os módulos nf-core seguem este padrão
- **Previsível**: Fácil saber como os arquivos de saída serão chamados

Muito bom, certo?
Bem, há mais uma mudança importante que precisamos fazer para melhorar nosso módulo para se adequar às diretrizes nf-core.

### 1.5. Centralizar a configuração de publicação

Você pode ter notado que temos publicado saídas em dois diretórios diferentes:

- **`results`** — O diretório de saída original que temos usado desde o início para nossos módulos locais, definido individualmente usando diretivas `publishDir` por módulo;
- **`core-hello-results`** — O diretório de saída definido com `--outdir` na linha de comando, que tem recebido os logs nf-core e os resultados publicados por `CAT_CAT`.

Isso é confuso e subótimo; seria melhor ter um local para tudo.
Claro, poderíamos entrar em cada um de nossos módulos locais e atualizar a diretiva `publishDir` manualmente para usar o diretório `core-hello-results`, mas e na próxima vez que decidirmos mudar o diretório de saída?

Ter módulos individuais tomando decisões de publicação claramente não é o caminho a seguir, especialmente em um mundo onde o mesmo módulo pode ser usado em muitos pipelines diferentes, por pessoas que têm necessidades ou preferências diferentes.
Queremos ser capazes de controlar onde as saídas são publicadas no nível da configuração do fluxo de trabalho.

"Ei," você pode dizer, "`CAT_CAT` está enviando suas saídas para o `--outdir`. Talvez devêssemos copiar sua diretiva `publishDir`?"

Sim, essa é uma ótima ideia.

Exceto que ele não tem uma diretiva `publishDir`. (Vá em frente, olhe o código do módulo.)

Isso porque pipelines nf-core centralizam o controle no nível do fluxo de trabalho configurando `publishDir` em `conf/modules.config` em vez de em módulos individuais.
Especificamente, o template nf-core declara uma diretiva `publishDir` padrão (com uma estrutura de diretório predefinida) que se aplica a todos os módulos a menos que uma diretiva de sobrescrita seja fornecida.

Isso não soa incrível? Poderia ser que para aproveitar esta diretiva padrão, tudo que precisamos fazer é remover a diretiva `publishDir` atual de nossos módulos locais?

Vamos tentar isso em `COWPY` para ver o que acontece, então olharemos o código para a configuração padrão para entender como funciona.

Finalmente, demonstraremos como sobrescrever o comportamento padrão se desejado.

#### 1.5.1. Remover a diretiva `publishDir` de `COWPY`

Vamos fazer isso.
Abra o arquivo do módulo `cowpy.nf` (em `core-hello/modules/local/`) e remova a diretiva `publishDir` como mostrado abaixo.

=== "Depois"

    ```groovy title="core-hello/modules/local/cowpy.nf (trecho)" linenums="1"
    #!/usr/bin/env nextflow

    // Generate ASCII art with cowpy (https://github.com/jeffbuttars/cowpy)
    process COWPY {

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
        conda 'conda-forge::cowpy==1.1.5'
    ```

=== "Antes"

    ```groovy title="core-hello/modules/local/cowpy.nf (trecho)" linenums="1" hl_lines="6"
    #!/usr/bin/env nextflow

    // Generate ASCII art with cowpy (https://github.com/jeffbuttars/cowpy)
    process COWPY {

        publishDir 'results', mode: 'copy'

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
        conda 'conda-forge::cowpy==1.1.5'

    ```

É isso!

#### 1.5.2. Executar o pipeline para testá-lo

Vamos dar uma olhada no que acontece se executarmos o pipeline agora.

```bash
nextflow run . --outdir core-hello-results -profile test,docker --validate_params false
```

??? success "Saída do comando"

    ```console
     N E X T F L O W   ~  version 25.04.3

    Launching `./main.nf` [silly_caravaggio] DSL2 - revision: b9e9b3b8de

    Input/output options
      input                     : /workspaces/training/hello-nf-core/core-hello/assets/greetings.csv
      outdir                    : core-hello-results

    Institutional config options
      config_profile_name       : Test profile
      config_profile_description: Minimal test dataset to check pipeline function

    Generic options
      validate_params           : false
      trace_report_suffix       : 2025-12-27_06-35-56

    Core Nextflow options
      runName                   : silly_caravaggio
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
    [db/39978e] CORE_HELLO:HELLO:sayHello (3)       [100%] 3 of 3 ✔
    [b5/bf6a8d] CORE_HELLO:HELLO:convertToUpper (3) [100%] 3 of 3 ✔
    [b7/c61842] CORE_HELLO:HELLO:CAT_CAT (test)     [100%] 1 of 1 ✔
    [46/5839d6] CORE_HELLO:HELLO:COWPY              [100%] 1 of 1 ✔
    -[core/hello] Pipeline completed successfully-
    ```

Dê uma olhada no seu diretório de trabalho atual.
Agora o `core-hello-results` também contém as saídas do módulo `COWPY`.

??? abstract "Conteúdo do diretório"

    ```console hl_lines="4-5"
    core-hello-results/
    ├── cat
    │   └── test.txt
    ├── cowpy
    │   └── cowpy-test.txt
    └── pipeline_info
        ├── execution_report_2025-12-27_06-16-55.html
        ├── execution_report_2025-12-27_06-23-13.html
        ├── execution_report_2025-12-27_06-29-02.html
        ├── execution_report_2025-12-27_06-35-56.html
        ├── execution_timeline_2025-12-27_06-16-55.html
        ├── execution_timeline_2025-12-27_06-23-13.html
        ├── execution_timeline_2025-12-27_06-29-02.html
        ├── execution_timeline_2025-12-27_06-35-56.html
        ├── execution_trace_2025-12-27_06-16-55.txt
        ├── execution_trace_2025-12-27_06-23-13.txt
        ├── execution_trace_2025-12-27_06-29-02.txt
        ├── execution_trace_2025-12-27_06-35-56.txt
        ├── hello_software_versions.yml
        ├── params_2025-12-27_06-17-00.json
        ├── params_2025-12-27_06-23-17.json
        ├── params_2025-12-27_06-29-07.json
        ├── params_2025-12-27_06-36-01.json
        ├── pipeline_dag_2025-12-27_06-16-55.html
        ├── pipeline_dag_2025-12-27_06-23-13.html
        ├── pipeline_dag_2025-12-27_06-29-02.html
        └── pipeline_dag_2025-12-27_06-35-56.html
    ```

Você pode ver que o Nextflow criou esta hierarquia de diretórios baseada nos nomes do fluxo de trabalho e do módulo.

O código responsável vive no arquivo `conf/modules.config`.
Esta é a configuração `publishDir` padrão que faz parte do template nf-core e se aplica a todos os processos:

```groovy
process {
    publishDir = [
        path: { "${params.outdir}/${task.process.tokenize(':')[-1].tokenize('_')[0].toLowerCase()}" },
        mode: params.publish_dir_mode,
        saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
    ]
}
```

Isso pode parecer complicado, então vamos olhar cada um dos três componentes:

- **`path:`** Determina o diretório de saída baseado no nome do processo.
  O nome completo de um processo contido em `task.process` inclui a hierarquia de importações de fluxo de trabalho e módulo (como `CORE_HELLO:HELLO:CAT_CAT`).
  As operações `tokenize` removem essa hierarquia para obter apenas o nome do processo, então pegam a primeira parte antes de qualquer underscore (se aplicável), e convertem para minúsculas.
  Isso é o que determina que os resultados de `CAT_CAT` sejam publicados em `${params.outdir}/cat/`.
- **`mode:`** Controla como os arquivos são publicados (cópia, link simbólico, etc.).
  Isso é configurável através do parâmetro `params.publish_dir_mode`.
- **`saveAs:`** Filtra quais arquivos publicar.
  Este exemplo exclui arquivos `versions.yml` retornando `null` para eles, impedindo que sejam publicados.

Isso fornece uma lógica consistente para organizar saídas.

A saída fica ainda melhor quando todos os módulos em um pipeline adotam esta convenção, então sinta-se à vontade para ir deletar as diretivas `publishDir` dos outros módulos em seu pipeline.
Este padrão será aplicado mesmo a módulos que não modificamos explicitamente para seguir as diretrizes nf-core.

Dito isso, você pode decidir que quer organizar suas entradas de forma diferente, e a boa notícia é que é fácil fazer isso.

#### 1.5.3. Sobrescrever o padrão

Para sobrescrever a diretiva `publishDir` padrão, você pode simplesmente adicionar suas próprias diretivas ao arquivo `conf/modules.config`.

Por exemplo, você poderia sobrescrever o padrão para um único processo usando o seletor `withName:`, como neste exemplo onde adicionamos uma diretiva `publishDir` personalizada para o processo 'COWPY'.

```groovy title="core-hello/conf/modules.config" linenums="13" hl_lines="8-10"
process {
    publishDir = [
        path: { "${params.outdir}/${task.process.tokenize(':')[-1].tokenize('_')[0].toLowerCase()}" },
        mode: params.publish_dir_mode,
        saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
    ]

    withName: 'COWPY' {
        ext.args = { "-c ${params.character}" }
        publishDir = [
            path: 'my_custom_results'
        ]
    }
}
```

Não vamos realmente fazer essa mudança, mas sinta-se à vontade para brincar com isso e ver que lógica você pode implementar.

O ponto é que este sistema permite dar a você o melhor dos dois mundos: consistência por padrão e a flexibilidade para personalizar a configuração sob demanda.

Para resumir, você obtém:

- **Fonte única de verdade**: Toda configuração de publicação vive em `modules.config`
- **Padrão útil**: Processos funcionam imediatamente sem configuração por módulo
- **Personalização fácil**: Sobrescreva o comportamento de publicação na configuração, não no código do módulo
- **Módulos portáveis**: Módulos não codificam locais de saída

Isso completa o conjunto de recursos de módulos nf-core que você deve absolutamente aprender a usar, mas há outros que você pode ler sobre nas [especificações de módulos nf-core](https://nf-co.re/docs/guidelines/components/modules).

### Conclusão

Você agora sabe como adaptar módulos locais para seguir as convenções nf-core:

- Projete seus módulos para aceitar e propagar tuplas de metadados;
- Use `ext.args` para manter as interfaces de módulos mínimas e portáveis;
- Use `ext.prefix` para nomenclatura de arquivo de saída configurável e padronizada;
- Adote a diretiva `publishDir` centralizada padrão para uma estrutura de diretório de resultados consistente.

### O que vem a seguir?

Aprenda como usar as ferramentas integradas baseadas em template do nf-core para criar módulos da maneira fácil.

---

## 2. Criar um módulo com as ferramentas nf-core

Agora que você aprendeu os padrões de módulos nf-core aplicando-os manualmente, vamos ver como você criaria módulos na prática.

### 2.1. Gerar um scaffold de módulo a partir de um template

Semelhante ao que existe para criar pipelines, o projeto nf-core fornece ferramentas para gerar módulos estruturados adequadamente baseados em um template, com todos esses padrões incorporados desde o início.

#### 2.1.1. Executar o comando de criação de módulo

O comando `nf-core modules create` gera um template de módulo que já segue todas as convenções que você aprendeu.

Vamos criar uma nova versão do módulo `COWPY` com um template mínimo executando este comando:

```bash
nf-core modules create --empty-template COWPY
```

A flag `--empty-template` cria um template inicial limpo sem código extra, facilitando ver a estrutura essencial.

O comando executa interativamente, guiando você através da configuração.
Ele automaticamente procura informações da ferramenta em repositórios de pacotes como Bioconda e bio.tools para pré-popular metadados.

Você será solicitado para várias opções de configuração:

- **Informações do autor**: Seu nome de usuário do GitHub para atribuição
- **Label de recurso**: Um conjunto predefinido de requisitos computacionais.
  O projeto nf-core fornece labels padrão como `process_single` para ferramentas leves e `process_high` para as exigentes.
  Essas labels ajudam a gerenciar alocação de recursos em diferentes ambientes de execução.
- **Requisito de metadados**: Se o módulo precisa de informações específicas de amostra via um map `meta` (geralmente sim para módulos de processamento de dados).

A ferramenta lida com a complexidade de encontrar informações de pacotes e configurar a estrutura, permitindo que você foque em implementar a lógica específica da ferramenta.

#### 2.1.2. Examinar o scaffold do módulo

A ferramenta cria uma estrutura de módulo completa em `modules/local/` (ou `modules/nf-core/` se você estiver no repositório nf-core/modules):

??? abstract "Conteúdo do diretório"

    ```console
    modules/local/cowpy
    ├── environment.yml
    ├── main.nf
    ├── meta.yml
    └── tests
        └── main.nf.test
    ```

Cada arquivo serve um propósito específico:

- **`main.nf`**: Definição do processo com todos os padrões nf-core incorporados
- **`meta.yml`**: Documentação do módulo descrevendo entradas, saídas e a ferramenta
- **`environment.yml`**: Especificação do ambiente Conda para dependências
- **`tests/main.nf.test`**: Casos de teste nf-test para validar que o módulo funciona

!!! tip "Saiba mais sobre testes"

    O arquivo de teste gerado usa nf-test, um framework de testes para pipelines e módulos Nextflow. Para aprender como escrever e executar esses testes, veja a [missão secundária nf-test](../side_quests/nf-test.md).

O `main.nf` gerado inclui todos os padrões que você acabou de aprender, mais alguns recursos adicionais:

```groovy title="modules/local/cowpy/main.nf" hl_lines="11 21 22"
process COWPY {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE':
        'biocontainers/YOUR-TOOL-HERE' }"

    input:
    tuple val(meta), path(input)        // Padrão 1: Tuplas de metadados ✓

    output:
    tuple val(meta), path("*"), emit: output
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''                // Padrão 2: ext.args ✓
    def prefix = task.ext.prefix ?: "${meta.id}"  // Padrão 3: ext.prefix ✓

    """
    // Add your tool command here

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        COWPY: \$(cowpy --version)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    echo $args
    touch ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        COWPY: \$(cowpy --version)
    END_VERSIONS
    """
}
```

Note como todos os padrões que você aplicou manualmente acima já estão lá!

O template também inclui várias convenções nf-core adicionais.
Algumas delas funcionam imediatamente, enquanto outras são placeholders que precisaremos preencher, conforme descrito abaixo.

**Recursos que funcionam como estão:**

- **`tag "$meta.id"`**: Adiciona ID de amostra aos nomes de processo nos logs para rastreamento mais fácil
- **`label 'process_single'`**: Label de recurso para configurar requisitos de CPU/memória
- **Bloco `when:`**: Permite execução condicional via configuração `task.ext.when`

Esses recursos já são funcionais e tornam os módulos mais fáceis de manter.

**Placeholders que personalizaremos abaixo:**

- **Blocos `input:` e `output:`**: Declarações genéricas que atualizaremos para corresponder à nossa ferramenta
- **Bloco `script:`**: Contém um comentário onde adicionaremos o comando `cowpy`
- **Bloco `stub:`**: Template que atualizaremos para produzir as saídas corretas
- **Container e ambiente**: Placeholders que preencheremos com informações de pacote

As próximas seções percorrem a conclusão dessas personalizações.

### 2.2. Configurar o contêiner e o ambiente conda

As diretrizes nf-core exigem que especifiquemos tanto um contêiner quanto um ambiente Conda como parte do módulo.

#### 2.2.1. Contêiner

Para o contêiner, você pode usar [Seqera Containers](https://seqera.io/containers/) para construir automaticamente um contêiner a partir de qualquer pacote Conda, incluindo pacotes conda-forge.
Neste caso estamos usando o mesmo contêiner pré-construído de antes.

O código padrão oferece alternar entre Docker e Singularity, mas vamos simplificar essa linha e apenas especificar o contêiner Docker que obtivemos do Seqera Containers acima.

=== "Depois"

```groovy title="modules/local/cowpy/main.nf" linenums="3" hl_lines="6"
process COWPY {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273"
```

=== "Antes"

```groovy title="modules/local/cowpy/main.nf" linenums="3" hl_lines="6"
process COWPY {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE':
        'biocontainers/YOUR-TOOL-HERE' }"
```

#### 2.2.2. Ambiente Conda

Para o ambiente Conda, o código do módulo especifica `conda "${moduleDir}/environment.yml"` o que significa que deve ser configurado no arquivo `environment.yml`.

A ferramenta de criação de módulo nos avisou que não conseguiu encontrar o pacote `cowpy` no Bioconda (o canal principal para ferramentas de bioinformática).
No entanto, `cowpy` está disponível no conda-forge, então você pode completar o `environment.yml` assim:

=== "Depois"

    ```yaml title="modules/local/cowpy/environment.yml"  linenums="1" hl_lines="1 3 5"
    name: COWPY
    channels:
      - conda-forge
    dependencies:
      - cowpy=1.1.5
    ```

=== "Antes"

    ```yaml title="modules/local/cowpy/environment.yml" linenums="1"
    ---
    # yaml-language-server: $schema=https://raw.githubusercontent.com/nf-core/modules/master/modules/environment-schema.json
    channels:
      - conda-forge
      - bioconda
    dependencies:
      # TODO nf-core: List required Conda package(s).
      #               Software MUST be pinned to channel (i.e. "bioconda"), version (i.e. "1.10").
      #               For Conda, the build (i.e. "h9402c20_2") must be EXCLUDED to support installation on different operating systems.
      - "YOUR-TOOL-HERE"
    ```

Para submissão ao nf-core, teríamos que seguir os padrões mais de perto, mas para nosso próprio uso podemos simplificar o código desta forma.

!!! tip "Pacotes Bioconda vs conda-forge"

    - **Pacotes Bioconda**: Automaticamente obtêm BioContainers construídos, fornecendo contêineres prontos para uso
    - **Pacotes conda-forge**: Podem usar Seqera Containers para construir contêineres sob demanda a partir da receita Conda

    A maioria das ferramentas de bioinformática está no Bioconda, mas para ferramentas conda-forge, Seqera Containers fornece uma solução fácil para containerização.

### 2.3. Conectar a lógica `COWPY`

Agora vamos atualizar os elementos de código que são específicos do que o processo `COWPY` faz: as entradas e saídas, e o bloco script.

#### 2.3.1. Entradas e saídas

O template gerado inclui declarações genéricas de entrada e saída que você precisará personalizar para sua ferramenta específica.
Olhando de volta para nosso módulo `COWPY` manual da seção 1, podemos usar isso como guia.

Atualize os blocos de entrada e saída:

=== "Depois"

    ```groovy title="modules/local/cowpy/main.nf" linenums="8" hl_lines="2 5"
    input:
    tuple val(meta), path(input_file)

    output:
    tuple val(meta), path("${prefix}.txt"), emit: cowpy_output
    path "versions.yml"           , emit: versions
    ```

=== "Antes"

    ```groovy title="modules/local/cowpy/main.nf" linenums="8" hl_lines="2 5"
    input:
    tuple val(meta), path(input)

    output:
    tuple val(meta), path("*"), emit: output
    path "versions.yml"           , emit: versions
    ```

Isso especifica:

- O nome do parâmetro do arquivo de entrada (`input_file` em vez de `input` genérico)
- O nome do arquivo de saída usando o padrão de prefixo configurável (`${prefix}.txt` em vez do curinga `*`)
- Um nome de emit descritivo (`cowpy_output` em vez de `output` genérico)

Se você estiver usando o servidor de linguagem Nextflow para validar sintaxe, a parte `${prefix}` será sinalizada como um erro neste estágio porque ainda não a adicionamos ao bloco script.
Vamos fazer isso agora.

#### 2.3.2. O bloco script

O template fornece um placeholder de comentário no bloco script onde você deve adicionar o comando real da ferramenta.

Baseado no módulo que escrevemos manualmente anteriormente, devemos fazer as seguintes edições:

=== "Depois"

    ```groovy title="modules/local/cowpy/main.nf" linenums="15" hl_lines="3 6"
    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    cat $input_file | cowpy $args > ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        COWPY: \$(cowpy --version)
    END_VERSIONS
    """
    ```

=== "Antes"

    ```groovy title="modules/local/cowpy/main.nf" linenums="15" hl_lines="6"
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    // Add your tool command here

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        COWPY: \$(cowpy --version)
    END_VERSIONS
    """
    ```

Mudanças principais:

- Mudar `def prefix` para apenas `prefix` (sem `def`) para torná-lo acessível no bloco output
- Substituir o comentário pelo comando `cowpy` real que usa tanto `$args` quanto `${prefix}.txt`

Note que se não tivéssemos já feito o trabalho de adicionar a configuração `ext.args` e `ext.prefix` para o processo `COWPY` ao arquivo `modules.config`, precisaríamos fazer isso agora.

#### 2.3.3. Implementar o bloco stub

No contexto Nextflow, um bloco [stub](https://www.nextflow.io/docs/latest/process.html#stub) permite que você defina um script leve e fictício usado para prototipagem rápida e teste da lógica de um pipeline sem executar o comando real.

Não se preocupe muito se isso parecer misterioso; incluímos isso para completude, mas você também pode simplesmente deletar a seção stub se não quiser lidar com isso, pois é completamente opcional.

=== "Depois"

    ```groovy title="modules/local/cowpy/main.nf" linenums="27" hl_lines="3 6"
    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        COWPY: \$(cowpy --version)
    END_VERSIONS
    """
    ```

=== "Antes"

    ```groovy title="modules/local/cowpy/main.nf" linenums="27" hl_lines="3 6"
    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    echo $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        COWPY: \$(cowpy --version)
    END_VERSIONS
    """
    ```

Mudanças principais:

- Mudar `def prefix` para apenas `prefix` para corresponder ao bloco script
- Remover a linha `echo $args` (que era apenas código placeholder do template)
- O stub cria um arquivo `${prefix}.txt` vazio correspondendo ao que o bloco script produz

Isso permite que você teste a lógica do fluxo de trabalho e manipulação de arquivos sem esperar a ferramenta real executar.

Uma vez que você tenha completado a configuração do ambiente (seção 2.2), entradas/saídas (seção 2.3.1), bloco script (seção 2.3.2) e bloco stub (seção 2.3.3), o módulo está pronto para testar!

### 2.4. Trocar o novo módulo `COWPY` e executar o pipeline

Tudo que precisamos fazer para experimentar esta nova versão do módulo `COWPY` é mudar a declaração de importação no arquivo de fluxo de trabalho `hello.nf` para apontar para o novo arquivo.

=== "Depois"

    ```groovy title="workflows/hello.nf" linenums="1" hl_lines="10"
    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    include { paramsSummaryMap       } from 'plugin/nf-schema'
    include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
    include { sayHello               } from '../modules/local/sayHello.nf'
    include { convertToUpper         } from '../modules/local/convertToUpper.nf'
    include { COWPY                  } from '../modules/local/cowpy/main.nf'
    include { CAT_CAT                } from '../modules/nf-core/cat/cat/main'
    ```

=== "Antes"

    ```groovy title="modules/local/cowpy/main.nf" linenums="1" hl_lines="10"
    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    include { paramsSummaryMap       } from 'plugin/nf-schema'
    include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
    include { sayHello               } from '../modules/local/sayHello.nf'
    include { convertToUpper         } from '../modules/local/convertToUpper.nf'
    include { COWPY                  } from '../modules/local/cowpy.nf'
    include { CAT_CAT                } from '../modules/nf-core/cat/cat/main'
    ```

Vamos executar o pipeline para testá-lo.

```bash
nextflow run . --outdir core-hello-results -profile test,docker --validate_params false
```

??? success "Saída do comando"

    ```console hl_lines="33"
      N E X T F L O W   ~  version 25.04.3

    Launching `./main.nf` [prickly_neumann] DSL2 - revision: b9e9b3b8de

    Input/output options
      input                     : /workspaces/training/hello-nf-core/core-hello/assets/greetings.csv
      outdir                    : core-hello-results

    Institutional config options
      config_profile_name       : Test profile
      config_profile_description: Minimal test dataset to check pipeline function

    Generic options
      validate_params           : false
      trace_report_suffix       : 2025-12-27_08-23-51

    Core Nextflow options
      runName                   : prickly_neumann
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
    [e9/008ede] CORE_HELLO:HELLO:sayHello (3)       [100%] 3 of 3 ✔
    [f0/d70cfe] CORE_HELLO:HELLO:convertToUpper (3) [100%] 3 of 3 ✔
    [be/0ecc58] CORE_HELLO:HELLO:CAT_CAT (test)     [100%] 1 of 1 ✔
    [11/8e082f] CORE_HELLO:HELLO:COWPY (test)       [100%] 1 of 1 ✔
    -[core/hello] Pipeline completed successfully-
    ```

Isso produz os mesmos resultados de antes.

### Conclusão

Você agora sabe como usar as ferramentas integradas nf-core para criar módulos eficientemente usando templates em vez de escrever tudo do zero.

### O que vem a seguir?

Aprenda quais são os benefícios de contribuir módulos para o nf-core e quais são as principais etapas e requisitos envolvidos.

---

## 3. Contribuir módulos de volta para o nf-core

O repositório [nf-core/modules](https://github.com/nf-core/modules) recebe contribuições de módulos bem testados e padronizados.

### 3.1. Por que contribuir?

Contribuir seus módulos para o nf-core:

- Torna suas ferramentas disponíveis para toda a comunidade nf-core através do catálogo de módulos em [nf-co.re/modules](https://nf-co.re/modules)
- Garante manutenção e melhorias contínuas da comunidade
- Fornece garantia de qualidade através de revisão de código e testes automatizados
- Dá ao seu trabalho visibilidade e reconhecimento

### 3.2. Lista de verificação do contribuidor

Para contribuir um módulo para o nf-core, você precisará passar pelas seguintes etapas:

1. Verificar se já existe em [nf-co.re/modules](https://nf-co.re/modules)
2. Fazer fork do repositório [nf-core/modules](https://github.com/nf-core/modules)
3. Usar `nf-core modules create` para gerar o template
4. Preencher a lógica e testes do módulo
5. Testar com `nf-core modules test tool/subtool`
6. Fazer lint com `nf-core modules lint tool/subtool`
7. Submeter um pull request

Para instruções detalhadas, veja o [tutorial de componentes nf-core](https://nf-co.re/docs/tutorials/nf-core_components/components).

### 3.3. Recursos

- **Tutorial de componentes**: [Guia completo para criar e contribuir módulos](https://nf-co.re/docs/tutorials/nf-core_components/components)
- **Especificações de módulos**: [Requisitos técnicos e diretrizes](https://nf-co.re/docs/guidelines/components/modules)
- **Suporte da comunidade**: [Slack nf-core](https://nf-co.re/join) - Junte-se ao canal `#modules`

### Conclusão

Você agora sabe como criar módulos nf-core! Você aprendeu os quatro padrões principais que tornam os módulos portáveis e de fácil manutenção:

- **Tuplas de metadados** propagam metadados através do fluxo de trabalho
- **`ext.args`** simplifica interfaces de módulos manipulando argumentos opcionais via configuração
- **`ext.prefix`** padroniza a nomenclatura de arquivos de saída
- **Publicação centralizada** via `publishDir` configurado em `modules.config` em vez de codificado em módulos

Ao transformar `COWPY` passo a passo, você desenvolveu uma compreensão profunda desses padrões, tornando-o equipado para trabalhar com, depurar e criar módulos nf-core.
Na prática, você usará `nf-core modules create` para gerar módulos estruturados adequadamente com esses padrões incorporados desde o início.

Finalmente, você aprendeu como contribuir módulos para a comunidade nf-core, tornando ferramentas disponíveis para pesquisadores em todo o mundo enquanto se beneficia da manutenção contínua da comunidade.

### O que vem a seguir?

Quando estiver pronto, continue para [Parte 5: Validação de entrada](./05_input_validation.md) para aprender como adicionar validação de entrada baseada em schema ao seu pipeline.
