# Parte 4: Criar um módulo nf-core

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tradução assistida por IA - [saiba mais e sugira melhorias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Nesta quarta parte do curso de treinamento Hello nf-core, mostramos como criar um módulo nf-core aplicando as convenções principais que tornam os módulos portáteis e de fácil manutenção.

O projeto nf-core fornece um comando (`nf-core modules create`) que gera templates de módulos estruturados corretamente de forma automática, semelhante ao que usamos para o fluxo de trabalho na Parte 2.
No entanto, para fins didáticos, vamos começar fazendo manualmente: transformar o módulo local `cowpy` em seu pipeline `core-hello` em um módulo no estilo nf-core passo a passo.
Depois disso, mostraremos como usar a criação de módulos baseada em template para trabalhar de forma mais eficiente no futuro.

??? info "Como começar desta seção"

    Esta seção pressupõe que você completou a [Parte 3: Usar um módulo nf-core](./03_use_module.md) e integrou o módulo `CAT_CAT` ao seu pipeline.

    Se você não completou a Parte 3 ou deseja começar do zero para esta parte, pode usar a solução `core-hello-part3` como ponto de partida.
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

// Gera arte ASCII com cowpy (https://github.com/jeffbuttars/cowpy)
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

Aplicaremos as seguintes convenções nf-core incrementalmente:

1. **Colocar o nome do processo em maiúsculas para `COWPY`** para seguir a convenção.
2. **Atualizar `COWPY` para usar tuplas de metadata** para propagar metadata da amostra através do fluxo de trabalho.
3. **Centralizar a configuração de argumentos da ferramenta com `ext.args`** para aumentar a versatilidade do módulo mantendo a interface mínima.
4. **Padronizar a nomeação de saída com `ext.prefix`** para promover consistência.
5. **Centralizar a configuração de publicação** para promover consistência.

Após cada passo, executaremos o pipeline para testar que tudo funciona como esperado.

!!! warning "Diretório de trabalho"

    Certifique-se de estar no diretório `core-hello` (a raiz do seu pipeline) para todas as edições de arquivo e execuções de comando nesta seção.

    ```bash
    cd core-hello
    ```

### 1.1. Colocar o nome do processo em maiúsculas

Esta é puramente uma convenção estilística (não há justificativa técnica), mas como é a norma para módulos nf-core, vamos seguir.

Precisamos fazer três conjuntos de alterações:

1. Atualizar o nome do processo no módulo
2. Atualizar a declaração de importação do módulo no cabeçalho do fluxo de trabalho
3. Atualizar a chamada do processo e a declaração de emit no corpo do fluxo de trabalho

Vamos começar!

#### 1.1.1. Atualizar o nome do processo no módulo

Abra o arquivo do módulo `cowpy.nf` (em `core-hello/modules/local/`) e modifique o nome do processo para maiúsculas:

=== "Depois"

    ```groovy title="core-hello/modules/local/cowpy.nf" linenums="3" hl_lines="2"
    // Gera arte ASCII com cowpy (https://github.com/jeffbuttars/cowpy)
    process COWPY {
    ```

=== "Antes"

    ```groovy title="core-hello/modules/local/cowpy.nf" linenums="3" hl_lines="2"
    // Gera arte ASCII com cowpy (https://github.com/jeffbuttars/cowpy)
    process cowpy {
    ```

Neste caso, a conversão para maiúsculas é completamente direta.

Se o nome do processo fosse composto por várias palavras, por exemplo se tivéssemos um processo chamado MyCowpyTool originalmente em camel case, a convenção nf-core seria usar underscores para separá-las, resultando em MY_COWPY_TOOL.

#### 1.1.2. Atualizar a declaração de importação do módulo

Os nomes dos processos diferenciam maiúsculas de minúsculas, então agora que mudamos o nome do processo, precisamos atualizar a declaração de importação do módulo adequadamente no cabeçalho do fluxo de trabalho de `hello.nf`:

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

Poderíamos usar um alias na declaração de importação para evitar ter que atualizar as chamadas ao processo, mas isso de certa forma frustraria o objetivo de adotar a convenção de maiúsculas.

#### 1.1.3. Atualizar a chamada do processo e a declaração de emit

Então agora vamos atualizar as duas referências ao processo no bloco workflow de `hello.nf`:

=== "Depois"

    ```groovy title="core-hello/workflows/hello.nf" linenums="43" hl_lines="2 17"
    // gera arte ASCII das saudações com cowpy
    COWPY(CAT_CAT.out.file_out)

    //
    // Agrupar e salvar versões de software
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
    // gera arte ASCII das saudações com cowpy
    cowpy(CAT_CAT.out.file_out)

    //
    // Agrupar e salvar versões de software
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

Certifique-se de fazer **ambas** as alterações, caso contrário você obterá um erro ao executar isto.

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

### 1.2. Atualizar `COWPY` para usar tuplas de metadata

Na versão atual do pipeline `core-hello`, estamos extraindo o arquivo da tupla de saída do `CAT_CAT` para passar ao `COWPY`, como mostrado na metade superior do diagrama abaixo.

<figure class="excalidraw">
    --8<-- "docs/en/docs/hello_nf-core/img/cowpy-inputs.svg"
</figure>

Seria melhor ter o `COWPY` aceitando tuplas de metadata diretamente, permitindo que os metadata fluam através do fluxo de trabalho, como mostrado na metade inferior do diagrama.

Para isso, precisaremos fazer as seguintes alterações:

1. Atualizar as definições de entrada e saída
2. Atualizar a chamada do processo no fluxo de trabalho
3. Atualizar o bloco emit no fluxo de trabalho

Uma vez que tivermos feito tudo isso, executaremos o pipeline para testar que tudo ainda funciona como antes.

#### 1.2.1. Atualizar as definições de entrada e saída

Retorne ao arquivo do módulo `cowpy.nf` e modifique-o para aceitar tuplas de metadata como mostrado abaixo.

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
Agora que a saída do `CAT_CAT` e a entrada do `COWPY` têm a mesma 'forma', ou seja, ambas consistem em uma estrutura `tuple val(meta), path(input_file)`, podemos simplesmente conectá-las diretamente em vez de ter que extrair o arquivo explicitamente da saída do processo `CAT_CAT`.

Abra o arquivo de fluxo de trabalho `hello.nf` (em `core-hello/workflows/`) e atualize a chamada para `COWPY` como mostrado abaixo.

=== "Depois"

    ```groovy title="core-hello/workflows/hello.nf" linenums="43" hl_lines="2"
        // gera arte ASCII das saudações com cowpy
        COWPY(CAT_CAT.out.file_out, params.character)
    ```

=== "Antes"

    ```groovy title="core-hello/workflows/hello.nf" linenums="43" hl_lines="1-2 5"
        // extract the file from the tuple since cowpy doesn't use metadata yet
        ch_for_cowpy = CAT_CAT.out.file_out.map{ meta, file -> file }

        // gera arte ASCII das saudações com cowpy
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

Isso não é tecnicamente necessário, mas é uma boa prática referir-se a saídas nomeadas sempre que possível.

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

O pipeline deve executar com sucesso, com os metadata agora fluindo do `CAT_CAT` através do `COWPY`.

Isso completa o que precisávamos fazer para que o `COWPY` manipule tuplas de metadata.
Agora, vamos ver o que mais podemos fazer para aproveitar os padrões de módulos nf-core.

### 1.3. Centralizar a configuração de argumentos da ferramenta com `ext.args`

Em seu estado atual, o processo `COWPY` espera receber um valor para o parâmetro `character`.
Como resultado, temos que fornecer um valor toda vez que chamamos o processo, mesmo se estivéssemos satisfeitos com os padrões definidos pela ferramenta.
Para o `COWPY` isso não é admitidamente um grande problema, mas para ferramentas com muitos parâmetros opcionais, pode se tornar bastante trabalhoso.

O projeto nf-core recomenda usar um recurso do Nextflow chamado [`ext.args`](https://www.nextflow.io/docs/latest/reference/process.html#ext) para gerenciar argumentos de ferramentas de forma mais conveniente via arquivos de configuração.

Em vez de declarar entradas de processo para cada opção da ferramenta, você escreve o módulo para referenciar `ext.args` na construção de sua linha de comando.
Então é só uma questão de configurar a variável `ext.args` para conter os argumentos e valores que você deseja usar no arquivo `modules.config`, que consolida detalhes de configuração para todos os módulos.
O Nextflow adicionará esses argumentos com seus valores à linha de comando da ferramenta em tempo de execução.

Vamos aplicar essa abordagem ao módulo `COWPY`.
Vamos precisar fazer as seguintes alterações:

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

    // Gera arte ASCII com cowpy (https://github.com/jeffbuttars/cowpy)
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

    // Gera arte ASCII com cowpy (https://github.com/jeffbuttars/cowpy)
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
   Daqui para frente, forneceremos esse argumento via a configuração `ext.args` conforme descrito abaixo.

2. **No bloco `script:`, adicionamos a linha `def args = task.ext.args ?: ''`.**
   Essa linha usa o operador `?:` para determinar o valor da variável `args`: o conteúdo de `task.ext.args` se não estiver vazio, ou uma string vazia se estiver.
   Note que embora geralmente nos referimos a `ext.args`, este código deve referenciar `task.ext.args` para extrair a configuração `ext.args` no nível do módulo.

3. **Na linha de comando, substituímos `-c "$character"` por `$args`.**
   É aqui que o Nextflow injetará quaisquer argumentos de ferramenta definidos em `ext.args` no arquivo `modules.config`.

Como resultado, a interface do módulo agora é mais simples: ela só espera as entradas essenciais de metadata e arquivo.

!!! note

    O operador `?:` é frequentemente chamado de 'operador Elvis' porque se parece com um rosto de Elvis Presley de lado, com o caractere `?` simbolizando a onda em seu cabelo.

#### 1.3.2. Configurar `ext.args` no arquivo `modules.config`

Agora que tiramos a declaração de `character` do módulo, temos que adicioná-la ao `ext.args` no arquivo de configuração `modules.config`.

Especificamente, vamos adicionar este pequeno pedaço de código ao bloco `process {}`:

```groovy title="Código a adicionar"
withName: 'COWPY' {
    ext.args = { "-c ${params.character}" }
}
```

A sintaxe `withName:` atribui essa configuração apenas ao processo `COWPY`, e `ext.args = { "-c ${params.character}" }` simplesmente compõe uma string que incluirá o valor do parâmetro `character`.
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

- A **interface do módulo permanece simples** - Ela só aceita as entradas essenciais de metadata e arquivo
- O **pipeline ainda expõe `params.character`** - Os usuários finais ainda podem configurá-lo como antes
- O **módulo agora é portátil** - Ele pode ser reutilizado em outros pipelines sem esperar um nome de parâmetro específico
- A configuração é **centralizada** em `modules.config`, mantendo a lógica do fluxo de trabalho limpa

Ao usar o arquivo `modules.config` como o lugar onde todos os pipelines centralizam a configuração por módulo, tornamos nossos módulos mais reutilizáveis em diferentes pipelines.

#### 1.3.3. Atualizar o fluxo de trabalho `hello.nf`

Como o módulo `COWPY` não requer mais o parâmetro `character` como entrada, precisamos atualizar a chamada do fluxo de trabalho adequadamente.

Abra o arquivo de fluxo de trabalho `hello.nf` (em `core-hello/workflows/`) e atualize a chamada para `COWPY` como mostrado abaixo.

=== "Depois"

    ```groovy title="core-hello/workflows/hello.nf" linenums="39" hl_lines="2"
        // gera arte ASCII das saudações com cowpy
        COWPY(CAT_CAT.out.file_out)
    ```

=== "Antes"

    ```groovy title="core-hello/workflows/hello.nf" linenums="39" hl_lines="2"
        // gera arte ASCII das saudações com cowpy
        COWPY(CAT_CAT.out.file_out, params.character)
    ```

O código do fluxo de trabalho agora está mais limpo: não precisamos passar `params.character` diretamente para o processo.
A interface do módulo é mantida mínima, tornando-a mais portátil, enquanto o pipeline ainda fornece a opção explícita através da configuração.

#### 1.3.4. Executar o pipeline para testá-lo

Vamos testar que o fluxo de trabalho ainda funciona como esperado, especificando um personagem diferente para verificar que a configuração `ext.args` está funcionando.

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
Encontre a saída no navegador de arquivos ou use o hash da tarefa (a parte `38/eb29ea` no exemplo acima) para visualizar o arquivo de saída:

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

Você deve ver a arte ASCII exibida com o personagem `kosh`, confirmando que a configuração `ext.args` funcionou!

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
Esta abordagem mantém a interface do módulo focada em dados essenciais (arquivos, metadata e quaisquer parâmetros obrigatórios por amostra), enquanto opções que controlam o comportamento da ferramenta são tratadas separadamente através da configuração.

Isso pode parecer desnecessário para uma ferramenta simples como `cowpy`, mas pode fazer uma grande diferença para ferramentas de análise de dados que têm muitos argumentos opcionais.

Para resumir os benefícios desta abordagem:

- **Interface limpa**: O módulo foca em entradas de dados essenciais (metadata e arquivos)
- **Flexibilidade**: Os usuários podem especificar argumentos de ferramentas via configuração, incluindo valores específicos por amostra
- **Consistência**: Todos os módulos nf-core seguem este padrão
- **Portabilidade**: Módulos podem ser reutilizados sem opções de ferramenta codificadas
- **Sem mudanças no fluxo de trabalho**: Adicionar ou alterar opções de ferramenta não requer atualizar o código do fluxo de trabalho

!!! note

    O sistema `ext.args` tem capacidades adicionais poderosas não cobertas aqui, incluindo alternar valores de argumentos dinamicamente com base em metadata. Veja as [especificações de módulos nf-core](https://nf-co.re/docs/guidelines/components/modules) para mais detalhes.

### 1.4. Padronizar a nomeação de saída com `ext.prefix`

Agora que demos ao processo `COWPY` acesso ao metamap, podemos começar a aproveitar outro padrão útil do nf-core: nomear arquivos de saída com base em metadata.

Aqui vamos usar um recurso do Nextflow chamado `ext.prefix` que nos permitirá padronizar a nomeação de arquivos de saída em todos os módulos usando `meta.id` (o identificador incluído no metamap), ainda sendo capaz de configurar módulos individualmente se desejado.

Isso será semelhante ao que fizemos com `ext.args`, com algumas diferenças que detalharemos conforme avançamos.

Vamos aplicar essa abordagem ao módulo `COWPY`.
Vamos precisar fazer as seguintes alterações:

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
   Note que embora geralmente nos referimos a `ext.prefix`, este código deve referenciar `task.ext.prefix` para extrair a configuração `ext.prefix` no nível do módulo.

2. **Na linha de comando, substituímos `cowpy-${input_file}` por `${prefix}.txt`.**
   É aqui que o Nextflow injetará o valor de `prefix` determinado pela linha acima.

3. **No bloco `output:`, substituímos `path("cowpy-${input_file}")` por `path("${prefix}.txt")`.**
   Isso simplesmente reitera qual será o caminho do arquivo de acordo com o que está escrito na linha de comando.

Como resultado, o nome do arquivo de saída agora é construído usando um padrão sensato (o identificador do metamap) combinado com a extensão de formato de arquivo apropriada.

#### 1.4.2. Configurar `ext.prefix` no arquivo `modules.config`

Neste caso, o padrão sensato não é suficientemente expressivo para nosso gosto; queremos usar um padrão de nomeação personalizado que inclua o nome da ferramenta, `cowpy-<id>.txt`, como tínhamos antes.

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

Caso você esteja se perguntando, a closure `ext.prefix` tem acesso ao pedaço correto de metadata porque a configuração é avaliada no contexto da execução do processo, onde os metadata estão disponíveis.

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
Você deve ver o arquivo de saída do cowpy com a mesma nomeação de antes: `cowpy-test.txt`, baseado no nome de lote padrão.

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

Sinta-se à vontade para alterar a configuração `ext.prefix` em `conf/modules.config` para satisfazer-se de que você pode alterar o padrão de nomeação sem ter que fazer nenhuma mudança no código do módulo ou do fluxo de trabalho.

Alternativamente, você também pode tentar executar isso novamente com um parâmetro `--batch` diferente especificado na linha de comando para satisfazer-se de que essa parte ainda é personalizável em tempo real.

Isso demonstra como `ext.prefix` permite que você mantenha sua convenção de nomeação preferida enquanto mantém a interface do módulo flexível.

Para resumir os benefícios desta abordagem:

- **Nomeação padronizada**: Os arquivos de saída são tipicamente nomeados usando IDs de amostra dos metadata
- **Configurável**: Os usuários podem sobrescrever a nomeação padrão se necessário
- **Consistente**: Todos os módulos nf-core seguem este padrão
- **Previsível**: Fácil saber como os arquivos de saída serão chamados

Muito bom, certo?
Bem, há mais uma mudança importante que precisamos fazer para melhorar nosso módulo para se ajustar às diretrizes nf-core.

### 1.5. Centralizar a configuração de publicação

Você pode ter notado que temos publicado saídas em dois diretórios diferentes:

- **`results`** — O diretório de saída original que temos usado desde o início para nossos módulos locais, definido individualmente usando diretivas `publishDir` por módulo;
- **`core-hello-results`** — O diretório de saída definido com `--outdir` na linha de comando, que tem recebido os logs nf-core e os resultados publicados pelo `CAT_CAT`.

Isso é confuso e subótimo; seria melhor ter um local para tudo.
Claro, poderíamos entrar em cada um de nossos módulos locais e atualizar a diretiva `publishDir` manualmente para usar o diretório `core-hello-results`, mas e na próxima vez que decidirmos mudar o diretório de saída?

Ter módulos individuais tomando decisões de publicação claramente não é o caminho a seguir, especialmente em um mundo onde o mesmo módulo pode ser usado em muitos pipelines diferentes, por pessoas que têm necessidades ou preferências diferentes.
Queremos ser capazes de controlar onde as saídas são publicadas no nível da configuração do fluxo de trabalho.

"Ei," você pode dizer, "`CAT_CAT` está enviando suas saídas para o `--outdir`. Talvez devêssemos copiar sua diretiva `publishDir`?"

Sim, essa é uma ótima ideia.

Exceto que ele não tem uma diretiva `publishDir`. (Vá em frente, olhe o código do módulo.)

Isso porque os pipelines nf-core centralizam o controle no nível do fluxo de trabalho configurando `publishDir` em `conf/modules.config` em vez de em módulos individuais.
Especificamente, o template nf-core declara uma diretiva `publishDir` padrão (com uma estrutura de diretórios predefinida) que se aplica a todos os módulos, a menos que uma diretiva sobrescrevendo seja fornecida.

Isso não parece incrível? Seria possível que para aproveitar essa diretiva padrão, tudo o que precisamos fazer é remover a diretiva `publishDir` atual de nossos módulos locais?

Vamos tentar isso no `COWPY` para ver o que acontece, então veremos o código para a configuração padrão para entender como funciona.

Finalmente, demonstraremos como sobrescrever o comportamento padrão se desejado.

#### 1.5.1. Remover a diretiva `publishDir` do `COWPY`

Vamos fazer isso.
Abra o arquivo do módulo `cowpy.nf` (em `core-hello/modules/local/`) e remova a diretiva `publishDir` como mostrado abaixo.

=== "Depois"

    ```groovy title="core-hello/modules/local/cowpy.nf (trecho)" linenums="1"
    #!/usr/bin/env nextflow

    // Gera arte ASCII com cowpy (https://github.com/jeffbuttars/cowpy)
    process COWPY {

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
        conda 'conda-forge::cowpy==1.1.5'
    ```

=== "Antes"

    ```groovy title="core-hello/modules/local/cowpy.nf (trecho)" linenums="1" hl_lines="6"
    #!/usr/bin/env nextflow

    // Gera arte ASCII com cowpy (https://github.com/jeffbuttars/cowpy)
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

Você pode ver que o Nextflow criou essa hierarquia de diretórios baseada nos nomes do fluxo de trabalho e do módulo.

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

- **`path:`** Determina o diretório de saída com base no nome do processo.
  O nome completo de um processo contido em `task.process` inclui a hierarquia de importações de fluxo de trabalho e módulo (como `CORE_HELLO:HELLO:CAT_CAT`).
  As operações `tokenize` removem essa hierarquia para obter apenas o nome do processo, então pegam a primeira parte antes de qualquer underscore (se aplicável), e convertem para minúsculas.
  Isso é o que determina que os resultados do `CAT_CAT` sejam publicados em `${params.outdir}/cat/`.
- **`mode:`** Controla como os arquivos são publicados (copy, symlink, etc.).
  Isso é configurável via o parâmetro `params.publish_dir_mode`.
- **`saveAs:`** Filtra quais arquivos publicar.
  Este exemplo exclui arquivos `versions.yml` retornando `null` para eles, impedindo que sejam publicados.

Isso fornece uma lógica consistente para organizar saídas.

A saída fica ainda melhor quando todos os módulos em um pipeline adotam essa convenção, então sinta-se à vontade para ir deletar as diretivas `publishDir` dos outros módulos em seu pipeline.
Esse padrão será aplicado mesmo a módulos que não modificamos explicitamente para seguir as diretrizes nf-core.

Dito isso, você pode decidir que deseja organizar suas entradas de forma diferente, e a boa notícia é que é fácil fazer isso.

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

Não vamos realmente fazer essa mudança, mas sinta-se à vontade para experimentar isso e ver que lógica você pode implementar.

O ponto é que esse sistema permite que você tenha o melhor de ambos os mundos: consistência por padrão e a flexibilidade para personalizar a configuração sob demanda.

Para resumir, você obtém:

- **Fonte única de verdade**: Toda configuração de publicação vive em `modules.config`
- **Padrão útil**: Processos funcionam prontos para uso sem configuração por módulo
- **Personalização fácil**: Sobrescreva o comportamento de publicação na configuração, não no código do módulo
- **Módulos portáteis**: Módulos não codificam locais de saída

Isso completa o conjunto de recursos de módulos nf-core que você deve absolutamente aprender a usar, mas há outros que você pode ler sobre nas [especificações de módulos nf-core](https://nf-co.re/docs/guidelines/components/modules).

### Conclusão

Você agora sabe como adaptar módulos locais para seguir as convenções nf-core:

- Projete seus módulos para aceitar e propagar tuplas de metadata;
- Use `ext.args` para manter as interfaces dos módulos mínimas e portáteis;
- Use `ext.prefix` para nomeação de arquivos de saída configurável e padronizada;
- Adote a diretiva `publishDir` centralizada padrão para uma estrutura consistente de diretório de resultados.

### O que vem a seguir?

Aprenda como usar as ferramentas integradas baseadas em template do nf-core para criar módulos de forma mais fácil.

---

## 2. Criar um módulo com as ferramentas nf-core

Agora que você aprendeu os padrões de módulos nf-core aplicando-os manualmente, vamos ver como você criaria módulos na prática.

### 2.1. Gerar um scaffold de módulo a partir de um template

Semelhante ao que existe para criar pipelines, o projeto nf-core fornece ferramentas para gerar módulos adequadamente estruturados baseados em um template, com todos esses padrões incorporados desde o início.

#### 2.1.1. Executar o comando de criação de módulo

O comando `nf-core modules create` gera um template de módulo que já segue todas as convenções que você aprendeu.

Vamos criar uma nova versão do módulo `COWPY` com um template mínimo executando este comando:

```bash
nf-core modules create --empty-template COWPY
```

A flag `--empty-template` cria um template inicial limpo sem código extra, facilitando ver a estrutura essencial.

O comando é executado interativamente, guiando você pela configuração.
Ele automaticamente busca informações da ferramenta de repositórios de pacotes como Bioconda e bio.tools para pré-popular metadata.

Você será solicitado para várias opções de configuração:

- **Informações do autor**: Seu nome de usuário do GitHub para atribuição
- **Label de recurso**: Um conjunto predefinido de requisitos computacionais.
  O projeto nf-core fornece labels padrão como `process_single` para ferramentas leves e `process_high` para as exigentes.
  Essas labels ajudam a gerenciar alocação de recursos em diferentes ambientes de execução.
- **Requisito de metadata**: Se o módulo precisa de informações específicas da amostra via um map `meta` (geralmente sim para módulos de processamento de dados).

A ferramenta lida com a complexidade de encontrar informações de pacotes e configurar a estrutura, permitindo que você se concentre em implementar a lógica específica da ferramenta.

#### 2.1.2. Examinar o scaffold do módulo

A ferramenta cria uma estrutura completa de módulo em `modules/local/` (ou `modules/nf-core/` se você estiver no repositório nf-core/modules):

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
- **`environment.yml`**: Especificação de ambiente Conda para dependências
- **`tests/main.nf.test`**: Casos de teste nf-test para validar que o módulo funciona

!!! tip "Saiba mais sobre testes"

    O arquivo de teste gerado usa nf-test, um framework de testes para pipelines e módulos Nextflow. Para aprender como escrever e executar esses testes, veja a [missão paralela nf-test](../side_quests/nf-test.md).

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
    tuple val(meta), path(input)        // Padrão 1: Tuplas de metadata ✓

    output:
    tuple val(meta), path("*"), emit: output
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''                // Padrão 2: ext.args ✓
    def prefix = task.ext.prefix ?: "${meta.id}"  // Padrão
```
