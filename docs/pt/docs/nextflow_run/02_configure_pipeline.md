# Parte 2: Configurar o pipeline

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tradução assistida por IA - [saiba mais e sugira melhorias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Na [Parte 1](./01_run_nextflow.md), você executou um pipeline completo de múltiplas etapas que processa várias entradas em paralelo usando contêineres.
Agora vamos ver como configurar o comportamento do pipeline usando `nextflow.config`: primeiro examinando o arquivo de configuração que já fornecemos, depois explorando algumas outras formas de fornecer configuração e, por fim, controlando como e onde as saídas são publicadas.

---

## 1. Examinar o arquivo de configuração principal

O Nextflow automaticamente detecta o `nextflow.config` no diretório de trabalho e aplica suas configurações a cada execução.

Fornecemos um arquivo de configuração que cobre quatro áreas: empacotamento de software, configurações de processo, parâmetros do pipeline e perfis de execução.

??? full-code "nextflow.config"

    ```groovy title="nextflow.config" linenums="1"
    /*
     * Empacotamento de software
     */
    docker.enabled = true

    /*
     * Configurações de processo
     */
    process {
        cpus = 1
        memory = 1.GB
    }

    /*
     * Parâmetros do pipeline
     */
    params {
        input = 'data/greetings.csv'
        batch = 'batch'
        character = 'turkey'
    }

    /*
     * Perfis
     */
    profiles {
        test {
            params.input = 'data/greetings.csv'
            params.batch = 'test'
            params.character = 'tux'
        }
        conda {
            docker.enabled = false
            conda.enabled = true
        }
    }
    ```

Vamos percorrer cada um deles e depois colocar os perfis em uso executando o pipeline com um deles.

!!! note "Nota"

    Esta configuração cobre a execução local em uma única máquina.
    O Nextflow também suporta schedulers HPC (SLURM, PBS, LSF) e executores em nuvem (AWS Batch, Google Cloud Batch, Azure Batch), todos configurados pelo mesmo mecanismo do `nextflow.config`.
    Consulte a [Parte 1: Adaptar ao seu ambiente de computação](../execution_config/01_packaging_and_execution.md) no curso [Execution Config](../execution_config/index.md) para um guia completo dessas opções.

### 1.1. Empacotamento de software

Empacotamento de software é a forma como o Nextflow fornece as ferramentas que seus processos precisam, seja uma imagem de contêiner, um ambiente Conda ou outra coisa.

```groovy title="nextflow.config" linenums="1"
/*
 * Empacotamento de software
 */
docker.enabled = true
```

Esta linha habilita o Docker para todos os processos.
Qualquer processo que declare uma diretiva `container` é executado dentro da imagem especificada.

### 1.2. Configurações de processo

Lembre-se de que um processo é uma única etapa do seu pipeline, como `sayHello` ou `cowpy`.
O Nextflow permite configurar várias coisas sobre como cada um é executado: quanto de CPU e memória ele recebe, qual contêiner ou ambiente Conda ele usa, e muito mais.

```groovy title="nextflow.config" linenums="6"
/*
 * Configurações de processo
 */
process {
    cpus = 1
    memory = 1.GB
}
```

Isso limita cada processo a uma única CPU e 1 GB de memória.

O Nextflow também permite definir valores diferentes para processos individuais nomeados ou grupos de processos; você aprenderá como fazer isso na [Parte 2: Gerenciar recursos de computação e falhas](../execution_config/02_resources_and_retries.md#12-set-resource-allocations-for-a-specific-process) do curso [Execution Config](../execution_config/index.md).

### 1.3. Parâmetros do pipeline

Os parâmetros são as entradas de linha de comando do pipeline — os mesmos flags `--input`, `--batch` e `--character` que você já vinha definindo diretamente na linha de comando.
Definir valores padrão para eles aqui significa que você não precisa digitá-los toda vez; porém, como verá mais adiante nesta parte, há algumas outras formas de fornecê-los também.

```groovy title="nextflow.config" linenums="14"
/*
 * Parâmetros do pipeline
 */
params {
    input = 'data/greetings.csv'
    batch = 'batch'
    character = 'turkey'
}
```

Esses valores padrão entram em ação sempre que um parâmetro não é fornecido na linha de comando, então executar `nextflow run main.nf` sem nenhum flag ainda funciona.

### 1.4. Perfis

Os perfis permitem agrupar um conjunto de configurações sob um único nome, para que você possa alternar entre configurações completas com um único flag em vez de alterar valores manualmente toda vez.

```groovy title="nextflow.config" linenums="23"
/*
 * Perfis
 */
profiles {
    test {
        params.input = 'data/greetings.csv'
        params.batch = 'test'
        params.character = 'tux'
    }
    conda {
        docker.enabled = false
        conda.enabled = true
    }
}
```

O perfil `test` substitui três parâmetros para executar o pipeline com um conjunto de entradas pequeno e bem definido; todo pipeline nf-core vem com um desses para validação rápida, e é uma convenção que vale a pena seguir nos seus próprios pipelines também.

O perfil `conda` troca o empacotamento de software de Docker para Conda.

Você ativa um perfil passando `-profile <nome>` na linha de comando.

Vamos colocar o perfil `test` em uso.

```bash
nextflow run main.nf -profile test
```

??? success "Saída do comando"

    ```console hl_lines="6 7 8 9"
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [reverent_heisenberg] revision: ce74f81996

    executor >  local (8)
    [3d/8a12c7] sayHello (3)       | 3 of 3 ✔
    [e7/e0934f] convertToUpper (2) | 3 of 3 ✔
    [1d/616569] collectGreetings   | 1 of 1 ✔
    [44/7d46cf] cowpy              | 1 of 1 ✔

    Outputs:

      /workspaces/training/nextflow-run/results

      first_output:
        - test/intermediates/Bonjour-output.txt
        - test/intermediates/Hello-output.txt
        - test/intermediates/Hola-output.txt

      uppercased:
        - test/intermediates/UPPER-Bonjour-output.txt
        - test/intermediates/UPPER-Hello-output.txt
        - test/intermediates/UPPER-Hola-output.txt

      collected: test/intermediates/COLLECTED-test-output.txt

      batch_report: test/test-report.txt

      cowpy_art: test/cowpy-COLLECTED-test-output.txt
    ```

O pipeline é executado com `batch = 'test'` e `character = 'tux'`.
Verifique `results/test/`: o nome do batch agora faz parte do próprio caminho do diretório, e a arte ASCII exibe o pinguim tux em vez de um peru.

!!! note "Nota"

    Você pode ativar vários perfis ao mesmo tempo e usar `nextflow config -profile <nome>,<nome>` para ver o resultado totalmente resolvido antes de executar qualquer coisa.
    A combinação de perfis e como o Nextflow resolve conflitos entre eles é abordada em profundidade na [Parte 3: Usar perfis para alternar configurações](../execution_config/03_profiles.md) do curso [Execution Config](../execution_config/index.md).

### Conclusão

Você sabe o que fazem os elementos mais comuns de um arquivo `nextflow.config` e como ativar um perfil.

### O que vem a seguir?

Aprenda algumas outras formas de fornecer valores de configuração sem modificar o `nextflow.config` principal — úteis para configurar execuções individuais e para compartilhar um conjunto exato de configurações com outra pessoa.

---

## 2. Fornecer configuração por meio de arquivos suplementares

Definir valores padrão no `nextflow.config` funciona bem para valores que raramente mudam.
O Nextflow também oferece dois mecanismos mais específicos: um arquivo de configuração específico para uma execução, para adaptar a execução a um ambiente particular, e um arquivo de parâmetros para compartilhar um conjunto exato de valores de entrada com um colaborador.

### 2.1. Usar um arquivo de configuração específico para uma execução

Digamos que você está movendo o pipeline para uma máquina que não tem Docker e quer dar a cada processo mais recursos para trabalhar.
Crie um novo arquivo de configuração com apenas as substituições que você precisa:

```groovy title="custom.config" linenums="1"
process {
    cpus = 2
    memory = 2.GB
}

docker.enabled = false
conda.enabled = true
```

Passe-o junto com seu pipeline principal usando `-c`:

```bash
nextflow run main.nf -c custom.config
```

??? success "Saída do comando"

    ```console hl_lines="6 7 8 9"
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [exotic_cray] revision: ce74f81996

    executor >  local (8)
    [77/e02315] sayHello (1)       | 3 of 3 ✔
    [a6/ccf44b] convertToUpper (3) | 3 of 3 ✔
    [56/fd1296] collectGreetings   | 1 of 1 ✔
    [2c/205a94] cowpy              | 1 of 1 ✔

    Outputs:

      /workspaces/training/nextflow-run/results

      first_output:
        - batch/intermediates/Bonjour-output.txt
        - batch/intermediates/Hello-output.txt
        - batch/intermediates/Hola-output.txt

      uppercased:
        - batch/intermediates/UPPER-Bonjour-output.txt
        - batch/intermediates/UPPER-Hello-output.txt
        - batch/intermediates/UPPER-Hola-output.txt

      collected: batch/intermediates/COLLECTED-batch-output.txt

      batch_report: batch/batch-report.txt

      cowpy_art: batch/cowpy-COLLECTED-batch-output.txt
    ```

O Nextflow mescla o `custom.config` sobre o `nextflow.config` do próprio pipeline, então cada processo agora recebe 2 CPUs e 2 GB de memória em vez dos valores padrão, e é executado via Conda em vez de Docker.
`cowpy` é o único processo com um pacote Conda declarado junto ao seu contêiner, então é o que você verá o Nextflow construir um ambiente para:

```console
Creating env using conda: conda-forge::cowpy==1.1.5 [cache /path/to/work/conda/env-898314d566668b6587ad714ae06b8520]
```

Um arquivo pequeno que apenas substitui a alocação de recursos e o empacotamento, sem tocar nos parâmetros do pipeline, é exatamente o padrão que os pipelines nf-core esperam de configurações institucionais.
Navegue pelo repositório [nf-core/configs](https://github.com/nf-core/configs) para ver exemplos do mundo real.

Isso oferece uma forma descartável de adaptar um pipeline a um novo ambiente sem tocar na sua configuração normal.

### 2.2. Usar um arquivo de parâmetros

Digamos que você precisa compartilhar um conjunto exato de parâmetros de execução com um colaborador, ou registrá-los para uma publicação.

O Nextflow permite fornecer [arquivos de parâmetros](https://nextflow.io/docs/latest/config.html#parameter-file) nos formatos YAML ou JSON, que são uma forma mais simples de distribuir um conjunto exato e reproduzível de valores.

Um arquivo de parâmetros chamado `test-params.yaml` já está disponível no seu diretório de trabalho:

```yaml title="test-params.yaml" linenums="1"
input: "data/greetings.csv"
batch: "yaml"
character: "stegosaurus"
```

A sintaxe usa dois-pontos (`:`) em vez dos sinais de igual (`=`) usados no `nextflow.config`, pois este arquivo é YAML puro em vez de Groovy.

!!! info "Info"

    Uma versão JSON, `test-params.json`, também está disponível. Fique à vontade para experimentá-la por conta própria; a sintaxe para passá-la é idêntica.

Passe o arquivo com `-params-file`:

```bash
nextflow run main.nf -params-file test-params.yaml
```

??? success "Saída do comando"

    ```console
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [sharp_faraday] revision: ce74f81996

    executor >  local (8)
    [1c/9ff63e] sayHello (1)       | 3 of 3 ✔
    [3b/bb5691] convertToUpper (2) | 3 of 3 ✔
    [cd/2c1f6e] collectGreetings   | 1 of 1 ✔
    [89/c333bc] cowpy              | 1 of 1 ✔

    Outputs:

      /workspaces/training/nextflow-run/results

      first_output:
        - yaml/intermediates/Bonjour-output.txt
        - yaml/intermediates/Hello-output.txt
        - yaml/intermediates/Hola-output.txt

      uppercased:
        - yaml/intermediates/UPPER-Bonjour-output.txt
        - yaml/intermediates/UPPER-Hello-output.txt
        - yaml/intermediates/UPPER-Hola-output.txt

      collected: yaml/intermediates/COLLECTED-yaml-output.txt

      batch_report: yaml/yaml-report.txt

      cowpy_art: yaml/cowpy-COLLECTED-yaml-output.txt
    ```

??? abstract "Conteúdo do arquivo"

    ```console title="results/yaml/cowpy-COLLECTED-yaml-output.txt"
     _________
    / BONJOUR \
    | HOLA    |
    \ HELLO   /
     ---------
    \                             .       .
     \                           / `.   .' "
      \                  .---.  <    > <    >  .---.
       \                 |    \  \ - ~ ~ - /  /    |
             _____          ..-~             ~-..-~
            |     |   \~~~\.'                    `./~~~/
           ---------   \__/                        \__/
          .'  O    \     /               /       \  "
         (_____,    `._.'               |         }  \/~~~/
          `----.          /       }     |        /    \__/
                `-.      |       /      |       /      `. ,~~|
                    ~-.__|      /_ - ~ ^|      /- _      `..-'
                        |     /        |     /     ~-.     `-. _  _  _
                        |_____|        |_____|         ~ - . _ _ _ _ _>
    ```

Um arquivo de parâmetros é especialmente valioso quando um pipeline tem mais do que alguns parâmetros: ele permite fornecê-los todos de uma vez, sem uma linha de comando extensa ou qualquer alteração no script do fluxo de trabalho, e é fácil de distribuir junto com seus resultados.

### Conclusão

Você conhece mais duas formas de fornecer configuração: um arquivo de configuração específico para uma execução, para adaptar a execução a um novo ambiente, e um arquivo de parâmetros para compartilhar valores de entrada exatos e reproduzíveis.

### O que vem a seguir?

Aprenda como controlar como e onde as saídas do seu pipeline são publicadas.

---

## 3. Gerenciar as saídas do pipeline

O autor do pipeline decide como as saídas são organizadas no código, mas você não precisa tocar nesse código para controlar onde elas acabam ou como chegam lá.
O Nextflow oferece formas de fazer isso no nível de configuração: definir um diretório base de saída e escolher se os arquivos são copiados ou vinculados simbolicamente.

### 3.1. Personalizar o diretório de saída

Por padrão, o Nextflow publica as saídas em `results/`.
Aponte para outro lugar com `-output-dir` (ou sua forma abreviada, `-o`):

```bash
nextflow run main.nf -output-dir outputs
```

??? success "Saída do comando"

    ```console hl_lines="6 7 8 9"
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [serene_kimura] revision: ce74f81996

    executor >  local (8)
    [31/df5c15] sayHello (3)       | 3 of 3 ✔
    [08/8bb2e5] convertToUpper (1) | 3 of 3 ✔
    [e5/5814da] collectGreetings   | 1 of 1 ✔
    [cf/8ab8c6] cowpy              | 1 of 1 ✔

    Outputs:

      /workspaces/training/nextflow-run/outputs

      first_output:
        - batch/intermediates/Bonjour-output.txt
        - batch/intermediates/Hello-output.txt
        - batch/intermediates/Hola-output.txt

      uppercased:
        - batch/intermediates/UPPER-Bonjour-output.txt
        - batch/intermediates/UPPER-Hello-output.txt
        - batch/intermediates/UPPER-Hola-output.txt

      collected: batch/intermediates/COLLECTED-batch-output.txt

      batch_report: batch/batch-report.txt

      cowpy_art: batch/cowpy-COLLECTED-batch-output.txt
    ```

??? abstract "Conteúdo do diretório"

    ```console
    outputs/batch
    ├── batch-report.txt
    ├── cowpy-COLLECTED-batch-output.txt
    └── intermediates
        ├── Bonjour-output.txt
        ├── COLLECTED-batch-output.txt
        ├── Hello-output.txt
        ├── Hola-output.txt
        ├── UPPER-Bonjour-output.txt
        ├── UPPER-Hello-output.txt
        └── UPPER-Hola-output.txt
    ```

As saídas agora ficam em `outputs/batch/` em vez do padrão `results/batch/`.
O próprio código do pipeline ainda decide a estrutura dentro desse diretório base, como os subdiretórios `batch/` e `intermediates/`; `-output-dir` controla apenas onde essa estrutura começa.

`-output-dir` é basicamente um atalho de linha de comando para a opção de configuração `outputDir`, então pode ser usado em qualquer lugar onde a configuração pode estar: diretamente no `nextflow.config`, dentro de um perfil ou em um arquivo de sobreposição `-c` como o que você usou anteriormente nesta parte.
Por exemplo, este trecho mostra a mesma configuração colocada diretamente no `nextflow.config` em vez de passada na linha de comando:

```groovy title="nextflow.config"
outputDir = 'outputs'
```

Consulte [Configuration file](https://nextflow.io/docs/latest/config.html) na referência do Nextflow para a lista completa de lugares onde uma opção de configuração como essa pode estar.

### 3.2. Escolher como as saídas são publicadas

Por padrão, o Nextflow publica as saídas como links simbólicos que apontam para os locais das saídas em `work/`, não como cópias reais:

```console
$ ls -l results/batch/intermediates/Hello-output.txt
lrwxr-xr-x  ...  Hello-output.txt -> /workspaces/training/nextflow-run/work/b7/b8c4d1.../Hello-output.txt
```

Os autores do pipeline podem definir o 'modo de publicação' como `'copy'` ou `'move'` para cada processo individual no código do fluxo de trabalho.
Eles normalmente fazem isso para as saídas finais do pipeline, deixando o comportamento padrão `'symlink'` para os arquivos intermediários que podem ser excluídos após a execução completa do pipeline.

Isso evita duplicar dados em disco, mas significa que você não pode excluir os diretórios de tarefa em `work/` sem quebrar o link, perdendo a capacidade de usar `-resume`.
Se você quiser que todos os arquivos de saída sejam copiados de verdade, defina [`workflow.output.mode`](https://nextflow.io/docs/latest/reference/config.html#workflow) como `'copy'` na configuração do seu pipeline. (Ao contrário de `-output-dir`, não há um flag de linha de comando para isso; é apenas via configuração.)

Tente definir isso no `nextflow.config`:

=== "Depois"

    ```groovy title="nextflow.config" hl_lines="7"
    params {
        input = 'data/greetings.csv'
        batch = 'batch'
        character = 'turkey'
    }

    workflow.output.mode = 'copy'
    ```

=== "Antes"

    ```groovy title="nextflow.config"
    params {
        input = 'data/greetings.csv'
        batch = 'batch'
        character = 'turkey'
    }
    ```

Em seguida, execute o pipeline alterando o nome do batch para que você possa ver a diferença nas saídas:

```bash
nextflow run main.nf --batch withmode
```

??? success "Saída do comando"

    ```console hl_lines="6 7 8 9"
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [angry_noether] revision: ce74f81996

    executor >  local (8)
    [41/2b478d] sayHello (3)       | 3 of 3 ✔
    [bf/dd2840] convertToUpper (2) | 3 of 3 ✔
    [ea/364e97] collectGreetings   | 1 of 1 ✔
    [10/76fe7b] cowpy              | 1 of 1 ✔

    Outputs:

      /workspaces/training/nextflow-run/results

      first_output:
        - withmode/intermediates/Bonjour-output.txt
        - withmode/intermediates/Hello-output.txt
        - withmode/intermediates/Hola-output.txt

      uppercased:
        - withmode/intermediates/UPPER-Bonjour-output.txt
        - withmode/intermediates/UPPER-Hello-output.txt
        - withmode/intermediates/UPPER-Hola-output.txt

      collected: withmode/intermediates/COLLECTED-withmode-output.txt

      batch_report: withmode/withmode-report.txt

      cowpy_art: withmode/cowpy-COLLECTED-withmode-output.txt
    ```

Dê uma olhada em um dos arquivos de saída como antes:

```console
$ ls -l results/withmode/intermediates/Hello-output.txt
-rw-r--r--  ...  Hello-output.txt
```

Agora é um arquivo real e independente que continuará disponível mesmo que `work/` seja limpo.

!!! warning "Aviso"

    A configuração `workflow.output.mode` apenas preenche um valor padrão para saídas que ainda não têm um modo definido no código do pipeline.
    Ela não pode substituir um modo que o autor definiu diretamente no código, independentemente do que você configurar.

### Conclusão

Você sabe como personalizar o diretório base de saída e escolher entre saídas copiadas e vinculadas simbolicamente, tudo sem tocar no código do pipeline.

### O que vem a seguir?

Siga para a [Parte 3](./03_manage_executions.md), onde você aprenderá como inspecionar o histórico de execuções anteriores, gerar relatórios de execução e limpar diretórios de trabalho antigos.

---

## Resumo

Nesta parte você aprendeu a:

- Configurar o comportamento do pipeline usando `nextflow.config` e perfis
- Fornecer configuração por meio de um arquivo de configuração específico para uma execução ou de um arquivo de parâmetros
- Personalizar o diretório de saída e escolher entre saídas copiadas e vinculadas simbolicamente
