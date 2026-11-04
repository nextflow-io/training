# Parte 2: Configurar a execução do pipeline

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tradução assistida por IA - [saiba mais e sugira melhorias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Em [Parte 1](./01_run_demo.md), você encontrou e executou o pipeline nf-core/demo usando seu perfil de teste.
Agora vamos ver como configurar a execução do pipeline: definir parâmetros, entender a validação e personalizar a alocação de recursos e os argumentos das ferramentas.

Como explicado em [Hello Config](../hello_nextflow/06_hello_config.md), queremos poder alterar em quais dados nosso pipeline será executado e como ele será executado sem modificar o código do pipeline em si.
Para isso, o Nextflow oferece várias formas de controlar a configuração do pipeline, o que pode ser um pouco avassalador.

O projeto nf-core especifica convenções para organizar os elementos de configuração, distinguindo dois tipos de configuração no nível superior: **parâmetros do pipeline** e **configuração** no sentido estrito.

- **Parâmetros do pipeline** (definidos através do sistema `params`) geralmente incluem coisas como arquivos de entrada, flags de comportamento de ferramentas e parâmetros de análise.
- **Configuração** no sentido estrito refere-se à logística de como o pipeline é executado, ou seja, o executor, alocações de recursos computacionais e assim por diante.

<figure class="excalidraw">
    --8<-- "docs/en/docs/nfcore_use/img/params_vs_config.excalidraw.svg"
</figure>

Vamos começar pelos parâmetros do pipeline e depois veremos a configuração no sentido estrito.

---

## 1. Parâmetros do pipeline

Para todos os pipelines nf-core, você pode obter uma lista completa dos parâmetros do pipeline diretamente pela linha de comando usando a flag `--help`, que é ela própria um parâmetro do pipeline.

### 1.1. Obter a lista de parâmetros com `--help`

Execute o comando de ajuda para o pipeline demo:

```bash
nextflow run nf-core/demo --help
```

??? success "Saída do comando"

    ```console
     N E X T F L O W   ~  version 26.04.4

    Launching `https://github.com/nf-core/demo` [adoring_meucci] revision: 32893afef8 [master]


    ------------------------------------------------------
                                            ,--./,-.
            ___     __   __   __   ___     /,-._.--~'
      |\ | |__  __ /  ` /  \ |__) |__         }  {
      | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                            `._,._,'
      nf-core/demo 1.2.0
    ------------------------------------------------------
    Typical pipeline command:

      nextflow run nf-core/demo -profile <docker/singularity/.../institute> --input samplesheet.csv --outdir <OUTDIR>


    Input/output options
      --input                       [string] Path to a metadata file containing information about the samples in the experiment.
      --outdir                      [string] The output directory where the results will be saved. You have to use absolute paths to storage on Cloud infrastructure.
      --email                       [string] Email address for completion summary.
      --multiqc_title               [string] MultiQC report title. Printed as page header, used for filename if not otherwise specified.

    Reference genome options
      --genome                      [string] Name of iGenomes reference.
      --fasta                       [string] Path to FASTA genome file.

    Process skipping options
      --skip_trim                   [boolean] Skip trimming fastq files with seqtk

    Generic options
      --multiqc_methods_description [string]          Custom MultiQC yaml file containing HTML including a methods description.
      --help                        [boolean, string] Display the help message.
      --help_full                   [boolean]         Display the full detailed help message.
      --show_hidden                 [boolean]         Display hidden parameters in the help message (only works when --help or --help_full are provided).
     !! Hiding 19 param(s), use the `--showHidden` parameter to show them !!
    ------------------------------------------------------

    * The pipeline
        https://doi.org/10.5281/zenodo.12192442

    * The nf-core framework
        https://doi.org/10.1038/s41587-020-0439-x

    * Software dependencies
        https://github.com/nf-core/demo/blob/master/CITATIONS.md
    ```

Como você pode ver, a saída agrupa os parâmetros em categorias (Input/output options, Reference genome options, etc.) com tipos e descrições para cada um.

Essa categorização é determinada por um arquivo de schema, que é abordado mais adiante.
Em pipelines Nextflow simples, `--help` só funciona se o desenvolvedor o implementou manualmente.

!!! tip "Dica"

    Use `--help --show_hidden` para ver parâmetros adicionais que ficam ocultos por padrão, como `--publish_dir_mode` ou `--monochrome_logs`.

### 1.2. Definir valores de parâmetros

Como abordado em [Hello Config](../hello_nextflow/06_hello_config.md), você pode definir valores de parâmetros na linha de comando com `--param_name` ou reunir um conjunto de parâmetros em um arquivo YAML e passá-lo com `-params-file`.
Ambas as abordagens funcionam da mesma forma com pipelines nf-core.

Por exemplo, para pular a etapa de trimagem, queremos definir o parâmetro booleano `skip_trim` como `true`.
Um arquivo de parâmetros chamado `my_params.yml` está disponível no seu diretório de trabalho com esse valor já definido:

```yaml title="my_params.yml"
skip_trim: true
```

Passe-o com `-params-file`:

```bash
nextflow run nf-core/demo -profile test,docker --outdir demo-results-notrim -params-file my_params.yml
```

??? success "Saída do comando"

    ```console
     N E X T F L O W   ~  version 26.04.4

    Launching `https://github.com/nf-core/demo` [focused_heisenberg] revision: 32893afef8 [master]


    ------------------------------------------------------
                                            ,--./,-.
            ___     __   __   __   ___     /,-._.--~'
      |\ | |__  __ /  ` /  \ |__) |__         }  {
      | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                            `._,._,'
      nf-core/demo 1.2.0
    ------------------------------------------------------

    Input/output options
      input                     : https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/samplesheet/samplesheet_test_illumina_amplicon.csv
      outdir                    : demo-results-notrim

    Process skipping options
      skip_trim                 : true

    Institutional config options
      config_profile_name       : Test profile
      config_profile_description: Minimal test dataset to check pipeline function

    Generic options
      trace_report_suffix       : 2026-07-03_22-08-47

    Core Nextflow options
      revision                  : master
      runName                   : focused_heisenberg
      containerEngine           : docker
      launchDir                 : /workspaces/training/nfcore-use
      workDir                   : /workspaces/training/nfcore-use/work
      projectDir                : /workspaces/.nextflow/assets/.repos/nf-core/demo/clones/32893afef8076a03a2767a020b3f0cab2e0b40b2
      userName                  : root
      profile                   : test,docker
      configFiles               : /workspaces/.nextflow/assets/.repos/nf-core/demo/clones/32893afef8076a03a2767a020b3f0cab2e0b40b2/nextflow.config

    !! Only displaying parameters that differ from the pipeline defaults !!
    ------------------------------------------------------

    * The pipeline
        https://doi.org/10.5281/zenodo.12192442

    * The nf-core framework
        https://doi.org/10.1038/s41587-020-0439-x

    * Software dependencies
        https://github.com/nf-core/demo/blob/master/CITATIONS.md

    executor >  local (5)
    [7a/f3599e] NFCORE_DEMO:DEMO:FASTQC (SAMPLE3_SE) | 3 of 3 ✔
    [b0/2f0bdc] NFCORE_DEMO:DEMO:COWPY               | 1 of 1 ✔
    [c3/3c2278] NFCORE_DEMO:DEMO:MULTIQC (demo)      | 1 of 1 ✔
    -[nf-core/demo] Pipeline completed successfully-
    ```

O processo `SEQTK_TRIM` não aparece mais na saída.

!!! warning "Aviso"

    **Limitações importantes sobre entradas de parâmetros**

    **Definindo parâmetros booleanos na linha de comando**

    A partir da versão 26.04 do Nextflow, todos os valores fornecidos na linha de comando são tipados como strings.
    Para um parâmetro booleano como `skip_trim`, passá-lo como uma flag simples (`--skip_trim`) ou como `--skip_trim true` é avaliado como a **string** `"true"`, o que falha na validação do schema:

    ```console
    * --skip_trim (true): Value is [string] but should be [boolean]
    ```

    Para definir um parâmetro booleano com um valor genuíno `true`/`false`, use um `-params-file` como mostrado acima, ou defina-o em um arquivo de configuração.
    Parâmetros do tipo string, integer e file-path não são afetados e ainda podem ser definidos diretamente na linha de comando.
    Este curso usa esse padrão ao longo de todo o material para parâmetros booleanos.

    **Usando arquivos de configuração personalizados**

    Embora seja tecnicamente possível definir parâmetros do pipeline em um arquivo de configuração personalizado passado com `-c`, isso pode não sobrescrever os valores padrão já definidos no `nextflow.config` do próprio pipeline, dependendo das regras de precedência de configuração do Nextflow.
    Usar `--param_name` na linha de comando ou `-params-file` é mais confiável, pois esses sempre têm precedência.

    Como regra geral: se aparece na saída do `--help`, defina-o via linha de comando ou arquivo de parâmetros, e não via arquivo de configuração.

### 1.3. Validação de parâmetros

Curiosidade: o comando `--help` funciona para todos os pipelines nf-core porque o projeto nf-core exige que os desenvolvedores definam formalmente todos os parâmetros do pipeline em um arquivo de schema JSON (`nextflow_schema.json`).
Esse schema registra o tipo, a descrição, o valor padrão e o agrupamento de cada parâmetro.

Além de alimentar a saída do `--help`, o arquivo de schema também permite a validação automatizada no momento da execução.
Isso significa que o Nextflow pode verificar se cada parâmetro que você passa existe e recebeu um valor adequado (do tipo correto, dentro do intervalo de valores permitidos, etc.).

Abordamos isso com mais detalhes na [seção de validação de entrada](../nfcore_build/04_input_validation.md), mas você já pode ver isso em ação fornecendo ao pipeline demo uma entrada de parâmetro inválida.

#### 1.3.1. Parâmetros não reconhecidos

Tente passar um parâmetro que não existe:

```bash
nextflow run nf-core/demo -profile test,docker --outdir demo-results --foobar "invalid"
```

A saída do console inclui um aviso:

```console
WARN: The following invalid input values have been detected:

* --foobar: invalid
```

O pipeline ainda é executado, mas o aviso alerta imediatamente que `--foobar` não é um parâmetro reconhecido.
Isso serve para chamar sua atenção para erros de digitação que não interrompem a execução, como usar `--outDir` em vez de `--outdir`, o que pode ajudá-lo a evitar desperdício de tempo e recursos computacionais.

#### 1.3.2. Valores de parâmetros inválidos

A validação também verifica os **valores** dos parâmetros.
O parâmetro `--skip_trim` é uma flag booleana, então passar um valor string faz o pipeline falhar imediatamente:

```bash
nextflow run nf-core/demo -profile test,docker --outdir demo-results --skip_trim yes
```

```console
ERROR ~ Validation of pipeline parameters failed!

 -- Check '.nextflow.log' file for details
The following invalid input values have been detected:

* --skip_trim (yes): Value is [string] but should be [boolean]

 -- Check script '/workspaces/.nextflow/assets/.repos/nf-core/demo/clones/32893afef8076a03a2767a020b3f0cab2e0b40b2/subworkflows/nf-core/utils_nfschema_plugin/main.nf' at line: 68 or see '.nextflow.log' file for more details
```

O pipeline para antes que qualquer processo seja executado, evitando uma execução com falha ou incorreta.
Como observado em [1.2](#12-set-parameter-values), parâmetros booleanos devem ser definidos com um valor genuíno `true`/`false` em um arquivo de parâmetros, em vez de serem passados na linha de comando, pois os valores da linha de comando são tipados como strings.

### 1.4. Validação de entrada

A mesma lógica de validação também pode ser usada para verificar a validade dos arquivos de entrada.
Por exemplo, se um pipeline espera uma samplesheet como sua principal entrada de dados (o que é o caso de muitos, senão da maioria dos pipelines nf-core), o desenvolvedor pode fornecer um schema de entrada (distinto do schema de parâmetros) descrevendo como o arquivo de entrada deve ser estruturado.

Então, em tempo de execução, o Nextflow pode verificar se o arquivo de entrada fornecido é válido.

Também abordamos isso com mais detalhes na [seção de validação de entrada](../nfcore_build/04_input_validation.md), mas você já pode ver isso em ação fornecendo ao pipeline demo uma samplesheet de entrada inválida.

O pipeline `nf-core/demo` espera um arquivo CSV com as colunas `sample`, `fastq_1` e `fastq_2`.
Isso é definido em um arquivo de schema (`assets/schema_input.json`) que especifica a estrutura esperada, os tipos de colunas e as restrições.

??? abstract "Arquivo de schema para entradas"

    ```json title="assets/schema_input.json"
    {
        "$schema": "https://json-schema.org/draft/2020-12/schema",
        "$id": "https://raw.githubusercontent.com/nf-core/demo/master/assets/schema_input.json",
        "title": "nf-core/demo pipeline - params.input schema",
        "description": "Schema for the file provided with params.input",
        "type": "array",
        "items": {
            "type": "object",
            "properties": {
                "sample": {
                    "type": "string",
                    "pattern": "^\\S+$",
                    "errorMessage": "Sample name must be provided and cannot contain spaces",
                    "meta": ["id"]
                },
                "fastq_1": {
                    "type": "string",
                    "format": "file-path",
                    "exists": true,
                    "pattern": "^([\\S\\s]*\\/)?[^\\s\\/]+\\.f(ast)?q\\.gz$",
                    "errorMessage": "FastQ file for reads 1 must be provided, cannot contain spaces and must have extension '.fq.gz' or '.fastq.gz'"
                },
                "fastq_2": {
                    "type": "string",
                    "format": "file-path",
                    "exists": true,
                    "pattern": "^([\\S\\s]*\\/)?[^\\s\\/]+\\.f(ast)?q\\.gz$",
                    "errorMessage": "FastQ file for reads 2 cannot contain spaces and must have extension '.fq.gz' or '.fastq.gz'"
                }
            },
            "required": ["sample", "fastq_1"]
        }
    }
    ```

O schema especifica que `sample` e `fastq_1` são obrigatórios, enquanto `fastq_2` é opcional (suportando dados paired-end e single-end).
Os caminhos de arquivo são validados quanto à existência e ao padrão de extensão.

Para demonstrar isso, fornecemos uma samplesheet malformada chamada `malformed_samplesheet.csv` no seu diretório de trabalho:

```csv title="malformed_samplesheet.csv"
sample,fastq_2
SAMPLE1,/not/a/real/file.fastq.gz
```

Essa samplesheet está faltando a coluna obrigatória `fastq_1` e tem um caminho de arquivo inexistente em `fastq_2`.

Execute o pipeline demo usando `malformed_samplesheet.csv` como entrada:

```bash
nextflow run nf-core/demo -profile test,docker --outdir demo-results --input malformed_samplesheet.csv
```

```console
ERROR ~ Validation of pipeline parameters failed!

 -- Check '.nextflow.log' file for details
The following invalid input values have been detected:

* --input (malformed_samplesheet.csv): Validation of file failed:
    -> Entry 1: Error for field 'fastq_2' (/not/a/real/file.fastq.gz): the file or directory '/not/a/real/file.fastq.gz' does not exist (FastQ file for reads 2 cannot contain spaces and must have extension '.fq.gz' or '.fastq.gz')
    -> Entry 1: Missing required field(s): fastq_1

 -- Check script '/workspaces/.nextflow/assets/.repos/nf-core/demo/clones/32893afef8076a03a2767a020b3f0cab2e0b40b2/subworkflows/nf-core/utils_nfschema_plugin/main.nf' at line: 68 or see '.nextflow.log' file for more details
```

Como você pode ver, o pipeline falha imediatamente e reporta **todos** os erros de validação de uma vez.
O nf-schema não para no primeiro erro — ele coleta todos os problemas e os lista juntos, para que você possa corrigir tudo de uma vez em vez de descobrir os problemas um por um.

Cada erro identifica exatamente a entrada e o campo que causou o problema, para que você possa corrigir sua samplesheet e relançar o pipeline com a confiança de que ele não vai falhar em algum momento posterior quando o Nextflow for acessar o caminho do arquivo.

Para desenvolvedores, tudo isso é abordado com mais detalhes na [Parte 4 de Build with nf-core](../nfcore_build/04_input_validation.md).

### Conclusão

Você sabe como obter uma lista completa dos parâmetros de um pipeline com `--help`, defini-los via linha de comando ou arquivo de parâmetros, e como o Nextflow valida tanto os valores dos parâmetros quanto os arquivos de entrada em relação aos schemas do pipeline.

### O que vem a seguir?

Aprenda sobre o outro tipo de configuração: como o pipeline é executado, abrangendo alocação de recursos e argumentos de ferramentas.

---

## 2. Configuração

A configuração no sentido estrito controla **como** o pipeline é executado: alocação de recursos, argumentos específicos de ferramentas, onde os jobs são executados e qual sistema de empacotamento de software usar.

Os pipelines nf-core incluem configuração padrão em `nextflow.config` e no diretório `conf/`.
Antes de sobrescrever qualquer coisa, é útil saber onde ficam os valores padrão.

### 2.1. Explorar os arquivos de configuração

Você já viu em [Parte 1](./01_run_demo.md) que o código-fonte do pipeline fica em `$NXF_HOME/assets`.
Usando o link simbólico `pipelines` que você criou em [Parte 1](./01_run_demo.md), liste os arquivos de configuração para ver o que está disponível:

```bash
ls pipelines/nf-core/demo/conf/
```

```console
base.config
containers_conda_lock_files_amd64.config
containers_conda_lock_files_arm64.config
containers_docker_amd64.config
containers_docker_arm64.config
containers_singularity_https_amd64.config
containers_singularity_https_arm64.config
containers_singularity_oras_amd64.config
containers_singularity_oras_arm64.config
igenomes.config
igenomes_ignored.config
modules.config
test.config
test_full.config
```

<figure class="excalidraw">
--8<-- "docs/en/docs/nfcore_use/img/nfcore_config_files.excalidraw.svg"
</figure>

Os arquivos de configuração mais importantes são:

- **`conf/base.config`**: Define labels de recursos (`process_low`, `process_medium`, `process_high`) que atribuem CPUs, memória e tempo aos processos. Quando você vê um processo usando mais recursos do que o esperado, é aqui que esses valores padrão estão definidos.
- **`conf/modules.config`**: Define os argumentos de ferramentas por processo (`ext.args`) e as configurações de publicação de saída (`publishDir`). Abra este arquivo para ver quais argumentos cada ferramenta recebe por padrão.
- **`conf/test.config`**: O perfil de teste que você usou em [Parte 1](./01_run_demo.md), que limita os recursos via `resourceLimits` e define uma samplesheet de teste. Ativado com `-profile test`.
  Há também um `conf/test_full.config` para executar com um conjunto de dados de teste completo, útil para benchmarking.

O `nextflow.config` central carrega todos os arquivos acima e define os valores padrão apropriados para tudo.

Se você quiser modificar qualquer uma das configurações especificadas nesses arquivos, não modifique nenhum deles diretamente.
Em vez disso, crie seu próprio arquivo de configuração e passe-o com `-c`.
Os valores que você especificar vão sobrescrever os valores padrão definidos nesses outros arquivos.

Vamos experimentar isso na prática.

### 2.2. Personalizar recursos de processos e argumentos de ferramentas

Os módulos nf-core suportam dois tipos comuns de substituição de configuração: **alocação de recursos** (CPUs, memória, tempo) e **argumentos de ferramentas** via `ext.args`.

Muitas ferramentas de linha de comando têm argumentos que não são usados com frequência suficiente para serem expostos como parâmetros do pipeline.
A convenção `ext.args` permite que você passe esses argumentos para a ferramenta subjacente através de um arquivo de configuração.

O arquivo `custom.config` fornecido no seu diretório de trabalho demonstra ambas as substituições:

```groovy title="custom.config" linenums="1"
process {
    withName: 'FASTQC' {
        cpus = 2
        memory = 4.GB
    }
    withName: 'SEQTK_TRIM' {
        ext.args = '-b 5'
    }
}
```

O primeiro bloco sobrescreve a alocação de recursos do `FASTQC`.
Por padrão, o `FASTQC` usa o label `process_medium` do `base.config`, que aloca 6 CPUs e 36 GB de memória; aqui limitamos a 2 CPUs e 4 GB.

O segundo bloco passa um argumento extra para o `SEQTK_TRIM` via `ext.args`.
A flag `-b 5` instrui o `seqtk trimfq` a remover 5 bases do início de cada read, além da trimagem por qualidade.

Execute o pipeline com esta configuração:

```bash
nextflow run nf-core/demo -profile test,docker --outdir demo-results-custom -c custom.config
```

??? success "Saída do comando"

    ```console
    executor >  local (8)
    [95/b32876] NFCORE_DEMO:DEMO:FASTQC (SAMPLE1_PE)     | 3 of 3 ✔
    [17/428668] NFCORE_DEMO:DEMO:SEQTK_TRIM (SAMPLE1_PE) | 3 of 3 ✔
    [cf/85991a] NFCORE_DEMO:DEMO:COWPY                   | 1 of 1 ✔
    [3c/94a7a0] NFCORE_DEMO:DEMO:MULTIQC (demo)          | 1 of 1 ✔
    -[nf-core/demo] Pipeline completed successfully-
    ```

A flag `-c` adiciona sua configuração sobre a configuração integrada do pipeline.

Para verificar se a substituição do `ext.args` teve efeito, encontre o hash do diretório de trabalho do `SEQTK_TRIM` na saída da execução (por exemplo, `work/17/428668...`) e verifique o arquivo `.command.sh` dentro dele:

```bash
cat work/17/428668/.command.sh
```

??? success "Saída do comando"

    ```console
    #!/usr/bin/env bash -e -u -o pipefail
    printf "%s\n" sample1_R1.fastq.gz sample1_R2.fastq.gz | while read f;
    do
        seqtk \
            trimfq \
            -b 5 \
            $f \
            | gzip --no-name > SAMPLE1_PE_$(basename $f)
    done
    ...
    ```

Você deve ver `-b 5` no comando `seqtk trimfq`.

Uma coisa importante a saber sobre `ext.args`: se um módulo já tem um valor padrão definido, seu valor vai **substituí-lo completamente** em vez de ser acrescentado a ele.
Por exemplo, o `FASTQC` tem `ext.args = '--quiet'` definido por padrão em `conf/modules.config`:

```groovy title="conf/modules.config" linenums="21" hl_lines="2"
    withName: FASTQC {
        ext.args   = '--quiet'
        publishDir = [
            path: { "${params.outdir}/fastqc/${meta.id}" },
            mode: params.publish_dir_mode,
            pattern: "*.{html,json}",
        ]
    }
```

Se você definir `ext.args = '--kmers 8'` para o `FASTQC`, a flag `--quiet` não será mais aplicada.
Para manter ambas, defina `ext.args = '--quiet --kmers 8'`.

Você deve sempre verificar a configuração padrão de um módulo antes de sobrescrever `ext.args`.

### Conclusão

Você sabe onde ficam os valores padrão de configuração dos pipelines nf-core e como sobrescrever alocações de recursos e argumentos de ferramentas com um arquivo de configuração personalizado.

### O que vem a seguir?

Siga para a [Parte 3](./03_run_production_pipeline.md), onde você aplicará o que aprendeu a um pipeline de produção real.

---

## Resumo

Nesta parte você aprendeu a:

- Obter ajuda, definir parâmetros e entender a validação de parâmetros e de entradas
- Personalizar a alocação de recursos e os argumentos de ferramentas através de arquivos de configuração
