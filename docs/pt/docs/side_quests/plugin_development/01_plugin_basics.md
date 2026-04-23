# Parte 1: Conceitos Básicos de Plugins

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tradução assistida por IA - [saiba mais e sugira melhorias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Nesta seção, você aprenderá como os plugins estendem o Nextflow e, em seguida, experimentará três plugins diferentes para vê-los em ação.

---

## 1. Como os plugins funcionam

Os plugins estendem o Nextflow por meio de vários tipos de extensão:

| Tipo de Extensão   | O que faz                                                        | Exemplo                                 |
| ------------------ | ---------------------------------------------------------------- | --------------------------------------- |
| Funções            | Adiciona funções personalizadas chamáveis nos fluxos de trabalho | `samplesheetToList()`                   |
| Monitores de fluxo | Respondem a eventos como conclusão de tarefas                    | Logging personalizado, alertas no Slack |
| Executors          | Adicionam backends de execução de tarefas                        | AWS Batch, Kubernetes                   |
| Filesystems        | Adicionam backends de armazenamento                              | S3, Azure Blob                          |

Funções e monitores de fluxo de trabalho (chamados de "trace observers" na API do Nextflow) são os tipos mais comuns para autores de plugins.
Executors e filesystems são tipicamente criados por fornecedores de plataformas.

Os próximos exercícios mostram plugins de função e um plugin observer, para que você possa ver os dois tipos em ação.

---

## 2. Usar plugins de função

Plugins de função adicionam funções chamáveis que você importa nos seus fluxos de trabalho.
Você experimentará dois: nf-hello (um exemplo simples) e nf-schema (um plugin real amplamente utilizado).
Ambos os exercícios modificam o mesmo pipeline `hello.nf`, para que você possa ver como os plugins aprimoram um fluxo de trabalho existente.

### 2.1. nf-hello: substituir código escrito manualmente

O plugin [nf-hello](https://github.com/nextflow-io/nf-hello) fornece uma função `randomString` que gera strings aleatórias.
O pipeline já define sua própria versão inline dessa função, que você substituirá pela versão do plugin.

#### 2.1.1. Ver o ponto de partida

Veja o pipeline:

```bash
cat hello.nf
```

```groovy title="Output"
#!/usr/bin/env nextflow

params.input = 'greetings.csv'

/**
 * Gera uma string alfanumérica aleatória
 */
def randomString(int length) {
    def chars = ('a'..'z') + ('A'..'Z') + ('0'..'9')
    def random = new Random()
    return (1..length).collect { chars[random.nextInt(chars.size())] }.join()
}

process SAY_HELLO {
    input:
        val greeting
    output:
        stdout
    script:
    """
    echo '$greeting'
    """
}

workflow {
    greeting_ch = channel.fromPath(params.input)
        .splitCsv(header: true)
        .map { row -> "${row.greeting}_${randomString(8)}" }
    SAY_HELLO(greeting_ch)
    SAY_HELLO.out.view { result -> "Output: ${result.trim()}" }
}
```

O pipeline define sua própria função `randomString` inline e a utiliza para adicionar um ID aleatório a cada saudação.

Execute-o:

```bash
nextflow run hello.nf
```

```console title="Output"
Output: Hello_aBcDeFgH
Output: Bonjour_xYzWvUtS
Output: Holà_qRsPdMnK
Output: Ciao_jLhGfEcB
Output: Hallo_tNwOiAuR
```

A ordem da saída e as strings aleatórias serão diferentes, e se você executar o script novamente obterá um conjunto diferente de saudações aleatórias.

#### 2.1.2. Configurar o plugin

Substitua a função inline por uma do plugin. Adicione isto ao seu `nextflow.config`:

```groovy title="nextflow.config"
// Configuração para os exercícios de desenvolvimento de plugins
plugins {
    id 'nf-hello@0.5.0'
}
```

Os plugins são declarados no `nextflow.config` usando o bloco `plugins {}`.
O Nextflow os baixa automaticamente do [Nextflow Plugin Registry](https://registry.nextflow.io/), um repositório central de plugins da comunidade e oficiais.

#### 2.1.3. Usar a função do plugin

Substitua a função `randomString` inline pela versão do plugin:

=== "Depois"

    ```groovy title="hello.nf" hl_lines="3"
    #!/usr/bin/env nextflow

    include { randomString } from 'plugin/nf-hello'

    params.input = 'greetings.csv'

    process SAY_HELLO {
        input:
            val greeting
        output:
            stdout
        script:
        """
        echo '$greeting'
        """
    }

    workflow {
        greeting_ch = channel.fromPath(params.input)
            .splitCsv(header: true)
            .map { row -> "${row.greeting}_${randomString(8)}" }
        SAY_HELLO(greeting_ch)
        SAY_HELLO.out.view { result -> "Output: ${result.trim()}" }
    }
    ```

=== "Antes"

    ```groovy title="hello.nf" hl_lines="5-12"
    #!/usr/bin/env nextflow

    params.input = 'greetings.csv'

    /**
     * Gera uma string alfanumérica aleatória
     */
    def randomString(int length) {
        def chars = ('a'..'z') + ('A'..'Z') + ('0'..'9')
        def random = new Random()
        return (1..length).collect { chars[random.nextInt(chars.size())] }.join()
    }

    process SAY_HELLO {
        input:
            val greeting
        output:
            stdout
        script:
        """
        echo '$greeting'
        """
    }

    workflow {
        greeting_ch = channel.fromPath(params.input)
            .splitCsv(header: true)
            .map { row -> "${row.greeting}_${randomString(8)}" }
        SAY_HELLO(greeting_ch)
        SAY_HELLO.out.view { result -> "Output: ${result.trim()}" }
    }
    ```

A instrução `include` importa `randomString` de uma biblioteca que é comprovada, testada e mantida por um grupo mais amplo de colaboradores que podem identificar e corrigir bugs.
Em vez de cada pipeline manter sua própria cópia da função, todo pipeline que usa o plugin obtém a mesma implementação verificada.
Isso reduz a duplicação de código e o esforço de manutenção que vem com ela.
A sintaxe `#!groovy include { function } from 'plugin/plugin-id'` é o mesmo `include` usado para módulos Nextflow, com o prefixo `plugin/`.
Você pode ver o [código-fonte de `randomString`](https://github.com/nextflow-io/nf-hello/blob/e67bddebfa589c7ae51f41bf780c92068dc09e93/plugins/nf-hello/src/main/nextflow/hello/HelloExtension.groovy#L110) no repositório nf-hello no GitHub.

#### 2.1.4. Executar

```bash
nextflow run hello.nf
```

```console title="Output"
Pipeline is starting! 🚀
Output: Hello_yqvtclcc
Output: Bonjour_vwwpyzcs
Output: Holà_wrghmgab
Output: Ciao_noniajuy
Output: Hallo_tvrtuxtp
Pipeline complete! 👋
```

(Suas strings aleatórias serão diferentes.)

A saída ainda tem sufixos aleatórios, mas agora `randomString` vem do plugin nf-hello em vez do código inline.
As mensagens "Pipeline is starting!" e "Pipeline complete!" são novas.
Elas vêm do componente observer do plugin, que você explorará na Parte 5.

O Nextflow baixa os plugins automaticamente na primeira vez que são usados, então qualquer pipeline que declare `nf-hello@0.5.0` obtém exatamente a mesma função `randomString` testada sem precisar copiar código entre projetos.

Você agora viu os três passos para usar um plugin de função: declará-lo no `nextflow.config`, importar a função com `include` e chamá-la no seu fluxo de trabalho.
O próximo exercício aplica esses mesmos passos a um plugin do mundo real.

### 2.2. nf-schema: análise de CSV com validação

O plugin [nf-schema](https://github.com/nextflow-io/nf-schema) é um dos plugins Nextflow mais amplamente utilizados.
Ele fornece `samplesheetToList`, uma função que analisa arquivos CSV/TSV usando um JSON schema que define as colunas e tipos esperados.

O pipeline atualmente lê `greetings.csv` usando `splitCsv` e um `map` manual, mas o nf-schema pode substituir isso por uma análise validada e orientada por schema.
Um arquivo JSON schema (`greetings_schema.json`) já está disponível no diretório do exercício.

??? info "O que é um schema?"

    Um schema é uma descrição formal de como dados válidos devem ser.
    Ele define coisas como quais colunas são esperadas, qual tipo cada valor deve ter (string, número, etc.) e quais campos são obrigatórios.

    Pense nele como um contrato: se os dados de entrada não corresponderem ao schema, a ferramenta pode detectar o problema cedo, em vez de deixá-lo causar erros confusos mais adiante no pipeline.

#### 2.2.1. Ver o schema

```bash
cat greetings_schema.json
```

```json title="Output"
{
  "$schema": "https://json-schema.org/draft/2020-12/schema",
  "type": "array",
  "items": {
    "type": "object",
    "properties": {
      "greeting": {
        "type": "string",
        "description": "The greeting text"
      },
      "language": {
        "type": "string",
        "description": "The language of the greeting"
      }
    },
    "required": ["greeting"]
  }
}
```

O schema define duas colunas (`greeting` e `language`) e marca `greeting` como obrigatória.
Se alguém passar um CSV sem a coluna `greeting`, o nf-schema detecta o erro antes que o pipeline seja executado.

#### 2.2.2. Adicionar nf-schema à configuração

Atualize o `nextflow.config` para incluir ambos os plugins:

=== "Depois"

    ```groovy title="nextflow.config" hl_lines="3"
    plugins {
        id 'nf-hello@0.5.0'
        id 'nf-schema@2.6.1'
    }
    ```

=== "Antes"

    ```groovy title="nextflow.config"
    plugins {
        id 'nf-hello@0.5.0'
    }
    ```

#### 2.2.3. Atualizar hello.nf para usar samplesheetToList

Substitua a entrada com `splitCsv` por `samplesheetToList`:

=== "Depois"

    ```groovy title="hello.nf" hl_lines="4 20 21 22"
    #!/usr/bin/env nextflow

    include { randomString } from 'plugin/nf-hello'
    include { samplesheetToList } from 'plugin/nf-schema'

    params.input = 'greetings.csv'

    process SAY_HELLO {
        input:
            val greeting
        output:
            stdout
        script:
        """
        echo '$greeting'
        """
    }

    workflow {
        def samplesheet_list = samplesheetToList(params.input, 'greetings_schema.json')
        greeting_ch = Channel.fromList(samplesheet_list)
            .map { row -> "${row[0]}_${randomString(8)}" }
        SAY_HELLO(greeting_ch)
        SAY_HELLO.out.view { result -> "Output: ${result.trim()}" }
    }
    ```

=== "Antes"

    ```groovy title="hello.nf" hl_lines="19 20 21"
    #!/usr/bin/env nextflow

    include { randomString } from 'plugin/nf-hello'

    params.input = 'greetings.csv'

    process SAY_HELLO {
        input:
            val greeting
        output:
            stdout
        script:
        """
        echo '$greeting'
        """
    }

    workflow {
        greeting_ch = channel.fromPath(params.input)
            .splitCsv(header: true)
            .map { row -> "${row.greeting}_${randomString(8)}" }
        SAY_HELLO(greeting_ch)
        SAY_HELLO.out.view { result -> "Output: ${result.trim()}" }
    }
    ```

O código personalizado de análise com `splitCsv` e `map` é substituído por `samplesheetToList`, uma função comprovada e testada que também valida o samplesheet em relação ao schema antes que o pipeline seja executado.
Isso reduz o esforço de manutenção da lógica de análise escrita manualmente, ao mesmo tempo que melhora a experiência para os usuários do pipeline, que recebem mensagens de erro claras quando a entrada não corresponde ao formato esperado.
Cada linha se torna uma lista de valores na ordem das colunas, então `row[0]` é a saudação e `row[1]` é o idioma.

#### 2.2.4. Executar

```bash
nextflow run hello.nf
```

```console title="Output"
Pipeline is starting! 🚀
Output: Hello_diozjdwm
Output: Bonjour_speathmm
Output: Holà_dllxnzap
Output: Ciao_wzueddzc
Output: Hallo_hsxwrjbh
Pipeline complete! 👋
```

(Suas strings aleatórias serão diferentes.)

A saída é a mesma, mas agora o schema valida a estrutura do CSV antes que o pipeline seja executado.
Em pipelines reais com sample sheets complexos e muitas colunas, esse tipo de validação previne erros que o `splitCsv` + `map` manual deixaria passar.

#### 2.2.5. Ver a validação em ação

Para ver o que a validação de schema detecta, tente introduzir erros no `greetings.csv`.

Renomeie a coluna obrigatória `greeting` para `message`:

```csv title="greetings.csv" hl_lines="1"
message,language
Hello,English
Bonjour,French
Holà,Spanish
Ciao,Italian
Hallo,German
```

Execute o pipeline:

```bash
nextflow run hello.nf
```

```console title="Output"
ERROR ~ Validation of samplesheet failed!

The following errors have been detected in greetings.csv:

-> Entry 1: Missing required field(s): greeting
-> Entry 2: Missing required field(s): greeting
-> Entry 3: Missing required field(s): greeting
-> Entry 4: Missing required field(s): greeting
-> Entry 5: Missing required field(s): greeting
```

O pipeline se recusa a executar porque o schema exige uma coluna `greeting` e não consegue encontrá-la.

Agora restaure a coluna obrigatória, mas renomeie a coluna opcional `language` para `lang`:

```csv title="greetings.csv" hl_lines="1"
greeting,lang
Hello,English
Bonjour,French
Holà,Spanish
Ciao,Italian
Hallo,German
```

```bash
nextflow run hello.nf
```

Desta vez o pipeline é executado, mas exibe um aviso:

```console title="Output (partial)"
WARN: Found the following unidentified headers in greetings.csv:
	- lang
```

Colunas obrigatórias causam erros graves; colunas opcionais causam avisos.
Esse é o tipo de feedback antecipado que economiza tempo de depuração em pipelines reais com dezenas de colunas.

#### 2.2.6. Configurar o comportamento de validação

O aviso sobre `lang` é útil, mas você pode controlar sua severidade por meio da configuração.
Os plugins podem incluir seus próprios escopos de configuração que controlam seu comportamento.
O plugin nf-schema inclui o escopo de configuração `validation`; modificando as configurações aqui você pode alterar como o nf-schema se comporta.

Adicione um bloco `validation` ao `nextflow.config` para fazer com que cabeçalhos não reconhecidos causem um erro em vez de um aviso:

=== "Depois"

    ```groovy title="nextflow.config" hl_lines="6-10"
    plugins {
        id 'nf-hello@0.5.0'
        id 'nf-schema@2.6.1'
    }

    validation {
        logging {
            unrecognisedHeaders = "error"
        }
    }
    ```

=== "Antes"

    ```groovy title="nextflow.config"
    plugins {
        id 'nf-hello@0.5.0'
        id 'nf-schema@2.6.1'
    }
    ```

Execute o pipeline novamente com a mesma coluna `lang` ainda presente:

```bash
nextflow run hello.nf
```

```console title="Output (partial)"
Found the following unidentified headers in greetings.csv:
	- lang
 -- Check script 'hello.nf' at line: 20 or see '.nextflow.log' file for more details
```

O pipeline agora falha em vez de apenas avisar.
O código do pipeline não mudou; apenas a configuração foi alterada.

Restaure o `greetings.csv` ao seu estado original e remova o bloco `validation` antes de continuar:

```csv title="greetings.csv"
greeting,language
Hello,English
Bonjour,French
Holà,Spanish
Ciao,Italian
Hallo,German
```

```groovy title="nextflow.config"
plugins {
    id 'nf-hello@0.5.0'
    id 'nf-schema@2.6.1'
}
```

Tanto o nf-hello quanto o nf-schema são plugins de função: eles fornecem funções que você importa com `include` e chama no seu código de fluxo de trabalho.
O próximo exercício mostra um tipo diferente de plugin que funciona sem nenhuma instrução `include`.

---

## 3. Usar um plugin observer: nf-co2footprint

Nem todos os plugins fornecem funções para importar.
O plugin [nf-co2footprint](https://github.com/nextflow-io/nf-co2footprint) usa um **trace observer** para monitorar o uso de recursos do seu pipeline e estimar sua pegada de carbono.
Você não precisa alterar nenhum código do pipeline; basta adicioná-lo à configuração.

### 3.1. Adicionar nf-co2footprint à configuração

Atualize o `nextflow.config`:

=== "Depois"

    ```groovy title="nextflow.config" hl_lines="4"
    plugins {
        id 'nf-hello@0.5.0'
        id 'nf-schema@2.6.1'
        id 'nf-co2footprint@1.2.0'
    }
    ```

=== "Antes"

    ```groovy title="nextflow.config"
    plugins {
        id 'nf-hello@0.5.0'
        id 'nf-schema@2.6.1'
    }
    ```

### 3.2. Executar o pipeline

```bash
nextflow run hello.nf
```

O plugin produz várias mensagens INFO e WARN durante a execução.
Isso é normal para um exemplo pequeno sendo executado em uma máquina local:

```console title="Output (partial)"
nf-co2footprint plugin  ~  version 1.2.0
WARN - [nf-co2footprint] Target zone null not found. Attempting to retrieve carbon intensity for fallback zone GLOBAL.
INFO - [nf-co2footprint] Using fallback carbon intensity from GLOBAL from CI table: 480.0 gCO₂eq/kWh.
WARN - [nf-co2footprint] Executor 'null' not mapped.
WARN - [nf-co2footprint] Fallback to: `machineType = null`, `pue = 1.0`. ...
...
WARN - [nf-co2footprint] No CPU model detected. Using default CPU power draw value (11.41 W).
WARN - [nf-co2footprint] 🔁 Requested memory is null for task 2. Using maximum consumed memory/`peak_rss` (0 GB) for CO₂e footprint computation.
```

Os avisos sobre zona, executor, modelo de CPU e memória aparecem porque o plugin não consegue detectar todos os detalhes de hardware de um ambiente de treinamento local.
Em um ambiente de produção (por exemplo, um cluster HPC ou nuvem), esses valores estariam disponíveis e as estimativas seriam mais precisas.

Ao final, procure uma linha como:

```console title="Output (partial)"
🌱 The workflow run used 126.76 uWh of electricity, resulting in the release of 60.84 ug of CO₂ equivalents into the atmosphere.
```

(Seus números serão diferentes.)

### 3.3. Ver o relatório

O plugin gera arquivos de saída no seu diretório de trabalho:

```bash
ls co2footprint_*
```

```console title="Output"
co2footprint_report_<timestamp>.html
co2footprint_summary_<timestamp>.txt
co2footprint_trace_<timestamp>.txt
```

Veja o resumo:

```bash
cat co2footprint_summary_*.txt
```

```console title="Output"
Total CO₂e footprint measures of this workflow run (including cached tasks):
  CO₂e emissions: 60.84 ug
  Energy consumption: 126.76 uWh
  CO₂e emissions (market): -

Which equals:
  - 3.48E-7 km travelled by car
  - It takes one tree 0.17s to sequester the equivalent amount of CO₂ from the atmosphere
  - 1.22E-7 % of a flight from Paris to London
```

(Seus números serão diferentes.)

A primeira seção mostra os valores brutos de energia e emissões.
A seção "Which equals" coloca esses números em perspectiva, convertendo-os em equivalentes familiares.
O resumo também inclui uma seção listando as opções de configuração do plugin e uma citação ao artigo de pesquisa [Green Algorithms](https://doi.org/10.1002/advs.202100707) no qual o método de cálculo é baseado.

### 3.4. Configurar o plugin

O aviso "Target zone null" da seção 3.2 apareceu porque o plugin não tinha nenhuma localização configurada.
O plugin nf-co2footprint define um escopo de configuração `co2footprint` onde você pode definir sua localização geográfica.

Adicione um bloco `co2footprint` ao `nextflow.config`:

=== "Depois"

    ```groovy title="nextflow.config" hl_lines="7-9"
    plugins {
        id 'nf-hello@0.5.0'
        id 'nf-schema@2.6.1'
        id 'nf-co2footprint@1.2.0'
    }

    co2footprint {
        location = 'GB'
    }
    ```

=== "Antes"

    ```groovy title="nextflow.config"
    plugins {
        id 'nf-hello@0.5.0'
        id 'nf-schema@2.6.1'
        id 'nf-co2footprint@1.2.0'
    }
    ```

!!! tip "Dica"

    Use o código do seu próprio país se preferir (por exemplo, `'US'`, `'DE'`, `'FR'`).

Execute o pipeline:

```bash
nextflow run hello.nf
```

```console title="Output (partial)"
INFO - [nf-co2footprint] Using fallback carbon intensity from GB from CI table: 163.92 gCO₂eq/kWh.
```

O aviso de zona desapareceu.
O plugin agora usa a intensidade de carbono específica do GB (163.92 gCO₂eq/kWh) em vez do valor global de fallback (480.0 gCO₂eq/kWh).

!!! note "Nota"

    Você também pode ver uma mensagem `WARN: Unrecognized config option 'co2footprint.location'`.
    Isso é apenas cosmético e pode ser ignorado com segurança; o plugin ainda lê o valor corretamente.

Na Parte 6, você criará um escopo de configuração para o seu próprio plugin.

Este plugin funciona inteiramente por meio do mecanismo observer, conectando-se aos eventos do ciclo de vida do fluxo de trabalho para coletar métricas de recursos e gerar seu relatório quando o pipeline é concluído.

Você agora experimentou plugins de função (importados com `include`) e um plugin observer (ativado apenas pela configuração).
Esses são os dois tipos de extensão mais comuns, mas como a tabela na seção 1 mostra, os plugins também podem adicionar executors e filesystems.

---

## 4. Descobrindo plugins

O [Nextflow Plugin Registry](https://registry.nextflow.io/) é o hub central para encontrar plugins disponíveis.

![A página do plugin nf-hello em registry.nextflow.io](img/plugin-registry-nf-hello.png)

Cada página de plugin mostra sua descrição, versões disponíveis, instruções de instalação e links para a documentação.

---

## 5. Preparar para o desenvolvimento de plugins

As seções a seguir (Partes 2-6) usam um arquivo de pipeline separado, `greet.nf`, que depende do nf-schema, mas não do nf-hello ou nf-co2footprint.

Atualize o `nextflow.config` para manter apenas o nf-schema:

```groovy title="nextflow.config"
// Configuração para os exercícios de desenvolvimento de plugins
plugins {
    id 'nf-schema@2.6.1'
}
```

Remova os arquivos de saída do co2footprint:

```bash
rm -f co2footprint_*
```

O arquivo `hello.nf` mantém o trabalho da Parte 1 para referência; daqui em diante, você trabalhará com `greet.nf`.

---

## Conclusão

Você usou três plugins diferentes:

- **nf-hello**: Um plugin de função que fornece `randomString`, importado com `include`
- **nf-schema**: Um plugin de função que fornece `samplesheetToList` para análise de CSV com validação por schema
- **nf-co2footprint**: Um plugin observer que monitora o uso de recursos automaticamente, sem necessidade de `include`

Padrões principais:

- Os plugins são declarados no `nextflow.config` com `#!groovy plugins { id 'plugin-name@version' }`
- Plugins de função requerem `#!groovy include { function } from 'plugin/plugin-id'`
- Plugins observer funcionam automaticamente uma vez declarados na configuração
- Os plugins podem definir escopos de configuração (por exemplo, `#!groovy validation {}`, `#!groovy co2footprint {}`) para personalizar o comportamento
- O [Nextflow Plugin Registry](https://registry.nextflow.io/) lista os plugins disponíveis

---

## O que vem a seguir?

As seções a seguir mostram como construir seu próprio plugin.
Se você não tem interesse no desenvolvimento de plugins, pode parar aqui ou pular para o [Resumo](summary.md).

[Continuar para a Parte 2 :material-arrow-right:](02_create_project.md){ .md-button .md-button--primary }
