---
description: Material de treinamento básico do Nextflow
---

# Canais

Canais são uma estrutura de dados chave do Nextflow que permite a implementação de fluxos de trabalho computacionais utilizando paradigmas funcional e reativo com base no paradigma de programação [Dataflow](https://en.wikipedia.org/wiki/Dataflow_programming).

Eles são usados para conectar logicamente tarefas entre si ou para implementar transformações de dados de estilo funcional.

<figure class="excalidraw">
--8<-- "docs/basic_training/img/channel-files.excalidraw.svg"
</figure>

## Tipos de canais

O Nextflow distingue dois tipos diferentes de canais: canais de **fila** e canais de **valor**.

### Canal de fila

Um canal de _fila_ é uma fila assíncrona undirecional FIFO (First-in-First-out, o primeiro a entrar, é o primeiro a sair) que conecta dois processos ou operadores.

-   _assíncrono_ significa que as operações ocorrem sem bloqueio.
-   _unidirecional_ significa que os dados fluem do gerador para o consumidor.
-   _FIFO_ significa que os dados são entregues na mesma ordem em que são produzidos. Primeiro a entrar, primeiro a sair.

Um canal de fila é criado implicitamente por definições de saída de um processo ou usando fábricas de canal, como o [Channel.of](https://www.nextflow.io/docs/latest/channel.html#of) ou [Channel.fromPath](https://www.nextflow.io/docs/latest/channel.html#frompath).

Tente os seguintes trechos de código:

!!! info ""

    Clique no ícone :material-plus-circle: no código para ver explicações.

```groovy linenums="1"
canal = Channel.of(1,2,3)
println(canal) // (1)!
canal.view() // (2)!
```

1. Use a função `println` embutida no Nextflow por padrão para imprimir o conteúdo do canal `canal`
2. Aplique o método `view` no canal `canal` para imprimir cada ítem emitido por esse canal

!!! exercise

    Tente executar este trecho de código. Você pode fazer isso criando um novo arquivo `.nf` ou editando um arquivo `.nf` já existente.

    ```groovy linenums="1"
    canal = Channel.of(1,2,3)
    canal.view()
    ```

### Canais de valor

Um canal de **valor** (também conhecido como canal singleton), por definição, está vinculado a um único valor e pode ser lido quantas vezes for necessário sem consumir seu conteúdo. Um canal de `valor` é criado usando a fábrica de canal [value](https://www.nextflow.io/docs/latest/channel.html#value) ou por operadores que retornam um valor apenas, como [first](https://www.nextflow.io/docs/latest/operator.html#first), [last](https://www.nextflow.io/docs/latest/operator.html#last), [collect](https://www.nextflow.io/docs/latest/operator.html#operator-collect), [count](https://www.nextflow.io/docs/latest/operator.html#operator-count), [min](https://www.nextflow.io/docs/latest/operator.html#operator-min), [max](https://www.nextflow.io/docs/latest/operator.html#operator-max), [reduce](https://www.nextflow.io/docs/latest/operator.html#operator-reduce), e [sum](https://www.nextflow.io/docs/latest/operator.html#operator-sum).

Para entender melhor a diferença entre canais de valor e de fila, salve o trecho abaixo como `example.nf`.

```groovy linenums="1" title="example.nf" linenums="1"
canal1 = Channel.of(1,2,3)
canal2 = Channel.of(1)

process SUM {
    input:
    val x
    val y

    output:
    stdout

    script:
    """
    echo \$(($x+$y))
    """
}

workflow {
    SUM(canal1,canal2).view()
}
```

Ao rodar o script, ele imprime apenas 2, como você pode ver abaixo:

```console
2
```

Para entender o motivo, podemos inspecionar o canal de fila executando o Nextflow com DSL1, o que nos dá uma compreensão mais explícita do que está por trás das cortinas.

```groovy linenums="1"
canal1 = Channel.of(1)
println canal1
```

```console
$ nextflow run example.nf -dsl1
...
DataflowQueue(queue=[DataflowVariable(value=1), DataflowVariable(value=groovyx.gpars.dataflow.operator.PoisonPill@34be065a)])
```

Temos o valor 1 como único elemento do nosso canal de fila e uma pílula de veneno, que vai dizer ao processo que não há mais nada para ser consumido. É por isso que temos apenas uma saída para o exemplo acima, que é 2. Vamos inspecionar um canal de valor agora.

```groovy linenums="1"
canal1 = Channel.value(1)
println canal1
```

```console
$ nextflow run example.nf -dsl1
...
DataflowVariable(value=1)
```

Não há pílula de veneno, e é por isso que obtemos uma saída diferente com o código abaixo, onde `canal2` é transformado em um canal de valor por meio do operador `first`.

```groovy linenums="1"
canal1 = Channel.of(1,2,3)
canal2 = Channel.of(1)

process SUM {
    input:
    val x
    val y

    output:
    stdout

    script:
    """
    echo \$(($x+$y))
    """
}

workflow {
    SUM(canal1,canal2.first()).view()
}
```

```console title="Output"
4

3

2
```

Além disso, em muitas situações, o Nextflow converterá implicitamente variáveis em canais de valor quando forem usadas em uma chamada de processo. Por exemplo, quando você chama um processo com um parâmetro de pipeline (`params.example`) que possui um valor de string, ele é automaticamente convertido em um canal de valor.

## Fábricas de canal

Estes são comandos do Nextflow para criar canais que possuem entradas e funções implícitas esperadas.

### `value()`

A fábrica de canal `value` é utilizada para criar um canal de _valor_. Um argumento opcional não `nulo` pode ser especificado para vincular o canal a um valor específico. Por exemplo:

```groovy linenums="1"
canal1 = Channel.value() // (1)!
canal2 = Channel.value( 'Olá, você!' ) // (2)!
canal3 = Channel.value( [1,2,3,4,5] ) // (3)!
```

1. Cria um canal de valor _vazio_
2. Cria um canal de valor e vincula uma string a ele
3. Cria um canal de valor e vincula a ele um objeto de lista que será emitido como uma única emissão

### `of()`

A fábrica `Channel.of` permite a criação de um canal de fila com os valores especificados como argumentos.

```groovy linenums="1"
canal = Channel.of( 1, 3, 5, 7 )
canal.view{ "numero: $it" }
```

A primeira linha neste exemplo cria uma variável `canal` que contém um objeto de canal. Este canal emite os valores especificados como parâmetro no método `of`. Assim, a segunda linha imprimirá o seguinte:

```console
numero: 1
numero: 3
numero: 5
numero: 7
```

O método `Channel.of` funciona de maneira semelhante ao `Channel.from` (que foi [descontinuado](https://www.nextflow.io/docs/latest/channel.html#of)), corrigindo alguns comportamentos inconsistentes do último e fornecendo um melhor manuseio quando um intervalo de valores é especificado. Por exemplo, o seguinte funciona com um intervalo de 1 a 23:

```groovy linenums="1"
Channel
  .of(1..23, 'X', 'Y')
  .view()
```

### `fromList()`

O método `Channel.fromList` cria um canal emitindo os elementos fornecidos por um objeto de lista especificado como um argumento:

```groovy linenums="1"
list = ['olá', 'mundo']

Channel
  .fromList(list)
  .view()
```

### `fromPath()`

A fábrica `fromPath` cria um canal de fila emitindo um ou mais arquivos correspondentes ao padrão glob especificado.

```groovy linenums="1"
Channel.fromPath( './data/meta/*.csv' )
```

Este exemplo cria um canal e emite tantos itens quanto arquivos com extensão `csv` existirem na pasta `/data/meta`. Cada elemento é um objeto de arquivo implementando a interface [Path](https://docs.oracle.com/javase/8/docs/api/java/nio/file/Paths.html) do Java.

!!! tip

    Dois asteriscos, ou seja, `**`, funcionam como `*`, mas cruzam os limites do diretório. Essa sintaxe geralmente é usada para percorrer caminhos completos. Os colchetes especificam uma coleção de subpadrões.

| Nome          | Descrição                                                                                                                                                    |
| ------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| glob          | Quando `true` interpreta caracteres `*`, `?`, `[]` e `{}` como glob wildcards, caso contrário, os trata como caracteres normais (padrão: `true`)             |
| type          | Tipo de caminho retornado, ou `file`, `dir` ou `any` (padrão: `file`)                                                                                        |
| hidden        | Quando `true` inclui arquivos ocultos nos caminhos resultantes (padrão: `false`)                                                                             |
| maxDepth      | Número máximo de níveis de diretório a serem visitados (padrão: `no limit`)                                                                                  |
| followLinks   | Quando `true` links simbólicos são seguidos durante a travessia da árvore de diretórios, caso contrário, eles são gerenciados como arquivos (padrão: `true`) |
| relative      | Quando `true` os caminhos de retorno são relativos ao diretório de topo mais comum (padrão: `false`)                                                         |
| checkIfExists | Quando `true` lança uma exceção quando o caminho especificado não existe no sistema de arquivos (padrão: `false`)                                            |

Saiba mais sobre a sintaxe dos padrões glob [neste link](https://docs.oracle.com/javase/tutorial/essential/io/fileOps.html#glob).

!!! exercise

    Use o método `Channel.fromPath` para criar um canal emitindo todos os arquivos com o sufixo `.fq` no diretório `data/ggal/` e qualquer subdiretório, além dos arquivos ocultos. Em seguida, imprima os nomes dos arquivos.

    ??? solution

        ```groovy linenums="1"
        Channel.fromPath( './data/ggal/**.fq' , hidden:true)
          .view()
        ```

### `fromFilePairs()`

O método `fromFilePairs` cria um canal emitindo os pares de arquivos correspondentes a um padrão glob fornecido pelo usuário. Os arquivos correspondentes são emitidos como tuplas, nas quais o primeiro elemento é a chave de agrupamento do par correspondente e o segundo elemento é a lista de arquivos (classificados em ordem lexicográfica).

```groovy linenums="1"
Channel
  .fromFilePairs('./data/ggal/*_{1,2}.fq')
  .view()
```

Ele produzirá uma saída semelhante à seguinte:

```groovy
[liver, [/workspace/gitpod/nf-training/data/ggal/liver_1.fq, /workspace/gitpod/nf-training/data/ggal/liver_2.fq]]
[gut, [/workspace/gitpod/nf-training/data/ggal/gut_1.fq, /workspace/gitpod/nf-training/data/ggal/gut_2.fq]]
[lung, [/workspace/gitpod/nf-training/data/ggal/lung_1.fq, /workspace/gitpod/nf-training/data/ggal/lung_2.fq]]
```

!!! warning

    O padrão glob _precisa_ conter pelo menos um caractere curinga de estrela (`*`).

| Nome          | Descrição                                                                                                                                                    |
| ------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| type          | Tipo de caminhos retornados, ou `file`, `dir` ou `any` (padrão: `file`)                                                                                      |
| hidden        | Quando `true` includes hidden files in the resulting paths (padrão: `false`)                                                                                 |
| maxDepth      | Número máximo de níveis de diretório a serem visitados (padrão: `no limit`)                                                                                  |
| followLinks   | Quando `true` links simbólicos são seguidos durante a travessia da árvore de diretórios, caso contrário, eles são gerenciados como arquivos (padrão: `true`) |
| size          | Define o número de arquivos que cada item emitido deve conter (padrão: 2). Use `-1` para qualquer número                                                     |
| flat          | Quando `true` os arquivos correspondentes são produzidos como únicos elementos nas tuplas emitidas (padrão: `false`)                                         |
| checkIfExists | Quando `true`, lança uma exceção do caminho especificado que não existe no sistema de arquivos (padrão: `false`)                                             |

!!! exercise

    Use o método `fromFilePairs` para criar um canal emitindo todos os pares de leituras em fastq no diretório `data/ggal/` e imprima-os. Em seguida, use a opção `flat:true` e compare a saída com a execução anterior.

    ??? solution

        Use o seguinte, com ou sem `flat:true`:

        ```groovy linenums="1"
        Channel.fromFilePairs( './data/ggal/*_{1,2}.fq', flat:true)
          .view()
        ```

        Em seguida, verifique os colchetes ao redor dos nomes dos arquivos, para ver a diferença com `flat`.

### `fromSRA()`

O método `Channel.fromSRA` permite consultar o banco de dados [NCBI SRA](https://www.ncbi.nlm.nih.gov/sra) e retorna um canal que emite os arquivos FASTQ correspondentes aos critérios de seleção especificados.

A consulta pode ser ID(s) de projeto(s) ou número(s) de acesso suportado(s) pela API do [NCBI ESearch](https://www.ncbi.nlm.nih.gov/books/NBK25499/#chapter4.ESearch).

!!! info

    Esta função agora requer uma chave de API que você só pode obter fazendo login em sua conta NCBI.

??? example "Instruções para login do NCBI e aquisição de chave"

    1. Vá para: <https://www.ncbi.nlm.nih.gov/>
    2. Clique no botão "Login" no canto superior direito para entrar no NCBI. Siga suas instruções.
    3. Uma vez em sua conta, clique no botão no canto superior direito, geralmente seu ID.
    4. Vá para Account settings
    5. Role para baixo até a seção API Key Management.
    6. Clique em "Criar uma chave de API".
    7. A página será atualizada e a chave será exibida onde estava o botão. Copie sua chave.

Por exemplo, o trecho a seguir imprimirá o conteúdo de um ID de projeto NCBI:

```groovy linenums="1"
params.ncbi_api_key = '<Sua chave da API aqui>'

Channel
  .fromSRA(['SRP073307'], apiKey: params.ncbi_api_key)
  .view()
```

!!! info ""

    :material-lightbulb: Substitua `<Your API key here>` com sua chave de API.

Isso deve imprimir:

```groovy
[SRR3383346, [/vol1/fastq/SRR338/006/SRR3383346/SRR3383346_1.fastq.gz, /vol1/fastq/SRR338/006/SRR3383346/SRR3383346_2.fastq.gz]]
[SRR3383347, [/vol1/fastq/SRR338/007/SRR3383347/SRR3383347_1.fastq.gz, /vol1/fastq/SRR338/007/SRR3383347/SRR3383347_2.fastq.gz]]
[SRR3383344, [/vol1/fastq/SRR338/004/SRR3383344/SRR3383344_1.fastq.gz, /vol1/fastq/SRR338/004/SRR3383344/SRR3383344_2.fastq.gz]]
[SRR3383345, [/vol1/fastq/SRR338/005/SRR3383345/SRR3383345_1.fastq.gz, /vol1/fastq/SRR338/005/SRR3383345/SRR3383345_2.fastq.gz]]
// (o restante foi omitido)
```

Vários IDs de acesso podem ser especificados usando um objeto lista:

```groovy linenums="1"
ids = ['ERR908507', 'ERR908506', 'ERR908505']
Channel
  .fromSRA(ids, apiKey: params.ncbi_api_key)
  .view()
```

```groovy
[ERR908507, [/vol1/fastq/ERR908/ERR908507/ERR908507_1.fastq.gz, /vol1/fastq/ERR908/ERR908507/ERR908507_2.fastq.gz]]
[ERR908506, [/vol1/fastq/ERR908/ERR908506/ERR908506_1.fastq.gz, /vol1/fastq/ERR908/ERR908506/ERR908506_2.fastq.gz]]
[ERR908505, [/vol1/fastq/ERR908/ERR908505/ERR908505_1.fastq.gz, /vol1/fastq/ERR908/ERR908505/ERR908505_2.fastq.gz]]
```

!!! info

    Os pares de leitura são gerenciados implicitamente e são retornados como uma lista de arquivos.

É fácil usar este canal como uma entrada usando a sintaxe usual do Nextflow. O código abaixo cria um canal contendo
duas amostras de um estudo SRA público e executa o FASTQC nos arquivos resultantes. Veja:

```groovy linenums="1"
params.ncbi_api_key = '<Your API key here>'

params.accession = ['ERR908507', 'ERR908506']

process fastqc {
  input:
  tuple val(sample_id), path(reads_file)

  output:
  path("fastqc_${sample_id}_logs")

  script:
  """
  mkdir fastqc_${sample_id}_logs
  fastqc -o fastqc_${sample_id}_logs -f fastq -q ${reads_file}
  """
}

workflow {
  reads = Channel.fromSRA(params.accession, apiKey: params.ncbi_api_key)
  fastqc(reads)
}
```

Se você deseja executar o pipeline acima e não possui o fastqc instalado em sua máquina, não esqueça o que aprendeu na seção anterior. Execute este pipeline com `-with-docker biocontainers/fastqc:v0.11.5`, por exemplo.

### Arquivos de texto

O operador `splitText` permite dividir strings de várias linhas ou itens de arquivo de texto, emitidos por um canal de origem em blocos contendo n linhas, que serão emitidos pelo canal resultante. Veja:

```groovy linenums="1"
Channel
  .fromPath('data/meta/random.txt') // (1)!
  .splitText() // (2)!
  .view() // (3)!
```

1. Instrui o Nextflow a criar um canal a partir do caminho `data/meta/random.txt`
2. O operador `splitText` divide cada item em pedaços de uma linha por padrão.
3. Veja o conteúdo do canal.

Você pode definir o número de linhas em cada bloco usando o parâmetro `by`, conforme mostrado no exemplo a seguir:

```groovy linenums="1"
Channel
  .fromPath('data/meta/random.txt')
  .splitText( by: 2 )
  .subscribe {
    print it;
    print "--- fim do bloco ---\n"
  }
```

!!! info

    O operador `subscribe` permite a execução de funções definidas pelo usuário cada vez que um novo valor é emitido pelo canal de origem.

Uma clausura opcional pode ser especificada para transformar os blocos de texto produzidos pelo operador. O exemplo a seguir mostra como dividir arquivos de texto em blocos de 10 linhas e transformá-los em letras maiúsculas:

```groovy linenums="1"
Channel
  .fromPath('data/meta/random.txt')
  .splitText( by: 10 ) { it.toUpperCase() }
  .view()
```

Você também pode fazer contagens para cada linha:

```groovy linenums="1"
count=0

Channel
  .fromPath('data/meta/random.txt')
  .splitText()
  .view { "${count++}: ${it.toUpperCase().trim()}" }
```

Por fim, você também pode usar o operador em arquivos simples (fora do contexto do canal):

```groovy linenums="1"
def f = file('data/meta/random.txt')
def lines = f.splitText()
def count=0
for( String row : lines ) {
  log.info "${count++} ${row.toUpperCase()}"
}
```

### Valores separados por vírgula (.csv)

O operador `splitCsv` permite analisar itens de texto formatados em CSV (Comma-separated value) emitidos por um canal.

Em seguida, ele os divide em registros ou os agrupa como uma lista de registros com um comprimento especificado.

No caso mais simples, basta aplicar o operador `splitCsv` a um canal que emite arquivos de texto ou entradas de texto no formato CSV. Por exemplo, para visualizar apenas a primeira e a quarta colunas:

```groovy linenums="1"
Channel
  .fromPath("data/meta/patients_1.csv")
  .splitCsv()
  // linha é um objeto de lista
  .view { linha -> "${linha[0]},${linha[3]}" }
```

Quando o CSV começa com uma linha de cabeçalho definindo os nomes das colunas, você pode especificar o parâmetro `header: true` que permite referenciar cada valor pelo nome da coluna, conforme mostrado no exemplo a seguir:

```groovy linenums="1"
Channel
  .fromPath("data/meta/patients_1.csv")
  .splitCsv(header: true)
  // linha é um objeto de lista
  .view { linha -> "${linha.patient_id},${linha.num_samples}" }
```

Como alternativa, você pode fornecer nomes de cabeçalho personalizados especificando uma lista de strings no parâmetro de cabeçalho, conforme mostrado abaixo:

```groovy linenums="1"
Channel
  .fromPath("data/meta/patients_1.csv")
  .splitCsv(header: ['col1', 'col2', 'col3', 'col4', 'col5'] )
  // linha é um objeto de lista
  .view { linha -> "${linha.col1},${linha.col4}" }
```

Você também pode processar vários arquivos CSV ao mesmo tempo:

```groovy linenums="1"
Channel
  .fromPath("data/meta/patients_*.csv") // <-- use um padrão de captura
  .splitCsv(header:true)
  .view { row -> "${row.patient_id}\t${row.num_samples}" }
```

!!! tip

    Observe que você pode alterar o formato de saída simplesmente adicionando um delimitador diferente.

Por fim, você também pode operar em arquivos CSV fora do contexto do canal:

```groovy linenums="1"
def f = file('data/meta/patients_1.csv')
def lines = f.splitCsv()
for( List row : lines ) {
  log.info "${row[0]} -- ${row[2]}"
}
```

!!! exercise

    Tente inserir leituras fastq no fluxo de trabalho do RNA-Seq anterior usando `.splitCsv`.

    ??? solution

        Adicione um arquivo de texto CSV contendo o seguinte, como uma entrada de exemplo com o nome "fastq.csv":

        ```csv
        gut,/workspace/gitpod/nf-training/data/ggal/gut_1.fq,/workspace/gitpod/nf-training/data/ggal/gut_2.fq
        ```

        Em seguida, substitua o canal de entrada para as leituras em `script7.nf`, alterando as seguintes linhas:

        ```groovy linenums="1"
        Channel
          .fromFilePairs( params.reads, checkIfExists: true )
          .set { read_pairs_ch }
        ```

        Para uma entrada de fábrica de canal splitCsv:

        ```groovy linenums="1" hl_lines="2 3 4"
        Channel
          .fromPath("fastq.csv")
          .splitCsv()
          .view () { row -> "${row[0]},${row[1]},${row[2]}" }
          .set { read_pairs_ch }
        ```

        Por fim, altere a cardinalidade dos processos que usam os dados de entrada. Por exemplo, para o processo de quantificação, mude de:

        ```groovy linenums="1"
        process QUANTIFICATION {
          tag "$sample_id"

          input:
          path salmon_index
          tuple val(sample_id), path(reads)

          output:
          path sample_id, emit: quant_ch

          script:
          """
          salmon quant --threads $task.cpus --libType=U -i $salmon_index -1 ${reads[0]} -2 ${reads[1]} -o $sample_id
          """
        }
        ```

        Para:

        ```groovy linenums="1" hl_lines="6 13"
        process QUANTIFICATION {
          tag "$sample_id"

          input:
          path salmon_index
          tuple val(sample_id), path(reads1), path(reads2)

          output:
          path sample_id, emit: quant_ch

          script:
          """
          salmon quant --threads $task.cpus --libType=U -i $salmon_index -1 ${reads1} -2 ${reads2} -o $sample_id
          """
        }
        ```

        Repita o procedimento acima para a etapa fastqc.

        ```groovy linenums="1"
        process FASTQC {
          tag "FASTQC on $sample_id"

          input:
          tuple val(sample_id), path(reads1), path(reads2)

          output:
          path "fastqc_${sample_id}_logs"

          script:
          """
          mkdir fastqc_${sample_id}_logs
          fastqc -o fastqc_${sample_id}_logs -f fastq -q ${reads1} ${reads2}
          """
        }
        ```

        Agora o fluxo de trabalho deve ser executado a partir de um arquivo CSV.

### Valores separados por tabulação (.tsv)

A análise de arquivos tsv funciona de maneira semelhante, basta adicionar a opção `sep:'\t'` no contexto `splitCsv`:

```groovy linenums="1"
Channel
  .fromPath("data/meta/regions.tsv", checkIfExists:true)
  // Use a opção `sep` para analisar arquivos com tabulação como separador
  .splitCsv(sep:'\t')
  // linha é um objeto de lista
  .view()
```

!!! exercise

    Tente usar a técnica de separação por tabulação no arquivo `data/meta/regions.tsv`, mas imprima apenas a primeira coluna e remova o cabeçalho.


    ??? solution

        ```groovy linenums="1"
        Channel
          .fromPath("data/meta/regions.tsv", checkIfExists:true)
          // Use a opção `sep` para analisar arquivos com tabulação como separador
          .splitCsv(sep:'\t', header:true )
          // linha é um objeto de lista
          .view { linha -> "${row.patient_id}" }
        ```

## Formatos de arquivo mais complexos

### JSON

Também podemos analisar facilmente o formato de arquivo JSON usando o seguinte esquema do Groovy:

```groovy linenums="1"
import groovy.json.JsonSlurper

def f = file('data/meta/regions.json')
def records = new JsonSlurper().parse(f)


for( def entry : records ) {
  log.info "$entry.patient_id -- $entry.feature"
}
```

!!! warning

    Ao usar uma versão JSON mais antiga, pode ser necessário substituir `parse(f)` por `parseText(f.text)`

### YAML

Isso também pode ser usado como uma forma de analisar arquivos YAML:

```groovy linenums="1"
import org.yaml.snakeyaml.Yaml

def f = file('data/meta/regions.yml')
def records = new Yaml().load(f)


for( def entry : records ) {
  log.info "$entry.patient_id -- $entry.feature"
}
```

### Armazenamento em módulos de analisadores sintáticos

A melhor maneira de armazenar scripts com analisadores é mantê-los em um arquivo de módulo Nextflow.

Veja o seguinte script Nextflow:

```groovy linenums="1"
include{ parseJsonFile } from './modules/parsers.nf'

process foo {
  input:
  tuple val(meta), path(data_file)

  """
  echo seu_comando $meta.region_id $data_file
  """
}

workflow {
  Channel.fromPath('data/meta/regions*.json') \
    | flatMap { parseJsonFile(it) } \
    | map { entry -> tuple(entry,"/some/data/${entry.patient_id}.txt") } \
    | foo
}
```

Para que este script funcione, um arquivo de módulo chamado `parsers.nf` precisa ser criado e armazenado em uma pasta de módulos no diretório atual.

O arquivo `parsers.nf` deve conter a função `parseJsonFile`.

O Nextflow usará isso como uma função personalizada dentro do escopo do fluxo de trabalho.

!!! tip

    Você aprenderá mais sobre arquivos de módulo posteriormente na [seção de Modularização](modules/) desse tutorial.
