---
description: Material de treinamento básico do Nextflow
---

# Processos

No Nextflow, um processo (`process`) é a primitiva de computação básica para executar funções estrangeiras (ou seja, scripts personalizados ou ferramentas).

A definição do processo começa com a palavra-chave `process`, seguida pelo nome do processo e, finalmente, o corpo do processo delimitado por chaves.

O nome do process é comumente escrito em letras maiúsculas por convenção.

Um processo básico, usando apenas o bloco de definição `script`, se parece com o seguinte:

```groovy linenums="1"
process DIGAOLA {

  script:
  """
  echo 'Olá mundo!'
  """
}
```

Em exemplos mais complexos, o corpo do processo pode conter até **cinco** blocos de definição:

1. **Diretivas** são declarações iniciais que definem configurações opcionais
2. **Input** (Bloco de entrada) define o(s) arquivo(s) de entrada esperado(s) e o canal onde encontrá-los
3. **Output** (Bloco de saída) define o(s) arquivo(s) de saída esperado(s) e o canal para enviar os dados
4. **When** é uma declaração de cláusula opcional para permitir processos condicionais
5. **Script** é uma string que define o comando a ser executado pelo processo

A sintaxe completa do processo é definida da seguinte forma:

!!! info ""

    Clique no ícone :material-plus-circle: no código para ver explicações.

```groovy linenums="1"
process < nome > {

  [ diretivas ] // (1)!

  input: // (2)!
  < entradas do processo >

  output: // (3)!
  < saídas do processo >

  when: // (4)!
  < condição >

  [script|shell|exec]: // (5)!
  """
  < script do usuário a ser executado >
  """
}
```

1. Zero, uma ou mais diretivas de processo
2. Zero, uma ou mais entradas para o processo
3. Zero, uma ou mais saídas para o processo
4. Uma condicional booleana opcional para acionar a execução do processo
5. O comando a ser executado

## Script

O bloco `script` é uma string que define o comando a ser executado pelo processo.

Um processo pode executar apenas um bloco `script`. Deve ser a última instrução quando o processo contém declarações de entrada e saída.

O bloco `script` pode ser uma string de uma ou várias linhas. A de várias linhas simplifica a escrita de scripts não triviais compostos por vários comandos abrangendo várias linhas. Por exemplo:

```groovy linenums="1"
process EXEMPLO {

  script:
  """
  echo 'Olá mundo!\nHola mundo!\nCiao mondo!\nHallo Welt!' > arquivo
  cat arquivo | head -n 1 | head -c 5 > pedaco_1.txt
  gzip -c pedaco_1.txt  > pedacos.gz
  """
}

workflow {
  EXEMPLO()
}
```

Por padrão, o comando `process` é interpretado como um script **Bash**. No entanto, qualquer outra linguagem de script pode ser usada simplesmente iniciando o script com a declaração [Shebang](<https://en.wikipedia.org/wiki/Shebang_(Unix)>) adequada. Por examplo:

```groovy linenums="1"
process CODIGOPYTHON {

  script:
  """
  #!/usr/bin/env python

  x = 'Olá'
  y = 'mundo!'
  print ("%s - %s" % (x,y))
  """
}

workflow {
  CODIGOPYTHON()
}
```

!!! tip

    Várias linguagens de programação podem ser usadas no mesmo script de fluxo de trabalho. No entanto, para grandes blocos de código, é melhor salvá-los em arquivos separados e invocá-los a partir do script do processo. Pode-se armazenar os scripts específicos na pasta `./bin/`.

### Parâmetros do script

Parâmetros de script (`params`) podem ser definidos dinamicamente usando valores variáveis. Por exemplo:

```groovy linenums="1"
params.data = 'Mundo'

process FOO {

  script:
  """
  echo Olá $params.data
  """
}

workflow {
  FOO()
}
```

!!! info

    Um script de processo pode conter qualquer formato de string suportado pela linguagem de programação Groovy. Isso nos permite usar a interpolação de strings como no script acima ou strings multilinha. Consulte [Interpolação de string](#groovy.adoc#_string_interpolation) para obter mais informações.

!!! warning

    Como o Nextflow usa a mesma sintaxe Bash para substituições de variáveis em strings, as variáveis de ambiente Bash precisam ser escapadas usando o caractere `\`.

```groovy linenums="1"
process FOO {

  script:
  """
  echo "O diretório atual é \$PWD"
  """
}

workflow {
  FOO()
}
```

Pode ser complicado escrever um script que usa muitas variáveis Bash. Uma alternativa possível é usar uma string de script delimitada por aspas simples

```groovy linenums="1"
process BAR {

  script:
  """
  echo $PATH | tr : '\\n'
  """
}

workflow {
  BAR()
}
```

No entanto, isso bloqueia o uso de variáveis Nextflow no script de comando.

Outra alternativa é usar uma instrução `shell` em vez de `script` e usar uma sintaxe diferente para variáveis do Nextflow, por exemplo, `!{..}`. Isso permite o uso das variáveis Nextflow e Bash no mesmo script.

```groovy linenums="1"
params.data = 'le monde'

process BAZ {

  shell:
  '''
  X='Bonjour'
  echo $X !{params.data}
  '''
}

workflow {
  BAZ()
}
```

### Scripts condicionais

O script do processo também pode ser definido de maneira completamente dinâmica usando uma instrução `if` ou qualquer outra expressão para avaliar um valor de string. Por exemplo:

```groovy linenums="1"
params.compressao = 'gzip'
params.arquivo_a_comprimir = "$baseDir/data/ggal/transcriptome.fa"

process FOO {

  input:
  path arquivo

  script:
  if( params.compressao == 'gzip' )
    """
    gzip -c $arquivo > ${arquivo}.gz
    """
  else if( params.compressao == 'bzip2' )
    """
    bzip2 -c $arquivo > ${arquivo}.bz2
    """
  else
    throw new IllegalArgumentException("Alinhador $params.compressao desconhecido")
}

workflow {
  FOO(params.arquivo_a_comprimir)
}
```

## Canais de entradas

Os processos Nextflow são isolados uns dos outros, mas podem se comunicar entre si enviando valores por meio de canais.

As entradas determinam implicitamente as dependências e a execução paralela do processo. A execução do processo é disparada cada vez que dados _novos_ estão prontos para serem consumidos do canal de entrada:

<figure class="excalidraw">
--8<-- "docs/basic_training/img/channel-process.excalidraw.svg"
</figure>

O bloco `input` define de quais canais o processo espera receber dados. Você só pode definir um bloco `input` por vez e deve conter uma ou mais declarações de entrada.

O bloco `input` segue a sintaxe mostrada abaixo:

```groovy linenums="1"
input:
  <qualificador da variável entrada> <nome da variável entrada>
```

### Valores de entrada

O qualificador `val` permite receber dados de qualquer tipo como entrada. Ele pode ser acessado no script do processo usando o nome de entrada especificado, conforme mostrado no exemplo a seguir:

```groovy linenums="1"
num = Channel.of( 1, 2, 3 )

process EXEMPLOBASICO {
  debug true

  input:
  val x

  script:
  """
  echo tarefa $x do processo
  """
}

workflow {
  myrun = EXEMPLOBASICO(num)
}
```

No exemplo acima, o processo é executado três vezes, cada vez que um valor é recebido do canal `num` e usado para processar o script. Assim, resulta em uma saída semelhante à mostrada abaixo:

```console
tarefa 3 do processo
tarefa 1 do processo
tarefa 2 do processo
```

!!! warning

    O canal garante que os itens sejam entregues na mesma ordem em que foram enviados - mas - como o processo é executado de forma paralela, não há garantia de que sejam processados na mesma ordem em que foram recebidos.

### Arquivo e caminhos de entrada

The `path` qualifier allows the handling of file values in the process execution context. This means that Nextflow will stage it in the process execution directory, and it can be accessed in the script by using the name specified in the input declaration.

```groovy linenums="1"
reads = Channel.fromPath( 'data/ggal/*.fq' )

process FOO {
  debug true

  input:
  path 'sample.fastq'

  script:
  """
  ls sample.fastq
  """
}

workflow {
  result = FOO(reads)
}
```

The input file name can also be defined using a variable reference as shown below:

```groovy linenums="1"
reads = Channel.fromPath( 'data/ggal/*.fq' )

process FOO {
  debug true

  input:
  path sample

  script:
  """
  ls  $sample
  """
}

workflow {
  result = FOO(reads)
}
```

The same syntax is also able to handle more than one input file in the same execution and only requires changing the channel composition.

```groovy linenums="1"
reads = Channel.fromPath( 'data/ggal/*.fq' )

process FOO {
  debug true

  input:
  path sample

  script:
  """
  ls -lh $sample
  """
}

workflow {
  FOO(reads.collect())
}
```

!!! warning

    In the past, the `file` qualifier was used for files, but the `path` qualifier should be preferred over file to handle process input files when using Nextflow 19.10.0 or later. When a process declares an input file, the corresponding channel elements must be **file** objects created with the file helper function from the file specific channel factories (e.g., `Channel.fromPath` or `Channel.fromFilePairs`).

!!! exercise

    Write a script that creates a channel containing all read files matching the pattern `data/ggal/*_1.fq` followed by a process that concatenates them into a single file and prints the first 20 lines.


    ??? solution

        ```groovy linenums="1"
        params.reads = "$baseDir/data/ggal/*_1.fq"

        Channel
          .fromPath( params.reads )
          .set { read_ch }

        process CONCATENATE {
          tag "Concat all files"

          input:
          path '*'

          output:
          path 'top_10_lines'

          script:
          """
          cat * > concatenated.txt
          head -n 20 concatenated.txt > top_10_lines
          """
        }

        workflow {
          concat_ch = CONCATENATE(read_ch.collect())
          concat_ch.view()
        }
        ```

### Combinando canais de entrada

A key feature of processes is the ability to handle inputs from multiple channels. However, it’s important to understand how channel contents and their semantics affect the execution of a process.

Consider the following example:

```groovy linenums="1"
ch1 = Channel.of(1,2,3)
ch2 = Channel.of('a','b','c')

process FOO {
  debug true

  input:
  val x
  val y

  script:
    """
    echo $x and $y
    """
}

workflow {
  FOO(ch1, ch2)
}
```

Both channels emit three values, therefore the process is executed three times, each time with a different pair:

-   `(1, a)`
-   `(2, b)`
-   `(3, c)`

What is happening is that the process waits until there’s a complete input configuration, i.e., it receives an input value from all the channels declared as input.

When this condition is verified, it consumes the input values coming from the respective channels, spawns a task execution, then repeats the same logic until one or more channels have no more content.

This means channel values are consumed serially one after another and the first empty channel causes the process execution to stop, even if there are other values in other channels.

**So what happens when channels do not have the same cardinality (i.e., they emit a different number of elements)?**

For example:

```groovy linenums="1"
input1 = Channel.of(1,2)
input2 = Channel.of('a','b','c','d')

process FOO {
  debug true

  input:
  val x
  val y

  script:
    """
    echo $x and $y
    """
}

workflow {
  FOO(input1, input2)
}
```

In the above example, the process is only executed twice because the process stops when a channel has no more data to be processed.

However, what happens if you replace value x with a `value` channel?

Compare the previous example with the following one :

```groovy linenums="1"
input1 = Channel.value(1)
input2 = Channel.of('a','b','c')

process BAR {
  debug true

  input:
  val x
  val y

  script:
    """
    echo $x and $y
    """
}

workflow {
  BAR(input1, input2)
}
```

```console title="Script output"
1 and b
1 and a
1 and c
```

This is because _value_ channels can be consumed multiple times and do not affect process termination.

!!! exercise

    Write a process that is executed for each read file matching the pattern `data/ggal/*_1.fq` and use the same `data/ggal/transcriptome.fa` in each execution.

    ??? solution

        ```groovy linenums="1"
        params.reads = "$baseDir/data/ggal/*_1.fq"
        params.transcriptome_file = "$baseDir/data/ggal/transcriptome.fa"

        Channel
            .fromPath( params.reads )
            .set { read_ch }

        process COMMAND {
          tag "Run_command"

          input:
          path reads
          path transcriptome

          output:
          path result

          script:
          """
          echo your_command $reads $transcriptome > result
          """
        }

        workflow {
          concat_ch = COMMAND(read_ch, params.transcriptome_file)
          concat_ch.view()
        }
        ```

### Repetidores de entradas

The `each` qualifier allows you to repeat the execution of a process for each item in a collection every time new data is received. For example:

```groovy linenums="1"
sequences = Channel.fromPath('data/prots/*.tfa')
methods = ['regular', 'espresso', 'psicoffee']

process ALIGNSEQUENCES {
  debug true

  input:
  path seq
  each mode

  script:
  """
  echo t_coffee -in $seq -mode $mode
  """
}

workflow {
  ALIGNSEQUENCES(sequences, methods)
}
```

In the above example, every time a file of sequences is received as an input by the process, it executes three tasks, each running a different alignment method set as a `mode` variable. This is useful when you need to repeat the same task for a given set of parameters.

!!! exercise

    Extend the previous example so a task is executed for each read file matching the pattern `data/ggal/*_1.fq` and repeat the same task with both `salmon` and `kallisto`.

    ??? solution

        ```groovy linenums="1"
        params.reads = "$baseDir/data/ggal/*_1.fq"
        params.transcriptome_file = "$baseDir/data/ggal/transcriptome.fa"
        methods= ['salmon', 'kallisto']

        Channel
            .fromPath( params.reads )
            .set { read_ch }

        process COMMAND {
          tag "Run_command"

          input:
          path reads
          path transcriptome
          each mode

          output:
          path result

          script:
          """
          echo $mode $reads $transcriptome > result
          """
        }

        workflow {
          concat_ch = COMMAND(read_ch , params.transcriptome_file, methods)
          concat_ch
              .view { "To run : ${it.text}" }
        }
        ```

## Canais de saída

The _output_ declaration block defines the channels used by the process to send out the results produced.

Only one output block, that can contain one or more output declaration, can be defined. The output block follows the syntax shown below:

```groovy linenums="1"
output:
  <output qualifier> <output name> , emit: <output channel>
```

### Valores de saída

The `val` qualifier specifies a defined _value_ in the script context. Values are frequently defined in the _input_ and/or _output_ declaration blocks, as shown in the following example:

```groovy linenums="1"
methods = ['prot','dna', 'rna']

process FOO {

  input:
  val x

  output:
  val x

  script:
  """
  echo $x > file
  """
}

workflow {
  receiver_ch = FOO(Channel.of(methods))
  receiver_ch.view { "Received: $it" }
}
```

### Caminhos e arquivos de saída

The `path` qualifier specifies one or more files produced by the process into the specified channel as an output.

```groovy linenums="1"
process RANDOMNUM {

    output:
    path 'result.txt'

    script:
    """
    echo $RANDOM > result.txt
    """
}


workflow {
  receiver_ch = RANDOMNUM()
  receiver_ch.view { "Received: " + it.text }
}
```

In the above example the process `RANDOMNUM` creates a file named `result.txt` containing a random number.

Since a file parameter using the same name is declared in the output block, the file is sent over the `receiver_ch` channel when the task is complete. A downstream `process` declaring the same channel as _input_ will be able to receive it.

### Múltiplos arquivos de saída

When an output file name contains a wildcard character (`*` or `?`) it is interpreted as a [glob](http://docs.oracle.com/javase/tutorial/essential/io/fileOps.html#glob) path matcher. This allows us to _capture_ multiple files into a list object and output them as a sole emission. For example:

```groovy linenums="1"
process SPLITLETTERS {

    output:
    path 'chunk_*'

    """
    printf 'Hola' | split -b 1 - chunk_
    """
}

workflow {
    letters = SPLITLETTERS()
    letters
        .flatMap()
        .view { "File: ${it.name} => ${it.text}" }
}
```

Prints the following:

```console
File: chunk_aa => H
File: chunk_ab => o
File: chunk_ac => l
File: chunk_ad => a
```

Some caveats on glob pattern behavior:

-   Input files are not included in the list of possible matches
-   Glob pattern matches both files and directory paths
-   When a two stars pattern `**` is used to recourse across directories, only file paths are matched i.e., directories are not included in the result list.

!!! exercise

    Remove the `flatMap` operator and see out the output change. The documentation for the `flatMap` operator is available at [this link](https://www.nextflow.io/docs/latest/operator.html#flatmap).

    ??? result

        ```groovy
        File: [chunk_aa, chunk_ab, chunk_ac, chunk_ad] => [H, o, l, a]
        ```

### Nomes dinâmicos de arquivos de saída

When an output file name needs to be expressed dynamically, it is possible to define it using a dynamic string that references values defined in the input declaration block or in the script global context. For example:

```groovy linenums="1"
species = ['cat','dog', 'sloth']
sequences = ['AGATAG','ATGCTCT', 'ATCCCAA']

Channel.fromList(species)
        .set { species_ch }

process ALIGN {

  input:
  val x
  val seq

  output:
  path "${x}.aln"

  script:
  """
  echo align -in $seq > ${x}.aln
  """
}

workflow {
  genomes = ALIGN( species_ch, sequences )
  genomes.view()
}
```

In the above example, each time the process is executed an alignment file is produced whose name depends on the actual value of the `x` input.

### Entradas e saídas compostas

So far we have seen how to declare multiple input and output channels that can handle one value at a time. However, Nextflow can also handle a _tuple_ of values.

The input and output declarations for tuples must be declared with a `tuple` qualifier followed by the definition of each element in the tuple.

```groovy linenums="1"
reads_ch = Channel.fromFilePairs('data/ggal/*_{1,2}.fq')

process FOO {

  input:
    tuple val(sample_id), path(sample_id)

  output:
    tuple val(sample_id), path('sample.bam')

  script:
  """
    echo your_command_here --reads $sample_id > sample.bam
  """
}

workflow {
  bam_ch = FOO(reads_ch)
  bam_ch.view()
}
```

!!! info

    In previous versions of Nextflow `tuple` was called `set` but it was used the same way with the same semantic.

!!! exercise

    Modify the script of the previous exercise so that the _bam_ file is named as the given `sample_id`.

    ??? solution

        ```groovy linenums="1"
        reads_ch = Channel.fromFilePairs('data/ggal/*_{1,2}.fq')

        process FOO {

          input:
            tuple val(sample_id), path(sample_files)

          output:
            tuple val(sample_id), path("${sample_id}.bam")

          script:
          """
            echo your_command_here --reads $sample_id > ${sample_id}.bam
          """
        }

        workflow {
          bam_ch = FOO(reads_ch)
          bam_ch.view()
        }
        ```

## Quando

The `when` declaration allows you to define a condition that must be verified in order to execute the process. This can be any expression that evaluates a boolean value.

It is useful to enable/disable the process execution depending on the state of various inputs and parameters. For example:

```groovy linenums="1"
params.dbtype = 'nr'
params.prot = 'data/prots/*.tfa'
proteins = Channel.fromPath(params.prot)

process FIND {
  debug true

  input:
  path fasta
  val type

  when:
  fasta.name =~ /^BB11.*/ && type == 'nr'

  script:
  """
  echo blastp -query $fasta -db nr
  """
}

workflow {
  result = FIND(proteins, params.dbtype)
}
```

## Diretivas

Directive declarations allow the definition of optional settings that affect the execution of the current process without affecting the _semantic_ of the task itself.

They must be entered at the top of the process body, before any other declaration blocks (i.e., `input`, `output`, etc.).

Directives are commonly used to define the amount of computing resources to be used or other meta directives that allow the definition of extra configuration of logging information. For example:

```groovy linenums="1"
process FOO {
  cpus 2
  memory 1.GB
  container 'image/name'

  script:
  """
  echo your_command --this --that
  """
}
```

!!! info ""

    :material-lightbulb: The complete list of directives is available [at this link](https://www.nextflow.io/docs/latest/process.html#directives).

| Name                                                                | Description                                                                                                                                          |
| ------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------- |
| [`cpus`](https://www.nextflow.io/docs/latest/process.html#cpus)     | Allows you to define the number of (logical) CPUs required by the process’ task.                                                                     |
| [`time`](https://www.nextflow.io/docs/latest/process.html#time)     | Allows you to define how long a process is allowed to run (e.g., time _1h_: 1 hour, _1s_ 1 second, _1m_ 1 minute, _1d_ 1 day).                       |
| [`memory`](https://www.nextflow.io/docs/latest/process.html#memory) | Allows you to define how much memory the process is allowed to use (e.g., _2 GB_ is 2 GB). Can also use B, KB,MB,GB and TB.                          |
| [`disk`](https://www.nextflow.io/docs/latest/process.html#disk)     | Allows you to define how much local disk storage the process is allowed to use.                                                                      |
| [`tag`](https://www.nextflow.io/docs/latest/process.html#tag)       | Allows you to associate each process execution with a custom label to make it easier to identify them in the log file or the trace execution report. |

## Organizando as saídas

### A diretiva PublishDir

Given each process is being executed in separate temporary `work/` folder (e.g., `work/f1/850698…`; `work/g3/239712…`; etc.), we may want to save important, non-intermediary, and/or final files in a results folder.

!!! tip

    Remember to delete the work folder from time to time to clear your intermediate files and stop them from filling your computer!

To store our workflow result files, we need to explicitly mark them using the directive [publishDir](https://www.nextflow.io/docs/latest/process.html#publishdir) in the process that’s creating the files. For example:

```groovy linenums="1"
params.outdir = 'my-results'
params.prot = 'data/prots/*.tfa'
proteins = Channel.fromPath(params.prot)


process BLASTSEQ {
    publishDir "$params.outdir/bam_files", mode: 'copy'

    input:
    path fasta

    output:
    path ('*.txt')

    script:
    """
    echo blastp $fasta > ${fasta}_result.txt
    """
}

workflow {
  blast_ch = BLASTSEQ(proteins)
  blast_ch.view()
}
```

The above example will copy all blast script files created by the `BLASTSEQ` task into the directory path `my-results`.

!!! tip

    The publish directory can be local or remote. For example, output files could be stored using an [AWS S3 bucket](https://aws.amazon.com/s3/) by using the `s3://` prefix in the target path.

### Gerenciar semântica de subdiretórios

You can use more than one `publishDir` to keep different outputs in separate directories. For example:

```groovy linenums="1"
params.reads = 'data/reads/*_{1,2}.fq.gz'
params.outdir = 'my-results'

samples_ch = Channel.fromFilePairs(params.reads, flat: true)

process FOO {
  publishDir "$params.outdir/$sampleId/", pattern: '*.fq'
  publishDir "$params.outdir/$sampleId/counts", pattern: "*_counts.txt"
  publishDir "$params.outdir/$sampleId/outlooks", pattern: '*_outlook.txt'

  input:
    tuple val(sampleId), path('sample1.fq.gz'), path('sample2.fq.gz')

  output:
    path "*"

  script:
  """
    < sample1.fq.gz zcat > sample1.fq
    < sample2.fq.gz zcat > sample2.fq

    awk '{s++}END{print s/4}' sample1.fq > sample1_counts.txt
    awk '{s++}END{print s/4}' sample2.fq > sample2_counts.txt

    head -n 50 sample1.fq > sample1_outlook.txt
    head -n 50 sample2.fq > sample2_outlook.txt
  """
}

workflow {
  out_channel = FOO(samples_ch)
}
```

The above example will create an output structure in the directory `my-results`, that contains a separate sub-directory for each given sample ID, each containing the folders `counts` and `outlooks`.
