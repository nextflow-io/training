---
title: Cache and resume
description: Basic Nextflow Training Workshop
---

# Execução de cache e resume

O mecanismo de caching do Nextflow funciona atribuindo uma ID única para cada tarefa que é usada para criar um diretório de execução separado onde as tarefas são executadas e os resultados guardados.

A ID de tarefa única é gerada como um número hash de 128-bit compondo os valores de entrada da tarefa, arquivos e linha de comando.

O diretório de trabalho do pipeline é organizado como mostrado abaixo:

```txt
work/
├── 12
│   └── 1adacb582d2198cd32db0e6f808bce
│       ├── genome.fa -> /data/../genome.fa
│       └── index
│           ├── hash.bin
│           ├── header.json
│           ├── indexing.log
│           ├── quasi_index.log
│           ├── refInfo.json
│           ├── rsd.bin
│           ├── sa.bin
│           ├── txpInfo.bin
│           └── versionInfo.json
├── 19
│   └── 663679d1d87bfeafacf30c1deaf81b
│       ├── ggal_gut
│       │   ├── aux_info
│       │   │   ├── ambig_info.tsv
│       │   │   ├── expected_bias.gz
│       │   │   ├── fld.gz
│       │   │   ├── meta_info.json
│       │   │   ├── observed_bias.gz
│       │   │   └── observed_bias_3p.gz
│       │   ├── cmd_info.json
│       │   ├── libParams
│       │   │   └── flenDist.txt
│       │   ├── lib_format_counts.json
│       │   ├── logs
│       │   │   └── salmon_quant.log
│       │   └── quant.sf
│       ├── ggal_gut_1.fq -> /data/../ggal_gut_1.fq
│       ├── ggal_gut_2.fq -> /data/../ggal_gut_2.fq
│       └── index -> /data/../asciidocs/day2/work/12/1adacb582d2198cd32db0e6f808bce/index
```

!!! info

    Você pode criar esse plot usando a função `tree` se você tiver instalado. No unix, simplesmente use `sudo apt install -y tree` ou com Homebrew: `brew install tree`

## Como funciona resume

A opção de linha de comando `-resume` permite a continuação da execução do pipeline pelo último passo que foi completado com sucesso:

```bash
nextflow run <script> -resume
```

Em termos práticos, o pipeline é executado do início. Entretanto, antes do lançamento da execução do `process`, Nextflow usa a ID única de tarefa para checar se o diretório de trabalho existe e se ele contém um estado de saída válido do comando com os esperados arquivos de saída.

Se a condição satisfeita a tarefa é ignorada e os resultado computados previamente são usados como resultados do `process`.

A primeira tarefa que a nova saída é computada invalida todas execuções posteriores no restante do DAG.

## Diretório de trabalho

Os diretórios de trabalho da tarefa são criados na pasta `work` no caminho de inicialização por padrão. Isso é suposto ser uma área de armazenamento **provisória** que pode ser limpada quando a computação termina.

!!! note

    Saídas finais do workflow são supostas a ser guardadas em uma localização diferente especificada usando uma ou mais diretivas [publishDir](https://www.nextflow.io/docs/latest/process.html#publishdir).

!!! warning

    Tenha certeza de deletar o diretório de trabalho ocasionalmente, se não sua máquina ou ambiente estará cheia de arquivos sem uso.

Uma diferente localização para o diretório de trabalho pode ser especificada usando `-w`. Por exemplo:

```bash
nextflow run <script> -w /some/scratch/dir
```

!!! warning

    Se você deletar ou mover o diretório de trabalho do pipeline, isso irá prevenir que você use o recurso resume nas execuções posteriores.

O código hash para os arquivos de entrada são computados usando:

-   O caminho completo do arquivo
-   O tamanho do arquivo
-   A última marca temporal modificada

Portanto, apenas usar **touch** em um arquivo irá invalidar a execução da tarefa relacionada.

## Como organizar experimentos _in-silico_

É uma boa prática organizar cada **experimento** em sua própria pasta. O experimento principal input parameters should be specified using a Nextflow config file. This makes it simple to track and replicate an experiment over time.

!!! note

    In the same experiment, the same pipeline can be executed multiple times, however, launching two (or more) Nextflow instances in the same directory concurrently should be avoided.

The `nextflow log` command lists the executions run in the current folder:

```console
$ nextflow log

TIMESTAMP            DURATION  RUN NAME          STATUS  REVISION ID  SESSION ID                            COMMAND
2019-05-06 12:07:32  1.2s      focused_carson    ERR     a9012339ce   7363b3f0-09ac-495b-a947-28cf430d0b85  nextflow run hello
2019-05-06 12:08:33  21.1s     mighty_boyd       OK      a9012339ce   7363b3f0-09ac-495b-a947-28cf430d0b85  nextflow run rnaseq-nf -with-docker
2019-05-06 12:31:15  1.2s      insane_celsius    ERR     b9aefc67b4   4dc656d2-c410-44c8-bc32-7dd0ea87bebf  nextflow run rnaseq-nf
2019-05-06 12:31:24  17s       stupefied_euclid  OK      b9aefc67b4   4dc656d2-c410-44c8-bc32-7dd0ea87bebf  nextflow run rnaseq-nf -resume -with-docker
```

You can use either the **session ID** or the **run name** to recover a specific execution. For example:

```bash
nextflow run rnaseq-nf -resume mighty_boyd
```

## Execution provenance

The `log` command, when provided with a **run name** or **session ID**, can return many useful bits of information about a pipeline execution that can be used to create a provenance report.

By default, it will list the work directories used to compute each task. For example:

```console
$ nextflow log tiny_fermat

/data/.../work/7b/3753ff13b1fa5348d2d9b6f512153a
/data/.../work/c1/56a36d8f498c99ac6cba31e85b3e0c
/data/.../work/f7/659c65ef60582d9713252bcfbcc310
/data/.../work/82/ba67e3175bd9e6479d4310e5a92f99
/data/.../work/e5/2816b9d4e7b402bfdd6597c2c2403d
/data/.../work/3b/3485d00b0115f89e4c202eacf82eba
```

The `-f` (fields) option can be used to specify which metadata should be printed by the `log` command. For example:

```console
$ nextflow log tiny_fermat -f 'process,exit,hash,duration'

index    0   7b/3753ff  2.0s
fastqc   0   c1/56a36d  9.3s
fastqc   0   f7/659c65  9.1s
quant    0   82/ba67e3  2.7s
quant    0   e5/2816b9  3.2s
multiqc  0   3b/3485d0  6.3s
```

The complete list of available fields can be retrieved with the command:

```bash
nextflow log -l
```

The `-F` option allows the specification of filtering criteria to print only a subset of tasks. For example:

```console
$ nextflow log tiny_fermat -F 'process =~ /fastqc/'

/data/.../work/c1/56a36d8f498c99ac6cba31e85b3e0c
/data/.../work/f7/659c65ef60582d9713252bcfbcc310
```

This can be useful to locate specific task work directories.

Finally, the `-t` option enables the creation of a basic custom provenance report, showing a template file in any format of your choice. For example:

```html
<div>
    <h2>${name}</h2>
    <div>
        Script:
        <pre>${script}</pre>
    </div>

    <ul>
        <li>Exit: ${exit}</li>
        <li>Status: ${status}</li>
        <li>Work dir: ${workdir}</li>
        <li>Container: ${container}</li>
    </ul>
</div>
```

!!! exercise

    Save the above snippet in a file named `template.html`. Then run this command (using the correct id for your run, e.g. not `tiny_fermat`):

    ```bash
    nextflow log tiny_fermat -t template.html > prov.html
    ```

    Finally, open the `prov.html` file with a browser.

## Resume troubleshooting

If your workflow execution is not resumed as expected with one or more tasks being unexpectedly re-executed each time, these may be the most likely causes:

#### Input file changed

Make sure that there’s no change in your input file(s). Don’t forget the task unique hash is computed by taking into account the complete file path, the last modified timestamp and the file size. If any of this information has changed, the workflow will be re-executed even if the input content is the same.

#### A process modifies an input

A process should never alter input files, otherwise the `resume` for future executions will be invalidated for the same reason explained in the previous point.

#### Inconsistent file attributes

Some shared file systems, such as [NFS](https://en.wikipedia.org/wiki/Network_File_System), may report an inconsistent file timestamp (i.e. a different timestamp for the same file) even if it has not been modified. To prevent this problem use the [lenient cache strategy](https://www.nextflow.io/docs/latest/process.html#cache).

#### Race condition in global variable

Nextflow is designed to simplify parallel programming without taking care about race conditions and the access to shared resources. One of the few cases in which a race condition can arise is when using a global variable with two (or more) operators.
For example:

```groovy linenums="1"
Channel
    .of(1,2,3)
    .map { it -> X=it; X+=2 }
    .view { "ch1 = $it" }

Channel
    .of(1,2,3)
    .map { it -> X=it; X*=2 }
    .view { "ch2 = $it" }
```

The problem in this snippet is that the `X` variable in the closure definition is defined in the global scope. Therefore, since operators are executed in parallel, the `X` value can be overwritten by the other `map` invocation.

The correct implementation requires the use of the `def` keyword to declare the variable **local**.

```groovy linenums="1"
Channel
    .of(1,2,3)
    .map { it -> def X=it; X+=2 }
    .println { "ch1 = $it" }

Channel
    .of(1,2,3)
    .map { it -> def X=it; X*=2 }
    .println { "ch2 = $it" }
```

#### Non-deterministic input channels

While dataflow channel ordering is guaranteed (i.e. data is read in the same order in which it’s written in the channel), a process can declare as input two or more channels each of which is the output of a **different** process, the overall input ordering is not consistent over different executions.

In practical terms, consider the following snippet:

```groovy linenums="1"
process foo {
    input:
    tuple val(pair), path(reads)

    output:
    tuple val(pair), path('*.bam'), emit: bam_ch

    script:
    """
    your_command --here
    """
}

process bar {
    input:
    tuple val(pair), path(reads)

    output:
    tuple val(pair), path('*.bai'), emit: bai_ch

    script:
    """
    other_command --here
    """
}

process gather {
    input:
    tuple val(pair), path(bam)
    tuple val(pair), path(bai)

    script:
    """
    merge_command $bam $bai
    """
}
```

The inputs declared at line 29 and 30 can be delivered in any order because the execution order of the process `foo` and `bar` are not deterministic due to their parallel execution.

Therefore the input of the third process needs to be synchronized using the [join](https://www.nextflow.io/docs/latest/operator.html#join) operator, or a similar approach. The third process should be written as:

```groovy
...

process gather {
    input:
    tuple val(pair), path(bam), path(bai)

    script:
    """
    merge_command $bam $bai
    """
}
```
