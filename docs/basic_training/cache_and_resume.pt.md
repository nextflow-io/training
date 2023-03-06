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

    Saídas finais do fluxo de trabalho são supostas a ser guardadas em uma localização diferente especificada usando uma ou mais diretivas [publishDir](https://www.nextflow.io/docs/latest/process.html#publishdir).

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

É uma boa prática organizar cada **experimento** em sua própria pasta. Os parâmetros de entrada do experimento principal devem ser especificados usando o arquivo de configuração do Nextflow. Isso deixa simples de acompanhar e replicar o experimento ao longo do tempo.

!!! note

    No mesmo experimento, o mesmo pipeline pode ser executado diversas vezes, entretanto, inciar duas (ou mais) instâncias do Nextflow no mesmo diretório atualmente deve ser evitado.

O comando `nextflow log` lista todas execuções na pasta atual:

```console
$ nextflow log

TIMESTAMP            DURATION  RUN NAME          STATUS  REVISION ID  SESSION ID                            COMMAND
2019-05-06 12:07:32  1.2s      focused_carson    ERR     a9012339ce   7363b3f0-09ac-495b-a947-28cf430d0b85  nextflow run hello
2019-05-06 12:08:33  21.1s     mighty_boyd       OK      a9012339ce   7363b3f0-09ac-495b-a947-28cf430d0b85  nextflow run rnaseq-nf -with-docker
2019-05-06 12:31:15  1.2s      insane_celsius    ERR     b9aefc67b4   4dc656d2-c410-44c8-bc32-7dd0ea87bebf  nextflow run rnaseq-nf
2019-05-06 12:31:24  17s       stupefied_euclid  OK      b9aefc67b4   4dc656d2-c410-44c8-bc32-7dd0ea87bebf  nextflow run rnaseq-nf -resume -with-docker
```

Você pode usar tanto o **ID da seção** ou o **nome da execução** para recuperar uma execução específica. Por exemplo:

```bash
nextflow run rnaseq-nf -resume mighty_boyd
```

## Proveniência da execução

O comando `log`, quando provido do **nome da execução** ou **ID da seção**, pode retornar bits de informações importantes sobre um pipeline em execução que pode ser usado para criar um reporte de proveniência.

Por padrão, irá listar todos diretórios de trabalho usados em cada tarefa. Por exemplo:

```console
$ nextflow log tiny_fermat

/data/.../work/7b/3753ff13b1fa5348d2d9b6f512153a
/data/.../work/c1/56a36d8f498c99ac6cba31e85b3e0c
/data/.../work/f7/659c65ef60582d9713252bcfbcc310
/data/.../work/82/ba67e3175bd9e6479d4310e5a92f99
/data/.../work/e5/2816b9d4e7b402bfdd6597c2c2403d
/data/.../work/3b/3485d00b0115f89e4c202eacf82eba
```

A opção `-f` (do inglês, fields) pode ser usada para especificar qual metadado deve ser impresso pelo comando `log`. Por exemplo:

```console
$ nextflow log tiny_fermat -f 'process,exit,hash,duration'

index    0   7b/3753ff  2.0s
fastqc   0   c1/56a36d  9.3s
fastqc   0   f7/659c65  9.1s
quant    0   82/ba67e3  2.7s
quant    0   e5/2816b9  3.2s
multiqc  0   3b/3485d0  6.3s
```

A lista completa dos campos que pode ser recuperada com o comando:

```bash
nextflow log -l
```

A opção `-F` permite a especificação de um critério de filtro para imprimir apenas um subconjunto de tarefas. Por exemplo:

```console
$ nextflow log tiny_fermat -F 'process =~ /fastqc/'

/data/.../work/c1/56a36d8f498c99ac6cba31e85b3e0c
/data/.../work/f7/659c65ef60582d9713252bcfbcc310
```

Isso pode ser útil para localizar um diretório de trabalho de uma específica tarefa.

Finalmente, a opção `-t` permite a criação de um reporte básico e customizável de providência, mostrando um modelo de arquivo em qualquer formato de sua escolha. Por exemplo:

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

    Salve o trecho acima em um arquivo chamado `template.html`. Então execute o comando (usando o ID correto para sua execução, ex. não `tiny_fermat`):

    ```bash
    nextflow log tiny_fermat -t template.html > prov.html
    ```

    Finalmente, abra o arquivo `prov.html` com um navegador.

## Resolução de problemas do resume

Se a execução do seu fluxo de trabalho não foi retomada como esperado com uma ou mais tarefas sendo inesperadamente re-executadas toda vez, essas são as causas mais prováveis:

#### Arquivos de entrada mudados

Tenha certeza que não há nenhuma mudança no(s) arquivo(s) de entrada. Não esqueça que cada tarefa tem seu hash único que é computado levando em conta o caminho completo do arquivo, a última marca temporal modificada e o tamanho do arquivo. Se alguma dessas informações foi alterada, o fluxo de trabalho deve ser re-executado mesmo que o conteúdo do arquivo é o mesmo.

#### Um processo modifica uma entrada

Um processo nunca deve alterar os arquivos de entrada, se não a função `resume` para execuções futuras será invalidado pela mesma razão explicada no ponto anterior.

#### Atributos de arquivos inconsistentes

Alguns sistemas de arquivos compartilhado, como [NFS](https://en.wikipedia.org/wiki/Network_File_System), devem reportar uma marca temporal de arquivo inconsistente (ex. uma diferente marca temporal para um mesmo arquivo) até mesmo quando não foi modificado. Para prevenir esse problema use a [estratégia do cache leniente](https://www.nextflow.io/docs/latest/process.html#cache).

#### Condição de corrida em uma variável global

Nextflow é desenvolvido para simplificar programação paralela sem ter que abdicar de condições de corrida e acessar recursos compartilhados. Um dos poucos casos que uma condição de corrida pode surgir é quando uma variável global com dois (ou mais) operadores. Por exemplo:

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

O problema desse trecho é que a variável `X` na definição fechada é definida no escopo global. Portanto, desde que operadores são executados em paralelo, o valor de `X` pode ser sobrescrito por outra invocação `map`.

A implementação correta requer o uso da palavra chave `def` para declarar a variável **local**.

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

#### Canais de entrada não determinísticos.

Por enquanto que o pedido dos canais dataflow são garantidos (ex. dados são lidos na mesma ordem que são escritos pelo canal), um processo pode declarar como entrada dois ou mais canais que cada um pode ter saídas de processos **diferentes**, a entrada geral não é consistente em várias execuções.

Em termos práticos, considere o trecho a seguir:

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

As entradas declaradas nas linhas 29 e 30 podem ser colocadas em qualquer ordem porque a execução do processo `foo` e `bar` não são determinísticos por causa de sua execução paralela.

Portanto a entrada do terceiro processo precisa está sincronizada usando o operador [join](https://www.nextflow.io/docs/latest/operator.html#join), ou com uma abordagem similar. O terceiro processo deve ser escrito assim:

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
