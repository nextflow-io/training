---
title: Cache e reentrância
description: Material de treinamento básico do Nextflow
---

# Execução de cache e de reentrância

O mecanismo de caching do Nextflow funciona atribuindo uma ID única para cada tarefa que é usada para criar um diretório de execução separado onde as tarefas são executadas e os resultados guardados.

A ID única de tarefa é gerada como uma hash de 128-bit compondo os valores de entrada da tarefa, arquivos e a string de comando.

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

    Você pode criar esse plot usando a função `tree` se você a tiver instalada. No unix, simplesmente use `sudo apt install -y tree` ou com Homebrew: `brew install tree`

## Como funciona a reentrância

A opção de linha de comando `-resume` permite a continuação da execução do pipeline pelo último passo que foi completado com sucesso:

```bash
nextflow run <script> -resume
```

Em termos práticos, o pipeline é executado do início. Entretanto, antes do lançamento da execução de um `process`, o Nextflow usa a ID única da tarefa para checar se o diretório de trabalho existe e se ele contém um estado de saída válido do comando com os esperados arquivos de saída.

Se esta condição for satisfeita a tarefa é ignorada e os resultados computados previamente são usados como resultados do `process`.

A primeira tarefa que tem uma nova saída computada invalida todas execuções posteriores no que resta do Grafo Acílico Direcionado (DAG, do inglês Directed Acyclic Graph).

## Diretório de trabalho

O diretório de trabalho das tarefas é criado por padrão na pasta `work` no mesmo diretório onde o pipeline foi executado. Essa localização é supostamente uma área de armazenamento **provisória** que pode ser limpada quando a execução do pipeline for finalizado.

!!! note

    As saídas finais do fluxo de trabalho geralmente são guardadas em uma localização diferente especificada usando uma ou mais diretivas [publishDir](https://www.nextflow.io/docs/latest/process.html#publishdir).

!!! warning

    Certifique-se de deletar o diretório de trabalho ocasionalmente, se não sua máquina ou ambiente estará cheia de arquivos sem uso.

Uma localização diferente para o diretório de trabalho pode ser especificada usando a opção `-w`. Por exemplo:

```bash
nextflow run <script> -w /algum/diretorio/de/scratch
```

!!! warning

    Se você deletar ou mover o diretório de trabalho do pipeline, isso irá impedir que você use o recurso de reentrância nas execuções posteriores.

O código hash para os arquivos de entrada são computados usando:

-   O caminho completo do arquivo
-   O tamanho do arquivo
-   A última marcação de tempo de modificação

Portanto, o simples uso do **touch** em um arquivo irá invalidar a execução da tarefa relacionada.

## Como organizar experimentos _in-silico_

É uma boa prática organizar cada **experimento** em sua própria pasta. Os parâmetros de entrada do experimento principal devem ser especificados usando o arquivo de configuração do Nextflow. Isso torna simples acompanhar e replicar o experimento ao longo do tempo.

!!! note

    No mesmo experimento, o mesmo pipeline pode ser executado diversas vezes, entretanto, iniciar duas (ou mais) instâncias do Nextflow no mesmo diretório ao mesmo tempo deve ser evitado.

O comando `nextflow log` lista todas as execuções na pasta atual:

```console
$ nextflow log

TIMESTAMP            DURATION  RUN NAME          STATUS  REVISION ID  SESSION ID                            COMMAND
2019-05-06 12:07:32  1.2s      focused_carson    ERR     a9012339ce   7363b3f0-09ac-495b-a947-28cf430d0b85  nextflow run hello
2019-05-06 12:08:33  21.1s     mighty_boyd       OK      a9012339ce   7363b3f0-09ac-495b-a947-28cf430d0b85  nextflow run rnaseq-nf -with-docker
2019-05-06 12:31:15  1.2s      insane_celsius    ERR     b9aefc67b4   4dc656d2-c410-44c8-bc32-7dd0ea87bebf  nextflow run rnaseq-nf
2019-05-06 12:31:24  17s       stupefied_euclid  OK      b9aefc67b4   4dc656d2-c410-44c8-bc32-7dd0ea87bebf  nextflow run rnaseq-nf -resume -with-docker
```

Você pode usar tanto o **ID da sessão** ou o **nome da execução** para recuperar uma execução específica. Por exemplo:

```bash
nextflow run rnaseq-nf -resume mighty_boyd
```

## Proveniência da execução

O comando `log`, quando provido do **nome da execução** ou **ID da sessão**, pode retornar algumas informações importantes sobre um pipeline em execução que podem ser usadas para criar um relatório de proveniência.

Por padrão, o comando irá listar todos diretórios de trabalho usados em cada tarefa. Por exemplo:

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

A lista completa dos campos disponíveis pode ser recuperada com o comando:

```bash
nextflow log -l
```

A opção `-F` permite a especificação de um critério de filtro para imprimir apenas um subconjunto de tarefas. Por exemplo:

```console
$ nextflow log tiny_fermat -F 'process =~ /fastqc/'

/data/.../work/c1/56a36d8f498c99ac6cba31e85b3e0c
/data/.../work/f7/659c65ef60582d9713252bcfbcc310
```

Isso pode ser útil para localizar um diretório de trabalho de uma tarefa específica.

Finalmente, a opção `-t` permite a criação de um relatório básico e customizável de proveniência, mostrando um modelo de arquivo em qualquer formato de sua escolha. Por exemplo:

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

## Resolução de problemas de reentrância

Se a execução do seu fluxo de trabalho não foi retomada como esperado com uma ou mais tarefas sendo inesperadamente re-executadas toda vez, essas são as causas mais prováveis:

#### Arquivos de entrada mudados

Tenha certeza que não há nenhuma mudança no(s) arquivo(s) de entrada. Não esqueça que cada tarefa tem seu hash único que é computado levando em conta o caminho completo do arquivo, a última marcação de tempo de modificação e o tamanho do arquivo. Se alguma dessas informações foi alterada, o fluxo de trabalho deve ser re-executado mesmo que o conteúdo do arquivo seja o mesmo.

#### Um processo modifica uma entrada

Um processo nunca deve alterar os arquivos de entrada, se não a função `resume` em execuções futuras será invalidada pela mesma razão explicada no ponto anterior.

#### Atributos de arquivos inconsistentes

Alguns sistemas de arquivos compartilhados, como o [NFS](https://en.wikipedia.org/wiki/Network_File_System), podem reportar uma marcação de tempo inconsistente para os arquivos (por exemplo, uma diferente marcação de tempo para um mesmo arquivo) até mesmo quando este não foi modificado. Para prevenir esse problema use a [estratégia de cache leniente](https://www.nextflow.io/docs/latest/process.html#cache).

#### Condição de corrida em uma variável global

O Nextflow é desenvolvido para simplificar programação paralela, de modo que você não precise se preocupar com condições de corrida e acesso a recursos compartilhados. Um dos poucos casos que uma condição de corrida pode surgir é quando uma variável global é utilizada com dois (ou mais) operadores. Por exemplo:

```groovy linenums="1"
Channel
    .of(1,2,3)
    .map { it -> X=it; X+=2 }
    .view { "canal1 = $it" }

Channel
    .of(1,2,3)
    .map { it -> X=it; X*=2 }
    .view { "canal2 = $it" }
```

O problema desse trecho é que a variável `X` na clausura é definida no escopo global. Portanto, como operadores são executados em paralelo, o valor de `X` pode ser sobrescrito pela outra invocação do `map`.

A implementação correta requer o uso da palavra chave `def` para declarar a variável **local**.

```groovy linenums="1"
Channel
    .of(1,2,3)
    .map { it -> def X=it; X+=2 }
    .println { "canal1 = $it" }

Channel
    .of(1,2,3)
    .map { it -> def X=it; X*=2 }
    .println { "canal2 = $it" }
```

#### Canais de entrada não determinísticos.

Embora a ordem de elementos em canais dataflow seja garantida (ou seja, os dados são lidos na mesma ordem que são escritos no canal), um processo pode declarar como entrada dois ou mais canais que são canais de saída de processos **diferentes**, fazendo com que a entrada não seja consistente entre diferentes execuções.

Em termos práticos, considere o trecho a seguir:

```groovy linenums="1"
process foo {
    input:
    tuple val(par), path(leituras)

    output:
    tuple val(par), path('*.bam'), emit: canal_bam

    script:
    """
    seu_comando --aqui
    """
}

process bar {
    input:
    tuple val(par), path(leituras)

    output:
    tuple val(par), path('*.bai'), emit: canal_bai

    script:
    """
    outro_comando --aqui
    """
}

process gather {
    input:
    tuple val(par), path(bam)
    tuple val(par), path(bai)

    script:
    """
    comando_de_mesclar $bam $bai
    """
}
```

As entradas declaradas nas linhas 29 e 30 podem ser colocadas em qualquer ordem porque a ordem da execução dos processos `foo` e `bar` não é determinística por causa de sua execução em paralelo.

Portanto a entrada do terceiro processo precisa estar sincronizada usando o operador [join](https://www.nextflow.io/docs/latest/operator.html#join), ou com uma abordagem similar. O terceiro processo deve ser escrito assim:

```groovy
...

process gather {
    input:
    tuple val(par), path(bam), path(bai)

    script:
    """
    comando_de_mesclar $bam $bai
    """
}
```
