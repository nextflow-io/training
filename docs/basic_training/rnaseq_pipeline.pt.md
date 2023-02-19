---
description: Material de treinamento básico do Nextflow
---

# Pipeline simples de RNA-Seq

Para demonstrar um cenário biomédico da vida real, nós iremos implementar uma prova de conceito de pipeline RNA-Seq que:

1. Cria arquivo de índice de transcriptoma
2. Realiza controles de qualidade
3. Realiza quantificação
4. Cria um relatório MultiQC

Isso será feito usando uma série de sete scripts, cada qual se baseia no anterior, para criar um fluxo de trabalho completo. Você poderá encontrá-los no diretório do tutorial (`script1.nf` - `script7.nf`). 

## Defina os parâmetros do pipeline

Parâmetros são entradas e opções que podem ser modificados quando um pipeline é executado.

O script `script1.nf` define os parâmetros de entrada do pipeline.

```groovy
params.reads = "$projectDir/data/ggal/gut_{1,2}.fq"
params.transcriptome_file = "$projectDir/data/ggal/transcriptome.fa"
params.multiqc = "$projectDir/multiqc"

println "reads: $params.reads"
```

Execute-o usando o comando a seguir:

```bash
nextflow run script1.nf
```

Tente especificar um parâmetro de entrada diferente no seu comando de execução, por exemplo:

```bash
nextflow run script1.nf --reads '/workspace/gitpod/nf-training/data/ggal/lung_{1,2}.fq'
```

### :material-progress-question: Exercícios

!!! exercise

    Modifique o `script1.nf` ao adicionar um quarto parâmetro chamado `outdir` e defina-o como um caminho padrão que será usado como o diretório de saída do pipeline.

    ??? result

        ```groovy
        params.reads = "$projectDir/data/ggal/gut_{1,2}.fq"
        params.transcriptome_file = "$projectDir/data/ggal/transcriptome.fa"
        params.multiqc = "$projectDir/multiqc"
        params.outdir = "results"
        ```

!!! exercise

    Modifique o `script1.nf` para imprimir todos os parâmetros do pipeline usando um único comando`log.info` como uma declaração de [string multilinha](https://www.nextflow.io/docs/latest/script.html#multi-line-strings) (multiline string).

    !!! tip ""

        :material-lightbulb: Veja um exemplo [aqui](https://github.com/nextflow-io/rnaseq-nf/blob/3b5b49f/main.nf#L41-L48).


    ??? result

        Adicione o código abaixo para seu arquivo de script:

        ```groovy
        log.info """\
                    R N A S E Q - N F   P I P E L I N E
                    ===================================
                    transcriptome: ${params.transcriptome_file}
                    reads        : ${params.reads}
                    outdir       : ${params.outdir}
                    """
                    .stripIndent()
        ```

### :material-check-all: Resumo

Nesta etapa você aprendeu:

1. Como definir parâmetros em seu script de pipeline
2. Como atribuir parâmetros usando a linha de comando
3. O uso de `$var` e `${var}` como espaço reservado para variáveis 
4. Como usar strings multilinhas
5. Como usar `log.info` para imprimir informações e salvá-las no arquivo de execução de log

## Criar um arquivo de índice de transcriptoma

Nextflow permite a execução de qualquer comando ou script usando uma definição de `processo`.

Um `processo` é definido ao fornecer três principais declarações: as [`entradas`](https://www.nextflow.io/docs/latest/process.html#inputs), [`saídas`](https://www.nextflow.io/docs/latest/process.html#outputs) e comandos de [`script`](https://www.nextflow.io/docs/latest/process.html#script) do processo.

Para adicionar uma etapa de processamento de índice do transcriptoma `INDEX`, tente adicionar os blocos de código a seguir no seu `script1.nf`. Como alternativa, esses blocos de código já foram adicionados ao `script2.nf`.

```groovy
/*
 * define the `index` process that creates a binary index
 * given the transcriptome file
 */
process INDEX {
  input:
  path transcriptome

  output:
  path 'salmon_index'

  script:
  """
  salmon index --threads $task.cpus -t $transcriptome -i salmon_index
  """
}
```

Além disso, adicione um escopo de fluxo de trabalho contendo uma definição de canal de entrada e o processo de índice:

```groovy
workflow {
    index_ch = INDEX(params.transcriptome_file)
}
```

Aqui, o parâmetro `params.transcriptome_file` é usado como entrada para o processo `INDEX`. O processo `INDEX` (usando a ferramenta `salmon`) cria `salmon_index`, um arquivo índice de transcriptoma que é passado como saída ao canal `index_ch`.

!!! info

    A declaração de `entrada` define a variável de caminho `transcriptoma` que é usada no `script` como uma referência (usando o símbolo de cifrão) na linha de comando Salmon.

!!! warning

    Os requisitos de recursos, como CPUs e limites de memória, podem mudar com diferentes execuções e plataformas de fluxo de trabalho. Nextflow pode usar `$task.cpus` como uma variável para o número de CPUs.. Veja a [documentação de diretivas de processo](https://www.nextflow.io/docs/latest/process.html#directives) para mais detalhes.

Execute-o usando o comando:

```bash
nextflow run script2.nf
```

A execução irá falhar porque `salmon` não está instalado em seu ambiente.

Adicione a opção de linha de comando `-with-docker` para iniciar a execução através do container Docker, como mostrado abaixo:

```bash
nextflow run script2.nf -with-docker
```

Dessa vez a execução vai funcionar porque usa o container Docker `nextflow/rnaseq-nf` que é definido no arquivo `nextflow.config` do seu diretório atual. Se você está executando esse script localmente, você precisará baixar o docker em seu computador, fazer log in e ativar o docker, e permitir que o script baixe o container contendo os scripts de execução. Você pode aprender mais sobre o docker [aqui](https://www.nextflow.io/docs/latest/docker.html).

Para evitar adicionar `-with-docker` cada vez que você executar o script, adicione a linha a seguir ao arquivo `nextflow.config`:

```groovy
docker.enabled = true
```

### :material-progress-question: Exercícios

!!! exercise

    Ative a execução do Docker por padrão adicionando a configuração acima no arquivo `nextflow.config`.

!!! exercise

    Imprima a saída do canal `index_ch` usando o operador [view](https://www.nextflow.io/docs/latest/operator.html#view).

    ??? result

        Adicione o código a seguir ao final do bloco de fluxo de trabalho em seu arquivo de script

        ```groovy
        index_ch.view()
        ```

!!! exercise

    Se você tiver mais CPUs disponíveis, tente alterar seu script para solicitar mais recursos para este processo. Por exemplo, consulte os [documentos de diretiva](https://www.nextflow.io/docs/latest/process.html#cpus). `$task.cpus` já está especificado no script, portanto definir o número de CPUs como uma diretiva informará ao Nextflow para executar este trabalho.

    ??? result

        Adicione `cpus 2` no top do processo de índice:

        ```groovy
        process INDEX {
            cpus 2
            input:
            ...
        ```

        Em seguida verifique se funcionou observando o script executado no diretório de trabalho. Procure pelo hexadecimal (por exemplo, `work/7f/f285b80022d9f61e82cd7f90436aa4/`), depois `cat` o arquivo `.command.sh`.

!!! exercise "Bonus Exercise"

    Use o comando `tree work` para observar como Nextflow organiza o diretório de trabalho do processo. Verifique [aqui](https://www.tecmint.com/linux-tree-command-examples/) se você precisa baixar `tree`.

    ??? result

        Deve ser algo assim:

        ```
        work
        ├── 17
        │   └── 263d3517b457de4525513ae5e34ea8
        │       ├── index
        │       │   ├── complete_ref_lens.bin
        │       │   ├── ctable.bin
        │       │   ├── ctg_offsets.bin
        │       │   ├── duplicate_clusters.tsv
        │       │   ├── eqtable.bin
        │       │   ├── info.json
        │       │   ├── mphf.bin
        │       │   ├── pos.bin
        │       │   ├── pre_indexing.log
        │       │   ├── rank.bin
        │       │   ├── refAccumLengths.bin
        │       │   ├── ref_indexing.log
        │       │   ├── reflengths.bin
        │       │   ├── refseq.bin
        │       │   ├── seq.bin
        │       │   └── versionInfo.json
        │       └── transcriptome.fa -> /workspace/Gitpod_test/data/ggal/transcriptome.fa
        ├── 7f
        ```

### :material-check-all: Resumo

Nesta etapa você aprendeu:

1. Como definir um processo executando um comando personalizado
2. Como as entradas do processo são declaradas
3. Como as saídas do processo são declaradas
4. Como imprimir o conteúdo de um canal
5. Como acessar o número de CPUs disponíveis

## Colete arquivos lidos por pares

Essa etapa mostra como combinar arquivos **lidos** em pares, para que eles possam ser pareados pelo **Salmon**.

Edite o script `script3.nf` adicionando a seguinte declaração como a última linha do arquivo:

```groovy
read_pairs_ch.view()
```

Salve-o e execute-o com o comando a seguir:

```bash
nextflow run script3.nf
```

Isso irá imprimir algo semelhante a isso:

```bash
[gut, [/.../data/ggal/gut_1.fq, /.../data/ggal/gut_2.fq]]
```

O exemplo acima mostra como o canal `read_pairs_ch` emite tuplas compostas de dois elementos, onde o primeiro é o prefixo do par de leitura e o segundo é uma lista que representa os arquivos reais.

Tente novamente especificando arquivos de leitura diferentes usando o padrão glob:

```bash
nextflow run script3.nf --reads 'data/ggal/*_{1,2}.fq'
```

!!! warning

    Caminhos de arquivo que incluem um ou mais caracteres especiais, como `*`, `?` etc., DEVEM ser colocados entre aspas simples para evitar que o Bash expanda o glob.

### :material-progress-question: Exercícios

!!! exercise

    Use o operador [set](https://www.nextflow.io/docs/latest/operator.html#set) no lugar da atribuição `=` para definir o canal `read_pairs_ch`.

    ??? result

        ```groovy
        Channel
            .fromFilePairs( params.reads )
            .set { read_pairs_ch }
        ```

!!! exercise

    Use a opção `checkIfExists` para o método [fromFilePairs](https://www.nextflow.io/docs/latest/channel.html#fromfilepairs) para checar se o caminho especificado contém os pares de arquivos.

    ??? result

        ```groovy
        Channel
            .fromFilePairs( params.reads, checkIfExists: true )
            .set { read_pairs_ch }
        ```

### :material-check-all: Resumo

Nessa etapa você aprendeu:

1. Como usar `fromFilePairs` para lidar com pares de arquivos de leitura
2. Como usar a opção `checkIfExists` para checar a existência de arquivos de entrada
3. Como usar o operador `set` para definir uma uma nova variável de canal

!!! info

    A declaração de um canal pode ser feita antes do escopo do fluxo de trabalho dentro dele. Desde que a declaração esteja acima do processo que requer o canal específico.

## Realizar a quantificação da expressão

`script4.nf` adiciona um processo `QUANTIFICATION` para quantificação de expressão e chama dentro do escopo do fluxo de trabalho. A quantificação requer o arquivo de índice de transcriptoma e os arquivos fastq do par de leitura de RNA-Seq.

No escopo do fluxo de trabalho, observe como o canal `index_ch` é designado como saída do processo `INDEX`.

A seguir, note que o primeiro canal de entrada para o processo de `QUANTIFICATION` é o `index_ch` declarado previamente, que contém o caminho para `salmon_index`.

Além disso, observe que o segundo canal de entrada para o processo `QUANTIFICATION` é o `read_pair_ch` que acabamos de criar. Este sendo uma `tupla` composta de dois elementos (um valor: `sample_id` e a lista de caminhos para os arquivos de leitura fastq: `reads`) para corresponder à estrutura dos itens emitidos pela fábrica de canais `fromFilePairs`. 

Execute-o usando o comando a seguir:

```bash
nextflow run script4.nf -resume
```

Você irá ver a execução do processo `QUANTIFICATION`.

Ao usar a opção `-resume`, qualquer etapa que já foi processada é ignorada.

Tente executar o mesmo script novamente com mais arquivos de leitura, como mostrado abaixo:

```bash
nextflow run script4.nf -resume --reads 'data/ggal/*_{1,2}.fq'
```

Você irá perceber que o processo `QUANTIFICATION` é executado múltiplas vezes.

Nextflow paraleliza a execução de seu pipeline simplesmente fornecendo vários conjuntos de dados de entrada para seu script.

!!! tip

    Pode ser útil aplicar configurações opcionais a um processo específico usando [diretivas](https://www.nextflow.io/docs/latest/process.html#directives) especificando-as no corpo do processo.

### :material-progress-question: Exercícios

!!! exercise

    Adicione uma diretiva de [tag](https://www.nextflow.io/docs/latest/process.html#tag), ao processo `QUANTIFICATION` para fornecer um log de execução mais legível.

    ??? result

        Adicione o código a seguir antes da declaração de entrada:

        ```groovy
        tag "Salmon on $sample_id"
        ```

!!! exercise

    Adicione a diretiva [publishDir](https://www.nextflow.io/docs/latest/process.html#publishdir) para o processo `QUANTIFICATION` para armazenar os resultados do processo em um diretório de sua escolha.

    ??? result

        Adicione o código a seguir antes da declaração de `entrada` no processo de `QUANTIFICATION`:

        ```groovy
        publishDir params.outdir, mode:'copy'
        ```

### :material-check-all: Resumo

Nessa etapa você aprendeu:

1. Como conectar dois processos juntos usando declarações de canal
2. Como `retomar` a execução de script e pular etapas em cache
3. Como usar a diretiva `tag` para fornecer uma saída de execução mais legível
4. Como usar a diretiva `publishDir` para armanezar os resultados do processo em um caminho da sua escolha

## Controle de qualidade

A seguir, nós implementamos uma etapa de controle de qualidade `FASTQC` para seus arquivos de leitura de entrada (usando a etiqueta `fastqc`). As entradas são as mesmas que os pares de arquivos de leitura na etapa `QUANTIFICATION`.

Você pode executá-lo usando o comando a seguir:

```bash
nextflow run script5.nf -resume
```

Nextflow DSL2 sabe como dividir `reads_pair_ch` em dois canais idênticos, já que eles são requiridos duas vezes como entrada para os processos `FASTQC` e `QUANTIFICATION`.

## Relatório MultiQC

Essa etapa coleta as saídas dos processos `QUANTIFICATION` e `FASTQC` para criar um relatório final usando a ferramenta [MultiQC](http://multiqc.info/).

Execute o próximo script com o comando a seguir:

```bash
nextflow run script6.nf -resume --reads 'data/ggal/*_{1,2}.fq'
```

Isso cria um relatório final no pasta de `resultados` no diretório de `trabalho` atual.

Neste script, observe que o uso dos operadores [mix](https://www.nextflow.io/docs/latest/operator.html#mix) e [collect](https://www.nextflow.io/docs/latest/operator.html#collect) de forma encadeada para reunir as saídas dos processos `QUANTIFICATION` e `FASTQC` como uma única entrada. [Operadores](https://www.nextflow.io/docs/latest/operator.html) podem ser usados para combinar e transformar canais.

```groovy
MULTIQC(quant_ch.mix(fastqc_ch).collect())
```

Queremos que apenas uma tarefa do MultiQC seja executada para produzir um relatório. Portanto, usamos o operador de canal `mix` para combinar os dois canais, seguido pelo operador `collect` para retornar os conteúdos completos do canal como um único elemento.

### :material-check-all: Resumo

Nessa etapa você aprendeu:

1. Como coletar muitas saídas para uma única entrada com o operador `collect`
2. Como combinar com `mix` dois canais em um único canal
3. Como encadear dois ou mais operadores juntos

## Handle completion event

This step shows how to execute an action when the pipeline completes the execution.

Note that Nextflow processes define the execution of **asynchronous** tasks i.e. they are not executed one after another as if they were written in the pipeline script in a common **imperative** programming language.

The script uses the `workflow.onComplete` event handler to print a confirmation message when the script completes.

Try to run it by using the following command:

```bash
nextflow run script7.nf -resume --reads 'data/ggal/*_{1,2}.fq'
```

## Email notifications

Send a notification email when the workflow execution completes using the `-N <email address>` command-line option.

Note: this requires the configuration of a SMTP server in the nextflow config file. Below is an example `nextflow.config` file showing the settings you would have to configure:

```groovy
mail {
    from = 'info@nextflow.io'
    smtp.host = 'email-smtp.eu-west-1.amazonaws.com'
    smtp.port = 587
    smtp.user = "xxxxx"
    smtp.password = "yyyyy"
    smtp.auth = true
    smtp.starttls.enable = true
    smtp.starttls.required = true
}
```

See [mail documentation](https://www.nextflow.io/docs/latest/mail.html#mail-configuration) for details.

## Custom scripts

Real-world pipelines use a lot of custom user scripts (BASH, R, Python, etc.) Nextflow allows you to consistently use and manage these scripts. Simply put them in a directory named `bin` in the pipeline project root. They will be automatically added to the pipeline execution `PATH`.

For example, create a file named `fastqc.sh` with the following content:

```bash
#!/bin/bash
set -e
set -u

sample_id=${1}
reads=${2}

mkdir fastqc_${sample_id}_logs
fastqc -o fastqc_${sample_id}_logs -f fastq -q ${reads}
```

Save it, give execute permission, and move it into the `bin` directory as shown below:

```bash
chmod +x fastqc.sh
mkdir -p bin
mv fastqc.sh bin
```

Then, open the `script7.nf` file and replace the `FASTQC` process’ script with the following code:

```groovy
script:
"""
fastqc.sh "$sample_id" "$reads"
"""
```

Run it as before:

```bash
nextflow run script7.nf -resume --reads 'data/ggal/*_{1,2}.fq'
```

### :material-check-all: Summary

In this step you have learned:

1. How to write or use existing custom scripts in your Nextflow pipeline.
2. How to avoid the use of absolute paths by having your scripts in the `bin/` folder.

## Metrics and reports

Nextflow can produce multiple reports and charts providing several runtime metrics and execution information.

Run the [rnaseq-nf](https://github.com/nextflow-io/rnaseq-nf) pipeline previously introduced as shown below:

```bash
nextflow run rnaseq-nf -with-docker -with-report -with-trace -with-timeline -with-dag dag.png
```

The `-with-docker` option launches each task of the execution as a Docker container run command.

The `-with-report` option enables the creation of the workflow execution report. Open the file `report.html` with a browser to see the report created with the above command.

The `-with-trace` option enables the creation of a tab separated file containing runtime information for each executed task. Check the `trace.txt` for an example.

The `-with-timeline` option enables the creation of the workflow timeline report showing how processes were executed over time. This may be useful to identify the most time consuming tasks and bottlenecks. See an example at [this link](https://www.nextflow.io/docs/latest/tracing.html#timeline-report).

Finally, the `-with-dag` option enables the rendering of the workflow execution direct acyclic graph representation. Note: This feature requires the installation of [Graphviz](http://www.graphviz.org/) on your computer. See [here](https://www.nextflow.io/docs/latest/tracing.html#dag-visualisation) for further details. Then try running :

```bash
dot -Tpng dag.dot > graph.png
open graph.png
```

!!! warning

    Run time metrics may be incomplete for runs with short running tasks, as in the case of this tutorial.

!!! info

    You view the HTML files by right-clicking on the file name in the left side-bar and choosing the **Preview** menu item.

## Run a project from GitHub

Nextflow allows the execution of a pipeline project directly from a GitHub repository (or similar services, e.g., BitBucket and GitLab).

This simplifies the sharing and deployment of complex projects and tracking changes in a consistent manner.

The following GitHub repository hosts a complete version of the workflow introduced in this tutorial: <https://github.com/nextflow-io/rnaseq-nf>

You can run it by specifying the project name and launching each task of the execution as a Docker container run command:

```bash
nextflow run nextflow-io/rnaseq-nf -with-docker
```

It automatically downloads the container and stores it in the `$HOME/.nextflow` folder.

Use the command `info` to show the project information:

```bash
nextflow info nextflow-io/rnaseq-nf
```

Nextflow allows the execution of a specific revision of your project by using the `-r` command line option. For example:

```bash
nextflow run nextflow-io/rnaseq-nf -r v2.1 -with-docker
```

Revision are defined by using Git tags or branches defined in the project repository.

Tags enable precise control of the changes in your project files and dependencies over time.

## More resources

-   [Nextflow documentation](http://docs.nextflow.io) - The Nextflow docs home.
-   [Nextflow patterns](https://github.com/nextflow-io/patterns) - A collection of Nextflow implementation patterns.
-   [CalliNGS-NF](https://github.com/CRG-CNAG/CalliNGS-NF) - A Variant calling pipeline implementing GATK best practices.
-   [nf-core](http://nf-co.re/) - A community collection of production ready genomic pipelines.
