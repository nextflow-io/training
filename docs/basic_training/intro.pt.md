---
description: Introdução ao Nextflow
---

# Introdução

## Conceitos básicos

O Nextflow é tanto um motor de orquestração de fluxo de trabalho quanto uma linguagem de domínio específico (Domain-Specific Language - DSL) que facilita a escrita de pipelines computacionais com uso intensivo de dados.

Ele foi projetado com base na ideia de que a plataforma Linux é a _língua franca_ da ciência de dados. O Linux fornece muitas ferramentas de linha de comando que, ainda que simples, são poderosas ferramentas de script que, quando encadeadas, facilitam manipulações complexas de dados.

O Nextflow estende essa abordagem, adicionando a capacidade de definir interações complexas entre programas e um ambiente de computação paralela de alto nível, baseado no modelo de programação Dataflow. Os principais recursos do Nextflow são:

-   Portabilidade e reprodutibilidade do fluxo de trabalho
-   Escalabilidade na paralelização e na implantação
-   Integração de ferramentas já existentes, sistemas e padrões da indústria

### Processos e Canais

Na prática, um pipeline Nextflow é feito juntando diferentes processos. Cada `processo` pode ser escrito em qualquer linguagem de script que possa ser executada pela plataforma Linux (Bash, Perl, Ruby, Python, etc.).

Os processos são executados de forma independente e isolados uns dos outros, ou seja, não compartilham um estado (gravável) comum. A única maneira de eles se comunicarem é por meio de filas assíncronas, chamadas de `canais`, onde o primeiro elemento a entrar, é o primeiro a sair (FIFO - First-in-First-out).

Qualquer `processo` pode definir um ou mais `canais` como uma `entrada` e `saída`. A interação entre esses processos e, em última análise, o próprio fluxo de execução do pipeline, é definido implicitamente por essas declarações de `entrada` e `saída`.

<figure class="excalidraw">
--8<-- "docs/basic_training/img/channel-process.excalidraw.svg"
</figure>

### Abstração de execução

Enquanto um `processo` define _qual_ comando ou `script` deve ser executado, o executor determina _como_ esse `script` é executado na plataforma alvo.

Se não for especificado de outra forma, os processos são executados no computador local. O executor local é muito útil para fins de desenvolvimento e teste de pipeline, no entanto, para pipelines computacionais do mundo real, uma plataforma de computação de alto desempenho (HPC) ou nuvem geralmente é necessária.

Em outras palavras, o Nextflow fornece uma abstração entre a lógica funcional do pipeline e o sistema de execução subjacente (ou sistema de tempo de execução). Assim, é possível escrever um pipeline que seja executado perfeitamente em seu computador, em um cluster ou na nuvem, sem ser modificado. Você simplesmente define a plataforma de execução alvo no arquivo de configuração.

<figure markdown>

![Abstração de execução](img/execution_abstraction.png)

</figure>

### Linguagem de script

O Nextflow implementa uma DSL declarativa que simplifica a escrita de fluxos de trabalho complexos de análise de dados como uma extensão de uma linguagem de programação de uso geral.

Essa abordagem torna o Nextflow flexível — ele fornece os benefícios de uma DSL concisa para lidar com casos de uso recorrentes com facilidade **e** a flexibilidade e o poder de uma linguagem de programação de propósito geral para lidar com casos extremos no mesmo ambiente de computação. Isso seria difícil de implementar usando uma abordagem puramente declarativa.

Em termos práticos, a linguagem de script Nextflow é uma extensão da [linguagem de programação Groovy](https://groovy-lang.org/) a qual, por sua vez, é um superconjunto da linguagem de programação Java. Groovy pode ser pensado como "Python para Java", pois simplifica a escrita do código e é mais acessível.

## Seu primeiro script

Aqui você executará seu primeiro script Nextflow (`hello.nf`), que veremos linha por linha.

Neste exemplo ilustrativo, o script recebe no primeiro processo uma string de entrada (um parâmetro chamado `params.saudacao`) e a divide em blocos de seis caracteres. O segundo processo converte os caracteres em maiúsculas. O resultado é então finalmente exibido na tela.

### Código em Nextflow

<!-- NOTE: (Phil, Jan 2023)
We can dynamically include external files using mkdocs, as follows:

```groovy title="nf-training/hello.nf" linenums="1"
--8<-- "nf-training/hello.nf"
```

This inserts a code snippet identical to the one below, and we don't have to worry about keeping the two in sync.

HOWEVER - currently the line annotations cannot be added for external files. So for now, we still need to copy the scripts.

TODO: Maybe either:
    - Rewrite docs to not use loads of annotations
    - Wait for future versions to allow annotations with external files
-->

!!! info

    Clique no ícone :material-plus-circle: no código para ver explicações.

```groovy title="nf-training/hello.nf" linenums="1"
#!/usr/bin/env nextflow
// (1)!

params.saudacao = 'Olá mundo!' // (2)!
canal_saudacao = Channel.of(params.saudacao) // (3)!

process SEPARELETRAS { // (4)!
    input: // (5)!
    val x // (6)!

    output: // (7)!
    path 'chunk_*' // (8)!

    // (9)!
    """
    printf '$x' | split -b 6 - chunk_
    """
} // (10)!

process CONVERTAEMMAIUSCULAS { // (11)!
    input: // (12)!
    path y // (13)!

    output: // (14)!
    stdout // (15)!

    // (16)!
    """
    cat $y | tr '[a-z]' '[A-Z]'
    """
} // (17)!

workflow { // (18)!
    canal_letras = SEPARELETRAS(canal_saudacao) // (19)!
    canal_resultados = CONVERTAEMMAIUSCULAS(canal_letras.flatten()) // (20)!
    canal_resultados.view{ it } // (21)!
} // (22)!
```

1. O código começa com um shebang, que declara o Nextflow como o interpretador.
2. Declara um parâmetro `saudacao` que é inicializado com o valor 'Olá mundo!'.
3. Inicializa um `canal` chamado `canal_saudacao`, que contém o valor de `params.saudacao`. Os canais são o tipo de entrada para processos no Nextflow.
4. Inicia o primeiro bloco do processo, definido como `SEPARELETRAS`.
5. Declaração de entrada para o processo `SEPARELETRAS`. As entradas podem ser valores (`val`), arquivos ou caminhos (`path`) ou outros qualificadores ([veja aqui](https://www.nextflow.io/docs/latest/process.html#inputs)).
6. Diz ao `processo` para esperar um valor de entrada (`val`), que atribuímos à variável 'x'.
7. Declaração de saída para o processo `SEPARELETRAS`.
8. Diz ao processo para esperar um ou mais arquivo de saída (`path`), com um nome de arquivo começando com 'chunk\_\*', como saída do script. O processo envia a saída como um canal.
9. Três aspas duplas iniciam e terminam o bloco de código para executar este `processo`.
   Dentro está o código a ser executado — imprimindo o valor de `entrada` x (chamado usando o prefixo do símbolo de dólar [$]), dividindo a string em pedaços com um comprimento de 6 caracteres ("Olá mu" e "ndo!") e salvando cada um para um arquivo (chunk_aa e chunk_ab).
10. Fim do primeiro bloco de processo.
11. Inicia o segundo bloco de processo, definido como `CONVERTAEMMAIUSCULAS`.
12. Declaração de entrada para o `processo` `CONVERTAEMMAIUSCULAS`.
13. Diz ao `processo` para esperar um ou mais arquivos de `entrada` (`path`; ou seja, chunk_aa e chunk_ab), que atribuímos à variável 'y'.
14. Declaração de saída para o processo `CONVERTAEMMAIUSCULAS`.
15. Diz ao processo para esperar a saída padrão (stdout) como saída e envia essa saída como um canal.
16. Três aspas duplas iniciam e terminam o bloco de código para executar este `processo`.
    Dentro do bloco, há um script para ler arquivos (cat) usando a variável de entrada '$y' e, em seguida, um pipe (|) para a conversão em maiúsculas, imprimindo na saída padrão.
17. Fim do segundo bloco de `processo`.
18. Início do bloco de fluxo de trabalho (`workflow`) onde cada processo pode ser chamado.
19. Execute o `processo` `SEPARELETRAS` no `canal_saudacao` (também conhecido como canal de saudação) e armazene a saída no canal `canal_letras`.
20. Execute o `processo` `CONVERTAEMMAIUSCULAS` no canal de letras `canal_letras`, que é achatado usando o operador `.flatten()`. Isso transforma o canal de entrada de forma que cada item seja um elemento separado. Armazenamos a saída no canal `canal_resultados`.
21. A saída final (no canal `canal_resultados`) é impressa na tela usando o operador `view` (aplicado ao nome do canal).
22. Fim do bloco do fluxo de trabalho (`workflow`).

O uso do operador `.flatten()` aqui é para dividir os dois arquivos em dois itens separados para serem colocados no próximo processo (caso contrário, eles seriam tratados como um único elemento).

### Hora de praticar

Agora copie o exemplo acima em seu editor de texto favorito e salve-o em um arquivo chamado `hello.nf`.

!!! warning

    Para o tutorial do Gitpod, verifique se você está na pasta chamada `nf-training`

Execute o script digitando o seguinte comando em seu terminal:

```bash
nextflow run hello.nf
```

A saída será semelhante ao texto mostrado abaixo:

```linenums="1"
N E X T F L O W  ~  version 22.04.5
Launching `hello.nf` [gigantic_poitras] DSL2 - revision: 197a0e289a
executor >  local (3)
[c8/c36893] process > SEPARELETRAS (1)   [100%] 1 of 1 ✔
[1a/3c54ed] process > CONVERTAEMMAIUSCULAS (2) [100%] 2 of 2 ✔
WORLD!
HELLO
```

A saída padrão mostra (linha por linha):

1. A versão do Nextflow que foi executada.
2. Os nomes do script e da versão.
3. O executor usado (no caso acima: local).
4. O primeiro `processo` é executado uma vez. A linha começa com um valor hexadecimal exclusivo (consulte a dica abaixo) e termina com as informações de porcentagem e conclusão do trabalho.
5. O segundo processo é executado duas vezes (uma vez para chunk_aa e outra para chunk_ab).
6. A string de resultado de stdout é impressa na tela.

!!! info

    Os números hexadecimais, como `c8/c36893`, identificam a execução do processo exclusivo. Esses números também são o prefixo dos diretórios onde cada processo é executado. Você pode inspecionar os arquivos produzidos mudando para o diretório `$PWD/work` e usando esses números para encontrar o caminho de execução específico do processo.

!!! tip

    O segundo processo é executado duas vezes, em dois diretórios de trabalho diferentes para cada arquivo de entrada. A saída de log [ANSI](https://en.wikipedia.org/wiki/ANSI_escape_code) do Nextflow é atualizada dinamicamente conforme o pipeline é executado; no exemplo anterior, o diretório de trabalho `[1a/3c54ed]` é o segundo dos dois diretórios que foram processados (sobrescrevendo o log com o primeiro). Para imprimir para a tela todos os caminhos relevantes, desative a saída de log ANSI usando o sinalizador `-ansi-log` (por exemplo, `nextflow run hello.nf -ansi-log false`).

Vale ressaltar que o processo `CONVERTAEMMAIUSCULAS` é executado em paralelo, portanto não há garantia de que a instância que processa a primeira divisão (o chunk _Olá mu') será executada antes daquela que processa o segundo split (o chunk 'ndo!_).

Assim, pode ser que seu resultado final seja impresso em uma ordem diferente:

```
ndo!
Olá mu
```

## Modifique e retome

Nextflow keeps track of all the processes executed in your pipeline. If you modify some parts of your script, only the processes that are changed will be re-executed. The execution of the processes that are not changed will be skipped and the cached result will be used instead.

This allows for testing or modifying part of your pipeline without having to re-execute it from scratch.

For the sake of this tutorial, modify the `CONVERTAEMMAIUSCULAS` process in the previous example, replacing the process script with the string `rev $y`, so that the process looks like this:

```groovy
process CONVERTAEMMAIUSCULAS {
    input:
    path y

    output:
    stdout

    """
    rev $y
    """
}
```

Then save the file with the same name, and execute it by adding the `-resume` option to the command line:

```console
$ nextflow run hello.nf -resume

N E X T F L O W  ~  version 22.04.5
Launching `hello.nf` [amazing_becquerel] DSL2 - revision: 525206806b
executor >  local (2)
[c8/c36893] process > SEPARELETRAS (1)   [100%] 1 of 1, cached: 1 ✔
[77/cf83b6] process > CONVERTAEMMAIUSCULAS (1) [100%] 2 of 2 ✔
!dlrow
 olleH
```

You will see that the execution of the process `SEPARELETRAS` is skipped (the process ID is the same as in the first output) — its results are retrieved from the cache. The second process is executed as expected, printing the reversed strings.

!!! info

    The pipeline results are cached by default in the directory `$PWD/work`. Depending on your script, this folder can take up a lot of disk space. If you are sure you won’t need to resume your pipeline execution, clean this folder periodically.

## Parâmetros do pipeline

Pipeline parameters are simply declared by prepending the prefix `params` to a variable name, separated by a dot character. Their value can be specified on the command line by prefixing the parameter name with a double dash character, i.e. `--paramName`.

Now, let’s try to execute the previous example specifying a different input string parameter, as shown below:

```bash
nextflow run hello.nf --saudacao 'Bonjour le monde!'
```

The string specified on the command line will override the default value of the parameter. The output will look like this:

```
N E X T F L O W  ~  version 22.04.5
Launching `hello.nf` [fervent_galileo] DSL2 - revision: 525206806b
executor >  local (4)
[e9/139d7d] process > SEPARELETRAS (1)   [100%] 1 of 1 ✔
[bb/fc8548] process > CONVERTAEMMAIUSCULAS (1) [100%] 3 of 3 ✔
m el r
!edno
uojnoB
```

### Em formato de DAG

Para entender melhor como o Nextflow está lidando com os dados neste pipeline, abaixo está uma figura tipo DAG para visualizar todas as `entradas`, `saídas`, `canais` e `processos`:

<figure markdown>

![Hello world diagram](img/helloworlddiagram.png)

</figure>
