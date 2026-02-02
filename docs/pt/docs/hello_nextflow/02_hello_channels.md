# Parte 2: Hello Channels

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tradução assistida por IA - [saiba mais e sugira melhorias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<!--
<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/lJ41WMMm44M?si=xCItHLiOQWqoqBB9&amp;list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

/// caption
:fontawesome-brands-youtube:{ .youtube } Veja [a playlist completa](https://www.youtube.com/playlist?list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik) no canal do Nextflow no YouTube.

:green_book: A transcrição do vídeo está disponível [aqui](./transcripts/02_hello_channels.md).
///
-->

Na Parte 1 deste curso (Hello World), mostramos como fornecer uma entrada variável para um processo fornecendo a entrada diretamente na chamada do processo: `sayHello(params.input)`.
Essa foi uma abordagem deliberadamente simplificada.
Na prática, essa abordagem tem grandes limitações; ou seja, ela só funciona para casos muito simples onde queremos executar o processo apenas uma vez, em um único valor.
Na maioria dos casos de uso realistas de fluxos de trabalho, queremos processar múltiplos valores (dados experimentais para múltiplas amostras, por exemplo), então precisamos de uma maneira mais sofisticada de lidar com entradas.

É para isso que servem os **canais** do Nextflow.
Canais são filas projetadas para lidar com entradas eficientemente e transportá-las de uma etapa para outra em fluxos de trabalho de múltiplas etapas, ao mesmo tempo que fornecem paralelismo integrado e muitos benefícios adicionais.

Nesta parte do curso, você aprenderá como usar um canal para lidar com múltiplas entradas de uma variedade de fontes diferentes.
Você também aprenderá a usar **operadores** para transformar o conteúdo dos canais conforme necessário.

??? info "Como começar a partir desta seção"

    Esta seção do curso pressupõe que você completou a Parte 1 do curso [Hello Nextflow](./index.md), mas se você está confortável com os conceitos básicos cobertos naquela seção, pode começar a partir daqui sem fazer nada especial.

---

## 0. Aquecimento: Execute `hello-channels.nf`

Vamos usar o script de fluxo de trabalho `hello-channels.nf` como ponto de partida.
Ele é equivalente ao script produzido ao trabalhar na Parte 1 deste curso de treinamento, exceto que mudamos o destino da saída:

```groovy title="hello-channels.nf" linenums="37" hl_lines="3"
output {
    first_output {
        path 'hello_channels'
        mode 'copy'
    }
}
```

Apenas para garantir que tudo está funcionando, execute o script uma vez antes de fazer qualquer alteração:

```bash
nextflow run hello-channels.nf --input 'Hello Channels!'
```

??? success "Saída do comando"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-channels.nf` [wise_jennings] DSL2 - revision: b24f4902d6

    executor >  local (1)
    [6f/824bc1] process > sayHello [100%] 1 of 1 ✔
    ```

Como anteriormente, você encontrará o arquivo de saída chamado `output.txt` no diretório `results/hello_channels` (como especificado no bloco `output` do script de fluxo de trabalho, mostrado acima).

??? abstract "Conteúdo do diretório"

    ```console title="results/hello_channels" hl_lines="2-3"
    results
    ├── hello_channels
    │   └── output.txt
    ├── hello_world
    │   └── output.txt
    └── output.txt -> /workspaces/training/hello-nextflow/work/8c/79499c11beea6e9d43605141f2817f/output.txt
    ```

??? abstract "Conteúdo do arquivo"

    ```console title="results/hello_channels/output.txt"
    Hello Channels!
    ```

Se isso funcionou para você, está pronto para aprender sobre canais.

---

## 1. Forneça entradas variáveis através de um canal explicitamente

Vamos criar um **canal** para passar a entrada variável para o processo `sayHello()` em vez de depender do tratamento implícito, que tem certas limitações.

### 1.1. Crie um canal de entrada

Existem várias **fábricas de canais** que podemos usar para configurar um canal.
Para manter as coisas simples por enquanto, vamos usar a fábrica de canais mais básica, chamada `channel.of`, que criará um canal contendo um único valor.
Funcionalmente, isso será semelhante a como tínhamos configurado antes, mas em vez de ter o Nextflow criando um canal implicitamente, estamos fazendo isso explicitamente agora.

Esta é a linha de código que vamos usar:

```console title="Sintaxe"
greeting_ch = channel.of('Hello Channels!')
```

Isso cria um canal chamado `greeting_ch` usando a fábrica de canais `channel.of()`, que configura um canal de fila simples, e carrega a string `'Hello Channels!'` para usar como o valor da saudação.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-pipeline-channel.svg"
</figure>

!!! note

    Estamos temporariamente voltando a usar strings codificadas em vez de usar um parâmetro CLI para fins de legibilidade. Voltaremos a usar parâmetros CLI assim que tivermos coberto o que está acontecendo no nível do canal.

No bloco workflow, adicione o código da fábrica de canais:

=== "Depois"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="4 5"
    workflow {

        main:
        // cria um canal para entradas
        greeting_ch = channel.of('Hello Channels!')
        // emite uma saudação
        sayHello(params.input)

        publish:
        first_output = sayHello.out
    }
    ```

=== "Antes"

    ```groovy title="hello-channels.nf" linenums="27"
    workflow {

        main:
        // emite uma saudação
        sayHello(params.input)

        publish:
        first_output = sayHello.out
    }
    ```

Isso ainda não está funcional, pois ainda não mudamos a entrada para a chamada do processo.

### 1.2. Adicione o canal como entrada para a chamada do processo

Agora precisamos realmente conectar nosso canal recém-criado à chamada do processo `sayHello()`, substituindo o parâmetro CLI que estávamos fornecendo diretamente antes.

No bloco workflow, faça a seguinte alteração de código:

=== "Depois"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="7"
    workflow {

        main:
        // cria um canal para entradas
        greeting_ch = channel.of('Hello Channels!')
        // emite uma saudação
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

=== "Antes"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="7"
    workflow {

        main:
        // cria um canal para entradas
        greeting_ch = channel.of('Hello Channels!')
        // emite uma saudação
        sayHello(params.input)

        publish:
        first_output = sayHello.out
    }
    ```

Isso informa ao Nextflow para executar o processo `sayHello` no conteúdo do canal `greeting_ch`.

Agora nosso fluxo de trabalho está devidamente funcional; é o equivalente explícito de escrever `sayHello('Hello Channels!')`.

### 1.3. Execute o fluxo de trabalho

Vamos executá-lo!

```bash
nextflow run hello-channels.nf
```

??? success "Saída do comando"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-channels.nf` [fabulous_crick] DSL2 - revision: 23e20f76e8

    executor >  local (1)
    [c0/4f1872] process > sayHello (1) [100%] 1 of 1 ✔
    ```

Se você fez ambas as edições corretamente, deve obter uma execução bem-sucedida.
Você pode verificar o diretório de resultados para se certificar de que o resultado ainda é o mesmo de antes.

??? abstract "Conteúdo do arquivo"

    ```console title="results/hello_channels/output.txt"
    Hello Channels!
    ```

Então aumentamos a flexibilidade do nosso fluxo de trabalho ao mesmo tempo que alcançamos o mesmo resultado final.
Isso pode parecer que estamos escrevendo mais código sem benefício tangível, mas o valor ficará claro assim que começarmos a lidar com mais entradas.

Como uma prévia disso, vamos olhar mais uma coisa antes de seguir em frente: um benefício pequeno mas conveniente de usar um canal explícito para gerenciar a entrada de dados.

### 1.4. Use `view()` para inspecionar o conteúdo do canal

Os canais do Nextflow são construídos de uma maneira que nos permite operar em seu conteúdo usando operadores, que cobriremos em detalhes mais tarde neste capítulo.

Por enquanto, vamos apenas mostrar como usar um operador super simples chamado [`view()`](https://www.nextflow.io/docs/latest/reference/operator.html#view) para inspecionar o conteúdo de um canal.
Você pode pensar em `view()` como uma ferramenta de depuração, como uma instrução `print()` em Python, ou seu equivalente em outras linguagens.

Adicione esta pequena linha ao bloco workflow:

=== "Depois"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="7"
    workflow {

        main:
        // cria um canal para entradas
        greeting_ch = channel.of('Hello Channels!')
                             .view()
        // emite uma saudação
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

=== "Antes"

    ```groovy title="hello-channels.nf" linenums="27"
    workflow {

        main:
        // cria um canal para entradas
        greeting_ch = channel.of('Hello Channels!')
        // emite uma saudação
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

A quantidade exata de espaços não importa, desde que seja um múltiplo de 4; estamos apenas tentando alinhar o início da instrução `.view()` com a parte `.of()` da construção do canal.

Agora execute o fluxo de trabalho novamente:

```bash
nextflow run hello-channels.nf
```

??? success "Saída do comando"

    ```console hl_lines="7"
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-channels.nf` [scruffy_shaw] DSL2 - revision: 2ede41e14a

    executor >  local (1)
    [ef/f7e40a] sayHello (1) [100%] 1 of 1 ✔
    Hello Channels!
    ```

Como você pode ver, isso exibe o conteúdo do canal no console.
Aqui temos apenas um elemento, mas quando começarmos a carregar múltiplos valores no canal na próxima seção, você verá que isso está configurado para exibir um elemento por linha.

### Conclusão

Você sabe como usar uma fábrica de canais básica para fornecer uma entrada a um processo.

### O que vem a seguir?

Aprenda como usar canais para fazer o fluxo de trabalho iterar sobre múltiplos valores de entrada.

---

## 2. Modifique o fluxo de trabalho para executar em múltiplos valores de entrada

Fluxos de trabalho normalmente executam em lotes de entradas que devem ser processadas em massa, então queremos atualizar o fluxo de trabalho para aceitar múltiplos valores de entrada.

### 2.1. Carregue múltiplas saudações no canal de entrada

Convenientemente, a fábrica de canais `channel.of()` que estivemos usando está bastante feliz em aceitar mais de um valor, então não precisamos modificar isso de jeito nenhum.
Podemos apenas carregar múltiplos valores no canal.

Vamos fazê-los `'Hello'`, `'Bonjour'` e `'Holà'`.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-pipeline-channel-multi.svg"
</figure>

_No diagrama, o canal é representado em verde, e a ordem dos elementos é representada como bolas de gude em um tubo: o primeiro carregado está à direita, depois o segundo no meio, depois o terceiro está à esquerda._

#### 2.1.1. Adicione mais saudações

Antes do bloco workflow, faça a seguinte alteração de código:

=== "Depois"

    ```groovy title="hello-channels.nf" linenums="30" hl_lines="2"
    // cria um canal para entradas
    greeting_ch = channel.of('Hello','Bonjour','Holà')
                         .view()
    ```

=== "Antes"

    ```groovy title="hello-channels.nf" linenums="30" hl_lines="2"
    // cria um canal para entradas
    greeting_ch = channel.of('Hello Channels')
                         .view()
    ```

A documentação nos diz que isso deve funcionar. Pode ser realmente tão simples?

#### 2.1.2. Execute o comando e observe a saída do log

Vamos tentar.

```bash
nextflow run hello-channels.nf
```

??? success "Saída do comando"

    ```console hl_lines="6"
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-channels.nf` [amazing_crick] DSL2 - revision: 59a9a5888a

    executor >  local (3)
    [f4/c9962c] process > sayHello (1) [100%] 3 of 3 ✔
    Hello
    Bonjour
    Holà
    ```

Certamente parece ter executado bem.
O monitor de execução mostra que `3 of 3` chamadas foram feitas para o processo `sayHello`, e vemos as três saudações enumeradas pela instrução `view()`, uma por linha como prometido.

No entanto, ainda há apenas uma saída no diretório de resultados:

??? abstract "Conteúdo do diretório"

    ```console title="results/hello_channels" hl_lines="3"
    results
    ├── hello_channels
    │   └── output.txt
    ├── hello_world
    │   └── output.txt
    └── output.txt -> /workspaces/training/hello-nextflow/work/8c/79499c11beea6e9d43605141f2817f/output.txt
    ```

??? abstract "Conteúdo do arquivo"

    ```console title="results/hello_channels/output.txt"
    Holà
    ```

Você deve ver uma das três saudações lá, mas a que você obteve pode ser diferente do que é mostrado aqui.
Consegue pensar por que isso pode ser?

Olhando de volta para o monitor de execução, ele nos deu apenas um caminho de subdiretório (`f4/c9962c`).
Vamos dar uma olhada lá.

??? abstract "Conteúdo do diretório"

    ```console hl_lines="9"
    work/f4/c9962ce91ef87480babcb86b2b9042/
    ├── .command.begin
    ├── .command.err
    ├── .command.log
    ├── .command.out
    ├── .command.run
    ├── .command.sh
    ├── .exitcode
    └── output.txt
    ```

??? abstract "Conteúdo do arquivo"

    ```console title="work/f4/c9962ce91ef87480babcb86b2b9042/output.txt"
    Hello
    ```

Essa nem é a mesma saudação que obtivemos no diretório de resultados! O que está acontecendo?

Neste ponto, precisamos informar que, por padrão, o sistema de log ANSI escreve o log de múltiplas chamadas ao mesmo processo na mesma linha.
Então o status de todas as três chamadas ao processo sayHello() estão chegando no mesmo lugar.

Felizmente, podemos desabilitar esse comportamento para ver a lista completa de chamadas de processo.

#### 2.1.3. Execute o comando novamente com a opção `-ansi-log false`

Para expandir o log para exibir uma linha por chamada de processo, adicione `-ansi-log false` ao comando.

```bash
nextflow run hello-channels.nf -ansi-log false
```

??? success "Saída do comando"

    ```console
     N E X T F L O W  ~  version 25.10.2
    Launching `hello-channels.nf` [desperate_monod] DSL2 - revision: 59a9a5888a
    Hello
    Bonjour
    Holà
    [23/871c7e] Submitted process > sayHello (2)
    [7f/21e2c2] Submitted process > sayHello (1)
    [f4/ea10a6] Submitted process > sayHello (3)
    ```

Desta vez vemos todas as três execuções de processo e seus subdiretórios de trabalho associados listados na saída.

Isso é muito melhor, pelo menos para um fluxo de trabalho simples.
Para um fluxo de trabalho complexo, ou um grande número de entradas, ter a lista completa exibida no terminal ficaria um pouco avassalador.
É por isso que `-ansi-log false` não é o comportamento padrão.

!!! tip

    A maneira como o status é relatado é um pouco diferente entre os dois modos de log.
    No modo condensado, o Nextflow relata se as chamadas foram concluídas com sucesso ou não.
    Neste modo expandido, ele apenas relata que foram submetidas.

De qualquer forma, agora que temos os subdiretórios de cada chamada de processo, podemos procurar seus logs e saídas.

??? abstract "Conteúdo do diretório"

    ```console
    work/23/871c7ec3642a898ecd5e6090d21300/
    ├── .command.begin
    ├── .command.err
    ├── .command.log
    ├── .command.out
    ├── .command.run
    ├── .command.sh
    ├── .exitcode
    └── output.txt
    ```

    ```console
    work/7f/21e2c2f3cc8833ef3858b236e5575c/
    ├── .command.begin
    ├── .command.err
    ├── .command.log
    ├── .command.out
    ├── .command.run
    ├── .command.sh
    ├── .exitcode
    └── output.txt
    ```

    ```console
    work/f4/ea10a680d5687596d3eaa3fcf69272/
    ├── .command.begin
    ├── .command.err
    ├── .command.log
    ├── .command.out
    ├── .command.run
    ├── .command.sh
    ├── .exitcode
    └── output.txt
    ```

??? abstract "Conteúdo do arquivo"

    ```txt title="work/23/871c7ec3642a898ecd5e6090d21300/output.txt"
    Bonjour
    ```

    ```txt title="work/7f/21e2c2f3cc8833ef3858b236e5575c/output.txt"
    Hello
    ```

    ```txt title="work/f4/ea10a680d5687596d3eaa3fcf69272/output.txt"
    Holà
    ```

Isso mostra que todos os três processos foram executados com sucesso (eba).

Dito isso, ainda temos o problema de que há apenas um arquivo de saída no diretório de resultados.

Você pode se lembrar de que codificamos o nome do arquivo de saída para o processo `sayHello`, então todas as três chamadas produziram um arquivo chamado `output.txt`.

Enquanto os arquivos de saída permanecerem nos subdiretórios de trabalho, isolados dos outros processos, isso está ok.
Mas quando são publicados no mesmo diretório de resultados, aquele que foi copiado primeiro é sobrescrito pelo próximo, e assim por diante.

### 2.2. Garanta que os nomes dos arquivos de saída sejam únicos

Podemos continuar publicando todas as saídas no mesmo diretório de resultados, mas precisamos garantir que terão nomes únicos.
Especificamente, precisamos modificar o primeiro processo para gerar um nome de arquivo dinamicamente para que os nomes de arquivo finais sejam únicos.

Então como tornamos os nomes de arquivo únicos?
Uma maneira comum de fazer isso é usar alguma informação única dos metadados das entradas (recebidos do canal de entrada) como parte do nome do arquivo de saída.
Aqui, por conveniência, vamos apenas usar a própria saudação, já que é apenas uma string curta, e prefixá-la ao nome base do arquivo de saída.

#### 2.2.1. Construa um nome de arquivo de saída dinâmico

No bloco process, faça as seguintes alterações de código:

=== "Depois"

    ```groovy title="hello-channels.nf" linenums="6" hl_lines="7 11"
    process sayHello {

        input:
        val greeting

        output:
        path "${greeting}-output.txt"

        script:
        """
        echo '${greeting}' > '${greeting}-output.txt'
        """
    }
    ```

=== "Antes"

    ```groovy title="hello-channels.nf" linenums="6" hl_lines="7 11"
    process sayHello {

        input:
        val greeting

        output:
        path 'output.txt'

        script:
        """
        echo '${greeting}' > output.txt
        """
    }
    ```

Certifique-se de substituir `output.txt` tanto na definição de saída quanto no bloco de comando `script:`.

!!! tip

    Na definição de saída, você DEVE usar aspas duplas em torno da expressão do nome do arquivo (NÃO aspas simples), caso contrário falhará.

Isso deve produzir um nome de arquivo de saída único toda vez que o processo for chamado, para que possa ser distinguido das saídas de outras chamadas ao mesmo processo no diretório de saída.

#### 2.2.2. Execute o fluxo de trabalho

Vamos executá-lo. Note que estamos de volta a executar com as configurações de log ANSI padrão.

```bash
nextflow run hello-channels.nf
```

??? success "Saída do comando"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-channels.nf` [sharp_minsky] DSL2 - revision: 16a291febe

    executor >  local (3)
    [e8/33ee64] sayHello (2) [100%] 3 of 3 ✔
    Hello
    Bonjour
    Holà
    ```

Voltando à visualização de resumo, a saída é resumida em uma linha novamente.
Dê uma olhada no diretório `results` para ver se todas as saudações de saída estão lá.

??? abstract "Conteúdo do diretório"

    ```console
    results/hello_channels/
    ├── Bonjour-output.txt
    ├── Hello-output.txt
    ├── Holà-output.txt
    └── output.txt
    ```

Sim! E cada uma tem o conteúdo esperado.

??? abstract "Conteúdo do arquivo"

    ```console title="Bonjour-output.txt"
    Bonjour
    ```

    ```console title="Hello-output.txt"
    Hello
    ```

    ```console title="Holà-output.txt"
    Holà
    ```

Sucesso! Agora podemos adicionar quantas saudações quisermos sem nos preocupar com arquivos de saída sendo sobrescritos.

!!! tip

    Na prática, nomear arquivos com base nos dados de entrada em si é quase sempre impraticável.
    A melhor maneira de gerar nomes de arquivo dinâmicos é passar metadados para um processo junto com os arquivos de entrada.
    Os metadados são normalmente fornecidos através de uma 'planilha de amostras' ou equivalentes.
    Você aprenderá como fazer isso mais tarde no seu treinamento de Nextflow (veja [Missão secundária de metadados](../side_quests/metadata.md)).

### Conclusão

Você sabe como alimentar múltiplos elementos de entrada através de um canal.

### O que vem a seguir?

Aprenda a usar um operador para transformar o conteúdo de um canal.

---

## 3. Forneça múltiplas entradas através de um array

Acabamos de mostrar como lidar com múltiplos elementos de entrada que foram codificados diretamente na fábrica de canais.
E se quisermos fornecer essas múltiplas entradas de uma maneira diferente?

Por exemplo, imagine que configuramos uma variável de entrada contendo um array de elementos assim:

`greetings_array = ['Hello','Bonjour','Holà']`

Podemos carregar isso em nosso canal de saída e esperar que funcione?

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-pipeline-multi-inputs-array.svg"
</figure>

Vamos descobrir.

### 3.1. Forneça um array de valores como entrada para o canal

O bom senso sugere que deveríamos ser capazes de simplesmente passar um array de valores em vez de um único valor.
Vamos tentar; precisaremos configurar a variável de entrada e carregá-la na fábrica de canais.

#### 3.1.1. Configure a variável de entrada

Vamos pegar a variável `greetings_array` que acabamos de imaginar e torná-la realidade adicionando-a ao bloco workflow:

=== "Depois"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="4 5"
    workflow {

        main:
        // declara um array de saudações de entrada
        greetings_array = ['Hello','Bonjour','Holà']
        // cria um canal para entradas
        greeting_ch = channel.of('Hello','Bonjour','Holà')
                             .view()
        // emite uma saudação
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

=== "Antes"

    ```groovy title="hello-channels.nf" linenums="27"
    workflow {

        main:
        // cria um canal para entradas
        greeting_ch = channel.of('Hello','Bonjour','Holà')
                             .view()
        // emite uma saudação
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

Isso ainda não está funcional, apenas adicionamos uma declaração para o array.

#### 3.1.2. Defina o array de saudações como entrada para a fábrica de canais

Agora vamos substituir os valores `'Hello','Bonjour','Holà'` atualmente codificados na fábrica de canais pelo `greetings_array` que acabamos de criar.

No bloco workflow, faça a seguinte alteração:

=== "Depois"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="7"
    workflow {

        main:
        // declara um array de saudações de entrada
        greetings_array = ['Hello','Bonjour','Holà']
        // cria um canal para entradas
        greeting_ch = channel.of(greetings_array)
                             .view()
        // emite uma saudação
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

=== "Antes"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="7"
    workflow {

        main:
        // declara um array de saudações de entrada
        greetings_array = ['Hello','Bonjour','Holà']
        // cria um canal para entradas
        greeting_ch = channel.of('Hello','Bonjour','Holà')
                             .view()
        // emite uma saudação
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

Isso deve estar funcional agora.

#### 3.1.3. Execute o fluxo de trabalho

Vamos tentar executá-lo:

```bash
nextflow run hello-channels.nf
```

??? failure "Saída do comando"

    ```console hl_lines="7 11 16"
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-channels.nf` [friendly_koch] DSL2 - revision: 97256837a7

    executor >  local (1)
    [a8/1f6ead] sayHello (1) | 0 of 1
    [Hello, Bonjour, Holà]
    ERROR ~ Error executing process > 'sayHello (1)'

    Caused by:
      Missing output file(s) `[Hello, Bonjour, Holà]-output.txt` expected by process `sayHello (1)`


    Command executed:

      echo '[Hello, Bonjour, Holà]' > '[Hello, Bonjour, Holà]-output.txt'

    Command exit status:
      0

    Command output:
      (empty)

    Work dir:
      /workspaces/training/hello-nextflow/work/a8/1f6ead5f3fa30a3c508e2e7cf83ffb

    Tip: you can replicate the issue by changing to the process work dir and entering the command `bash .command.run`

    -- Check '.nextflow.log' file for details
    ```

Oh não! Há um erro!

Observe a saída de `view()` e as mensagens de erro.

Parece que o Nextflow tentou executar uma única chamada de processo, usando `[Hello, Bonjour, Holà]` como um valor de string, em vez de usar as três strings no array como valores separados.

Então é o 'empacotamento' que está causando o problema.
Como fazemos o Nextflow desempacotar o array e carregar as strings individuais no canal?

### 3.2. Use um operador para transformar o conteúdo do canal

É aqui que os **[operadores](https://www.nextflow.io/docs/latest/reference/operator.html)** entram em jogo.
Você já usou o operador `.view()`, que apenas observa o que está lá.
Agora vamos olhar para operadores que nos permitem agir sobre o conteúdo de um canal.

Se você examinar a [lista de operadores](https://www.nextflow.io/docs/latest/reference/operator.html) na documentação do Nextflow, encontrará [`flatten()`](https://www.nextflow.io/docs/latest/reference/operator.html#flatten), que faz exatamente o que precisamos: desempacotar o conteúdo de um array e emiti-los como itens individuais.

#### 3.2.1. Adicione o operador `flatten()`

Para aplicar o operador `flatten()` ao nosso canal de entrada, anexamos ele à declaração da fábrica de canais.

No bloco workflow, faça a seguinte alteração de código:

=== "Depois"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="9"
    workflow {

        main:
        // declara um array de saudações de entrada
        greetings_array = ['Hello','Bonjour','Holà']
        // cria um canal para entradas
        greeting_ch = channel.of(greetings_array)
                             .view()
                             .flatten()
        // emite uma saudação
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

=== "Antes"

    ```groovy title="hello-channels.nf" linenums="27"
    workflow {

        main:
        // declara um array de saudações de entrada
        greetings_array = ['Hello','Bonjour','Holà']
        // cria um canal para entradas
        greeting_ch = channel.of(greetings_array)
                             .view()
        // emite uma saudação
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

Aqui adicionamos o operador na próxima linha para legibilidade, mas você pode adicionar operadores na mesma linha que a fábrica de canais se preferir, assim:
`greeting_ch = channel.of(greetings_array).view().flatten()`

#### 3.2.2. Refine a(s) instrução(ões) `view()`

Poderíamos executar isso imediatamente para testar se funciona, mas enquanto estamos nisso, vamos refinar como inspecionamos o conteúdo do canal.

Queremos ser capazes de contrastar como o conteúdo se parece antes e depois que o operador `flatten()` é aplicado, então vamos adicionar um segundo, E vamos adicionar um pouco de código para obtê-los rotulados mais claramente na saída.

No bloco workflow, faça a seguinte alteração de código:

=== "Depois"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="8-10"
    workflow {

        main:
        // declara um array de saudações de entrada
        greetings_array = ['Hello','Bonjour','Holà']
        // cria um canal para entradas
        greeting_ch = channel.of(greetings_array)
                             .view { greeting -> "Before flatten: $greeting" }
                             .flatten()
                             .view { greeting -> "After flatten: $greeting" }
        // emite uma saudação
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

=== "Antes"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="8-9"
    workflow {

        main:
        // declara um array de saudações de entrada
        greetings_array = ['Hello','Bonjour','Holà']
        // cria um canal para entradas
        greeting_ch = channel.of(greetings_array)
                             .view()
                             .flatten()
        // emite uma saudação
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

Você vê que adicionamos uma segunda instrução `.view`, e para cada uma delas, substituímos os parênteses vazios (`()`) por chaves contendo algum código, como `{ greeting -> "Before flatten: $greeting" }`.

Essas são chamadas _closures_. O código que elas contêm será executado para cada item no canal.
Definimos uma variável temporária para o valor interno, aqui chamada `greeting` (mas poderia ser qualquer nome arbitrário), que é usada apenas dentro do escopo dessa closure.

Neste exemplo, `$greeting` representa cada item individual carregado no canal.
Isso resultará em uma saída de console bem rotulada.

!!! info

    Em alguns pipelines você pode ver uma variável especial chamada `$it` usada dentro de closures de operadores.
    Esta é uma variável _implícita_ que permite um acesso de forma abreviada à variável interna,
    sem precisar defini-la com um `->`.

    Preferimos ser explícitos para ajudar na clareza do código, como tal a sintaxe `$it` é desencorajada e será lentamente eliminada da linguagem Nextflow.

#### 3.2.3. Execute o fluxo de trabalho

Finalmente, você pode tentar executar o fluxo de trabalho novamente!

```bash
nextflow run hello-channels.nf
```

??? success "Saída do comando"

    ```console hl_lines="7-10"
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-channels.nf` [sleepy_gutenberg] DSL2 - revision: 1db4f760ee

    executor >  local (3)
    [b1/6a1e15] sayHello (2) [100%] 3 of 3 ✔
    Before flatten: [Hello, Bonjour, Holà]
    After flatten: Hello
    After flatten: Bonjour
    After flatten: Holà
    ```

Desta vez funciona E nos dá a percepção adicional do que o conteúdo do canal parece antes e depois de executarmos o operador `flatten()`.

- Você vê que obtemos uma única instrução `Before flatten:` porque nesse ponto o canal contém um item, o array original.
  Então obtemos três instruções `After flatten:` separadas, uma para cada saudação, que agora são itens individuais no canal.

Importante, isso significa que cada item agora pode ser processado separadamente pelo fluxo de trabalho.

!!! tip

    É tecnicamente possível alcançar os mesmos resultados usando uma fábrica de canais diferente, [`channel.fromList`](https://nextflow.io/docs/latest/reference/channel.html#fromlist), que inclui uma etapa de mapeamento implícita em sua operação.
    Aqui escolhemos não usar isso para demonstrar o uso de um operador em um caso de uso simples.

### Conclusão

Você sabe como usar um operador como `flatten()` para transformar o conteúdo de um canal, e como usar o operador `view()` para inspecionar o conteúdo do canal antes e depois de aplicar um operador.

### O que vem a seguir?

Aprenda como fazer o fluxo de trabalho receber um arquivo como sua fonte de valores de entrada.

---

## 4. Leia valores de entrada de um arquivo CSV

Realisticamente, raramente se alguma vez vamos começar de um array de valores.
Muito provavelmente, teremos um ou mais arquivos contendo os dados que precisam ser processados, em algum tipo de formato estruturado.

Preparamos um arquivo CSV chamado `greetings.csv` que contém várias saudações de entrada, imitando o tipo de dados colunar que você pode querer processar em uma análise de dados real, armazenado em `data/`.
(Os números não são significativos, eles estão lá apenas para fins ilustrativos.)

```csv title="data/greetings.csv" linenums="1"
Hello,English,123
Bonjour,French,456
Holà,Spanish,789
```

Nossa próxima tarefa é adaptar nosso fluxo de trabalho para ler os valores deste arquivo.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-pipeline-multi-inputs-csv.svg"
</figure>

Vamos ver como podemos fazer isso acontecer.

### 4.1. Modifique o script para esperar um arquivo CSV como fonte de saudações

Para começar, vamos precisar fazer duas alterações principais no script:

- Mudar o parâmetro de entrada para apontar para o arquivo CSV
- Mudar a fábrica de canais para uma projetada para lidar com um arquivo

#### 4.1.1. Mude o parâmetro de entrada para apontar para o arquivo CSV

Lembra do parâmetro `params.input` que configuramos na Parte 1?
Vamos atualizá-lo para apontar para o arquivo CSV contendo nossas saudações.

Faça a seguinte edição na declaração do parâmetro:

=== "Depois"

    ```groovy title="hello-channels.nf" linenums="20" hl_lines="5"
    /*
    * Pipeline parameters
    */
    params {
        input: Path = 'data/greetings.csv'
    }
    ```

=== "Antes"

    ```groovy title="hello-channels.nf" linenums="20" hl_lines="5"
    /*
     * Pipeline parameters
     */
    input: String = 'Holà mundo!'
    ```

Isso assume que o arquivo está localizado junto com o código do fluxo de trabalho.
Você aprenderá como lidar com outros locais de dados mais tarde em sua jornada com Nextflow.

#### 4.1.2. Mude para uma fábrica de canais projetada para lidar com um arquivo

Como agora queremos usar um arquivo em vez de strings simples como entrada, não podemos usar a fábrica de canais `channel.of()` de antes.
Precisamos mudar para usar uma nova fábrica de canais, [`channel.fromPath()`](https://www.nextflow.io/docs/latest/reference/channel.html#channel-path), que tem alguma funcionalidade integrada para lidar com caminhos de arquivo.

No bloco workflow, faça a seguinte alteração de código:

=== "Depois"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="4-8"
    workflow {

        main:
        // cria um canal para entradas de um arquivo CSV
        greeting_ch = channel.fromPath(params.input)
                             .view { greeting -> "Before flatten: $greeting" }
                             // .flatten()
                             // .view { greeting -> "After flatten: $greeting" }
        // emite uma saudação
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

=== "Antes"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="4-8"
    workflow {

        main:
        // declara um array de saudações de entrada
        greetings_array = ['Hello','Bonjour','Holà']
        // cria um canal para entradas
        greeting_ch = channel.of(greetings_array)
                             .view { greeting -> "Before flatten: $greeting" }
                             .flatten()
                             .view { greeting -> "After flatten: $greeting" }
        // emite uma saudação
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

Você notará que mudamos a entrada do canal de volta para `param.input`, e deletamos a declaração `greetings_array` já que não precisaremos mais dela.
Também comentamos o `flatten()` e a segunda instrução `view()`.

#### 4.1.3. Execute o fluxo de trabalho

Vamos tentar executar o fluxo de trabalho com a nova fábrica de canais e o arquivo de entrada.

```bash
nextflow run hello-channels.nf
```

??? failure "Saída do comando"

    ```console hl_lines="5 6 9 14"
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-channels.nf` [peaceful_poisson] DSL2 - revision: a286c08ad5

    [-        ] sayHello [  0%] 0 of 1
    Before flatten: /workspaces/training/hello-nextflow/data/greetings.csv
    ERROR ~ Error executing process > 'sayHello (1)'

    Caused by:
      File `/workspaces/training/hello-nextflow/data/greetings.csv-output.txt` is outside the scope of the process work directory: /workspaces/training/hello-nextflow/work/30/e610cb4ea5ae8693f456ac3329c92f


    Command executed:

      echo '/workspaces/training/hello-nextflow/data/greetings.csv' > '/workspaces/training/hello-nextflow/data/greetings.csv-output.txt'

    Command exit status:
      -

    Command output:
      (empty)

    Work dir:
      /workspaces/training/hello-nextflow/work/30/e610cb4ea5ae8693f456ac3329c92f

    Tip: when you have fixed the problem you can continue the execution adding the option `-resume` to the run command line

    -- Check '.nextflow.log' file for details
    ```

Oh não, não funciona. Dê uma olhada no início da saída do console e na mensagem de erro.
A parte `Command executed:` é especialmente útil aqui.

Isso pode parecer um pouco familiar.
Parece que o Nextflow tentou executar uma única chamada de processo usando o próprio caminho do arquivo como um valor de string.
Então ele resolveu o caminho do arquivo corretamente, mas não analisou realmente seu conteúdo, que é o que queríamos.

Como fazemos o Nextflow abrir o arquivo e carregar seu conteúdo no canal?

Parece que precisamos de outro [operador](https://www.nextflow.io/docs/latest/reference/operator.html)!

### 4.2. Use o operador `splitCsv()` para analisar o arquivo

Olhando através da lista de operadores novamente, encontramos [`splitCsv()`](https://www.nextflow.io/docs/latest/reference/operator.html#splitCsv), que é projetado para analisar e dividir texto formatado em CSV.

#### 4.2.1. Aplique `splitCsv()` ao canal

Para aplicar o operador, anexamos ele à linha da fábrica de canais como anteriormente.

No bloco workflow, faça a seguinte alteração de código para substituir `flatten()` por `splitcsv()` (descomentado):

=== "Depois"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="6-8"
    workflow {

        main:
        // cria um canal para entradas de um arquivo CSV
        greeting_ch = channel.fromPath(params.input)
                             .view { csv -> "Before splitCsv: $csv" }
                             .splitCsv()
                             .view { csv -> "After splitCsv: $csv" }
        // emite uma saudação
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

=== "Antes"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="6-8"
    workflow {

        main:
        // cria um canal para entradas de um arquivo CSV
        greeting_ch = channel.fromPath(params.input)
                             .view { greeting -> "Before flatten: $greeting" }
                             // .flatten()
                             // .view { greeting -> "After flatten: $greeting" }
        // emite uma saudação
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

Como você pode ver, também atualizamos as instruções `view()` de antes/depois.
Tecnicamente poderíamos ter usado o mesmo nome de variável (`greeting`) mas atualizamos para algo mais apropriado (`csv`) para tornar o código mais legível por outros.

#### 4.2.2. Execute o fluxo de trabalho novamente

Vamos tentar executar o fluxo de trabalho com a lógica de análise de CSV adicionada.

```bash
nextflow run hello-channels.nf
```

??? failure "Saída do comando"

    ```console hl_lines="7-11 14 19"
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-channels.nf` [insane_fermat] DSL2 - revision: 8e62fcbeb1

    executor >  local (3)
    [24/76da2f] sayHello (2) [  0%] 0 of 3 ✘
    Before splitCsv: /workspaces/training/hello-nextflow/data/greetings.csv
    After splitCsv: [Hello, English, 123]
    After splitCsv: [Bonjour, French, 456]
    After splitCsv: [Holà, Spanish, 789]
    ERROR ~ Error executing process > 'sayHello (2)'

    Caused by:
      Missing output file(s) `[Bonjour, French, 456]-output.txt` expected by process `sayHello (2)`


    Command executed:

      echo '[Bonjour, French, 456]' > '[Bonjour, French, 456]-output.txt'

    Command exit status:
      0

    Command output:
      (empty)

    Work dir:
      /workspaces/training/hello-nextflow/work/24/76da2fcc4876b61632749f99e26a50

    Tip: you can try to figure out what's wrong by changing to the process work dir and showing the script file named `.command.sh`

    -- Check '.nextflow.log' file for details
    ```

Interessantemente, isso também falha, mas com um erro diferente.
Desta vez o Nextflow analisou o conteúdo do arquivo (eba!) mas carregou cada linha como um array, e cada array é um elemento no canal.

Precisamos dizer a ele para pegar apenas a primeira coluna em cada linha.
Então como desempacotamos isso?

Usamos anteriormente `flatten()` para desempacotar o conteúdo de um canal, mas isso não funcionaria aqui porque flatten desempacota _tudo_ (sinta-se livre para tentar se quiser ver por si mesmo).

Em vez disso, usaremos outro operador chamado `map()` que é realmente útil e aparece muito em pipelines Nextflow.

### 4.3. Use o operador `map()` para extrair as saudações

O operador [`map()`](https://www.nextflow.io/docs/latest/reference/operator.html#map) é uma ferramenta muito útil que nos permite fazer todos os tipos de mapeamentos para o conteúdo de um canal.

Neste caso, vamos usá-lo para extrair aquele único elemento que queremos de cada linha em nosso arquivo de dados.
Esta é a aparência da sintaxe:

```groovy title="Sintaxe"
.map { row -> row[0] }
```

Isso significa 'para cada linha no canal, pegue o 0º (primeiro) item que ela contém'.

Então vamos aplicar isso à nossa análise de CSV.

#### 4.3.1. Aplique `map()` ao canal

No bloco workflow, faça a seguinte alteração de código:

=== "Depois"

    ```groovy title="hello-channels.nf" linenums="27" hl_lines="9 10"
    workflow {

        main:
        // cria um canal para entradas de um arquivo CSV
        greeting_ch = channel.fromPath(params.input)
                             .view { csv -> "Before splitCsv: $csv" }
                             .splitCsv()
                             .view { csv -> "After splitCsv: $csv" }
                             .map { item -> item[0] }
                             .view { csv -> "After map: $csv" }
        // emite uma saudação
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

=== "Antes"

    ```groovy title="hello-channels.nf" linenums="27"
    workflow {

        main:
        // cria um canal para entradas de um arquivo CSV
        greeting_ch = channel.fromPath(params.input)
                             .view { csv -> "Before splitCsv: $csv" }
                             .splitCsv()
                             .view { csv -> "After splitCsv: $csv" }
        // emite uma saudação
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }
    ```

Você vê que adicionamos outra chamada `view()` para confirmar que o operador faz o que esperamos.

#### 4.3.2. Execute o fluxo de trabalho

Vamos executar isso mais uma vez:

```bash
nextflow run hello-channels.nf
```

??? success "Saída do comando"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-channels.nf` [focused_volhard] DSL2 - revision: de435e45be

    executor >  local (3)
    [54/6eebe3] sayHello (3) [100%] 3 of 3 ✔
    Before splitCsv: /workspaces/training/hello-nextflow/data/greetings.csv
    After splitCsv: [Hello, English, 123]
    After splitCsv: [Bonjour, French, 456]
    After splitCsv: [Holà, Spanish, 789]
    After map: Hello
    After map: Bonjour
    After map: Holà
    ```

Desta vez deve executar sem erros.

Olhando para a saída das instruções `view()`, você vê o seguinte:

- Uma única instrução `Before splitCsv:`: nesse ponto o canal contém um item, o caminho do arquivo original.
- Três instruções `After splitCsv:` separadas: uma para cada saudação, mas cada uma está contida dentro de um array que corresponde àquela linha no arquivo.
- Três instruções `After map:` separadas: uma para cada saudação, que agora são elementos individuais no canal.

Note que as linhas podem aparecer em uma ordem diferente na sua saída.

Você também pode olhar os arquivos de saída para verificar que cada saudação foi corretamente extraída e processada através do fluxo de trabalho.

Alcançamos o mesmo resultado de antes, mas agora temos muito mais flexibilidade para adicionar mais elementos ao canal de saudações que queremos processar modificando um arquivo de entrada, sem modificar nenhum código.
Você aprenderá abordagens mais sofisticadas para lidar com entradas complexas em um treinamento posterior.

### Conclusão

Você sabe como usar o construtor de canal `.fromPath()` e os operadores `splitCsv()` e `map()` para ler um arquivo de valores de entrada e lidar com eles apropriadamente.

De forma mais geral, você tem uma compreensão básica de como o Nextflow usa **canais** para gerenciar entradas para processos e **operadores** para transformar seu conteúdo.

### O que vem a seguir?

Faça uma grande pausa, você trabalhou duro nesta seção!

Quando estiver pronto, prossiga para [**Parte 3: Olá Fluxo de Trabalho**](./03_hello_workflow.md) para aprender como adicionar mais etapas e conectá-las em um fluxo de trabalho adequado.

---

## Quiz

<quiz>
O que é um canal no Nextflow?
- [ ] Uma especificação de caminho de arquivo
- [ ] Uma definição de processo
- [x] Uma estrutura tipo fila para passar dados entre processos
- [ ] Uma configuração de definição

Saiba mais: [1.1. Crie um canal de entrada](#11-crie-um-canal-de-entrada)
</quiz>

<quiz>
O que este código produzirá como saída?

```groovy
channel.of('Hello', 'Bonjour', 'Hola')
    .view()
```

- [ ] `['Hello', 'Bonjour', 'Hola']` (uma única lista)
- [x] Cada elemento em uma linha separada: `Hello`, `Bonjour`, `Hola`
- [ ] Nada (canais não imprimem por padrão)
- [ ] Um erro (sintaxe inválida)

Saiba mais: [1.1. Crie um canal de entrada](#11-crie-um-canal-de-entrada)
</quiz>

<quiz>
Quando um canal contém múltiplos valores, como o Nextflow lida com a execução do processo?
- [ ] O processo executa uma vez com todos os valores
- [x] O processo executa uma vez para cada valor no canal
- [ ] O processo executa apenas com o primeiro valor
- [ ] O processo executa apenas com o último valor

Saiba mais: [2. Modifique o fluxo de trabalho para executar em múltiplos valores de entrada](#2-modifique-o-fluxo-de-trabalho-para-executar-em-multiplos-valores-de-entrada)
</quiz>

<quiz>
O que o operador `flatten()` faz?
- [ ] Combina múltiplos canais em um
- [ ] Ordena elementos do canal
- [x] Desempacota arrays em elementos individuais
- [ ] Remove elementos duplicados

Saiba mais: [3.2.1. Adicione o operador `flatten()`](#321-adicione-o-operador-flatten)
</quiz>

<quiz>
Qual é o propósito do operador `view()`?
- [ ] Para filtrar o conteúdo do canal
- [ ] Para transformar elementos do canal
- [x] Para inspecionar e depurar o conteúdo do canal
- [ ] Para salvar o conteúdo do canal em um arquivo

Saiba mais: [1.4. Use `view()` para inspecionar o conteúdo do canal](#14-use-view-para-inspecionar-o-conteudo-do-canal)
</quiz>

<quiz>
O que `splitCsv()` faz?
- [ ] Cria um arquivo CSV a partir do conteúdo do canal
- [ ] Divide uma string por vírgulas
- [x] Analisa um arquivo CSV em arrays representando cada linha
- [ ] Mescla múltiplos arquivos CSV

Saiba mais: [4.2. Use o operador `splitCsv()` para analisar o arquivo](#42-use-o-operador-splitcsv-para-analisar-o-arquivo)
</quiz>

<quiz>
Qual é o propósito do operador `map()`?
- [ ] Para filtrar elementos de um canal
- [ ] Para combinar múltiplos canais
- [x] Para transformar cada elemento em um canal
- [ ] Para contar elementos em um canal

Saiba mais: [4.3. Use o operador `map()` para extrair as saudações](#43-use-o-operador-map-para-extrair-as-saudacoes)
</quiz>

<quiz>
Por que é importante usar nomes de arquivo de saída dinâmicos ao processar múltiplas entradas?
- [ ] Para melhorar o desempenho
- [ ] Para reduzir o espaço em disco
- [x] Para evitar que arquivos de saída se sobrescrevam
- [ ] Para habilitar a funcionalidade de resume

Saiba mais: [2.2. Garanta que os nomes dos arquivos de saída sejam únicos](#22-garanta-que-os-nomes-dos-arquivos-de-saida-sejam-unicos)
</quiz>
