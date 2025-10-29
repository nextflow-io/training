# Part 1: Hello World

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/8X2hHI-9vms?si=F0t9LFYLjAWoyRXj&amp;list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

/// caption
:fontawesome-brands-youtube:{ .youtube } Veja a [playlist completa](https://www.youtube.com/playlist?list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik) no canal de YouTube do Nextflow.

:green_book: A transcrição desse vídeo está disponível [aqui](./transcripts/01_hello_world.md).
///

Um exemplo "Hello World!" é um modelo minimalista projetado para demonstrar a sintaxe básica e a estrutura de uma linguagem de programação ou framework de software. Esse exemplo normalmente consiste em exibir a frase "Hello, World!" no dispositivo de saída, como o console ou terminal, ou gravá-la em um arquivo.

Nesta primeira parte do curso de treinamento Hello Nextflow, introduzimos o tópico com um exemplo Hello World muito simples e independente de domínio, que será gradualmente expandido para demonstrar o uso dos componentes e da lógica fundamental do Nextflow.

---

## 0. Aquecimento: Execute Hello World diretamente

Vamos demonstrar isso com um comando simples executado diretamente no terminal, para entender seu funcionamento antes de encapsulá-lo no Nextflow.

### 0.1. Faça o terminal dizer hello

```bash
echo 'Hello World!'
```

### 0.2. Agora faça-o escrever a saída em um arquivo

```bash
echo 'Hello World!' > output.txt
```

### 0.3. Verifique se o arquivo de saída está presente usando o comando `ls`

```bash
ls
```

### 0.4. Exiba o conteúdo do arquivo

```bash
cat output.txt
```

!!! dica

    No ambiente Gitpod, você também pode encontrar o arquivo de saída no explorador de arquivos e visualizar seu conteúdo clicando nele. Alternativamente, pode usar o comando `code` para abrir o arquivo para visualização.


    ```bash
    code output.txt
    ```

### Conclusão

Agora você sabe como executar um comando simples no terminal que exibe um texto e, opcionalmente, como gravar essa saída em um arquivo.

### O que vem a seguir?

Descubra como isso seria estruturado em um fluxo de trabalho do Nextflow.

---

## 1. Experimente o script inicial do workflow Hello World

Conforme mencionado na introdução, fornecemos um script de workflow funcional, embora minimalista, chamado hello-world.nf. Ele realiza a mesma tarefa de antes (escrever "Hello World!"), mas utilizando o Nextflow.

Para começar, primeiro abriremos o script do workflow para entender sua estrutura, e só depois o executaremos (antes de tentar qualquer modificação) para verificar se ele se comporta conforme esperado.

### 1.1. Decifrando a estrutura do código

Let's open the `hello-world.nf` script in the editor pane.

!!! nota

    O arquivo está no diretório hello-nextflow, que deve ser seu diretório de trabalho atual. Você pode abrir o arquivo clicando no mesmo no explorador de arquivo ou digitar `ls` no terminal e Cmd+Click (MacOS) ou Ctrl+Click (PC).

```groovy title="hello-world.nf" linenums="1"
#!/usr/bin/env nextflow

/*
 * Use echo para exibir 'Hello World!'
 */
process sayHello {

    output:
        stdout

    script:
    """
    echo 'Hello World!'
    """
}

workflow {

    // emita uma saudação
    sayHello()
}
```

Como você pode ver, um script Nextflow envolve dois tipos de componentes principais: um ou mais **processos** e o **fluxo de trabalho** propriamente dito.
Cada **processo** descreve a(s) operação(ões) que a etapa correspondente no pipeline deve realizar, enquanto o **fluxo de trabalho** descreve a lógica do fluxo de dados que conecta as várias etapas.
Vamos dar uma olhada mais de perto no bloco **process** primeiro e, depois, no bloco **workflow**.

#### 1.1.1. A definição `process`

O primeiro bloco de código descreve um **processo**. A definição do processo começa com a palavra-chave `process`, seguida do nome do processo e, finalmente, do corpo do processo delimitado por chaves.
O corpo do processo deve conter um bloco de script que especifica o comando a ser executado, que pode ser qualquer coisa que você executaria em um terminal de linha de comando.

Aqui temos um **processo** chamado `sayHello` que grava sua **saída** em `stdout`.

```groovy title="hello-world.nf" linenums="3"
/*
 * Use echo para imprimir "Hello World!" na saída padrão
 */
processo sayHello {
    output:
        stdout
    script:
    """
    echo 'Hello World!'
    """
}
```

Essa é uma definição mínima de processo que contém apenas uma definição de saída e o próprio script.
Em um pipeline do mundo real, um processo geralmente contém blocos adicionais, como diretivas, entradas e cláusulas condicionais, que apresentaremos mais adiante neste curso de treinamento.

!!! nota
A definição de saída não _determina_ qual saída será criada.
Ela simplesmente _declara_ qual é a saída esperada, de modo que o Nextflow possa procurá-la quando a execução estiver concluída.
Isso é necessário para verificar se o comando foi executado com êxito e para passar a saída para processos posteriores, se necessário.

#### 1.1.2. A definição `workflow`

O segundo bloco de código descreve o **fluxo de trabalho** em si.
A definição do fluxo de trabalho começa com a palavra-chave `workflow`, seguida de um nome opcional e, em seguida, o corpo do fluxo de trabalho delimitado por chaves.
Aqui temos um **fluxo de trabalho** que consiste em uma chamada para o processo `sayHello`.

```groovy title="hello-world.nf" linenums="16"
workflow {
    // emite uma saudação
    sayHello()
}
```

Essa é uma definição mínima de **fluxo de trabalho**.
Em um pipeline do mundo real, o fluxo de trabalho normalmente contém várias chamadas para **processos** conectados por **canais**.
Você aprenderá como adicionar mais processos e conectá-los por canais daqui a pouco.

### 1.2. Executando o fluxo de trabalho

Olhar para o código não é tão divertido quanto executá-lo, portanto, vamos testar isso na prática.

```bash
nextflow run hello-world.nf
```

A saída do console deve ser semelhante a esta:

```console title="Output"
 N E X T F L O W   ~  version 24.02.0-edge

 ┃ Launching `hello-world.nf` [reverent_carson] DSL2 - revision: 463b611a35

executor >  local (1)
[1c/7d08e6] sayHello [100%] 1 of 1 ✔
```

Parabéns, você acabou de executar seu primeiro fluxo de trabalho Nextflow!

O resultado mais importante aqui é a última linha (linha 6), que informa que o processo `sayHello` foi executado com sucesso uma vez.

Ok, isso é ótimo, mas onde podemos encontrar o resultado?
A definição do processo `sayHello` dizia que a saída seria enviada para a saída padrão, mas nada foi impresso no console, não é mesmo?

### 1.3. Localizando a saída e os logs no diretório `work`

Quando você executa o Nextflow pela primeira vez em um determinado diretório, ele cria um diretório chamado `work` no qual gravará todos os arquivos (e links simbólicos) gerados durante a execução.

Dê uma olhada lá dentro; você encontrará um subdiretório nomeado com um hash (para torná-lo exclusivo; discutiremos o motivo daqui a pouco), aninhado em dois níveis de profundidade e contendo alguns arquivos de log.

!!! dica
Se você procurar o conteúdo do subdiretório de tarefas no explorador de arquivos VSCode do Gitpod, verá todos esses arquivos imediatamente.
No entanto, esses arquivos estão configurados para serem ocultados no terminal. Portanto, se quiser usar o `ls` ou o `tree` para visualizá-los, será necessário definir a opção relevante para exibir arquivos ocultados.

    ```bash
    tree -a work
    ```

    Você deverá ver algo parecido com isto, embora os nomes exatos dos subdiretórios sejam diferentes em seu sistema.

    ```console title="Directory contents"
    work
    └── 1c
        └── 7d08e685a7aa7060b9c21667924824
            ├── .command.begin
            ├── .command.err
            ├── .command.log
            ├── .command.out
            ├── .command.run
            ├── .command.sh
            └── .exitcode
    ```

Você deve ter notado que os nomes dos subdiretórios apareceram (de forma truncada) na saída da execução do fluxo de trabalho, na linha que diz:

```console title="Output"
[1c/7d08e6] sayHello [100%] 1 of 1 ✔
```

Isso informa qual é o caminho do subdiretório para essa chamada de processo específica (às vezes chamada de tarefa).

!!! nota
O Nextflow cria um subdiretório exclusivo e separado para cada chamada de processo.
Ele prepara os arquivos de entrada relevantes, o script e outros arquivos auxiliares lá, e grava todos os arquivos de saída e registros lá também.

Se olharmos dentro do subdiretório, encontraremos os seguintes arquivos de log:

- **`.command.begin`**: Metadados relacionados ao início da execução da tarefa do processo.
- **`.command.err`**: Mensagens de erro (stderr) emitidas pela tarefa de processo.
- **`.command.log`**: Saída de registro completa emitida pela tarefa de processo.
- **`.command.out`**: Saída regular (stdout) pela tarefa de processo.
- **`.command.sh`**: O comando que foi executado pela chamada da tarefa de processo
- **`.exitcode`**: O código de saída resultante do comando.

Nesse caso, você pode procurar sua saída no arquivo `.command.out`, pois é nele que a saída stdout é capturada.
Se você abri-lo, encontrará a saudação `Hello World!`, que era o resultado esperado de nosso fluxo de trabalho minimalista.

Também vale a pena dar uma olhada no arquivo `.command.sh`, que informa qual comando o Nextflow realmente executou. Nesse caso, é muito simples, mas mais adiante no curso você verá comandos que envolvem alguma interpolação de variáveis. Quando estiver lidando com isso, você precisará ser capaz de verificar exatamente o que foi executado, especialmente ao solucionar um problema.

### Conclusão

Você agora sabe como decifrar um script simples do Nextflow, executá-lo e encontrar a saída e os logs no diretório de trabalho.

### O que vem a seguir?

Saiba como fazer com que o script produza um arquivo nomeado.

---

## 2. Enviando a saída para um arquivo

Em vez de imprimir "Hello World!" na saída padrão, seria melhor salvar essa saída em um arquivo específico, exatamente como fizemos ao executar no terminal anteriormente.
É assim que a maioria das ferramentas que você executará como parte dos pipelines do mundo real normalmente se comporta; veremos exemplos disso mais tarde.
Para obter esse resultado, o script e os blocos de definição de saída precisam ser atualizados.

### 2.1. Alterando o comando `process` para gerar um arquivo nomeado

Essa é a mesma alteração que fizemos quando executamos o comando diretamente no terminal anteriormente.

_Antes:_

```groovy title="hello-world.nf" linenums="11"
'''
echo 'Hello World!'
'''
```

_Depois:_

```groovy title="hello-world.nf" linenums="11"
'''
echo 'Hello World!' > output.txt
'''
```

### 2.2. Alterando a declaração de saída no processo `sayHello`

Precisamos informar ao Nextflow que agora ele deve procurar um arquivo específico a ser produzido pela execução do processo.

_Antes:_

```groovy title="hello-world.nf" linenums="8"
output:
    stdout
```

_Depois:_

```groovy title="hello-world.nf" linenums="8"
output:
    path 'output.txt'
```

!!! nota
As entradas e saídas nos blocos de processo normalmente exigem um qualificador e um nome de variável:

    ```
    <qualificador de entrada/saída> <nome da entrada/saída>
    ```

    O qualificador define o tipo de dados a serem recebidos.
    Essas informações são usadas pelo Nextflow para aplicar as regras semânticas associadas a cada qualificador e tratá-las adequadamente.

    Os qualificadores comuns incluem `val` e `path`.

    No exemplo acima, `stdout` é uma exceção, pois não está associado a um nome.

### 2.3. Executando o fluxo de trabalho novamente

```bash
nextflow run hello-world.nf
```

A saída do log deve ser muito semelhante à primeira vez que você executou o fluxo de trabalho:

```console title="Output"
 N E X T F L O W   ~  version 24.02.0-edge

 ┃ Launching `hello-world.nf` [cranky_sinoussi] DSL2 - revision: 30b437bb96

executor >  local (1)
[7a/6bd54c] sayHello [100%] 1 of 1 ✔
```

Como fez anteriormente, localize o diretório `work` no explorador de arquivos.
Lá, localize o arquivo de saída `output.txt`, clique nele para abri-lo e verifique se ele contém a saudação conforme esperado.

!!! aviso
Este exemplo é frágil porque codificamos o nome do arquivo de saída em dois lugares diferentes (nos blocos `script` e `output`).
Se alterarmos um, mas não o outro, o script será interrompido.
Mais tarde, você aprenderá a usar variáveis para evitar esse problema.

### 2.4. Adicionando uma diretiva `publishDir` ao processo

Você deve ter notado que a saída está enterrada em um diretório de trabalho com várias camadas de profundidade.
O Nextflow está no controle desse diretório e não devemos interagir com ele.
Para tornar o arquivo de saída mais acessível, podemos utilizar a diretiva `publishDir`.

Ao especificar essa diretiva, estamos dizendo ao Nextflow para copiar automaticamente o arquivo de saída para um diretório de saída designado.
Isso nos permite deixar o diretório de trabalho sozinho e, ao mesmo tempo, ter acesso fácil ao arquivo de saída desejado.

_Antes:_

```groovy title="hello-world.nf" linenums="6"
process sayHello {

    output:
        path 'output.txt'
```

_Depois:_

```groovy title="hello-world.nf" linenums="6"
process sayHello {

    publishDir 'results', mode: 'copy'

    output:
        path 'output.txt'
```

!!! nota
Há uma nova opção de sintaxe que possibilita declarar e publicar saídas no nível do fluxo de trabalho, documentada [aqui] (https://www.nextflow.io/docs/latest/workflow.html#publishing-outputs), o que torna redundante o uso do `publishDir` no nível do processo quando o pipeline estiver totalmente operacional.
No entanto, o `publishDir` ainda é muito útil durante o desenvolvimento do pipeline; é por isso que o incluímos nesta série de treinamento.
Isso também garantirá que você possa ler e entender o grande número de pipelines que já foram escritos com o `publishDir`.
Você aprenderá a usar a sintaxe de saídas em nível de fluxo de trabalho mais adiante nesta série de treinamento.

### 2.5. Executando o fluxo de trabalho novamente

```bash
nextflow run hello-world.nf
```

A saída do log deve começar a parecer muito familiar:

```console title="Output"
 N E X T F L O W   ~  version 24.02.0-edge

 ┃ Launching `hello-world.nf` [mighty_lovelace] DSL2 - revision: 6654bc1327

executor >  local (1)
[10/15498d] sayHello [100%] 1 of 1 ✔
```

Desta vez, o Nextflow terá criado um novo diretório chamado `results/`.
Nesse diretório, está o nosso arquivo `output.txt`.
Se você verificar o conteúdo, ele deverá corresponder à saída em nosso diretório work/task.
É assim que movemos os arquivos de resultados para fora dos diretórios de trabalho.

### Conclusão

Você agora já sabe como enviar resultados para um arquivo nomeado específico e usar a diretiva `publishDir` para mover arquivos para fora do diretório de trabalho do Nextflow.

### O que vem a seguir?

Saiba como fazer com que o Nextflow retome a execução de um pipeline usando resultados em cache de uma execução anterior para ignorar quaisquer etapas que já tenham sido concluídas com êxito.

---

## 3. Usando o recurso `resume` do Nextflow

O Nextflow tem uma opção chamada `resume`, que permite que você execute novamente um pipeline que já tenha sido iniciado anteriormente.
Quando iniciado com `-resume`, qualquer processo que já tenha sido executado exatamente com o mesmo código, configurações e entradas será ignorado.
Usar esse modo significa que o Nextflow só executará processos novos, que tenham sido modificados ou que estejam recebendo novas configurações ou entradas.

Há duas vantagens principais em fazer isso:

- Se você estiver no meio do desenvolvimento do pipeline, poderá iterar mais rapidamente, pois só precisará executar efetivamente o(s) processo(s) em que estiver trabalhando ativamente para testar as alterações.

- Se estiver executando um pipeline em produção e algo der errado, em muitos casos, você poderá corrigir o problema e reiniciar o pipeline, e ele voltará a ser executado a partir do ponto de falha, o que pode economizar muito tempo e computação.

### 3.1. Executando o fluxo de trabalho novamente com `-resume`

```bash
nextflow run hello-world.nf -resume
```

A saída do console deve ser semelhante.

```console title="Output"
 N E X T F L O W   ~  version 24.02.0-edge

 ┃ Launching `hello-world.nf` [thirsty_gautier] DSL2 - revision: 6654bc1327

[10/15498d] sayHello [100%] 1 of 1, cached: 1 ✔
```

Observe o bit adicional `cached:` na linha de status do processo, o que significa que o Nextflow reconheceu que já havia feito esse trabalho e simplesmente reutilizou o resultado da última execução.

!!! nota
Quando você executa novamente um pipeline com `resume`, o Nextflow não sobrescreve nenhum arquivo gravado em um diretório publishDir por qualquer chamada de processo que tenha sido executada com sucesso anteriormente.

### Conclusão

Agora você sabe como reiniciar um pipeline sem repetir etapas que já foram executadas de maneira idêntica.

### O que vem a seguir?

Saiba como adicionar entradas variáveis.

---

## 4. Adicionando entradas variáveis usando um canal

Até agora, estamos emitindo uma saudação codificada no comando `process`.
Agora, vamos adicionar alguma flexibilidade usando uma variável de entrada, para que possamos alterar facilmente a saudação.

Para isso, precisamos fazer uma série de alterações inter-relacionadas:

1. Informar o processo sobre as entradas de variáveis esperadas usando o bloco `input:`
2. Editar o processo para usar a entrada
3. Criar um **canal** para passar a entrada para o processo (falaremos mais sobre isso em um minuto)
4. Adicionar o canal como entrada à chamada do processo

### 4.1. Adicionando uma definição de `input` ao bloco do processo

Primeiro, precisamos adaptar a definição do processo para aceitar uma entrada.

_Antes:_

```groovy title="hello-world.nf" linenums="6"
process sayHello {

    publishDir 'results', mode: 'copy'

    output:
        path "output.txt"
```

_Depois:_

```groovy title="hello-world.nf" linenums="6"
process sayHello {

    publishDir 'results', mode: 'copy'

    input:
        val greeting

    output:
        path "output.txt"
```

### 4.2. Editando o comando do processo para usar a variável de entrada

Agora, trocamos o valor original codificado para a variável de entrada.

_Antes:_

```groovy title="hello-world.nf" linenums="16"
"""
echo 'Hello World!' > output.txt
"""
```

_Depois:_

```groovy title="hello-world.nf" linenums="16"
"""
echo '$greeting' > output.txt
"""
```

### 4.3. Criando um canal de entrada

Agora que o nosso processo espera uma entrada, precisamos configurar essa entrada no corpo do fluxo de trabalho.

É aqui que entram os canais: o Nextflow usa canais para alimentar as entradas dos processos e transportar dados entre os processos que estão conectados entre si.
Há várias maneiras de fazer isso, mas, por enquanto, usaremos apenas o canal mais simples possível, contendo um único valor.

Vamos criar o canal usando a fábrica `channel.of()`, que configura um canal de valor simples, e fornecer a ele uma string codificada para ser usada como saudação, declarando `greeting_ch = channel.of('Hello world!')`.

_Antes:_

```groovy title="hello-world.nf" linenums="21"
workflow {

    // emite uma saudação
    sayHello()
}
```

_Depois:_

```groovy title="hello-world.nf" linenums="21"
workflow {

    // cria um canal para os inputs
    greeting_ch = channel.of('Hello world!')

    // emite uma saudação
    sayHello()
}
```

### 4.4. Adicionando o canal como entrada à chamada de processo

Agora precisamos realmente conectar nosso canal recém-criado à chamada de processo `sayHello()`.

_Antes:_

```groovy title="hello-world.nf" linenums="26"
// emit a greeting
sayHello()
```

_Depois:_

```groovy title="hello-world.nf" linenums="26"
// emit a greeting
sayHello(greeting_ch)
```

### 4.5. Executando o comando do fluxo de trabalho novamente

Vamos executá-lo!

```bash
nextflow run hello-world.nf
```

Se você fez todas as quatro edições corretamente, deverá obter outra execução bem-sucedida:

```console title="Output"
 N E X T F L O W   ~  version 24.02.0-edge

 ┃ Launching `hello-world.nf` [prickly_avogadro] DSL2 - revision: b58b6ab94b

executor >  local (1)
[1f/50efd5] sayHello (1) [100%] 1 of 1 ✔
```

Sinta-se à vontade para verificar o diretório de resultados para se certificar de que o resultado ainda é o mesmo de antes; até agora, estamos apenas ajustando progressivamente o encanamento interno para aumentar a flexibilidade do nosso fluxo de trabalho e, ao mesmo tempo, obter o mesmo resultado final.

### Conclusão

Você agora sabe como usar um canal simples para fornecer uma entrada a um processo.

### O que vem a seguir?

Saiba como passar entradas da linha de comando.

---

## 5. Usando parâmetros da CLI para entradas

Queremos poder especificar a entrada da linha de comando, pois essa é a parte que quase sempre será diferente nas execuções subsequentes do fluxo de trabalho.

Boas notícias: O Nextflow tem um sistema de parâmetros de fluxo de trabalho integrado chamado `params`, que facilita a declaração e o uso de parâmetros da CLI.

### 5.1. Editando a declaração do canal de entrada para usar um parâmetro

Aqui, trocamos a string codificada por `params.greeting` na linha de criação do canal.

_Antes:_

```groovy title="hello-world.nf" linenums="23"
// criar um canal para os inputs
greeting_ch = channel.of('Hello world!')
```

_Depois:_

```groovy title="hello-world.nf" linenums="23"
// criar um canal para os inputs
greeting_ch = channel.of(params.greeting)
```

Isso cria automaticamente um parâmetro chamado `greeting` que você pode usar para fornecer um valor na linha de comando.

### 5.2. Executar o fluxo de trabalho novamente com o parâmetro `--greeting`

Para fornecer um valor para esse parâmetro, basta adicionar `--greeting <value>` à sua linha de comando.

```bash
nextflow run hello-world.nf --greeting 'Bonjour le monde!'
```

A execução desse comando já deve lhe parecer extremamente familiar.

```console title="Output"
 N E X T F L O W   ~  version 24.02.0-edge

 ┃ Launching `hello-world.nf` [cheesy_engelbart] DSL2 - revision: b58b6ab94b

executor >  local (1)
[1c/9b6dc9] sayHello (1) [100%] 1 of 1 ✔
```

Não se esqueça de abrir o arquivo de saída para verificar se agora você tem a nova versão da saudação. Pronto!

!!! dica
É útil distinguir os parâmetros em nível de Nextflow dos parâmetros em nível de pipeline.
Para parâmetros que se aplicam a um pipeline, usamos um hífen duplo (`--`), enquanto usamos um único hífen (`-`) para parâmetros que modificam uma configuração específica do Nextflow, por exemplo, o recurso `-resume` que usamos anteriormente.

### 5.3. Definindo um valor padrão para um parâmetro de linha de comando

Em muitos casos, faz sentido fornecer um valor padrão para um determinado parâmetro para que você não precise especificá-lo em cada execução.
Vamos inicializar o parâmetro `greeting` com um valor padrão, adicionando a declaração do parâmetro na parte superior do script (com um bloco de comentários como bônus gratuito).

```groovy title="hello-world.nf" linenums="3"
/*
 * Parâmetros do Pipeline
 */
params.greeting = "Olá mundo!"
```

### 5.4. Executando o fluxo de trabalho novamente sem especificar o parâmetro

Agora que você tem um valor padrão definido, pode executar o fluxo de trabalho novamente sem precisar especificar um valor na linha de comando.

```bash
nextflow run hello-world.nf
```

A saída deve ser a mesma.

```console title="Output"
 N E X T F L O W   ~  version 24.02.0-edge

 ┃ Launching `hello-world.nf` [wise_waddington] DSL2 - revision: 988fc779cf

executor >  local (1)
[c0/8b8332] sayHello (1) [100%] 1 of 1 ✔
```

Verifique a saída no diretório de resultados e... Tcharam! Funciona! O Nextflow usou o valor padrão para nomear a saída. Mas espere aí, o que acontece agora se fornecermos o parâmetro na linha de comando?

### 5.5. Executando o fluxo de trabalho novamente com o parâmetro `--greeting` na linha de comando usando uma saudação diferente

```bash
nextflow run hello-world.nf --greeting 'Hola Mundo!'
```

O Nextflow não está reclamando, o que é um bom sinal:

```console title="Output"
 N E X T F L O W   ~  version 24.02.0-edge

 ┃ Launching `hello-world.nf` [prickly_miescher] DSL2 - revision: 988fc779cf

executor >  local (1)
[56/f88a56] sayHello (1) [100%] 1 of 1 ✔
```

Verifique o diretório de resultados e veja o conteúdo de `output.txt`. Tcharam novamente!
O valor do parâmetro que passamos na linha de comando substituiu o valor que demos à variável no script. De fato, os parâmetros podem ser definidos de várias maneiras diferentes; se o mesmo parâmetro for definido em vários locais, seu valor será determinado com base na ordem de precedência descrita [aqui](https://www.nextflow.io/docs/latest/config.html).

!!! dica
Você pode colocar a declaração do parâmetro dentro do bloco de fluxo de trabalho, se preferir. Seja qual for a sua escolha, tente agrupar coisas semelhantes no mesmo lugar para não acabar com declarações espalhadas por toda parte.

### Conclusão

Até agora você já sabe como configurar uma variável de entrada para um processo e fornecer um valor na linha de comando.

### O que vem a seguir?

Saiba como adicionar um segundo processo e encadeá-los.

---

## 6. Adicionando uma segunda etapa ao fluxo de trabalho

A maioria dos fluxos de trabalho do mundo real envolve mais de uma etapa. Aqui apresentamos um segundo processo que converte o texto em maiúsculas (all-caps), usando o clássico UNIX one-liner:

```bash
tr '[a-z]' '[A-Z]'
```

Primeiro, vamos executar o comando sozinho no terminal para verificar se ele funciona conforme o esperado, sem que nenhum código de fluxo de trabalho atrapalhe a clareza, assim como fizemos no início com `echo 'Hello World'`. Em seguida, escreveremos um processo que faz a mesma coisa e, finalmente, conectaremos os dois processos para que a saída do primeiro sirva de entrada para o segundo.

### 6.1. Executando o comando no terminal por si só

```bash
echo 'Hello World' | tr '[a-z]' '[A-Z]'
```

A saída é simplesmente a versão em maiúsculas da string de texto:

```console title="Output"
HELLO WORLD
```

!!! nota
Esse é um one-liner de substituição de texto muito ingênuo que não leva em conta as letras acentuadas, portanto, por exemplo, 'olá' se tornará 'OLà'. Isso é esperado.

### 6.2. Fazendo com que o comando receba um arquivo como entrada e grave a saída em um arquivo

Como anteriormente, queremos enviar os resultados para um arquivo dedicado, que nomeamos prefixando o nome do arquivo original com `UPPER-`.

```bash
cat output.txt | tr '[a-z]' '[A-Z]' > UPPER-output.txt
```

Agora, a saída do `HELLO WORLD` está no novo arquivo de saída, `UPPER-output.txt`.

### 6.3. Envolvendo o comando em uma nova definição de processo do Nextflow

Podemos modelar nosso novo processo com base no primeiro, já que queremos usar todos os mesmos componentes.

```groovy title="hello-world.nf" linenums="26"
/*
 * Use um utilitário de substituição de texto para converter a saudação em maiúsculas
 */
process convertToUpper {

    publishDir 'results', mode: 'copy'

    input:
        path input_file

    output:
        path "UPPER-${input_file}"

    script:
    """
    cat '$input_file' | tr '[a-z]' '[A-Z]' > UPPER-${input_file}
    """
}
```

Como um pequeno bônus, aqui compomos o segundo nome de arquivo de saída com base no primeiro.

!!! dica
É muito importante lembrar: é necessário usar aspas duplas ao redor da expressão do nome do arquivo de saída (NÃO aspas simples) ou haverá falha.

### 6.4. Adicionando uma chamada ao novo processo no corpo do fluxo de trabalho

Não se esqueça de que precisamos dizer ao Nextflow para realmente chamar o processo que acabamos de criar! Para fazer isso, nós o adicionamos ao corpo do fluxo de trabalho.

```groovy title="hello-world.nf" linenums="44"
workflow {
    // criar um canal para entradas
    greeting_ch = channel.of(params.greeting)
    // emite uma saudação
    sayHello(greeting_ch)
    // converter a saudação em maiúsculas
    convertToUpper()
}
```

Parece bom! Mas ainda precisamos conectar a chamada do processo `convertToUpper` para ser executada na saída do `sayHello`.

### 6.5. Passando a saída do primeiro processo para o segundo processo

A saída do processo `sayHello` é automaticamente empacotada como um canal chamado `sayHello.out`, portanto, tudo o que precisamos fazer é passá-la como entrada para o processo `convertToUpper`.

```groovy title="hello-world.nf" linenums="52"
// converter a saudação em letras maiúsculas
convertToUpper(sayHello.out)
```

Para um caso simples como esse, isso é tudo o que precisamos fazer para conectar dois processos!

### 6.6. Executar o mesmo comando de fluxo de trabalho anterior

Vamos nos certificar de que isso funcione:

```bash
nextflow run hello-world.nf --greeting 'Hello World!'
```

Que emocionante! Agora há uma linha extra na saída de registro, que corresponde ao novo processo que acabamos de adicionar:

```console title="Output"
 N E X T F L O W   ~  version 24.02.0-edge

 ┃ Launching `hello-world.nf` [magical_brenner] DSL2 - revision: 0e18f34798

executor >  local (2)
[57/3836c0] sayHello (1)       [100%] 1 of 1 ✔
[ee/bb3cc8] convertToUpper (1) [100%] 1 of 1 ✔
```

Você notará que, desta vez, o fluxo de trabalho produziu dois novos subdiretórios de trabalho, um por chamada de processo.
Verifique o diretório de trabalho da chamada para o segundo processo, onde você deve encontrar dois arquivos de saída diferentes listados. Se olhar com atenção, você notará que um deles (a saída do primeiro processo) tem um pequeno ícone de seta à direita, o que significa que é um link simbólico.
Ele aponta para o local onde esse arquivo está no diretório de trabalho do primeiro processo.
Por padrão, o Nextflow usa links simbólicos para encenar arquivos de entrada sempre que possível, para evitar fazer cópias duplicadas.

!!! nota
Tudo o que fizemos foi conectar a saída do `sayHello` à entrada do `convertToUpper` e os dois processos puderam ser executados em série.
O Nextflow fez o trabalho pesado de manipular os arquivos de entrada e saída e passá-los entre os dois comandos para nós.
Esse é o poder dos canais no Nextflow, fazendo o trabalho pesado de conectar as etapas do nosso pipeline.
Além disso, o Nextflow determinará automaticamente qual chamada precisa ser executada primeiro com base em como elas estão conectadas, de modo que a ordem em que estão escritas no corpo do fluxo de trabalho não importa.
No entanto, recomendamos que você seja gentil com seus colaboradores e com seu futuro eu, e tente escrevê-las em uma ordem lógica!

### Conclusão

Você já sabe como adicionar uma segunda etapa que usa a saída da primeira etapa como entrada.

### O que vem a seguir?

Saiba como fazer com que o fluxo de trabalho seja executado em um lote de valores de entrada.

---

## 7. Modificando o fluxo de trabalho para ser executado em um lote de valores de entrada

Os fluxos de trabalho normalmente são executados em lotes de entradas que devem ser processados em massa, portanto, queremos atualizar o fluxo de trabalho para aceitar vários valores de entrada.

Convenientemente, a fábrica `channel.of()` que estamos usando aceita de bom grado mais de um valor, portanto, não precisamos modificá-la; basta carregar mais valores no canal.

### 7.1. Carregando várias saudações no canal de entrada

Para manter as coisas simples, voltamos a codificar as saudações na fábrica do canal em vez de usar um parâmetro para a entrada, mas melhoraremos isso em breve.
_Antes:_

```groovy title="hello-world.nf" linenums="46"
// criar um canal para entradas
greeting_ch = channel.of(params.greeting)
```

_Depois:_

```groovy title="hello-world.nf" linenums="46"
// criar um canal para entradas
greeting_ch = channel.of('Hello','Bonjour','Holà')
```

A documentação nos diz que isso deve funcionar. Será que realmente pode ser tão simples?

### 7.2. Executando o comando e vendo a saída do registro

Vamos tentar.

```bash
nextflow run hello-world.nf
```

Bem, ele parece estar funcionando bem.

```console title="Output"
 N E X T F L O W   ~  version 24.02.0-edge

 ┃ Launching `hello-world.nf` [lonely_pare] DSL2 - revision: b9f1d96905

executor >  local (6)
[3d/1fe62c] sayHello (2)       [100%] 3 of 3 ✔
[86/695813] convertToUpper (3) [100%] 3 of 3 ✔
```

No entanto... Isso parece indicar que foram feitas “3 de 3” chamadas para cada processo, o que é animador, mas isso nos dá apenas um caminho de subdiretório para cada um. O que está acontecendo?

Por padrão, o sistema de registro ANSI grava o registro de várias chamadas para o mesmo processo na mesma linha. Felizmente, podemos desativar esse comportamento.

### 7.3. Execute o comando novamente com a opção `-ansi-log false`

Para expandir o registro em log para exibir uma linha por chamada de processo, basta adicionar `-ansi-log false` ao comando.

```bash
nextflow run hello-world.nf -ansi-log false
```

Desta vez, vemos todos os seis subdiretórios de trabalho listados na saída:

```console title="Output"
N E X T F L O W  ~  version 24.02.0-edge
Launching `hello-world.nf` [big_woese] DSL2 - revision: 53f20aeb70
[62/d81e63] Submitted process > sayHello (1)
[19/507af3] Submitted process > sayHello (2)
[8a/3126e6] Submitted process > sayHello (3)
[12/48a5c6] Submitted process > convertToUpper (1)
[73/e6e746] Submitted process > convertToUpper (2)
[c5/4fedda] Submitted process > convertToUpper (3)
```

Isso é muito melhor, pelo menos para esse número de processos.
No caso de um fluxo de trabalho complexo ou de um grande número de entradas, ter a lista completa de saída para o terminal pode ser um pouco cansativo.

Dito isso, temos outro problema. Se você olhar o diretório `results`, há apenas dois arquivos: `output.txt` e `UPPER-output.txt`!

```console title="Conteúdo do diretório”
resultados
├── output.txt
└── UPPER-output.txt
```

O que está acontecendo com isso? Não deveríamos estar esperando dois arquivos por saudação de entrada, ou seja, seis arquivos no total?
Você deve se lembrar que codificamos o nome do arquivo de saída para o primeiro processo.
Isso foi bom desde que houvesse apenas uma única chamada feita por processo, mas quando começamos a processar vários valores de entrada e a publicar as saídas no mesmo diretório de resultados, isso se torna um problema.
Para um determinado processo, cada chamada produz uma saída com o mesmo nome de arquivo, portanto, o Nextflow simplesmente substitui o arquivo de saída anterior sempre que um novo é produzido.

### 7.4. Certifique-se de que os nomes dos arquivos de saída sejam exclusivos

Como publicaremos todos os resultados no mesmo diretório de resultados, precisamos garantir que eles tenham nomes exclusivos.
Especificamente, precisamos modificar o primeiro processo para gerar um nome de arquivo dinamicamente, de modo que os nomes dos arquivos finais sejam exclusivos.
Então, como tornar os nomes dos arquivos exclusivos? Uma maneira comum de fazer isso é usar alguma parte exclusiva de metadados como parte do nome do arquivo.
Aqui, por conveniência, usaremos apenas a saudação em si.

_Antes:_

```groovy title="hello-world.nf" linenums="11"
process sayHello {

    publishDir 'results', mode: 'copy'

    input:
        val greeting

    output:
        path "output.txt"

    script:
    """
    echo '$greeting' > "output.txt"
    """
}
```

_Depois:_

```groovy title="hello-world.nf" linenums="11"
process sayHello {

    publishDir 'results', mode: 'copy'

    input:
        val greeting

    output:
        path "${greeting}-output.txt"

    script:
    """
    echo '$greeting' > '$greeting-output.txt'
    """
}
```

Isso deve produzir um nome de arquivo de saída exclusivo para cada chamada de cada processo.

### 7.5. Execute o fluxo de trabalho e veja o diretório de resultados

Vamos executá-lo e verificar se ele funciona.

```bash
nextflow run hello-world.nf
```

Ao voltar para a visualização de resumo, a saída tem a seguinte aparência novamente:

```console title="Output"
 N E X T F L O W   ~  version 24.02.0-edge

 ┃ Launching `hello-world.nf` [jovial_mccarthy] DSL2 - revision: 53f20aeb70

executor >  local (6)
[03/f007f2] sayHello (1)       [100%] 3 of 3 ✔
[e5/dd2890] convertToUpper (3) [100%] 3 of 3 ✔
```

Mas o mais importante é que agora temos seis novos arquivos, além dos dois que já tínhamos no diretório `results`:

```console title="Directory contents"
results
├── Bonjour-output.txt
├── Hello-output.txt
├── Holà-output.txt
├── output.txt
├── UPPER-Bonjour-output.txt
├── UPPER-Hello-output.txt
├── UPPER-Holà-output.txt
└── UPPER-output.txt
```

Sucesso! Agora podemos adicionar quantas saudações quisermos sem nos preocuparmos com o fato de os arquivos de saída serem sobrescritos.

!!! nota
Na prática, nomear arquivos com base nos próprios dados de entrada é quase sempre impraticável. A melhor maneira de gerar nomes de arquivos dinâmicos é usar uma planilha de amostras que contenha metadados relevantes (como IDs de amostra exclusivos) e criar uma estrutura de dados chamada “mapa”, que passamos para os processos e da qual podemos obter um identificador apropriado para gerar os nomes de arquivos.
Mostraremos a você como fazer isso mais adiante neste curso de treinamento.

### Conclusão

Você já sabe como alimentar um lote de vários elementos de entrada por meio de um canal.

### O que vem a seguir?

Saiba como fazer com que o fluxo de trabalho use um arquivo como sua fonte de valores de entrada.

---

## 8. Modificando o fluxo de trabalho para usar um arquivo como fonte de valores de entrada

Muitas vezes, quando queremos executar um lote de vários elementos de entrada, os valores de entrada podem estar contidos em um arquivo.

Como exemplo, fornecemos a você um arquivo CSV chamado `greetings.csv` no diretório `data/`, contendo várias saudações separadas por vírgulas.

```csv title="greetings.csv"
Hello,Bonjour,Holà
```

Portanto, só precisamos modificar nosso fluxo de trabalho para ler os valores de um arquivo como esse.

### 8.1. Configurando um parâmetro da CLI com um valor padrão apontando para um arquivo de entrada

Primeiro, vamos usar o argumento `params` para configurar um novo parâmetro chamado `input_file`, substituindo o parâmetro `greeting`, agora inútil, com um valor padrão que aponta para o arquivo `greetings.csv`.

_Antes:_

```groovy title="hello-world.nf" linenums="6"
/*
 * Parâmetros do pipeline
 */
params.greeting = "Bonjour le monde!"
```

_Depois:_

```groovy title="hello-world.nf" linenums="6"
/*
 * Parâmetros do pipeline
 */
params.input_file = "data/greetings.csv"
```

### 8.2. Atualizando a declaração do canal para lidar com o arquivo de entrada

Neste ponto, apresentamos uma nova fábrica de canais, `channel.fromPath()`, que tem algumas funcionalidades integradas para lidar com caminhos de arquivos.

Vamos usá-la em vez da fábrica `channel.of()` que usamos anteriormente; a sintaxe básica é a seguinte:

```groovy title="channel construction syntax"
channel.fromPath(input_file)
```

Agora, vamos implantar um novo conceito, um “operador” para transformar esse arquivo CSV em conteúdo de canal. Você aprenderá mais sobre operadores mais tarde, mas, por enquanto, basta entendê-los como formas de transformar canais de várias maneiras.

Como nosso objetivo é ler o conteúdo de um arquivo `.csv`, adicionaremos o operador `.splitCsv()` para fazer com que o Nextflow analise o conteúdo do arquivo adequadamente, bem como o operador `.flatten()` para transformar o elemento de matriz produzido por `.splitCsv()` em um canal de elementos individuais.
Portanto, a instrução de construção do canal passa a ser:

```groovy title="channel construction syntax"
channel.fromPath(input_file)
       .splitCsv()
       .flatten()
```

E aqui está ele no contexto do corpo do fluxo de trabalho:

_Antes:_

```groovy title="hello-world.nf" linenums="46"
// criar um canal para entradas
greeting_ch = channel.of('Hello','Bonjour','Holà')
```

_Depois:_

```groovy title="hello-world.nf" linenums="46"
// criar um canal para entradas de um arquivo CSV
greeting_ch = channel.fromPath(params.input_file)
                     .splitCsv()
                     .flatten()
```

### 8.3. Executando o fluxo de trabalho (uma última vez!)

```bash
nextflow run hello-world.nf
```

Mais uma vez, vemos que cada processo é executado três vezes:

```console title="Output"
 N E X T F L O W   ~  version 24.02.0-edge

 ┃ Launching `hello-world.nf` [angry_spence] DSL2 - revision: d171cc0193

executor >  local (6)
[0e/ceb175] sayHello (2)       [100%] 3 of 3 ✔
[01/046714] convertToUpper (3) [100%] 3 of 3 ✔
```

Observando os resultados, vemos que cada saudação foi extraída e processada corretamente pelo fluxo de trabalho. Obtivemos o mesmo resultado da etapa anterior, mas agora temos muito mais flexibilidade para adicionar mais elementos ao canal de saudações que queremos processar.

!!! dica
Durante o desenvolvimento do pipeline, você pode inspecionar o conteúdo de qualquer canal adicionando o operador `.view()` ao nome do canal.
Por exemplo, se você adicionar `greeting_ch.view()` em qualquer lugar do corpo do fluxo de trabalho, ao executar o script, o Nextflow imprimirá o conteúdo do canal na saída padrão.
Você também pode usar isso para inspecionar o efeito dos operadores.

    Por exemplo, a saída de `channel.fromPath(params.input_file).splitCsv().view()` terá a seguinte aparência:

    ```console title="Output"
    [Hello, Bonjour, Holà]
    ```

    Enquanto a saída de `channel.fromPath(params.input_file).splitCsv().flatten().view()` terá a seguinte aparência:

    ```console title="Output”
    Olá
    Bom dia
    Olá
    ```

### Conclusão

Você agora sabe como fornecer os valores de entrada para o fluxo de trabalho por meio de um arquivo.
De modo mais geral, você aprendeu a usar os componentes essenciais do Nextflow e tem uma compreensão básica da lógica de como criar um fluxo de trabalho e gerenciar entradas e saídas.

### O que vem a seguir?

Comemore seu sucesso e faça uma pausa!
Não se preocupe se os tipos de canais e operadores parecerem muito difíceis de lidar na primeira vez que você os encontrar.
Você terá mais oportunidades de praticar o uso desses componentes em várias configurações ao longo deste curso de treinamento.
Quando estiver pronto, vá para a Parte 2 para aprender sobre outro conceito importante: o provisionamento do software necessário para cada processo.
