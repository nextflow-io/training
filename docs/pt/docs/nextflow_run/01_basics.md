# Parte 1: Operações básicas de execução

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tradução assistida por IA - [saiba mais e sugira melhorias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Nesta primeira parte do curso de treinamento Nextflow Run, começamos suavemente o tópico com um exemplo muito básico Hello World independente de domínio, que usaremos para demonstrar operações essenciais e apontar os componentes de código Nextflow correspondentes.

??? info "O que é um exemplo Hello World?"

    Um "Hello World!" é um exemplo minimalista que visa demonstrar a sintaxe básica e estrutura de uma linguagem de programação ou framework de software.
    O exemplo tipicamente consiste em imprimir a frase "Hello, World!" para o dispositivo de saída, como o console ou terminal, ou escrevê-la em um arquivo.

---

## 1. Execute um Hello World diretamente

Vamos demonstrar este conceito com um comando simples que executamos diretamente no terminal, para mostrar o que ele faz antes de envolvê-lo no Nextflow.

!!! tip "Dica"

    Lembre-se de que você agora deve estar dentro do diretório `nextflow-run/` conforme descrito na página [Primeiros Passos](00_orientation.md).

### 1.1. Faça o terminal dizer olá

Execute o seguinte comando no seu terminal.

```bash
echo 'Hello World!'
```

??? success "Saída do comando"

    ```console
    Hello World!
    ```

Isso imprime o texto 'Hello World' ali mesmo no terminal.

### 1.2. Escreva a saída em um arquivo

Executar pipelines geralmente envolve ler dados de arquivos e escrever resultados em outros arquivos, então vamos modificar o comando para escrever a saída de texto em um arquivo para tornar o exemplo um pouco mais relevante.

```bash
echo 'Hello World!' > output.txt
```

??? success "Saída do comando"

    ```console

    ```

Isso não imprime nada no terminal.

### 1.3. Encontre a saída

O texto 'Hello World' agora deve estar no arquivo de saída que especificamos, chamado `output.txt`.
Você pode abri-lo no explorador de arquivos ou pela linha de comando usando o utilitário `cat`, por exemplo.

??? abstract "Conteúdo do arquivo"

    ```console title="output.txt" linenums="1"
    Hello World!
    ```

Isso é o que vamos tentar replicar com nosso primeiro fluxo de trabalho Nextflow.

### Conclusão

Agora você sabe como executar um comando simples no terminal que produz algum texto e, opcionalmente, como fazê-lo escrever a saída em um arquivo.

### O que vem a seguir?

Descubra o que é necessário para executar um fluxo de trabalho Nextflow que alcança o mesmo resultado.

---

## 2. Execute o fluxo de trabalho

Fornecemos a você um script de fluxo de trabalho chamado `1-hello.nf` que recebe uma saudação de entrada através de um argumento de linha de comando chamado `--input` e produz um arquivo de texto contendo essa saudação.

Não vamos olhar o código ainda; primeiro vamos ver como é executá-lo.

### 2.1. Inicie o fluxo de trabalho e monitore a execução

No terminal, execute o seguinte comando:

```bash
nextflow run 1-hello.nf --input 'Hello World!'
```

??? success "Saída do comando"

    ```console hl_lines="6"
    N E X T F L O W   ~  version 25.10.2

    Launching `1-hello.nf` [goofy_torvalds] DSL2 - revision: c33d41f479

    executor >  local (1)
    [a3/7be2fa] sayHello | 1 of 1 ✔
    ```

Se sua saída do console se parece com isso, então parabéns, você acabou de executar seu primeiro fluxo de trabalho Nextflow!

A saída mais importante aqui é a última linha, que está destacada na saída acima:

```console
[a3/7be2fa] sayHello | 1 of 1 ✔
```

Isso nos diz que o processo `sayHello` foi executado com sucesso uma vez (`1 of 1 ✔`).

Ótimo, mas você pode estar se perguntando: onde está a saída?

### 2.2. Encontre o arquivo de saída no diretório `results`

Este fluxo de trabalho está configurado para publicar sua saída em um diretório de resultados.
Se você olhar para o seu diretório atual, verá que quando executou o fluxo de trabalho, o Nextflow criou um novo diretório chamado `results`, bem como um subdiretório chamado `1-hello` dentro dele, contendo um arquivo chamado `output.txt`.

```console title="results/"
results
└── 1-hello
    └── output.txt
```

Abra o arquivo; o conteúdo deve corresponder à string que você especificou na linha de comando.

```console title="results/1-hello/output.txt" linenums="1"
Hello World!
```

Ótimo, nosso fluxo de trabalho fez o que deveria fazer!

### 2.3. Salve os resultados em um diretório diferente

Por padrão, o Nextflow salvará as saídas do pipeline em um diretório chamado `results` no seu caminho atual.
Para mudar onde seus arquivos são publicados, use a flag CLI `-output-dir` (ou `-o` de forma abreviada)

!!! danger "Aviso"

    Note que `--input` tem dois traços e `-output-dir` tem um!
    Isso ocorre porque `--input` é um _parâmetro_ do pipeline e `-output-dir` é uma flag CLI central do Nextflow.
    Mais sobre isso adiante.

```bash
nextflow run 1-hello.nf --input 'Hello World!' -output-dir hello_results
```

??? success "Saída do comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `1-hello.nf` [hungry_celsius] DSL2 - revision: f048d6ea78

    executor >  local (1)
    [a3/1e1535] sayHello [100%] 1 of 1 ✔
    ```

Você deve ver que suas saídas agora são publicadas em um diretório chamado `hello_results` em vez de `results`:

```console title="hello_results/"
hello_results
└── 1-hello
    └── output.txt
```

Os arquivos dentro deste diretório são os mesmos de antes, é apenas o diretório de nível superior que é diferente.
No entanto, esteja ciente em ambos os casos de que o resultado 'publicado' é uma cópia (ou em alguns casos um link simbólico) da saída real produzida pelo Nextflow quando ele executou o fluxo de trabalho.

Então agora, vamos espiar sob o capô para ver onde o Nextflow realmente executou o trabalho.

!!! Warning "Aviso"

    Nem todos os fluxos de trabalho serão configurados para publicar saídas em um diretório de resultados, e/ou os nomes de diretórios e estrutura podem ser diferentes.
    Um pouco mais adiante nesta seção, mostraremos como descobrir onde esse comportamento é especificado.

### 2.4. Encontre a saída original e os logs no diretório `work/`

Quando você executa um fluxo de trabalho, o Nextflow cria um 'diretório de tarefa' distinto para cada invocação individual de cada processo no fluxo de trabalho (=cada etapa no pipeline).
Para cada um, ele prepara as entradas necessárias, executa a(s) instrução(ões) relevante(s) e escreve saídas e arquivos de log dentro desse diretório, que é nomeado automaticamente usando um hash para torná-lo único.

Todos esses diretórios de tarefa ficarão sob um diretório chamado `work` dentro do seu diretório atual (onde você está executando o comando).

Isso pode parecer confuso, então vamos ver como isso se parece na prática.

Voltando à saída do console para o fluxo de trabalho que executamos anteriormente, tínhamos esta linha:

```console
[a3/1e1535] sayHello [100%] 1 of 1 ✔
```

Vê como a linha começa com `[a3/1e1535]`?
Essa é uma forma truncada do caminho do diretório de tarefa para aquela chamada de processo, e diz onde encontrar a saída da chamada de processo `sayHello` dentro do caminho do diretório `work/`.

Você pode encontrar o caminho completo digitando o seguinte comando (substituindo `a3/1e1535` pelo que você vê no seu próprio terminal) e pressionando a tecla tab para autocompletar o caminho ou adicionando um asterisco:

```bash
ls work/a3/1e1535*
```

Isso deve retornar o caminho completo do diretório: `work/a3/1e153543b0a7f9d2c4735ddb4ab231`

Vamos dar uma olhada no que está lá dentro.

??? abstract "Conteúdo do diretório"

    ```console
    work
    ├── a3
    │   └── 1e153543b0a7f9d2c4735ddb4ab231
    │       ├── .command.begin
    │       ├── .command.err
    │       ├── .command.log
    │       ├── .command.out
    │       ├── .command.run
    │       ├── .command.sh
    │       ├── .exitcode
    │       └── output.txt
    └── a4
        └── aa3694b8808bdcc1135ef4a1187a4d
            ├── .command.begin
            ├── .command.err
            ├── .command.log
            ├── .command.out
            ├── .command.run
            ├── .command.sh
            ├── .exitcode
            └── output.txt
    ```

??? question "Não está vendo a mesma coisa?"

    Os nomes exatos dos subdiretórios serão diferentes no seu sistema.

    Se você navegar pelo conteúdo do subdiretório de tarefa no explorador de arquivos do VSCode, verá todos os arquivos imediatamente.
    No entanto, os arquivos de log estão configurados para serem invisíveis no terminal, então se você quiser usar `ls` ou `tree` para visualizá-los, precisará definir a opção relevante para exibir arquivos invisíveis.

    ```bash
    tree -a work
    ```

Há dois conjuntos de diretórios em `work/`, das duas execuções de pipeline diferentes que fizemos.
Cada execução de tarefa obtém seu próprio diretório isolado para trabalhar.
Neste caso, o pipeline fez a mesma coisa ambas as vezes, então o conteúdo de cada diretório de tarefa é idêntico

Você deve reconhecer imediatamente o arquivo `output.txt`, que é de fato a saída original do processo `sayHello` que foi publicada no diretório `results`.
Se você abri-lo, encontrará a saudação `Hello World!` novamente.

```console title="work/a3/1e153543b0a7f9d2c4735ddb4ab231/output.txt"
Hello World!
```

E quanto a todos aqueles outros arquivos?

Esses são os arquivos auxiliares e de log que o Nextflow escreveu como parte da execução da tarefa:

- **`.command.begin`**: Arquivo sentinela criado assim que a tarefa é iniciada.
- **`.command.err`**: Mensagens de erro (`stderr`) emitidas pela chamada de processo
- **`.command.log`**: Saída de log completa emitida pela chamada de processo
- **`.command.out`**: Saída regular (`stdout`) pela chamada de processo
- **`.command.run`**: Script completo executado pelo Nextflow para executar a chamada de processo
- **`.command.sh`**: O comando que foi realmente executado pela chamada de processo
- **`.exitcode`**: O código de saída resultante do comando

O arquivo `.command.sh` é especialmente útil porque mostra o comando principal que o Nextflow executou, não incluindo toda a contabilidade e configuração de tarefa/ambiente.

```console title="work/a3/1e153543b0a7f9d2c4735ddb4ab231/.command.sh"
#!/bin/bash -ue
echo 'Hello World!' > output.txt

```

Isso confirma que o fluxo de trabalho compôs o mesmo comando que executamos diretamente na linha de comando anteriormente.

Quando algo dá errado e você precisa solucionar o que aconteceu, pode ser útil olhar o script `command.sh` para verificar exatamente qual comando o Nextflow compôs com base nas instruções do fluxo de trabalho, interpolação de variáveis e assim por diante.

### 2.5. Re-execute o fluxo de trabalho com diferentes saudações

Tente re-executar o fluxo de trabalho algumas vezes com diferentes valores para o argumento `--input`, depois olhe os diretórios de tarefa.

??? abstract "Conteúdo do diretório"

    ```console
    work/
    ├── 09
    │   └── 5ea8665939daf6f04724286c9b3c8a
    │       ├── .command.begin
    │       ├── .command.err
    │       ├── .command.log
    │       ├── .command.out
    │       ├── .command.run
    │       ├── .command.sh
    │       ├── .exitcode
    │       └── output.txt
    ├── 92
    │   └── ceb95e05d87621c92a399da9bd2067
    │       ├── .command.begin
    │       ├── .command.err
    │       ├── .command.log
    │       ├── .command.out
    │       ├── .command.run
    │       ├── .command.sh
    │       ├── .exitcode
    │       └── output.txt
    ├── 93
    │   └── 6708dbc20c7efdc6769cbe477061ec
    │       ├── .command.begin
    │       ├── .command.err
    │       ├── .command.log
    │       ├── .command.out
    │       ├── .command.run
    │       ├── .command.sh
    │       ├── .exitcode
    │       └── output.txt
    ├── a3
    │   └── 1e153543b0a7f9d2c4735ddb4ab231
    │       ├── .command.begin
    │       ├── .command.err
    │       ├── .command.log
    │       ├── .command.out
    │       ├── .command.run
    │       ├── .command.sh
    │       ├── .exitcode
    │       └── output.txt
    └── a4
        └── aa3694b8808bdcc1135ef4a1187a4d
            ├── .command.begin
            ├── .command.err
            ├── .command.log
            ├── .command.out
            ├── .command.run
            ├── .command.sh
            ├── .exitcode
            └── output.txt
    ```

Você vê que um novo subdiretório com um conjunto completo de arquivos de saída e log foi criado para cada execução.

Em contraste, se você olhar o diretório `results`, ainda há apenas um conjunto de resultados, e o conteúdo do arquivo de saída corresponde ao que você executou por último.

??? abstract "Conteúdo do diretório"

    ```console title="results/"
    results
    └── 1-hello
        └── output.txt
    ```

Isso mostra que os resultados publicados serão sobrescritos por execuções subsequentes, enquanto os diretórios de tarefa em `work/` são preservados.

### Conclusão

Você sabe como executar um script Nextflow simples, monitorar sua execução e encontrar suas saídas.

### O que vem a seguir?

Aprenda a ler um script Nextflow básico e identificar como seus componentes se relacionam com sua funcionalidade.

---

## 3. Examine o script inicial do fluxo de trabalho Hello World

O que fizemos ali foi basicamente tratar o script de fluxo de trabalho como uma caixa preta.
Agora que vimos o que ele faz, vamos abrir a caixa e olhar dentro.

Nosso objetivo aqui não é memorizar a sintaxe do código Nextflow, mas formar alguma intuição básica sobre quais são os principais componentes e como eles estão organizados.

### 3.1. Examine a estrutura geral do código

Você encontrará o script `1-hello.nf` no seu diretório atual, que deve ser `nextflow-run`. Abra-o no painel do editor.

??? full-code "Arquivo de código completo"

    ```groovy title="1-hello.nf" linenums="1"
    #!/usr/bin/env nextflow

    /*
    * Usa echo para imprimir 'Hello World!' em um arquivo
    */
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

    /*
    * Pipeline parameters
    */
    params {
        input: String
    }

    workflow {

        main:
        // emite uma saudação
        sayHello(params.input)

        publish:
        first_output = sayHello.out
    }

    output {
        first_output {
            path '1-hello'
            mode 'copy'
        }
    }
    ```

Um script de fluxo de trabalho Nextflow tipicamente inclui uma ou mais definições de **processo**, o **fluxo de trabalho** em si, e alguns blocos opcionais como **params** e **output**.

Cada **processo** descreve qual(is) operação(ões) a etapa correspondente no pipeline deve realizar, enquanto o **fluxo de trabalho** descreve a lógica de fluxo de dados que conecta as várias etapas.

Vamos dar uma olhada mais de perto no bloco **processo** primeiro, depois olharemos o bloco **workflow**.

### 3.2. A definição do `process`

O primeiro bloco de código descreve um [**processo**](https://nextflow.io/docs/latest/process.html).
A definição do processo começa com a palavra-chave `process`, seguida pelo nome do processo e finalmente o corpo do processo delimitado por chaves.
O corpo do processo deve conter um bloco script que especifica o comando a executar, que pode ser qualquer coisa que você seria capaz de executar em um terminal de linha de comando.

```groovy title="1-hello.nf" linenums="3"
/*
* Usa echo para imprimir uma saudação em um arquivo
*/
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

Aqui temos um **processo** chamado `sayHello` que recebe uma variável de **entrada** chamada `greeting` e escreve sua **saída** em um arquivo chamado `output.txt`.

<figure class="excalidraw">
--8<-- "docs/en/docs/nextflow_run/img/sayhello_with_input.svg"
</figure>

Esta é uma definição de processo muito mínima que contém apenas uma definição de `input`, uma definição de `output` e o `script` a executar.

A definição de `input` inclui o qualificador `val`, que diz ao Nextflow para esperar um valor de algum tipo (pode ser uma string, um número, qualquer coisa).

A definição de `output` inclui o qualificador `path`, que diz ao Nextflow que isso deve ser tratado como um caminho (inclui tanto caminhos de diretório quanto arquivos).

### 3.3. A definição do `workflow`

O segundo bloco de código descreve o próprio [**fluxo de trabalho**](https://nextflow.io/docs/latest/workflow.html).
A definição do fluxo de trabalho começa com a palavra-chave `workflow`, seguida por um nome opcional, depois o corpo do fluxo de trabalho delimitado por chaves.

Aqui temos um **fluxo de trabalho** que consiste em um bloco `main:` e um bloco `publish:`.
O bloco `main:` é o corpo principal do fluxo de trabalho e o bloco `publish:` lista as saídas que devem ser publicadas no diretório `results`.

```groovy title="1-hello.nf" linenums="27"
workflow {

    main:
    // emite uma saudação
    sayHello(params.input)

    publish:
    first_output = sayHello.out
}
```

Neste caso, o bloco `main:` contém uma chamada ao processo `sayHello` e fornece a ele uma entrada chamada `params.input` para usar como saudação.

Como discutiremos com mais detalhes em um momento, `params.input` contém o valor que demos ao parâmetro `--input` em nossa linha de comando.

O bloco `publish:` lista a saída da chamada de processo `sayHello()`, que ele se refere como `sayHello.out` e dá o nome `first_output` (isso pode ser qualquer coisa que o autor do fluxo de trabalho quiser).

Esta é uma definição de **fluxo de trabalho** muito mínima.
Em um pipeline do mundo real, o fluxo de trabalho tipicamente contém múltiplas chamadas a **processos** conectados por **canais**, e pode haver valores padrão configurados para as entradas variáveis.

Entraremos nisso na Parte 2 do curso.
Por agora, vamos dar uma olhada mais de perto em como nosso fluxo de trabalho está lidando com entradas e saídas.

### 3.4. O sistema `params` de parâmetros de linha de comando

O `params.input` que fornecemos à chamada do processo `sayHello()` é um pedaço elegante de código Nextflow e vale a pena gastar um minuto extra nele.

Como mencionado acima, é assim que passamos o valor do parâmetro de linha de comando `--input` para a chamada do processo `sayHello()`.
Na verdade, simplesmente declarar `params.someParameterName` é suficiente para dar ao fluxo de trabalho um parâmetro chamado `--someParameterName` da linha de comando.

Aqui formalizamos essa declaração de parâmetro configurando um bloco `params` que especifica o tipo de entrada que o fluxo de trabalho espera (Nextflow 25.10.2 e posterior).

```groovy title="1-hello.nf" linenums="20"
/*
 * Pipeline parameters
 */
params {
    input: String
}
```

Os tipos suportados incluem `String`, `Integer`, `Float`, `Boolean` e `Path`.
Para saber mais, veja [Workflow parameters](https://nextflow.io/docs/latest/config.html#workflow-parameters) na documentação de referência do Nextflow.

!!! tip "Dica"

    Lembre-se de que parâmetros de _fluxo de trabalho_ declarados usando o sistema `params` sempre levam dois traços na linha de comando (`--`).
    Isso os distingue de flags CLI de _nível Nextflow_, que levam apenas um traço (`-`).

### 3.5. A diretiva `publish`

Na outra ponta do fluxo de trabalho, já demos uma olhada no bloco `publish:`.
Essa é uma metade do sistema de tratamento de saída; a outra metade é o bloco `output` localizado abaixo.

```groovy title="1-hello.nf" linenums="37"
output {
    first_output {
        path '1-hello'
        mode 'copy'
    }
}
```

Isso especifica que a saída `first_output` listada no bloco `publish:` deve ser copiada para um subdiretório chamado `1-hello` sob o diretório de saída padrão `results`.

A linha `mode 'copy'` substitui o comportamento padrão do sistema, que é criar um link simbólico (ou symlink) para o arquivo original no diretório `work/` em vez de uma cópia apropriada.

Há mais opções do que as exibidas aqui para controlar o comportamento de publicação; cobriremos algumas mais tarde.
Você também verá que quando um fluxo de trabalho gera múltiplas saídas, cada uma é listada dessa forma no bloco `output`.

Para saber mais, veja [Publishing outputs](https://nextflow.io/docs/latest/workflow.html#publishing-outputs) na documentação de referência do Nextflow.

??? info "Sintaxe mais antiga para publicar saídas usando `publishDir`"

    Até muito recentemente, a forma estabelecida de publicar saídas era fazer isso no nível de cada processo individual usando uma diretiva `publishDir`.

    Você ainda encontrará esse padrão de código por todo lugar em pipelines e módulos de processo Nextflow mais antigos, então é importante estar ciente disso.

    Em vez de ter um bloco `publish:` no fluxo de trabalho e um bloco `output` no nível superior, você veria uma linha `publishDir` na definição do processo `sayHello`:

    ```groovy title="Exemplo de sintaxe" linenums="1" hl_lines="3"
    process sayHello {

        publishDir 'results/1-hello', mode: 'copy'

        output:
        path 'output.txt'

        script:
        """
        echo 'Hello World!' > output.txt
        """
    }
    ```

    No entanto, não recomendamos usar isso em nenhum trabalho novo, pois eventualmente será proibido em versões futuras da linguagem Nextflow.

### Conclusão

Agora você sabe como um fluxo de trabalho Nextflow simples é estruturado e como os componentes básicos se relacionam com sua funcionalidade.

### O que vem a seguir?

Aprenda a gerenciar suas execuções de fluxo de trabalho convenientemente.

---

## 4. Gerencie execuções de fluxo de trabalho

Saber como iniciar fluxos de trabalho e recuperar saídas é ótimo, mas você rapidamente descobrirá que há alguns outros aspectos do gerenciamento de fluxo de trabalho que facilitarão sua vida.

Aqui mostramos como aproveitar o recurso `resume` para quando você precisar relançar o mesmo fluxo de trabalho, como inspecionar os logs de execução com `nextflow log`, e como excluir diretórios de trabalho mais antigos com `nextflow clean`.

### 4.1. Relance um fluxo de trabalho com `-resume`

Às vezes, você vai querer re-executar um pipeline que já executou anteriormente sem refazer qualquer trabalho que já foi concluído com sucesso.

O Nextflow tem uma opção chamada `-resume` que permite fazer isso.
Especificamente, neste modo, quaisquer processos que já foram executados com exatamente o mesmo código, configurações e entradas serão pulados.
Isso significa que o Nextflow executará apenas processos que você adicionou ou modificou desde a última execução, ou aos quais você está fornecendo novas configurações ou entradas.

Há duas vantagens principais em fazer isso:

- Se você está no meio do desenvolvimento de um pipeline, pode iterar mais rapidamente, pois só precisa executar o(s) processo(s) em que está trabalhando ativamente para testar suas mudanças.
- Se você está executando um pipeline em produção e algo dá errado, em muitos casos você pode corrigir o problema e relançar o pipeline, e ele retomará a execução do ponto de falha, o que pode economizar muito tempo e computação.

Para usá-lo, simplesmente adicione `-resume` ao seu comando e execute:

```bash
nextflow run 1-hello.nf --input 'Hello World!' -resume
```

??? success "Saída do comando"

    ```console linenums="1"
    N E X T F L O W   ~  version 25.10.2

    Launching `1-hello.nf` [tiny_noyce] DSL2 - revision: c33d41f479

    [a3/7be2fa] sayHello | 1 of 1, cached: 1 ✔
    ```

A saída do console deve parecer familiar, mas há uma coisa um pouco diferente comparado a antes.

Procure pelo bit `cached:` que foi adicionado na linha de status do processo (linha 5), que significa que o Nextflow reconheceu que já fez este trabalho e simplesmente reutilizou o resultado da execução bem-sucedida anterior.

Você também pode ver que o hash do subdiretório de trabalho é o mesmo da execução anterior.
O Nextflow está literalmente apontando para a execução anterior e dizendo "Eu já fiz isso ali."

!!! tip "Dica"

    Quando você re-executa um pipeline com `resume`, o Nextflow não sobrescreve nenhum arquivo publicado fora do diretório de trabalho por quaisquer execuções que foram executadas com sucesso anteriormente.

    Para saber mais, veja [Cache and resume](https://nextflow.io/docs/latest/cache-and-resume.html) na documentação de referência do Nextflow.

### 4.2. Inspecione o log de execuções passadas

Sempre que você inicia um fluxo de trabalho Nextflow, uma linha é escrita em um arquivo de log chamado `history`, sob um diretório oculto chamado `.nextflow` no diretório de trabalho atual.

??? abstract "Conteúdo do arquivo"

    ```txt title=".nextflow/history" linenums="1"
    2025-07-04 19:27:09	1.8s	wise_watson	OK	3539118582ccde68dde471cc2c66295c	a02c9c46-c3c7-4085-9139-d1b9b5b194c8	nextflow run 1-hello.nf --input 'Hello World'
    2025-07-04 19:27:20	2.9s	spontaneous_blackwell	OK	3539118582ccde68dde471cc2c66295c	59a5db23-d83c-4c02-a54e-37ddb73a337e	nextflow run 1-hello.nf --input Bonjour
    2025-07-04 19:27:31	1.8s	gigantic_yonath	OK	3539118582ccde68dde471cc2c66295c	5acaa83a-6ad6-4509-bebc-cb25d5d7ddd0	nextflow run 1-hello.nf --input 'Dobry den'
    2025-07-04 19:27:45	2.4s	backstabbing_swartz	OK	3539118582ccde68dde471cc2c66295c	5f4b3269-5b53-404a-956c-cac915fbb74e	nextflow run 1-hello.nf --input Konnichiwa
    2025-07-04 19:27:57	2.1s	goofy_wilson	OK	3539118582ccde68dde471cc2c66295c	5f4b3269-5b53-404a-956c-cac915fbb74e	nextflow run 1-hello.nf --input Konnichiwa -resume
    ```

Este arquivo dá a você o timestamp, nome da execução, status, ID de revisão, ID de sessão e linha de comando completa para cada execução Nextflow que foi iniciada dentro do diretório de trabalho atual.

Uma forma mais conveniente de acessar esta informação é usar o comando [`nextflow log`](https://nextflow.io/docs/latest/reference/cli.html#log).

```bash
nextflow log
```

??? success "Saída do comando"

    ```console linenums="1"
    TIMESTAMP               DURATION        RUN NAME                STATUS  REVISION ID     SESSION ID                              COMMAND
    2025-07-04 19:27:09     1.8s            wise_watson             OK       3539118582     a02c9c46-c3c7-4085-9139-d1b9b5b194c8    nextflow run 1-hello.nf --input 'Hello World'
    2025-07-04 19:27:20     2.9s            spontaneous_blackwell   OK       3539118582     59a5db23-d83c-4c02-a54e-37ddb73a337e    nextflow run 1-hello.nf --input Bonjour
    2025-07-04 19:27:31     1.8s            gigantic_yonath         OK       3539118582     5acaa83a-6ad6-4509-bebc-cb25d5d7ddd0    nextflow run 1-hello.nf --input 'Dobry den'
    2025-07-04 19:27:45     2.4s            backstabbing_swartz     OK       3539118582     5f4b3269-5b53-404a-956c-cac915fbb74e    nextflow run 1-hello.nf --input Konnichiwa
    2025-07-04 19:27:57     2.1s            goofy_wilson            OK       3539118582     5f4b3269-5b53-404a-956c-cac915fbb74e    nextflow run 1-hello.nf --input Konnichiwa -resume
    ```

Isso exibirá o conteúdo do arquivo de log no terminal, aumentado com uma linha de cabeçalho.

Você notará que o ID de sessão muda sempre que você executa um novo comando `nextflow run`, EXCETO se você está usando a opção `-resume`.
Nesse caso, o ID de sessão permanece o mesmo.

O Nextflow usa o ID de sessão para agrupar informações de cache de execução sob o diretório `cache`, também localizado em `.nextflow`.

### 4.3. Exclua diretórios de trabalho mais antigos

Se você executar muitos pipelines, pode acabar acumulando muitos arquivos em muitos subdiretórios.
Como os subdiretórios são nomeados aleatoriamente, é difícil dizer pelos nomes quais são execuções mais antigas vs. mais recentes.

Felizmente o Nextflow inclui um comando útil chamado [`nextflow clean`](https://www.nextflow.io/docs/latest/reference/cli.html#clean) que pode excluir automaticamente os subdiretórios de trabalho para execuções passadas que você não se importa mais.

#### 4.3.1. Determine os critérios de exclusão

Há múltiplas opções para determinar o que excluir, que você pode explorar na documentação vinculada acima.
Aqui mostramos um exemplo que exclui todos os subdiretórios de execuções antes de uma determinada execução, especificada usando seu nome de execução.

Procure a execução bem-sucedida mais recente onde você não usou `-resume`; no nosso caso o nome da execução era `backstabbing_swartz`.

O nome da execução é a string de duas partes gerada pela máquina mostrada entre colchetes na linha de saída do console `Launching (...)`.
Você também pode usar o log do Nextflow para procurar uma execução com base em seu timestamp e/ou linha de comando.

#### 4.3.2. Faça uma execução de teste

Primeiro usamos a flag de execução de teste `-n` para verificar o que será excluído dado o comando:

```bash
nextflow clean -before backstabbing_swartz -n
```

??? success "Saída do comando"

    ```console
    Would remove /workspaces/training/hello-nextflow/work/eb/1a5de36637b475afd88fca7f79e024
    Would remove /workspaces/training/hello-nextflow/work/6b/19b0e002ea13486d3a0344c336c1d0
    Would remove /workspaces/training/hello-nextflow/work/45/9a6dd7ab771f93003d040956282883
    ```

Sua saída terá nomes de diretórios de tarefa diferentes e pode ter um número diferente de linhas, mas deve parecer similar ao exemplo.

Se você não vir nenhuma linha de saída, você ou não forneceu um nome de execução válido ou não há execuções passadas para excluir. Certifique-se de mudar `backstabbing_swartz` no comando de exemplo para qualquer que seja o nome de execução mais recente correspondente no seu log.

#### 4.3.3. Prossiga com a exclusão

Se a saída parecer como esperado e você quiser prosseguir com a exclusão, re-execute o comando com a flag `-f` em vez de `-n`:

```bash
nextflow clean -before backstabbing_swartz -f
```

??? success "Saída do comando"

    ```console
    Removed /workspaces/training/hello-nextflow/work/eb/1a5de36637b475afd88fca7f79e024
    Removed /workspaces/training/hello-nextflow/work/6b/19b0e002ea13486d3a0344c336c1d0
    Removed /workspaces/training/hello-nextflow/work/45/9a6dd7ab771f93003d040956282883
    ```

A saída deve ser similar a antes, mas agora dizendo 'Removed' em vez de 'Would remove'.
Note que isso não remove os subdiretórios de dois caracteres (como `eb/` acima), mas esvazia seu conteúdo.

!!! Warning "Aviso"

    Excluir subdiretórios de trabalho de execuções passadas os remove do cache do Nextflow e exclui quaisquer saídas que estavam armazenadas nesses diretórios.
    Isso significa que quebra a capacidade do Nextflow de retomar a execução sem re-executar os processos correspondentes.

    Você é responsável por salvar quaisquer saídas que você se importe! Essa é a principal razão pela qual preferimos usar o modo `copy` em vez do modo `symlink` para a diretiva `publish`.

### Conclusão

Você sabe como relançar um pipeline sem repetir etapas que já foram executadas de forma idêntica, inspecionar o log de execução e usar o comando `nextflow clean` para limpar diretórios de trabalho antigos.

### O que vem a seguir?

Faça uma pequena pausa! Você acabou de absorver os blocos de construção da sintaxe Nextflow e instruções básicas de uso.

Na próxima seção deste treinamento, vamos olhar quatro versões sucessivamente mais realistas do pipeline Hello World que demonstrarão como o Nextflow permite processar múltiplas entradas eficientemente, executar fluxos de trabalho compostos de múltiplas etapas conectadas, aproveitar componentes de código modulares e utilizar contêineres para maior reprodutibilidade e portabilidade.

---

## Quiz

<quiz>
Na linha de saída do console `[a3/7be2fa] SAYHELLO | 1 of 1 ✔`, o que `[a3/7be2fa]` representa?
- [ ] O número da versão do processo
- [ ] Um identificador único de execução
- [x] O caminho truncado para o diretório de trabalho da tarefa
- [ ] O checksum do arquivo de saída

Saiba mais: [2.4. Encontre a saída original e os logs no diretório `work/`](#24-encontre-a-saida-original-e-os-logs-no-diretorio-work)
</quiz>

<quiz>
Qual é o propósito do arquivo `.command.sh` em um diretório de tarefa?
- [ ] Ele armazena as configurações da tarefa
- [x] Ele mostra o comando real que foi executado pelo processo
- [ ] Ele contém mensagens de erro de tarefas que falharam
- [ ] Ele lista os arquivos de entrada preparados para a tarefa

Saiba mais: [2.4. Encontre a saída original e os logs no diretório `work/`](#24-encontre-a-saida-original-e-os-logs-no-diretorio-work)
</quiz>

<quiz>
O que acontece com os resultados publicados quando você re-executa um fluxo de trabalho sem `-resume`?
- [ ] Eles são preservados em diretórios separados com timestamp
- [x] Eles são sobrescritos pela nova execução
- [ ] O Nextflow impede a sobrescrita e falha
- [ ] Eles são automaticamente copiados como backup

Saiba mais: [2.5. Re-execute o fluxo de trabalho com diferentes saudações](#25-re-execute-o-fluxo-de-trabalho-com-diferentes-saudacoes)
</quiz>

<quiz>
O que esta saída do console indica?

```console
[skipped  ] process > sayHello (1) [100%] 1 of 1, cached: 1 ✔
```

- [ ] A tarefa falhou e foi pulada
- [ ] A tarefa está esperando em uma fila
- [x] O Nextflow reutilizou resultados de uma execução idêntica anterior
- [ ] A tarefa foi cancelada manualmente

Saiba mais: [4.1. Relance um fluxo de trabalho com `-resume`](#41-relance-um-fluxo-de-trabalho-com--resume)
</quiz>

<quiz>
Onde o Nextflow armazena o histórico de execução que o comando `nextflow log` exibe?
- [ ] No diretório results
- [ ] No diretório work
- [x] No arquivo `.nextflow/history`
- [ ] No `nextflow.config`

Saiba mais: [4.2. Inspecione o log de execuções passadas](#42-inspecione-o-log-de-execucoes-passadas)
</quiz>

<quiz>
Qual é o propósito do bloco `params` em um arquivo de fluxo de trabalho?
- [ ] Definir requisitos de recursos do processo
- [ ] Configurar o executor
- [x] Declarar e tipar parâmetros de entrada do fluxo de trabalho
- [ ] Especificar opções de publicação de saída

Saiba mais: [3.4. O sistema params de parâmetros de linha de comando](#34-o-sistema-params-de-parametros-de-linha-de-comando)
</quiz>

<quiz>
No bloco `output` do fluxo de trabalho, o que `mode 'copy'` faz?
- [ ] Cria um backup do diretório de trabalho
- [x] Faz uma cópia completa dos arquivos em vez de links simbólicos
- [ ] Copia o script do fluxo de trabalho para results
- [ ] Habilita cópia incremental de arquivos

Saiba mais: [3.5. A diretiva publish](#35-a-diretiva-publish)
</quiz>

<quiz>
Qual é a flag recomendada para usar com o comando `nextflow clean` antes de realmente excluir arquivos?
- [x] `-n` (dry run) para visualizar o que seria excluído
- [ ] `-v` (verbose) para ver saída detalhada
- [ ] `-a` (all) para selecionar todos os diretórios
- [ ] `-q` (quiet) para suprimir avisos

Saiba mais: [4.3. Exclua diretórios de trabalho mais antigos](#43-exclua-diretorios-de-trabalho-mais-antigos)
</quiz>
