# Parte 1: Executar operações básicas

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tradução assistida por IA - [saiba mais e sugira melhorias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Nesta primeira parte do curso de treinamento Nextflow para Bioimagem, usaremos um exemplo básico Hello World agnóstico de domínio para demonstrar operações essenciais e apontar os componentes de código Nextflow correspondentes.

## 1. Executar o fluxo de trabalho

Fornecemos um script de fluxo de trabalho chamado `hello-world.nf` que recebe uma entrada via argumento de linha de comando chamado `--greeting` e produz um arquivo de texto contendo essa saudação.
Ainda não vamos olhar o código; primeiro vamos ver como é executá-lo.

### 1.1. Iniciar o fluxo de trabalho e monitorar a execução

No terminal, execute o seguinte comando:

```bash
nextflow run hello-world.nf --greeting 'Hello World!'
```

A saída do console deve se parecer com isto:

```console title="Saída" linenums="1"
 N E X T F L O W   ~  version 25.04.3

Launching `hello-world.nf` [goofy_torvalds] DSL2 - revision: c33d41f479

executor >  local (1)
[a3/7be2fa] sayHello | 1 of 1 ✔
```

Parabéns, você acabou de executar seu primeiro fluxo de trabalho Nextflow!

A saída mais importante aqui é a última linha (linha 6):

```console title="Saída" linenums="6"
[a3/7be2fa] sayHello | 1 of 1 ✔
```

Isso nos diz que o processo `sayHello` foi executado com sucesso uma vez (`1 of 1 ✔`).

Isso é ótimo, mas você pode estar se perguntando: onde está a saída?

### 1.2. Encontrar o arquivo de saída no diretório `results`

Este fluxo de trabalho está configurado para publicar sua saída em um diretório chamado `results`.
Se você olhar seu diretório atual, verá que quando você executou o fluxo de trabalho, o Nextflow criou um novo diretório chamado `results`, que contém um arquivo chamado `output.txt`.

```console title="results/" linenums="1"
results
└── output.txt
```

Abra o arquivo; o conteúdo deve corresponder à saudação que você especificou na linha de comando.

<details>
  <summary>Conteúdo do arquivo</summary>

```console title="results/output.txt" linenums="1"
Hello World!
```

</details>

Isso é ótimo, nosso fluxo de trabalho fez o que deveria fazer!

No entanto, esteja ciente de que o resultado 'publicado' é uma cópia (ou em alguns casos um symlink) da saída real produzida pelo Nextflow quando executou o fluxo de trabalho.

Então agora, vamos olhar por baixo do capô para ver onde o Nextflow realmente executou o trabalho.

!!! warning "Aviso"

    Nem todos os fluxos de trabalho estarão configurados para publicar saídas em um diretório results, e/ou o nome do diretório pode ser diferente.
    Um pouco mais adiante nesta seção, mostraremos como descobrir onde esse comportamento é especificado.

### 1.3. Encontrar a saída original e logs no diretório `work/`

Quando você executa um fluxo de trabalho, o Nextflow cria um 'diretório de tarefa' distinto para cada invocação de cada processo no fluxo de trabalho (=cada etapa no pipeline).
Para cada um, ele preparará as entradas necessárias, executará a(s) instrução(ões) relevante(s) e gravará saídas e arquivos de log dentro daquele único diretório, que é nomeado automaticamente usando um hash para torná-lo único.

Todos esses diretórios de tarefa ficarão em um diretório chamado `work` dentro do seu diretório atual (onde você está executando o comando).

Isso pode parecer confuso, então vamos ver como isso se parece na prática.

Voltando à saída do console para o fluxo de trabalho que executamos anteriormente, tínhamos esta linha:

```console title="Trecho da saída do comando" linenums="6"
[a3/7be2fa] sayHello | 1 of 1 ✔
```

Vê como a linha começa com `[a3/7be2fa]`?
Essa é uma forma truncada do caminho do diretório de tarefa para aquela chamada de processo, e diz onde encontrar a saída da chamada do processo `sayHello` dentro do caminho do diretório `work/`.

Você pode encontrar o caminho completo digitando o seguinte comando (substituindo `a3/7be2fa` pelo que você vê em seu próprio terminal) e pressionando a tecla tab para autocompletar o caminho ou adicionando um asterisco:

```bash
tree work/a3/7be2fa*
```

Isso deve gerar o caminho completo do diretório: `work/a3/7be2fa7be2fad5e71e5f49998f795677fd68`

Vamos dar uma olhada no que há lá dentro.

!!! Tip "Dica"

    Se você navegar pelo conteúdo do subdiretório de tarefa no explorador de arquivos do VSCode, verá todos os arquivos imediatamente.
    No entanto, os arquivos de log estão configurados para serem invisíveis no terminal, então se você quiser usar `ls` ou `tree` para visualizá-los, precisará definir a opção relevante para exibir arquivos invisíveis.

    ```bash
    tree -a work
    ```

Os nomes exatos dos subdiretórios serão diferentes no seu sistema.

<details>
  <summary>Conteúdo do diretório</summary>

```console title="work/"
work
└── a3
    └── 7be2fad5e71e5f49998f795677fd68
        ├── .command.begin
        ├── .command.err
        ├── .command.log
        ├── .command.out
        ├── .command.run
        ├── .command.sh
        ├── .exitcode
        └── output.txt
```

</details>

Você deve reconhecer imediatamente o arquivo `output.txt`, que é de fato a saída original do processo `sayHello` que foi publicada no diretório `results`.
Se você abri-lo, encontrará a saudação `Hello World!` novamente.

<details>
  <summary>Conteúdo do arquivo output.txt</summary>

```console title="work/a3/7be2fa7be2fad5e71e5f49998f795677fd68/output.txt" linenums="1"
Hello World!
```

</details>

E quanto a todos aqueles outros arquivos?

Estes são os arquivos auxiliares e de log que o Nextflow escreveu como parte da execução da tarefa:

- **`.command.begin`**: Arquivo sentinela criado assim que a tarefa é lançada.
- **`.command.err`**: Mensagens de erro (`stderr`) emitidas pela chamada do processo
- **`.command.log`**: Saída de log completa emitida pela chamada do processo
- **`.command.out`**: Saída regular (`stdout`) pela chamada do processo
- **`.command.run`**: Script completo executado pelo Nextflow para executar a chamada do processo
- **`.command.sh`**: O comando que foi realmente executado pela chamada do processo
- **`.exitcode`**: O código de saída resultante do comando

O arquivo `.command.sh` é especialmente útil porque mostra o comando principal que o Nextflow executou, não incluindo toda a contabilidade e configuração de tarefa/ambiente.

<details>
  <summary>Conteúdo do arquivo</summary>

```console title="work/a3/7be2fa7be2fad5e71e5f49998f795677fd68/.command.sh" linenums="1"
#!/bin/bash -ue
echo 'Hello World!' > output.txt

```

</details>

!!! Tip "Dica"

    Quando algo dá errado e você precisa solucionar o que aconteceu, pode ser útil olhar o script `command.sh` para verificar exatamente qual comando o Nextflow compôs com base nas instruções do fluxo de trabalho, interpolação de variáveis e assim por diante.

### 1.4. Exercício opcional: executar novamente com saudações diferentes

Tente executar novamente o fluxo de trabalho algumas vezes com valores diferentes para o argumento `--greeting`, depois olhe tanto o conteúdo do diretório `results/` quanto os diretórios de tarefa.

Observe como as saídas e logs de diretórios de tarefa isolados são preservados, enquanto o conteúdo do diretório `results` é sobrescrito pela saída de execuções subsequentes.

### Conclusão

Você sabe como executar um script Nextflow simples, monitorar sua execução e encontrar suas saídas.

### O que vem a seguir?

Aprenda a ler um script Nextflow básico e identificar como seus componentes se relacionam com sua funcionalidade.

---

## 2. Examinar o script inicial do fluxo de trabalho Hello World

O que fizemos lá foi basicamente tratar o script do fluxo de trabalho como uma caixa preta.
Agora que vimos o que ele faz, vamos abrir a caixa e olhar dentro.

_O objetivo aqui não é memorizar a sintaxe do código Nextflow, mas formar alguma intuição básica sobre quais são os componentes principais e como eles são organizados._

### 2.1. Examinar a estrutura geral do código

Vamos abrir o script `hello-world.nf` no painel do editor.

<details>
  <summary>Código</summary>

```groovy title="hello-world.nf" linenums="1"
#!/usr/bin/env nextflow

/*
 * Usa echo para imprimir uma saudação em um arquivo
 */
process sayHello {

    publishDir 'results', mode: 'copy'

    input:
    val greeting

    output:
    path 'output.txt'

    script:
    """
    echo '$greeting' > output.txt
    """
}

workflow {

    // emite uma saudação
    sayHello(params.greeting)
}
```

</details>

Um script Nextflow envolve dois tipos principais de componentes centrais: um ou mais **processos**, e o **fluxo de trabalho** em si.
Cada **processo** descreve quais operação(ões) a etapa correspondente no pipeline deve realizar, enquanto o **fluxo de trabalho** descreve a lógica de fluxo de dados que conecta as várias etapas.

Vamos dar uma olhada mais de perto no bloco **process** primeiro, depois veremos o bloco **workflow**.

### 2.2. A definição de `process`

O primeiro bloco de código descreve um **processo**.
A definição do processo começa com a palavra-chave `process`, seguida pelo nome do processo e finalmente o corpo do processo delimitado por chaves.
O corpo do processo deve conter um bloco script que especifica o comando a ser executado, que pode ser qualquer coisa que você seria capaz de executar em um terminal de linha de comando.

Aqui temos um **processo** chamado `sayHello` que recebe uma variável de **entrada** chamada `greeting` e escreve sua **saída** em um arquivo chamado `output.txt`.

<details>
  <summary>Código</summary>

```groovy title="hello-world.nf" linenums="3"
/*
 * Usa echo para imprimir uma saudação em um arquivo
 */
process sayHello {

    publishDir 'results', mode: 'copy'

    input:
    val greeting

    output:
    path 'output.txt'

    script:
    """
    echo '$greeting' > output.txt
    """
}
```

</details>

Esta é uma definição de processo muito mínima que contém apenas uma definição de `input`, uma definição de `output` e o `script` a executar.

A definição de `input` inclui o qualificador `val`, que diz ao Nextflow para esperar um valor de algum tipo (pode ser uma string, um número, qualquer coisa).

A definição de `output` inclui o qualificador `path`, que diz ao Nextflow que isso deve ser tratado como um caminho (inclui tanto caminhos de diretório quanto arquivos).

!!! Tip "Dica"

    A definição de saída não _determina_ qual saída será criada.
    Ela simplesmente _declara_ onde encontrar o(s) arquivo(s) de saída esperado(s), para que o Nextflow possa procurá-lo uma vez que a execução esteja completa.

    Isso é necessário para verificar se o comando foi executado com sucesso e para passar a saída para processos posteriores, se necessário.
    A saída produzida que não corresponder ao que é declarado no bloco de saída não será passada para processos posteriores.

Em um pipeline do mundo real, um processo geralmente contém informações adicionais, como diretivas de processo, que introduziremos daqui a pouco.

### 2.3. A definição de `workflow`

O segundo bloco de código descreve o **fluxo de trabalho** em si.
A definição do fluxo de trabalho começa com a palavra-chave `workflow`, seguida por um nome opcional, depois o corpo do fluxo de trabalho delimitado por chaves.

Aqui temos um **fluxo de trabalho** que consiste em uma chamada ao processo `sayHello`, que recebe uma entrada, `params.greeting`, que contém o valor que demos ao parâmetro `--greeting`.

```groovy title="hello-world.nf" linenums="22"
workflow {

    // emite uma saudação
    sayHello(params.greeting)
}
```

Esta é uma definição de **fluxo de trabalho** muito mínima.
Em um pipeline do mundo real, o fluxo de trabalho normalmente contém múltiplas chamadas a **processos** conectados por **canais**, e pode haver valores padrão configurados para as entradas de variáveis.

Veremos isso em ação quando executarmos nf-core/molkart na Parte 2 do curso.

### 2.4. O sistema `params` de parâmetros de linha de comando

O `params.greeting` que fornecemos à chamada do processo `sayHello()` é um pedaço interessante de código Nextflow e vale a pena gastar um minuto extra nele.

Como mencionado acima, é assim que passamos o valor do parâmetro de linha de comando `--greeting` para a chamada do processo `sayHello()`.
Na verdade, simplesmente declarar `params.someParameterName` nos permitirá dar ao fluxo de trabalho um parâmetro chamado `--someParameterName` a partir da linha de comando.

!!! Tip "Dica"

    Esses parâmetros de fluxo de trabalho declarados usando o sistema `params` sempre levam dois traços (`--`).
    Isso os distingue dos parâmetros de nível Nextflow, que levam apenas um traço (`-`).

### Conclusão

Você agora sabe como um fluxo de trabalho Nextflow simples é estruturado, e como os componentes básicos se relacionam com sua funcionalidade.

### O que vem a seguir?

Aprenda a gerenciar suas execuções de fluxo de trabalho de forma conveniente.

---

## 3. Gerenciar execuções de fluxo de trabalho

Saber como iniciar fluxos de trabalho e recuperar saídas é ótimo, mas você rapidamente descobrirá que existem alguns outros aspectos do gerenciamento de fluxo de trabalho que tornarão sua vida mais fácil.

Aqui mostramos como aproveitar o recurso `resume` para quando você precisar relançar o mesmo fluxo de trabalho, como inspecionar os logs de execução com `nextflow log`, e como excluir diretórios de trabalho mais antigos com `nextflow clean`.

### 3.1. Relançar um fluxo de trabalho com `-resume`

Às vezes, você vai querer executar novamente um pipeline que você já lançou anteriormente sem refazer nenhum trabalho que já foi concluído com sucesso.

O Nextflow tem uma opção chamada `-resume` que permite fazer isso.
Especificamente, neste modo, quaisquer processos que já foram executados com exatamente o mesmo código, configurações e entradas serão ignorados.
Isso significa que o Nextflow só executará processos que você adicionou ou modificou desde a última execução, ou aos quais você está fornecendo novas configurações ou entradas.

Existem duas vantagens principais em fazer isso:

- Se você está no meio do desenvolvimento de um pipeline, pode iterar mais rapidamente, já que só precisa executar o(s) processo(s) em que está trabalhando ativamente para testar suas alterações.
- Se você está executando um pipeline em produção e algo dá errado, em muitos casos você pode corrigir o problema e relançar o pipeline, e ele retomará a execução do ponto de falha, o que pode economizar muito tempo e computação.

Para usá-lo, simplesmente adicione `-resume` ao seu comando e execute-o:

```bash
nextflow run hello-world.nf --greeting 'Hello World!' -resume
```

??? success "Saída do comando"

    ```console
    N E X T F L O W   ~  version 25.04.3

    Launching `hello-world.nf` [tiny_noyce] DSL2 - revision: c33d41f479

    [a3/7be2fa] process > sayHello [100%] 1 of 1, cached: 1 ✔
    ```

Procure pelo trecho `cached:` que foi adicionado na linha de status do processo (linha 5), o que significa que o Nextflow reconheceu que já fez este trabalho e simplesmente reutilizou o resultado da execução bem-sucedida anterior.

Você também pode ver que o hash do subdiretório de trabalho é o mesmo da execução anterior.
O Nextflow está literalmente apontando para a execução anterior e dizendo "Eu já fiz isso ali."

!!! Tip "Dica"

    Quando você executa novamente um pipeline com `resume`, o Nextflow não sobrescreve nenhum arquivo escrito em um diretório `publishDir` por qualquer chamada de processo que foi executada com sucesso anteriormente.

### 3.2. Inspecionar o log de execuções passadas

Sempre que você inicia um fluxo de trabalho nextflow, uma linha é escrita em um arquivo de log chamado `history`, em um diretório oculto chamado `.nextflow` no diretório de trabalho atual.

Uma maneira mais conveniente de acessar essas informações é usar o comando `nextflow log`.

```bash
nextflow log
```

Isso exibirá o conteúdo do arquivo de log no terminal, mostrando o timestamp, nome da execução, status e linha de comando completa para cada execução Nextflow que foi lançada de dentro do diretório de trabalho atual.

### 3.3. Excluir diretórios de trabalho mais antigos

Durante o processo de desenvolvimento, você normalmente executará seus rascunhos de pipelines um grande número de vezes, o que pode levar a uma acumulação de muitos arquivos em muitos subdiretórios.
Como os subdiretórios são nomeados aleatoriamente, é difícil dizer pelos nomes quais são execuções mais antigas vs. mais recentes.

O Nextflow inclui um subcomando `clean` conveniente que pode excluir automaticamente os subdiretórios de trabalho de execuções passadas que você não se importa mais, com várias [opções](https://www.nextflow.io/docs/latest/reference/cli.html#clean) para controlar o que será excluído.

Você pode usar o log do Nextflow para procurar uma execução com base em seu timestamp e/ou linha de comando, depois usar `nextflow clean -before <run_name> -f` para excluir diretórios de trabalho de execuções anteriores.

!!! Warning "Aviso"

    Excluir subdiretórios de trabalho de execuções passadas os remove do cache do Nextflow e exclui quaisquer saídas que estavam armazenadas nesses diretórios.
    Isso significa que quebra a capacidade do Nextflow de retomar a execução sem executar novamente os processos correspondentes.

    Você é responsável por salvar quaisquer saídas que você se importa ou planeja usar! Se você está usando a diretiva `publishDir` para esse propósito, certifique-se de usar o modo `copy`, não o modo `symlink`.

### Conclusão

Você sabe como relançar um pipeline sem repetir etapas que já foram executadas de forma idêntica, inspecionar o log de execução e usar o comando `nextflow clean` para limpar diretórios de trabalho antigos.

### O que vem a seguir?

Agora que você entende operações básicas do Nextflow, está pronto para executar um pipeline de bioimagem real com nf-core/molkart.
