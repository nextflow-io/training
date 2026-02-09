# Testando com nf-test

Ser capaz de testar sistematicamente que cada parte do seu fluxo de trabalho está fazendo o que deveria fazer é fundamental para reprodutibilidade e manutenção a longo prazo, e pode ser uma grande ajuda durante o processo de desenvolvimento.

Vamos dedicar um minuto para falar sobre por que testar é tão importante. Se você está desenvolvendo um fluxo de trabalho, uma das primeiras coisas que você fará é pegar alguns dados de teste que você sabe que são válidos e devem produzir um resultado. Você adiciona o primeiro processo ao pipeline e o conecta às suas entradas para fazê-lo funcionar. Então, para verificar se está tudo funcionando, você o executa com os dados de teste. Assumindo que funcione, você passa para o próximo processo e executa os dados de teste novamente. Você repete esse processo até ter um pipeline com o qual está satisfeito.

Então, talvez você adicione um parâmetro simples verdadeiro ou falso como `--skip_process`. Agora você deve executar o pipeline duas vezes, uma com cada parâmetro para garantir que funcione como esperado. Mas espere, como verificamos se o `--skip_process` realmente pula o processo? Temos que investigar as saídas ou verificar os arquivos de log! Isso é trabalhoso e propenso a erros.

À medida que você desenvolve seu pipeline, ele rapidamente se tornará tão complexo que testar manualmente cada iteração é lento e propenso a erros. Além disso, se você encontrar um erro, será muito difícil identificar exatamente de onde no seu pipeline o erro está vindo. É aqui que os testes entram.

Testar permite que você verifique sistematicamente que cada parte do seu pipeline está funcionando como esperado. Os benefícios para um desenvolvedor de testes bem escritos são enormes:

- **Confiança**: Como os testes cobrem todo o pipeline, você pode ter confiança de que mudar algo não afeta nada mais
- **Confiabilidade**: Quando múltiplos desenvolvedores trabalham no pipeline, eles sabem que os outros desenvolvedores não quebraram o pipeline e cada componente.
- **Transparência**: Os testes mostram onde um pipeline está falando e facilitam rastrear o problema. Eles também funcionam como uma forma de documentação, mostrando como executar um processo ou fluxo de trabalho.
- **Velocidade**: Como os testes são automatizados, eles podem ser executados muito rapidamente e repetidamente. Você pode iterar rapidamente com menos medo de introduzir novos bugs.

Existem muitos tipos diferentes de testes que podemos escrever:

1. **Testes em nível de módulo**: Para processos individuais
2. **Testes em nível de fluxo de trabalho**: Para um único fluxo de trabalho
3. **Testes em nível de pipeline**: Para o pipeline como um todo
4. **Testes de desempenho**: Para a velocidade e eficiência do pipeline
5. **Testes de estresse**: Avaliando o desempenho do pipeline sob condições extremas para determinar seus limites

Testar processos individuais é análogo a testes unitários em outras linguagens. Testar o fluxo de trabalho ou o pipeline inteiro é análogo ao que é chamado de testes de integração em outras linguagens, onde testamos as interações dos componentes.

[**nf-test**](https://www.nf-test.com/) é uma ferramenta que permite escrever testes em nível de módulo, fluxo de trabalho e pipeline. Em resumo, permite que você verifique sistematicamente que cada parte individual do pipeline está funcionando como esperado, _isoladamente_.

### Objetivos de aprendizado

Nesta missão secundária, você aprenderá a usar nf-test para escrever um teste em nível de fluxo de trabalho para o pipeline, bem como testes em nível de módulo para os três processos que ele chama.

Ao final desta missão secundária, você será capaz de usar as seguintes técnicas efetivamente:

- Inicializar nf-test no seu projeto
- Gerar testes em nível de módulo e fluxo de trabalho
- Adicionar tipos comuns de asserções
- Entender quando usar snapshots vs. asserções de conteúdo
- Executar testes para um projeto inteiro

Essas habilidades ajudarão você a implementar uma estratégia de testes abrangente em seus projetos de pipeline, garantindo que sejam mais robustos e fáceis de manter.

### Pré-requisitos

Antes de assumir esta missão secundária, você deve:

- Ter completado o tutorial [Hello Nextflow](../hello_nextflow/README.md) ou curso equivalente para iniciantes.
- Estar confortável usando conceitos e mecanismos básicos do Nextflow (processos, canais, operadores, trabalhando com arquivos, metadados)

---

## 0. Primeiros passos

#### Abra o codespace de treinamento

Se você ainda não o fez, certifique-se de abrir o ambiente de treinamento conforme descrito na [Configuração do Ambiente](../envsetup/index.md).

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### Mova-se para o diretório do projeto

Vamos mover para o diretório onde os arquivos para este tutorial estão localizados.

```bash
cd side-quests/nf-test
```

Você pode configurar o VSCode para focar neste diretório:

```bash
code .
```

#### Revise os materiais

Você encontrará um arquivo de fluxo de trabalho principal e um arquivo CSV chamado `greetings.csv` que contém a entrada para o pipeline.

```console title="Directory contents"
.
├── greetings.csv
└── main.nf
```

Para uma descrição detalhada dos arquivos, consulte o [aquecimento do Hello Nextflow](../hello_nextflow/00_orientation.md).

O fluxo de trabalho que testaremos é um subconjunto do fluxo de trabalho Hello construído em [Hello Workflow](../hello_nextflow/03_hello_workflow.md).

??? example "O que o fluxo de trabalho Hello Nextflow faz?"

    Se você não fez o treinamento [Hello Nextflow](../hello_nextflow/index.md), aqui está uma visão geral rápida do que este fluxo de trabalho simples faz.

    O fluxo de trabalho pega um arquivo CSV contendo saudações, executa quatro etapas consecutivas de transformação nelas e gera um único arquivo de texto contendo uma imagem ASCII de um personagem divertido dizendo as saudações.

    As quatro etapas são implementadas como processos Nextflow (`sayHello`, `convertToUpper`, `collectGreetings` e `cowpy`) armazenados em arquivos de módulo separados.

    1. **`sayHello`:** Escreve cada saudação em seu próprio arquivo de saída (por exemplo, "Hello-output.txt")
    2. **`convertToUpper`:** Converte cada saudação para maiúsculas (por exemplo, "HELLO")
    3. **`collectGreetings`:** Coleta todas as saudações em maiúsculas em um único arquivo em lote
    4. **`cowpy`:** Gera arte ASCII usando a ferramenta `cowpy`

    Os resultados são publicados em um diretório chamado `results/`, e a saída final do pipeline (quando executado com parâmetros padrão) é um arquivo de texto simples contendo arte ASCII de um personagem dizendo as saudações em maiúsculas.

    Nesta missão secundária, usamos uma forma intermediária do fluxo de trabalho Hello que contém apenas os dois primeiros processos. <!-- TODO: change this to use the full finished workflow as suggested in https://github.com/nextflow-io/training/issues/735 -->

O subconjunto com o qual trabalharemos é composto por dois processos: `sayHello` e `convertToUpper`.
Você pode ver o código completo do fluxo de trabalho abaixo.

??? example "Código do fluxo de trabalho"

    ```groovy title="main.nf"
    /*
    * Parâmetros do pipeline
    */
    params.input_file = "greetings.csv"

    /*
    * Use echo para imprimir 'Hello World!' na saída padrão
    */
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

    /*
    * Use um utilitário de substituição de texto para converter a saudação para maiúsculas
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

    workflow {

        // create a channel for inputs from a CSV file
        greeting_ch = channel.fromPath(params.input_file).splitCsv().flatten()

        // emit a greeting
        sayHello(greeting_ch)

        // convert the greeting to uppercase
        convertToUpper(sayHello.out)
    }
    ```

#### Execute o fluxo de trabalho

Vamos executar o fluxo de trabalho para garantir que está funcionando como esperado.

```bash
nextflow run main.nf
```

```console title="Result of running the workflow"
 N E X T F L O W   ~  version 24.10.2

Launching `main.nf` [soggy_linnaeus] DSL2 - revision: bbf79d5c31

executor >  local (6)
[f7/c3be66] sayHello (3)       | 3 of 3 ✔
[cd/e15303] convertToUpper (3) | 3 of 3 ✔
```

PARABÉNS! Você acabou de executar um teste!

"Espere, o quê? Eu apenas executei o fluxo de trabalho e funcionou! Como isso é um teste?"

Boa pergunta!

Vamos analisar o que acabou de acontecer.

Você executou o fluxo de trabalho com os parâmetros padrão, confirmou que funcionou e está satisfeito com os resultados. Esta é a essência dos testes. Se você trabalhou no curso de treinamento Hello Nextflow, você terá notado que sempre começamos cada seção executando o fluxo de trabalho que estávamos usando como ponto de partida, para confirmar que tudo está configurado corretamente.

Testar software essencialmente faz esse processo por nós.

#### Revise a tarefa

Seu desafio é adicionar testes padronizados a este fluxo de trabalho usando nf-test, a fim de facilitar a verificação de que cada parte continua funcionando como esperado caso quaisquer alterações adicionais sejam feitas.

<!-- TODO: give a bit more details, similar to how it's done in the Metadata side quest -->

#### Lista de verificação de prontidão

Acha que está pronto para mergulhar?

- [ ] Eu entendo o objetivo deste curso e seus pré-requisitos
- [ ] Meu codespace está funcionando
- [ ] Defini meu diretório de trabalho apropriadamente
- [ ] Executei o fluxo de trabalho com sucesso
- [ ] Eu entendo a tarefa

Se você pode marcar todas as caixas, está pronto para começar.

---

## 1. Inicialize `nf-test`

O pacote `nf-test` fornece um comando de inicialização que configura algumas coisas para que possamos começar a desenvolver testes para nosso projeto.

```bash
nf-test init
```

Isso deve produzir a seguinte saída:

```bash
🚀 nf-test 0.9.3
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr

Project configured. Configuration is stored in nf-test.config
```

Ele também cria um diretório `tests` contendo um esboço de arquivo de configuração.

### 1.1. Gere um nf-test

`nf-test` vem com um conjunto de ferramentas para construir arquivos nf-test, economizando a maior parte do trabalho. Elas vêm sob o subcomando `generate`. Vamos gerar um teste para o pipeline:

```bash
nf-test generate pipeline main.nf
```

```console title="Output"
> nf-test generate pipeline main.nf

Load source file '/workspaces/training/side-quests/nf-test/main.nf'
Wrote pipeline test file '/workspaces/training/side-quests/nf-test/tests/main.nf.test

SUCCESS: Generated 1 test files.
```

Isso criará um arquivo `main.nf.test` dentro do diretório `tests`. Este é nosso arquivo de teste em nível de pipeline. Se você executar `tree tests/` você deve ver algo assim:

```console title="Test directory contents"
tests/
├── main.nf.test
└── nextflow.config
```

O arquivo `main.nf.test` é nosso arquivo de teste em nível de pipeline. Vamos abri-lo e dar uma olhada no conteúdo.

```groovy title="tests/main.nf.test"
nextflow_pipeline {

    name "Test Workflow main.nf"
    script "main.nf"

    test("Should run without failures") {

        when {
            params {
                // define parameters here. Example:
                // outdir = "tests/results"
            }
        }

        then {
            assert workflow.success
        }

    }

}
```

Vamos dedicar um segundo para entender a estrutura do arquivo de teste.

O bloco `nextflow_pipeline` é o ponto de entrada para todos os testes em nível de pipeline. Ele contém o seguinte:

- `name`: O nome do teste.
- `script`: O caminho para o script do pipeline.

O bloco `test` é o teste real. Ele contém o seguinte:

- `when`: As condições sob as quais o teste deve ser executado. Isso inclui os parâmetros que serão usados para executar o pipeline.
- `then`: As asserções que devem ser feitas. Isso inclui os resultados esperados do pipeline.

Em português simples, a lógica do teste se lê da seguinte forma:
"**Quando** estes _parâmetros_ são fornecidos a este _pipeline_, **então** esperamos ver estes resultados."

Este não é um teste funcional, demonstraremos como transformá-lo em um na próxima seção.

### Uma Nota sobre Nomes de Testes

No exemplo acima, usamos o nome padrão "Should run without failures" que é apropriado para um teste básico que apenas verifica se o pipeline é executado com sucesso. No entanto, à medida que adicionamos casos de teste mais específicos, devemos usar nomes mais descritivos que indiquem o que estamos realmente testando. Por exemplo:

- "Should convert input to uppercase" - ao testar funcionalidade específica
- "Should handle empty input gracefully" - ao testar casos extremos
- "Should respect max memory parameter" - ao testar restrições de recursos
- "Should create expected output files" - ao testar geração de arquivos

Bons nomes de teste devem:

1. Começar com "Should" para deixar claro qual é o comportamento esperado
2. Descrever a funcionalidade ou cenário específico sendo testado
3. Ser claro o suficiente para que, se o teste falhar, você saiba qual funcionalidade está quebrada

À medida que adicionarmos mais asserções e casos de teste específicos mais tarde, usaremos esses nomes mais descritivos para deixar claro o que cada teste está verificando.

### 1.2. Execute o teste

Vamos executar o teste para ver o que acontece.

```bash
nf-test test tests/main.nf.test
```

```console title="nf-test pipeline fail"
> nf-test test tests/main.nf.test

🚀 nf-test 0.9.3
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Workflow main.nf

  Test [693ba951] 'Should run without failures' FAILED (4.652s)

  Assertion failed:

  assert workflow.success
         |        |
         workflow false

  Nextflow stdout:

  ERROR ~ No such file or directory: /workspaces/training/side-quests/nf-test/.nf-test/tests/693ba951a20fec36a5a9292ed1cc8a9f/greetings.csv

   -- Check '/workspaces/training/side-quests/nf-test/.nf-test/tests/693ba951a20fec36a5a9292ed1cc8a9f/meta/nextflow.log' file for details
  Nextflow stderr:

FAILURE: Executed 1 tests in 4.679s (1 failed)
```

O teste falha! O que aconteceu?

1. nf-test tentou executar o pipeline como está, usando as configurações no bloco `when`:

```groovy title="tests/main.nf.test"
when {
    params {
        // define parameters here. Example:
        // outdir = "tests/results"
    }
}
```

2. nf-test verificou o status do pipeline e o comparou com o bloco `when`:

```groovy title="tests/main.nf.test"
then {
    assert workflow.success
}
```

Note como nf-test relatou que o pipeline falhou e forneceu a mensagem de erro do Nextflow:

```console title="Error"
ERROR ~ No such file or directory: /workspaces/training/side-quests/nf-test/.nf-test/tests/693ba951a20fec36a5a9292ed1cc8a9f/greetings.csv
```

Então qual foi o problema? Lembre-se de que o pipeline tem um arquivo greetings.csv no diretório do projeto. Quando nf-test executa o pipeline, ele procurará este arquivo, mas não consegue encontrá-lo. O arquivo está lá, o que está acontecendo? Bem, se olharmos para o caminho, podemos ver que o teste está ocorrendo no caminho `./nf-test/tests/longHashString/`. Assim como o Nextflow, nf-test cria um novo diretório para cada teste para manter tudo isolado. O arquivo de dados não está localizado lá, então devemos corrigir o caminho para o arquivo no teste original.

Vamos voltar ao arquivo de teste e alterar o caminho para o arquivo no bloco `when`.

Você pode estar se perguntando como vamos apontar para a raiz do pipeline no teste. Como esta é uma situação comum, nf-test tem uma variedade de variáveis globais que podemos usar para facilitar nossas vidas. Você pode encontrar a lista completa [aqui](https://www.nf-test.com/docs/testcases/global_variables/) mas, enquanto isso, usaremos a variável `projectDir`, que significa a raiz do projeto do pipeline.

_Antes:_

```groovy title="tests/main.nf.test" linenums="1" hl_lines="3 4"
when {
    params {
        // define parameters here. Example:
        // outdir = "tests/results"
    }
}
```

_Depois:_

```groovy title="tests/main.nf.test" linenums="1" hl_lines="3"
when {
    params {
        input_file = "${projectDir}/greetings.csv"
    }
}
```

Vamos executar o teste novamente para ver se funciona.

```bash title="nf-test pipeline pass"
nf-test test tests/main.nf.test
```

```console title="Pipeline passes"
> nf-test test tests/main.nf.test

🚀 nf-test 0.9.3
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Workflow main.nf

  Test [1d4aaf12] 'Should run without failures' PASSED (1.619s)


SUCCESS: Executed 1 tests in 1.626s
```

Sucesso! O pipeline é executado com sucesso e o teste passa. Execute quantas vezes quiser e você sempre obterá o mesmo resultado!

Por padrão, a saída do Nextflow está oculta, mas para se convencer de que nf-test está definitivamente executando o fluxo de trabalho, você pode usar a flag `--verbose`:

```bash
nf-test test tests/main.nf.test --verbose
```

```console title="Pipeline runs all processes"
> nf-test test tests/main.nf.test

🚀 nf-test 0.9.3
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Workflow main.nf

  Test [693ba951] 'Should run without failures'
    > Nextflow 24.10.4 is available - Please consider updating your version to it
    > N E X T F L O W  ~  version 24.10.0
    > Launching `/workspaces/training/side-quests/nf-test/main.nf` [zen_ampere] DSL2 - revision: bbf79d5c31
    > [2b/61e453] Submitted process > sayHello (2)
    > [31/4e1606] Submitted process > sayHello (1)
    > [bb/5209ee] Submitted process > sayHello (3)
    > [83/83db6f] Submitted process > convertToUpper (2)
    > [9b/3428b1] Submitted process > convertToUpper (1)
    > [ca/0ba51b] Submitted process > convertToUpper (3)
    PASSED (5.206s)


SUCCESS: Executed 1 tests in 5.239s
```

### 1.3. Adicione asserções

Uma verificação simples é garantir que nosso pipeline esteja executando todos os processos que esperamos e não pulando nenhum silenciosamente. Lembre-se de que nosso pipeline executa 6 processos, um chamado `sayHello` e um chamado `convertToUpper` para cada uma das 3 saudações.

Vamos adicionar uma asserção ao nosso teste para verificar se o pipeline executa o número esperado de processos. Também atualizaremos o nome do nosso teste para refletir melhor o que estamos testando.

**Antes:**

```groovy title="tests/main.nf.test" linenums="1" hl_lines="1"
    test("Should run without failures") {

        when {
            params {
                input_file = "${projectDir}/greetings.csv"
            }
        }

        then {
            assert workflow.success
        }

    }
```

**Depois:**

```groovy title="tests/main.nf.test" linenums="1" hl_lines="1 11"
    test("Should run successfully with correct number of processes") {

        when {
            params {
                input_file = "${projectDir}/greetings.csv"
            }
        }

        then {
            assert workflow.success
            assert workflow.trace.tasks().size() == 6
        }

    }
```

O nome do teste agora reflete melhor o que estamos realmente verificando - não apenas que o pipeline é executado sem falhar, mas que ele executa o número esperado de processos.

Vamos executar o teste novamente para ver se funciona.

```bash title="nf-test pipeline pass"
nf-test test tests/main.nf.test
```

```console title="Pipeline passes with assertions"
🚀 nf-test 0.9.3
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Workflow main.nf

  Test [1d4aaf12] 'Should run successfully with correct number of processes' PASSED (1.567s)


SUCCESS: Executed 1 tests in 1.588s
```

Sucesso! O pipeline é executado com sucesso e o teste passa. Agora começamos a testar os detalhes do pipeline, bem como o status geral.

### 1.4. Teste a saída

Vamos adicionar uma asserção ao nosso teste para verificar se o arquivo de saída foi criado. Vamos adicioná-lo como um teste separado, com um nome informativo, para tornar os resultados mais fáceis de interpretar.

**Antes:**

```groovy title="tests/main.nf.test" linenums="1" hl_lines="14"
    test("Should run successfully with correct number of processes") {

        when {
            params {
                input_file = "${projectDir}/greetings.csv"
            }
        }

        then {
            assert workflow.success
            assert workflow.trace.tasks().size() == 6
        }

    }
```

**Depois:**

```groovy title="tests/main.nf.test" linenums="1" hl_lines="14-33"
    test("Should run successfully with correct number of processes") {

        when {
            params {
                input_file = "${projectDir}/greetings.csv"
            }
        }

        then {
            assert workflow.success
            assert workflow.trace.tasks().size() == 6
        }

    }

    test("Should produce correct output files") {

        when {
            params {
                input_file = "${projectDir}/greetings.csv"
            }
        }

        then {
            assert file("$launchDir/results/Bonjour-output.txt").exists()
            assert file("$launchDir/results/Hello-output.txt").exists()
            assert file("$launchDir/results/Holà-output.txt").exists()
            assert file("$launchDir/results/UPPER-Bonjour-output.txt").exists()
            assert file("$launchDir/results/UPPER-Hello-output.txt").exists()
            assert file("$launchDir/results/UPPER-Holà-output.txt").exists()
        }

    }
```

Execute o teste novamente para ver se funciona.

```bash title="nf-test pipeline pass"
nf-test test tests/main.nf.test
```

```console title="Pipeline passes with file assertions"
> nf-test test tests/main.nf.test

🚀 nf-test 0.9.3
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Workflow main.nf

  Test [f0e08a68] 'Should run successfully with correct number of processes' PASSED (8.144s)
  Test [d7e32a32] 'Should produce correct output files' PASSED (6.994s)


SUCCESS: Executed 2 tests in 15.165s
```

Sucesso! Os testes passam porque o pipeline foi concluído com sucesso, o número correto de processos foi executado e os arquivos de saída foram criados. Isso também deve mostrar a você como é útil fornecer esses nomes informativos para seus testes.

Isso é apenas a superfície, podemos continuar escrevendo asserções para verificar os detalhes do pipeline, mas por enquanto vamos passar para testar os componentes internos do pipeline.

### Conclusão

Você sabe como escrever um nf-test para um pipeline.

### O que vem a seguir?

Aprenda como testar um processo Nextflow.

---

## 2. Teste um processo Nextflow

Não precisamos escrever testes para cada parte do pipeline, mas quanto mais testes tivermos, mais abrangentes podemos ser sobre o pipeline e mais confiantes podemos estar de que está funcionando como esperado. Nesta seção, vamos testar ambos os processos no pipeline como unidades individuais.

### 2.1. Teste o processo `sayHello`

Vamos começar com o processo `sayHello`.

Vamos usar o comando `nf-test generate` novamente para gerar testes para o processo.

```bash
nf-test generate process main.nf
```

```console title="Output"
> nf-test generate process main.nf

Load source file '/workspaces/training/side-quests/nf-test/main.nf'
Wrote process test file '/workspaces/training/side-quests/nf-test/tests/main.sayhello.nf.test
Wrote process test file '/workspaces/training/side-quests/nf-test/tests/main.converttoupper.nf.test

SUCCESS: Generated 2 test files.
```

Vamos focar por enquanto no processo `sayhello` no arquivo `main.sayhello.nf.test`.

Vamos abrir o arquivo e dar uma olhada no conteúdo.

```groovy title="tests/main.sayhello.nf.test"
nextflow_process {

    name "Test Process sayHello"
    script "main.nf"
    process "sayHello"

    test("Should run without failures") {

        when {
            params {
                // define parameters here. Example:
                // outdir = "tests/results"
            }
            process {
                """
                // define inputs of the process here. Example:
                // input[0] = file("test-file.txt")
                """
            }
        }

        then {
            assert process.success
            assert snapshot(process.out).match()
        }

    }

}
```

Como antes, começamos com os detalhes do teste, seguidos pelos blocos `when` e `then`. No entanto, também temos um bloco `process` adicional que nos permite definir as entradas para o processo.

Vamos executar o teste para ver se funciona.

```bash title="nf-test pipeline pass"
nf-test test tests/main.sayhello.nf.test
```

```console title="Process test fails"
> nf-test test tests/main.sayhello.nf.test

🚀 nf-test 0.9.3
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Process sayHello

  Test [1eaad118] 'Should run without failures' FAILED (4.876s)

  Assertion failed:

  assert process.success
         |       |
         |       false
         sayHello

  Nextflow stdout:

  Process `sayHello` declares 1 input but was called with 0 arguments
  Nextflow stderr:

FAILURE: Executed 1 tests in 4.884s (1 failed)
```

O teste falha porque o processo `sayHello` declara 1 entrada mas foi chamado com 0 argumentos. Vamos corrigir isso adicionando uma entrada ao processo. Lembre-se de [Hello Workflow](../hello_nextflow/03_hello_workflow.md) (e da seção de aquecimento acima) que nosso processo `sayHello` recebe uma única entrada de valor, que precisaremos fornecer. Também devemos corrigir o nome do teste para refletir melhor o que estamos testando.

**Antes:**

```groovy title="tests/main.sayhello.nf.test" linenums="1" hl_lines="1 10 11"
    test("Should run without failures") {

        when {
            params {
                // define parameters here. Example:
                // outdir = "tests/results"
            }
            process {
                """
                // define inputs of the process here. Example:
                // input[0] = file("test-file.txt")
                """
            }
        }

        then {
            assert process.success
            assert snapshot(process.out).match()
        }

    }
```

**Depois:**

```groovy title="tests/main.sayhello.nf.test" linenums="1" hl_lines="1 10"
    test("Should run without failures and produce correct output") {

        when {
            params {
                // define parameters here. Example:
                // outdir = "tests/results"
            }
            process {
                """
                input[0] = "hello"
                """
            }
        }

        then {
            assert process.success
            assert snapshot(process.out).match()
        }

    }
```

Vamos executar o teste novamente para ver se funciona.

```console title="nf-test pipeline pass"
> nf-test test tests/main.sayhello.nf.test

🚀 nf-test 0.9.3
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Process sayHello

  Test [f91a1bcd] 'Should run without failures and produce correct output' PASSED (1.604s)
  Snapshots:
    1 created [Should run without failures and produce correct output]


Snapshot Summary:
  1 created

SUCCESS: Executed 1 tests in 1.611s
```

Sucesso! O teste passa porque o processo `sayHello` foi executado com sucesso e a saída foi criada.

### 2.2. Confira o snapshot criado pelo teste

Se olharmos para o arquivo `tests/main.sayhello.nf.test`, podemos ver que ele usa um método `snapshot()` no bloco de asserção:

```groovy title="tests/main.sayhello.nf.test"
assert snapshot(process.out).match()
```

Isso está dizendo ao nf-test para criar um snapshot da saída do processo `sayHello`. Vamos dar uma olhada no conteúdo do arquivo de snapshot.

```console title="Snapshot file contents"
code tests/main.sayhello.nf.test.snap
```

Não vamos imprimi-lo aqui, mas você deve ver um arquivo JSON contendo detalhes do processo e saídas do processo. Em particular, podemos ver uma linha que se parece com isso:

```json title="Snapshot file contents"
"0": [
    "hello-output.txt:md5,b1946ac92492d2347c6235b4d2611184"
]
```

Isso representa as saídas criadas pelo processo `sayHello`, que estamos testando explicitamente. Se executarmos novamente o teste, o programa verificará se a nova saída corresponde à saída que foi originalmente registrada. Esta é uma maneira rápida e simples de testar que as saídas do processo não mudam, e é por isso que nf-test fornece isso como padrão.

!!!warning

    Isso significa que temos que ter certeza de que a saída que registramos na execução original está correta!

Se, no curso do desenvolvimento futuro, algo no código mudar que faça com que a saída seja diferente, o teste falhará e teremos que determinar se a mudança é esperada ou não.

- Se acontecer de algo no código ter quebrado, teremos que consertá-lo, com a expectativa de que o código corrigido passará no teste.
- Se for uma mudança esperada (por exemplo, a ferramenta foi melhorada e os resultados são melhores), então precisaremos atualizar o snapshot para aceitar a nova saída como a referência a corresponder. nf-test tem um parâmetro `--update-snapshot` para este propósito.

Podemos executar o teste novamente e ver que o teste deve passar:

```console title="nf-test process pass with snapshot"
> nf-test test tests/main.sayhello.nf.test

🚀 nf-test 0.9.3
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Process sayHello

  Test [f91a1bcd] 'Should run without failures and produce correct output' PASSED (1.675s)


SUCCESS: Executed 1 tests in 1.685s
```

Sucesso! O teste passa porque o processo `sayHello` foi executado com sucesso e a saída correspondeu ao snapshot.

### 2.3. Alternativa aos Snapshots: Asserções Diretas de Conteúdo

Embora os snapshots sejam ótimos para capturar quaisquer mudanças na saída, às vezes você quer verificar conteúdo específico sem ser tão rigoroso sobre o arquivo inteiro corresponder. Por exemplo:

- Quando partes da saída podem mudar (timestamps, IDs aleatórios, etc.) mas certo conteúdo chave deve estar presente
- Quando você quer verificar padrões ou valores específicos na saída
- Quando você quer tornar o teste mais explícito sobre o que constitui sucesso

Aqui está como poderíamos modificar nosso teste para verificar conteúdo específico:

**Antes:**

```groovy title="tests/main.sayhello.nf.test" linenums="1" hl_lines="1 5 6 17"
    test("Should run without failures and produce correct output") {

        when {
            params {
                // define parameters here. Example:
                // outdir = "tests/results"
            }
            process {
                """
                input[0] = "hello"
                """
            }
        }

        then {
            assert process.success
            assert snapshot(process.out).match()
        }

    }
```

**Depois:**

```groovy title="tests/main.sayhello.nf.test" linenums="1" hl_lines="1 5 16 17"
    test("Should run without failures and contain expected greeting") {

        when {
            params {
                // define parameters here
            }
            process {
                """
                input[0] = "hello"
                """
            }
        }

        then {
            assert process.success
            assert path(process.out[0][0]).readLines().contains('hello')
            assert !path(process.out[0][0]).readLines().contains('HELLO')
        }

    }
```

Note que nf-test vê as saídas do processo como uma lista de listas, então `process.out[0][0]` está buscando a primeira parte do primeiro item do canal (ou 'emissão') deste processo.

Esta abordagem:

- Deixa claro exatamente o que esperamos na saída
- É mais resiliente a mudanças irrelevantes na saída
- Fornece melhores mensagens de erro quando os testes falham
- Permite validações mais complexas (padrões regex, comparações numéricas, etc.)

Vamos executar o teste para ver se funciona.

```bash title="nf-test pipeline pass"
nf-test test tests/main.sayhello.nf.test
```

```console title="Process test fails"
> nf-test test tests/main.sayhello.nf.test

🚀 nf-test 0.9.3
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Process sayHello

  Test [58df4e4b] 'Should run without failures and contain expected greeting' PASSED (7.196s)


SUCCESS: Executed 1 tests in 7.208s
```

### 2.4. Teste o processo `convertToUpper`

Vamos abrir o arquivo `tests/main.converttoupper.nf.test` e dar uma olhada no conteúdo:

```groovy title="tests/main.converttoupper.nf.test"
nextflow_process {

    name "Test Process convertToUpper"
    script "main.nf"
    process "convertToUpper"

    test("Should run without failures") {

        when {
            params {
                // define parameters here. Example:
                // outdir = "tests/results"
            }
            process {
                """
                // define inputs of the process here. Example:
                // input[0] = file("test-file.txt")
                """
            }
        }

        then {
            assert process.success
            assert snapshot(process.out).match()
        }

    }

}
```

Este é um teste semelhante ao processo `sayHello`, mas está testando o processo `convertToUpper`. Sabemos que este falhará porque, assim como com `sayHello`, o processo `convertToUpper` recebe uma única entrada de caminho, mas não especificamos uma.

Agora precisamos fornecer um único arquivo de entrada para o processo convertToUpper, que inclui algum texto que queremos converter para maiúsculas. Existem muitas maneiras de fazer isso:

- Poderíamos criar um arquivo dedicado para testar
- Poderíamos reutilizar o arquivo data/greetings.csv existente
- Poderíamos criá-lo dinamicamente dentro do teste

Por enquanto, vamos reutilizar o arquivo data/greetings.csv existente usando o exemplo que usamos com o teste em nível de pipeline. Como antes, podemos nomear o teste para refletir melhor o que estamos testando, mas desta vez vamos deixá-lo fazer um 'snapshot' do conteúdo em vez de verificar strings específicas (como fizemos no outro processo).

**Antes:**

```groovy title="tests/main.converttoupper.nf.test" linenums="1" hl_lines="1 10 11"
    test("Should run without failures") {

        when {
            params {
                // define parameters here. Example:
                // outdir = "tests/results"
            }
            process {
                """
                // define inputs of the process here. Example:
                // input[0] = file("test-file.txt")
                """
            }
        }

        then {
            assert process.success
            assert snapshot(process.out).match()
        }

    }
```

**Depois:**

```groovy title="tests/main.converttoupper.nf.test" linenums="1" hl_lines="1 10"
    test("Should run without failures and produce correct output") {

        when {
            params {
                // define parameters here. Example:
                // outdir = "tests/results"
            }
            process {
                """
                input[0] = "${projectDir}/greetings.csv"
                """
            }
        }

        then {
            assert process.success
            assert snapshot(process.out).match()
        }

    }
```

E execute o teste!

```bash title="nf-test pipeline pass"
nf-test test tests/main.converttoupper.nf.test
```

```console title="nf-test process convertToUpper pass"
> nf-test test tests/main.converttoupper.nf.test

🚀 nf-test 0.9.3
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Process convertToUpper

  Test [c59b6044] 'Should run without failures and produce correct output' PASSED (1.755s)
  Snapshots:
    1 created [Should run without failures and produce correct output]


Snapshot Summary:
  1 created

SUCCESS: Executed 1 tests in 1.764s
```

Note, criamos um arquivo de snapshot para o processo `convertToUpper` em `tests/main.converttoupper.nf.test.snap`. Se executarmos o teste novamente, devemos ver o nf-test passar novamente.

```bash title="nf-test process convertToUpper pass"
nf-test test tests/main.converttoupper.nf.test
```

```console title="nf-test process convertToUpper pass"
> nf-test test tests/main.converttoupper.nf.test

🚀 nf-test 0.9.3
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Process convertToUpper

  Test [c59b6044] 'Should run without failures and produce correct output' PASSED (1.798s)


SUCCESS: Executed 1 tests in 1.811s
```

### Conclusão

Você sabe como escrever testes para um processo Nextflow e executá-los.

### O que vem a seguir?

Aprenda como executar testes para tudo de uma vez!

## 3. Execute testes para todo o repositório

Executar nf-test em cada componente é bom, mas trabalhoso e propenso a erros. Não podemos simplesmente testar tudo de uma vez?

Sim, podemos!

Vamos executar nf-test em todo o repositório.

### 3.1. Execute nf-test em todo o repositório

Podemos executar nf-test em todo o repositório executando o comando `nf-test test`.

```bash
nf-test test .
```

Note, estamos apenas usando o `.` para executar tudo do nosso diretório atual. Isso incluirá todos os testes!

```console title="nf-test repo pass"
> nf-test test .

🚀 nf-test 0.9.3
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Process convertToUpper

  Test [3d26d9af] 'Should run without failures and produce correct output' PASSED (4.155s)

Test Workflow main.nf

  Test [f183df37] 'Should run successfully with correct number of processes' PASSED (3.33s)
  Test [d7e32a32] 'Should produce correct output files' PASSED (3.102s)

Test Process sayHello

  Test [58df4e4b] 'Should run without failures and contain expected greeting' PASSED (2.614s)


SUCCESS: Executed 4 tests in 13.481s
```

Olha só! Executamos 4 testes, 1 para cada processo e 2 para todo o pipeline com um único comando. Imagine quão poderoso isso é em uma base de código grande!

---

## Resumo

Nesta missão secundária, você aprendeu a aproveitar os recursos do nf-test para criar e executar testes para processos individuais, bem como testes de ponta a ponta para todo o pipeline.
Agora você está ciente das duas principais abordagens para validação de saída, snapshots e asserções diretas de conteúdo, e quando usar cada uma.
Você também sabe como executar testes um por um ou para um projeto inteiro.

Aplicar essas técnicas em seu próprio trabalho permitirá que você garanta que:

- Seu código funciona como esperado
- Mudanças não quebram funcionalidades existentes
- Outros desenvolvedores podem contribuir com confiança
- Problemas podem ser identificados e corrigidos rapidamente
- O conteúdo da saída corresponde às expectativas

### Padrões principais

<!-- TODO: Can we add snippets of code below to illustrate? -->

1. Testes em nível de pipeline:
   - Teste básico de sucesso
   - Verificação de contagem de processos
   - Verificações de existência de arquivos de saída
2. Testes em nível de processo
3. Duas abordagens para validação de saída:
   - Usando snapshots para verificação completa de saída
   - Usando asserções diretas de conteúdo para verificações de conteúdo específico
4. Executando todos os testes em um repositório com um único comando

### Recursos adicionais

Confira a [documentação do nf-test](https://www.nf-test.com/) para recursos de teste mais avançados e melhores práticas. Você pode querer:

- Adicionar asserções mais abrangentes aos seus testes
- Escrever testes para casos extremos e condições de erro
- Configurar integração contínua para executar testes automaticamente
- Aprender sobre outros tipos de testes como testes de fluxo de trabalho e módulo
- Explorar técnicas mais avançadas de validação de conteúdo

**Lembre-se:** Testes são documentação viva de como seu código deve se comportar. Quanto mais testes você escrever, e quanto mais específicas forem suas asserções, mais confiante você pode estar na confiabilidade do seu pipeline.

---

## O que vem a seguir?

Retorne ao [menu de Missões Secundárias](./index.md) ou clique no botão no canto inferior direito da página para passar para o próximo tópico da lista.
