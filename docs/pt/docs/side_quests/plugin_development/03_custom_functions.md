# Parte 3: Funções Personalizadas

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tradução assistida por IA - [saiba mais e sugira melhorias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Ao final desta seção, você terá funções personalizadas no seu plugin, compiladas e instaladas localmente, sendo executadas em um fluxo de trabalho real.

!!! tip "Começando por aqui?"

    Se você está entrando nesta parte, copie a solução da Parte 2 para usar como ponto de partida:

    ```bash
    cp -r solutions/2-create-project/* .
    ```

---

## 1. Veja o que o template gerou

Antes de escrever suas próprias funções, observe a função de exemplo que o template criou para entender o padrão.

Entre no diretório do plugin:

```bash
cd nf-greeting
```

O template criou um arquivo chamado `GreetingExtension.groovy` onde as funções do plugin são definidas.
Abra-o para ver o ponto de partida:

```bash
cat src/main/groovy/training/plugin/GreetingExtension.groovy
```

```groovy title="Output" hl_lines="29 40-43"
/*
 * Copyright 2025, Seqera Labs
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

package training.plugin

import groovy.transform.CompileStatic
import nextflow.Session
import nextflow.plugin.extension.Function
import nextflow.plugin.extension.PluginExtensionPoint

/**
 * Implementa uma função personalizada que pode ser importada por
 * scripts Nextflow.
 */
@CompileStatic
class GreetingExtension extends PluginExtensionPoint { // (1)!

    @Override
    protected void init(Session session) {             // (2)!
    }

    /**
     * Diz olá para o alvo especificado.
     *
     * @param target
     */
    @Function                                          // (3)!
    void sayHello(String target) {
        println "Hello, ${target}!"
    }

}
```

1. A classe na qual sua extensão se baseia. O Nextflow exige isso para reconhecer suas funções.
2. Chamado quando o plugin é carregado; use para inicialização
3. Torna este método chamável a partir de fluxos de trabalho via `include`

O template inclui uma função `sayHello` de exemplo.
A anotação `@Function` é o que torna um método chamável a partir de fluxos de trabalho Nextflow.
Sem ela, o método existe apenas dentro do código do plugin.

Em Groovy (e Java), os métodos declaram o tipo que retornam e os tipos dos seus parâmetros.
Por exemplo, `String reverseGreeting(String greeting)` declara um método que recebe um parâmetro `String` e retorna um `String`.
A palavra-chave `void` significa que o método não retorna nada, como em `sayHello` acima.
Isso é diferente de Python ou R, onde os tipos não precisam ser declarados explicitamente.

---

## 2. Substitua sayHello por reverseGreeting

A função `sayHello` do template é um placeholder.
Substitua-a pela sua própria função para ver o ciclo completo de escrever, compilar e usar uma função de plugin.

Edite `src/main/groovy/training/plugin/GreetingExtension.groovy` para substituir o método `sayHello`:

=== "Depois"

    ```groovy title="GreetingExtension.groovy" linenums="24" hl_lines="8-14"
    @CompileStatic
    class GreetingExtension extends PluginExtensionPoint {

        @Override
        protected void init(Session session) {
        }

        /**
         * Inverte uma string de saudação
         */
        @Function                                // (1)!
        String reverseGreeting(String greeting) { // (2)!
            return greeting.reverse()             // (3)!
        }

    }
    ```

    1. Torna o método chamável a partir de fluxos de trabalho Nextflow
    2. Recebe um String, retorna um String
    3. Método de inversão de string nativo do Groovy

=== "Antes"

    ```groovy title="GreetingExtension.groovy" linenums="24" hl_lines="12-20"
    /**
     * Implementa uma função personalizada que pode ser importada por
     * scripts Nextflow.
     */
    @CompileStatic
    class GreetingExtension extends PluginExtensionPoint {

        @Override
        protected void init(Session session) {
        }

        /**
         * Diz olá para o alvo especificado.
         *
         * @param target
         */
        @Function
        void sayHello(String target) {
            println "Hello, ${target}!"
        }

    }
    ```

Partes principais desta função:

- **`@Function`**: Torna o método chamável a partir de fluxos de trabalho Nextflow
- **`String reverseGreeting(String greeting)`**: Recebe um String, retorna um String
- **`greeting.reverse()`**: Método de inversão de string nativo do Groovy

!!! tip "Métodos públicos e privados"

    Métodos sem `@Function` não são expostos aos fluxos de trabalho Nextflow.
    Você pode adicionar métodos auxiliares à sua classe sem se preocupar com eles vazando para o namespace do fluxo de trabalho.

---

## 3. Compile e instale seu plugin

Compile e instale o plugin:

```bash
make install
```

!!! tip "Se a compilação falhar"

    Leia a mensagem de erro com atenção; ela geralmente inclui um número de linha e descreve o problema.
    Causas comuns são erros de sintaxe (colchete ou aspas faltando), nomes de classes com erro de digitação e incompatibilidade de tipos.
    Se você estiver travado, compare seu código caractere por caractere com os exemplos.

---

## 4. Use sua função em um fluxo de trabalho

O plugin está compilado e instalado.
O próximo passo é usar `reverseGreeting` em um fluxo de trabalho para verificar se funciona de ponta a ponta.

Volte para o diretório do pipeline:

```bash
cd ..
```

Edite `greet.nf` para importar e usar `reverseGreeting`:

=== "Depois"

    ```groovy title="greet.nf" hl_lines="4 23-25" linenums="1"
    #!/usr/bin/env nextflow

    include { samplesheetToList } from 'plugin/nf-schema'
    include { reverseGreeting } from 'plugin/nf-greeting'

    params.input = 'greetings.csv'

    process SAY_HELLO {
        input:
            val greeting
        output:
            stdout
        script:
        """
        echo '$greeting'
        """
    }

    workflow {
        greeting_ch = Channel.fromList(samplesheetToList(params.input, 'greetings_schema.json'))
            .map { row -> row[0] }

        greeting_ch
            .map { greeting -> reverseGreeting(greeting) }
            .view { reversed -> "Reversed: $reversed" }

        SAY_HELLO(greeting_ch)
        SAY_HELLO.out.view { result -> "Output: ${result.trim()}" }
    }
    ```

=== "Antes"

    ```groovy title="greet.nf" linenums="1"
    #!/usr/bin/env nextflow

    include { samplesheetToList } from 'plugin/nf-schema'

    params.input = 'greetings.csv'

    process SAY_HELLO {
        input:
            val greeting
        output:
            stdout
        script:
        """
        echo '$greeting'
        """
    }

    workflow {
        greeting_ch = Channel.fromList(samplesheetToList(params.input, 'greetings_schema.json'))
            .map { row -> row[0] }
        SAY_HELLO(greeting_ch)
        SAY_HELLO.out.view { result -> "Output: ${result.trim()}" }
    }
    ```

Execute o pipeline:

```bash
nextflow run greet.nf
```

??? example "Output"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `greet.nf` [elated_marconi] DSL2 - revision: cd8d52c97c

    Pipeline is starting! 🚀
    executor >  local (5)
    [fe/109754] process > SAY_HELLO (5) [100%] 5 of 5 ✔
    Reversed: olleH
    Reversed: ruojnoB
    Reversed: àloH
    Reversed: oaiC
    Reversed: ollaH
    Output: Hello
    Output: Bonjour
    Output: Holà
    Output: Ciao
    Output: Hallo
    Pipeline complete! 👋
    ```

Sua primeira função de plugin personalizada está funcionando em um fluxo de trabalho real.
O mesmo padrão `include { ... } from 'plugin/...'` que você usou com nf-hello e nf-schema na Parte 1 funciona com o seu próprio plugin.

---

## 5. Adicione decorateGreeting

Um plugin pode fornecer múltiplas funções.
Adicione uma segunda que envolve uma saudação com marcadores decorativos; você a tornará configurável na Parte 6.

Edite `GreetingExtension.groovy` para adicionar `decorateGreeting` após `reverseGreeting`, antes da chave de fechamento da classe:

=== "Depois"

    ```groovy title="GreetingExtension.groovy" linenums="24" hl_lines="16-22"
    @CompileStatic
    class GreetingExtension extends PluginExtensionPoint {

        @Override
        protected void init(Session session) {
        }

        /**
         * Inverte uma string de saudação
         */
        @Function
        String reverseGreeting(String greeting) {
            return greeting.reverse()
        }

        /**
         * Decora uma saudação com marcadores comemorativos
         */
        @Function
        String decorateGreeting(String greeting) {
            return "*** ${greeting} ***"             // (1)!
        }

    }
    ```

    1. Interpolação de string do Groovy: `#!groovy ${...}` insere o valor da variável na string

=== "Antes"

    ```groovy title="GreetingExtension.groovy" linenums="24"
    @CompileStatic
    class GreetingExtension extends PluginExtensionPoint {

        @Override
        protected void init(Session session) {
        }

        /**
         * Inverte uma string de saudação
         */
        @Function
        String reverseGreeting(String greeting) {
            return greeting.reverse()
        }

    }
    ```

Esta função usa interpolação de string do Groovy (`"*** ${greeting} ***"`) para incorporar a variável de saudação dentro de uma string.

Compile, instale e atualize o fluxo de trabalho:

```bash
cd nf-greeting && make install && cd ..
```

Atualize `greet.nf` para também importar e usar `decorateGreeting`:

=== "Depois"

    ```groovy title="greet.nf" hl_lines="4-6 14 16-17 19 33" linenums="1"
    #!/usr/bin/env nextflow

    include { samplesheetToList } from 'plugin/nf-schema'
    // Importa funções personalizadas do nosso plugin
    include { reverseGreeting } from 'plugin/nf-greeting'
    include { decorateGreeting } from 'plugin/nf-greeting'  // (1)!

    params.input = 'greetings.csv'

    process SAY_HELLO {
        input:
            val greeting
        output:
            path 'greeting.txt'
        script:
        // Usa nossa função de plugin personalizada para decorar a saudação
        def decorated = decorateGreeting(greeting)  // (2)!
        """
        echo '$decorated' > greeting.txt
        """
    }

    workflow {
        greeting_ch = Channel.fromList(samplesheetToList(params.input, 'greetings_schema.json'))
            .map { row -> row[0] }

        // Demonstra o uso da função reverseGreeting
        greeting_ch
            .map { greeting -> reverseGreeting(greeting) }
            .view { reversed -> "Reversed: $reversed" }

        SAY_HELLO(greeting_ch)
        SAY_HELLO.out.view { file -> "Decorated: ${file.text.trim()}" }
    }
    ```

    1. Múltiplas funções do mesmo plugin precisam de declarações `include` separadas
    2. Funções de plugin também funcionam dentro de blocos `script:` de processos

=== "Antes"

    ```groovy title="greet.nf" linenums="1" hl_lines="4 12 15 28"
    #!/usr/bin/env nextflow

    include { samplesheetToList } from 'plugin/nf-schema'
    include { reverseGreeting } from 'plugin/nf-greeting'

    params.input = 'greetings.csv'

    process SAY_HELLO {
        input:
            val greeting
        output:
            stdout
        script:
        """
        echo '$greeting'
        """
    }

    workflow {
        greeting_ch = Channel.fromList(samplesheetToList(params.input, 'greetings_schema.json'))
            .map { row -> row[0] }

        greeting_ch
            .map { greeting -> reverseGreeting(greeting) }
            .view { reversed -> "Reversed: $reversed" }

        SAY_HELLO(greeting_ch)
        SAY_HELLO.out.view { result -> "Output: ${result.trim()}" }
    }
    ```

```bash
nextflow run greet.nf
```

??? example "Output"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `greet.nf` [elated_marconi] DSL2 - revision: cd8d52c97c

    Pipeline is starting! 🚀
    executor >  local (5)
    [fe/109754] process > SAY_HELLO (5) [100%] 5 of 5 ✔
    Reversed: olleH
    Reversed: ruojnoB
    Reversed: àloH
    Reversed: oaiC
    Reversed: ollaH
    Decorated: *** Hello ***
    Decorated: *** Bonjour ***
    Decorated: *** Holà ***
    Decorated: *** Ciao ***
    Decorated: *** Hallo ***
    Pipeline complete! 👋
    ```

As funções de plugin funcionam tanto em scripts de processos (como `decorateGreeting` dentro de `SAY_HELLO`) quanto em operações de fluxo de trabalho (como `reverseGreeting` em um `map`).

---

## Conclusão

Você aprendeu que:

- As funções são definidas com a anotação `@Function` em subclasses de `PluginExtensionPoint`
- Funções de plugin importadas com `include` funcionam de forma idêntica, seja do seu próprio plugin ou de um existente
- Funções de plugin funcionam tanto em scripts de processos quanto em operações de fluxo de trabalho

---

## O que vem a seguir?

Suas funções estão funcionando, mas até agora você só verificou isso executando o pipeline completo e conferindo a saída visualmente.
Essa abordagem não escala: à medida que você adiciona mais funções, precisa de uma maneira mais rápida de verificar que cada uma se comporta corretamente, especialmente após fazer alterações.
A próxima seção apresenta testes unitários, que permitem verificar funções individuais automaticamente sem executar um pipeline.

[Continuar para a Parte 4 :material-arrow-right:](04_build_and_test.md){ .md-button .md-button--primary }
