# Parte 4: Testes

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tradução assistida por IA - [saiba mais e sugira melhorias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Plugins são softwares independentes nos quais os desenvolvedores de pipelines precisam confiar.
Testar cada funcionalidade de forma independente, fora de um pipeline, garante que o plugin funcione corretamente antes que alguém o integre a um fluxo de trabalho.
Nesta seção, você vai escrever e executar testes usando o framework de testes Spock.

!!! tip "Começando a partir daqui?"

    Se você está entrando nesta parte, copie a solução da Parte 3 para usar como ponto de partida:

    ```bash
    cp -r solutions/3-custom-functions/* .
    ```

    Em seguida, entre no diretório do plugin:

    ```bash
    cd nf-greeting
    ```

Certifique-se de que você está no diretório do plugin:

```bash
cd nf-greeting
```

---

## 1. Por que testar?

Um build bem-sucedido significa que o código compila, mas não verifica se ele funciona como esperado.
Testes unitários são pequenos trechos de código que verificam automaticamente se suas funções produzem a saída correta para uma determinada entrada.
Por exemplo, um teste pode verificar que `#!groovy reverseGreeting("Hello")` retorna `"olleH"`.

Os testes são valiosos porque:

- Eles detectam bugs antes que os usuários o façam
- Eles dão confiança para fazer alterações sem quebrar o que já funciona
- Eles servem como documentação mostrando como as funções devem ser usadas

---

## 2. Entendendo os testes Spock

O template do plugin usa o [Spock](https://spockframework.org/), um framework de testes para Groovy.
O Spock já está configurado no projeto (via `build.gradle`), então você não precisa adicionar nada.

Se você já usou ferramentas de teste antes (como `pytest` em Python ou `testthat` em R), o Spock cumpre o mesmo papel: você escreve pequenas funções que chamam seu código com entradas conhecidas e verificam as saídas.
A diferença é que o Spock usa blocos com rótulos (`given:`, `expect:`, `when:`, `then:`) que são semelhantes a um processo ou fluxo de trabalho do Nextflow.

Aqui está a estrutura básica:

```groovy
def 'should reverse a greeting'() {   // (1)!
    given:                             // (2)!
    def ext = new GreetingExtension()

    expect:                            // (3)!
    ext.reverseGreeting('Hello') == 'olleH'
}
```

1. **Nome do teste entre aspas**: Descreve o que o teste verifica. Use linguagem simples.
2. **Bloco `given:`**: Configure o que você precisa para o teste (criar objetos, preparar dados)
3. **Bloco `expect:`**: As verificações em si. Cada linha deve ser `true` para o teste passar

Essa estrutura torna os testes legíveis: "Dado um objeto de extensão, espera-se que `reverseGreeting('Hello')` seja igual a `'olleH'`."

---

## 3. Escreva os testes

Escreva testes para as duas funções que você criou na Parte 3: `reverseGreeting` e `decorateGreeting`.

### 3.1. Crie a classe de teste

```bash
touch src/test/groovy/training/plugin/GreetingExtensionTest.groovy
```

Abra o arquivo no seu editor e adicione o esqueleto vazio da classe de teste:

```groovy title="src/test/groovy/training/plugin/GreetingExtensionTest.groovy" linenums="1"
package training.plugin

import spock.lang.Specification

/**
 * Testes para as funções de extensão de saudação
 */
class GreetingExtensionTest extends Specification {  // (1)!

}
```

1. Todas as classes de teste Spock estendem `Specification`. Este é o ponto de partida para qualquer arquivo de teste Spock.

### 3.2. Teste reverseGreeting

Adicione um método de teste dentro do corpo da classe.
O bloco `given:` cria uma instância de `GreetingExtension`, e o bloco `expect:` verifica que `reverseGreeting` inverte corretamente duas entradas diferentes.
Isso testa a função diretamente, sem executar um pipeline.

=== "Depois"

    ```groovy title="GreetingExtensionTest.groovy" linenums="1" hl_lines="10-17"
    package training.plugin

    import spock.lang.Specification

    /**
     * Testes para as funções de extensão de saudação
     */
    class GreetingExtensionTest extends Specification {

        def 'should reverse a greeting'() {
            given:
            def ext = new GreetingExtension()            // (1)!

            expect:
            ext.reverseGreeting('Hello') == 'olleH'     // (2)!
            ext.reverseGreeting('Bonjour') == 'ruojnoB'
        }
    }
    ```

    1. Crie uma instância da sua extensão para testar diretamente, sem executar um pipeline
    2. Cada linha em `expect:` é uma asserção; o teste passa somente se todas forem `true`

=== "Antes"

    ```groovy title="GreetingExtensionTest.groovy" linenums="1"
    package training.plugin

    import spock.lang.Specification

    /**
     * Testes para as funções de extensão de saudação
     */
    class GreetingExtensionTest extends Specification {

    }
    ```

### 3.3. Teste decorateGreeting

Adicione um segundo método de teste após o primeiro.
Este verifica que `decorateGreeting` envolve a string de entrada com `***` em cada lado.

=== "Depois"

    ```groovy title="GreetingExtensionTest.groovy" linenums="1" hl_lines="18-25"
    package training.plugin

    import spock.lang.Specification

    /**
     * Testes para as funções de extensão de saudação
     */
    class GreetingExtensionTest extends Specification {

        def 'should reverse a greeting'() {
            given:
            def ext = new GreetingExtension()

            expect:
            ext.reverseGreeting('Hello') == 'olleH'
            ext.reverseGreeting('Bonjour') == 'ruojnoB'
        }

        def 'should decorate a greeting'() {
            given:
            def ext = new GreetingExtension()

            expect:
            ext.decorateGreeting('Hello') == '*** Hello ***'
        }
    }
    ```

=== "Antes"

    ```groovy title="GreetingExtensionTest.groovy" linenums="1" hl_lines="18"
    package training.plugin

    import spock.lang.Specification

    /**
     * Testes para as funções de extensão de saudação
     */
    class GreetingExtensionTest extends Specification {

        def 'should reverse a greeting'() {
            given:
            def ext = new GreetingExtension()

            expect:
            ext.reverseGreeting('Hello') == 'olleH'
            ext.reverseGreeting('Bonjour') == 'ruojnoB'
        }
    }
    ```

---

## 4. Execute os testes

```bash
make test
```

??? example "Saída dos testes"

    ```console
    BUILD SUCCESSFUL in 5s
    6 actionable tasks: 6 executed
    ```

    **Onde estão os resultados dos testes?** O Gradle oculta a saída detalhada quando todos os testes passam.
    "BUILD SUCCESSFUL" significa que tudo funcionou.
    Se algum teste falhar, você verá mensagens de erro detalhadas.

??? exercise "Adicione um teste de caso extremo"

    Adicione um teste que verifique se `reverseGreeting` lida com uma string vazia.
    O que `reverseGreeting('')` deve retornar?
    Adicione o teste, execute `make test` e verifique se ele passa.

    ??? solution "Solução"

        Adicione este método de teste ao `GreetingExtensionTest.groovy`:

        ```groovy
        def 'should handle empty string'() {
            given:
            def ext = new GreetingExtension()

            expect:
            ext.reverseGreeting('') == ''
        }
        ```

        Uma string vazia invertida continua sendo uma string vazia.

---

## 5. Visualize o relatório de testes

O Gradle gera um relatório de testes em HTML com resultados detalhados para cada teste.
Inicie um servidor web no diretório do relatório:

```bash
pushd build/reports/tests/test
python -m http.server
```

O VS Code vai solicitar que você abra a aplicação no seu navegador.
Clique até chegar à sua classe de teste para ver os resultados individuais de cada teste:

![Relatório de testes mostrando que todos os testes passaram](./img/test_report.png)

O relatório mostra cada método de teste e se ele passou ou falhou.

Pressione ++ctrl+c++ para parar o servidor e, em seguida, retorne ao diretório anterior:

```bash
popd
```

Volte ao diretório principal do projeto:

```bash
cd ..
```

---

## Conclusão

Você aprendeu que:

- Os testes Spock usam uma estrutura legível com `given:`/`expect:`
- Use `make test` para executar os testes e `build/reports/tests/test/` para o relatório HTML
- Os testes verificam o comportamento e servem como documentação de como as funções devem ser usadas

---

## O que vem a seguir?

Até agora, seu plugin adiciona funções personalizadas que os pipelines podem chamar.
Plugins também podem reagir a eventos do fluxo de trabalho (uma tarefa sendo concluída, um arquivo sendo publicado, o pipeline terminando) usando trace observers.
Na próxima seção, você vai construir um observer que conta as tarefas concluídas e imprime um resumo quando o pipeline termina.

[Continuar para a Parte 5 :material-arrow-right:](05_observers.md){ .md-button .md-button--primary }
