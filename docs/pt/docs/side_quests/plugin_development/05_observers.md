# Parte 5: Trace Observers

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Tradução assistida por IA - [saiba mais e sugira melhorias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Trace observers permitem que seu plugin responda a eventos do fluxo de trabalho, como a conclusão de uma tarefa, a publicação de um arquivo ou o término do pipeline.
Isso possibilita casos de uso como relatórios personalizados, notificações no Slack, coleta de métricas ou integração com sistemas externos de monitoramento.
Nesta seção, você vai construir um observer que conta as tarefas concluídas e imprime um resumo.

!!! tip "Começando por aqui?"

    Se você está entrando nesta parte agora, copie a solução da Parte 4 para usar como ponto de partida:

    ```bash
    cp -r solutions/4-build-and-test/* .
    ```

---

## 1. Entendendo o trace observer existente

A mensagem "Pipeline is starting!" exibida quando você executou o pipeline veio da classe `GreetingObserver` do seu plugin.

Veja o código do observer:

```bash
cat nf-greeting/src/main/groovy/training/plugin/GreetingObserver.groovy
```

```groovy title="Output" hl_lines="30 32-34 37-39"
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
import groovy.util.logging.Slf4j
import nextflow.Session
import nextflow.trace.TraceObserver

/**
 * Implementa um observer que permite adicionar lógica personalizada
 * aos eventos de execução do Nextflow.
 */
@Slf4j
@CompileStatic
class GreetingObserver implements TraceObserver {    // (1)!

    @Override
    void onFlowCreate(Session session) {            // (2)!
        println "Pipeline is starting! 🚀"
    }

    @Override
    void onFlowComplete() {                         // (3)!
        println "Pipeline complete! 👋"
    }
}
```

1. Interface para se conectar aos eventos do ciclo de vida do fluxo de trabalho
2. Chamado quando o fluxo de trabalho inicia; recebe a sessão para acessar a configuração
3. Chamado quando o fluxo de trabalho termina com sucesso

Há dois pontos importantes a observar:

1. **`class GreetingObserver implements TraceObserver`**: `TraceObserver` é uma interface definida pelo Nextflow. Se sua classe implementa essa interface, o Nextflow pode se conectar a ela e chamar seus métodos quando os eventos ocorrem.
2. **`@Override`**: A interface `TraceObserver` define métodos como `onFlowCreate` e `onFlowComplete`. Quando você escreve métodos com esses nomes e adiciona a anotação `@Override`, o Nextflow os chama no momento apropriado. Quaisquer métodos que você não sobrescrever são ignorados.

O conjunto completo de eventos do ciclo de vida aos quais você pode se conectar no momento da escrita são:

| Método              | Quando é chamado                    |
| ------------------- | ----------------------------------- |
| `onFlowCreate`      | O fluxo de trabalho inicia          |
| `onFlowComplete`    | O fluxo de trabalho termina         |
| `onProcessStart`    | Uma tarefa começa a ser executada   |
| `onProcessComplete` | Uma tarefa termina                  |
| `onProcessCached`   | Uma tarefa em cache é reutilizada   |
| `onFilePublish`     | Um arquivo é publicado              |

Para uma lista completa, consulte a [interface TraceObserver](https://github.com/nextflow-io/nextflow/blob/master/modules/nextflow/src/main/groovy/nextflow/trace/TraceObserver.groovy) no código-fonte do Nextflow.

---

## 2. Adicionar um observer contador de tarefas

O objetivo é construir um observer que conta as tarefas concluídas e imprime um resumo ao final.
Adicionar um novo observer a um plugin requer duas coisas: escrever a classe do observer e registrá-la na factory para que o Nextflow a carregue.

### 2.1. Criar um observer mínimo

Crie um novo arquivo:

```bash
touch nf-greeting/src/main/groovy/training/plugin/TaskCounterObserver.groovy
```

Comece com o observer mais simples possível, que imprime uma mensagem quando qualquer tarefa é concluída:

```groovy title="nf-greeting/src/main/groovy/training/plugin/TaskCounterObserver.groovy" linenums="1"
package training.plugin

import groovy.transform.CompileStatic
import nextflow.processor.TaskHandler       // (1)!
import nextflow.trace.TraceObserver
import nextflow.trace.TraceRecord

/**
 * Observer que responde à conclusão de tarefas
 */
@CompileStatic
class TaskCounterObserver implements TraceObserver {  // (2)!

    @Override
    void onProcessComplete(TaskHandler handler, TraceRecord trace) {  // (3)!
        println "✓ Task completed!"
    }
}
```

1. Importe as classes necessárias: `TraceObserver`, `TaskHandler` e `TraceRecord`
2. Crie uma classe que `implements TraceObserver`
3. Sobrescreva `onProcessComplete` para executar código quando uma tarefa terminar

Isso é o mínimo necessário:

- Importar as classes necessárias (`TraceObserver`, `TaskHandler`, `TraceRecord`)
- Criar uma classe que `implements TraceObserver`
- Sobrescrever `onProcessComplete` para fazer algo quando uma tarefa terminar

### 2.2. Registrar o observer

O `GreetingFactory` cria os observers.
Veja como ele está:

```bash
cat nf-greeting/src/main/groovy/training/plugin/GreetingFactory.groovy
```

```groovy title="Output" hl_lines="25 27-29"
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
import nextflow.trace.TraceObserver
import nextflow.trace.TraceObserverFactory

@CompileStatic
class GreetingFactory implements TraceObserverFactory {

    @Override
    Collection<TraceObserver> create(Session session) {
        return List.<TraceObserver>of(new GreetingObserver())
    }

}
```

Edite o `GreetingFactory.groovy` para adicionar o novo observer:

=== "Depois"

    ```groovy title="GreetingFactory.groovy" linenums="31" hl_lines="3-6"
    @Override
    Collection<TraceObserver> create(Session session) {
        return [
            new GreetingObserver(),
            new TaskCounterObserver()
        ]
    }
    ```

=== "Antes"

    ```groovy title="GreetingFactory.groovy" linenums="31" hl_lines="3"
    @Override
    Collection<TraceObserver> create(Session session) {
        return List.<TraceObserver>of(new GreetingObserver())
    }
    ```

!!! note "Sintaxe de lista em Groovy"

    Substituímos o estilo Java `List.<TraceObserver>of(...)` pelo literal de lista mais simples do Groovy `[...]`.
    Ambos retornam uma `Collection`, mas a sintaxe Groovy é mais legível ao adicionar múltiplos itens.

### 2.3. Compilar, instalar e testar

```bash
cd nf-greeting && make install && cd ..
nextflow run greet.nf -ansi-log false
```

!!! tip "Por que `-ansi-log false`?"

    Por padrão, o display de progresso ANSI do Nextflow sobrescreve as linhas anteriores para mostrar uma visão limpa e atualizada do progresso.
    Isso significa que você veria apenas a *contagem final* de tarefas, e não as mensagens intermediárias.

    Usar `-ansi-log false` desativa esse comportamento e exibe toda a saída sequencialmente, o que é essencial ao testar observers que imprimem mensagens durante a execução.

Você deve ver "✓ Task completed!" impresso cinco vezes (uma vez por tarefa), intercalado com a saída existente do pipeline:

```console title="Output (partial)"
...
[9b/df7630] Submitted process > SAY_HELLO (4)
Decorated: *** Hello ***
✓ Task completed!
✓ Task completed!
Decorated: *** Holà ***
✓ Task completed!
...
Pipeline complete! 👋
```

O observer está funcionando.
Cada vez que uma tarefa termina, o Nextflow chama `onProcessComplete`, e nossa implementação imprime uma mensagem.

??? exercise "Personalizar a mensagem"

    Tente alterar a mensagem em `onProcessComplete` para algo de sua escolha, recompile e execute novamente.
    Isso confirma que o ciclo completo de edição-compilação-execução funciona para observers.

### 2.4. Adicionar lógica de contagem

O observer mínimo prova que o hook funciona, mas não rastreia nada.

Uma classe pode conter variáveis (chamadas de campos ou variáveis de instância) que persistem durante o tempo de vida do objeto.
Isso significa que um observer pode acumular estado ao longo de múltiplos eventos durante uma execução do pipeline.

A próxima versão adiciona uma variável contadora (`taskCount`) que começa em zero.
Cada vez que uma tarefa é concluída, o contador aumenta em um.
Quando todo o fluxo de trabalho termina, o observer imprime o total final.

Atualize o `TaskCounterObserver.groovy` com as alterações destacadas:

```groovy title="nf-greeting/src/main/groovy/training/plugin/TaskCounterObserver.groovy" linenums="1" hl_lines="14 18-19 22-24"
package training.plugin

import groovy.transform.CompileStatic
import nextflow.processor.TaskHandler
import nextflow.trace.TraceObserver
import nextflow.trace.TraceRecord

/**
 * Observer que conta as tarefas concluídas
 */
@CompileStatic
class TaskCounterObserver implements TraceObserver {

    private int taskCount = 0                // (1)!

    @Override
    void onProcessComplete(TaskHandler handler, TraceRecord trace) {
        taskCount++                          // (2)!
        println "📊 Tasks completed so far: ${taskCount}"
    }

    @Override
    void onFlowComplete() {                  // (3)!
        println "📈 Final task count: ${taskCount}"
    }
}
```

1. `taskCount` é uma variável que pertence ao objeto observer. Ela mantém seu valor entre as chamadas de método, permitindo acumular uma contagem ao longo de toda a execução do fluxo de trabalho. `private` significa que apenas esta classe pode acessá-la.
2. `taskCount++` adiciona um ao contador. Esta linha é executada toda vez que uma tarefa é concluída, então a contagem cresce conforme o fluxo de trabalho avança.
3. `onFlowComplete` é um segundo hook do ciclo de vida. Ele é executado uma vez quando o fluxo de trabalho termina, tornando-o um bom lugar para imprimir um resumo.

Em resumo:

- `taskCount` persiste entre as chamadas de método, acumulando uma contagem ao longo de toda a execução
- `onProcessComplete` incrementa o contador e imprime o total acumulado cada vez que uma tarefa termina
- `onFlowComplete` é executado uma vez ao final, imprimindo a contagem final

Recompile e teste:

```bash
cd nf-greeting && make install && cd ..
nextflow run greet.nf -ansi-log false
```

??? example "Saída"

    ```console
    N E X T F L O W  ~  version 25.10.2
    Launching `greet.nf` [pensive_engelbart] DSL2 - revision: 85fefd90d0
    Pipeline is starting! 🚀
    Reversed: olleH
    Reversed: ruojnoB
    Reversed: àloH
    Reversed: oaiC
    Reversed: ollaH
    [be/bd8e72] Submitted process > SAY_HELLO (2)
    [5b/d24c2b] Submitted process > SAY_HELLO (1)
    [14/1f9dbe] Submitted process > SAY_HELLO (3)
    Decorated: *** Bonjour ***
    Decorated: *** Hello ***
    [85/a6b3ad] Submitted process > SAY_HELLO (4)
    📊 Tasks completed so far: 1
    📊 Tasks completed so far: 2
    Decorated: *** Holà ***
    📊 Tasks completed so far: 3
    Decorated: *** Ciao ***
    [3c/be6686] Submitted process > SAY_HELLO (5)
    📊 Tasks completed so far: 4
    Decorated: *** Hallo ***
    📊 Tasks completed so far: 5
    Pipeline complete! 👋
    📈 Final task count: 5
    ```

    As mensagens do contador são intercaladas com os envios de tarefas porque os observers são executados conforme as tarefas são concluídas.

---

## 3. Rastrear arquivos publicados

O observer também pode responder quando arquivos são publicados.
O método `onFilePublish` recebe os caminhos de destino e de origem, que você pode usar para registrar, validar ou processar as saídas publicadas.

### 3.1. Adicionar um diretório de publicação

Primeiro, atualize o `greet.nf` para que o processo `SAY_HELLO` publique seus arquivos de saída:

=== "Depois"

    ```groovy title="greet.nf" linenums="10" hl_lines="2"
    process SAY_HELLO {
        publishDir 'results'
        input:
            val greeting
        output:
            path 'greeting.txt'
        script:
        // Use nossa função personalizada do plugin para decorar a saudação
        def decorated = decorateGreeting(greeting)
        """
        echo '$decorated' > greeting.txt
        """
    }
    ```

=== "Antes"

    ```groovy title="greet.nf" linenums="10"
    process SAY_HELLO {
        input:
            val greeting
        output:
            path 'greeting.txt'
        script:
        // Use nossa função personalizada do plugin para decorar a saudação
        def decorated = decorateGreeting(greeting)
        """
        echo '$decorated' > greeting.txt
        """
    }
    ```

### 3.2. Adicionar o método onFilePublish

Adicione um método `onFilePublish` e o import necessário ao `TaskCounterObserver.groovy`:

```groovy title="nf-greeting/src/main/groovy/training/plugin/TaskCounterObserver.groovy" linenums="1" hl_lines="5 23-26"
package training.plugin

import groovy.transform.CompileStatic
import nextflow.processor.TaskHandler
import java.nio.file.Path
import nextflow.trace.TraceObserver
import nextflow.trace.TraceRecord

/**
 * Observer que conta as tarefas concluídas
 */
@CompileStatic
class TaskCounterObserver implements TraceObserver {

    private int taskCount = 0

    @Override
    void onProcessComplete(TaskHandler handler, TraceRecord trace) {
        taskCount++
        println "📊 Tasks completed so far: ${taskCount}"
    }

    @Override
    void onFilePublish(Path destination, Path source) {
        println "📁 Published: ${destination.fileName}"
    }

    @Override
    void onFlowComplete() {
        println "📈 Final task count: ${taskCount}"
    }
}
```

### 3.3. Compilar e testar

```bash
cd nf-greeting && make install && cd ..
nextflow run greet.nf -ansi-log false
```

Você deve ver mensagens "Published:" para cada arquivo de saída, junto com a saída do contador de tarefas:

```console title="Output (partial)"
...
📊 Tasks completed so far: 1
📁 Published: greeting.txt
📊 Tasks completed so far: 2
📁 Published: greeting.txt
...
📈 Final task count: 5
Pipeline complete! 👋
```

O método `onFilePublish` é acionado cada vez que o Nextflow publica um arquivo no diretório `results`.
Esse padrão é útil para construir logs de auditoria, acionar ações subsequentes ou validar saídas conforme são produzidas.

---

## Conclusão

Você aprendeu que:

- Trace observers se conectam a eventos do ciclo de vida do fluxo de trabalho como `onFlowCreate`, `onProcessComplete`, `onFilePublish` e `onFlowComplete`
- Crie observers implementando `TraceObserver` e registrando-os em uma Factory
- Observers podem conter variáveis de instância para acumular estado ao longo dos eventos
- Observers são úteis para registro personalizado, coleta de métricas, notificações e relatórios

---

## O que vem a seguir?

O contador de tarefas funciona, mas está sempre ativo.
Em um plugin real, os usuários devem poder habilitar ou desabilitar funcionalidades, ou ajustar o comportamento, a partir do `nextflow.config` sem precisar editar o código-fonte do plugin.
A próxima seção mostra como tornar seu observer configurável e como compartilhar seu plugin finalizado com outras pessoas.

[Continuar para a Parte 6 :material-arrow-right:](06_configuration.md){ .md-button .md-button--primary }
