# Parte 3: Funzioni Personalizzate

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduzione assistita da IA - [scopri di più e suggerisci miglioramenti](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Alla fine di questa sezione, avrete funzioni personalizzate nel vostro plugin, compilate e installate localmente, in esecuzione in un flusso di lavoro reale.

!!! tip "Iniziate da qui?"

    Se vi unite a partire da questa parte, copiate la soluzione dalla Parte 2 da usare come punto di partenza:

    ```bash
    cp -r solutions/2-create-project/* .
    ```

---

## 1. Vedere cosa ha generato il template

Prima di scrivere le vostre funzioni, osservate la funzione di esempio che il template ha creato per capire il pattern.

Spostatevi nella directory del plugin:

```bash
cd nf-greeting
```

Il template ha creato un file chiamato `GreetingExtension.groovy` dove vengono definite le funzioni del plugin.
Apritelo per vedere il punto di partenza:

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
 * Implementa una funzione personalizzata che può essere importata da
 * script Nextflow.
 */
@CompileStatic
class GreetingExtension extends PluginExtensionPoint { // (1)!

    @Override
    protected void init(Session session) {             // (2)!
    }

    /**
     * Saluta il destinatario specificato.
     *
     * @param target
     */
    @Function                                          // (3)!
    void sayHello(String target) {
        println "Hello, ${target}!"
    }

}
```

1. La classe su cui si basa la vostra estensione. Nextflow richiede questa classe per riconoscere le vostre funzioni.
2. Chiamata quando il plugin viene caricato; da usare per l'inizializzazione
3. Rende questo metodo richiamabile dai flussi di lavoro tramite `include`

Il template include una funzione `sayHello` di esempio.
L'annotazione `@Function` è ciò che rende un metodo richiamabile dai flussi di lavoro Nextflow.
Senza di essa, il metodo esiste solo all'interno del codice del plugin.

In Groovy (e Java), i metodi dichiarano il tipo che restituiscono e i tipi dei loro parametri.
Ad esempio, `String reverseGreeting(String greeting)` dichiara un metodo che accetta un parametro `String` e restituisce una `String`.
La parola chiave `void` significa che il metodo non restituisce nulla, come nel caso di `sayHello` sopra.
Questo è diverso da Python o R, dove i tipi non devono essere dichiarati esplicitamente.

---

## 2. Sostituire sayHello con reverseGreeting

La funzione `sayHello` del template è un segnaposto.
Sostituitela con la vostra funzione per vedere il ciclo completo di scrittura, compilazione e utilizzo di una funzione del plugin.

Modificate `src/main/groovy/training/plugin/GreetingExtension.groovy` per sostituire il metodo `sayHello`:

=== "Dopo"

    ```groovy title="GreetingExtension.groovy" linenums="24" hl_lines="8-14"
    @CompileStatic
    class GreetingExtension extends PluginExtensionPoint {

        @Override
        protected void init(Session session) {
        }

        /**
         * Inverte una stringa di saluto
         */
        @Function                                // (1)!
        String reverseGreeting(String greeting) { // (2)!
            return greeting.reverse()             // (3)!
        }

    }
    ```

    1. Rende il metodo richiamabile dai flussi di lavoro Nextflow
    2. Accetta una String, restituisce una String
    3. Il metodo di inversione delle stringhe integrato in Groovy

=== "Prima"

    ```groovy title="GreetingExtension.groovy" linenums="24" hl_lines="12-20"
    /**
     * Implementa una funzione personalizzata che può essere importata da
     * script Nextflow.
     */
    @CompileStatic
    class GreetingExtension extends PluginExtensionPoint {

        @Override
        protected void init(Session session) {
        }

        /**
         * Saluta il destinatario specificato.
         *
         * @param target
         */
        @Function
        void sayHello(String target) {
            println "Hello, ${target}!"
        }

    }
    ```

Elementi chiave di questa funzione:

- **`@Function`**: Rende il metodo richiamabile dai flussi di lavoro Nextflow
- **`String reverseGreeting(String greeting)`**: Accetta una String, restituisce una String
- **`greeting.reverse()`**: Il metodo di inversione delle stringhe integrato in Groovy

!!! tip "Metodi pubblici e privati"

    I metodi senza `@Function` non sono esposti ai flussi di lavoro Nextflow.
    Potete aggiungere metodi di supporto alla vostra classe senza preoccuparvi che vengano esposti nel namespace del flusso di lavoro.

---

## 3. Compilare e installare il plugin

Compilate e installate il plugin:

```bash
make install
```

!!! tip "Se la compilazione fallisce"

    Leggete attentamente il messaggio di errore; di solito include un numero di riga e descrive il problema.
    Le cause più comuni sono errori di sintassi (parentesi o virgolette mancanti), nomi di classe scritti in modo errato e mancata corrispondenza dei tipi.
    Se siete bloccati, confrontate il vostro codice carattere per carattere con gli esempi.

---

## 4. Usare la funzione in un flusso di lavoro

Il plugin è compilato e installato.
Il passo successivo è usare `reverseGreeting` in un flusso di lavoro per verificare che funzioni dall'inizio alla fine.

Tornate alla directory della pipeline:

```bash
cd ..
```

Modificate `greet.nf` per importare e usare `reverseGreeting`:

=== "Dopo"

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

=== "Prima"

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

Eseguite la pipeline:

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

La vostra prima funzione personalizzata del plugin funziona in un flusso di lavoro reale.
Lo stesso pattern `include { ... } from 'plugin/...'` che avete usato con nf-hello e nf-schema nella Parte 1 funziona anche con il vostro plugin.

---

## 5. Aggiungere decorateGreeting

Un plugin può fornire più funzioni.
Aggiungete una seconda funzione che racchiude un saluto con marcatori decorativi; la renderete configurabile nella Parte 6.

Modificate `GreetingExtension.groovy` per aggiungere `decorateGreeting` dopo `reverseGreeting`, prima della parentesi graffa di chiusura della classe:

=== "Dopo"

    ```groovy title="GreetingExtension.groovy" linenums="24" hl_lines="16-22"
    @CompileStatic
    class GreetingExtension extends PluginExtensionPoint {

        @Override
        protected void init(Session session) {
        }

        /**
         * Inverte una stringa di saluto
         */
        @Function
        String reverseGreeting(String greeting) {
            return greeting.reverse()
        }

        /**
         * Decora un saluto con marcatori celebrativi
         */
        @Function
        String decorateGreeting(String greeting) {
            return "*** ${greeting} ***"             // (1)!
        }

    }
    ```

    1. Interpolazione di stringhe in Groovy: `#!groovy ${...}` inserisce il valore della variabile nella stringa

=== "Prima"

    ```groovy title="GreetingExtension.groovy" linenums="24"
    @CompileStatic
    class GreetingExtension extends PluginExtensionPoint {

        @Override
        protected void init(Session session) {
        }

        /**
         * Inverte una stringa di saluto
         */
        @Function
        String reverseGreeting(String greeting) {
            return greeting.reverse()
        }

    }
    ```

Questa funzione usa l'interpolazione di stringhe di Groovy (`"*** ${greeting} ***"`) per incorporare la variabile del saluto all'interno di una stringa.

Compilate, installate e aggiornate il flusso di lavoro:

```bash
cd nf-greeting && make install && cd ..
```

Aggiornate `greet.nf` per importare e usare anche `decorateGreeting`:

=== "Dopo"

    ```groovy title="greet.nf" hl_lines="4-6 14 16-17 19 33" linenums="1"
    #!/usr/bin/env nextflow

    include { samplesheetToList } from 'plugin/nf-schema'
    // Importa le funzioni personalizzate dal nostro plugin
    include { reverseGreeting } from 'plugin/nf-greeting'
    include { decorateGreeting } from 'plugin/nf-greeting'  // (1)!

    params.input = 'greetings.csv'

    process SAY_HELLO {
        input:
            val greeting
        output:
            path 'greeting.txt'
        script:
        // Usa la nostra funzione personalizzata del plugin per decorare il saluto
        def decorated = decorateGreeting(greeting)  // (2)!
        """
        echo '$decorated' > greeting.txt
        """
    }

    workflow {
        greeting_ch = Channel.fromList(samplesheetToList(params.input, 'greetings_schema.json'))
            .map { row -> row[0] }

        // Dimostra l'uso della funzione reverseGreeting
        greeting_ch
            .map { greeting -> reverseGreeting(greeting) }
            .view { reversed -> "Reversed: $reversed" }

        SAY_HELLO(greeting_ch)
        SAY_HELLO.out.view { file -> "Decorated: ${file.text.trim()}" }
    }
    ```

    1. Più funzioni dallo stesso plugin richiedono istruzioni `include` separate
    2. Le funzioni del plugin funzionano anche all'interno dei blocchi `script:` dei processi

=== "Prima"

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

Le funzioni del plugin funzionano sia negli script dei processi (come `decorateGreeting` all'interno di `SAY_HELLO`) che nelle operazioni del flusso di lavoro (come `reverseGreeting` in un `map`).

---

## Takeaway

Avete imparato che:

- Le funzioni sono definite con l'annotazione `@Function` nelle sottoclassi di `PluginExtensionPoint`
- Le funzioni del plugin importate con `include` funzionano in modo identico sia che provengano dal vostro plugin che da uno esistente
- Le funzioni del plugin funzionano sia negli script dei processi che nelle operazioni del flusso di lavoro

---

## Cosa c'è dopo?

Le vostre funzioni funzionano, ma finora le avete verificate solo eseguendo l'intera pipeline e controllando l'output a occhio.
Questo approccio non è scalabile: man mano che aggiungete più funzioni, avete bisogno di un modo più rapido per verificare che ciascuna si comporti correttamente, specialmente dopo aver apportato modifiche.
La sezione successiva introduce i test unitari, che vi permettono di verificare le singole funzioni automaticamente senza eseguire una pipeline.

[Continua alla Parte 4 :material-arrow-right:](04_build_and_test.md){ .md-button .md-button--primary }
