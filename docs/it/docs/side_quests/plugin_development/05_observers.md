# Parte 5: Trace Observer

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduzione assistita da IA - [scopri di più e suggerisci miglioramenti](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

I trace observer permettono al tuo plugin di rispondere agli eventi del flusso di lavoro, come il completamento di un'attività, la pubblicazione di un file o il termine della pipeline.
Questo abilita casi d'uso come report personalizzati, notifiche Slack, raccolta di metriche o integrazione con sistemi di monitoraggio esterni.
In questa sezione, costruiremo un observer che conta le attività completate e stampa un riepilogo.

!!! tip "Parti da qui?"

    Se ti unisci a questa parte, copia la soluzione dalla Parte 4 da usare come punto di partenza:

    ```bash
    cp -r solutions/4-build-and-test/* .
    ```

---

## 1. Comprendere il trace observer esistente

Il messaggio "Pipeline is starting!" che hai visto quando hai eseguito la pipeline proveniva dalla classe `GreetingObserver` nel tuo plugin.

Esaminiamo il codice dell'observer:

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
 * Implementa un observer che permette di definire logica personalizzata
 * sugli eventi di esecuzione di Nextflow.
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

1. Interfaccia per agganciarsi agli eventi del ciclo di vita del flusso di lavoro
2. Chiamato quando il flusso di lavoro si avvia; riceve la sessione per accedere alla configurazione
3. Chiamato quando il flusso di lavoro termina con successo

Ci sono due cose da notare:

1. **`class GreetingObserver implements TraceObserver`**: `TraceObserver` è un'interfaccia definita da Nextflow. Se la tua classe implementa questa interfaccia, Nextflow può agganciarsi ad essa e chiamare i tuoi metodi quando si verificano gli eventi.
2. **`@Override`**: L'interfaccia `TraceObserver` definisce metodi come `onFlowCreate` e `onFlowComplete`. Quando scrivi metodi con questi nomi e aggiungi l'annotazione `@Override`, Nextflow li chiama al momento opportuno. I metodi che non sovrascrivi vengono ignorati.

L'insieme completo degli eventi del ciclo di vita a cui puoi agganciarti al momento della stesura è:

| Metodo              | Quando viene chiamato                   |
| ------------------- | --------------------------------------- |
| `onFlowCreate`      | Il flusso di lavoro si avvia            |
| `onFlowComplete`    | Il flusso di lavoro termina             |
| `onProcessStart`    | Un'attività inizia l'esecuzione         |
| `onProcessComplete` | Un'attività termina                     |
| `onProcessCached`   | Un'attività in cache viene riutilizzata |
| `onFilePublish`     | Un file viene pubblicato                |

Per un elenco completo, consulta l'[interfaccia TraceObserver](https://github.com/nextflow-io/nextflow/blob/master/modules/nextflow/src/main/groovy/nextflow/trace/TraceObserver.groovy) nel codice sorgente di Nextflow.

---

## 2. Aggiungere un observer con contatore di attività

L'obiettivo è costruire un observer che conti le attività completate e stampi un riepilogo alla fine.
Aggiungere un nuovo observer a un plugin richiede due cose: scrivere la classe dell'observer e registrarla nella factory in modo che Nextflow la carichi.

### 2.1. Creare un observer minimale

Creiamo un nuovo file:

```bash
touch nf-greeting/src/main/groovy/training/plugin/TaskCounterObserver.groovy
```

Iniziamo con l'observer più semplice possibile, che stampa un messaggio quando un'attività si completa:

```groovy title="nf-greeting/src/main/groovy/training/plugin/TaskCounterObserver.groovy" linenums="1"
package training.plugin

import groovy.transform.CompileStatic
import nextflow.processor.TaskHandler       // (1)!
import nextflow.trace.TraceObserver
import nextflow.trace.TraceRecord

/**
 * Observer che risponde al completamento delle attività
 */
@CompileStatic
class TaskCounterObserver implements TraceObserver {  // (2)!

    @Override
    void onProcessComplete(TaskHandler handler, TraceRecord trace) {  // (3)!
        println "✓ Task completed!"
    }
}
```

1. Importa le classi necessarie: `TraceObserver`, `TaskHandler` e `TraceRecord`
2. Crea una classe che `implements TraceObserver`
3. Sovrascrive `onProcessComplete` per eseguire codice quando un'attività termina

Questo è il minimo necessario:

- Importare le classi richieste (`TraceObserver`, `TaskHandler`, `TraceRecord`)
- Creare una classe che `implements TraceObserver`
- Sovrascrivere `onProcessComplete` per fare qualcosa quando un'attività termina

### 2.2. Registrare l'observer

La `GreetingFactory` crea gli observer.
Diamo un'occhiata:

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

Modifichiamo `GreetingFactory.groovy` per aggiungere il nuovo observer:

=== "Dopo"

    ```groovy title="GreetingFactory.groovy" linenums="31" hl_lines="3-6"
    @Override
    Collection<TraceObserver> create(Session session) {
        return [
            new GreetingObserver(),
            new TaskCounterObserver()
        ]
    }
    ```

=== "Prima"

    ```groovy title="GreetingFactory.groovy" linenums="31" hl_lines="3"
    @Override
    Collection<TraceObserver> create(Session session) {
        return List.<TraceObserver>of(new GreetingObserver())
    }
    ```

!!! note "Sintassi delle liste in Groovy"

    Abbiamo sostituito lo stile Java `List.<TraceObserver>of(...)` con il più semplice letterale di lista di Groovy `[...]`.
    Entrambi restituiscono una `Collection`, ma la sintassi Groovy è più leggibile quando si aggiungono più elementi.

### 2.3. Build, installazione e test

```bash
cd nf-greeting && make install && cd ..
nextflow run greet.nf -ansi-log false
```

!!! tip "Perché `-ansi-log false`?"

    Per impostazione predefinita, la visualizzazione ANSI di Nextflow sovrascrive le righe precedenti per mostrare una vista aggiornata e pulita del progresso.
    Questo significa che vedresti solo il conteggio *finale* delle attività, non i messaggi intermedi.

    Usare `-ansi-log false` disabilita questo comportamento e mostra tutto l'output in sequenza, il che è essenziale quando si testano observer che stampano messaggi durante l'esecuzione.

Dovresti vedere "✓ Task completed!" stampato cinque volte (una per ogni attività), intercalato con l'output esistente della pipeline:

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

L'observer funziona.
Ogni volta che un'attività termina, Nextflow chiama `onProcessComplete` e la nostra implementazione stampa un messaggio.

??? exercise "Personalizza il messaggio"

    Prova a cambiare il messaggio in `onProcessComplete` con qualcosa di tuo, ricompila ed esegui di nuovo.
    Questo conferma che il ciclo completo modifica-build-esecuzione funziona per gli observer.

### 2.4. Aggiungere la logica di conteggio

L'observer minimale dimostra che l'hook funziona, ma non tiene traccia di nulla.

Una classe può contenere variabili (chiamate campi o variabili di istanza) che persistono per tutta la durata dell'oggetto.
Questo significa che un observer può accumulare stato attraverso più eventi durante l'esecuzione di una pipeline.

La versione successiva aggiunge una variabile contatore (`taskCount`) che parte da zero.
Ogni volta che un'attività si completa, il contatore aumenta di uno.
Quando l'intero flusso di lavoro termina, l'observer stampa il totale finale.

Aggiorniamo `TaskCounterObserver.groovy` con le modifiche evidenziate:

```groovy title="nf-greeting/src/main/groovy/training/plugin/TaskCounterObserver.groovy" linenums="1" hl_lines="14 18-19 22-24"
package training.plugin

import groovy.transform.CompileStatic
import nextflow.processor.TaskHandler
import nextflow.trace.TraceObserver
import nextflow.trace.TraceRecord

/**
 * Observer che conta le attività completate
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

1. `taskCount` è una variabile che appartiene all'oggetto observer. Mantiene il suo valore tra le chiamate ai metodi, così può accumulare un conteggio per l'intera esecuzione del flusso di lavoro. `private` significa che solo questa classe può accedervi.
2. `taskCount++` aggiunge uno al contatore. Questa riga viene eseguita ogni volta che un'attività si completa, quindi il conteggio cresce man mano che il flusso di lavoro avanza.
3. `onFlowComplete` è un secondo hook del ciclo di vita. Viene eseguito una volta quando il flusso di lavoro termina, rendendolo un buon posto per stampare un riepilogo.

In sintesi:

- `taskCount` persiste tra le chiamate ai metodi, accumulando un conteggio per l'intera esecuzione
- `onProcessComplete` incrementa il contatore e stampa il totale aggiornato ogni volta che un'attività termina
- `onFlowComplete` viene eseguito una volta alla fine, stampando il conteggio finale

Ricompiliamo e testiamo:

```bash
cd nf-greeting && make install && cd ..
nextflow run greet.nf -ansi-log false
```

??? example "Output"

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

    I messaggi del contatore sono intercalati con le sottomissioni delle attività perché gli observer vengono eseguiti al completamento delle attività.

---

## 3. Tracciare i file pubblicati

L'observer può anche rispondere quando i file vengono pubblicati.
Il metodo `onFilePublish` riceve i percorsi di destinazione e sorgente, che puoi usare per registrare, validare o elaborare gli output pubblicati.

### 3.1. Aggiungere una directory di pubblicazione

Prima di tutto, aggiorniamo `greet.nf` in modo che il processo `SAY_HELLO` pubblichi i suoi file di output:

=== "Dopo"

    ```groovy title="greet.nf" linenums="10" hl_lines="2"
    process SAY_HELLO {
        publishDir 'results'
        input:
            val greeting
        output:
            path 'greeting.txt'
        script:
        // Usa la nostra funzione personalizzata del plugin per decorare il saluto
        def decorated = decorateGreeting(greeting)
        """
        echo '$decorated' > greeting.txt
        """
    }
    ```

=== "Prima"

    ```groovy title="greet.nf" linenums="10"
    process SAY_HELLO {
        input:
            val greeting
        output:
            path 'greeting.txt'
        script:
        // Usa la nostra funzione personalizzata del plugin per decorare il saluto
        def decorated = decorateGreeting(greeting)
        """
        echo '$decorated' > greeting.txt
        """
    }
    ```

### 3.2. Aggiungere il metodo onFilePublish

Aggiungiamo un metodo `onFilePublish` e l'import necessario a `TaskCounterObserver.groovy`:

```groovy title="nf-greeting/src/main/groovy/training/plugin/TaskCounterObserver.groovy" linenums="1" hl_lines="5 23-26"
package training.plugin

import groovy.transform.CompileStatic
import nextflow.processor.TaskHandler
import java.nio.file.Path
import nextflow.trace.TraceObserver
import nextflow.trace.TraceRecord

/**
 * Observer che conta le attività completate
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

### 3.3. Build e test

```bash
cd nf-greeting && make install && cd ..
nextflow run greet.nf -ansi-log false
```

Dovresti vedere i messaggi "Published:" per ogni file di output insieme all'output del contatore di attività:

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

Il metodo `onFilePublish` si attiva ogni volta che Nextflow pubblica un file nella directory `results`.
Questo schema è utile per costruire log di audit, attivare azioni downstream o validare gli output man mano che vengono prodotti.

---

## Takeaway

Abbiamo imparato che:

- I trace observer si agganciano agli eventi del ciclo di vita del flusso di lavoro come `onFlowCreate`, `onProcessComplete`, `onFilePublish` e `onFlowComplete`
- Si creano observer implementando `TraceObserver` e registrandoli in una Factory
- Gli observer possono contenere variabili di istanza per accumulare stato attraverso gli eventi
- Gli observer sono utili per logging personalizzato, raccolta di metriche, notifiche e report

---

## Cosa c'è dopo?

Il contatore di attività funziona, ma è sempre attivo.
In un plugin reale, gli utenti dovrebbero poter abilitare o disabilitare funzionalità, o modificare il comportamento, da `nextflow.config` senza dover modificare il codice sorgente del plugin.
La prossima sezione mostra come rendere il tuo observer configurabile e come condividere il plugin finito con altri.

[Continua alla Parte 6 :material-arrow-right:](06_configuration.md){ .md-button .md-button--primary }
