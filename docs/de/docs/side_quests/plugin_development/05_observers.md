# Teil 5: Trace Observers

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } KI-gestützte Übersetzung - [mehr erfahren & Verbesserungen vorschlagen](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Trace Observers ermöglichen es deinem Plugin, auf Workflow-Ereignisse zu reagieren, z. B. wenn eine Aufgabe abgeschlossen wird, eine Datei veröffentlicht wird oder die Pipeline beendet ist.
Das ermöglicht Anwendungsfälle wie benutzerdefinierte Berichte, Slack-Benachrichtigungen, Metrikenerfassung oder die Integration mit externen Monitoring-Systemen.
In diesem Abschnitt baust du einen Observer, der abgeschlossene Aufgaben zählt und eine Zusammenfassung ausgibt.

!!! tip "Hier eingestiegen?"

    Wenn du erst ab diesem Teil mitmachst, kopiere die Lösung aus Teil 4 als Ausgangspunkt:

    ```bash
    cp -r solutions/4-build-and-test/* .
    ```

---

## 1. Den vorhandenen Trace Observer verstehen

Die Meldung „Pipeline is starting!" beim Ausführen der Pipeline stammt aus der Klasse `GreetingObserver` in deinem Plugin.

Schau dir den Observer-Code an:

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
 * Implementiert einen Observer, der benutzerdefinierte
 * Logik bei Nextflow-Ausführungsereignissen ermöglicht.
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

1. Interface zum Einhängen in Workflow-Lifecycle-Ereignisse
2. Wird beim Start des Workflows aufgerufen; erhält die Session für den Zugriff auf die Konfiguration
3. Wird aufgerufen, wenn der Workflow erfolgreich abgeschlossen wird

Zwei Dinge sind hier wichtig:

1. **`class GreetingObserver implements TraceObserver`**: `TraceObserver` ist ein von Nextflow definiertes Interface. Wenn deine Klasse dieses Interface implementiert, kann Nextflow sich einhängen und deine Methoden aufrufen, wenn Ereignisse eintreten.
2. **`@Override`**: Das `TraceObserver`-Interface definiert Methoden wie `onFlowCreate` und `onFlowComplete`. Wenn du Methoden mit diesen Namen schreibst und die Annotation `@Override` hinzufügst, ruft Nextflow sie zum richtigen Zeitpunkt auf. Methoden, die du nicht überschreibst, werden ignoriert.

Die vollständige Liste der Lifecycle-Ereignisse, in die du dich zum Zeitpunkt der Erstellung einhängen kannst:

| Methode             | Wann sie aufgerufen wird                   |
| ------------------- | ------------------------------------------ |
| `onFlowCreate`      | Workflow startet                           |
| `onFlowComplete`    | Workflow wird beendet                      |
| `onProcessStart`    | Eine Aufgabe beginnt die Ausführung        |
| `onProcessComplete` | Eine Aufgabe wird abgeschlossen            |
| `onProcessCached`   | Eine gecachte Aufgabe wird wiederverwendet |
| `onFilePublish`     | Eine Datei wird veröffentlicht             |

Eine vollständige Liste findest du im [TraceObserver-Interface](https://github.com/nextflow-io/nextflow/blob/master/modules/nextflow/src/main/groovy/nextflow/trace/TraceObserver.groovy) im Nextflow-Quellcode.

---

## 2. Einen Aufgabenzähler-Observer hinzufügen

Das Ziel ist, einen Observer zu bauen, der abgeschlossene Aufgaben zählt und am Ende eine Zusammenfassung ausgibt.
Um einen neuen Observer zu einem Plugin hinzuzufügen, sind zwei Dinge nötig: die Observer-Klasse schreiben und sie in der Factory registrieren, damit Nextflow sie lädt.

### 2.1. Einen minimalen Observer erstellen

Erstelle eine neue Datei:

```bash
touch nf-greeting/src/main/groovy/training/plugin/TaskCounterObserver.groovy
```

Beginne mit dem einfachstmöglichen Observer, der eine Meldung ausgibt, wenn eine Aufgabe abgeschlossen wird:

```groovy title="nf-greeting/src/main/groovy/training/plugin/TaskCounterObserver.groovy" linenums="1"
package training.plugin

import groovy.transform.CompileStatic
import nextflow.processor.TaskHandler       // (1)!
import nextflow.trace.TraceObserver
import nextflow.trace.TraceRecord

/**
 * Observer, der auf den Abschluss von Aufgaben reagiert
 */
@CompileStatic
class TaskCounterObserver implements TraceObserver {  // (2)!

    @Override
    void onProcessComplete(TaskHandler handler, TraceRecord trace) {  // (3)!
        println "✓ Task completed!"
    }
}
```

1. Die benötigten Klassen importieren: `TraceObserver`, `TaskHandler` und `TraceRecord`
2. Eine Klasse erstellen, die `TraceObserver` implementiert
3. `onProcessComplete` überschreiben, um Code auszuführen, wenn eine Aufgabe abgeschlossen wird

Das ist das Minimum:

- Die benötigten Klassen importieren (`TraceObserver`, `TaskHandler`, `TraceRecord`)
- Eine Klasse erstellen, die `TraceObserver` implementiert
- `onProcessComplete` überschreiben, um etwas zu tun, wenn eine Aufgabe abgeschlossen wird

### 2.2. Den Observer registrieren

Die `GreetingFactory` erstellt Observers.
Schau sie dir an:

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

Bearbeite `GreetingFactory.groovy`, um den neuen Observer hinzuzufügen:

=== "Danach"

    ```groovy title="GreetingFactory.groovy" linenums="31" hl_lines="3-6"
    @Override
    Collection<TraceObserver> create(Session session) {
        return [
            new GreetingObserver(),
            new TaskCounterObserver()
        ]
    }
    ```

=== "Vorher"

    ```groovy title="GreetingFactory.groovy" linenums="31" hl_lines="3"
    @Override
    Collection<TraceObserver> create(Session session) {
        return List.<TraceObserver>of(new GreetingObserver())
    }
    ```

!!! note "Groovy-Listen-Syntax"

    Wir haben das Java-ähnliche `List.<TraceObserver>of(...)` durch Groovys einfachere Listen-Literal-Syntax `[...]` ersetzt.
    Beide geben eine `Collection` zurück, aber die Groovy-Syntax ist lesbarer, wenn mehrere Elemente hinzugefügt werden.

### 2.3. Bauen, installieren und testen

```bash
cd nf-greeting && make install && cd ..
nextflow run greet.nf -ansi-log false
```

!!! tip "Warum `-ansi-log false`?"

    Standardmäßig überschreibt Nextflows ANSI-Fortschrittsanzeige vorherige Zeilen, um eine saubere, aktualisierte Ansicht des Fortschritts zu zeigen.
    Das bedeutet, du würdest nur den *letzten* Aufgabenzähler sehen, nicht die Zwischenmeldungen.

    Mit `-ansi-log false` wird dieses Verhalten deaktiviert und alle Ausgaben werden sequenziell angezeigt – das ist wichtig beim Testen von Observers, die während der Ausführung Meldungen ausgeben.

Du solltest „✓ Task completed!" fünfmal sehen (einmal pro Aufgabe), vermischt mit der bestehenden Pipeline-Ausgabe:

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

Der Observer funktioniert.
Jedes Mal, wenn eine Aufgabe abgeschlossen wird, ruft Nextflow `onProcessComplete` auf, und unsere Implementierung gibt eine Meldung aus.

??? exercise "Die Meldung anpassen"

    Ändere die Meldung in `onProcessComplete` nach deinen Wünschen, baue neu und führe die Pipeline erneut aus.
    Das bestätigt, dass der vollständige Bearbeiten-Bauen-Ausführen-Zyklus für Observers funktioniert.

### 2.4. Zähllogik hinzufügen

Der minimale Observer beweist, dass der Hook funktioniert, verfolgt aber nichts.

Eine Klasse kann Variablen (sogenannte Felder oder Instanzvariablen) enthalten, die für die Lebensdauer des Objekts erhalten bleiben.
Das bedeutet, ein Observer kann über mehrere Ereignisse hinweg während eines Pipeline-Laufs Zustand ansammeln.

Die nächste Version fügt eine Zählervariable (`taskCount`) hinzu, die bei null beginnt.
Jedes Mal, wenn eine Aufgabe abgeschlossen wird, erhöht sich der Zähler um eins.
Wenn der gesamte Workflow abgeschlossen ist, gibt der Observer die Gesamtsumme aus.

Aktualisiere `TaskCounterObserver.groovy` mit den hervorgehobenen Änderungen:

```groovy title="nf-greeting/src/main/groovy/training/plugin/TaskCounterObserver.groovy" linenums="1" hl_lines="14 18-19 22-24"
package training.plugin

import groovy.transform.CompileStatic
import nextflow.processor.TaskHandler
import nextflow.trace.TraceObserver
import nextflow.trace.TraceRecord

/**
 * Observer, der abgeschlossene Aufgaben zählt
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

1. `taskCount` ist eine Variable, die zum Observer-Objekt gehört. Sie behält ihren Wert zwischen Methodenaufrufen, sodass sie über den gesamten Workflow-Lauf hinweg einen Zähler ansammeln kann. `private` bedeutet, dass nur diese Klasse darauf zugreifen kann.
2. `taskCount++` erhöht den Zähler um eins. Diese Zeile wird jedes Mal ausgeführt, wenn eine Aufgabe abgeschlossen wird, sodass der Zähler mit dem Fortschritt des Workflows wächst.
3. `onFlowComplete` ist ein zweiter Lifecycle-Hook. Er wird einmal ausgeführt, wenn der Workflow abgeschlossen ist – ein guter Ort, um eine Zusammenfassung auszugeben.

Zusammenfassend:

- `taskCount` bleibt über Methodenaufrufe hinweg erhalten und sammelt einen Zähler über den gesamten Lauf an
- `onProcessComplete` erhöht den Zähler und gibt die laufende Summe aus, jedes Mal wenn eine Aufgabe abgeschlossen wird
- `onFlowComplete` wird einmal am Ende ausgeführt und gibt den endgültigen Zählerstand aus

Neu bauen und testen:

```bash
cd nf-greeting && make install && cd ..
nextflow run greet.nf -ansi-log false
```

??? example "Ausgabe"

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

    Die Zählermeldungen sind mit den Aufgaben-Einreichungen vermischt, weil Observers ausgeführt werden, wenn Aufgaben abgeschlossen werden.

---

## 3. Veröffentlichte Dateien verfolgen

Der Observer kann auch reagieren, wenn Dateien veröffentlicht werden.
Die Methode `onFilePublish` erhält die Ziel- und Quellpfade, die du zum Protokollieren, Validieren oder Verarbeiten veröffentlichter Ausgaben verwenden kannst.

### 3.1. Ein Ausgabeverzeichnis hinzufügen

Aktualisiere zunächst `greet.nf`, damit der Prozess `SAY_HELLO` seine Ausgabedateien veröffentlicht:

=== "Danach"

    ```groovy title="greet.nf" linenums="10" hl_lines="2"
    process SAY_HELLO {
        publishDir 'results'
        input:
            val greeting
        output:
            path 'greeting.txt'
        script:
        // Unsere benutzerdefinierte Plugin-Funktion verwenden, um die Begrüßung zu dekorieren
        def decorated = decorateGreeting(greeting)
        """
        echo '$decorated' > greeting.txt
        """
    }
    ```

=== "Vorher"

    ```groovy title="greet.nf" linenums="10"
    process SAY_HELLO {
        input:
            val greeting
        output:
            path 'greeting.txt'
        script:
        // Unsere benutzerdefinierte Plugin-Funktion verwenden, um die Begrüßung zu dekorieren
        def decorated = decorateGreeting(greeting)
        """
        echo '$decorated' > greeting.txt
        """
    }
    ```

### 3.2. Die Methode onFilePublish hinzufügen

Füge eine `onFilePublish`-Methode und den benötigten Import zu `TaskCounterObserver.groovy` hinzu:

```groovy title="nf-greeting/src/main/groovy/training/plugin/TaskCounterObserver.groovy" linenums="1" hl_lines="5 23-26"
package training.plugin

import groovy.transform.CompileStatic
import nextflow.processor.TaskHandler
import java.nio.file.Path
import nextflow.trace.TraceObserver
import nextflow.trace.TraceRecord

/**
 * Observer, der abgeschlossene Aufgaben zählt
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

### 3.3. Bauen und testen

```bash
cd nf-greeting && make install && cd ..
nextflow run greet.nf -ansi-log false
```

Du solltest „Published:"-Meldungen für jede Ausgabedatei zusammen mit der Aufgabenzähler-Ausgabe sehen:

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

Die Methode `onFilePublish` wird jedes Mal ausgelöst, wenn Nextflow eine Datei im Verzeichnis `results` veröffentlicht.
Dieses Muster ist nützlich für das Erstellen von Audit-Logs, das Auslösen nachgelagerter Aktionen oder das Validieren von Ausgaben während ihrer Erstellung.

---

## Fazit

Du hast gelernt, dass:

- Trace Observers sich in Workflow-Lifecycle-Ereignisse wie `onFlowCreate`, `onProcessComplete`, `onFilePublish` und `onFlowComplete` einhängen
- Observers erstellt werden, indem `TraceObserver` implementiert und in einer Factory registriert wird
- Observers Instanzvariablen halten können, um Zustand über Ereignisse hinweg anzusammeln
- Observers nützlich für benutzerdefiniertes Logging, Metrikenerfassung, Benachrichtigungen und Berichte sind

---

## Wie geht es weiter?

Der Aufgabenzähler funktioniert, ist aber immer aktiv.
In einem echten Plugin sollten Benutzer\*innen Funktionen aktivieren oder deaktivieren oder das Verhalten über `nextflow.config` anpassen können, ohne den Plugin-Quellcode zu bearbeiten.
Der nächste Abschnitt zeigt, wie du deinen Observer konfigurierbar machst und wie du dein fertiges Plugin mit anderen teilst.

[Weiter zu Teil 6 :material-arrow-right:](06_configuration.md){ .md-button .md-button--primary }
