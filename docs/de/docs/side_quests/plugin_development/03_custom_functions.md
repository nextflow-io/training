# Teil 3: Eigene Funktionen

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } KI-gestützte Übersetzung - [mehr erfahren & Verbesserungen vorschlagen](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Am Ende dieses Abschnitts hast du eigene Funktionen in deinem Plugin, die lokal gebaut und installiert sind und in einem echten Workflow laufen.

!!! tip "Hier eingestiegen?"

    Wenn du erst ab diesem Teil mitmachst, kopiere die Lösung aus Teil 2 als Ausgangspunkt:

    ```bash
    cp -r solutions/2-create-project/* .
    ```

---

## 1. Sieh dir an, was das Template generiert hat

Bevor du eigene Funktionen schreibst, schau dir die Beispielfunktion an, die das Template erstellt hat, um das Muster zu verstehen.

Wechsle in das Plugin-Verzeichnis:

```bash
cd nf-greeting
```

Das Template hat eine Datei namens `GreetingExtension.groovy` erstellt, in der Plugin-Funktionen definiert werden.
Öffne sie, um den Ausgangspunkt zu sehen:

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
 * Implementiert eine eigene Funktion, die von
 * Nextflow-Skripten importiert werden kann.
 */
@CompileStatic
class GreetingExtension extends PluginExtensionPoint { // (1)!

    @Override
    protected void init(Session session) {             // (2)!
    }

    /**
     * Begrüßt das angegebene Ziel.
     *
     * @param target
     */
    @Function                                          // (3)!
    void sayHello(String target) {
        println "Hello, ${target}!"
    }

}
```

1. Die Klasse, auf der deine Extension aufbaut. Nextflow benötigt sie, um deine Funktionen zu erkennen.
2. Wird beim Laden des Plugins aufgerufen; hier kannst du Initialisierungen vornehmen
3. Macht diese Methode über `include` aus Workflows heraus aufrufbar

Das Template enthält eine Beispielfunktion `sayHello`.
Die Annotation `@Function` macht eine Methode aus Nextflow-Workflows heraus aufrufbar.
Ohne sie existiert die Methode nur innerhalb des Plugin-Codes.

In Groovy (und Java) deklarieren Methoden, welchen Typ sie zurückgeben und welche Typen ihre Parameter haben.
Zum Beispiel deklariert `String reverseGreeting(String greeting)` eine Methode, die einen `String`-Parameter entgegennimmt und einen `String` zurückgibt.
Das Schlüsselwort `void` bedeutet, dass die Methode nichts zurückgibt, wie bei `sayHello` oben.
Das unterscheidet sich von Python oder R, wo Typen nicht explizit deklariert werden müssen.

---

## 2. sayHello durch reverseGreeting ersetzen

Die `sayHello`-Funktion des Templates ist ein Platzhalter.
Ersetze sie durch deine eigene Funktion, um den vollständigen Zyklus aus Schreiben, Bauen und Verwenden einer Plugin-Funktion zu durchlaufen.

Bearbeite `src/main/groovy/training/plugin/GreetingExtension.groovy`, um die `sayHello`-Methode zu ersetzen:

=== "Danach"

    ```groovy title="GreetingExtension.groovy" linenums="24" hl_lines="8-14"
    @CompileStatic
    class GreetingExtension extends PluginExtensionPoint {

        @Override
        protected void init(Session session) {
        }

        /**
         * Einen Begrüßungsstring umkehren
         */
        @Function                                // (1)!
        String reverseGreeting(String greeting) { // (2)!
            return greeting.reverse()             // (3)!
        }

    }
    ```

    1. Macht die Methode aus Nextflow-Workflows heraus aufrufbar
    2. Nimmt einen String entgegen, gibt einen String zurück
    3. Grooovys eingebaute Methode zur String-Umkehrung

=== "Vorher"

    ```groovy title="GreetingExtension.groovy" linenums="24" hl_lines="12-20"
    /**
     * Implementiert eine eigene Funktion, die von
     * Nextflow-Skripten importiert werden kann.
     */
    @CompileStatic
    class GreetingExtension extends PluginExtensionPoint {

        @Override
        protected void init(Session session) {
        }

        /**
         * Begrüßt das angegebene Ziel.
         *
         * @param target
         */
        @Function
        void sayHello(String target) {
            println "Hello, ${target}!"
        }

    }
    ```

Die wichtigsten Teile dieser Funktion:

- **`@Function`**: Macht die Methode aus Nextflow-Workflows heraus aufrufbar
- **`String reverseGreeting(String greeting)`**: Nimmt einen String entgegen, gibt einen String zurück
- **`greeting.reverse()`**: Grooovys eingebaute Methode zur String-Umkehrung

!!! tip "Öffentliche und private Methoden"

    Methoden ohne `@Function` werden Nextflow-Workflows nicht zugänglich gemacht.
    Du kannst deiner Klasse Hilfsmethoden hinzufügen, ohne dir Sorgen zu machen, dass sie in den Workflow-Namespace gelangen.

---

## 3. Dein Plugin bauen und installieren

Baue und installiere das Plugin:

```bash
make install
```

!!! tip "Falls der Build fehlschlägt"

    Lies die Fehlermeldung sorgfältig; sie enthält meist eine Zeilennummer und beschreibt das Problem.
    Häufige Ursachen sind Syntaxfehler (fehlende Klammer oder Anführungszeichen), falsch geschriebene Klassennamen und Typkonflikte.
    Wenn du nicht weiterkommst, vergleiche deinen Code Zeichen für Zeichen mit den Beispielen.

---

## 4. Deine Funktion in einem Workflow verwenden

Das Plugin ist gebaut und installiert.
Der nächste Schritt ist, `reverseGreeting` in einem Workflow zu verwenden, um zu überprüfen, ob es von Anfang bis Ende funktioniert.

Wechsle zurück in das Pipeline-Verzeichnis:

```bash
cd ..
```

Bearbeite `greet.nf`, um `reverseGreeting` zu importieren und zu verwenden:

=== "Danach"

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

=== "Vorher"

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

Führe die Pipeline aus:

```bash
nextflow run greet.nf
```

??? example "Ausgabe"

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

Deine erste eigene Plugin-Funktion funktioniert in einem echten Workflow.
Das gleiche `include { ... } from 'plugin/...'`-Muster, das du in Teil 1 mit nf-hello und nf-schema verwendet hast, funktioniert auch mit deinem eigenen Plugin.

---

## 5. decorateGreeting hinzufügen

Ein Plugin kann mehrere Funktionen bereitstellen.
Füge eine zweite hinzu, die eine Begrüßung mit dekorativen Markierungen umschließt; in Teil 6 wirst du sie konfigurierbar machen.

Bearbeite `GreetingExtension.groovy`, um `decorateGreeting` nach `reverseGreeting` und vor der schließenden Klammer der Klasse hinzuzufügen:

=== "Danach"

    ```groovy title="GreetingExtension.groovy" linenums="24" hl_lines="16-22"
    @CompileStatic
    class GreetingExtension extends PluginExtensionPoint {

        @Override
        protected void init(Session session) {
        }

        /**
         * Einen Begrüßungsstring umkehren
         */
        @Function
        String reverseGreeting(String greeting) {
            return greeting.reverse()
        }

        /**
         * Eine Begrüßung mit festlichen Markierungen versehen
         */
        @Function
        String decorateGreeting(String greeting) {
            return "*** ${greeting} ***"             // (1)!
        }

    }
    ```

    1. Groovy-String-Interpolation: `#!groovy ${...}` fügt den Wert der Variable in den String ein

=== "Vorher"

    ```groovy title="GreetingExtension.groovy" linenums="24"
    @CompileStatic
    class GreetingExtension extends PluginExtensionPoint {

        @Override
        protected void init(Session session) {
        }

        /**
         * Einen Begrüßungsstring umkehren
         */
        @Function
        String reverseGreeting(String greeting) {
            return greeting.reverse()
        }

    }
    ```

Diese Funktion verwendet Groovy-String-Interpolation (`"*** ${greeting} ***"`), um die Begrüßungsvariable in einen String einzubetten.

Baue und installiere das Plugin und aktualisiere den Workflow:

```bash
cd nf-greeting && make install && cd ..
```

Aktualisiere `greet.nf`, um auch `decorateGreeting` zu importieren und zu verwenden:

=== "Danach"

    ```groovy title="greet.nf" hl_lines="4-6 14 16-17 19 33" linenums="1"
    #!/usr/bin/env nextflow

    include { samplesheetToList } from 'plugin/nf-schema'
    // Eigene Funktionen aus unserem Plugin importieren
    include { reverseGreeting } from 'plugin/nf-greeting'
    include { decorateGreeting } from 'plugin/nf-greeting'  // (1)!

    params.input = 'greetings.csv'

    process SAY_HELLO {
        input:
            val greeting
        output:
            path 'greeting.txt'
        script:
        // Unsere eigene Plugin-Funktion verwenden, um die Begrüßung zu dekorieren
        def decorated = decorateGreeting(greeting)  // (2)!
        """
        echo '$decorated' > greeting.txt
        """
    }

    workflow {
        greeting_ch = Channel.fromList(samplesheetToList(params.input, 'greetings_schema.json'))
            .map { row -> row[0] }

        // Verwendung der reverseGreeting-Funktion demonstrieren
        greeting_ch
            .map { greeting -> reverseGreeting(greeting) }
            .view { reversed -> "Reversed: $reversed" }

        SAY_HELLO(greeting_ch)
        SAY_HELLO.out.view { file -> "Decorated: ${file.text.trim()}" }
    }
    ```

    1. Mehrere Funktionen aus demselben Plugin benötigen separate `include`-Anweisungen
    2. Plugin-Funktionen funktionieren auch innerhalb von `script:`-Blöcken eines Prozesses

=== "Vorher"

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

??? example "Ausgabe"

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

Plugin-Funktionen funktionieren sowohl in Prozess-Skripten (wie `decorateGreeting` innerhalb von `SAY_HELLO`) als auch in Workflow-Operationen (wie `reverseGreeting` in einem `map`).

---

## Fazit

Du hast gelernt, dass:

- Funktionen mit der Annotation `@Function` in `PluginExtensionPoint`-Unterklassen definiert werden
- Mit `include` importierte Plugin-Funktionen identisch funktionieren, egal ob sie aus deinem eigenen Plugin oder einem bestehenden stammen
- Plugin-Funktionen sowohl in Prozess-Skripten als auch in Workflow-Operationen funktionieren

---

## Wie geht es weiter?

Deine Funktionen funktionieren, aber bisher hast du das nur überprüft, indem du die vollständige Pipeline ausgeführt und die Ausgabe manuell kontrolliert hast.
Dieser Ansatz skaliert nicht: Wenn du mehr Funktionen hinzufügst, brauchst du eine schnellere Möglichkeit zu prüfen, dass jede einzelne korrekt funktioniert – besonders nach Änderungen.
Der nächste Abschnitt stellt Unit-Tests vor, mit denen du einzelne Funktionen automatisch überprüfen kannst, ohne eine Pipeline auszuführen.

[Weiter zu Teil 4 :material-arrow-right:](04_build_and_test.md){ .md-button .md-button--primary }
