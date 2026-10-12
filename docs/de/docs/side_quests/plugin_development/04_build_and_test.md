# Teil 4: Testing

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } KI-gestützte Übersetzung - [mehr erfahren & Verbesserungen vorschlagen](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Plugins sind eigenständige Software, der Pipeline-Entwickler\*innen vertrauen müssen.
Jedes Feature unabhängig zu testen – außerhalb einer Pipeline – stellt sicher, dass das Plugin korrekt funktioniert, bevor es jemand in einen Workflow integriert.
In diesem Abschnitt schreibst und führst du Tests mit dem Spock-Testing-Framework aus.

!!! tip "Hier eingestiegen?"

    Wenn du mit diesem Teil beginnst, kopiere die Lösung aus Teil 3 als Ausgangspunkt:

    ```bash
    cp -r solutions/3-custom-functions/* .
    ```

    Wechsle dann in das Plugin-Verzeichnis:

    ```bash
    cd nf-greeting
    ```

Stelle sicher, dass du dich im Plugin-Verzeichnis befindest:

```bash
cd nf-greeting
```

---

## 1. Warum testen?

Ein erfolgreicher Build bedeutet, dass der Code kompiliert – aber nicht, dass er wie erwartet funktioniert.
Unit-Tests sind kleine Code-Stücke, die automatisch prüfen, ob deine Funktionen für eine bestimmte Eingabe die richtige Ausgabe liefern.
Ein Test könnte zum Beispiel prüfen, dass `#!groovy reverseGreeting("Hello")` den Wert `"olleH"` zurückgibt.

Tests sind wertvoll, weil sie:

- Fehler finden, bevor Nutzer\*innen sie entdecken
- dir die Sicherheit geben, Änderungen vorzunehmen, ohne etwas kaputtzumachen
- als Dokumentation dienen und zeigen, wie Funktionen verwendet werden sollen

---

## 2. Spock-Tests verstehen

Das Plugin-Template verwendet [Spock](https://spockframework.org/), ein Testing-Framework für Groovy.
Spock ist bereits im Projekt konfiguriert (über `build.gradle`), du musst also nichts hinzufügen.

Wenn du bereits Testing-Tools kennst (wie `pytest` in Python oder `testthat` in R), erfüllt Spock dieselbe Rolle: Du schreibst kleine Funktionen, die deinen Code mit bekannten Eingaben aufrufen und die Ausgaben prüfen.
Der Unterschied ist, dass Spock beschriftete Blöcke verwendet (`given:`, `expect:`, `when:`, `then:`), die einem Nextflow-Prozess oder -Workflow ähneln.

Hier ist die grundlegende Struktur:

```groovy
def 'should reverse a greeting'() {   // (1)!
    given:                             // (2)!
    def ext = new GreetingExtension()

    expect:                            // (3)!
    ext.reverseGreeting('Hello') == 'olleH'
}
```

1. **Testname in Anführungszeichen**: Beschreibt, was der Test prüft. Verwende klares Englisch.
2. **`given:`-Block**: Richte alles ein, was du für den Test brauchst (Objekte erstellen, Daten vorbereiten)
3. **`expect:`-Block**: Die eigentlichen Prüfungen. Jede Zeile muss `true` sein, damit der Test besteht

Diese Struktur macht Tests lesbar: „Gegeben ein Extension-Objekt, erwarte, dass `reverseGreeting('Hello')` gleich `'olleH'` ist."

---

## 3. Tests schreiben

Schreibe Tests für die beiden Funktionen, die du in Teil 3 erstellt hast: `reverseGreeting` und `decorateGreeting`.

### 3.1. Die Test-Klasse erstellen

```bash
touch src/test/groovy/training/plugin/GreetingExtensionTest.groovy
```

Öffne die Datei in deinem Editor und füge das leere Test-Klassen-Gerüst hinzu:

```groovy title="src/test/groovy/training/plugin/GreetingExtensionTest.groovy" linenums="1"
package training.plugin

import spock.lang.Specification

/**
 * Tests für die Greeting-Extension-Funktionen
 */
class GreetingExtensionTest extends Specification {  // (1)!

}
```

1. Alle Spock-Test-Klassen erweitern `Specification`. Das ist der Ausgangspunkt für jede Spock-Testdatei.

### 3.2. reverseGreeting testen

Füge eine Testmethode innerhalb des Klassen-Rumpfs hinzu.
Der `given:`-Block erstellt eine `GreetingExtension`-Instanz, und der `expect:`-Block prüft, dass `reverseGreeting` zwei verschiedene Eingaben korrekt umkehrt.
Damit wird die Funktion direkt getestet, ohne eine Pipeline auszuführen.

=== "Danach"

    ```groovy title="GreetingExtensionTest.groovy" linenums="1" hl_lines="10-17"
    package training.plugin

    import spock.lang.Specification

    /**
     * Tests für die Greeting-Extension-Funktionen
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

    1. Erstelle eine Instanz deiner Extension, um sie direkt zu testen – ohne eine Pipeline auszuführen
    2. Jede Zeile in `expect:` ist eine Assertion; der Test besteht nur, wenn alle `true` sind

=== "Vorher"

    ```groovy title="GreetingExtensionTest.groovy" linenums="1"
    package training.plugin

    import spock.lang.Specification

    /**
     * Tests für die Greeting-Extension-Funktionen
     */
    class GreetingExtensionTest extends Specification {

    }
    ```

### 3.3. decorateGreeting testen

Füge eine zweite Testmethode nach der ersten hinzu.
Diese prüft, dass `decorateGreeting` den Eingabe-String auf jeder Seite mit `***` umschließt.

=== "Danach"

    ```groovy title="GreetingExtensionTest.groovy" linenums="1" hl_lines="18-25"
    package training.plugin

    import spock.lang.Specification

    /**
     * Tests für die Greeting-Extension-Funktionen
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

=== "Vorher"

    ```groovy title="GreetingExtensionTest.groovy" linenums="1" hl_lines="18"
    package training.plugin

    import spock.lang.Specification

    /**
     * Tests für die Greeting-Extension-Funktionen
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

## 4. Tests ausführen

```bash
make test
```

??? example "Test-Ausgabe"

    ```console
    BUILD SUCCESSFUL in 5s
    7 actionable tasks: 7 executed
    ```

    **Wo sind die Testergebnisse?** Gradle blendet die detaillierte Ausgabe aus, wenn alle Tests bestehen.
    „BUILD SUCCESSFUL" bedeutet, dass alles funktioniert hat.
    Schlägt ein Test fehl, werden detaillierte Fehlermeldungen angezeigt.

??? exercise "Einen Edge-Case-Test hinzufügen"

    Füge einen Test hinzu, der prüft, ob `reverseGreeting` einen leeren String verarbeitet.
    Was sollte `reverseGreeting('')` zurückgeben?
    Füge den Test hinzu, führe `make test` aus und überprüfe, ob er besteht.

    ??? solution "Lösung"

        Füge diese Testmethode zu `GreetingExtensionTest.groovy` hinzu:

        ```groovy
        def 'should handle empty string'() {
            given:
            def ext = new GreetingExtension()

            expect:
            ext.reverseGreeting('') == ''
        }
        ```

        Ein leerer String bleibt auch nach dem Umkehren ein leerer String.

---

## 5. Den Testbericht anzeigen

Gradle erstellt einen HTML-Testbericht mit detaillierten Ergebnissen für jeden Test.
Starte einen Webserver im Berichtsverzeichnis:

```bash
pushd build/reports/tests/test
python -m http.server
```

VS Code fordert dich auf, die Anwendung in deinem Browser zu öffnen.
Klicke dich zu deiner Test-Klasse durch, um die einzelnen Testergebnisse zu sehen:

![Testbericht, der zeigt, dass alle Tests bestanden haben](./img/test_report.png)

Der Bericht zeigt jede Testmethode und ob sie bestanden hat oder fehlgeschlagen ist.

Drücke ++ctrl+c++, um den Server zu stoppen, und kehre dann zum vorherigen Verzeichnis zurück:

```bash
popd
```

Wechsle zurück in das Hauptprojektverzeichnis:

```bash
cd ..
```

---

## Fazit

Du hast gelernt, dass:

- Spock-Tests eine lesbare `given:`/`expect:`-Struktur verwenden
- du mit `make test` Tests ausführst und den HTML-Bericht unter `build/reports/tests/test/` findest
- Tests das Verhalten prüfen und als Dokumentation dienen, die zeigt, wie Funktionen verwendet werden sollen

---

## Wie geht es weiter?

Bisher fügt dein Plugin benutzerdefinierte Funktionen hinzu, die Pipelines aufrufen können.
Plugins können auch auf Workflow-Ereignisse reagieren (eine abgeschlossene Aufgabe, eine veröffentlichte Datei, das Ende der Pipeline) – mithilfe von Trace-Observern.
Im nächsten Abschnitt baust du einen Observer, der abgeschlossene Aufgaben zählt und am Ende der Pipeline eine Zusammenfassung ausgibt.

[Weiter zu Teil 5 :material-arrow-right:](05_observers.md){ .md-button .md-button--primary }
