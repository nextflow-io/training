# Teil 3: Ausführungskonfiguration

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } KI-gestützte Übersetzung - [mehr erfahren & Verbesserungen vorschlagen](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Dieser Abschnitt wird erkunden, wie man die Konfiguration einer Nextflow Pipeline verwaltet, um ihr Verhalten anzupassen, sie an verschiedene Umgebungen anzupassen und die Ressourcennutzung zu optimieren - _ohne eine einzige Zeile des Workflow-Codes selbst zu ändern_.

Es gibt mehrere Möglichkeiten, dies zu tun, die kombiniert werden können und gemäß der [hier](https://www.nextflow.io/docs/latest/config.html) beschriebenen Rangfolge interpretiert werden.

In diesem Teil des Kurses werden wir dir den einfachsten und gebräuchlichsten Konfigurationsdatei-Mechanismus zeigen, die `nextflow.config`-Datei, die du bereits im Abschnitt über Container in Teil 2 kennengelernt hast.

Wir werden wesentliche Komponenten der Nextflow-Konfiguration durchgehen, wie process-Direktiven, executors, profiles und Parameter-Dateien.
Indem du lernst, diese Konfigurationsoptionen effektiv zu nutzen, kannst du die Flexibilität, Skalierbarkeit und Leistung von Nextflow Pipelines voll ausschöpfen.

Um diese Konfigurationselemente zu üben, werden wir eine frische Kopie des Workflows ausführen, den wir zuletzt am Ende von Teil 2 dieses Trainingskurses ausgeführt haben, umbenannt in `3-main.nf`.

Wenn du mit der Hello Pipeline nicht vertraut bist oder eine Auffrischung gebrauchen könntest, siehe [diese Info-Seite](../info/hello_pipeline.md).

---

## 1. Workflow-Eingabe-Parameter verwalten

??? example "Szenario"

    Du hast eine Pipeline heruntergeladen und möchtest sie wiederholt mit denselben Eingabedateien und Einstellungen ausführen, aber du möchtest nicht jedes Mal alle Parameter eintippen.
    Oder vielleicht richtest du die Pipeline für einen Kollegen ein, der sich nicht wohl mit Kommandozeilen-Argumenten fühlt.

Wir beginnen mit einem Aspekt der Konfiguration, der einfach eine Erweiterung dessen ist, womit wir bisher gearbeitet haben: die Verwaltung von Eingabe-Parametern.

Derzeit ist unser Workflow so eingerichtet, dass er mehrere Parameter-Werte über die Kommandozeile akzeptiert, die in einem `params`-Block im Workflow-Script selbst deklariert sind.
Einer hat einen Standardwert, der als Teil seiner Deklaration gesetzt ist.

Allerdings möchtest du vielleicht Standardwerte für alle setzen oder den bestehenden Standard überschreiben, ohne entweder Parameter auf der Kommandozeile angeben oder die Originaldatei ändern zu müssen.

Es gibt mehrere Möglichkeiten, das zu tun; wir werden dir drei grundlegende Wege zeigen, die sehr häufig verwendet werden.

### 1.1. Werte in `nextflow.config` einrichten

Das ist der einfachste Ansatz, obwohl er möglicherweise am wenigsten flexibel ist, da die Haupt-`nextflow.config`-Datei nichts ist, was du für jeden Lauf bearbeiten möchtest.
Aber es hat den Vorteil, die Anliegen der _Deklaration_ der Parameter im Workflow (was definitiv dorthin gehört) von der Bereitstellung von _Standardwerten_ zu trennen, die eher in einer Konfigurationsdatei zu Hause sind.

Lass uns das in zwei Schritten machen.

#### 1.1.1. Einen `params`-Block in der Konfigurationsdatei erstellen

Nimm die folgenden Code-Änderungen in der `nextflow.config`-Datei vor:

=== "Nachher"

    ```groovy title="nextflow.config" linenums="1" hl_lines="3-10"
    docker.enabled = true

    /*
    * Pipeline-Parameter
    */
    params {
        input = 'data/greetings.csv'
        batch = 'batch'
        character = 'turkey'
    }
    ```

=== "Vorher"

    ```groovy title="nextflow.config" linenums="1"
    docker.enabled = true
    ```

Beachte, dass wir nicht einfach den `params`-Block vom Workflow in die Konfigurationsdatei kopiert haben.
Für den `batch`-Parameter, der bereits einen Standardwert deklariert hatte, ist die Syntax etwas anders.
In der Workflow-Datei ist das eine typisierte Deklaration.
In der Konfiguration sind das Wertzuweisungen.

Technisch reicht das aus, um die noch in der Workflow-Datei angegebenen Standardwerte zu überschreiben.
Du könntest den Standardwert für `batch` ändern und den Workflow ausführen, um dich zu vergewissern, dass der in der Konfigurationsdatei gesetzte Wert den in der Workflow-Datei gesetzten überschreibt.

Aber im Sinne davon, die Konfiguration vollständig in die Konfigurationsdatei zu verschieben, lass uns diesen Standardwert ganz aus der Workflow-Datei entfernen.

#### 1.1.2. Den Standardwert für `batch` in der Workflow-Datei entfernen

Nimm die folgende Code-Änderung an der `3-main.nf` Workflow-Datei vor:

=== "Nachher"

    ```groovy title="3-main.nf" linenums="9" hl_lines="6"
    /*
    * Pipeline-Parameter
    */
    params {
        input: Path
        batch: String
        character: String
    }
    ```

=== "Vorher"

    ```groovy title="3-main.nf" linenums="9" hl_lines="6"
    /*
    * Pipeline-Parameter
    */
    params {
        input: Path
        batch: String = 'batch'
        character: String
    }
    ```

Jetzt setzt die Workflow-Datei selbst keine Standardwerte mehr für diese Parameter.

#### 1.1.3. Die Pipeline ausführen

Lass uns testen, dass es korrekt funktioniert, ohne irgendwelche Parameter in der Kommandozeile anzugeben.

```bash
nextflow run 3-main.nf
```

??? success "Befehlsausgabe"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `3-main.nf` [disturbed_einstein] DSL2 - revision: ede9037d02

    executor >  local (8)
    [f0/35723c] sayHello (2)       | 3 of 3 ✔
    [40/3efd1a] convertToUpper (3) | 3 of 3 ✔
    [17/e97d32] collectGreetings   | 1 of 1 ✔
    [98/c6b57b] cowpy              | 1 of 1 ✔
    ```

Das produziert immer noch dieselbe Ausgabe wie zuvor.

Die finale ASCII-Kunst-Ausgabe befindet sich im `results/3-main/`-Verzeichnis unter dem Namen `cowpy-COLLECTED-batch-output.txt`, wie zuvor.

??? abstract "Dateiinhalte"

    ```console title="results/3-main/cowpy-COLLECTED-batch-output.txt"
    _________
    / HOLà    \
    | HELLO   |
    \ BONJOUR /
    ---------
      \                                  ,+*^^*+___+++_
      \                           ,*^^^^              )
        \                       _+*                     ^**+_
        \                    +^       _ _++*+_+++_,         )
                  _+^^*+_    (     ,+*^ ^          \+_        )
                {       )  (    ,(    ,_+--+--,      ^)      ^\
                { (\@)    } f   ,(  ,+-^ __*_*_  ^^\_   ^\       )
              {:;-/    (_+*-+^^^^^+*+*<_ _++_)_    )    )      /
              ( /  (    (        ,___    ^*+_+* )   <    <      \
              U _/     )    *--<  ) ^\-----++__)   )    )       )
                (      )  _(^)^^))  )  )\^^^^^))^*+/    /       /
              (      /  (_))_^)) )  )  ))^^^^^))^^^)__/     +^^
            (     ,/    (^))^))  )  ) ))^^^^^^^))^^)       _)
              *+__+*       (_))^)  ) ) ))^^^^^^))^^^^^)____*^
              \             \_)^)_)) ))^^^^^^^^^^))^^^^)
              (_             ^\__^^^^^^^^^^^^))^^^^^^^)
                ^\___            ^\__^^^^^^))^^^^^^^^)\\
                      ^^^^^\uuu/^^\uuu/^^^^\^\^\^\^\^\^\^\
                        ___) >____) >___   ^\_\_\_\_\_\_\)
                        ^^^//\\_^^//\\_^       ^(\_\_\_\)
                          ^^^ ^^ ^^^ ^
    ```

Funktional hat diese Verschiebung nichts geändert, aber konzeptionell ist es etwas sauberer, die Standardwerte in der Konfigurationsdatei zu setzen.

### 1.2. Eine lauf-spezifische Konfigurationsdatei verwenden

??? example "Szenario"

    Du möchtest mit verschiedenen Einstellungen experimentieren, ohne deine Haupt-Konfigurationsdatei zu ändern.

Du kannst das tun, indem du eine neue `nextflow.config`-Datei in einem Unterverzeichnis erstellst, das du als Arbeitsverzeichnis für deine Experimente verwenden wirst.

#### 1.2.1. Das Arbeitsverzeichnis mit einer leeren Konfiguration erstellen

Lass uns beginnen, indem wir ein neues Verzeichnis erstellen und hineinwechseln:

```bash
mkdir -p tux-run
cd tux-run
```

Dann erstelle eine leere Konfigurationsdatei in diesem Verzeichnis:

```bash
touch nextflow.config
```

Das erzeugt eine leere Datei.

#### 1.2.2. Die experimentelle Konfiguration einrichten

Öffne jetzt die neue Datei und füge die Parameter hinzu, die du anpassen möchtest:

```groovy title="tux-run/nextflow.config" linenums="1"
params {
    input = '../data/greetings.csv'
    batch = 'experiment'
    character = 'tux'
}
```

Beachte, dass der Pfad zur Eingabedatei die Verzeichnisstruktur widerspiegeln muss.

#### 1.2.3. Die Pipeline ausführen

Wir können jetzt unsere Pipeline von innerhalb unseres neuen Arbeitsverzeichnisses ausführen.
Stelle sicher, den Pfad entsprechend anzupassen!

```bash
nextflow run ../3-main.nf
```

??? success "Befehlsausgabe"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `../3-main.nf` [trusting_escher] DSL2 - revision: 356df0818d

    executor >  local (8)
    [59/b66913] sayHello (2)       [100%] 3 of 3 ✔
    [ad/f06364] convertToUpper (3) [100%] 3 of 3 ✔
    [10/714895] collectGreetings   [100%] 1 of 1 ✔
    [88/3ece98] cowpy              [100%] 1 of 1 ✔
    ```

Das wird neue Verzeichnisse unter `tux-run/` erstellen, einschließlich `tux-run/work/` und `tux-run/results/`.

In diesem Lauf kombiniert Nextflow die `nextflow.config` in unserem aktuellen Verzeichnis mit der `nextflow.config` im Stammverzeichnis der Pipeline und überschreibt dadurch den Standard-Charakter (turkey) mit dem tux-Charakter.

Die finale Ausgabedatei sollte den tux-Charakter enthalten, der die Grüße sagt.

??? abstract "Dateiinhalte"

    ```console title="tux-run/results/3-main/cowpy-COLLECTED-experiment-output.txt"
    _________
    / HELLO   \
    | BONJOUR |
    \ HOLà    /
    ---------
      \
        \
            .--.
          |o_o |
          |:_/ |
          //   \ \
        (|     | )
        /'\_   _/`\
        \___)=(___/

    ```

Das war's; jetzt hast du einen Raum zum Experimentieren, ohne deine 'normale' Konfiguration zu ändern.

!!! warning "Warnung"

    Stelle sicher, dass du zum vorherigen Verzeichnis zurückwechselst, bevor du zum nächsten Abschnitt übergehst!

    ```bash
    cd ..
    ```

Jetzt schauen wir uns eine weitere nützliche Möglichkeit an, Parameter-Werte zu setzen.

### 1.3. Eine Parameter-Datei verwenden

??? example "Szenario"

    Du musst exakte Ausführungs-Parameter mit einem Mitarbeiter teilen oder sie für eine Publikation aufzeichnen.

Der Unterverzeichnis-Ansatz funktioniert großartig zum Experimentieren, aber er erfordert etwas Einrichtung und verlangt, dass du Pfade entsprechend anpasst.
Es gibt einen einfacheren Ansatz für den Fall, dass du deine Pipeline mit einem bestimmten Satz von Werten ausführen oder jemand anderem ermöglichen möchtest, dies mit minimalem Aufwand zu tun.

Nextflow erlaubt uns, Parameter über eine Parameter-Datei im YAML- oder JSON-Format anzugeben, was es sehr bequem macht, alternative Sätze von Standardwerten sowie lauf-spezifische Parameter-Werte zu verwalten und zu verteilen.

#### 1.3.1. Die Beispiel-Parameter-Datei untersuchen

Um das zu demonstrieren, stellen wir eine Beispiel-Parameter-Datei im aktuellen Verzeichnis bereit, genannt `test-params.yaml`:

```yaml title="test-params.yaml" linenums="1"
input: "data/greetings.csv"
batch: "yaml"
character: "stegosaurus"
```

Diese Parameter-Datei enthält ein Schlüssel-Wert-Paar für jede der Eingaben, die wir angeben möchten.
Beachte die Verwendung von Doppelpunkten (`:`) anstelle von Gleichheitszeichen (`=`), wenn du die Syntax mit der Konfigurationsdatei vergleichst.
Die config-Datei ist in Groovy geschrieben, während die Parameter-Datei in YAML geschrieben ist.

!!! info "Hinweis"

    Wir stellen auch eine JSON-Version der Parameter-Datei als Beispiel bereit, aber wir werden hier nicht damit ausführen.
    Fühl dich frei, diese auf eigene Faust auszuprobieren.

#### 1.3.2. Die Pipeline ausführen

Um den Workflow mit dieser Parameter-Datei auszuführen, füge einfach `-params-file <dateiname>` zum Basis-Befehl hinzu.

```bash
nextflow run 3-main.nf -params-file test-params.yaml
```

??? success "Befehlsausgabe"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `3-main.nf` [disturbed_sammet] DSL2 - revision: ede9037d02

    executor >  local (8)
    [f0/35723c] sayHello (2)       | 3 of 3 ✔
    [40/3efd1a] convertToUpper (3) | 3 of 3 ✔
    [17/e97d32] collectGreetings   | 1 of 1 ✔
    [98/c6b57b] cowpy              | 1 of 1 ✔
    ```

Die finale Ausgabedatei sollte den stegosaurus-Charakter enthalten, der die Grüße sagt.

??? abstract "Dateiinhalte"

    ```console title="results/3-main/cowpy-COLLECTED-yaml-output.txt"
    _________
    / HELLO   \
    | HOLà    |
    \ BONJOUR /
    ---------
    \                             .       .
    \                           / `.   .' "
      \                  .---.  <    > <    >  .---.
      \                 |    \  \ - ~ ~ - /  /    |
            _____          ..-~             ~-..-~
            |     |   \~~~\.'                    `./~~~/
          ---------   \__/                        \__/
          .'  O    \     /               /       \  "
        (_____,    `._.'               |         }  \/~~~/
          `----.          /       }     |        /    \__/
                `-.      |       /      |       /      `. ,~~|
                    ~-.__|      /_ - ~ ^|      /- _      `..-'
                        |     /        |     /     ~-.     `-. _  _  _
                        |_____|        |_____|         ~ - . _ _ _ _ _>
    ```

Die Verwendung einer Parameter-Datei mag übertrieben erscheinen, wenn du nur ein paar Parameter angeben musst, aber einige Pipelines erwarten Dutzende von Parametern.
In diesen Fällen erlaubt uns die Verwendung einer Parameter-Datei, Parameter-Werte zur Laufzeit bereitzustellen, ohne massive Kommandozeilen eintippen zu müssen und ohne das Workflow-Script zu ändern.

Es macht es auch einfacher, Sätze von Parametern an Mitarbeiter zu verteilen oder als Zusatzinformationen für eine Publikation zum Beispiel.
Das macht deine Arbeit für andere reproduzierbarer.

### Zusammenfassung

Du weißt, wie du wichtige Konfigurationsoptionen für die Verwaltung von Workflow-Eingaben nutzen kannst.

### Was kommt als Nächstes?

Lerne, wie du verwalten kannst, wo und wie deine Workflow-Ausgaben veröffentlicht werden.

---

## 2. Workflow-Ausgaben verwalten

??? example "Szenario"

    Deine Pipeline veröffentlicht Ausgaben in ein fest codiertes Verzeichnis, aber du möchtest die Ergebnisse nach Projekt- oder Experimentnamen organisieren, ohne jedes Mal den Workflow-Code bearbeiten zu müssen.

Der Workflow, den wir geerbt haben, verwendet Pfade für Workflow-Level-Ausgabe-Deklarationen, was nicht besonders flexibel ist und viel Wiederholung beinhaltet.

Schauen wir uns ein paar gängige Möglichkeiten an, wie du das konfigurieren könntest, um flexibler zu sein.

### 2.1. Den `outputDir`-Verzeichnisnamen anpassen

Jede Version des Workflows, die wir bisher ausgeführt haben, hat ihre Ausgaben in ein anderes Unterverzeichnis veröffentlicht, das in den Ausgabe-Definitionen fest codiert ist.

Lass uns das ändern, um einen benutzer-konfigurierbaren Parameter zu verwenden.
Wir könnten einen ganz neuen Parameter dafür erstellen, aber lass uns den `batch`-Parameter verwenden, da er direkt verfügbar ist.

#### 2.1.1. Einen Wert für `outputDir` in der Konfigurationsdatei setzen

Der Pfad, den Nextflow zum Veröffentlichen von Ausgaben verwendet, wird durch die `outputDir`-Option gesteuert.
Um den Pfad für alle Ausgaben zu ändern, kannst du in der `nextflow.config`-Konfigurationsdatei einen Wert für diese Option setzen.

Füge den folgenden Code zur `nextflow.config`-Datei hinzu:

=== "Nachher"

    ```groovy title="nextflow.config" linenums="9" hl_lines="10-13"
    /*
    * Pipeline-Parameter
    */
    params {
        input = 'data/greetings.csv'
        batch = 'batch'
        character = 'turkey'
    }

    /*
    * Output settings
    */
    outputDir = "results/${params.batch}"
    ```

=== "Vorher"

    ```groovy title="nextflow.config" linenums="9"
    /*
    * Pipeline-Parameter
    */
    params {
        input = 'data/greetings.csv'
        batch = 'batch'
        character = 'turkey'
    }
    ```

Das wird den eingebauten Standardpfad `results/` durch `results/` plus dem Wert des `batch`-Parameters als Unterverzeichnis ersetzen.
Du könntest auch den `results`-Teil ändern, wenn du möchtest.

Für eine temporäre Änderung könntest du diese Option von der Kommandozeile aus setzen, indem du den `-output-dir`-Parameter in deinem Befehl verwendest (aber dann könntest du nicht den `batch`-Parameter-Wert verwenden).

#### 2.1.2. Den wiederholten Teil des fest codierten Pfades entfernen

Wir haben immer noch ein Unterverzeichnis fest codiert in den Ausgabe-Optionen, also lass uns das jetzt loswerden.

Nimm die folgenden Code-Änderungen in der Workflow-Datei vor:

=== "Nachher"

    ```groovy title="3-main.nf" linenums="42" hl_lines="3 7 11 15 19"
    output {
        first_output {
            path 'intermediates'
            mode 'copy'
        }
        uppercased {
            path 'intermediates'
            mode 'copy'
        }
        collected {
            path 'intermediates'
            mode 'copy'
        }
        batch_report {
            path ''
            mode 'copy'
        }
        cowpy_art {
            path ''
            mode 'copy'
        }
    }
    ```

=== "Vorher"

    ```groovy title="3-main.nf" linenums="42" hl_lines="3 7 11 15 19"
    output {
        first_output {
            path '3-main/intermediates'
            mode 'copy'
        }
        uppercased {
            path '3-main/intermediates'
            mode 'copy'
        }
        collected {
            path '3-main/intermediates'
            mode 'copy'
        }
        batch_report {
            path '3-main'
            mode 'copy'
        }
        cowpy_art {
            path '3-main'
            mode 'copy'
        }
    }
    ```

Wir hätten auch einfach `${params.batch}` zu jedem Pfad hinzufügen können, anstatt den `outputDir`-Standard zu ändern, aber das ist prägnanter.

#### 2.1.3. Die Pipeline ausführen

Lass uns testen, dass es korrekt funktioniert, indem wir den Batch-Namen auf `outdir` von der Kommandozeile setzen.

```bash
nextflow run 3-main.nf --batch outdir
```

??? success "Befehlsausgabe"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `3-main.nf` [disturbed_einstein] DSL2 - revision: ede9037d02

    executor >  local (8)
    [f0/35723c] sayHello (2)       | 3 of 3 ✔
    [40/3efd1a] convertToUpper (3) | 3 of 3 ✔
    [17/e97d32] collectGreetings   | 1 of 1 ✔
    [98/c6b57b] cowpy              | 1 of 1 ✔
    ```

Das produziert immer noch dieselbe Ausgabe wie zuvor, außer dass wir unsere Ausgaben diesmal unter `results/outdir/` finden.

??? abstract "Verzeichnisinhalte"

    ```console
    results/outdir/
    ├── cowpy-COLLECTED-outdir-output.txt
    ├── intermediates
    │   ├── Bonjour-output.txt
    │   ├── COLLECTED-outdir-output.txt
    │   ├── Hello-output.txt
    │   ├── Holà-output.txt
    │   ├── UPPER-Bonjour-output.txt
    │   ├── UPPER-Hello-output.txt
    │   └── UPPER-Holà-output.txt
    └── outdir-report.txt
    ```

Du kannst diesen Ansatz mit benutzerdefinierten Pfad-Definitionen kombinieren, um jede gewünschte Verzeichnishierarchie zu konstruieren.

### 2.2. Ausgaben nach process organisieren

Eine beliebte Möglichkeit, Ausgaben weiter zu organisieren, ist es nach process zu tun, _d.h._ Unterverzeichnisse für jeden in der Pipeline ausgeführten process zu erstellen.

#### 2.2.1. Die Ausgabe-Pfade durch einen Verweis auf process-Namen ersetzen

Alles, was du tun musst, ist den Namen des process als `<task>.name` in der Ausgabe-Pfad-Deklaration zu referenzieren.

Nimm die folgenden Änderungen in der Workflow-Datei vor:

=== "Nachher"

    ```groovy title="3-main.nf" linenums="42" hl_lines="3 7 11 15 19"
    output {
        first_output {
            path { sayHello.name }
            mode 'copy'
        }
        uppercased {
            path { convertToUpper.name }
            mode 'copy'
        }
        collected {
            path { collectGreetings.name }
            mode 'copy'
        }
        batch_report {
            path { collectGreetings.name }
            mode 'copy'
        }
        cowpy_art {
            path { cowpy.name }
            mode 'copy'
        }
    }
    ```

=== "Vorher"

    ```groovy title="3-main.nf" linenums="42" hl_lines="3 7 11 15 19"
    output {
        first_output {
            path 'intermediates'
            mode 'copy'
        }
        uppercased {
            path 'intermediates'
            mode 'copy'
        }
        collected {
            path 'intermediates'
            mode 'copy'
        }
        batch_report {
            path ''
            mode 'copy'
        }
        cowpy_art {
            path ''
            mode 'copy'
        }
    }
    ```

Das entfernt die verbleibenden fest codierten Elemente aus der Ausgabe-Pfad-Konfiguration.

#### 2.2.2. Die Pipeline ausführen

Lass uns testen, dass es korrekt funktioniert, indem wir den Batch-Namen auf `pnames` von der Kommandozeile setzen.

```bash
nextflow run 3-main.nf --batch pnames
```

??? success "Befehlsausgabe"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `3-main.nf` [jovial_mcclintock] DSL2 - revision: ede9037d02

    executor >  local (8)
    [f0/35723c] sayHello (2)       | 3 of 3 ✔
    [40/3efd1a] convertToUpper (3) | 3 of 3 ✔
    [17/e97d32] collectGreetings   | 1 of 1 ✔
    [98/c6b57b] cowpy              | 1 of 1 ✔
    ```

Das produziert immer noch dieselbe Ausgabe wie zuvor, außer dass wir unsere Ausgaben diesmal unter `results/pnames/` finden, und sie sind nach process gruppiert.

??? abstract "Verzeichnisinhalte"

    ```console
    results/pnames/
    ├── collectGreetings
    │   ├── COLLECTED-pnames-output.txt
    │   └── pnames-report.txt
    ├── convertToUpper
    │   ├── UPPER-Bonjour-output.txt
    │   ├── UPPER-Hello-output.txt
    │   └── UPPER-Holà-output.txt
    ├── cowpy
    │   └── cowpy-COLLECTED-pnames-output.txt
    └── sayHello
        ├── Bonjour-output.txt
        ├── Hello-output.txt
        └── Holà-output.txt
    ```

Beachte, dass wir hier die Unterscheidung zwischen `intermediates` versus finalen Ausgaben auf der obersten Ebene gelöscht haben.
Du könntest natürlich diese Ansätze mischen und kombinieren, zum Beispiel indem du den Pfad der ersten Ausgabe als `intermediates/${sayHello.name}` setzt.

### 2.3. Den publish-Modus auf Workflow-Ebene setzen

Schließlich können wir im Sinne der Reduzierung der Menge an sich wiederholendem Code die pro-Ausgabe `mode`-Deklarationen durch eine einzelne Zeile in der Konfiguration ersetzen.

#### 2.3.1. `workflow.output.mode` zur Konfigurationsdatei hinzufügen

Füge den folgenden Code zur `nextflow.config`-Datei hinzu:

=== "Nachher"

    ```groovy title="nextflow.config" linenums="2" hl_lines="5"
    /*
    * Output settings
    */
    outputDir = "results/${params.batch}"
    workflow.output.mode = 'copy'
    ```

=== "Vorher"

    ```groovy title="nextflow.config" linenums="12"
    /*
    * Output settings
    */
    outputDir = "results/${params.batch}"
    ```

Genau wie bei der `outputDir`-Option würde es ausreichen, `workflow.output.mode` einen Wert in der Konfigurationsdatei zu geben, um zu überschreiben, was in der Workflow-Datei gesetzt ist, aber lass uns den unnötigen Code trotzdem entfernen.

#### 2.3.2. Den output-Modus aus der Workflow-Datei entfernen

Nimm die folgenden Änderungen in der Workflow-Datei vor:

=== "Nachher"

    ```groovy title="3-main.nf" linenums="42"
    output {
        first_output {
            path { sayHello.name }
        }
        uppercased {
            path { convertToUpper.name }
        }
        collected {
            path { collectGreetings.name }
        }
        batch_report {
            path { collectGreetings.name }
        }
        cowpy_art {
            path { cowpy.name }
        }
    }
    ```

=== "Vorher"

    ```groovy title="3-main.nf" linenums="42" hl_lines="3 7 11 15 19"
    output {
        first_output {
            path { sayHello.name }
            mode 'copy'
        }
        uppercased {
            path { convertToUpper.name }
            mode 'copy'
        }
        collected {
            path { collectGreetings.name }
            mode 'copy'
        }
        batch_report {
            path { collectGreetings.name }
            mode 'copy'
        }
        cowpy_art {
            path { cowpy.name }
            mode 'copy'
        }
    }
    ```

Das ist prägnanter, nicht wahr?

#### 2.3.3. Die Pipeline ausführen

Lass uns testen, dass es korrekt funktioniert, indem wir den Batch-Namen auf `outmode` von der Kommandozeile setzen.

```bash
nextflow run 3-main.nf --batch outmode
```

??? success "Befehlsausgabe"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `3-main.nf` [rowdy_sagan] DSL2 - revision: ede9037d02

    executor >  local (8)
    [f0/35723c] sayHello (2)       | 3 of 3 ✔
    [40/3efd1a] convertToUpper (3) | 3 of 3 ✔
    [17/e97d32] collectGreetings   | 1 of 1 ✔
    [98/c6b57b] cowpy              | 1 of 1 ✔
    ```

Das produziert immer noch dieselbe Ausgabe wie zuvor, außer dass wir unsere Ausgaben diesmal unter `results/outmode/` finden.
Es sind immer noch alles echte Kopien, keine Symlinks.

??? abstract "Verzeichnisinhalte"

    ```console
    results/outmode/
    ├── collectGreetings
    │   ├── COLLECTED-outmode-output.txt
    │   └── outmode-report.txt
    ├── convertToUpper
    │   ├── UPPER-Bonjour-output.txt
    │   ├── UPPER-Hello-output.txt
    │   └── UPPER-Holà-output.txt
    ├── cowpy
    │   └── cowpy-COLLECTED-outmode-output.txt
    └── sayHello
        ├── Bonjour-output.txt
        ├── Hello-output.txt
        └── Holà-output.txt
    ```

Der Hauptgrund, warum du vielleicht immer noch die pro-Ausgabe-Art des Setzens des mode verwenden möchtest, ist, wenn du innerhalb desselben Workflows mischen und kombinieren möchtest, _d.h._ einige Ausgaben kopiert und einige per Symlink verlinkt haben möchtest.

Es gibt viele andere Optionen, die du auf diese Weise anpassen kannst, aber hoffentlich gibt dir das einen Eindruck von der Bandbreite der Optionen und wie du sie effektiv nutzen kannst, um deinen Präferenzen zu entsprechen.

### Zusammenfassung

Du weißt, wie du die Benennung und Struktur der Verzeichnisse kontrollieren kannst, in denen deine Ausgaben veröffentlicht werden, sowie den Workflow-Ausgabe-Veröffentlichungsmodus.

### Was kommt als Nächstes?

Lerne, wie du deine Workflow-Konfiguration an deine Rechenumgebung anpassen kannst, beginnend mit der Software-Paketierungstechnologie.

---

## 3. Eine Software-Paketierungstechnologie auswählen

Bisher haben wir Konfigurationselemente betrachtet, die steuern, wie Eingaben hineingehen und wo Ausgaben herauskommen. Jetzt ist es an der Zeit, uns spezieller darauf zu konzentrieren, deine Workflow-Konfiguration an deine Rechenumgebung anzupassen.

Der erste Schritt auf diesem Weg ist anzugeben, woher die Software-Pakete kommen werden, die in jedem Schritt ausgeführt werden.
Sind sie bereits in der lokalen Rechenumgebung installiert?
Müssen wir Images abrufen und sie über ein Container-System ausführen?
Oder müssen wir Conda-Pakete abrufen und eine lokale Conda-Umgebung erstellen?

Im allerersten Teil dieses Trainingskurses (Teile 1-4) haben wir in unserem Workflow einfach lokal installierte Software verwendet.
Dann haben wir in Teil 5 Docker-Container und die `nextflow.config`-Datei eingeführt, die wir verwendet haben, um die Verwendung von Docker-Containern zu aktivieren.

Jetzt schauen wir uns an, wie wir eine alternative Software-Paketierungsoption über die `nextflow.config`-Datei konfigurieren können.

### 3.1. Docker deaktivieren und Conda in der config-Datei aktivieren

??? example "Szenario"

    Du verschiebst deine Pipeline auf einen HPC-Cluster, wo Docker aus Sicherheitsgründen nicht erlaubt ist.
    Der Cluster unterstützt Singularity und Conda, also musst du deine Konfiguration entsprechend umstellen.

Nextflow unterstützt mehrere Container-Technologien einschließlich Singularity (das auf HPC weiter verbreitet ist), sowie Software-Paketmanager wie Conda.

Wir können unsere Konfigurationsdatei ändern, um Conda anstelle von Docker zu verwenden.
Dazu setzen wir den Wert von `docker.enabled` auf `false` und fügen eine Direktive hinzu, die die Verwendung von Conda aktiviert:

=== "Nachher"

    ```groovy title="nextflow.config" linenums="1" hl_lines="1-2"
    docker.enabled = false
    conda.enabled = true
    ```

=== "Vorher"

    ```groovy title="nextflow.config" linenums="1" hl_lines="1"
    docker.enabled = true
    ```

Das wird Nextflow erlauben, Conda-Umgebungen für processes zu erstellen und zu nutzen, die Conda-Pakete angegeben haben.
Das bedeutet, dass wir jetzt eines davon zu unserem `cowpy` process hinzufügen müssen!

### 3.2. Ein Conda-Paket in der process-Definition angeben

Wir haben bereits die URI für ein Conda-Paket abgerufen, das das `cowpy`-Tool enthält: `conda-forge::cowpy==1.1.5`

Jetzt fügen wir die URI zur `cowpy` process-Definition hinzu, indem wir die `conda`-Direktive verwenden:

=== "Nachher"

    ```groovy title="modules/cowpy.nf" linenums="4" hl_lines="4"
    process cowpy {

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
        conda 'conda-forge::cowpy==1.1.5'

        input:
    ```

=== "Vorher"

    ```groovy title="modules/cowpy.nf" linenums="4"
    process cowpy {

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'

        input:
    ```

Um das klarzustellen, wir _ersetzen_ nicht die `docker`-Direktive, wir _fügen_ eine alternative Option hinzu.

!!! tip "Tipp"

    Es gibt einige verschiedene Möglichkeiten, die URI für ein bestimmtes Conda-Paket zu bekommen.
    Wir empfehlen, die [Seqera Containers](https://seqera.io/containers/) Suchabfrage zu verwenden, die dir eine URI gibt, die du kopieren und einfügen kannst, selbst wenn du nicht planst, einen Container daraus zu erstellen.

### 3.3. Den Workflow ausführen, um zu verifizieren, dass er Conda verwenden kann

Lass es uns ausprobieren.

```bash
nextflow run 3-main.nf --batch conda
```

??? success "Befehlsausgabe"

    ```console title="Ausgabe"
    N E X T F L O W   ~  version 25.10.2

    Launching `3-main.nf` [trusting_lovelace] DSL2 - revision: 028a841db1

    executor >  local (8)
    [ee/4ca1f2] sayHello (3)       | 3 of 3 ✔
    [20/2596a7] convertToUpper (1) | 3 of 3 ✔
    [b3/e15de5] collectGreetings   | 1 of 1 ✔
    [c5/af5f88] cowpy              | 1 of 1 ✔
    ```

Das sollte ohne Probleme funktionieren und dieselben Ausgaben wie zuvor unter `results/conda` produzieren.

Hinter den Kulissen hat Nextflow die Conda-Pakete abgerufen und die Umgebung erstellt, was normalerweise etwas Arbeit erfordert; also ist es schön, dass wir nichts davon selbst machen müssen!

!!! info "Hinweis"

    Das läuft schnell, weil das `cowpy`-Paket ziemlich klein ist, aber wenn du mit großen Paketen arbeitest, kann es beim ersten Mal etwas länger dauern als üblich, und du könntest sehen, dass die Konsolenausgabe für eine Minute oder so 'stecken bleibt', bevor sie abgeschlossen ist.
    Das ist normal und liegt an der zusätzlichen Arbeit, die Nextflow beim ersten Mal macht, wenn du ein neues Paket verwendest.

Aus unserer Sicht sieht es so aus, als ob es genau so funktioniert wie mit Docker, obwohl die Mechanik im Hintergrund etwas anders ist.

Das bedeutet, dass wir bereit sind, bei Bedarf mit Conda-Umgebungen auszuführen.

??? info "Docker und Conda mischen und kombinieren"

    Da diese Direktiven pro process zugewiesen werden, ist es möglich, zu 'mischen und kombinieren', _d.h._ einige der processes in deinem Workflow so zu konfigurieren, dass sie mit Docker laufen und andere mit Conda, zum Beispiel, wenn die von dir verwendete Recheninfrastruktur beides unterstützt.
    In diesem Fall würdest du sowohl Docker als auch Conda in deiner Konfigurationsdatei aktivieren.
    Wenn beide für einen bestimmten process verfügbar sind, wird Nextflow Container priorisieren.

    Und wie bereits erwähnt, unterstützt Nextflow mehrere andere Software-Paketierungs- und Container-Technologien, sodass du nicht auf nur diese beiden beschränkt bist.

### Zusammenfassung

Du weißt, wie du konfigurierst, welches Software-Paket jeder process verwenden soll, und wie du zwischen Technologien wechselst.

### Was kommt als Nächstes?

Lerne, wie du die Ausführungsplattform änderst, die Nextflow verwendet, um die eigentliche Arbeit zu erledigen.

---

## 4. Eine Ausführungsplattform auswählen

??? example "Szenario"

    Du hast deine Pipeline auf deinem Laptop entwickelt und getestet, aber jetzt musst du sie auf Tausenden von Proben ausführen.
    Deine Institution hat einen HPC-Cluster mit einem Slurm-Scheduler, den du stattdessen verwenden möchtest.

Bis jetzt haben wir unsere Pipeline mit dem local executor ausgeführt.
Dieser führt jede Aufgabe auf der Maschine aus, auf der Nextflow läuft.
Wenn Nextflow startet, schaut es sich die verfügbaren CPUs und den Speicher an.
Wenn die Ressourcen der zur Ausführung bereiten Aufgaben die verfügbaren Ressourcen überschreiten, wird Nextflow die letzten Aufgaben von der Ausführung zurückhalten, bis eine oder mehrere der früheren Aufgaben beendet sind und die notwendigen Ressourcen freigeben.

Der local executor ist bequem und effizient, aber er ist auf diese einzelne Maschine beschränkt. Bei sehr großen Workloads könntest du feststellen, dass deine lokale Maschine ein Engpass ist, entweder weil du eine einzelne Aufgabe hast, die mehr Ressourcen erfordert als du verfügbar hast, oder weil du so viele Aufgaben hast, dass das Warten darauf, dass eine einzelne Maschine sie ausführt, zu lange dauern würde.

Nextflow unterstützt [viele verschiedene Ausführungs-Backends](https://www.nextflow.io/docs/latest/executor.html), einschließlich HPC-Scheduler (Slurm, LSF, SGE, PBS, Moab, OAR, Bridge, HTCondor und andere) sowie Cloud-Ausführungs-Backends (AWS Batch, Google Cloud Batch, Azure Batch, Kubernetes und mehr).

### 4.1. Ein anderes Backend anvisieren

Die Wahl des executor wird durch eine process-Direktive namens `executor` gesetzt.
Standardmäßig ist sie auf `local` gesetzt, also ist die folgende Konfiguration impliziert:

```groovy title="Eingebaute Konfiguration"
process {
    executor = 'local'
}
```

Um den executor auf ein anderes Backend zu setzen, würdest du einfach den gewünschten executor mit ähnlicher Syntax angeben, wie oben für Ressourcenzuweisungen beschrieben (siehe [Dokumentation](https://www.nextflow.io/docs/latest/executor.html) für alle Optionen).

```groovy title="nextflow.config"
process {
    executor = 'slurm'
}
```

!!! warning "Warnung"

    Wir können das in der Trainingsumgebung nicht wirklich testen, weil sie nicht eingerichtet ist, um sich mit einem HPC zu verbinden.

### 4.2. Mit Backend-spezifischer Syntax für Ausführungsparameter umgehen

Die meisten Hochleistungs-Rechenplattformen erlauben (und erfordern manchmal), dass du bestimmte Parameter wie Ressourcenzuweisungsanfragen und -beschränkungen (z.B. Anzahl der CPUs und Speicher) und den Namen der zu verwendenden Job-Queue angibst.

Leider verwendet jedes dieser Systeme unterschiedliche Technologien, Syntaxen und Konfigurationen, um zu definieren, wie ein Job definiert und an den entsprechenden Scheduler übermittelt werden soll.

??? abstract "Beispiele"

    Zum Beispiel muss derselbe Job, der 8 CPUs und 4GB RAM benötigt und auf der Queue "my-science-work" ausgeführt werden soll, je nach Backend auf unterschiedliche Weise ausgedrückt werden.

    ```bash title="Config für SLURM / übermitteln mit sbatch"
    #SBATCH -o /path/to/my/task/directory/my-task-1.log
    #SBATCH --no-requeue
    #SBATCH -c 8
    #SBATCH --mem 4096M
    #SBATCH -p my-science-work
    ```

    ```bash title="Config für PBS / übermitteln mit qsub"
    #PBS -o /path/to/my/task/directory/my-task-1.log
    #PBS -j oe
    #PBS -q my-science-work
    #PBS -l nodes=1:ppn=5
    #PBS -l mem=4gb
    ```

    ```bash title="Config für SGE / übermitteln mit qsub"
    #$ -o /path/to/my/task/directory/my-task-1.log
    #$ -j y
    #$ -terse
    #$ -notify
    #$ -q my-science-work
    #$ -l slots=5
    #$ -l h_rss=4096M,mem_free=4096M
    ```

Glücklicherweise vereinfacht Nextflow all das.
Es bietet eine standardisierte Syntax, sodass du die relevanten Eigenschaften wie `cpus`, `memory` und `queue` (siehe Dokumentation für andere Eigenschaften) nur einmal angeben kannst.
Dann wird Nextflow zur Laufzeit diese Einstellungen verwenden, um die entsprechenden Backend-spezifischen Scripts basierend auf der executor-Einstellung zu generieren.

Wir werden diese standardisierte Syntax im nächsten Abschnitt behandeln.

### Zusammenfassung

Du weißt jetzt, wie du den executor änderst, um verschiedene Arten von Recheninfrastruktur zu verwenden.

### Was kommt als Nächstes?

Lerne, wie du Rechenressourcenzuweisungen und -beschränkungen in Nextflow ausdrückst.

---

## 5. Rechenressourcenzuweisungen steuern

??? example "Szenario"

    Deine Pipeline schlägt auf dem Cluster immer fehl, weil Aufgaben wegen Überschreitung von Speicherlimits beendet werden.
    Oder vielleicht werden dir Ressourcen berechnet, die du nicht nutzt, und du möchtest die Kosten optimieren.

Die meisten Hochleistungs-Rechenplattformen erlauben (und erfordern manchmal), dass du bestimmte Ressourcenzuweisungsparameter wie Anzahl der CPUs und Speicher angibst.

Standardmäßig wird Nextflow eine einzelne CPU und 2GB Speicher für jeden process verwenden.
Die entsprechenden process-Direktiven heißen `cpus` und `memory`, also ist die folgende Konfiguration impliziert:

```groovy title="Eingebaute Konfiguration" linenums="1"
process {
    cpus = 1
    memory = 2.GB
}
```

Du kannst diese Werte ändern, entweder für alle processes oder für bestimmte benannte processes, indem du zusätzliche process-Direktiven in deiner Konfigurationsdatei verwendest.
Nextflow wird sie in die entsprechenden Anweisungen für den gewählten executor übersetzen.

Aber woher weißt du, welche Werte du verwenden sollst?

### 5.1. Den Workflow ausführen, um einen Ressourcennutzungsbericht zu generieren

??? example "Szenario"

    Du weißt nicht, wie viel Speicher oder CPU deine processes brauchen und möchtest vermeiden, Ressourcen zu verschwenden oder Jobs beenden zu lassen.

Wenn du nicht im Voraus weißt, wie viel CPU und Speicher deine processes wahrscheinlich brauchen, kannst du etwas Ressourcen-Profiling machen, was bedeutet, dass du den Workflow mit einigen Standardzuweisungen ausführst, aufzeichnest, wie viel jeder process verwendet hat, und von dort aus schätzt, wie du die Basiszuweisungen anpassen kannst.

Praktischerweise enthält Nextflow eingebaute Tools dafür und wird auf Anfrage gerne einen Bericht für dich generieren.

Um das zu tun, füge `-with-report <dateiname>.html` zu deiner Kommandozeile hinzu.

```bash
nextflow run 3-main.nf -with-report report-config-1.html
```

Der Bericht ist eine HTML-Datei, die du herunterladen und in deinem Browser öffnen kannst. Du kannst auch mit der rechten Maustaste darauf im Datei-Explorer auf der linken Seite klicken und auf `Show preview` klicken, um sie in der Trainingsumgebung anzuzeigen.

Nimm dir ein paar Minuten Zeit, um den Bericht durchzusehen und zu schauen, ob du einige Möglichkeiten zur Anpassung der Ressourcen identifizieren kannst.
Stelle sicher, dass du auf die Tabs klickst, die die Nutzungsergebnisse als Prozentsatz dessen zeigen, was zugewiesen wurde.
Es gibt [Dokumentation](https://www.nextflow.io/docs/latest/reports.html), die alle verfügbaren Funktionen beschreibt.

### 5.2. Ressourcenzuweisungen für alle processes setzen

Das Profiling zeigt, dass die processes in unserem Trainings-Workflow sehr leichtgewichtig sind, also lass uns die Standard-Speicherzuweisung auf 1GB pro process reduzieren.

Füge das Folgende zu deiner `nextflow.config`-Datei hinzu, vor dem Abschnitt mit den Pipeline-Parametern:

```groovy title="nextflow.config" linenums="4"
/*
* Process settings
*/
process {
    memory = 1.GB
}
```

Das wird helfen, die Menge an Rechenleistung zu reduzieren, die wir verbrauchen.

### 5.3. Ressourcenzuweisungen für einen bestimmten process setzen

Gleichzeitig werden wir so tun, als ob der `cowpy` process mehr Ressourcen benötigt als die anderen, nur um zu demonstrieren, wie man Zuweisungen für einen einzelnen process anpasst.

=== "Nachher"

    ```groovy title="nextflow.config" linenums="4" hl_lines="6-9"
    /*
    * Process settings
    */
    process {
        memory = 1.GB
        withName: 'cowpy' {
            memory = 2.GB
            cpus = 2
        }
    }
    ```

=== "Vorher"

    ```groovy title="nextflow.config" linenums="4"
    /*
    * Process settings
    */
    process {
        memory = 1.GB
    }
    ```

Mit dieser Konfiguration werden alle processes 1GB Speicher und eine einzelne CPU (der implizierte Standard) anfordern, außer dem `cowpy` process, der 2GB und 2 CPUs anfordern wird.

!!! info "Hinweis"

    Wenn du eine Maschine mit wenigen CPUs hast und du eine hohe Anzahl pro process zuweist, könntest du sehen, dass process-Aufrufe hintereinander in die Warteschlange gestellt werden.
    Das liegt daran, dass Nextflow sicherstellt, dass wir nicht mehr CPUs anfordern als verfügbar sind.

### 5.4. Den Workflow mit der aktualisierten Konfiguration ausführen

Lass uns das ausprobieren, indem wir einen anderen Dateinamen für den Profiling-Bericht angeben, damit wir die Leistung vor und nach den Konfigurationsänderungen vergleichen können.

```bash
nextflow run 3-main.nf -with-report report-config-2.html
```

Du wirst wahrscheinlich keinen echten Unterschied bemerken, da dies ein so kleiner Workload ist, aber das ist der Ansatz, den du verwenden würdest, um die Leistung und Ressourcenanforderungen eines realen Workflows zu analysieren.

Es ist sehr nützlich, wenn deine processes unterschiedliche Ressourcenanforderungen haben. Es befähigt dich, die Ressourcenzuweisungen, die du für jeden process einrichtest, basierend auf tatsächlichen Daten richtig zu dimensionieren, nicht auf Vermutungen.

!!! tip "Tipp"

    Das ist nur ein kleiner Vorgeschmack auf das, was du tun kannst, um deine Ressourcennutzung zu optimieren.
    Nextflow selbst hat wirklich raffinierte [dynamische Wiederholungslogik](https://www.nextflow.io/docs/latest/process.html#dynamic-task-resources) eingebaut, um Jobs, die aufgrund von Ressourcenbeschränkungen fehlschlagen, erneut zu versuchen.
    Zusätzlich bietet die Seqera Platform auch KI-gesteuerte Tools zur automatischen Optimierung deiner Ressourcenzuweisungen.

### 5.5. Ressourcenlimits hinzufügen

Abhängig davon, welchen Rechen-executor und welche Recheninfrastruktur du verwendest, kann es einige Einschränkungen geben, was du zuweisen kannst (oder musst).
Zum Beispiel kann dein Cluster verlangen, dass du innerhalb bestimmter Limits bleibst.

Du kannst die `resourceLimits`-Direktive verwenden, um die relevanten Beschränkungen zu setzen. Die Syntax sieht so aus, wenn sie allein in einem process-Block steht:

```groovy title="Syntax-Beispiel"
process {
    resourceLimits = [
        memory: 750.GB,
        cpus: 200,
        time: 30.d
    ]
}
```

Nextflow wird diese Werte in die entsprechenden Anweisungen übersetzen, abhängig vom executor, den du angegeben hast.

Wir werden das nicht ausführen, da wir keinen Zugang zu relevanter Infrastruktur in der Trainingsumgebung haben.
Wenn du jedoch versuchen würdest, den Workflow mit Ressourcenzuweisungen auszuführen, die diese Limits überschreiten, und dann den `sbatch`-Befehl in der `.command.run`-Script-Datei nachschlagen würdest, würdest du sehen, dass die Anforderungen, die tatsächlich an den executor gesendet werden, auf die durch `resourceLimits` angegebenen Werte begrenzt sind.

??? info "Institutionelle Referenzkonfigurationen"

    Das nf-core-Projekt hat eine [Sammlung von Konfigurationsdateien](https://nf-co.re/configs/) zusammengestellt, die von verschiedenen Institutionen auf der ganzen Welt geteilt werden und eine breite Palette von HPC- und Cloud-executors abdecken.

    Diese geteilten Configs sind sowohl für Leute wertvoll, die dort arbeiten und daher einfach die Konfiguration ihrer Institution sofort nutzen können, als auch als Modell für Leute, die eine Konfiguration für ihre eigene Infrastruktur entwickeln möchten.

### Zusammenfassung

Du weißt, wie du einen Profiling-Bericht generierst, um die Ressourcennutzung zu bewerten, und wie du Ressourcenzuweisungen für alle processes und/oder für einzelne processes änderst, sowie Ressourcenbeschränkungen für die Ausführung auf HPC setzt.

### Was kommt als Nächstes?

Lerne, wie du voreingestellte Konfigurationsprofile einrichtest und zur Laufzeit zwischen ihnen wechselst.

---

## 6. Profile verwenden, um zwischen voreingestellten Konfigurationen zu wechseln

??? example "Szenario"

    Du wechselst regelmäßig zwischen dem Ausführen von Pipelines auf deinem Laptop für die Entwicklung und auf dem HPC deiner Institution für Produktionsläufe.
    Du bist es leid, jedes Mal manuell Konfigurationseinstellungen zu ändern, wenn du die Umgebung wechselst.

Wir haben dir eine Reihe von Möglichkeiten gezeigt, wie du deine Pipeline-Konfiguration anpassen kannst, abhängig vom Projekt, an dem du arbeitest, oder der Rechenumgebung, die du verwendest.

Du möchtest vielleicht zwischen alternativen Einstellungen wechseln, abhängig davon, welche Recheninfrastruktur du verwendest. Zum Beispiel möchtest du vielleicht lokal auf deinem Laptop entwickeln und kleine Tests ausführen, dann vollständige Workloads auf HPC oder Cloud ausführen.

Nextflow ermöglicht es dir, beliebig viele Profile einzurichten, die verschiedene Konfigurationen beschreiben, die du dann zur Laufzeit mit einem Kommandozeilen-Argument auswählen kannst, anstatt die Konfigurationsdatei selbst ändern zu müssen.

### 6.1. Profile für das Wechseln zwischen lokaler Entwicklung und Ausführung auf HPC erstellen

Lass uns zwei alternative Profile einrichten; eines für das Ausführen kleiner Lasten auf einem normalen Computer, wo wir Docker-Container verwenden werden, und eines für das Ausführen auf einem Universitäts-HPC mit einem Slurm-Scheduler, wo wir Conda-Pakete verwenden werden.

#### 6.1.1. Die Profile einrichten

Füge das Folgende zu deiner `nextflow.config`-Datei hinzu, nach dem Abschnitt mit den Pipeline-Parametern, aber vor den Ausgabe-Einstellungen:

```groovy title="nextflow.config" linenums="24"
/*
* Profiles
*/
profiles {
    my_laptop {
        process.executor = 'local'
        docker.enabled = true
    }
    univ_hpc {
        process.executor = 'slurm'
        conda.enabled = true
        process.resourceLimits = [
            memory: 750.GB,
            cpus: 200,
            time: 30.d
        ]
    }
}
```

Du siehst, dass wir für das Universitäts-HPC auch Ressourcenbeschränkungen angeben.

#### 6.1.2. Den Workflow mit einem Profil ausführen

Um ein Profil in unserer Nextflow-Kommandozeile anzugeben, verwenden wir das `-profile`-Argument.

Lass uns versuchen, den Workflow mit der `my_laptop`-Konfiguration auszuführen.

```bash
nextflow run 3-main.nf -profile my_laptop
```

??? success "Befehlsausgabe"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `3-main.nf` [gigantic_brazil] DSL2 - revision: ede9037d02

    executor >  local (8)
    [58/da9437] sayHello (3)       | 3 of 3 ✔
    [35/9cbe77] convertToUpper (2) | 3 of 3 ✔
    [67/857d05] collectGreetings   | 1 of 1 ✔
    [37/7b51b5] cowpy              | 1 of 1 ✔
    ```

Wie du sehen kannst, erlaubt uns das, sehr bequem zur Laufzeit zwischen Konfigurationen zu wechseln.

!!! warning "Warnung"

    Das `univ_hpc`-Profil wird in der Trainingsumgebung nicht richtig laufen, da wir keinen Zugang zu einem Slurm-Scheduler haben.

Wenn wir in Zukunft andere Konfigurationselemente finden, die immer gemeinsam mit diesen auftreten, können wir sie einfach zum entsprechenden Profil/zu den entsprechenden Profilen hinzufügen.
Wir können auch zusätzliche Profile erstellen, wenn es andere Konfigurationselemente gibt, die wir zusammenfassen möchten.

### 6.2. Ein Profil mit Test-Parametern erstellen

??? example "Szenario"

    Du möchtest, dass andere deine Pipeline schnell ausprobieren können, ohne ihre eigenen Eingabedaten sammeln zu müssen.

Profile sind nicht nur für Infrastruktur-Konfiguration.
Wir können sie auch verwenden, um Standardwerte für Workflow-Parameter zu setzen, um es anderen einfacher zu machen, den Workflow auszuprobieren, ohne selbst geeignete Eingabewerte sammeln zu müssen.
Du kannst das als Alternative zur Verwendung einer Parameter-Datei betrachten.

#### 6.2.1. Das Profil einrichten

Die Syntax für das Ausdrücken von Standardwerten in diesem Kontext sieht so aus, für ein Profil, das wir `test` nennen:

```groovy title="Syntax-Beispiel"
    test {
        params.<parameter1>
        params.<parameter2>
        ...
    }
```

Wenn wir ein test-Profil für unseren Workflow hinzufügen, wird der `profiles`-Block:

```groovy title="nextflow.config" linenums="24"
/*
* Profiles
*/
profiles {
    my_laptop {
        process.executor = 'local'
        docker.enabled = true
    }
    univ_hpc {
        process.executor = 'slurm'
        conda.enabled = true
        process.resourceLimits = [
            memory: 750.GB,
            cpus: 200,
            time: 30.d
        ]
    }
    test {
        params.input = 'data/greetings.csv'
        params.batch = 'test'
        params.character = 'dragonandcow'
    }
}
```

Genau wie für technische Konfigurationsprofile kannst du mehrere verschiedene Profile einrichten, die Parameter unter einem beliebigen Namen angeben, den du magst.

#### 6.2.2. Den Workflow lokal mit dem test-Profil ausführen

Praktischerweise schließen sich Profile nicht gegenseitig aus, sodass wir mehrere Profile in unserer Kommandozeile mit der folgenden Syntax angeben können `-profile <profil1>,<profil2>` (für eine beliebige Anzahl von Profilen).

Wenn du Profile kombinierst, die Werte für dieselben Konfigurationselemente setzen und in derselben Konfigurationsdatei beschrieben sind, wird Nextflow den Konflikt lösen, indem es den Wert verwendet, den es zuletzt eingelesen hat (_d.h._ was auch immer später in der Datei kommt).
Wenn die widersprüchlichen Einstellungen in verschiedenen Konfigurationsquellen gesetzt sind, gilt die Standard-[Rangfolge](https://www.nextflow.io/docs/latest/config.html).

Lass uns versuchen, das test-Profil zu unserem vorherigen Befehl hinzuzufügen:

```bash
nextflow run 3-main.nf -profile my_laptop,test
```

??? success "Befehlsausgabe"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `3-main.nf` [jovial_coulomb] DSL2 - revision: 46a6763141

    executor >  local (8)
    [9b/687cdc] sayHello (2)       | 3 of 3 ✔
    [ca/552187] convertToUpper (3) | 3 of 3 ✔
    [e8/83e306] collectGreetings   | 1 of 1 ✔
    [fd/e84fa9] cowpy              | 1 of 1 ✔
    ```

Das wird Docker verwenden, wo möglich, und Ausgaben unter `results/test` produzieren, und diesmal ist der Charakter das komische Duo `dragonandcow`.

??? abstract "Dateiinhalte"

    ```console title="results/test/"
     _________
    / HOLà    \
    | HELLO   |
    \ BONJOUR /
    ---------
                \                    ^    /^
                  \                  / \  // \
                  \   |\___/|      /   \//  .\
                    \  /O  O  \__  /    //  | \ \           *----*
                      /     /  \/_/    //   |  \  \          \   |
                      \@___\@`    \/_   //    |   \   \         \/\ \
                    0/0/|       \/_ //     |    \    \         \ \
                0/0/0/0/|        \///      |     \     \       | |
              0/0/0/0/0/_|_ /   (  //       |      \     _\     |  /
          0/0/0/0/0/0/`/,_ _ _/  ) ; -.    |    _ _\.-~       /   /
                      ,-}        _      *-.|.-~-.           .~    ~
      \     \__/        `/\      /                 ~-. _ .-~      /
      \____(oo)           *.   }            {                   /
      (    (--)          .----~-.\        \-`                 .~
      //__\\  \__ Ack!   ///.----..<        \             _ -~
      //    \\               ///-._ _ _ _ _ _ _{^ - - - - ~
    ```

Das bedeutet, solange wir alle Testdatendateien mit dem Workflow-Code verteilen, kann jeder den Workflow schnell ausprobieren, ohne eigene Eingaben über die Kommandozeile oder eine Parameter-Datei bereitstellen zu müssen.

!!! tip "Tipp"

    Wir können für größere Dateien, die extern gespeichert sind, auf URLs verweisen.
    Nextflow wird sie automatisch herunterladen, solange eine offene Verbindung besteht.

    Für mehr Details, siehe die Side Quest [Working with Files](../side_quests/working_with_files.md)

### 6.3. `nextflow config` verwenden, um die aufgelöste Konfiguration zu sehen

Wie oben erwähnt, kann manchmal derselbe Parameter auf verschiedene Werte in Profilen gesetzt sein, die du kombinieren möchtest.
Und allgemeiner gibt es zahlreiche Orte, an denen Konfigurationselemente gespeichert werden können, und manchmal können dieselben Eigenschaften an verschiedenen Orten auf verschiedene Werte gesetzt sein.

Nextflow wendet eine festgelegte [Rangfolge](https://www.nextflow.io/docs/latest/config.html) an, um Konflikte zu lösen, aber das kann schwierig sein, selbst zu bestimmen.
Und selbst wenn nichts im Konflikt steht, kann es mühsam sein, alle möglichen Orte nachzuschlagen, an denen Dinge konfiguriert sein könnten.

Glücklicherweise enthält Nextflow ein praktisches Hilfsprogramm namens `config`, das diesen gesamten Prozess für dich automatisieren kann.

Das `config`-Tool wird alle Inhalte in deinem aktuellen Arbeitsverzeichnis erkunden, alle Konfigurationsdateien aufsaugen und die vollständig aufgelöste Konfiguration produzieren, die Nextflow verwenden würde, um den Workflow auszuführen.
Das ermöglicht es dir herauszufinden, welche Einstellungen verwendet werden, ohne etwas starten zu müssen.

#### 6.3.1. Die Standardkonfiguration auflösen

Führe diesen Befehl aus, um die Konfiguration aufzulösen, die standardmäßig angewendet würde.

```bash
nextflow config
```

??? success "Befehlsausgabe"

    ```groovy
    docker {
      enabled = false
    }

    conda {
      enabled = true
    }

    process {
      memory = '1 GB'
      withName:cowpy {
          memory = '2 GB'
          cpus = 2
      }
    }

    params {
      input = 'greetings.csv'
      batch = 'batch'
      character = 'turkey'
    }
    ```

Das zeigt dir die Basiskonfiguration, die du bekommst, wenn du nichts Extra in der Kommandozeile angibst.

#### 6.3.2. Die Konfiguration mit spezifischen aktivierten Einstellungen auflösen

Wenn du Kommandozeilen-Parameter angibst, z.B. ein oder mehrere Profile aktivierst oder eine Parameter-Datei lädst, wird der Befehl diese zusätzlich berücksichtigen.

```bash
nextflow config -profile my_laptop,test
```

??? success "Befehlsausgabe"

    ```groovy
    docker {
      enabled = true
    }

    conda {
      enabled = true
    }

    process {
      memory = '1 GB'
      withName:cowpy {
          memory = '2 GB'
          cpus = 2
      }
      executor = 'local'
    }

    params {
      input = 'greetings.csv'
      batch = 'test'
      character = 'dragonandcow'
    }
    ```

Das wird besonders nützlich für komplexe Projekte, die mehrere Konfigurationsebenen beinhalten.

### Zusammenfassung

Du weißt, wie du Profile verwendest, um zur Laufzeit mit minimalem Aufwand eine voreingestellte Konfiguration auszuwählen.
Allgemeiner weißt du, wie du deine Workflow-Ausführungen konfigurierst, um verschiedenen Rechenplattformen zu entsprechen und die Reproduzierbarkeit deiner Analysen zu verbessern.

### Was kommt als Nächstes?

Lerne, wie du Pipelines direkt aus Remote-Repositories wie GitHub ausführst.

---

## 7. Pipelines aus Remote-Repositories ausführen

??? example "Szenario"

    Du möchtest eine etablierte Pipeline wie die von nf-core ausführen, ohne den Code selbst herunterladen und verwalten zu müssen.

Bisher haben wir Workflow-Scripts ausgeführt, die sich im aktuellen Verzeichnis befinden.
In der Praxis wirst du oft Pipelines ausführen wollen, die in Remote-Repositories gespeichert sind, wie GitHub.

Nextflow macht das einfach: Du kannst jede Pipeline direkt von einer Git-Repository-URL ausführen, ohne sie vorher manuell herunterladen zu müssen.

### 7.1. Eine Pipeline von GitHub ausführen

Die grundlegende Syntax für das Ausführen einer Remote-Pipeline ist `nextflow run <repository>`, wobei `<repository>` ein GitHub-Repository-Pfad wie `nextflow-io/hello`, eine vollständige URL oder ein Pfad zu GitLab, Bitbucket oder anderen Git-Hosting-Diensten sein kann.

Versuche, die offizielle Nextflow "hello" Demo-Pipeline auszuführen:

```bash
nextflow run nextflow-io/hello
```

Das erste Mal, wenn du eine Remote-Pipeline ausführst, lädt Nextflow sie herunter und speichert sie lokal im cache.
Nachfolgende Ausführungen verwenden die gecachte Version, es sei denn, du forderst explizit ein Update an.

### 7.2. Eine Version für Reproduzierbarkeit angeben

Standardmäßig führt Nextflow die neueste Version vom Standard-Branch aus.
Du kannst eine bestimmte Version, einen Branch oder einen Commit mit dem `-r`-Flag angeben:

```bash
nextflow run nextflow-io/hello -r v1.1
```

Das Angeben exakter Versionen ist essentiell für die Reproduzierbarkeit.

### Zusammenfassung

Du weißt, wie du Pipelines direkt von GitHub und anderen Remote-Repositories ausführst, und wie du Versionen für die Reproduzierbarkeit angibst.

### Was kommt als Nächstes?

Klopf dir kräftig auf die Schulter!
Du weißt alles, was du wissen musst, um mit dem Ausführen und Verwalten von Nextflow Pipelines zu beginnen.

Das schließt diesen Kurs ab, aber wenn du eifrig weiterlernen möchtest, haben wir zwei Hauptempfehlungen:

- Wenn du tiefer in die Entwicklung eigener Pipelines eintauchen möchtest, schau dir [Hello Nextflow](../hello_nextflow/index.md) an, einen Einsteigerkurs, der denselben allgemeinen Verlauf wie dieser abdeckt, aber viel mehr ins Detail über channels und Operatoren geht.
- Wenn du weiterhin lernen möchtest, wie man Nextflow Pipelines ausführt, ohne tiefer in den Code einzusteigen, schau dir den ersten Teil von [Hello nf-core](../hello_nf-core/index.md) an, der die Tools zum Finden und Ausführen von Pipelines aus dem äußerst beliebten [nf-core](https://nf-co.re/) Projekt vorstellt.

Viel Spaß!

---

## Quiz

<quiz>
Wenn Parameter-Werte sowohl in der Workflow-Datei als auch in `nextflow.config` gesetzt sind, welcher hat Vorrang?
- [ ] Der Wert aus der Workflow-Datei
- [x] Der Wert aus der Konfigurationsdatei
- [ ] Der zuerst gefundene Wert
- [ ] Es verursacht einen Fehler

Mehr erfahren: [1.1. Set up values in `nextflow.config`](#11-set-up-values-in-nextflowconfig)
</quiz>

<quiz>
Was ist der Syntax-Unterschied zwischen dem Setzen eines Parameter-Standards in einer Workflow-Datei vs. einer config-Datei?
- [ ] Sie verwenden dieselbe Syntax
- [x] Workflow verwendet typisierte Deklaration (`#!groovy param: Type = value`), config verwendet Zuweisung (`#!groovy param = value`)
- [ ] Config verwendet typisierte Deklaration, Workflow verwendet Zuweisung
- [ ] Nur config-Dateien können Standardwerte setzen

Mehr erfahren: [1.1. Set up values in `nextflow.config`](#11-set-up-values-in-nextflowconfig)
</quiz>

<quiz>
Wie gibst du eine Parameter-Datei an, wenn du einen Workflow ausführst?
- [ ] `--params params.yaml`
- [ ] `-config params.yaml`
- [x] `-params-file params.yaml`
- [ ] `--input-params params.yaml`

Mehr erfahren: [1.3. Use a parameter file](#13-use-a-parameter-file)
</quiz>

<quiz>
Was steuert die `outputDir`-Konfigurationsoption?
- [ ] Den Ort des work directory
- [x] Den Basispfad, wo Workflow-Ausgaben veröffentlicht werden
- [ ] Das Verzeichnis für Log-Dateien
- [ ] Den Ort der Modul-Dateien

Mehr erfahren: [2.1. Customize the outputDir directory name](#21-customize-the-outputdir-directory-name)
</quiz>

<quiz>
Wie referenzierst du einen process-Namen dynamisch in der Ausgabe-Pfad-Konfiguration?
- [ ] `#!groovy ${processName}`
- [ ] `process.name`
- [x] `#!groovy { meta.id }`
- [ ] `@processName`

Mehr erfahren: [2.2. Organize outputs by process](#22-organize-outputs-by-process)
</quiz>

<quiz>
Wenn sowohl Docker als auch Conda aktiviert sind und ein process beide Direktiven hat, welche wird priorisiert?
- [x] Docker (containers)
- [ ] Conda
- [ ] Die zuerst im process definierte
- [ ] Es verursacht einen Fehler

Mehr erfahren: [3. Select a software packaging technology](#3-select-a-software-packaging-technology)
</quiz>

<quiz>
Was ist der Standard-executor in Nextflow?
- [x] `local`
- [ ] `slurm`
- [ ] `kubernetes`
- [ ] `aws`

Mehr erfahren: [4. Select an execution platform](#4-select-an-execution-platform)
</quiz>

<quiz>
Welcher Befehl generiert einen Ressourcennutzungsbericht?
- [ ] `nextflow run workflow.nf -with-metrics`
- [ ] `nextflow run workflow.nf -with-stats`
- [x] `nextflow run workflow.nf -with-report report.html`
- [ ] `nextflow run workflow.nf -profile report`

Mehr erfahren: [5.1. Run the workflow to generate a resource utilization report](#51-run-the-workflow-to-generate-a-resource-utilization-report)
</quiz>

<quiz>
Wie setzt du Ressourcenanforderungen für einen bestimmten process namens `cowpy` in der config-Datei?
- [ ] `#!groovy cowpy.memory = '2.GB'`
- [ ] `#!groovy process.cowpy.memory = '2.GB'`
- [x] `#!groovy process { withName: 'cowpy' { memory = '2.GB' } }`
- [ ] `#!groovy resources.cowpy.memory = '2.GB'`

Mehr erfahren: [5.3. Set resource allocations for a specific process](#53-set-resource-allocations-for-a-specific-process)
</quiz>

<quiz>
Was macht die `resourceLimits`-Direktive?
- [ ] Setzt minimale Ressourcenanforderungen
- [ ] Weist Ressourcen zu processes zu
- [x] Begrenzt die maximalen Ressourcen, die angefordert werden können
- [ ] Überwacht die Ressourcennutzung in Echtzeit

Mehr erfahren: [5.5. Add resource limits](#55-add-resource-limits)
</quiz>

<quiz>
Wie gibst du mehrere Profile in einem einzelnen Befehl an?
- [ ] `-profile profile1 -profile profile2`
- [ ] `-profiles profile1,profile2`
- [x] `-profile profile1,profile2`
- [ ] `--profile profile1 --profile profile2`

Mehr erfahren: [6. Use profiles to switch between preset configurations](#6-use-profiles-to-switch-between-preset-configurations)
</quiz>

<quiz>
Welcher Befehl zeigt die vollständig aufgelöste Konfiguration, die Nextflow verwenden würde?
- [ ] `nextflow show-config`
- [ ] `nextflow settings`
- [x] `nextflow config`
- [ ] `nextflow resolve`

Mehr erfahren: [6.3. Use `nextflow config` to see the resolved configuration](#63-use-nextflow-config-to-see-the-resolved-configuration)
</quiz>

<quiz>
Wofür können Profile verwendet werden? (Wähle alle zutreffenden aus)
- [x] Definieren von infrastrukturspezifischen Einstellungen (executors, containers)
- [x] Setzen von Ressourcenlimits für verschiedene Umgebungen
- [x] Bereitstellen von Test-Parametern für einfaches Workflow-Testen
- [ ] Definieren neuer processes

Mehr erfahren: [6. Use profiles to switch between preset configurations](#6-use-profiles-to-switch-between-preset-configurations)
</quiz>
