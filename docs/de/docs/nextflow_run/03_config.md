# Teil 3: Konfiguration der Ausführung

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } KI-gestützte Übersetzung - [mehr erfahren & Verbesserungen vorschlagen](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Dieser Abschnitt zeigt dir, wie du die Konfiguration einer Nextflow-Pipeline verwaltest, um ihr Verhalten anzupassen, sie an verschiedene Umgebungen anzupassen und die Ressourcennutzung zu optimieren – _ohne eine einzige Zeile des Workflow-Codes selbst zu ändern_.

Es gibt mehrere Möglichkeiten, dies zu tun, die kombiniert werden können und gemäß der in der [Konfigurations](https://nextflow.io/docs/latest/config.html)-Dokumentation beschriebenen Rangfolge interpretiert werden.

In diesem Teil des Kurses zeigen wir dir den einfachsten und gebräuchlichsten Mechanismus für Konfigurationsdateien: die `nextflow.config`-Datei, die du bereits im Abschnitt über Container in Teil 2 kennengelernt hast.

Wir werden wesentliche Komponenten der Nextflow-Konfiguration wie Prozess-Direktiven, Executors, Profile und Parameterdateien durchgehen.
Indem du lernst, diese Konfigurationsoptionen effektiv zu nutzen, kannst du die Flexibilität, Skalierbarkeit und Leistung von Nextflow-Pipelines voll ausschöpfen.

Um diese Konfigurationselemente zu üben, werden wir eine neue Kopie des Workflows ausführen, den wir zuletzt am Ende von Teil 2 dieses Trainingskurses ausgeführt haben, umbenannt in `3-main.nf`.

Falls du mit der Hello-Pipeline nicht vertraut bist oder eine Auffrischung gebrauchen könntest, siehe [diese Infoseite](../info/hello_pipeline.md).

---

## 1. Workflow-Eingabeparameter verwalten

??? example "Szenario"

    Du hast eine Pipeline heruntergeladen und möchtest sie wiederholt mit denselben Eingabedateien und Einstellungen ausführen, aber du möchtest nicht jedes Mal alle Parameter eingeben.
    Oder vielleicht richtest du die Pipeline für eine\*n Kolleg\*in ein, der\*die sich mit Kommandozeilenargumenten nicht wohl fühlt.

Wir beginnen mit einem Aspekt der Konfiguration, der einfach eine Erweiterung dessen ist, womit wir bisher gearbeitet haben: die Verwaltung von Eingabeparametern.

Derzeit ist unser Workflow so eingerichtet, dass er mehrere Parameterwerte über die Kommandozeile akzeptiert, die in einem `params`-Block im Workflow-Skript selbst deklariert sind.
Einer hat einen Standardwert, der als Teil seiner Deklaration festgelegt ist.

Du möchtest jedoch möglicherweise Standardwerte für alle festlegen oder den vorhandenen Standardwert überschreiben, ohne Parameter in der Kommandozeile angeben oder die ursprüngliche Skriptdatei ändern zu müssen.

Es gibt mehrere Möglichkeiten, dies zu tun; wir zeigen dir drei grundlegende Methoden, die sehr häufig verwendet werden.

### 1.1. Werte in `nextflow.config` einrichten

Dies ist der einfachste Ansatz, obwohl er möglicherweise der am wenigsten flexible ist, da die Haupt-`nextflow.config`-Datei nichts ist, was du für jeden Lauf bearbeiten möchtest.
Aber es hat den Vorteil, die Belange der _Deklaration_ der Parameter im Workflow (die definitiv dorthin gehört) von der Bereitstellung von _Standardwerten_ zu trennen, die eher in einer Konfigurationsdatei zu Hause sind.

Lass uns dies in zwei Schritten tun.

#### 1.1.1. Einen `params`-Block in der Konfigurationsdatei erstellen

Nimm die folgenden Codeänderungen in der `nextflow.config`-Datei vor:

=== "Danach"

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

Beachte, dass wir den `params`-Block nicht einfach vom Workflow in die Konfigurationsdatei kopiert haben.
Für den `batch`-Parameter, der bereits einen Standardwert deklariert hatte, ist die Syntax etwas anders.
In der Workflow-Datei ist das eine typisierte Deklaration.
In der Konfiguration sind das Wertzuweisungen.

Technisch gesehen reicht dies aus, um die Standardwerte zu überschreiben, die noch in der Workflow-Datei angegeben sind.
Du könntest den Standardwert für `batch` ändern und den Workflow ausführen, um dich zu überzeugen, dass der in der Konfigurationsdatei festgelegte Wert den in der Workflow-Datei festgelegten überschreibt.

Aber im Sinne der vollständigen Verlagerung der Konfiguration in die Konfigurationsdatei, lass uns diesen Standardwert vollständig aus der Workflow-Datei entfernen.

#### 1.1.2. Den Standardwert für `batch` in der Workflow-Datei entfernen

Nimm die folgende Codeänderung in der `3-main.nf`-Workflow-Datei vor:

=== "Danach"

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

Jetzt setzt die Workflow-Datei selbst keine Standardwerte für diese Parameter.

#### 1.1.3. Die Pipeline ausführen

Lass uns testen, dass es korrekt funktioniert, ohne Parameter in der Kommandozeile anzugeben.

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

Dies erzeugt immer noch dieselbe Ausgabe wie zuvor.

Die finale ASCII-Art-Ausgabe befindet sich im Verzeichnis `results/3-main/` unter dem Namen `cowpy-COLLECTED-batch-output.txt`, wie zuvor.

??? abstract "Dateiinhalt"

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

Funktional hat diese Änderung nichts geändert, aber konzeptionell ist es etwas sauberer, die Standardwerte in der Konfigurationsdatei festzulegen.

### 1.2. Eine laufspezifische Konfigurationsdatei verwenden

??? example "Szenario"

    Du möchtest mit verschiedenen Einstellungen experimentieren, ohne deine Hauptkonfigurationsdatei zu ändern.

Du kannst dies tun, indem du eine neue `nextflow.config`-Datei in einem Unterverzeichnis erstellst, das du als Arbeitsverzeichnis für deine Experimente verwenden wirst.

#### 1.2.1. Das Arbeitsverzeichnis mit einer leeren Konfiguration erstellen

Lass uns zunächst ein neues Verzeichnis erstellen und hineinwechseln:

```bash
mkdir -p tux-run
cd tux-run
```

Erstelle dann eine leere Konfigurationsdatei in diesem Verzeichnis:

```bash
touch nextflow.config
```

Dies erzeugt eine leere Datei.

#### 1.2.2. Die experimentelle Konfiguration einrichten

Öffne nun die neue Datei und füge die Parameter hinzu, die du anpassen möchtest:

```groovy title="tux-run/nextflow.config" linenums="1"
params {
    input = '../data/greetings.csv'
    batch = 'experiment'
    character = 'tux'
}
```

Beachte, dass der Pfad zur Eingabedatei die Verzeichnisstruktur widerspiegeln muss.

#### 1.2.3. Die Pipeline ausführen

Wir können jetzt unsere Pipeline aus unserem neuen Arbeitsverzeichnis heraus ausführen.
Achte darauf, den Pfad entsprechend anzupassen!

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

Dies erstellt einen neuen Satz von Verzeichnissen unter `tux-run/`, einschließlich `tux-run/work/` und `tux-run/results/`.

In diesem Lauf kombiniert Nextflow die `nextflow.config` in unserem aktuellen Verzeichnis mit der `nextflow.config` im Stammverzeichnis der Pipeline und überschreibt dadurch den Standardcharakter (turkey) mit dem tux-Charakter.

Die finale Ausgabedatei sollte den tux-Charakter enthalten, der die Grüße sagt.

??? abstract "Dateiinhalt"

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

!!! warning

    Stelle sicher, dass du zum vorherigen Verzeichnis zurückwechselst, bevor du zum nächsten Abschnitt übergehst!

    ```bash
    cd ..
    ```

Schauen wir uns nun eine weitere nützliche Möglichkeit an, Parameterwerte festzulegen.

### 1.3. Eine Parameterdatei verwenden

??? example "Szenario"

    Du musst genaue Laufparameter mit einer\*m Kolleg\*in teilen oder sie für eine Veröffentlichung aufzeichnen.

Der Unterverzeichnis-Ansatz funktioniert hervorragend zum Experimentieren, erfordert aber etwas Einrichtung und dass du Pfade entsprechend anpasst.
Es gibt einen einfacheren Ansatz, wenn du deine Pipeline mit einem bestimmten Satz von Werten ausführen möchtest oder es jemand anderem mit minimalem Aufwand ermöglichen möchtest.

Nextflow ermöglicht es uns, Parameter über eine [Parameterdatei](https://nextflow.io/docs/latest/config.html#parameter-file) im YAML- oder JSON-Format anzugeben, was es sehr bequem macht, alternative Sätze von Standardwerten zu verwalten und zu verteilen, zum Beispiel, sowie laufspezifische Parameterwerte.

#### 1.3.1. Die Beispiel-Parameterdatei untersuchen

Um dies zu demonstrieren, stellen wir eine Beispiel-Parameterdatei im aktuellen Verzeichnis bereit, genannt `test-params.yaml`:

```yaml title="test-params.yaml" linenums="1"
input: "data/greetings.csv"
batch: "yaml"
character: "stegosaurus"
```

Diese Parameterdatei enthält ein Schlüssel-Wert-Paar für jede der Eingaben, die wir angeben möchten.
Beachte die Verwendung von Doppelpunkten (`:`) anstelle von Gleichheitszeichen (`=`), wenn du die Syntax mit der Konfigurationsdatei vergleichst.
Die Config-Datei ist in Groovy geschrieben, während die Parameterdatei in YAML geschrieben ist.

!!! info

    Wir stellen auch eine JSON-Version der Parameterdatei als Beispiel bereit, aber wir werden sie hier nicht ausführen.
    Probiere diese gerne selbst aus.

#### 1.3.2. Die Pipeline ausführen

Um den Workflow mit dieser Parameterdatei auszuführen, füge einfach `-params-file <dateiname>` zum Basisbefehl hinzu.

```bash
nextflow run 3-main.nf -params-file test-params.yaml
```

??? success "Befehlsausgabe"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `3-main.nf` [disturbed_sammet] DSL2 - revision: ede9037d02

    executor >  local (8)
    [2b/9a7d1e] sayHello (2)       | 3 of 3 ✔
    [5c/8f3b2a] convertToUpper (3) | 3 of 3 ✔
    [a3/29d8fb] collectGreetings   | 1 of 1 ✔
    [b7/83ef12] cowpy              | 1 of 1 ✔
    ```

Die finale Ausgabedatei sollte den Stegosaurus-Charakter enthalten, der die Grüße sagt.

??? abstract "Dateiinhalt"

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

Die Verwendung einer Parameterdatei mag übertrieben erscheinen, wenn du nur wenige Parameter anzugeben hast, aber einige Pipelines erwarten Dutzende von Parametern.
In diesen Fällen ermöglicht uns die Verwendung einer Parameterdatei, Parameterwerte zur Laufzeit bereitzustellen, ohne massive Kommandozeilen eingeben und ohne das Workflow-Skript ändern zu müssen.

Es erleichtert auch die Verteilung von Parametersätzen an Kolleg\*innen oder als unterstützende Informationen für eine Veröffentlichung, zum Beispiel.
Dies macht deine Arbeit reproduzierbarer für andere.

### Fazit

Du weißt, wie du wichtige Konfigurationsoptionen zur Verwaltung von Workflow-Eingaben nutzen kannst.

### Wie geht es weiter?

Lerne, wie du verwaltest, wo und wie deine Workflow-Ausgaben veröffentlicht werden.

---

## 2. Workflow-Ausgaben verwalten

??? example "Szenario"

    Deine Pipeline veröffentlicht Ausgaben in ein fest codiertes Verzeichnis, aber du möchtest Ergebnisse nach Projekt- oder Experimentnamen organisieren, ohne jedes Mal den Workflow-Code zu bearbeiten.

Der Workflow, den wir geerbt haben, verwendet Pfade für Ausgabedeklarationen auf Workflow-Ebene, was nicht besonders flexibel ist und viel Wiederholung beinhaltet.

Schauen wir uns einige gängige Möglichkeiten an, wie du dies flexibler konfigurieren kannst.

### 2.1. Den `outputDir`-Verzeichnisnamen anpassen

Jede Version des Workflows, die wir bisher ausgeführt haben, hat ihre Ausgaben in ein anderes Unterverzeichnis veröffentlicht, das fest in die Ausgabedefinitionen codiert ist.

Wir haben geändert, wo sich dieses Unterverzeichnis befand, indem wir in Teil 1 das CLI-Flag `-output-dir` verwendet haben, aber das ist immer noch nur ein statischer String.
Lass uns dies stattdessen in einer Config-Datei konfigurieren, wo wir komplexere dynamische Pfade definieren können.
Wir könnten einen ganz neuen Parameter dafür erstellen, aber lass uns den `batch`-Parameter verwenden, da er direkt verfügbar ist.

#### 2.1.1. Einen Wert für `outputDir` in der Konfigurationsdatei festlegen

Der Pfad, den Nextflow zum Veröffentlichen von Ausgaben verwendet, wird durch die `outputDir`-Option gesteuert.
Um den Pfad für alle Ausgaben zu ändern, kannst du einen Wert für diese Option in der `nextflow.config`-Konfigurationsdatei festlegen.

Füge den folgenden Code zur `nextflow.config`-Datei hinzu:

=== "Danach"

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
    * Ausgabeeinstellungen
    */
    outputDir = "results_config/${params.batch}"
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

Dies ersetzt den eingebauten Standardpfad `results/` durch `results_config/` plus den Wert des `batch`-Parameters als Unterverzeichnis.

Denke daran, dass du diese Option auch über die Kommandozeile mit dem Parameter `-output-dir` in deinem Befehl (`-o` als Kurzform) setzen kannst, aber dann könntest du den `batch`-Parameterwert nicht verwenden.
Die Verwendung des CLI-Flags überschreibt `outputDir` in der Config, falls es gesetzt ist.

#### 2.1.2. Den wiederholten Teil des fest codierten Pfads entfernen

Wir haben immer noch ein Unterverzeichnis fest in den Ausgabeoptionen codiert, also lass uns das jetzt entfernen.

Nimm die folgenden Codeänderungen in der Workflow-Datei vor:

=== "Danach"

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

Wir hätten auch einfach `${params.batch}` zu jedem Pfad hinzufügen können, anstatt den `outputDir`-Standard zu ändern, aber dies ist prägnanter.

#### 2.1.3. Die Pipeline ausführen

Lass uns testen, dass es korrekt funktioniert, indem wir den Batch-Namen von der Kommandozeile aus auf `outdir` setzen.

```bash
nextflow run 3-main.nf --batch outdir
```

??? success "Befehlsausgabe"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `3-main.nf` [amazing_church] DSL2 - revision: 6e18cd130e

    executor >  local (8)
    [9c/6a03ea] sayHello (2)       [100%] 3 of 3 ✔
    [11/9e58a6] convertToUpper (3) [100%] 3 of 3 ✔
    [c8/1977e5] collectGreetings   [100%] 1 of 1 ✔
    [38/f01eda] cowpy              [100%] 1 of 1 ✔
    ```

Dies erzeugt immer noch dieselbe Ausgabe wie zuvor, außer dass wir diesmal unsere Ausgaben unter `results_config/outdir/` finden.

??? abstract "Verzeichnisinhalt"

    ```console
    results_config/outdir
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

Du kannst diesen Ansatz mit benutzerdefinierten Pfaddefinitionen kombinieren, um jede gewünschte Verzeichnishierarchie zu konstruieren.

### 2.2. Ausgaben nach Prozess organisieren

Eine beliebte Möglichkeit, Ausgaben weiter zu organisieren, besteht darin, dies nach Prozess zu tun, _d.h._ Unterverzeichnisse für jeden in der Pipeline ausgeführten Prozess zu erstellen.

#### 2.2.1. Die Ausgabepfade durch einen Verweis auf Prozessnamen ersetzen

Alles, was du tun musst, ist den Namen des Prozesses als `<prozess>.name` in der Ausgabepfaddeklaration zu referenzieren.

Nimm die folgenden Änderungen in der Workflow-Datei vor:

=== "Danach"

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

Dies entfernt die verbleibenden fest codierten Elemente aus der Ausgabepfadkonfiguration.

#### 2.2.2. Die Pipeline ausführen

Lass uns testen, dass es korrekt funktioniert, indem wir den Batch-Namen von der Kommandozeile aus auf `pnames` setzen.

```bash
nextflow run 3-main.nf --batch pnames
```

??? success "Befehlsausgabe"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `3-main.nf` [jovial_mcclintock] DSL2 - revision: ede9037d02

    executor >  local (8)
    [4a/c2e6b8] sayHello (2)       | 3 of 3 ✔
    [6f/d4a172] convertToUpper (3) | 3 of 3 ✔
    [e8/4f19d7] collectGreetings   | 1 of 1 ✔
    [f2/a85c36] cowpy              | 1 of 1 ✔
    ```

Dies erzeugt immer noch dieselbe Ausgabe wie zuvor, außer dass wir diesmal unsere Ausgaben unter `results_config/pnames/` finden und sie nach Prozess gruppiert sind.

??? abstract "Verzeichnisinhalt"

    ```console
    results_config/pnames/
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

!!! note

    Beachte, dass wir hier die Unterscheidung zwischen `intermediates` und finalen Ausgaben auf der obersten Ebene gelöscht haben.
    Du kannst diese Ansätze mischen und sogar mehrere Variablen einbeziehen, zum Beispiel indem du den Pfad der ersten Ausgabe als `#!groovy "${params.batch}/intermediates/${sayHello.name}"` festlegst

### 2.3. Den Veröffentlichungsmodus auf Workflow-Ebene festlegen

Schließlich können wir im Sinne der Reduzierung von sich wiederholendem Code die `mode`-Deklarationen pro Ausgabe durch eine einzige Zeile in der Konfiguration ersetzen.

#### 2.3.1. `workflow.output.mode` zur Konfigurationsdatei hinzufügen

Füge den folgenden Code zur `nextflow.config`-Datei hinzu:

=== "Danach"

    ```groovy title="nextflow.config" linenums="2" hl_lines="5"
    /*
    * Ausgabeeinstellungen
    */
    outputDir = "results_config/${params.batch}"
    workflow.output.mode = 'copy'
    ```

=== "Vorher"

    ```groovy title="nextflow.config" linenums="12"
    /*
    * Ausgabeeinstellungen
    */
    outputDir = "results_config/${params.batch}"
    ```

Genau wie bei der `outputDir`-Option würde es ausreichen, `workflow.output.mode` einen Wert in der Konfigurationsdatei zu geben, um zu überschreiben, was in der Workflow-Datei festgelegt ist, aber lass uns den unnötigen Code trotzdem entfernen.

#### 2.3.2. Ausgabemodus aus der Workflow-Datei entfernen

Nimm die folgenden Änderungen in der Workflow-Datei vor:

=== "Danach"

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

Lass uns testen, dass es korrekt funktioniert, indem wir den Batch-Namen von der Kommandozeile aus auf `outmode` setzen.

```bash
nextflow run 3-main.nf --batch outmode
```

??? success "Befehlsausgabe"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `3-main.nf` [rowdy_sagan] DSL2 - revision: ede9037d02

    executor >  local (8)
    [5b/d91e3c] sayHello (2)       | 3 of 3 ✔
    [8a/f6c241] convertToUpper (3) | 3 of 3 ✔
    [89/cd3a48] collectGreetings   | 1 of 1 ✔
    [9e/71fb52] cowpy              | 1 of 1 ✔
    ```

Dies erzeugt immer noch dieselbe Ausgabe wie zuvor, außer dass wir diesmal unsere Ausgaben unter `results_config/outmode/` finden.
Sie sind immer noch alle echte Kopien, keine Symlinks.

??? abstract "Verzeichnisinhalt"

    ```console
    results_config/outmode/
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

Der Hauptgrund, warum du möglicherweise immer noch die Methode pro Ausgabe zum Festlegen des Modus verwenden möchtest, ist, wenn du innerhalb desselben Workflows mischen und anpassen möchtest, _d.h._ einige Ausgaben kopiert und einige als Symlinks haben möchtest.

Es gibt viele andere Optionen, die du auf diese Weise anpassen kannst, aber hoffentlich gibt dir dies ein Gefühl für die Bandbreite der Optionen und wie du sie effektiv nutzen kannst, um deinen Vorlieben zu entsprechen.

### Fazit

Du weißt, wie du die Benennung und Struktur der Verzeichnisse steuerst, in denen deine Ausgaben veröffentlicht werden, sowie den Workflow-Ausgabeveröffentlichungsmodus.

### Wie geht es weiter?

Lerne, wie du deine Workflow-Konfiguration an deine Rechenumgebung anpasst, beginnend mit der Software-Packaging-Technologie.

---

## 3. Eine Software-Packaging-Technologie auswählen

Bisher haben wir uns Konfigurationselemente angesehen, die steuern, wie Eingaben hineingehen und wo Ausgaben herauskommen. Jetzt ist es an der Zeit, sich spezifischer darauf zu konzentrieren, deine Workflow-Konfiguration an deine Rechenumgebung anzupassen.

Der erste Schritt auf diesem Weg ist die Angabe, woher die Softwarepakete kommen, die in jedem Schritt ausgeführt werden.
Sind sie bereits in der lokalen Rechenumgebung installiert?
Müssen wir Images abrufen und sie über ein Container-System ausführen?
Oder müssen wir Conda-Pakete abrufen und eine lokale Conda-Umgebung erstellen?

Im allerersten Teil dieses Trainingskurses (Teile 1-4) haben wir einfach lokal installierte Software in unserem Workflow verwendet.
Dann haben wir in Teil 5 Docker-Container und die `nextflow.config`-Datei eingeführt, die wir verwendet haben, um die Verwendung von Docker-Containern zu aktivieren.

Schauen wir uns nun an, wie wir eine alternative Software-Packaging-Option über die `nextflow.config`-Datei konfigurieren können.

### 3.1. Docker deaktivieren und Conda in der Config-Datei aktivieren

??? example "Szenario"

    Du verschiebst deine Pipeline auf einen HPC-Cluster, wo Docker aus Sicherheitsgründen nicht erlaubt ist.
    Der Cluster unterstützt Singularity und Conda, also musst du deine Konfiguration entsprechend umstellen.

Wie bereits erwähnt, unterstützt Nextflow mehrere Container-Technologien, einschließlich Singularity (das auf HPC weiter verbreitet ist), sowie Software-Paketmanager wie Conda.

Wir können unsere Konfigurationsdatei ändern, um Conda anstelle von Docker zu verwenden.
Dazu ändern wir den Wert von `docker.enabled` auf `false` und fügen eine Direktive hinzu, die die Verwendung von Conda aktiviert:

=== "Danach"

    ```groovy title="nextflow.config" linenums="1" hl_lines="1-2"
    docker.enabled = false
    conda.enabled = true
    ```

=== "Vorher"

    ```groovy title="nextflow.config" linenums="1" hl_lines="1"
    docker.enabled = true
    ```

Dies ermöglicht es Nextflow, Conda-Umgebungen für Prozesse zu erstellen und zu nutzen, die Conda-Pakete angegeben haben.
Was bedeutet, dass wir jetzt eines davon zu unserem `cowpy`-Prozess hinzufügen müssen!

### 3.2. Ein Conda-Paket in der Prozessdefinition angeben

Wir haben bereits die URI für ein Conda-Paket abgerufen, das das `cowpy`-Tool enthält: `conda-forge::cowpy==1.1.5`

Jetzt fügen wir die URI zur `cowpy`-Prozessdefinition mit der `conda`-Direktive hinzu:

=== "Danach"

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

Um es klar zu sagen, wir _ersetzen_ die `docker`-Direktive nicht, wir _fügen_ eine alternative Option hinzu.

!!! tip

    Es gibt einige verschiedene Möglichkeiten, die URI für ein bestimmtes Conda-Paket zu erhalten.
    Wir empfehlen die Verwendung der [Seqera Containers](https://seqera.io/containers/)-Suchabfrage, die dir eine URI gibt, die du kopieren und einfügen kannst, auch wenn du nicht planst, einen Container daraus zu erstellen.

### 3.3. Den Workflow ausführen, um zu überprüfen, dass er Conda verwenden kann

Lass es uns ausprobieren.

```bash
nextflow run 3-main.nf --batch conda
```

??? success "Befehlsausgabe"

    ```console title="Output"
    N E X T F L O W   ~  version 25.10.2

    Launching `3-main.nf` [trusting_lovelace] DSL2 - revision: 028a841db1

    executor >  local (8)
    [ee/4ca1f2] sayHello (3)       | 3 of 3 ✔
    [20/2596a7] convertToUpper (1) | 3 of 3 ✔
    [b3/e15de5] collectGreetings   | 1 of 1 ✔
    [c5/af5f88] cowpy              | 1 of 1 ✔
    ```

Dies sollte problemlos funktionieren und dieselben Ausgaben wie zuvor unter `results_config/conda` erzeugen.

Hinter den Kulissen hat Nextflow die Conda-Pakete abgerufen und die Umgebung erstellt, was normalerweise etwas Arbeit erfordert; es ist also schön, dass wir das nicht selbst machen müssen!

!!! info

    Dies läuft schnell, weil das `cowpy`-Paket recht klein ist, aber wenn du mit großen Paketen arbeitest, kann es beim ersten Mal etwas länger als üblich dauern, und du siehst möglicherweise, dass die Konsolenausgabe für eine Minute oder so 'hängen bleibt', bevor sie abgeschlossen wird.
    Dies ist normal und liegt an der zusätzlichen Arbeit, die Nextflow beim ersten Mal leistet, wenn du ein neues Paket verwendest.

Aus unserer Sicht sieht es so aus, als würde es genauso funktionieren wie die Ausführung mit Docker, obwohl die Mechanik im Backend etwas anders ist.

Das bedeutet, dass wir bereit sind, bei Bedarf mit Conda-Umgebungen zu arbeiten.

??? info "Docker und Conda mischen"

    Da diese Direktiven pro Prozess zugewiesen werden, ist es möglich, zu 'mischen', _d.h._ einige der Prozesse in deinem Workflow so zu konfigurieren, dass sie mit Docker laufen, und andere mit Conda, zum Beispiel, wenn die von dir verwendete Recheninfrastruktur beide unterstützt.
    In diesem Fall würdest du sowohl Docker als auch Conda in deiner Konfigurationsdatei aktivieren.
    Wenn beide für einen bestimmten Prozess verfügbar sind, priorisiert Nextflow Container.

    Und wie bereits erwähnt, unterstützt Nextflow mehrere andere Software-Packaging- und Container-Technologien, sodass du nicht nur auf diese beiden beschränkt bist.

### Fazit

Du weißt, wie du konfigurierst, welches Softwarepaket jeder Prozess verwenden soll, und wie du zwischen Technologien wechselst.

### Wie geht es weiter?

Lerne, wie du die von Nextflow verwendete Ausführungsplattform änderst, um die Arbeit tatsächlich zu erledigen.

---

## 4. Eine Ausführungsplattform auswählen

??? example "Szenario"

    Du hast deine Pipeline auf deinem Laptop entwickelt und getestet, aber jetzt musst du sie auf Tausenden von Proben ausführen.
    Deine Institution hat einen HPC-Cluster mit einem Slurm-Scheduler, den du stattdessen verwenden möchtest.

Bis jetzt haben wir unsere Pipeline mit dem lokalen Executor ausgeführt.
Dieser führt jede Aufgabe auf der Maschine aus, auf der Nextflow läuft.
Wenn Nextflow startet, schaut es sich die verfügbaren CPUs und den Speicher an.
Wenn die Ressourcen der zur Ausführung bereiten Aufgaben die verfügbaren Ressourcen überschreiten, hält Nextflow die letzten Aufgaben von der Ausführung zurück, bis eine oder mehrere der früheren Aufgaben beendet sind und die notwendigen Ressourcen freigeben.

Der lokale Executor ist bequem und effizient, aber er ist auf diese einzelne Maschine beschränkt. Für sehr große Workloads kannst du feststellen, dass deine lokale Maschine ein Engpass ist, entweder weil du eine einzelne Aufgabe hast, die mehr Ressourcen erfordert, als du verfügbar hast, oder weil du so viele Aufgaben hast, dass das Warten darauf, dass eine einzelne Maschine sie ausführt, zu lange dauern würde.

Nextflow unterstützt [viele verschiedene Ausführungs-Backends](https://nextflow.io/docs/latest/executor.html), einschließlich HPC-Scheduler (Slurm, LSF, SGE, PBS, Moab, OAR, Bridge, HTCondor und andere) sowie Cloud-Ausführungs-Backends (AWS Batch, Google Cloud Batch, Azure Batch, Kubernetes und mehr).

### 4.1. Ein anderes Backend anvisieren

Die Wahl des Executors wird durch eine Prozess-Direktive namens `executor` festgelegt.
Standardmäßig ist sie auf `local` gesetzt, sodass die folgende Konfiguration impliziert ist:

```groovy title="Eingebaute Konfiguration"
process {
    executor = 'local'
}
```

Um den Executor so einzustellen, dass er ein anderes Backend anvisiert, würdest du einfach den gewünschten Executor mit ähnlicher Syntax angeben, wie oben für Ressourcenzuweisungen beschrieben (siehe [Executors](https://nextflow.io/docs/latest/executor.html) für alle Optionen).

```groovy title="nextflow.config"
process {
    executor = 'slurm'
}
```

!!! warning

    Wir können dies in der Trainingsumgebung nicht tatsächlich testen, da sie nicht für die Verbindung zu einem HPC eingerichtet ist.

### 4.2. Umgang mit backend-spezifischer Syntax für Ausführungsparameter

Die meisten Hochleistungsrechenplattformen erlauben (und erfordern manchmal), dass du bestimmte Parameter wie Ressourcenzuweisungsanfragen und -beschränkungen (z.B. Anzahl der CPUs und Speicher) und den Namen der zu verwendenden Job-Warteschlange angibst.

Leider verwendet jedes dieser Systeme unterschiedliche Technologien, Syntaxen und Konfigurationen, um zu definieren, wie ein Job definiert und an den entsprechenden Scheduler übermittelt werden soll.

??? abstract "Beispiele"

    Zum Beispiel muss derselbe Job, der 8 CPUs und 4 GB RAM erfordert und in der Warteschlange "my-science-work" ausgeführt werden soll, je nach Backend auf folgende unterschiedliche Weise ausgedrückt werden.

    ```bash title="Config für SLURM / Übermittlung mit sbatch"
    #SBATCH -o /path/to/my/task/directory/my-task-1.log
    #SBATCH --no-requeue
    #SBATCH -c 8
    #SBATCH --mem 4096M
    #SBATCH -p my-science-work
    ```

    ```bash title="Config für PBS / Übermittlung mit qsub"
    #PBS -o /path/to/my/task/directory/my-task-1.log
    #PBS -j oe
    #PBS -q my-science-work
    #PBS -l nodes=1:ppn=5
    #PBS -l mem=4gb
    ```

    ```bash title="Config für SGE / Übermittlung mit qsub"
    #$ -o /path/to/my/task/directory/my-task-1.log
    #$ -j y
    #$ -terse
    #$ -notify
    #$ -q my-science-work
    #$ -l slots=5
    #$ -l h_rss=4096M,mem_free=4096M
    ```

Glücklicherweise vereinfacht Nextflow all dies.
Es bietet eine standardisierte Syntax, sodass du die relevanten Eigenschaften wie `cpus`, `memory` und `queue` nur einmal angeben kannst (siehe [Prozess-Direktiven](https://nextflow.io/docs/latest/reference/process.html#process-directives) für alle verfügbaren Optionen).
Dann generiert Nextflow zur Laufzeit die entsprechenden backend-spezifischen Skripte basierend auf der Executor-Einstellung.

Wir werden diese standardisierte Syntax im nächsten Abschnitt behandeln.

### Fazit

Du weißt jetzt, wie du den Executor änderst, um verschiedene Arten von Recheninfrastruktur zu verwenden.

### Wie geht es weiter?

Lerne, wie du Ressourcenzuweisungen und -beschränkungen in Nextflow bewertest und ausdrückst.

---

## 5. Rechenressourcenzuweisungen steuern

??? example "Szenario"

    Deine Pipeline schlägt auf dem Cluster ständig fehl, weil Aufgaben wegen Überschreitung der Speichergrenzen beendet werden.
    Oder vielleicht werden dir Ressourcen in Rechnung gestellt, die du nicht nutzt, und du möchtest die Kosten optimieren.

Die meisten Hochleistungsrechenplattformen erlauben (und erfordern manchmal), dass du bestimmte Ressourcenzuweisungsparameter wie Anzahl der CPUs und Speicher angibst.

Standardmäßig verwendet Nextflow eine einzelne CPU und 2 GB Speicher für jeden Prozess.
Die entsprechenden Prozess-Direktiven heißen `cpus` und `memory`, sodass die folgende Konfiguration impliziert ist:

```groovy title="Eingebaute Konfiguration" linenums="1"
process {
    cpus = 1
    memory = 2.GB
}
```

Du kannst diese Werte ändern, entweder für alle Prozesse oder für bestimmte benannte Prozesse, indem du zusätzliche Prozess-Direktiven in deiner Konfigurationsdatei verwendest.
Nextflow übersetzt sie in die entsprechenden Anweisungen für den gewählten Executor.

Aber woher weißt du, welche Werte du verwenden sollst?

### 5.1. Den Workflow ausführen, um einen Ressourcennutzungsbericht zu generieren

??? example "Szenario"

    Du weißt nicht im Voraus, wie viel Speicher oder CPU deine Prozesse wahrscheinlich benötigen, und möchtest vermeiden, Ressourcen zu verschwenden oder Jobs beendet zu bekommen.

Wenn du nicht im Voraus weißt, wie viel CPU und Speicher deine Prozesse wahrscheinlich benötigen, kannst du ein Ressourcenprofiling durchführen, was bedeutet, dass du den Workflow mit einigen Standardzuweisungen ausführst, aufzeichnest, wie viel jeder Prozess verwendet hat, und von dort aus abschätzt, wie die Basiszuweisungen angepasst werden sollten.

Praktischerweise enthält Nextflow eingebaute Tools dafür und generiert auf Anfrage gerne einen Bericht für dich.

Füge dazu `-with-report <dateiname>.html` zu deiner Kommandozeile hinzu.

```bash
nextflow run 3-main.nf -with-report report-config-1.html
```

Der Bericht ist eine HTML-Datei, die du herunterladen und in deinem Browser öffnen kannst. Du kannst auch mit der rechten Maustaste darauf im Datei-Explorer auf der linken Seite klicken und auf `Show preview` klicken, um ihn in der Trainingsumgebung anzuzeigen.

Nimm dir ein paar Minuten Zeit, um den Bericht durchzusehen und zu sehen, ob du einige Möglichkeiten zur Anpassung von Ressourcen identifizieren kannst.
Stelle sicher, dass du auf die Registerkarten klickst, die die Nutzungsergebnisse als Prozentsatz dessen zeigen, was zugewiesen wurde.

Siehe [Reports](https://nextflow.io/docs/latest/reports.html) für Dokumentation zu allen verfügbaren Funktionen.

### 5.2. Ressourcenzuweisungen für alle Prozesse festlegen

Das Profiling zeigt, dass die Prozesse in unserem Trainings-Workflow sehr leichtgewichtig sind, also lass uns die Standard-Speicherzuweisung auf 1 GB pro Prozess reduzieren.

Füge das Folgende zu deiner `nextflow.config`-Datei hinzu, vor dem Abschnitt Pipeline-Parameter:

```groovy title="nextflow.config" linenums="4"
/*
* Prozesseinstellungen
*/
process {
    memory = 1.GB
}
```

Das wird helfen, die Menge an Rechenleistung zu reduzieren, die wir verbrauchen.

### 5.3. Ressourcenzuweisungen für einen bestimmten Prozess festlegen

Gleichzeitig werden wir so tun, als ob der `cowpy`-Prozess mehr Ressourcen als die anderen benötigt, nur damit wir demonstrieren können, wie man Zuweisungen für einen einzelnen Prozess anpasst.

=== "Danach"

    ```groovy title="nextflow.config" linenums="4" hl_lines="6-9"
    /*
    * Prozesseinstellungen
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
    * Prozesseinstellungen
    */
    process {
        memory = 1.GB
    }
    ```

Mit dieser Konfiguration werden alle Prozesse 1 GB Speicher und eine einzelne CPU (der implizierte Standard) anfordern, außer dem `cowpy`-Prozess, der 2 GB und 2 CPUs anfordern wird.

!!! info

    Wenn du eine Maschine mit wenigen CPUs hast und eine hohe Anzahl pro Prozess zuweist, siehst du möglicherweise, dass Prozessaufrufe hintereinander in die Warteschlange gestellt werden.
    Dies liegt daran, dass Nextflow sicherstellt, dass wir nicht mehr CPUs anfordern, als verfügbar sind.

### 5.4. Den Workflow mit der aktualisierten Konfiguration ausführen

Lass uns das ausprobieren und einen anderen Dateinamen für den Profiling-Bericht angeben, damit wir die Leistung vor und nach den Konfigurationsänderungen vergleichen können.

```bash
nextflow run 3-main.nf -with-report report-config-2.html
```

Du wirst wahrscheinlich keinen wirklichen Unterschied bemerken, da dies eine so kleine Workload ist, aber dies ist der Ansatz, den du verwenden würdest, um die Leistung und Ressourcenanforderungen eines realen Workflows zu analysieren.

Es ist sehr nützlich, wenn deine Prozesse unterschiedliche Ressourcenanforderungen haben. Es ermöglicht dir, die Ressourcenzuweisungen, die du für jeden Prozess einrichtest, basierend auf tatsächlichen Daten richtig zu dimensionieren, nicht auf Vermutungen.

!!! tip

    Dies ist nur ein kleiner Vorgeschmack darauf, was du tun kannst, um deine Ressourcennutzung zu optimieren.
    Nextflow selbst hat einige wirklich nette [dynamische Wiederholungslogik](https://nextflow.io/docs/latest/process.html#dynamic-task-resources) eingebaut, um Jobs, die aufgrund von Ressourcenbeschränkungen fehlschlagen, erneut zu versuchen.
    Darüber hinaus bietet die Seqera Platform KI-gesteuerte Tools zur automatischen Optimierung deiner Ressourcenzuweisungen.

### 5.5. Ressourcengrenzen hinzufügen

Abhängig davon, welchen Computing-Executor und welche Recheninfrastruktur du verwendest, kann es einige Einschränkungen geben, was du zuweisen kannst (oder musst).
Zum Beispiel kann dein Cluster verlangen, dass du innerhalb bestimmter Grenzen bleibst.

Du kannst die `resourceLimits`-Direktive verwenden, um die relevanten Beschränkungen festzulegen. Die Syntax sieht so aus, wenn sie allein in einem Prozessblock steht:

```groovy title="Syntaxbeispiel"
process {
    resourceLimits = [
        memory: 750.GB,
        cpus: 200,
        time: 30.d
    ]
}
```

Nextflow übersetzt diese Werte in die entsprechenden Anweisungen, abhängig vom Executor, den du angegeben hast.

Wir werden dies nicht ausführen, da wir in der Trainingsumgebung keinen Zugriff auf relevante Infrastruktur haben.
Wenn du jedoch versuchen würdest, den Workflow mit Ressourcenzuweisungen auszuführen, die diese Grenzen überschreiten, und dann den `sbatch`-Befehl in der `.command.run`-Skriptdatei nachschlagen würdest, würdest du sehen, dass die Anfragen, die tatsächlich an den Executor gesendet werden, auf die durch `resourceLimits` angegebenen Werte begrenzt sind.

??? info "Institutionelle Referenzkonfigurationen"

    Das nf-core-Projekt hat eine [Sammlung von Konfigurationsdateien](https://nf-co.re/configs/) zusammengestellt, die von verschiedenen Institutionen auf der ganzen Welt geteilt werden und eine breite Palette von HPC- und Cloud-Executors abdecken.

    Diese geteilten Configs sind sowohl für Personen wertvoll, die dort arbeiten und daher einfach die Konfiguration ihrer Institution sofort nutzen können, als auch als Modell für Personen, die eine Konfiguration für ihre eigene Infrastruktur entwickeln möchten.

### Fazit

Du weißt, wie du einen Profiling-Bericht generierst, um die Ressourcennutzung zu bewerten, und wie du Ressourcenzuweisungen für alle Prozesse und/oder für einzelne Prozesse änderst sowie Ressourcenbeschränkungen für die Ausführung auf HPC festlegst.

### Wie geht es weiter?

Lerne, wie du voreingestellte Konfigurationsprofile einrichtest und zur Laufzeit zwischen ihnen wechselst.

---

## 6. Profile verwenden, um zwischen voreingestellten Konfigurationen zu wechseln

??? example "Szenario"

    Du wechselst regelmäßig zwischen der Ausführung von Pipelines auf deinem Laptop für die Entwicklung und auf dem HPC deiner Institution für Produktionsläufe.
    Du bist es leid, jedes Mal manuell Konfigurationseinstellungen zu ändern, wenn du die Umgebung wechselst.

Wir haben dir eine Reihe von Möglichkeiten gezeigt, wie du deine Pipeline-Konfiguration je nach Projekt, an dem du arbeitest, oder der von dir verwendeten Rechenumgebung anpassen kannst.

Du möchtest möglicherweise zwischen alternativen Einstellungen wechseln, je nachdem, welche Recheninfrastruktur du verwendest. Zum Beispiel möchtest du möglicherweise lokal auf deinem Laptop entwickeln und kleine Tests durchführen und dann vollständige Workloads auf HPC oder Cloud ausführen.

Nextflow ermöglicht es dir, eine beliebige Anzahl von [**Profilen**](https://nextflow.io/docs/latest/config.html#profiles) einzurichten, die verschiedene Konfigurationen beschreiben, die du dann zur Laufzeit mit einem Kommandozeilenargument auswählen kannst, anstatt die Konfigurationsdatei selbst ändern zu müssen.

### 6.1. Profile für den Wechsel zwischen lokaler Entwicklung und Ausführung auf HPC erstellen

Lass uns zwei alternative Profile einrichten; eines für die Ausführung kleiner Lasten auf einem normalen Computer, wo wir Docker-Container verwenden, und eines für die Ausführung auf einem Universitäts-HPC mit einem Slurm-Scheduler, wo wir Conda-Pakete verwenden.

#### 6.1.1. Die Profile einrichten

Füge das Folgende zu deiner `nextflow.config`-Datei hinzu, nach dem Abschnitt Pipeline-Parameter, aber vor den Ausgabeeinstellungen:

```groovy title="nextflow.config" linenums="24"
/*
* Profile
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

Du siehst, dass wir für den Universitäts-HPC auch Ressourcenbeschränkungen angeben.

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

Wie du sehen kannst, ermöglicht uns dies, zur Laufzeit sehr bequem zwischen Konfigurationen umzuschalten.

!!! warning

    Das `univ_hpc`-Profil wird in der Trainingsumgebung nicht ordnungsgemäß ausgeführt, da wir keinen Zugriff auf einen Slurm-Scheduler haben.

Wenn wir in Zukunft andere Konfigurationselemente finden, die immer mit diesen zusammen auftreten, können wir sie einfach zum entsprechenden Profil hinzufügen.
Wir können auch zusätzliche Profile erstellen, wenn es andere Konfigurationselemente gibt, die wir zusammen gruppieren möchten.

### 6.2. Ein Profil von Testparametern erstellen

??? example "Szenario"

    Du möchtest, dass andere deine Pipeline schnell ausprobieren können, ohne ihre eigenen Eingabedaten sammeln zu müssen.

Profile sind nicht nur für Infrastrukturkonfiguration.
Wir können sie auch verwenden, um Standardwerte für Workflow-Parameter festzulegen, um es anderen zu erleichtern, den Workflow auszuprobieren, ohne selbst geeignete Eingabewerte sammeln zu müssen.
Du kannst dies als Alternative zur Verwendung einer Parameterdatei betrachten.

#### 6.2.1. Das Profil einrichten

Die Syntax zum Ausdrücken von Standardwerten in diesem Kontext sieht so aus, für ein Profil, das wir `test` nennen:

```groovy title="Syntaxbeispiel"
    test {
        params.<parameter1>
        params.<parameter2>
        ...
    }
```

Wenn wir ein Testprofil für unseren Workflow hinzufügen, wird der `profiles`-Block zu:

```groovy title="nextflow.config" linenums="24"
/*
* Profile
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

Genau wie bei technischen Konfigurationsprofilen kannst du mehrere verschiedene Profile einrichten, die Parameter unter einem beliebigen Namen angeben, den du magst.

#### 6.2.2. Den Workflow lokal mit dem Testprofil ausführen

Praktischerweise schließen sich Profile nicht gegenseitig aus, sodass wir mehrere Profile in unserer Kommandozeile mit der folgenden Syntax angeben können: `-profile <profil1>,<profil2>` (für eine beliebige Anzahl von Profilen).

Wenn du Profile kombinierst, die Werte für dieselben Konfigurationselemente setzen und in derselben Konfigurationsdatei beschrieben sind, löst Nextflow den Konflikt, indem es den Wert verwendet, den es zuletzt gelesen hat (_d.h._ was später in der Datei kommt).
Wenn die widersprüchlichen Einstellungen in verschiedenen Konfigurationsquellen festgelegt sind, gilt die Standard-[Rangfolge](https://www.nextflow.io/docs/latest/config.html).

Lass uns versuchen, das Testprofil zu unserem vorherigen Befehl hinzuzufügen:

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

Dies wird Docker verwenden, wo möglich, und Ausgaben unter `results_config/test` erzeugen, und diesmal ist der Charakter das komödiantische Duo `dragonandcow`.

??? abstract "Dateiinhalt"

    ```console title="results_config/test/"
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

Das bedeutet, dass jeder, solange wir Testdatendateien mit dem Workflow-Code verteilen, den Workflow schnell ausprobieren kann, ohne eigene Eingaben über die Kommandozeile oder eine Parameterdatei bereitstellen zu müssen.

!!! tip

    Wir können auf URLs für größere Dateien verweisen, die extern gespeichert sind.
    Nextflow lädt sie automatisch herunter, solange eine offene Verbindung besteht.

    Für weitere Details siehe die Side Quest [Working with Files](../side_quests/working_with_files.md)

### 6.3. `nextflow config` verwenden, um die aufgelöste Konfiguration zu sehen

Wie oben erwähnt, kann manchmal derselbe Parameter in Profilen, die du kombinieren möchtest, auf unterschiedliche Werte gesetzt werden.
Und allgemeiner gibt es zahlreiche Orte, an denen Konfigurationselemente gespeichert werden können, und manchmal können dieselben Eigenschaften an verschiedenen Orten auf unterschiedliche Werte gesetzt werden.

Nextflow wendet eine festgelegte [Rangfolge](https://nextflow.io/docs/latest/config.html#configuration-file) an, um Konflikte zu lösen, aber das kann schwierig sein, selbst zu bestimmen.
Und selbst wenn nichts in Konflikt steht, kann es mühsam sein, alle möglichen Orte nachzuschlagen, an denen Dinge konfiguriert werden könnten.

Glücklicherweise enthält Nextflow ein praktisches Dienstprogramm namens `config`, das diesen gesamten Prozess für dich automatisieren kann.

Das `config`-Tool durchsucht alle Inhalte in deinem aktuellen Arbeitsverzeichnis, sammelt alle Konfigurationsdateien ein und erzeugt die vollständig aufgelöste Konfiguration, die Nextflow zum Ausführen des Workflows verwenden würde.
Dies ermöglicht es dir, herauszufinden, welche Einstellungen verwendet werden, ohne etwas starten zu müssen.

#### 6.3.1. Die Standardkonfiguration auflösen

Führe diesen Befehl aus, um die Konfiguration aufzulösen, die standardmäßig angewendet würde.

```bash
nextflow config
```

??? success "Befehlsausgabe"

    ```groovy
    params {
      input = 'data/greetings.csv'
      batch = 'batch'
      character = 'turkey'
    }

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

    outputDir = 'results_config/batch'

    workflow {
      output {
          mode = 'copy'
      }
    }
    ```

Dies zeigt dir die Basiskonfiguration, die du erhältst, wenn du nichts Zusätzliches in der Kommandozeile angibst.

#### 6.3.2. Die Konfiguration mit aktivierten spezifischen Einstellungen auflösen

Wenn du Kommandozeilenparameter bereitstellst, z.B. ein oder mehrere Profile aktivierst oder eine Parameterdatei lädst, berücksichtigt der Befehl diese zusätzlich.

```bash
nextflow config -profile my_laptop,test
```

??? success "Befehlsausgabe"

    ```groovy
    params {
      input = 'data/greetings.csv'
      batch = 'test'
      character = 'dragonandcow'
    }

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

    outputDir = 'results_config/test'

    workflow {
      output {
          mode = 'copy'
      }
    }
    ```

Dies wird besonders nützlich für komplexe Projekte, die mehrere Konfigurationsebenen beinhalten.

### Fazit

Du weißt, wie du Profile verwendest, um zur Laufzeit mit minimalem Aufwand eine voreingestellte Konfiguration auszuwählen.
Allgemeiner weißt du, wie du deine Workflow-Ausführungen konfigurierst, um verschiedenen Rechenplattformen zu entsprechen und die Reproduzierbarkeit deiner Analysen zu verbessern.

### Wie geht es weiter?

Lerne, wie du Pipelines direkt aus Remote-Repositories wie GitHub ausführst.

---

## 7. Pipelines aus Remote-Repositories ausführen

??? example "Szenario"

    Du möchtest eine etablierte Pipeline wie die von nf-core ausführen, ohne den Code selbst herunterladen und verwalten zu müssen.

Bisher haben wir Workflow-Skripte ausgeführt, die sich im aktuellen Verzeichnis befinden.
In der Praxis möchtest du oft Pipelines ausführen, die in Remote-Repositories wie GitHub gespeichert sind.

Nextflow macht dies unkompliziert: Du kannst jede Pipeline direkt von einer Git-Repository-URL ausführen, ohne sie zuerst manuell herunterladen zu müssen.

### 7.1. Eine Pipeline von GitHub ausführen

Die grundlegende Syntax zum Ausführen einer Remote-Pipeline ist `nextflow run <repository>`, wobei `<repository>` ein GitHub-Repository-Pfad wie `nextflow-io/hello`, eine vollständige URL oder ein Pfad zu GitLab, Bitbucket oder anderen Git-Hosting-Diensten sein kann.

Versuche, die offizielle Nextflow-"hello"-Demo-Pipeline auszuführen:

```bash
nextflow run nextflow-io/hello
```

??? success "Befehlsausgabe"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Pulling nextflow-io/hello ...
     downloaded from https://github.com/nextflow-io/hello.git
    Launching `https://github.com/nextflow-io/hello` [sleepy_swanson] DSL2 - revision: 2ce0b0e294 [master]

    executor >  local (4)
    [ba/08236d] sayHello (4) [100%] 4 of 4 ✔
    Ciao world!

    Hello world!

    Bonjour world!

    Hola world!
    ```

Wenn du zum ersten Mal eine Remote-Pipeline ausführst, lädt Nextflow sie herunter und speichert sie lokal zwischen.
Nachfolgende Läufe verwenden die zwischengespeicherte Version, es sei denn, du forderst explizit ein Update an.

### 7.2. Eine Version für Reproduzierbarkeit angeben

Standardmäßig führt Nextflow die neueste Version vom Standard-Branch aus.
Du kannst eine bestimmte Version (Tag), einen Branch oder einen Commit mit dem `-r`-Flag angeben:

```bash
nextflow run nextflow-io/hello -r v1.3
```

??? success "Befehlsausgabe"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `https://github.com/nextflow-io/hello` [sick_carson] DSL2 - revision: 2ce0b0e294 [v1.3]

    executor >  local (4)
    [61/e11f77] sayHello (4) [100%] 4 of 4 ✔
    Ciao world!

    Bonjour world!

    Hello world!

    Hola world!
    ```

Die Angabe genauer Versionen ist für die Reproduzierbarkeit unerlässlich.

### Fazit

Du weißt, wie du Pipelines direkt von GitHub und anderen Remote-Repositories ausführst und wie du Versionen für Reproduzierbarkeit angibst.

### Wie geht es weiter?

Klopf dir selbst auf die Schulter!
Du weißt alles, was du wissen musst, um mit dem Ausführen und Verwalten von Nextflow-Pipelines zu beginnen.

Das schließt diesen Kurs ab, aber wenn du weiter lernen möchtest, haben wir zwei Hauptempfehlungen:

- Wenn du tiefer in die Entwicklung deiner eigenen Pipelines eintauchen möchtest, schau dir [Hello Nextflow](../hello_nextflow/index.md) an, einen Kurs für Anfänger\*innen, der dieselbe allgemeine Progression wie dieser abdeckt, aber viel detaillierter auf Kanäle und Operatoren eingeht.
- Wenn du weiter lernen möchtest, wie man Nextflow-Pipelines ausführt, ohne tiefer in den Code einzusteigen, schau dir den ersten Teil von [Hello nf-core](../hello_nf-core/index.md) an, der die Tools zum Finden und Ausführen von Pipelines aus dem äußerst beliebten [nf-core](https://nf-co.re/)-Projekt vorstellt.

Viel Spaß!

---

## Quiz

<quiz>
Wenn Parameterwerte sowohl in der Workflow-Datei als auch in `nextflow.config` gesetzt sind, welcher hat Vorrang?
- [ ] Der Wert in der Workflow-Datei
- [x] Der Wert in der Konfigurationsdatei
- [ ] Der zuerst angetroffene Wert
- [ ] Es verursacht einen Fehler

Mehr erfahren: [1.1. Werte in `nextflow.config` einrichten](#11-werte-in-nextflowconfig-einrichten)
</quiz>

<quiz>
Was ist der Syntaxunterschied zwischen dem Setzen eines Parameter-Standardwerts in einer Workflow-Datei vs. einer Config-Datei?
- [ ] Sie verwenden dieselbe Syntax
- [x] Workflow verwendet typisierte Deklaration (`#!groovy param: Type = value`), Config verwendet Zuweisung (`#!groovy param = value`)
- [ ] Config verwendet typisierte Deklaration, Workflow verwendet Zuweisung
- [ ] Nur Config-Dateien können Standardwerte setzen

Mehr erfahren: [1.1. Werte in `nextflow.config` einrichten](#11-werte-in-nextflowconfig-einrichten)
</quiz>

<quiz>
Wie gibst du eine Parameterdatei beim Ausführen eines Workflows an?
- [ ] `--params params.yaml`
- [ ] `-config params.yaml`
- [x] `-params-file params.yaml`
- [ ] `--input-params params.yaml`

Mehr erfahren: [1.3. Eine Parameterdatei verwenden](#13-eine-parameterdatei-verwenden)
</quiz>

<quiz>
Was steuert die `outputDir`-Konfigurationsoption?
- [ ] Den Speicherort des Arbeitsverzeichnisses
- [x] Den Basispfad, wo Workflow-Ausgaben veröffentlicht werden
- [ ] Das Verzeichnis für Log-Dateien
- [ ] Den Speicherort von Moduldateien

Mehr erfahren: [2.1. Den `outputDir`-Verzeichnisnamen anpassen](#21-den-outputdir-verzeichnisnamen-anpassen)
</quiz>

<quiz>
Wie referenzierst du einen Prozessnamen dynamisch in der Ausgabepfadkonfiguration?
- [ ] `#!groovy ${processName}`
- [ ] `#!groovy path "<process>.name"`
- [x] `#!groovy path { <process>.name }`
- [ ] `@processName`

Mehr erfahren: [2.2. Ausgaben nach Prozess organisieren](#22-ausgaben-nach-prozess-organisieren)
</quiz>

<quiz>
Wenn sowohl Docker als auch Conda aktiviert sind und ein Prozess beide Direktiven hat, welche wird priorisiert?
- [x] Docker (Container)
- [ ] Conda
- [ ] Die zuerst im Prozess definierte
- [ ] Es verursacht einen Fehler

Mehr erfahren: [3. Eine Software-Packaging-Technologie auswählen](#3-eine-software-packaging-technologie-auswahlen)
</quiz>

<quiz>
Was ist der Standard-Executor in Nextflow?
- [x] `local`
- [ ] `slurm`
- [ ] `kubernetes`
- [ ] `aws`

Mehr erfahren: [4. Eine Ausführungsplattform auswählen](#4-eine-ausfuhrungsplattform-auswahlen)
</quiz>

<quiz>
Welcher Befehl generiert einen Ressourcennutzungsbericht?
- [ ] `nextflow run workflow.nf -with-metrics`
- [ ] `nextflow run workflow.nf -with-stats`
- [x] `nextflow run workflow.nf -with-report report.html`
- [ ] `nextflow run workflow.nf -profile report`

Mehr erfahren: [5.1. Den Workflow ausführen, um einen Ressourcennutzungsbericht zu generieren](#51-den-workflow-ausfuhren-um-einen-ressourcennutzungsbericht-zu-generieren)
</quiz>

<quiz>
Wie setzt du Ressourcenanforderungen für einen bestimmten Prozess namens `cowpy` in der Config-Datei?
- [ ] `#!groovy cowpy.memory = '2.GB'`
- [ ] `#!groovy process.cowpy.memory = '2.GB'`
- [x] `#!groovy process { withName: 'cowpy' { memory = '2.GB' } }`
- [ ] `#!groovy resources.cowpy.memory = '2.GB'`

Mehr erfahren: [5.3. Ressourcenzuweisungen für einen bestimmten Prozess festlegen](#53-ressourcenzuweisungen-fur-einen-bestimmten-prozess-festlegen)
</quiz>

<quiz>
Was macht die `resourceLimits`-Direktive?
- [ ] Setzt Mindestressourcenanforderungen
- [ ] Weist Prozessen Ressourcen zu
- [x] Begrenzt die maximalen Ressourcen, die angefordert werden können
- [ ] Überwacht die Ressourcennutzung in Echtzeit

Mehr erfahren: [5.5. Ressourcengrenzen hinzufügen](#55-ressourcengrenzen-hinzufugen)
</quiz>

<quiz>
Wie gibst du mehrere Profile in einem einzigen Befehl an?
- [ ] `-profile profile1 -profile profile2`
- [ ] `-profiles profile1,profile2`
- [x] `-profile profile1,profile2`
- [ ] `--profile profile1 --profile profile2`

Mehr erfahren: [6. Profile verwenden, um zwischen voreingestellten Konfigurationen zu wechseln](#6-profile-verwenden-um-zwischen-voreingestellten-konfigurationen-zu-wechseln)
</quiz>

<quiz>
Welcher Befehl zeigt die vollständig aufgelöste Konfiguration, die Nextflow verwenden würde?
- [ ] `nextflow show-config`
- [ ] `nextflow settings`
- [x] `nextflow config`
- [ ] `nextflow resolve`

Mehr erfahren: [6.3. `nextflow config` verwenden, um die aufgelöste Konfiguration zu sehen](#63-nextflow-config-verwenden-um-die-aufgeloste-konfiguration-zu-sehen)
</quiz>

<quiz>
Wofür können Profile verwendet werden? (Wähle alle zutreffenden aus)
- [x] Definieren infrastrukturspezifischer Einstellungen (Executors, Container)
- [x] Festlegen von Ressourcengrenzen für verschiedene Umgebungen
- [x] Bereitstellen von Testparametern für einfaches Workflow-Testen
- [ ] Definieren neuer Prozesse

Mehr erfahren: [6. Profile verwenden, um zwischen voreingestellten Konfigurationen zu wechseln](#6-profile-verwenden-um-zwischen-voreingestellten-konfigurationen-zu-wechseln)
</quiz>
