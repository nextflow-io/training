# Teil 6: Hello Config

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } KI-gestützte Übersetzung - [mehr erfahren & Verbesserungen vorschlagen](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/FcZTiE25TeA?si=y8lAedhEHWaTV4zd&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

/// caption
:fontawesome-brands-youtube:{ .youtube } Siehe [die gesamte Playlist](https://youtube.com/playlist?list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&si=eF7cLR62goy-lc6n) auf dem Nextflow YouTube-Kanal.

:green_book: Das Video-Transkript ist [hier](./transcripts/06_hello_config.md) verfügbar.
///

In diesem Abschnitt erfährst du, wie du die Konfiguration deiner Nextflow-Pipeline einrichtest und verwaltest, sodass du ihr Verhalten anpassen, sie an verschiedene Umgebungen anpassen und die Ressourcennutzung optimieren kannst, _ohne eine einzige Zeile des Workflow-Codes selbst zu ändern_.

Es gibt mehrere Möglichkeiten, dies zu tun, die kombiniert werden können und gemäß der in der Konfigurationsdokumentation beschriebenen [Rangfolge](https://nextflow.io/docs/latest/config.html) interpretiert werden.

In diesem Teil des Kurses zeigen wir dir den einfachsten und gebräuchlichsten Mechanismus für Konfigurationsdateien, die [`nextflow.config`](https://nextflow.io/docs/latest/config.html)-Datei, die du bereits in Teil 5: Hello Containers kennengelernt hast.

Wir werden wesentliche Komponenten der Nextflow-Konfiguration wie Prozess-Direktiven, Executors, Profile und Parameterdateien durchgehen.
Indem du lernst, diese Konfigurationsoptionen effektiv zu nutzen, kannst du die Flexibilität, Skalierbarkeit und Leistung deiner Pipelines verbessern.

??? info "Wie du von diesem Abschnitt aus beginnst"

    Dieser Abschnitt des Kurses setzt voraus, dass du die Teile 1-5 des [Hello Nextflow](./index.md)-Kurses abgeschlossen hast und eine vollständige funktionierende Pipeline besitzt.

    Wenn du den Kurs von diesem Punkt aus beginnst, musst du das `modules`-Verzeichnis und die `nextflow.config`-Datei aus den Lösungen kopieren:

    ```bash
    cp -r solutions/5-hello-containers/modules .
    cp solutions/5-hello-containers/nextflow.config .
    ```

    Die `nextflow.config`-Datei enthält die Zeile `docker.enabled = true`, die die Verwendung von Docker-Containern aktiviert.

    Wenn du mit der Hello-Pipeline nicht vertraut bist oder eine Auffrischung gebrauchen könntest, siehe [diese Infoseite](../info/hello_pipeline.md).

---

## 0. Aufwärmen: Führe `hello-config.nf` aus

Wir werden das Workflow-Skript `hello-config.nf` als Ausgangspunkt verwenden.
Es entspricht dem Skript, das durch die Bearbeitung von Teil 5 dieses Trainingskurses erstellt wurde, außer dass wir die Ausgabeziele geändert haben:

```groovy title="hello-config.nf" linenums="37" hl_lines="3 7 11 15"
output {
    first_output {
        path 'hello_config/intermediates'
        mode 'copy'
    }
    uppercased {
        path 'hello_config/intermediates'
        mode 'copy'
    }
    collected {
        path 'hello_config/intermediates'
        mode 'copy'
    }
    batch_report {
        path 'hello_config'
        mode 'copy'
    }
    cowpy_art {
        path 'hello_config'
        mode 'copy'
    }
}
```

Um sicherzustellen, dass alles funktioniert, führe das Skript einmal aus, bevor du Änderungen vornimmst:

```bash
nextflow run hello-config.nf
```

??? success "Befehlsausgabe"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-config.nf` [nice_escher] DSL2 - revision: d5dfdc9872

    executor > local (7)
    [6a/bc46a6] sayHello (2) [100%] 3 of 3 ✔
    [33/67bc48] convertToUpper (3) [100%] 3 of 3 ✔
    [b5/de03ba] collectGreetings [100%] 1 of 1 ✔
    [98/c6b57b] cowpy              | 1 of 1 ✔
    ```

Wie zuvor findest du die Ausgabedateien im Verzeichnis, das im `output`-Block angegeben ist (`results/hello_config/`).

??? abstract "Verzeichnisinhalt"

    ```console
    results/hello_config/
    ├── cowpy-COLLECTED-batch-output.txt
    ├── intermediates
    │   ├── Bonjour-output.txt
    │   ├── COLLECTED-batch-output.txt
    │   ├── Hello-output.txt
    │   ├── Holà-output.txt
    │   ├── UPPER-Bonjour-output.txt
    │   ├── UPPER-Hello-output.txt
    │   └── UPPER-Holà-output.txt
    └── batch-report.txt
    ```

Die finale ASCII-Art-Ausgabe befindet sich im Verzeichnis `results/hello_config/` unter dem Namen `cowpy-COLLECTED-batch-output.txt`.

??? abstract "Dateiinhalt"

    ```console title="results/hello_config/cowpy-COLLECTED-batch-output.txt"
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

Wenn das für dich funktioniert hat, bist du bereit zu lernen, wie du deine Pipelines konfigurierst.

---

## 1. Workflow-Eingabeparameter verwalten

Wir beginnen mit einem Aspekt der Konfiguration, der einfach eine Erweiterung dessen ist, womit wir bisher gearbeitet haben: die Verwaltung von Eingabeparametern.

Derzeit ist unser Workflow so eingerichtet, dass er mehrere Parameterwerte über die Befehlszeile akzeptiert, mit Standardwerten, die in einem `params`-Block im Workflow-Skript selbst festgelegt sind.
Möglicherweise möchtest du diese Standardwerte jedoch überschreiben, ohne entweder Parameter in der Befehlszeile angeben oder die ursprüngliche Skriptdatei ändern zu müssen.

Es gibt mehrere Möglichkeiten, dies zu tun; wir zeigen dir drei grundlegende Methoden, die sehr häufig verwendet werden.

### 1.1. Standardwerte in `nextflow.config` verschieben

Dies ist der einfachste Ansatz, obwohl er möglicherweise der am wenigsten flexible ist, da die Haupt-`nextflow.config`-Datei nichts ist, was du für jeden Lauf bearbeiten möchtest.
Aber es hat den Vorteil, die Belange der _Deklaration_ der Parameter im Workflow (die definitiv dorthin gehört) von der Bereitstellung von _Standardwerten_ zu trennen, die eher in einer Konfigurationsdatei zu Hause sind.

Lass uns dies in zwei Schritten tun.

#### 1.1.1. Erstelle einen `params`-Block in der Konfigurationsdatei

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
Die Syntax ist etwas anders.
In der Workflow-Datei sind das typisierte Deklarationen.
In der Konfiguration sind das Wertzuweisungen.

Technisch gesehen reicht dies aus, um die Standardwerte zu überschreiben, die noch in der Workflow-Datei angegeben sind.
Du könntest das Zeichen ändern, zum Beispiel, und den Workflow ausführen, um dich zu überzeugen, dass der in der Konfigurationsdatei festgelegte Wert den in der Workflow-Datei festgelegten überschreibt.

Aber im Sinne der vollständigen Verlagerung der Konfiguration in die Konfigurationsdatei, lass uns diese Werte vollständig aus der Workflow-Datei entfernen.

#### 1.1.2. Entferne die Werte aus dem `params`-Block in der Workflow-Datei

Nimm die folgenden Codeänderungen in der `hello-config.nf`-Workflow-Datei vor:

=== "Danach"

    ```groovy title="hello-config.nf" linenums="9" hl_lines="5-7"
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

    ```groovy title="hello-config.nf" linenums="9" hl_lines="5-7"
    /*
    * Pipeline-Parameter
    */
    params {
        input: Path = 'data/greetings.csv'
        batch: String = 'batch'
        character: String = 'turkey'
    }
    ```

Jetzt setzt die Workflow-Datei selbst keine Standardwerte für diese Parameter.

#### 1.1.3. Führe die Pipeline aus

Lass uns testen, dass es korrekt funktioniert.

```bash
nextflow run hello-config.nf
```

??? success "Befehlsausgabe"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-config.nf` [disturbed_einstein] DSL2 - revision: ede9037d02

    executor >  local (8)
    [f0/35723c] sayHello (2)       | 3 of 3 ✔
    [40/3efd1a] convertToUpper (3) | 3 of 3 ✔
    [17/e97d32] collectGreetings   | 1 of 1 ✔
    [98/c6b57b] cowpy              | 1 of 1 ✔
    ```

Dies erzeugt immer noch die gleiche Ausgabe wie zuvor.

Die finale ASCII-Art-Ausgabe befindet sich im Verzeichnis `results/hello_config/` unter dem Namen `cowpy-COLLECTED-batch-output.txt`, wie zuvor.

??? abstract "Dateiinhalt"

    ```console title="results/hello_config/cowpy-COLLECTED-batch-output.txt"
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

Funktional hat diese Verschiebung nichts geändert, aber konzeptionell ist es etwas sauberer, die Standardwerte in der Konfigurationsdatei festzulegen.

### 1.2. Verwende eine laufspezifische Konfigurationsdatei

Das ist großartig, aber manchmal möchtest du vielleicht einige temporäre Experimente mit verschiedenen Standardwerten durchführen, ohne die Haupt-Konfigurationsdatei zu verändern.
Du kannst dies tun, indem du eine neue `nextflow.config`-Datei in einem Unterverzeichnis erstellst, das du als Arbeitsverzeichnis für deine Experimente verwenden wirst.

#### 1.2.1. Erstelle das Arbeitsverzeichnis mit einer leeren Konfiguration

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

#### 1.2.2. Richte die experimentelle Konfiguration ein

Öffne nun die neue Datei und füge die Parameter hinzu, die du anpassen möchtest:

```groovy title="tux-run/nextflow.config" linenums="1"
params {
    input = '../data/greetings.csv'
    batch = 'experiment'
    character = 'tux'
}
```

Beachte, dass der Pfad zur Eingabedatei die Verzeichnisstruktur widerspiegeln muss.

#### 1.2.3. Führe die Pipeline aus

Wir können nun unsere Pipeline aus unserem neuen Arbeitsverzeichnis heraus ausführen.
Stelle sicher, den Pfad entsprechend anzupassen!

```bash
nextflow run ../hello-config.nf
```

??? success "Befehlsausgabe"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `../hello-config.nf` [trusting_escher] DSL2 - revision: 356df0818d

    executor >  local (8)
    [59/b66913] sayHello (2)       [100%] 3 of 3 ✔
    [ad/f06364] convertToUpper (3) [100%] 3 of 3 ✔
    [10/714895] collectGreetings   [100%] 1 of 1 ✔
    [88/3ece98] cowpy              [100%] 1 of 1 ✔
    ```

Dies erstellt einen neuen Satz von Verzeichnissen unter `tux-run/`, einschließlich `tux-run/work/` und `tux-run/results/`.

In diesem Lauf kombiniert Nextflow die `nextflow.config` in unserem aktuellen Verzeichnis mit der `nextflow.config` im Stammverzeichnis der Pipeline und überschreibt dadurch das Standardzeichen (turkey) mit dem tux-Zeichen.

Die finale Ausgabedatei sollte das tux-Zeichen enthalten, das die Begrüßungen sagt.

??? abstract "Dateiinhalt"

    ```console title="tux-run/results/hello_config/cowpy-COLLECTED-experiment-output.txt"
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

### 1.3. Verwende eine Parameterdatei

Der Unterverzeichnis-Ansatz funktioniert hervorragend zum Experimentieren, erfordert aber etwas Einrichtung und verlangt, dass du Pfade entsprechend anpasst.
Es gibt einen einfacheren Ansatz, wenn du deine Pipeline mit einem bestimmten Satz von Werten ausführen möchtest oder es jemand anderem mit minimalem Aufwand ermöglichen möchtest.

Nextflow ermöglicht es uns, Parameter über eine [Parameterdatei](https://nextflow.io/docs/latest/config.html#params-file) im YAML- oder JSON-Format anzugeben, was es sehr bequem macht, alternative Sätze von Standardwerten zu verwalten und zu verteilen, zum Beispiel, sowie laufspezifische Parameterwerte.

#### 1.3.1. Untersuche die Beispiel-Parameterdatei

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

#### 1.3.2. Führe die Pipeline aus

Um den Workflow mit dieser Parameterdatei auszuführen, füge einfach `-params-file <dateiname>` zum Basisbefehl hinzu.

```bash
nextflow run hello-config.nf -params-file test-params.yaml
```

??? success "Befehlsausgabe"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-config.nf` [disturbed_sammet] DSL2 - revision: ede9037d02

    executor >  local (8)
    [f0/35723c] sayHello (2)       | 3 of 3 ✔
    [40/3efd1a] convertToUpper (3) | 3 of 3 ✔
    [17/e97d32] collectGreetings   | 1 of 1 ✔
    [98/c6b57b] cowpy              | 1 of 1 ✔
    ```

Die finale Ausgabedatei sollte das Stegosaurus-Zeichen enthalten, das die Begrüßungen sagt.

??? abstract "Dateiinhalt"

    ```console title="results/hello_config/cowpy-COLLECTED-yaml-output.txt"
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
In diesen Fällen ermöglicht uns die Verwendung einer Parameterdatei, Parameterwerte zur Laufzeit bereitzustellen, ohne massive Befehlszeilen eingeben und ohne das Workflow-Skript ändern zu müssen.

Es erleichtert auch die Verteilung von Parametersätzen an Mitarbeiter\*innen oder als unterstützende Informationen für eine Veröffentlichung, zum Beispiel.
Dies macht deine Arbeit für andere reproduzierbarer.

### Fazit

Du weißt, wie du wichtige Konfigurationsoptionen zur Verwaltung von Workflow-Eingaben nutzen kannst.

### Wie geht es weiter?

Lerne, wie du verwaltest, wo und wie deine Workflow-Ausgaben veröffentlicht werden.

---

## 2. Workflow-Ausgaben verwalten

Bisher haben wir alle Pfade für Workflow-Level-Ausgabedeklarationen fest codiert, und wie wir festgestellt haben, als wir begannen, mehrere Ausgaben hinzuzufügen, kann dabei etwas Wiederholung auftreten.

Schauen wir uns einige gängige Möglichkeiten an, wie du dies flexibler konfigurieren kannst.

### 2.1. Passe das Ausgabeverzeichnis mit `-output-dir` an

Wenn wir kontrollieren, wie unsere 'veröffentlichten' Ausgaben organisiert werden, haben wir zwei unterschiedliche Prioritäten:

- Das oberste Ausgabeverzeichnis
- Wie Dateien innerhalb dieses Verzeichnisses organisiert werden

Wir haben bisher das Standard-Oberverzeichnis verwendet: `results`.
Lass uns damit beginnen, dies anzupassen, indem wir die CLI-Option `-output-dir` verwenden.

#### 2.1.1. Führe die Pipeline mit `-output-dir` aus

Die Option `-output-dir` (Kurzform: `-o`) überschreibt das Standard-Ausgabeverzeichnis (`results/`) für alle Workflow-Ausgaben.
Dies ist die empfohlene Methode, um den Stammpfad zu steuern, wo Ausgaben veröffentlicht werden.

```bash
nextflow run hello-config.nf -output-dir custom-outdir-cli/
```

??? success "Befehlsausgabe"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-config.nf` [prickly_kay] DSL2 - revision: 32ecc4fba2

    executor >  local (8)
    [9f/332636] sayHello (1)       [100%] 3 of 3 ✔
    [03/a55991] convertToUpper (3) [100%] 3 of 3 ✔
    [e5/ab7893] collectGreetings   [100%] 1 of 1 ✔
    [a8/97338e] cowpy              [100%] 1 of 1 ✔
    ```

Dies veröffentlicht Ausgaben nach `custom-outdir-cli/` anstelle von `results/`:

??? abstract "Verzeichnisinhalt"

    ```console
    custom-outdir-cli/
    └── hello_config
        ├── batch-report.txt
        ├── cowpy-COLLECTED-batch-output.txt
        └── intermediates
            ├── Bonjour-output.txt
            ├── COLLECTED-batch-output.txt
            ├── Hello-output.txt
            ├── Holà-output.txt
            ├── UPPER-Bonjour-output.txt
            ├── UPPER-Hello-output.txt
            └── UPPER-Holà-output.txt
    ```

Beachte, dass wir immer noch das `hello_config`-Unterverzeichnis aus den `path`-Deklarationen im output-Block haben.
Lass uns das aufräumen.

#### 2.1.2. Entferne fest codierte Pfade aus dem output-Block

Das `hello_config/`-Präfix wurde in früheren Kapiteln fest codiert, aber da wir jetzt lernen, Ausgabepfade flexibel zu konfigurieren, können wir diese Fest-Codierung entfernen.
Für Ausgaben, die kein Unterverzeichnis benötigen, können wir die `path`-Direktive auf einen leeren String setzen oder sie vollständig entfernen.

Nimm die folgenden Codeänderungen in der Workflow-Datei vor:

=== "Danach"

    ```groovy title="hello-config.nf" linenums="42" hl_lines="3 7 11 15 19"
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

    ```groovy title="hello-config.nf" linenums="42" hl_lines="3 7 11 15 19"
    output {
        first_output {
            path 'hello_config/intermediates'
            mode 'copy'
        }
        uppercased {
            path 'hello_config/intermediates'
            mode 'copy'
        }
        collected {
            path 'hello_config/intermediates'
            mode 'copy'
        }
        batch_report {
            path 'hello_config'
            mode 'copy'
        }
        cowpy_art {
            path 'hello_config'
            mode 'copy'
        }
    }
    ```

Führe die Pipeline erneut aus:

```bash
nextflow run hello-config.nf -output-dir custom-outdir-cli-2/
```

Jetzt werden die Ausgaben direkt unter `custom-outdir-cli-2/` veröffentlicht, ohne das `hello_config`-Unterverzeichnis:

??? abstract "Verzeichnisinhalt"

    ```console
    custom-outdir-cli-2/
    ├── batch-report.txt
    ├── cowpy-COLLECTED-batch-output.txt
    └── intermediates
        ├── Bonjour-output.txt
        ├── COLLECTED-batch-output.txt
        ├── Hello-output.txt
        ├── Holà-output.txt
        ├── UPPER-Bonjour-output.txt
        ├── UPPER-Hello-output.txt
        └── UPPER-Holà-output.txt
    ```

!!! tip

    Die Option `-output-dir` wird verwendet, um zu steuern, _wohin_ Ausgaben gehen, während die `path`-Direktive im output-Block die _Unterverzeichnisstruktur_ steuert.

### 2.2. Dynamische Ausgabepfade

Zusätzlich zur Änderung des Ausgabeverzeichnisses über die CLI können wir auch einen benutzerdefinierten Standardwert in der Config-Datei mit `outputDir` setzen.
Dies ermöglicht es uns, den Verzeichnispfad dynamisch festzulegen - nicht nur mit statischen Strings.

#### 2.2.1. Setze `outputDir` in der Konfigurationsdatei

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
    outputDir = "custom-outdir-config/${params.batch}"
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

Dies setzt das Ausgabeverzeichnis auf `custom-outdir-config/` plus den Wert des `batch`-Parameters als Unterverzeichnis.
Jetzt kannst du den Ausgabeort ändern, indem du den `--batch`-Parameter setzt:

```bash
nextflow run hello-config.nf --batch my_run
```

Dies veröffentlicht Ausgaben nach `custom-outdir-config/my_run/`.

!!! note

    Die CLI-Option `-output-dir` hat Vorrang vor der `outputDir`-Konfigurationseinstellung.
    Wenn sie gesetzt ist, wird die Config-Option vollständig ignoriert.

#### 2.2.2. Unterverzeichnisse mit Batch- und Prozessnamen

Wir können auch Unterverzeichnis-Ausgabe-`path`-Deklarationen dynamisch setzen, auf Basis pro Ausgabe.

Zum Beispiel können wir unsere Ausgaben nach Prozess organisieren, indem wir `<process>.name` in der Ausgabepfad-Deklaration referenzieren:

=== "Danach"

    ```groovy title="hello-config.nf" linenums="42" hl_lines="3 7 11 15 19"
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

    ```groovy title="hello-config.nf" linenums="42" hl_lines="3 7 11 15 19"
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

Wir können noch weiter gehen und komplexere Unterverzeichnispfade zusammensetzen.

In der obigen Bearbeitung haben wir die Unterscheidung zwischen `intermediates` versus finalen Ausgaben auf der obersten Ebene gelöscht.
Lass uns das zurückbringen und auch die Dateien in ein `params.batch`-Unterverzeichnis legen.

!!! tip

    Das Einbeziehen von `params.batch` im output-Block `path`, anstatt im `outputDir`-Config, bedeutet, dass es nicht mit `-output-dir` in der CLI überschrieben wird.

Aktualisiere zunächst die Config-Datei, um `${params.batch}` aus `outputDir` zu entfernen (da wir es zu den path-Deklarationen verschieben):

=== "Danach"

    ```groovy title="nextflow.config" linenums="12" hl_lines="4"
    /*
    * Ausgabeeinstellungen
    */
    outputDir = "custom-outdir-config/"
    ```

=== "Vorher"

    ```groovy title="nextflow.config" linenums="12" hl_lines="4"
    /*
    * Ausgabeeinstellungen
    */
    outputDir = "custom-outdir-config/${params.batch}"
    ```

Nimm dann die folgenden Änderungen in der Workflow-Datei vor:

=== "Danach"

    ```groovy title="hello-config.nf" linenums="42" hl_lines="3 7 11 15 19"
    output {
        first_output {
            path { "${params.batch}/intermediates/${sayHello.name}" }
            mode 'copy'
        }
        uppercased {
            path { "${params.batch}/intermediates/${convertToUpper.name}" }
            mode 'copy'
        }
        collected {
            path { "${params.batch}/intermediates/${collectGreetings.name}" }
            mode 'copy'
        }
        batch_report {
            path { "${params.batch}/${collectGreetings.name}" }
            mode 'copy'
        }
        cowpy_art {
            path { "${params.batch}/${cowpy.name}" }
            mode 'copy'
        }
    }
    ```

=== "Vorher"

    ```groovy title="hello-config.nf" linenums="42" hl_lines="3 7 11 15 19"
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

#### 2.2.3. Führe die Pipeline aus

Schauen wir uns an, wie das in der Praxis funktioniert, indem wir sowohl `-output-dir` (oder kurz `-o`) auf `custom-outdir-config-2` als auch den Batch-Namen auf `rep2` von der Befehlszeile aus setzen:

```bash
nextflow run hello-config.nf -output-dir custom-outdir-config-2 --batch rep2
```

??? success "Befehlsausgabe"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-config.nf` [mad_curry] DSL2 - revision: 668a98ccb9

    executor >  local (8)
    [9e/6095e0] sayHello (1)       [100%] 3 of 3 ✔
    [05/454d52] convertToUpper (3) [100%] 3 of 3 ✔
    [ed/e3ddfb] collectGreetings   [100%] 1 of 1 ✔
    [39/5e063a] cowpy              [100%] 1 of 1 ✔
    ```

Dies veröffentlicht Ausgaben nach `custom-outdir-config-2/rep2/`, mit dem angegebenen Basispfad _und_ dem Batch-Namen-Unterverzeichnis _und_ Ergebnissen gruppiert nach Prozess:

??? abstract "Verzeichnisinhalt"

    ```console
    custom-outdir-config-2
    └── rep2
        ├── collectGreetings
        │   └── rep2-report.txt
        ├── cowpy
        │   └── cowpy-COLLECTED-rep2-output.txt
        └── intermediates
            ├── collectGreetings
            │   └── COLLECTED-rep2-output.txt
            ├── convertToUpper
            │   ├── UPPER-Bonjour-output.txt
            │   ├── UPPER-Hello-output.txt
            │   └── UPPER-Holà-output.txt
            └── sayHello
                ├── Bonjour-output.txt
                ├── Hello-output.txt
                └── Holà-output.txt
    ```

### 2.3. Setze den Veröffentlichungsmodus auf Workflow-Ebene

Schließlich können wir im Sinne der Reduzierung der Menge an sich wiederholendem Code die pro-Ausgabe-`mode`-Deklarationen durch eine einzige Zeile in der Konfiguration ersetzen.

#### 2.3.1. Füge `workflow.output.mode` zur Konfigurationsdatei hinzu

Füge den folgenden Code zur `nextflow.config`-Datei hinzu:

=== "Danach"

    ```groovy title="nextflow.config" linenums="12" hl_lines="5"
    /*
    * Ausgabeeinstellungen
    */
    outputDir = "custom-outdir-config/"
    workflow.output.mode = 'copy'
    ```

=== "Vorher"

    ```groovy title="nextflow.config" linenums="12"
    /*
    * Ausgabeeinstellungen
    */
    outputDir = "custom-outdir-config/"
    ```

Das Setzen von `workflow.output.mode` in der Konfigurationsdatei reicht aus, um zu überschreiben, was in der Workflow-Datei gesetzt ist, aber lass uns den unnötigen Code trotzdem entfernen.

#### 2.3.2. Entferne den Ausgabemodus aus der Workflow-Datei

Nimm die folgenden Änderungen in der Workflow-Datei vor:

=== "Danach"

    ```groovy title="hello-config.nf" linenums="42"
    output {
        first_output {
            path { "${params.batch}/intermediates/${sayHello.name}" }
        }
        uppercased {
            path { "${params.batch}/intermediates/${convertToUpper.name}" }
        }
        collected {
            path { "${params.batch}/intermediates/${collectGreetings.name}" }
        }
        batch_report {
            path { "${params.batch}/${collectGreetings.name}" }
        }
        cowpy_art {
            path { "${params.batch}/${cowpy.name}" }
        }
    }
    ```

=== "Vorher"

    ```groovy title="hello-config.nf" linenums="42" hl_lines="4 8 12 16 20"
    output {
        first_output {
            path { "${params.batch}/intermediates/${sayHello.name}" }
            mode 'copy'
        }
        uppercased {
            path { "${params.batch}/intermediates/${convertToUpper.name}" }
            mode 'copy'
        }
        collected {
            path { "${params.batch}/intermediates/${collectGreetings.name}" }
            mode 'copy'
        }
        batch_report {
            path { "${params.batch}/${collectGreetings.name}" }
            mode 'copy'
        }
        cowpy_art {
            path { "${params.batch}/${cowpy.name}" }
            mode 'copy'
        }
    }
    ```

Das ist prägnanter, nicht wahr?

#### 2.3.3. Führe die Pipeline aus

Lass uns testen, dass es korrekt funktioniert:

```bash
nextflow run hello-config.nf -output-dir config-output-mode
```

??? success "Befehlsausgabe"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-config.nf` [small_stone] DSL2 - revision: 024d6361b5

    executor >  local (8)
    [e8/a0e93e] sayHello (1)       [100%] 3 of 3 ✔
    [14/176c9d] convertToUpper (3) [100%] 3 of 3 ✔
    [23/d667ca] collectGreetings   [100%] 1 of 1 ✔
    [e6/1dc80e] cowpy              [100%] 1 of 1 ✔
    ```

Dies veröffentlicht Ausgaben nach `config-output-mode/`, und sie sind immer noch alle echte Kopien, keine Symlinks.

??? abstract "Verzeichnisinhalt"

    ```console
    config-output-mode
    └── batch
        ├── collectGreetings
        │   └── batch-report.txt
        ├── cowpy
        │   └── cowpy-COLLECTED-batch-output.txt
        └── intermediates
            ├── collectGreetings
            │   └── COLLECTED-batch-output.txt
            ├── convertToUpper
            │   ├── UPPER-Bonjour-output.txt
            │   ├── UPPER-Hello-output.txt
            │   └── UPPER-Holà-output.txt
            └── sayHello
                ├── Bonjour-output.txt
                ├── Hello-output.txt
                └── Holà-output.txt
    ```

Der Hauptgrund, warum du möglicherweise immer noch die pro-Ausgabe-Methode zum Setzen des Modus verwenden möchtest, ist, wenn du innerhalb desselben Workflows mischen und anpassen möchtest, _d.h._ einige Ausgaben kopiert und einige als Symlinks haben möchtest.

Es gibt viele andere Optionen, die du auf diese Weise anpassen kannst, aber hoffentlich gibt dir dies ein Gefühl für die Bandbreite der Optionen und wie du sie effektiv nutzen kannst, um deinen Vorlieben zu entsprechen.

### Fazit

Du weißt, wie du die Benennung und Struktur der Verzeichnisse steuerst, in denen deine Ausgaben veröffentlicht werden, sowie den Workflow-Ausgabe-Veröffentlichungsmodus.

### Wie geht es weiter?

Lerne, wie du deine Workflow-Konfiguration an deine Rechenumgebung anpasst, beginnend mit der Software-Packaging-Technologie.

---

## 3. Wähle eine Software-Packaging-Technologie

Bisher haben wir uns Konfigurationselemente angesehen, die steuern, wie Eingaben hineingehen und wo Ausgaben herauskommen. Jetzt ist es an der Zeit, uns spezifischer darauf zu konzentrieren, deine Workflow-Konfiguration an deine Rechenumgebung anzupassen.

Der erste Schritt auf diesem Weg ist die Angabe, woher die Softwarepakete kommen, die in jedem Schritt ausgeführt werden.
Sind sie bereits in der lokalen Rechenumgebung installiert?
Müssen wir Images abrufen und sie über ein Container-System ausführen?
Oder müssen wir Conda-Pakete abrufen und eine lokale Conda-Umgebung erstellen?

Im allerersten Teil dieses Trainingskurses (Teile 1-4) haben wir einfach lokal installierte Software in unserem Workflow verwendet.
Dann haben wir in Teil 5 Docker-Container und die `nextflow.config`-Datei eingeführt, die wir verwendet haben, um die Verwendung von Docker-Containern zu aktivieren.

Schauen wir uns nun an, wie wir eine alternative Software-Packaging-Option über die `nextflow.config`-Datei konfigurieren können.

### 3.1. Deaktiviere Docker und aktiviere Conda in der Config-Datei

Lass uns so tun, als würden wir auf einem HPC-Cluster arbeiten und der Admin erlaubt aus Sicherheitsgründen nicht die Verwendung von Docker.
Glücklicherweise für uns unterstützt Nextflow mehrere andere Container-Technologien wie Singularity (das auf HPC weiter verbreitet ist) und Software-Paketmanager wie Conda.

Wir können unsere Konfigurationsdatei ändern, um [Conda](https://nextflow.io/docs/latest/conda.html) anstelle von Docker zu verwenden.
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

### 3.2. Gib ein Conda-Paket in der Prozessdefinition an

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

### 3.3. Führe den Workflow aus, um zu überprüfen, dass er Conda verwenden kann

Lass es uns ausprobieren.

```bash
nextflow run hello-config.nf --batch conda
```

??? success "Befehlsausgabe"

    ```console title="Ausgabe"
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-config.nf` [friendly_lamport] DSL2 - revision: 024d6361b5

    executor >  local (8)
    [e8/91c116] sayHello (2)       [100%] 3 of 3 ✔
    [fe/6a70ce] convertToUpper (3) [100%] 3 of 3 ✔
    [99/7cc493] collectGreetings   [100%] 1 of 1 ✔
    [3c/09fb59] cowpy              [100%] 1 of 1 ✔
    ```

Dies sollte problemlos funktionieren und die gleichen Ausgaben wie zuvor unter `custom-outdir-config/conda` erzeugen.

Hinter den Kulissen hat Nextflow die Conda-Pakete abgerufen und die Umgebung erstellt, was normalerweise etwas Arbeit erfordert; es ist also schön, dass wir das nicht selbst tun müssen!

!!! note

    Dies läuft schnell, weil das `cowpy`-Paket recht klein ist, aber wenn du mit großen Paketen arbeitest, kann es beim ersten Mal etwas länger als üblich dauern, und du siehst möglicherweise, dass die Konsolenausgabe für eine Minute oder so 'hängen bleibt', bevor sie abgeschlossen ist.
    Dies ist normal und liegt an der zusätzlichen Arbeit, die Nextflow beim ersten Mal leistet, wenn du ein neues Paket verwendest.

Aus unserer Sicht sieht es so aus, als würde es genauso funktionieren wie die Ausführung mit Docker, obwohl die Mechanik im Backend etwas anders ist.

Das bedeutet, dass wir bereit sind, bei Bedarf mit Conda-Umgebungen zu arbeiten.

??? info "Docker und Conda mischen und anpassen"

    Da diese Direktiven pro Prozess zugewiesen werden, ist es möglich, zu 'mischen und anzupassen', _d.h._ einige der Prozesse in deinem Workflow so zu konfigurieren, dass sie mit Docker laufen, und andere mit Conda, zum Beispiel, wenn die von dir verwendete Recheninfrastruktur beides unterstützt.
    In diesem Fall würdest du sowohl Docker als auch Conda in deiner Konfigurationsdatei aktivieren.
    Wenn beide für einen bestimmten Prozess verfügbar sind, priorisiert Nextflow Container.

    Und wie bereits erwähnt, unterstützt Nextflow mehrere andere Software-Packaging- und Container-Technologien, sodass du nicht nur auf diese beiden beschränkt bist.

### Fazit

Du weißt, wie du konfigurierst, welches Softwarepaket jeder Prozess verwenden soll, und wie du zwischen Technologien wechselst.

### Wie geht es weiter?

Lerne, wie du die von Nextflow verwendete Ausführungsplattform änderst, um die Arbeit tatsächlich zu erledigen.

---

## 4. Wähle eine Ausführungsplattform

Bis jetzt haben wir unsere Pipeline mit dem lokalen Executor ausgeführt.
Dieser führt jede Aufgabe auf der Maschine aus, auf der Nextflow läuft.
Wenn Nextflow beginnt, schaut es sich die verfügbaren CPUs und den Speicher an.
Wenn die Ressourcen der zur Ausführung bereiten Aufgaben die verfügbaren Ressourcen überschreiten, hält Nextflow die letzten Aufgaben von der Ausführung zurück, bis eine oder mehrere der früheren Aufgaben abgeschlossen sind und die notwendigen Ressourcen freigeben.

Der lokale Executor ist bequem und effizient, aber er ist auf diese einzelne Maschine beschränkt. Für sehr große Workloads stellst du möglicherweise fest, dass deine lokale Maschine ein Engpass ist, entweder weil du eine einzelne Aufgabe hast, die mehr Ressourcen erfordert, als du verfügbar hast, oder weil du so viele Aufgaben hast, dass das Warten darauf, dass eine einzelne Maschine sie ausführt, zu lange dauern würde.

Nextflow unterstützt [viele verschiedene Executors](https://nextflow.io/docs/latest/executor.html), einschließlich HPC-Scheduler (Slurm, LSF, SGE, PBS, Moab, OAR, Bridge, HTCondor und andere) sowie Cloud-Ausführungs-Backends wie (AWS Batch, Google Cloud Batch, Azure Batch, Kubernetes und mehr).

### 4.1. Auf ein anderes Backend abzielen

Die Wahl des Executors wird durch eine Prozess-Direktive namens `executor` festgelegt.
Standardmäßig ist sie auf `local` gesetzt, sodass die folgende Konfiguration impliziert ist:

```groovy title="Eingebaute Konfiguration"
process {
    executor = 'local'
}
```

Um den Executor so einzustellen, dass er auf ein anderes Backend abzielt, würdest du einfach den gewünschten Executor mit ähnlicher Syntax angeben, wie oben für Ressourcenzuweisungen beschrieben (siehe [Executor-Dokumentation](https://nextflow.io/docs/latest/executor.html) für alle Optionen).

```groovy title="nextflow.config"
process {
    executor = 'slurm'
}
```

!!! warning

    Wir können dies in der Trainingsumgebung nicht tatsächlich testen, da sie nicht für die Verbindung zu einem HPC eingerichtet ist.

### 4.2. Umgang mit backend-spezifischer Syntax für Ausführungsparameter

Die meisten Hochleistungsrechner-Plattformen erlauben (und erfordern manchmal), dass du bestimmte Parameter wie Ressourcenzuweisungsanfragen und -beschränkungen (z.B. Anzahl der CPUs und Speicher) und den Namen der zu verwendenden Job-Warteschlange angibst.

Leider verwendet jedes dieser Systeme unterschiedliche Technologien, Syntaxen und Konfigurationen, um zu definieren, wie ein Job definiert und an den entsprechenden Scheduler übermittelt werden soll.

??? abstract "Beispiele"

    Zum Beispiel muss derselbe Job, der 8 CPUs und 4 GB RAM erfordert, um in der Warteschlange "my-science-work" ausgeführt zu werden, je nach Backend auf unterschiedliche Weise ausgedrückt werden.

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
Es bietet eine standardisierte Syntax, sodass du die relevanten Eigenschaften wie [`cpus`](https://nextflow.io/docs/latest/reference/process.html#cpus), [`memory`](https://nextflow.io/docs/latest/reference/process.html#memory) und [`queue`](https://nextflow.io/docs/latest/reference/process.html#queue) (siehe [Prozess-Direktiven](https://nextflow.io/docs/latest/reference/process.html#process-directives) für andere Eigenschaften) nur einmal angeben kannst.
Dann generiert Nextflow zur Laufzeit die entsprechenden backend-spezifischen Skripte basierend auf der Executor-Einstellung.

Wir werden diese standardisierte Syntax im nächsten Abschnitt behandeln.

### Fazit

Du weißt jetzt, wie du den Executor änderst, um verschiedene Arten von Recheninfrastruktur zu verwenden.

### Wie geht es weiter?

Lerne, wie du Ressourcenzuweisungen und -beschränkungen in Nextflow bewertest und ausdrückst.

---

## 5. Steuere Rechenressourcenzuweisungen

Die meisten Hochleistungsrechner-Plattformen erlauben (und erfordern manchmal), dass du bestimmte Ressourcenzuweisungsparameter wie Anzahl der CPUs und Speicher angibst.

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

### 5.1. Führe den Workflow aus, um einen Ressourcennutzungsbericht zu generieren

Wenn du im Voraus nicht weißt, wie viel CPU und Speicher deine Prozesse wahrscheinlich benötigen, kannst du ein Ressourcen-Profiling durchführen, was bedeutet, dass du den Workflow mit einigen Standardzuweisungen ausführst, aufzeichnest, wie viel jeder Prozess verwendet hat, und von dort aus abschätzt, wie die Basiszuweisungen angepasst werden sollten.

Praktischerweise enthält Nextflow eingebaute Tools dafür und erstellt gerne auf Anfrage einen Bericht für dich.

Füge dazu `-with-report <dateiname>.html` zu deiner Befehlszeile hinzu.

```bash
nextflow run hello-config.nf -with-report report-config-1.html
```

Der Bericht ist eine HTML-Datei, die du herunterladen und in deinem Browser öffnen kannst. Du kannst auch mit der rechten Maustaste darauf im Datei-Explorer auf der linken Seite klicken und auf `Show preview` klicken, um ihn in der Trainingsumgebung anzuzeigen.

Nimm dir ein paar Minuten Zeit, um den Bericht durchzusehen und zu sehen, ob du einige Möglichkeiten zur Anpassung von Ressourcen identifizieren kannst.
Stelle sicher, auf die Registerkarten zu klicken, die die Nutzungsergebnisse als Prozentsatz dessen zeigen, was zugewiesen wurde.

Siehe [Reports](https://nextflow.io/docs/latest/reports.html) für Dokumentation zu allen verfügbaren Funktionen.

### 5.2. Setze Ressourcenzuweisungen für alle Prozesse

Das Profiling zeigt, dass die Prozesse in unserem Trainings-Workflow sehr leichtgewichtig sind, also lass uns die Standard-Speicherzuweisung auf 1 GB pro Prozess reduzieren.

Füge das Folgende zu deiner `nextflow.config`-Datei hinzu, vor dem Abschnitt Pipeline-Parameter:

=== "Danach"

    ```groovy title="nextflow.config" linenums="1" hl_lines="4-9"
    docker.enabled = false
    conda.enabled = true

    /*
    * Prozesseinstellungen
    */
    process {
        memory = 1.GB
    }

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
    docker.enabled = false
    conda.enabled = true

    /*
    * Pipeline-Parameter
    */
    params {
        input = 'data/greetings.csv'
        batch = 'batch'
        character = 'turkey'
    }
    ```

Das wird helfen, die Menge an Rechenleistung zu reduzieren, die wir verbrauchen.

### 5.3. Setze Ressourcenzuweisungen für einen bestimmten Prozess

Gleichzeitig werden wir so tun, als würde der `cowpy`-Prozess mehr Ressourcen als die anderen benötigen, nur damit wir demonstrieren können, wie man Zuweisungen für einen einzelnen Prozess anpasst.

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

!!! tip

    Wenn du eine Maschine mit wenigen CPUs hast und eine hohe Anzahl pro Prozess zuweist, siehst du möglicherweise, dass Prozessaufrufe hintereinander in die Warteschlange gestellt werden.
    Dies liegt daran, dass Nextflow sicherstellt, dass wir nicht mehr CPUs anfordern, als verfügbar sind.

### 5.4. Führe den Workflow mit der aktualisierten Konfiguration aus

Lass uns das ausprobieren und einen anderen Dateinamen für den Profiling-Bericht angeben, damit wir die Leistung vor und nach den Konfigurationsänderungen vergleichen können.

```bash
nextflow run hello-config.nf -with-report report-config-2.html
```

Du wirst wahrscheinlich keinen wirklichen Unterschied bemerken, da dies eine so kleine Arbeitslast ist, aber dies ist der Ansatz, den du verwenden würdest, um die Leistung und Ressourcenanforderungen eines realen Workflows zu analysieren.

Es ist sehr nützlich, wenn deine Prozesse unterschiedliche Ressourcenanforderungen haben. Es ermöglicht dir, die Ressourcenzuweisungen, die du für jeden Prozess einrichtest, basierend auf tatsächlichen Daten richtig zu dimensionieren, nicht auf Vermutungen.

!!! tip

    Dies ist nur ein kleiner Vorgeschmack darauf, was du tun kannst, um deine Ressourcennutzung zu optimieren.
    Nextflow selbst hat einige wirklich nette [dynamische Wiederholungslogik](https://nextflow.io/docs/latest/process.html#dynamic-task-resources) eingebaut, um Jobs, die aufgrund von Ressourcenbeschränkungen fehlschlagen, erneut zu versuchen.
    Zusätzlich bietet die Seqera Platform KI-gesteuerte Tools zur automatischen Optimierung deiner Ressourcenzuweisungen.

### 5.5. Füge Ressourcenlimits hinzu

Abhängig davon, welchen Computing-Executor und welche Recheninfrastruktur du verwendest, kann es einige Einschränkungen geben, was du zuweisen kannst (oder musst).
Zum Beispiel kann dein Cluster verlangen, dass du innerhalb bestimmter Grenzen bleibst.

Du kannst die `resourceLimits`-Direktive verwenden, um die relevanten Beschränkungen festzulegen. Die Syntax sieht so aus, wenn sie allein in einem Prozess-Block steht:

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
Wenn du jedoch versuchen würdest, den Workflow mit Ressourcenzuweisungen auszuführen, die diese Grenzen überschreiten, und dann den `sbatch`-Befehl in der `.command.run`-Skriptdatei nachschlagen würdest, würdest du sehen, dass die Anfragen, die tatsächlich an den Executor gesendet werden, auf die von `resourceLimits` angegebenen Werte begrenzt sind.

??? info "Institutionelle Referenzkonfigurationen"

    Das nf-core-Projekt hat eine [Sammlung von Konfigurationsdateien](https://nf-co.re/configs/) zusammengestellt, die von verschiedenen Institutionen auf der ganzen Welt geteilt werden und eine breite Palette von HPC- und Cloud-Executors abdecken.

    Diese geteilten Configs sind sowohl für Personen wertvoll, die dort arbeiten und daher einfach die Konfiguration ihrer Institution sofort nutzen können, als auch als Modell für Personen, die eine Konfiguration für ihre eigene Infrastruktur entwickeln möchten.

### Fazit

Du weißt, wie du einen Profiling-Bericht generierst, um die Ressourcennutzung zu bewerten, und wie du Ressourcenzuweisungen für alle Prozesse und/oder für einzelne Prozesse änderst sowie Ressourcenbeschränkungen für die Ausführung auf HPC festlegst.

### Wie geht es weiter?

Lerne, wie du voreingestellte Konfigurationsprofile einrichtest und zur Laufzeit zwischen ihnen wechselst.

---

## 6. Verwende Profile, um zwischen voreingestellten Konfigurationen zu wechseln

Wir haben dir eine Reihe von Möglichkeiten gezeigt, wie du deine Pipeline-Konfiguration je nach Projekt, an dem du arbeitest, oder der Rechenumgebung, die du verwendest, anpassen kannst.

Du möchtest möglicherweise zwischen alternativen Einstellungen wechseln, je nachdem, welche Recheninfrastruktur du verwendest. Zum Beispiel möchtest du vielleicht lokal auf deinem Laptop entwickeln und kleine Tests durchführen und dann vollständige Workloads auf HPC oder Cloud ausführen.

Nextflow ermöglicht es dir, eine beliebige Anzahl von [Profilen](https://nextflow.io/docs/latest/config.html#config-profiles) einzurichten, die verschiedene Konfigurationen beschreiben, die du dann zur Laufzeit mit einem Befehlszeilenargument auswählen kannst, anstatt die Konfigurationsdatei selbst ändern zu müssen.

### 6.1. Erstelle Profile zum Wechseln zwischen lokaler Entwicklung und Ausführung auf HPC

Lass uns zwei alternative Profile einrichten; eines für die Ausführung kleiner Lasten auf einem normalen Computer, wo wir Docker-Container verwenden, und eines für die Ausführung auf einem Universitäts-HPC mit einem Slurm-Scheduler, wo wir Conda-Pakete verwenden.

#### 6.1.1. Richte die Profile ein

Füge das Folgende zu deiner `nextflow.config`-Datei hinzu, nach dem Abschnitt Pipeline-Parameter, aber vor den Ausgabeeinstellungen:

=== "Danach"

    ```groovy title="nextflow.config" linenums="15" hl_lines="10-27"
    /*
    * Pipeline-Parameter
    */
    params {
        input = 'data/greetings.csv'
        batch = 'batch'
        character = 'turkey'
    }

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

    /*
    * Ausgabeeinstellungen
    */
    outputDir = "custom-outdir-config/"
    workflow.output.mode = 'copy'
    ```

=== "Vorher"

    ```groovy title="nextflow.config" linenums="15"
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
    outputDir = "custom-outdir-config/"
    workflow.output.mode = 'copy'
    ```

Du siehst, dass wir für den Universitäts-HPC auch Ressourcenbeschränkungen angeben.

#### 6.1.2. Führe den Workflow mit einem Profil aus

Um ein Profil in unserer Nextflow-Befehlszeile anzugeben, verwenden wir das Argument `-profile`.

Lass uns versuchen, den Workflow mit der `my_laptop`-Konfiguration auszuführen.

```bash
nextflow run hello-config.nf -profile my_laptop
```

??? success "Befehlsausgabe"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-config.nf` [hungry_sanger] DSL2 - revision: 024d6361b5

    executor >  local (8)
    [b0/fb2ec9] sayHello (3)       [100%] 3 of 3 ✔
    [4a/e039f0] convertToUpper (3) [100%] 3 of 3 ✔
    [6f/408fa9] collectGreetings   [100%] 1 of 1 ✔
    [f1/fd6520] cowpy              [100%] 1 of 1 ✔
    ```

Wie du sehen kannst, ermöglicht uns dies, sehr bequem zur Laufzeit zwischen Konfigurationen zu wechseln.

!!! warning

    Das `univ_hpc`-Profil wird in der Trainingsumgebung nicht ordnungsgemäß ausgeführt, da wir keinen Zugriff auf einen Slurm-Scheduler haben.

Wenn wir in Zukunft andere Konfigurationselemente finden, die immer mit diesen zusammen auftreten, können wir sie einfach zum entsprechenden Profil(en) hinzufügen.
Wir können auch zusätzliche Profile erstellen, wenn es andere Konfigurationselemente gibt, die wir zusammen gruppieren möchten.

### 6.2. Erstelle ein Profil von Testparametern

Profile sind nicht nur für Infrastrukturkonfiguration. Wir können sie auch verwenden, um Standardwerte für Workflow-Parameter festzulegen, um es anderen zu erleichtern, den Workflow auszuprobieren, ohne selbst geeignete Eingabewerte sammeln zu müssen.
Du kannst dies als Alternative zur Verwendung einer Parameterdatei betrachten.

#### 6.2.1. Richte das Profil ein

Die Syntax zum Ausdrücken von Standardwerten in diesem Kontext sieht so aus, für ein Profil, das wir `test` nennen:

```groovy title="Syntaxbeispiel"
    test {
        params.<parameter1>
        params.<parameter2>
        ...
    }
```

Wenn wir ein Testprofil für unseren Workflow hinzufügen, wird der `profiles`-Block zu:

```groovy title="nextflow.config" linenums="24" hl_lines="18-22"
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

#### 6.2.2. Führe den Workflow lokal mit dem Testprofil aus

Praktischerweise schließen sich Profile nicht gegenseitig aus, sodass wir mehrere Profile in unserer Befehlszeile mit der folgenden Syntax angeben können: `-profile <profil1>,<profil2>` (für eine beliebige Anzahl von Profilen).

Wenn du Profile kombinierst, die Werte für dieselben Konfigurationselemente setzen und in derselben Konfigurationsdatei beschrieben sind, löst Nextflow den Konflikt, indem es den Wert verwendet, den es zuletzt gelesen hat (_d.h._ was später in der Datei kommt).
Wenn die widersprüchlichen Einstellungen in verschiedenen Konfigurationsquellen gesetzt sind, gilt die Standard-[Rangfolge](https://nextflow.io/docs/latest/config.html).

Lass uns versuchen, das Testprofil zu unserem vorherigen Befehl hinzuzufügen:

```bash
nextflow run hello-config.nf -profile my_laptop,test
```

??? success "Befehlsausgabe"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-config.nf` [modest_becquerel] DSL2 - revision: 024d6361b5

    executor >  local (8)
    [4c/fe2580] sayHello (1)       [100%] 3 of 3 ✔
    [fd/7d9017] convertToUpper (3) [100%] 3 of 3 ✔
    [13/1523bd] collectGreetings   [100%] 1 of 1 ✔
    [06/a1ee14] cowpy              [100%] 1 of 1 ✔
    ```

Dies verwendet Docker, wo möglich, und erzeugt Ausgaben unter `custom-outdir-config/test`, und diesmal ist das Zeichen das komödiantische Duo `dragonandcow`.

??? abstract "Dateiinhalt"

    ```console title="custom-outdir-config/test/cowpy/cowpy-COLLECTED-test-output.txt"
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

Dies bedeutet, dass jeder, solange wir alle Testdatendateien mit dem Workflow-Code verteilen, den Workflow schnell ausprobieren kann, ohne eigene Eingaben über die Befehlszeile oder eine Parameterdatei bereitstellen zu müssen.

!!! tip

    Wir können auf URLs für größere Dateien verweisen, die extern gespeichert sind.
    Nextflow lädt sie automatisch herunter, solange eine offene Verbindung besteht.

    Für weitere Details siehe die Side Quest [Working with Files](../side_quests/working_with_files.md)

### 6.3. Verwende `nextflow config`, um die aufgelöste Konfiguration zu sehen

Wie oben erwähnt, kann manchmal derselbe Parameter in Profilen, die du kombinieren möchtest, auf unterschiedliche Werte gesetzt werden.
Und allgemeiner gibt es zahlreiche Orte, an denen Konfigurationselemente gespeichert werden können, und manchmal können dieselben Eigenschaften an verschiedenen Orten auf unterschiedliche Werte gesetzt werden.

Nextflow wendet eine festgelegte [Rangfolge](https://nextflow.io/docs/latest/config.html) an, um Konflikte zu lösen, aber das kann schwierig sein, selbst zu bestimmen.
Und selbst wenn nichts in Konflikt steht, kann es mühsam sein, alle möglichen Orte nachzuschlagen, an denen Dinge konfiguriert werden könnten.

Glücklicherweise enthält Nextflow ein praktisches Dienstprogramm-Tool namens `config`, das diesen gesamten Prozess für dich automatisieren kann.

Das `config`-Tool durchsucht alle Inhalte in deinem aktuellen Arbeitsverzeichnis, sammelt alle Konfigurationsdateien ein und erstellt die vollständig aufgelöste Konfiguration, die Nextflow zum Ausführen des Workflows verwenden würde.
Dies ermöglicht es dir, herauszufinden, welche Einstellungen verwendet werden, ohne etwas starten zu müssen.

#### 6.3.1. Löse die Standardkonfiguration auf

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

    outputDir = 'custom-outdir-config/'

    workflow {
      output {
          mode = 'copy'
      }
    }
    ```

Dies zeigt dir die Basiskonfiguration, die du erhältst, wenn du nichts Zusätzliches in der Befehlszeile angibst.

#### 6.3.2. Löse die Konfiguration mit aktivierten spezifischen Einstellungen auf

Wenn du Befehlszeilenparameter bereitstellst, z.B. ein oder mehrere Profile aktivierst oder eine Parameterdatei lädst, berücksichtigt der Befehl diese zusätzlich.

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

    outputDir = 'custom-outdir-config/'

    workflow {
      output {
          mode = 'copy'
      }
    }
    ```

Dies wird besonders nützlich für komplexe Projekte, die mehrere Konfigurationsebenen beinhalten.

### Fazit

Du weißt, wie du Profile verwendest, um eine voreingestellte Konfiguration zur Laufzeit mit minimalem Aufwand auszuwählen.
Allgemeiner weißt du, wie du deine Workflow-Ausführungen konfigurierst, um verschiedenen Rechenplattformen zu entsprechen und die Reproduzierbarkeit deiner Analysen zu verbessern.

### Wie geht es weiter?

Feiere und klopfe dir selbst auf die Schulter! Du hast deinen allerersten Nextflow-Entwicklerkurs abgeschlossen.

Gehe zur finalen [Kurszusammenfassung](./next_steps.md), um zu überprüfen, was du gelernt hast, und herauszufinden, was als Nächstes kommt.

---

## Quiz

<quiz>
Wie heißt die Konfigurationsdatei, die Nextflow automatisch lädt?
- [ ] `config.nf`
- [ ] `pipeline.config`
- [x] `nextflow.config`
- [ ] `workflow.config`
</quiz>

<quiz>
Was hat Vorrang, wenn derselbe Parameter sowohl in der Config-Datei als auch in der Befehlszeile gesetzt ist?
- [ ] Der Wert aus der Config-Datei
- [x] Der Wert aus der Befehlszeile
- [ ] Der zuerst angetroffene Wert
- [ ] Keiner; es verursacht einen Fehler

Mehr erfahren: [1.1. Standardwerte in `nextflow.config` verschieben](#11-standardwerte-in-nextflowconfig-verschieben)
</quiz>

<quiz>
Kannst du sowohl Docker als auch Conda in derselben Konfiguration aktiviert haben?
- [x] Ja, Nextflow kann beide je nach Prozess-Direktiven verwenden
- [ ] Nein, nur eines kann gleichzeitig aktiviert sein
- [ ] Ja, aber nur in Profilen
- [ ] Nein, sie schließen sich gegenseitig aus
</quiz>

<quiz>
Wenn sowohl Docker als auch Conda aktiviert sind und ein Prozess beide Direktiven hat, welche wird priorisiert?
- [x] Docker (Container)
- [ ] Conda
- [ ] Die zuerst definierte
- [ ] Es verursacht einen Fehler

Mehr erfahren: [3. Wähle eine Software-Packaging-Technologie](#3-wahle-eine-software-packaging-technologie)
</quiz>

<quiz>
Was ist die Standard-Speicherzuweisung für Nextflow-Prozesse?
- [ ] 1 GB
- [x] 2 GB
- [ ] 4 GB
- [ ] Kein Limit
</quiz>

<quiz>
Wie setzt du Ressourcenanforderungen für einen bestimmten Prozess in der Config-Datei?
- [ ] `#!groovy processName.memory = '4 GB'`
- [ ] `#!groovy process.memory.processName = '4 GB'`
- [x] `#!groovy process { withName: 'processName' { memory = '4 GB' } }`
- [ ] `#!groovy resources.processName.memory = '4 GB'`

Mehr erfahren: [5.3. Setze Ressourcenzuweisungen für einen bestimmten Prozess](#53-setze-ressourcenzuweisungen-fur-einen-bestimmten-prozess)
</quiz>

<quiz>
Welche Befehlszeilenoption generiert einen Ressourcennutzungsbericht?
- [ ] `-with-metrics`
- [ ] `-with-stats`
- [x] `-with-report`
- [ ] `-with-profile`

Mehr erfahren: [5.1. Führe den Workflow aus, um einen Ressourcennutzungsbericht zu generieren](#51-fuhre-den-workflow-aus-um-einen-ressourcennutzungsbericht-zu-generieren)
</quiz>

<quiz>
Was macht die `resourceLimits`-Direktive?
- [ ] Setzt Mindestressourcenanforderungen
- [ ] Weist Prozessen Ressourcen zu
- [x] Begrenzt die maximalen Ressourcen, die angefordert werden können
- [ ] Überwacht die Ressourcennutzung

Mehr erfahren: [5.5. Füge Ressourcenlimits hinzu](#55-fuge-ressourcenlimits-hinzu)
</quiz>

<quiz>
Was ist der Standard-Executor in Nextflow?
- [x] `local`
- [ ] `slurm`
- [ ] `kubernetes`
- [ ] `aws`

Mehr erfahren: [4. Wähle eine Ausführungsplattform](#4-wahle-eine-ausfuhrungsplattform)
</quiz>

<quiz>
Wie gibst du eine Parameterdatei beim Ausführen von Nextflow an?
- [ ] `--params params.json`
- [ ] `-config params.json`
- [x] `-params-file params.json`
- [ ] `--input params.json`

Mehr erfahren: [1.3. Verwende eine Parameterdatei](#13-verwende-eine-parameterdatei)
</quiz>

<quiz>
Wofür können Profile verwendet werden? (Wähle alle zutreffenden aus)
- [x] Definieren infrastrukturspezifischer Einstellungen
- [x] Setzen von Ressourcenlimits für verschiedene Umgebungen
- [x] Bereitstellen von Testparametern
- [ ] Definieren neuer Prozesse

Mehr erfahren: [6. Verwende Profile, um zwischen voreingestellten Konfigurationen zu wechseln](#6-verwende-profile-um-zwischen-voreingestellten-konfigurationen-zu-wechseln)
</quiz>

<quiz>
Wie gibst du mehrere Profile in einem einzigen Befehl an?
- [ ] `-profile profil1 -profile profil2`
- [ ] `-profiles profil1,profil2`
- [x] `-profile profil1,profil2`
- [ ] `--profile profil1 --profile profil2`

Mehr erfahren: [6. Verwende Profile, um zwischen voreingestellten Konfigurationen zu wechseln](#6-verwende-profile-um-zwischen-voreingestellten-konfigurationen-zu-wechseln)
</quiz>
