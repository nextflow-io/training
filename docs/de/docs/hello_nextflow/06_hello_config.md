# Teil 6: Hello Config

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } KI-gestützte Übersetzung - [mehr erfahren & Verbesserungen vorschlagen](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/FcZTiE25TeA?si=tnXTi6mRkITY0zW_&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1&amp;cc_lang_pref=de" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

/// caption
:fontawesome-brands-youtube:{ .youtube } Sieh dir [die ganze Playlist](https://youtube.com/playlist?list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&si=eF7cLR62goy-lc6n) auf dem Nextflow YouTube-Kanal an.

:green_book: Das Videotranskript ist [hier](./transcripts/06_hello_config.md) verfügbar.
///

In diesem Abschnitt lernst du, wie du die Konfiguration deiner Nextflow-Pipeline einrichtest und verwaltest. Du wirst in der Lage sein, ihr Verhalten anzupassen, sie an verschiedene Umgebungen anzupassen und die Ressourcennutzung zu optimieren, _ohne eine einzige Zeile des Workflow-Codes selbst zu ändern_.

Es gibt mehrere Möglichkeiten, dies zu tun, die kombiniert werden können und gemäß der in der Konfigurationsdokumentation beschriebenen [Vorrangordnung](https://nextflow.io/docs/latest/config.html) interpretiert werden.

In diesem Teil des Kurses zeigen wir dir den einfachsten und häufigsten Konfigurationsdatei-Mechanismus, die [`nextflow.config`](https://nextflow.io/docs/latest/config.html)-Datei, die du bereits in Teil 5: Hello Containers kennengelernt hast.

Wir behandeln wesentliche Komponenten der Nextflow-Konfiguration wie Prozess-Direktiven, Executors, Profile und Parameterdateien.
Indem du lernst, diese Konfigurationsoptionen effektiv zu nutzen, kannst du die Flexibilität, Skalierbarkeit und Leistung deiner Pipelines verbessern.

??? info "Wie du von diesem Abschnitt aus beginnst"

    Dieser Abschnitt des Kurses setzt voraus, dass du die Teile 1-5 des [Hello Nextflow](./index.md)-Kurses abgeschlossen hast und eine vollständig funktionierende Pipeline hast.

    Wenn du den Kurs von diesem Punkt aus beginnst, musst du das `modules`-Verzeichnis und die `nextflow.config`-Datei aus den Lösungen kopieren:

    ```bash
    cp -r solutions/5-hello-containers/modules .
    cp solutions/5-hello-containers/nextflow.config .
    ```

    Die `nextflow.config`-Datei enthält die Zeile `docker.enabled = true`, die die Verwendung von Docker-Containern aktiviert.

    Wenn du mit der Hello-Pipeline nicht vertraut bist oder eine Auffrischung brauchst, sieh dir [diese Infoseite](../info/hello_pipeline.md) an.

---

## 0. Aufwärmen: `hello-config.nf` ausführen

Wir werden das Workflow-Script `hello-config.nf` als Ausgangspunkt verwenden.
Es entspricht dem Script, das durch Durcharbeiten von Teil 5 dieses Kurses entstanden ist, außer dass wir die Ausgabeziele geändert haben:

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

Um sicherzustellen, dass alles funktioniert, führe das Script einmal aus, bevor du Änderungen vornimmst:

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

Wie zuvor findest du die Ausgabedateien in dem im `output`-Block angegebenen Verzeichnis (`results/hello_config/`).

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

Die finale ASCII-Kunst-Ausgabe befindet sich im Verzeichnis `results/hello_config/`, unter dem Namen `cowpy-COLLECTED-batch-output.txt`.

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

Wenn das bei dir funktioniert hat, bist du bereit zu lernen, wie man Pipelines konfiguriert.

---

## 1. Workflow-Eingabeparameter verwalten

Wir beginnen mit einem Aspekt der Konfiguration, der einfach eine Erweiterung dessen ist, womit wir bisher gearbeitet haben: die Verwaltung von Eingabeparametern.

Derzeit ist unser Workflow so eingerichtet, dass er mehrere Parameterwerte über die Befehlszeile akzeptiert, wobei Standardwerte in einem `params`-Block im Workflow-Script selbst festgelegt sind.
Möglicherweise möchtest du diese Standardwerte jedoch überschreiben, ohne entweder Parameter in der Befehlszeile angeben oder die ursprüngliche Script-Datei ändern zu müssen.

Es gibt mehrere Möglichkeiten, das zu tun; wir zeigen dir drei grundlegende Wege, die sehr häufig verwendet werden.

### 1.1. Standardwerte in `nextflow.config` verschieben

Dies ist der einfachste Ansatz, obwohl er möglicherweise am wenigsten flexibel ist, da die Haupt-`nextflow.config`-Datei nicht etwas ist, das du für jeden Lauf bearbeiten möchtest.
Aber es hat den Vorteil, die Belange des _Deklarierens_ der Parameter im Workflow (was definitiv dorthin gehört) vom Bereitstellen von _Standardwerten_ zu trennen, die eher in einer Konfigurationsdatei zu Hause sind.

Lass uns das in zwei Schritten tun.

#### 1.1.1. Einen `params`-Block in der Konfigurationsdatei erstellen

Nimm folgende Code-Änderungen in der `nextflow.config`-Datei vor:

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

Beachte, dass wir den `params`-Block nicht einfach vom Workflow in die Konfigurationsdatei kopiert haben.
Die Syntax ist etwas anders.
In der Workflow-Datei sind das typisierte Deklarationen.
In der Konfiguration sind das Wertzuweisungen.

Technisch gesehen reicht dies aus, um die Standardwerte zu überschreiben, die noch in der Workflow-Datei angegeben sind.
Du könntest zum Beispiel den Charakter ändern und den Workflow ausführen, um dich zu vergewissern, dass der in der Konfigurationsdatei festgelegte Wert den in der Workflow-Datei festgelegten überschreibt.

Aber im Sinne der vollständigen Verlagerung der Konfiguration in die Konfigurationsdatei, entfernen wir diese Werte vollständig aus der Workflow-Datei.

#### 1.1.2. Die Werte aus dem `params`-Block in der Workflow-Datei entfernen

Nimm folgende Code-Änderungen an der `hello-config.nf`-Workflow-Datei vor:

=== "Nachher"

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

#### 1.1.3. Die Pipeline ausführen

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

Dies produziert dieselbe Ausgabe wie zuvor.

Die finale ASCII-Kunst-Ausgabe befindet sich im Verzeichnis `results/hello_config/`, unter dem Namen `cowpy-COLLECTED-batch-output.txt`, wie zuvor.

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

Funktional hat diese Verschiebung nichts geändert, aber konzeptionell ist es etwas sauberer, die Standardwerte in der Konfigurationsdatei zu haben.

### 1.2. Eine laufspezifische Konfigurationsdatei verwenden

Das ist großartig, aber manchmal möchtest du vielleicht einige temporäre Experimente mit anderen Standardwerten durchführen, ohne die Hauptkonfigurationsdatei zu ändern.
Du kannst das tun, indem du eine neue `nextflow.config`-Datei in einem Unterverzeichnis erstellst, das du als Arbeitsverzeichnis für deine Experimente verwendest.

#### 1.2.1. Das Arbeitsverzeichnis mit einer leeren Konfiguration erstellen

Beginnen wir mit dem Erstellen eines neuen Verzeichnisses und wechseln hinein:

```bash
mkdir -p tux-run
cd tux-run
```

Dann erstelle eine leere Konfigurationsdatei in diesem Verzeichnis:

```bash
touch nextflow.config
```

Dies erzeugt eine leere Datei.

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

Wir können jetzt unsere Pipeline aus unserem neuen Arbeitsverzeichnis ausführen.
Stelle sicher, dass du den Pfad entsprechend anpasst!

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

Dies erstellt einen neuen Satz von Verzeichnissen unter `tux-run/` einschließlich `tux-run/work/` und `tux-run/results/`.

Bei diesem Lauf kombiniert Nextflow die `nextflow.config` in unserem aktuellen Verzeichnis mit der `nextflow.config` im Stammverzeichnis der Pipeline und überschreibt dadurch den Standard-Charakter (turkey) mit dem tux-Charakter.

Die finale Ausgabedatei sollte den tux-Charakter enthalten, der die Grüße sagt.

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

Das war's; jetzt hast du einen Bereich zum Experimentieren, ohne deine 'normale' Konfiguration zu ändern.

!!! warning "Warnung"

    Stelle sicher, dass du zum vorherigen Verzeichnis zurückwechselst, bevor du zum nächsten Abschnitt übergehst!

    ```bash
    cd ..
    ```

Schauen wir uns jetzt eine weitere nützliche Möglichkeit an, Parameterwerte festzulegen.

### 1.3. Eine Parameterdatei verwenden

Der Unterverzeichnis-Ansatz funktioniert gut zum Experimentieren, erfordert aber etwas Einrichtung und dass du Pfade entsprechend anpasst.
Es gibt einen einfacheren Ansatz, wenn du deine Pipeline mit einem bestimmten Satz von Werten ausführen möchtest, oder es jemand anderem mit minimalem Aufwand ermöglichen möchtest.

Nextflow ermöglicht es uns, Parameter über eine [Parameterdatei](https://nextflow.io/docs/latest/config.html#params-file) im YAML- oder JSON-Format anzugeben, was es sehr bequem macht, alternative Sätze von Standardwerten sowie laufspezifische Parameterwerte zu verwalten und zu verteilen.

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

    Wir stellen auch eine JSON-Version der Parameterdatei als Beispiel bereit, aber wir werden hier nicht damit arbeiten.
    Probiere sie gerne selbst aus.

#### 1.3.2. Die Pipeline ausführen

Um den Workflow mit dieser Parameterdatei auszuführen, füge einfach `-params-file <filename>` zum Basisbefehl hinzu.

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

Die finale Ausgabedatei sollte den stegosaurus-Charakter enthalten, der die Grüße sagt.

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

Die Verwendung einer Parameterdatei mag übertrieben erscheinen, wenn du nur wenige Parameter angeben musst, aber einige Pipelines erwarten Dutzende von Parametern.
In diesen Fällen ermöglicht uns die Verwendung einer Parameterdatei, Parameterwerte zur Laufzeit bereitzustellen, ohne massive Befehlszeilen eingeben oder das Workflow-Script ändern zu müssen.

Es macht es auch einfacher, Parametersätze an Kollegen zu verteilen oder als ergänzende Informationen für eine Veröffentlichung bereitzustellen.
Das macht deine Arbeit reproduzierbarer für andere.

### Fazit

Du weißt, wie du die wichtigsten Konfigurationsoptionen für die Verwaltung von Workflow-Eingaben nutzen kannst.

### Wie geht es weiter?

Lerne, wie du verwaltest, wo und wie deine Workflow-Ausgaben veröffentlicht werden.

---

## 2. Workflow-Ausgaben verwalten

Bisher haben wir alle Pfade für Workflow-Level-Ausgabedeklarationen hartcodiert, und wie wir beim Hinzufügen mehrerer Ausgaben festgestellt haben, kann dabei etwas Wiederholung auftreten.

Schauen wir uns einige gängige Möglichkeiten an, wie du dies flexibler konfigurieren kannst.

### 2.1. Das Ausgabeverzeichnis mit `-output-dir` anpassen

Wenn wir steuern, wie unsere 'veröffentlichten' Ausgaben organisiert werden, haben wir zwei unterschiedliche Prioritäten:

- Das oberste Ausgabeverzeichnis
- Wie Dateien innerhalb dieses Verzeichnisses organisiert werden

Wir haben bisher das Standard-Ausgabeverzeichnis verwendet: `results`.
Lass uns zunächst damit beginnen, dies mit der CLI-Option `-output-dir` anzupassen.

#### 2.1.1. Die Pipeline mit `-output-dir` ausführen

Die Option `-output-dir` (Kurzform: `-o`) überschreibt das Standard-Ausgabeverzeichnis (`results/`) für alle Workflow-Ausgaben.
Dies ist der empfohlene Weg, um den Basispfad zu steuern, in dem Ausgaben veröffentlicht werden.

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

Dies veröffentlicht Ausgaben nach `custom-outdir-cli/` statt nach `results/`:

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

Beachte, dass wir immer noch das Unterverzeichnis `hello_config` aus den `path`-Deklarationen im Output-Block haben.
Lass uns das aufräumen.

#### 2.1.2. Hartcodierte Pfade aus dem Output-Block entfernen

Das Präfix `hello_config/` wurde in früheren Kapiteln hartcodiert, aber da wir jetzt lernen, Ausgabepfade flexibel zu konfigurieren, können wir diese Hartcodierung entfernen.
Für Ausgaben, die kein Unterverzeichnis benötigen, können wir die `path`-Direktive auf einen leeren String setzen oder sie ganz entfernen.

Nimm folgende Code-Änderungen in der Workflow-Datei vor:

=== "Nachher"

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

Jetzt werden die Ausgaben direkt unter `custom-outdir-cli-2/` veröffentlicht, ohne das Unterverzeichnis `hello_config`:

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

!!! tip "Tipp"

    Die Option `-output-dir` wird verwendet, um zu steuern, _wohin_ Ausgaben gehen, während die `path`-Direktive im Output-Block die _Unterverzeichnisstruktur_ steuert.

### 2.2. Dynamische Ausgabepfade

Zusätzlich zum Ändern des Ausgabeverzeichnisses über die CLI können wir auch einen benutzerdefinierten Standardwert in der Config-Datei mit `outputDir` setzen.
Dies ermöglicht es uns, den Verzeichnispfad dynamisch zu setzen - nicht nur mit statischen Strings.

#### 2.2.1. `outputDir` in der Konfigurationsdatei setzen

Füge folgenden Code zur `nextflow.config`-Datei hinzu:

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
    * Ausgabe-Einstellungen
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

!!! note "Hinweis"

    Die CLI-Option `-output-dir` hat Vorrang vor der `outputDir`-Konfigurationseinstellung.
    Wenn sie gesetzt ist, wird die Config-Option vollständig ignoriert.

#### 2.2.2. Unterverzeichnisse mit Batch- und Prozessnamen

Wir können auch Unterverzeichnis-Ausgabe-`path`-Deklarationen dynamisch setzen, auf einer pro-Ausgabe-Basis.

Zum Beispiel können wir unsere Ausgaben nach Prozess organisieren, indem wir `<process>.name` in der Ausgabepfad-Deklaration referenzieren:

=== "Nachher"

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

Wir können weitergehen und komplexere Unterverzeichnispfade zusammenstellen.

In der obigen Bearbeitung haben wir die Unterscheidung zwischen `intermediates` gegenüber finalen Ausgaben auf der obersten Ebene aufgehoben.
Lass uns das zurückholen und die Dateien auch in ein `params.batch`-Unterverzeichnis legen.

!!! tip "Tipp"

    Das Einbeziehen von `params.batch` im Output-Block-`path` statt im `outputDir`-Config bedeutet, dass es nicht mit `-output-dir` auf der CLI überschrieben wird.

Aktualisiere zuerst die Config-Datei, um `${params.batch}` aus `outputDir` zu entfernen (da wir es in die Pfaddeklarationen verschieben):

=== "Nachher"

    ```groovy title="nextflow.config" linenums="12" hl_lines="4"
    /*
    * Ausgabe-Einstellungen
    */
    outputDir = "custom-outdir-config/"
    ```

=== "Vorher"

    ```groovy title="nextflow.config" linenums="12" hl_lines="4"
    /*
    * Ausgabe-Einstellungen
    */
    outputDir = "custom-outdir-config/${params.batch}"
    ```

Nimm dann folgende Änderungen in der Workflow-Datei vor:

=== "Nachher"

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

#### 2.2.3. Die Pipeline ausführen

Lass uns sehen, wie das in der Praxis funktioniert, indem wir sowohl `-output-dir` (oder `-o` kurz) auf `custom-outdir-config-2` als auch den Batch-Namen auf `rep2` von der Befehlszeile aus setzen:

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

### 2.3. Den Veröffentlichungsmodus auf Workflow-Ebene festlegen

Schließlich können wir im Sinne der Reduzierung von repetitivem Code die pro-Ausgabe-`mode`-Deklarationen durch eine einzelne Zeile in der Konfiguration ersetzen.

#### 2.3.1. `workflow.output.mode` zur Konfigurationsdatei hinzufügen

Füge folgenden Code zur `nextflow.config`-Datei hinzu:

=== "Nachher"

    ```groovy title="nextflow.config" linenums="12" hl_lines="5"
    /*
    * Ausgabe-Einstellungen
    */
    outputDir = "custom-outdir-config/"
    workflow.output.mode = 'copy'
    ```

=== "Vorher"

    ```groovy title="nextflow.config" linenums="12"
    /*
    * Ausgabe-Einstellungen
    */
    outputDir = "custom-outdir-config/"
    ```

Das Setzen von `workflow.output.mode` in der Konfigurationsdatei reicht aus, um das zu überschreiben, was in der Workflow-Datei gesetzt ist, aber entfernen wir den unnötigen Code trotzdem.

#### 2.3.2. Den Ausgabemodus aus der Workflow-Datei entfernen

Nimm folgende Änderungen in der Workflow-Datei vor:

=== "Nachher"

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

Das ist doch prägnanter, oder?

#### 2.3.3. Die Pipeline ausführen

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

Dies veröffentlicht Ausgaben nach `config-output-mode/`, und sie sind alle noch echte Kopien, keine Symlinks.

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

Der Hauptgrund, warum du möglicherweise immer noch die pro-Ausgabe-Methode zum Festlegen des Modus verwenden möchtest, ist, wenn du innerhalb desselben Workflows mischen möchtest, _d.h._ einige Ausgaben kopieren und einige als Symlinks haben möchtest.

Es gibt viele andere Optionen, die du auf diese Weise anpassen kannst, aber hoffentlich gibt dir dies einen Eindruck von der Bandbreite der Optionen und wie du sie effektiv nutzen kannst, um deinen Präferenzen gerecht zu werden.

### Fazit

Du weißt, wie du die Benennung und Struktur der Verzeichnisse steuerst, in denen deine Ausgaben veröffentlicht werden, sowie den Workflow-Ausgabe-Veröffentlichungsmodus.

### Wie geht es weiter?

Lerne, wie du deine Workflow-Konfiguration an deine Rechenumgebung anpasst, beginnend mit der Software-Paketierungstechnologie.

---

## 3. Eine Software-Paketierungstechnologie auswählen

Bisher haben wir uns Konfigurationselemente angesehen, die steuern, wie Eingaben hineingehen und wo Ausgaben herauskommen. Jetzt konzentrieren wir uns speziell darauf, deine Workflow-Konfiguration an deine Rechenumgebung anzupassen.

Der erste Schritt auf diesem Weg ist die Angabe, woher die Softwarepakete kommen, die in jedem Schritt ausgeführt werden.
Sind sie bereits in der lokalen Rechenumgebung installiert?
Müssen wir Images abrufen und über ein Container-System ausführen?
Oder müssen wir Conda-Pakete abrufen und eine lokale Conda-Umgebung erstellen?

Im allerersten Teil dieses Kurses (Teile 1-4) haben wir einfach lokal installierte Software in unserem Workflow verwendet.
Dann haben wir in Teil 5 Docker-Container und die `nextflow.config`-Datei eingeführt, die wir verwendet haben, um die Verwendung von Docker-Containern zu aktivieren.

Lass uns jetzt sehen, wie wir eine alternative Software-Paketierungsoption über die `nextflow.config`-Datei konfigurieren können.

### 3.1. Docker deaktivieren und Conda in der Config-Datei aktivieren

Stellen wir uns vor, wir arbeiten auf einem HPC-Cluster und die Administration erlaubt die Verwendung von Docker aus Sicherheitsgründen nicht.
Glücklicherweise unterstützt Nextflow mehrere andere Container-Technologien wie Singularity (das auf HPC weiter verbreitet ist) sowie Software-Paketmanager wie Conda.

Wir können unsere Konfigurationsdatei ändern, um [Conda](https://nextflow.io/docs/latest/conda.html) anstelle von Docker zu verwenden.
Dazu ändern wir den Wert von `docker.enabled` auf `false` und fügen eine Direktive hinzu, die die Verwendung von Conda aktiviert:

=== "Nachher"

    ```groovy title="nextflow.config" linenums="1" hl_lines="1-2"
    docker.enabled = false
    conda.enabled = true
    ```

=== "Vorher"

    ```groovy title="nextflow.config" linenums="1" hl_lines="1"
    docker.enabled = true
    ```

Dies ermöglicht es Nextflow, Conda-Umgebungen für Prozesse zu erstellen und zu nutzen, die Conda-Pakete angegeben haben.
Das bedeutet, dass wir jetzt eines davon zu unserem `cowpy`-Prozess hinzufügen müssen!

### 3.2. Ein Conda-Paket in der Prozessdefinition angeben

Wir haben bereits die URI für ein Conda-Paket abgerufen, das das `cowpy`-Tool enthält: `conda-forge::cowpy==1.1.5`

Jetzt fügen wir die URI zur `cowpy`-Prozessdefinition hinzu, indem wir die `conda`-Direktive verwenden:

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

Um klar zu sein, wir _ersetzen_ nicht die `container`-Direktive, wir _fügen_ eine alternative Option hinzu.

!!! tip "Tipp"

    Es gibt einige verschiedene Möglichkeiten, die URI für ein bestimmtes Conda-Paket zu erhalten.
    Wir empfehlen die Verwendung der [Seqera Containers](https://seqera.io/containers/)-Suchabfrage, die dir eine URI gibt, die du kopieren und einfügen kannst, auch wenn du nicht vorhast, einen Container daraus zu erstellen.

### 3.3. Den Workflow ausführen, um zu überprüfen, dass er Conda verwenden kann

Probieren wir es aus.

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

Dies sollte ohne Probleme funktionieren und dieselben Ausgaben wie zuvor unter `custom-outdir-config/conda` produzieren.

Hinter den Kulissen hat Nextflow die Conda-Pakete abgerufen und die Umgebung erstellt, was normalerweise etwas Arbeit erfordert; es ist also schön, dass wir das nicht selbst tun müssen!

!!! note "Hinweis"

    Dies läuft schnell, weil das `cowpy`-Paket ziemlich klein ist, aber wenn du mit großen Paketen arbeitest, kann es beim ersten Mal etwas länger dauern als üblich, und du könntest sehen, dass die Konsolenausgabe für eine Minute oder so 'hängen bleibt', bevor sie abgeschlossen wird.
    Das ist normal und liegt an der zusätzlichen Arbeit, die Nextflow beim ersten Mal erledigt, wenn du ein neues Paket verwendest.

Aus unserer Sicht sieht es so aus, als würde es genauso funktionieren wie mit Docker, obwohl die Mechanik im Backend etwas anders ist.

Das bedeutet, dass wir bereit sind, mit Conda-Umgebungen zu arbeiten, wenn nötig.

??? info "Docker und Conda mischen"

    Da diese Direktiven pro Prozess zugewiesen werden, ist es möglich, 'zu mischen', _d.h._ einige Prozesse in deinem Workflow mit Docker und andere mit Conda auszuführen, zum Beispiel wenn die Recheninfrastruktur, die du verwendest, beides unterstützt.
    In diesem Fall würdest du sowohl Docker als auch Conda in deiner Konfigurationsdatei aktivieren.
    Wenn beide für einen bestimmten Prozess verfügbar sind, wird Nextflow Container priorisieren.

    Und wie bereits erwähnt, unterstützt Nextflow mehrere andere Software-Paketierungs- und Container-Technologien, du bist also nicht nur auf diese beiden beschränkt.

### Fazit

Du weißt, wie du konfigurierst, welches Softwarepaket jeder Prozess verwenden soll, und wie du zwischen Technologien wechseln kannst.

### Wie geht es weiter?

Lerne, wie du die Ausführungsplattform änderst, die Nextflow verwendet, um die eigentliche Arbeit zu erledigen.

---

## 4. Eine Ausführungsplattform auswählen

Bis jetzt haben wir unsere Pipeline mit dem lokalen Executor ausgeführt.
Dieser führt jede Aufgabe auf der Maschine aus, auf der Nextflow läuft.
Wenn Nextflow startet, schaut es sich die verfügbaren CPUs und den Speicher an.
Wenn die Ressourcen der zur Ausführung bereiten Aufgaben die verfügbaren Ressourcen übersteigen, hält Nextflow die letzten Aufgaben von der Ausführung zurück, bis eine oder mehrere der früheren Aufgaben abgeschlossen sind und die notwendigen Ressourcen freigeben.

Der lokale Executor ist bequem und effizient, aber er ist auf diese einzelne Maschine beschränkt. Bei sehr großen Arbeitslasten kann deine lokale Maschine zum Engpass werden, entweder weil du eine einzelne Aufgabe hast, die mehr Ressourcen benötigt als verfügbar, oder weil du so viele Aufgaben hast, dass das Warten auf eine einzelne Maschine zu lange dauern würde.

Nextflow unterstützt [viele verschiedene Executors](https://nextflow.io/docs/latest/executor.html), einschließlich HPC-Scheduler (Slurm, LSF, SGE, PBS, Moab, OAR, Bridge, HTCondor und andere) sowie Cloud-Ausführungs-Backends wie (AWS Batch, Google Cloud Batch, Azure Batch, Kubernetes und mehr).

### 4.1. Ein anderes Backend anvisieren

Die Wahl des Executors wird durch eine Prozess-Direktive namens `executor` festgelegt.
Standardmäßig ist sie auf `local` gesetzt, daher ist die folgende Konfiguration impliziert:

```groovy title="Eingebaute Konfiguration"
process {
    executor = 'local'
}
```

Um den Executor auf ein anderes Backend zu setzen, gibst du einfach den gewünschten Executor mit einer ähnlichen Syntax an, wie oben für Ressourcenzuweisungen beschrieben (siehe [Executor-Dokumentation](https://nextflow.io/docs/latest/executor.html) für alle Optionen).

```groovy title="nextflow.config"
process {
    executor = 'slurm'
}
```

!!! warning "Warnung"

    Wir können dies in der Trainingsumgebung nicht tatsächlich testen, weil sie nicht für die Verbindung mit einem HPC eingerichtet ist.

### 4.2. Mit backend-spezifischer Syntax für Ausführungsparameter umgehen

Die meisten Hochleistungsrechenplattformen erlauben (und erfordern manchmal), dass du bestimmte Parameter wie Ressourcenzuweisungsanfragen und -beschränkungen (z.B. Anzahl der CPUs und Speicher) und den Namen der zu verwendenden Job-Warteschlange angibst.

Leider verwendet jedes dieser Systeme unterschiedliche Technologien, Syntaxen und Konfigurationen, um zu definieren, wie ein Job definiert und an den entsprechenden Scheduler übermittelt werden soll.

??? abstract "Beispiele"

    Zum Beispiel muss derselbe Job, der 8 CPUs und 4GB RAM benötigt und in der Warteschlange "my-science-work" ausgeführt werden soll, je nach Backend auf unterschiedliche Weise ausgedrückt werden.

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

Glücklicherweise vereinfacht Nextflow all das.
Es bietet eine standardisierte Syntax, damit du die relevanten Eigenschaften wie [`cpus`](https://nextflow.io/docs/latest/reference/process.html#cpus), [`memory`](https://nextflow.io/docs/latest/reference/process.html#memory) und [`queue`](https://nextflow.io/docs/latest/reference/process.html#queue) (siehe [Prozess-Direktiven](https://nextflow.io/docs/latest/reference/process.html#process-directives) für andere Eigenschaften) nur einmal angeben kannst.
Dann wird Nextflow zur Laufzeit diese Einstellungen verwenden, um die entsprechenden backend-spezifischen Scripts basierend auf der Executor-Einstellung zu generieren.

Wir werden diese standardisierte Syntax im nächsten Abschnitt behandeln.

### Fazit

Du weißt jetzt, wie du den Executor änderst, um verschiedene Arten von Recheninfrastruktur zu nutzen.

### Wie geht es weiter?

Lerne, wie du Ressourcenzuweisungen und -beschränkungen in Nextflow evaluierst und ausdrückst.

---

## 5. Rechenressourcen-Zuweisungen steuern

Die meisten Hochleistungsrechenplattformen erlauben (und erfordern manchmal), dass du bestimmte Ressourcenzuweisungsparameter wie Anzahl der CPUs und Speicher angibst.

Standardmäßig verwendet Nextflow eine einzelne CPU und 2GB Speicher für jeden Prozess.
Die entsprechenden Prozess-Direktiven heißen `cpus` und `memory`, daher ist die folgende Konfiguration impliziert:

```groovy title="Eingebaute Konfiguration" linenums="1"
process {
    cpus = 1
    memory = 2.GB
}
```

Du kannst diese Werte ändern, entweder für alle Prozesse oder für bestimmte benannte Prozesse, indem du zusätzliche Prozess-Direktiven in deiner Konfigurationsdatei verwendest.
Nextflow wird sie in die entsprechenden Anweisungen für den gewählten Executor übersetzen.

Aber woher weißt du, welche Werte du verwenden sollst?

### 5.1. Den Workflow ausführen, um einen Ressourcennutzungsbericht zu generieren

Wenn du nicht im Voraus weißt, wie viel CPU und Speicher deine Prozesse wahrscheinlich benötigen werden, kannst du ein Ressourcen-Profiling durchführen, das heißt, du führst den Workflow mit einigen Standard-Zuweisungen aus, zeichnest auf, wie viel jeder Prozess verwendet hat, und schätzt von dort aus, wie du die Basis-Zuweisungen anpassen solltest.

Praktischerweise enthält Nextflow eingebaute Tools dafür und wird auf Anfrage gerne einen Bericht für dich erstellen.

Dazu füge `-with-report <filename>.html` zu deiner Befehlszeile hinzu.

```bash
nextflow run hello-config.nf -with-report report-config-1.html
```

Der Bericht ist eine HTML-Datei, die du herunterladen und in deinem Browser öffnen kannst. Du kannst auch mit der rechten Maustaste darauf klicken im Datei-Explorer auf der linken Seite und auf `Show preview` klicken, um ihn in der Trainingsumgebung anzuzeigen.

Nimm dir ein paar Minuten, um den Bericht durchzusehen und zu schauen, ob du einige Möglichkeiten zur Anpassung der Ressourcen identifizieren kannst.
Klicke unbedingt auf die Tabs, die die Nutzungsergebnisse als Prozentsatz dessen zeigen, was zugewiesen wurde.

Siehe [Reports](https://nextflow.io/docs/latest/reports.html) für Dokumentation über alle verfügbaren Funktionen.

### 5.2. Ressourcenzuweisungen für alle Prozesse festlegen

Das Profiling zeigt, dass die Prozesse in unserem Trainings-Workflow sehr leichtgewichtig sind, also reduzieren wir die Standard-Speicherzuweisung auf 1GB pro Prozess.

Füge Folgendes zu deiner `nextflow.config`-Datei hinzu, vor dem Pipeline-Parameter-Abschnitt:

=== "Nachher"

    ```groovy title="nextflow.config" linenums="1" hl_lines="4-9"
    docker.enabled = false
    conda.enabled = true

    /*
    * Prozess-Einstellungen
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

Das wird helfen, den Rechenverbrauch zu reduzieren.

### 5.3. Ressourcenzuweisungen für einen bestimmten Prozess festlegen

Gleichzeitig werden wir so tun, als ob der `cowpy`-Prozess mehr Ressourcen benötigt als die anderen, nur um zu demonstrieren, wie man Zuweisungen für einen einzelnen Prozess anpasst.

=== "Nachher"

    ```groovy title="nextflow.config" linenums="4" hl_lines="6-9"
    /*
    * Prozess-Einstellungen
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
    * Prozess-Einstellungen
    */
    process {
        memory = 1.GB
    }
    ```

Mit dieser Konfiguration werden alle Prozesse 1GB Speicher und eine einzelne CPU (der implizierte Standard) anfordern, außer dem `cowpy`-Prozess, der 2GB und 2 CPUs anfordern wird.

!!! tip "Tipp"

    Wenn du eine Maschine mit wenigen CPUs hast und eine hohe Anzahl pro Prozess zuweist, könntest du sehen, dass Prozessaufrufe hintereinander in die Warteschlange gestellt werden.
    Das liegt daran, dass Nextflow sicherstellt, dass wir nicht mehr CPUs anfordern als verfügbar sind.

### 5.4. Den Workflow mit der aktualisierten Konfiguration ausführen

Probieren wir das aus, indem wir einen anderen Dateinamen für den Profiling-Bericht angeben, damit wir die Leistung vor und nach den Konfigurationsänderungen vergleichen können.

```bash
nextflow run hello-config.nf -with-report report-config-2.html
```

Du wirst wahrscheinlich keinen echten Unterschied bemerken, da dies eine so kleine Arbeitslast ist, aber das ist der Ansatz, den du verwenden würdest, um die Leistung und Ressourcenanforderungen eines realen Workflows zu analysieren.

Es ist sehr nützlich, wenn deine Prozesse unterschiedliche Ressourcenanforderungen haben. Es ermöglicht dir, die Ressourcenzuweisungen, die du für jeden Prozess einrichtest, basierend auf tatsächlichen Daten richtig zu dimensionieren, nicht auf Vermutungen.

!!! tip "Tipp"

    Dies ist nur ein kleiner Vorgeschmack dessen, was du tun kannst, um deine Ressourcennutzung zu optimieren.
    Nextflow selbst hat eine wirklich raffinierte [dynamische Wiederholungslogik](https://nextflow.io/docs/latest/process.html#dynamic-task-resources) eingebaut, um Jobs, die aufgrund von Ressourcenbeschränkungen fehlschlagen, automatisch zu wiederholen.
    Zusätzlich bietet die Seqera Platform KI-gestützte Tools zur automatischen Optimierung deiner Ressourcenzuweisungen.

### 5.5. Ressourcenlimits hinzufügen

Abhängig davon, welchen Rechen-Executor und welche Recheninfrastruktur du verwendest, kann es Einschränkungen geben, was du zuweisen kannst (oder musst).
Zum Beispiel kann dein Cluster erfordern, dass du innerhalb bestimmter Grenzen bleibst.

Du kannst die `resourceLimits`-Direktive verwenden, um die relevanten Beschränkungen festzulegen. Die Syntax sieht so aus, wenn sie allein in einem Process-Block steht:

```groovy title="Syntax-Beispiel"
process {
    resourceLimits = [
        memory: 750.GB,
        cpus: 200,
        time: 30.d
    ]
}
```

Nextflow wird diese Werte in die entsprechenden Anweisungen übersetzen, abhängig vom Executor, den du angegeben hast.

Wir werden das nicht ausführen, da wir in der Trainingsumgebung keinen Zugang zu relevanter Infrastruktur haben.
Wenn du jedoch versuchen würdest, den Workflow mit Ressourcenzuweisungen auszuführen, die diese Limits überschreiten, und dann den `sbatch`-Befehl in der `.command.run`-Script-Datei nachschlagen würdest, würdest du sehen, dass die Anfragen, die tatsächlich an den Executor gesendet werden, bei den durch `resourceLimits` angegebenen Werten gedeckelt sind.

??? info "Institutionelle Referenzkonfigurationen"

    Das nf-core-Projekt hat eine [Sammlung von Konfigurationsdateien](https://nf-co.re/configs/) zusammengestellt, die von verschiedenen Institutionen weltweit geteilt werden und eine breite Palette von HPC- und Cloud-Executors abdecken.

    Diese geteilten Configs sind sowohl für Personen wertvoll, die dort arbeiten und daher einfach die Konfiguration ihrer Institution direkt nutzen können, als auch als Modell für Personen, die eine Konfiguration für ihre eigene Infrastruktur entwickeln möchten.

### Fazit

Du weißt, wie du einen Profiling-Bericht erstellst, um die Ressourcennutzung zu bewerten, und wie du Ressourcenzuweisungen für alle Prozesse und/oder für einzelne Prozesse änderst sowie Ressourcenbeschränkungen für die Ausführung auf HPC festlegst.

### Wie geht es weiter?

Lerne, wie du voreingestellte Konfigurationsprofile einrichtest und zur Laufzeit zwischen ihnen wechselst.

---

## 6. Profile verwenden, um zwischen voreingestellten Konfigurationen zu wechseln

Wir haben dir eine Reihe von Möglichkeiten gezeigt, wie du deine Pipeline-Konfiguration anpassen kannst, abhängig vom Projekt, an dem du arbeitest, oder der Rechenumgebung, die du verwendest.

Möglicherweise möchtest du zwischen alternativen Einstellungen wechseln, abhängig davon, welche Recheninfrastruktur du verwendest. Zum Beispiel könntest du lokal auf deinem Laptop entwickeln und kleine Tests durchführen, dann vollständige Arbeitslasten auf HPC oder Cloud ausführen.

Nextflow ermöglicht es dir, beliebig viele [Profile](https://nextflow.io/docs/latest/config.html#config-profiles) einzurichten, die verschiedene Konfigurationen beschreiben, die du dann zur Laufzeit mit einem Befehlszeilenargument auswählen kannst, anstatt die Konfigurationsdatei selbst ändern zu müssen.

### 6.1. Profile für den Wechsel zwischen lokaler Entwicklung und Ausführung auf HPC erstellen

Lass uns zwei alternative Profile einrichten; eines für die Ausführung kleiner Lasten auf einem normalen Computer, wo wir Docker-Container verwenden werden, und eines für die Ausführung auf einem Universitäts-HPC mit einem Slurm-Scheduler, wo wir Conda-Pakete verwenden werden.

#### 6.1.1. Die Profile einrichten

Füge Folgendes zu deiner `nextflow.config`-Datei hinzu, nach dem Pipeline-Parameter-Abschnitt, aber vor den Ausgabe-Einstellungen:

=== "Nachher"

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
    * Ausgabe-Einstellungen
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
    * Ausgabe-Einstellungen
    */
    outputDir = "custom-outdir-config/"
    workflow.output.mode = 'copy'
    ```

Du siehst, dass wir für das Universitäts-HPC auch Ressourcenbeschränkungen angeben.

#### 6.1.2. Den Workflow mit einem Profil ausführen

Um ein Profil in unserer Nextflow-Befehlszeile anzugeben, verwenden wir das Argument `-profile`.

Versuchen wir, den Workflow mit der `my_laptop`-Konfiguration auszuführen.

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

Wie du sehen kannst, ermöglicht uns das, zur Laufzeit sehr bequem zwischen Konfigurationen zu wechseln.

!!! warning "Warnung"

    Das `univ_hpc`-Profil wird in der Trainingsumgebung nicht richtig funktionieren, da wir keinen Zugang zu einem Slurm-Scheduler haben.

Wenn wir in Zukunft andere Konfigurationselemente finden, die immer mit diesen zusammen auftreten, können wir sie einfach zu den entsprechenden Profilen hinzufügen.
Wir können auch zusätzliche Profile erstellen, wenn es andere Konfigurationselemente gibt, die wir zusammenfassen möchten.

### 6.2. Ein Profil mit Testparametern erstellen

Profile sind nicht nur für Infrastrukturkonfiguration.
Wir können sie auch verwenden, um Standardwerte für Workflow-Parameter festzulegen, um es anderen zu erleichtern, den Workflow auszuprobieren, ohne selbst geeignete Eingabewerte sammeln zu müssen.
Du kannst dies als Alternative zur Verwendung einer Parameterdatei betrachten.

#### 6.2.1. Das Profil einrichten

Die Syntax zum Ausdrücken von Standardwerten in diesem Kontext sieht so aus, für ein Profil, das wir `test` nennen:

```groovy title="Syntax-Beispiel"
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

Genau wie bei technischen Konfigurationsprofilen kannst du mehrere verschiedene Profile einrichten, die Parameter unter beliebigen Namen angeben.

#### 6.2.2. Den Workflow lokal mit dem Testprofil ausführen

Praktischerweise schließen sich Profile nicht gegenseitig aus, sodass wir mehrere Profile in unserer Befehlszeile mit der folgenden Syntax angeben können: `-profile <profile1>,<profile2>` (für beliebig viele Profile).

Wenn du Profile kombinierst, die Werte für dieselben Konfigurationselemente festlegen und in derselben Konfigurationsdatei beschrieben sind, wird Nextflow den Konflikt lösen, indem es den Wert verwendet, den es zuletzt gelesen hat (_d.h._ was später in der Datei kommt).
Wenn die widersprüchlichen Einstellungen in verschiedenen Konfigurationsquellen festgelegt sind, gilt die Standard-[Vorrangordnung](https://nextflow.io/docs/latest/config.html).

Versuchen wir, das Testprofil zu unserem vorherigen Befehl hinzuzufügen:

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

Dies wird Docker wo möglich verwenden und Ausgaben unter `custom-outdir-config/test` produzieren, und diesmal ist der Charakter das komische Duo `dragonandcow`.

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

Das bedeutet, dass wir, solange wir Testdatendateien mit dem Workflow-Code verteilen, jedem ermöglichen können, den Workflow schnell auszuprobieren, ohne eigene Eingaben über die Befehlszeile oder eine Parameterdatei bereitstellen zu müssen.

!!! tip "Tipp"

    Wir können auf URLs für größere Dateien verweisen, die extern gespeichert sind.
    Nextflow wird sie automatisch herunterladen, solange eine offene Verbindung besteht.

    Weitere Details findest du im Side Quest [Arbeiten mit Dateien](../side_quests/working_with_files.md)

### 6.3. `nextflow config` verwenden, um die aufgelöste Konfiguration zu sehen

Wie oben erwähnt, kann derselbe Parameter manchmal in Profilen, die du kombinieren möchtest, auf unterschiedliche Werte gesetzt werden.
Und allgemeiner gibt es zahlreiche Orte, an denen Konfigurationselemente gespeichert werden können, und manchmal können dieselben Eigenschaften an verschiedenen Orten auf unterschiedliche Werte gesetzt werden.

Nextflow wendet eine festgelegte [Vorrangordnung](https://nextflow.io/docs/latest/config.html) an, um Konflikte zu lösen, aber das kann schwierig selbst zu bestimmen sein.
Und selbst wenn nichts in Konflikt steht, kann es mühsam sein, alle möglichen Orte nachzuschlagen, an denen Dinge konfiguriert werden könnten.

Glücklicherweise enthält Nextflow ein praktisches Dienstprogramm-Tool namens `config`, das diesen gesamten Prozess für dich automatisieren kann.

Das `config`-Tool durchsucht alle Inhalte in deinem aktuellen Arbeitsverzeichnis, sammelt alle Konfigurationsdateien ein und produziert die vollständig aufgelöste Konfiguration, die Nextflow verwenden würde, um den Workflow auszuführen.
Dadurch kannst du herausfinden, welche Einstellungen verwendet werden, ohne etwas starten zu müssen.

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

    outputDir = 'custom-outdir-config/'

    workflow {
      output {
          mode = 'copy'
      }
    }
    ```

Dies zeigt dir die Basis-Konfiguration, die du erhältst, wenn du nichts Zusätzliches in der Befehlszeile angibst.

#### 6.3.2. Die Konfiguration mit aktivierten spezifischen Einstellungen auflösen

Wenn du Befehlszeilenparameter bereitstellst, z.B. ein oder mehrere Profile aktivierst oder eine Parameterdatei lädst, wird der Befehl diese zusätzlich berücksichtigen.

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

Das wird besonders nützlich für komplexe Projekte, die mehrere Konfigurationsebenen beinhalten.

### Fazit

Du weißt, wie du Profile verwendest, um eine voreingestellte Konfiguration zur Laufzeit mit minimalem Aufwand auszuwählen.
Allgemeiner weißt du, wie du deine Workflow-Ausführungen konfigurierst, um verschiedenen Rechenplattformen gerecht zu werden und die Reproduzierbarkeit deiner Analysen zu verbessern.

### Wie geht es weiter?

Feiere und klopfe dir kräftig auf die Schulter! Du hast deinen allerersten Nextflow-Entwicklerkurs abgeschlossen.

Gehe zur abschließenden [Kurszusammenfassung](./next_steps.md) weiter, um zu überprüfen, was du gelernt hast, und herauszufinden, was als Nächstes kommt.

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
Was hat Vorrang, wenn derselbe Parameter sowohl in der Config-Datei als auch in der Befehlszeile gesetzt wird?
- [ ] Der Wert in der Config-Datei
- [x] Der Wert in der Befehlszeile
- [ ] Der zuerst angetroffene Wert
- [ ] Keiner; es verursacht einen Fehler

Mehr erfahren: [1.1. Standardwerte in `nextflow.config` verschieben](#11-standardwerte-in-nextflowconfig-verschieben)
</quiz>

<quiz>
Können sowohl Docker als auch Conda in derselben Konfiguration aktiviert sein?
- [x] Ja, Nextflow kann beide je nach Prozess-Direktiven verwenden
- [ ] Nein, nur eines kann zur Zeit aktiviert sein
- [ ] Ja, aber nur in Profilen
- [ ] Nein, sie schließen sich gegenseitig aus
</quiz>

<quiz>
Wenn sowohl Docker als auch Conda aktiviert sind und ein Prozess beide Direktiven hat, was wird priorisiert?
- [x] Docker (Container)
- [ ] Conda
- [ ] Das zuerst definierte
- [ ] Es verursacht einen Fehler

Mehr erfahren: [3. Eine Software-Paketierungstechnologie auswählen](#3-eine-software-paketierungstechnologie-auswahlen)
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

Mehr erfahren: [5.3. Ressourcenzuweisungen für einen bestimmten Prozess festlegen](#53-ressourcenzuweisungen-fur-einen-bestimmten-prozess-festlegen)
</quiz>

<quiz>
Welche Befehlszeilenoption generiert einen Ressourcennutzungsbericht?
- [ ] `-with-metrics`
- [ ] `-with-stats`
- [x] `-with-report`
- [ ] `-with-profile`

Mehr erfahren: [5.1. Den Workflow ausführen, um einen Ressourcennutzungsbericht zu generieren](#51-den-workflow-ausfuhren-um-einen-ressourcennutzungsbericht-zu-generieren)
</quiz>

<quiz>
Was macht die `resourceLimits`-Direktive?
- [ ] Setzt Mindest-Ressourcenanforderungen
- [ ] Weist Prozessen Ressourcen zu
- [x] Deckelt die maximalen Ressourcen, die angefordert werden können
- [ ] Überwacht die Ressourcennutzung

Mehr erfahren: [5.5. Ressourcenlimits hinzufügen](#55-ressourcenlimits-hinzufugen)
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
Wie gibst du eine Parameterdatei beim Ausführen von Nextflow an?
- [ ] `--params params.json`
- [ ] `-config params.json`
- [x] `-params-file params.json`
- [ ] `--input params.json`

Mehr erfahren: [1.3. Eine Parameterdatei verwenden](#13-eine-parameterdatei-verwenden)
</quiz>

<quiz>
Wofür können Profile verwendet werden? (Wähle alle zutreffenden aus)
- [x] Definieren infrastrukturspezifischer Einstellungen
- [x] Festlegen von Ressourcenlimits für verschiedene Umgebungen
- [x] Bereitstellen von Testparametern
- [ ] Definieren neuer Prozesse

Mehr erfahren: [6. Profile verwenden, um zwischen voreingestellten Konfigurationen zu wechseln](#6-profile-verwenden-um-zwischen-voreingestellten-konfigurationen-zu-wechseln)
</quiz>

<quiz>
Wie gibst du mehrere Profile in einem einzigen Befehl an?
- [ ] `-profile profile1 -profile profile2`
- [ ] `-profiles profile1,profile2`
- [x] `-profile profile1,profile2`
- [ ] `--profile profile1 --profile profile2`

Mehr erfahren: [6. Profile verwenden, um zwischen voreingestellten Konfigurationen zu wechseln](#6-profile-verwenden-um-zwischen-voreingestellten-konfigurationen-zu-wechseln)
</quiz>
