# Teil 5: Hello Containers

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } KI-gestützte Übersetzung - [mehr erfahren & Verbesserungen vorschlagen](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/Xqr--bKEN9U?si=y8lAedhEHWaTV4zd&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

/// caption
:fontawesome-brands-youtube:{ .youtube } Siehe [die gesamte Playlist](https://youtube.com/playlist?list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&si=eF7cLR62goy-lc6n) auf dem Nextflow YouTube-Kanal.

:green_book: Das Video-Transkript ist [hier](./transcripts/05_hello_containers.md) verfügbar.
///

In den Teilen 1-4 dieses Trainings hast du gelernt, wie du die grundlegenden Bausteine von Nextflow verwendest, um einen einfachen Workflow zu erstellen, der Text verarbeiten, die Ausführung bei mehreren Eingaben parallelisieren und die Ergebnisse für die weitere Verarbeitung sammeln kann.

Allerdings warst du auf grundlegende UNIX-Tools beschränkt, die in deiner Umgebung verfügbar sind.
Aufgaben aus der Praxis erfordern oft verschiedene Tools und Pakete, die nicht standardmäßig enthalten sind.
Normalerweise müsstest du diese Tools installieren, ihre Abhängigkeiten verwalten und eventuelle Konflikte lösen.

Das ist alles sehr mühsam und nervig, deshalb zeigen wir dir, wie du **Container** verwendest, um dieses Problem viel bequemer zu lösen.

Ein **Container** ist eine leichtgewichtige, eigenständige, ausführbare Softwareeinheit, die aus einem Container-**Image** erstellt wird und alles enthält, was zum Ausführen einer Anwendung benötigt wird, einschließlich Code, Systembibliotheken und Einstellungen.
Wie du dir vorstellen kannst, wird das sehr hilfreich sein, um deine Pipelines reproduzierbarer zu machen.

Beachte, dass wir dies mit [Docker](https://www.docker.com/get-started/) lehren, aber denk daran, dass Nextflow auch [mehrere andere Container-Technologien](https://nextflow.io/docs/latest/container.html) unterstützt.

??? info "Wie du von diesem Abschnitt aus beginnst"

    Dieser Abschnitt des Kurses setzt voraus, dass du die Teile 1-4 des [Hello Nextflow](./index.md)-Kurses abgeschlossen hast und eine vollständige funktionierende Pipeline besitzt.

    Wenn du den Kurs von diesem Punkt aus beginnst, musst du das `modules`-Verzeichnis aus den Lösungen kopieren:

    ```bash
    cp -r solutions/4-hello-modules/modules .
    ```

---

## 0. Aufwärmen: Führe `hello-containers.nf` aus

Wir werden das Workflow-Skript `hello-containers.nf` als Ausgangspunkt verwenden.
Es entspricht dem Skript, das durch die Bearbeitung von Teil 4 dieses Trainings erstellt wurde, außer dass wir die Ausgabeziele geändert haben:

```groovy title="hello-containers.nf" linenums="37" hl_lines="3 7 11 15"
output {
    first_output {
        path 'hello_containers'
        mode 'copy'
    }
    uppercased {
        path 'hello_containers'
        mode 'copy'
    }
    collected {
        path 'hello_containers'
        mode 'copy'
    }
    batch_report {
        path 'hello_containers'
        mode 'copy'
    }
}
```

Um sicherzustellen, dass alles funktioniert, führe das Skript einmal aus, bevor du Änderungen vornimmst:

```bash
nextflow run hello-containers.nf
```

??? success "Befehlsausgabe"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-containers.nf` [nice_escher] DSL2 - revision: d5dfdc9872

    executor > local (7)
    [5a/ec1fa1] sayHello (2) [100%] 3 of 3 ✔
    [30/32b5b8] convertToUpper (3) [100%] 3 of 3 ✔
    [d3/be01bc] collectGreetings [100%] 1 of 1 ✔

    ```

Wie zuvor findest du die Ausgabedateien im Verzeichnis, das im `output`-Block angegeben ist (`results/hello_containers/`).

??? abstract "Verzeichnisinhalt"

    ```console
    results/hello_containers/
    ├── Bonjour-output.txt
    ├── COLLECTED-batch-output.txt
    ├── Hello-output.txt
    ├── Holà-output.txt
    ├── batch-report.txt
    ├── UPPER-Bonjour-output.txt
    ├── UPPER-Hello-output.txt
    └── UPPER-Holà-output.txt
    ```

Wenn das bei dir funktioniert hat, bist du bereit zu lernen, wie man Container verwendet.

---

## 1. Verwende einen Container 'manuell'

Was wir tun möchten, ist einen Schritt zu unserem Workflow hinzuzufügen, der einen Container für die Ausführung verwendet.

Allerdings werden wir zunächst einige grundlegende Konzepte und Operationen durchgehen, um dein Verständnis davon zu festigen, was Container sind, bevor wir sie in Nextflow verwenden.

### 1.1. Lade das Container-Image herunter

Um einen Container zu verwenden, lädst du normalerweise ein Container-Image aus einer Container-Registry herunter oder _pullst_ es, und führst dann das Container-Image aus, um eine Container-Instanz zu erstellen.

Die allgemeine Syntax lautet wie folgt:

```bash title="Syntax"
docker pull '<container>'
```

Der `docker pull`-Teil ist die Anweisung an das Container-System, ein Container-Image aus einem Repository zu pullen.

Der `'<container>'`-Teil ist die URI-Adresse des Container-Images.

Als Beispiel pullen wir ein Container-Image, das [cowpy](https://github.com/jeffbuttars/cowpy) enthält, eine Python-Implementierung eines Tools namens `cowsay`, das ASCII-Art generiert, um beliebige Texteingaben auf unterhaltsame Weise anzuzeigen.

```txt title="Example"
 ________________________
< Are we having fun yet? >
 ------------------------
    \                                  ___-------___
     \                             _-~~             ~~-_
      \                         _-~                    /~-_
             /^\__/^\         /~  \                   /    \
           /|  O|| O|        /      \_______________/        \
          | |___||__|      /       /                \          \
          |          \    /      /                    \          \
          |   (_______) /______/                        \_________ \
          |         / /         \                      /            \
           \         \^\\         \                  /               \     /
             \         ||           \______________/      _-_       //\__//
               \       ||------_-~~-_ ------------- \ --/~   ~\    || __/
                 ~-----||====/~     |==================|       |/~~~~~
                  (_(__/  ./     /                    \_\      \.
                         (_(___/                         \_____)_)
```

Es gibt verschiedene Repositories, in denen du veröffentlichte Container finden kannst.
Wir haben den [Seqera Containers](https://seqera.io/containers/)-Service verwendet, um dieses Docker-Container-Image aus dem `cowpy` Conda-Paket zu generieren: `'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'`.

Führe den vollständigen Pull-Befehl aus:

```bash
docker pull 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
```

??? success "Befehlsausgabe"

    ```console
    1.1.5--3db457ae1977a273: Pulling from library/cowpy
    dafa2b0c44d2: Pull complete
    dec6b097362e: Pull complete
    f88da01cff0b: Pull complete
    4f4fb700ef54: Pull complete
    92dc97a3ef36: Pull complete
    403f74b0f85e: Pull complete
    10b8c00c10a5: Pull complete
    17dc7ea432cc: Pull complete
    bb36d6c3110d: Pull complete
    0ea1a16bbe82: Pull complete
    030a47592a0a: Pull complete
    c23bdb422167: Pull complete
    e1686ff32a11: Pull complete
    Digest: sha256:1ebc0043e8cafa61203bf42d29fd05bd14e7b4298e5e8cf986504c15f5aa4160
    Status: Downloaded newer image for community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273
    community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273
    ```

Wenn du das Image noch nie heruntergeladen hast, kann dies eine Minute dauern.
Sobald es fertig ist, hast du eine lokale Kopie des Container-Images.

### 1.2. Verwende den Container, um `cowpy` als einmaligen Befehl auszuführen

Eine sehr häufige Art, wie Leute Container verwenden, ist, sie direkt auszuführen, _d.h._ nicht-interaktiv.
Das ist großartig für die Ausführung einmaliger Befehle.

Die allgemeine Syntax lautet wie folgt:

```bash title="Syntax"
docker run --rm '<container>' [tool command]
```

Der `docker run --rm '<container>'`-Teil ist die Anweisung an das Container-System, eine Container-Instanz aus einem Container-Image zu starten und einen Befehl darin auszuführen.
Das `--rm`-Flag weist das System an, die Container-Instanz herunterzufahren, nachdem der Befehl abgeschlossen wurde.

Die `[tool command]`-Syntax hängt vom Tool ab, das du verwendest, und davon, wie der Container eingerichtet ist.
Beginnen wir einfach mit `cowpy`.

Vollständig zusammengesetzt sieht der Container-Ausführungsbefehl so aus; führe ihn aus.

```bash
docker run --rm 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273' cowpy
```

??? success "Befehlsausgabe"

    ```console
    ______________________________________________________
    < Cowacter, eyes:default, tongue:False, thoughts:False >
    ------------------------------------------------------
        \   ^__^
          \  (oo)\_______
            (__)\       )\/\
              ||----w |
              ||     ||
    ```

Das System hat den Container gestartet, den `cowpy`-Befehl mit seinen Parametern ausgeführt, die Ausgabe an die Konsole gesendet und schließlich die Container-Instanz heruntergefahren.

### 1.3. Verwende den Container, um `cowpy` interaktiv auszuführen

Du kannst einen Container auch interaktiv ausführen, was dir eine Shell-Eingabeaufforderung innerhalb des Containers gibt und es dir ermöglicht, mit dem Befehl zu experimentieren.

#### 1.3.1. Starte den Container

Um interaktiv auszuführen, fügen wir einfach `-it` zum `docker run`-Befehl hinzu.
Optional können wir die Shell angeben, die wir innerhalb des Containers verwenden möchten, indem wir _z.B._ `/bin/bash` an den Befehl anhängen.

```bash
docker run --rm -it 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273' /bin/bash
```

Beachte, dass sich deine Eingabeaufforderung zu etwas wie `(base) root@b645838b3314:/tmp#` ändert, was anzeigt, dass du dich jetzt innerhalb des Containers befindest.

Du kannst dies überprüfen, indem du `ls /` ausführst, um den Verzeichnisinhalt vom Root des Dateisystems aufzulisten:

```bash
ls /
```

??? abstract "Befehlsausgabe"

    ```console
    bin  boot  dev  etc  home  lib  lib64  media  mnt  opt  proc  root  run  sbin  srv  sys  tmp  usr  var
    ```

Wir verwenden hier `ls` anstelle von `tree`, weil das `tree`-Dienstprogramm in diesem Container nicht verfügbar ist.
Du kannst sehen, dass das Dateisystem innerhalb des Containers sich vom Dateisystem auf deinem Host-System unterscheidet.

Eine Einschränkung dessen, was wir gerade getan haben, ist, dass der Container standardmäßig vollständig vom Host-System isoliert ist.
Das bedeutet, dass der Container nicht auf Dateien auf dem Host-System zugreifen kann, es sei denn, du erlaubst es ihm ausdrücklich.

Wir zeigen dir gleich, wie das geht.

#### 1.3.2. Führe den gewünschten Tool-Befehl aus

Jetzt, da du dich innerhalb des Containers befindest, kannst du den `cowpy`-Befehl direkt ausführen und ihm einige Parameter geben.
Zum Beispiel sagt die Tool-Dokumentation, dass wir den Charakter ('cowacter') mit `-c` ändern können.

```bash
cowpy "Hello Containers" -c tux
```

??? success "Befehlsausgabe"

    ```console
    __________________
    < Hello Containers >
    ------------------
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

Jetzt zeigt die Ausgabe den Linux-Pinguin Tux anstelle der Standard-Kuh, weil wir den Parameter `-c tux` angegeben haben.

Da du dich innerhalb des Containers befindest, kannst du den `cowpy`-Befehl so oft ausführen, wie du möchtest, und die Eingabeparameter variieren, ohne dich mit Docker-Befehlen herumschlagen zu müssen.

!!! Tip "Tipp"

    Verwende das '-c'-Flag, um einen anderen Charakter auszuwählen, einschließlich:
    `beavis`, `cheese`, `daemon`, `dragonandcow`, `ghostbusters`, `kitty`, `moose`, `milk`, `stegosaurus`, `turkey`, `turtle`, `tux`

Das ist nett. Was noch netter wäre, ist, wenn wir unsere `greetings.csv` als Eingabe dafür verwenden könnten.
Aber da wir keinen Zugriff auf das Dateisystem haben, können wir das nicht.

Lass uns das beheben.

#### 1.3.3. Verlasse den Container

Um den Container zu verlassen, kannst du `exit` an der Eingabeaufforderung eingeben oder die Tastenkombination ++ctrl+d++ verwenden.

```bash
exit
```

Deine Eingabeaufforderung sollte jetzt wieder so sein wie vor dem Start des Containers.

#### 1.3.4. Mounte Daten in den Container

Wie bereits erwähnt, ist der Container standardmäßig vom Host-System isoliert.

Um dem Container den Zugriff auf das Host-Dateisystem zu ermöglichen, kannst du ein **Volume** vom Host-System in den Container **mounten**, indem du die folgende Syntax verwendest:

```bash title="Syntax"
-v <outside_path>:<inside_path>
```

In unserem Fall wird `<outside_path>` das aktuelle Arbeitsverzeichnis sein, also können wir einfach einen Punkt (`.`) verwenden, und `<inside_path>` ist nur ein Alias, den wir uns ausdenken; nennen wir es `/my_project` (der innere Pfad muss absolut sein).

Um ein Volume zu mounten, ersetzen wir die Pfade und fügen das Volume-Mounting-Argument zum docker run-Befehl wie folgt hinzu:

```bash
docker run --rm -it -v .:/my_project 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273' /bin/bash
```

Dies mountet das aktuelle Arbeitsverzeichnis als Volume, das unter `/my_project` innerhalb des Containers zugänglich sein wird.

Du kannst überprüfen, ob es funktioniert, indem du den Inhalt von `/my_project` auflistest:

```bash
ls /my_project
```

??? success "Befehlsausgabe"

    ```console
    data               hello-config.nf      hello-modules.nf   hello-world.nf  nextflow.config  solutions         work
    hello-channels.nf  hello-containers.nf  hello-workflow.nf  modules         results          test-params.json
    ```

Du kannst jetzt den Inhalt des Arbeitsverzeichnisses von innerhalb des Containers sehen, einschließlich der `greetings.csv`-Datei unter `data/`.

Dies hat effektiv einen Tunnel durch die Container-Wand etabliert, den du verwenden kannst, um auf diesen Teil deines Dateisystems zuzugreifen.

#### 1.3.5. Verwende die gemounteten Daten

Jetzt, da wir das Arbeitsverzeichnis in den Container gemountet haben, können wir den `cowpy`-Befehl verwenden, um den Inhalt der `greetings.csv`-Datei anzuzeigen.

Dazu verwenden wir `cat /my_project/data/greetings.csv | `, um den Inhalt der CSV-Datei in den `cowpy`-Befehl zu pipen.

```bash
cat /my_project/data/greetings.csv | cowpy -c turkey
```

??? success "Befehlsausgabe"

    ```console title="data/greetings.csv"
     ____________________
    / Hello,English,123  \
    | Bonjour,French,456 |
    \ Holà,Spanish,789   /
    --------------------
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

Dies erzeugt die gewünschte ASCII-Art eines Truthahns, der unsere Beispiel-Begrüßungen herunterspult!
Außer dass der Truthahn hier die vollständigen Zeilen wiederholt, anstatt nur die Begrüßungen.
Wir wissen bereits, dass unser Nextflow-Workflow einen besseren Job machen wird!

Experimentiere ruhig mit diesem Befehl.
Wenn du fertig bist, verlasse den Container wie zuvor:

```bash
exit
```

Du befindest dich wieder in deiner normalen Shell.

### Fazit

Du weißt, wie man einen Container pullt und ihn entweder als einmaligen Befehl oder interaktiv ausführt. Du weißt auch, wie du deine Daten von innerhalb deines Containers zugänglich machst, was es dir ermöglicht, jedes Tool, an dem du interessiert bist, mit echten Daten auszuprobieren, ohne Software auf deinem System installieren zu müssen.

### Wie geht es weiter?

Lerne, wie man Container für die Ausführung von Nextflow-Prozessen verwendet.

---

## 2. Verwende Container in Nextflow

Nextflow hat eingebaute Unterstützung für die Ausführung von Prozessen innerhalb von Containern, damit du Tools ausführen kannst, die du nicht in deiner Rechenumgebung installiert hast.
Das bedeutet, dass du jedes beliebige Container-Image verwenden kannst, um deine Prozesse auszuführen, und Nextflow kümmert sich um das Pullen des Images, das Mounten der Daten und die Ausführung des Prozesses darin.

Um dies zu demonstrieren, werden wir einen `cowpy`-Schritt zur Pipeline hinzufügen, die wir entwickelt haben, nach dem `collectGreetings`-Schritt.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-pipeline-cowpy.svg"
</figure>

### 2.1. Schreibe ein `cowpy`-Modul

Erstellen wir zunächst das `cowpy`-Prozessmodul.

#### 2.1.1. Erstelle eine Datei-Vorlage für das neue Modul

Erstelle eine leere Datei für das Modul namens `cowpy.nf`.

```bash
touch modules/cowpy.nf
```

Dies gibt uns einen Ort, an dem wir den Prozesscode ablegen können.

#### 2.1.2. Kopiere den `cowpy`-Prozesscode in die Moduldatei

Wir können unseren `cowpy`-Prozess nach den anderen Prozessen modellieren, die wir zuvor geschrieben haben.

```groovy title="modules/cowpy.nf" linenums="1"
#!/usr/bin/env nextflow

// ASCII-Art mit cowpy generieren
process cowpy {

    input:
    path input_file
    val character

    output:
    path "cowpy-${input_file}"

    script:
    """
    cat ${input_file} | cowpy -c "${character}" > cowpy-${input_file}
    """

}
```

Der Prozess erwartet eine `input_file`, die die Begrüßungen enthält, sowie einen `character`-Wert.

Die Ausgabe wird eine neue Textdatei sein, die die vom `cowpy`-Tool generierte ASCII-Art enthält.

### 2.2. Füge cowpy zum Workflow hinzu

Jetzt müssen wir das Modul importieren und den Prozess aufrufen.

#### 2.2.1. Importiere den `cowpy`-Prozess in `hello-containers.nf`

Füge die Import-Deklaration über dem Workflow-Block ein und fülle sie entsprechend aus.

=== "Danach"

    ```groovy title="hello-containers.nf" linenums="3" hl_lines="5"
    // Module einbinden
    include { sayHello } from './modules/sayHello.nf'
    include { convertToUpper } from './modules/convertToUpper.nf'
    include { collectGreetings } from './modules/collectGreetings.nf'
    include { cowpy } from './modules/cowpy.nf'
    ```

=== "Vorher"

    ```groovy title="hello-containers.nf" linenums="3"
    // Module einbinden
    include { sayHello } from './modules/sayHello.nf'
    include { convertToUpper } from './modules/convertToUpper.nf'
    include { collectGreetings } from './modules/collectGreetings.nf'
    ```

Jetzt ist das `cowpy`-Modul im Workflow verfügbar.

#### 2.2.2. Füge einen Aufruf des `cowpy`-Prozesses im Workflow hinzu

Verbinden wir den `cowpy()`-Prozess mit der Ausgabe des `collectGreetings()`-Prozesses, der, wie du dich vielleicht erinnerst, zwei Ausgaben produziert:

- `collectGreetings.out.outfile` enthält die Ausgabedatei <--_was wir wollen_
- `collectGreetings.out.report` enthält die Berichtsdatei mit der Anzahl der Begrüßungen pro Batch

Nimm im Workflow-Block die folgende Codeänderung vor:

=== "Danach"

    ```groovy title="hello-containers.nf" linenums="19" hl_lines="12-13"
        main:
        // Einen Kanal für Eingaben aus einer CSV-Datei erstellen
        greeting_ch = channel.fromPath(params.input)
                            .splitCsv()
                            .map { line -> line[0] }
        // Eine Begrüßung ausgeben
        sayHello(greeting_ch)
        // Die Begrüßung in Großbuchstaben umwandeln
        convertToUpper(sayHello.out)
        // Alle Begrüßungen in einer Datei sammeln
        collectGreetings(convertToUpper.out.collect(), params.batch)
        // ASCII-Art der Begrüßungen mit cowpy generieren
        cowpy(collectGreetings.out.outfile, params.character)
    ```

=== "Vorher"

    ```groovy title="hello-containers.nf" linenums="19"
        main:
        // Einen Kanal für Eingaben aus einer CSV-Datei erstellen
        greeting_ch = channel.fromPath(params.input)
                            .splitCsv()
                            .map { line -> line[0] }
        // Eine Begrüßung ausgeben
        sayHello(greeting_ch)
        // Die Begrüßung in Großbuchstaben umwandeln
        convertToUpper(sayHello.out)
        // Alle Begrüßungen in einer Datei sammeln
        collectGreetings(convertToUpper.out.collect(), params.batch)
    ```

Beachte, dass wir einen neuen CLI-Parameter, `params.character`, deklariert haben, um anzugeben, welchen Charakter wir die Begrüßungen sagen lassen möchten.

#### 2.2.3. Füge den `character`-Parameter zum `params`-Block hinzu

Dies ist technisch optional, aber es ist die empfohlene Praxis und es ist eine Gelegenheit, einen Standardwert für den Charakter festzulegen.

=== "Danach"

    ```groovy title="hello-containers.nf" linenums="9" hl_lines="7"
    /*
    * Pipeline-Parameter
    */
    params {
        input: Path = 'data/greetings.csv'
        batch: String = 'batch'
        character: String = 'turkey'
    }
    ```

=== "Vorher"

    ```groovy title="hello-containers.nf" linenums="9"
    /*
    * Pipeline-Parameter
    */
    params {
        input: Path = 'data/greetings.csv'
        batch: String = 'batch'
    }
    ```

Jetzt können wir faul sein und das Eingeben des Charakter-Parameters in unseren Befehlszeilen überspringen.

#### 2.2.4. Aktualisiere die Workflow-Ausgaben

Wir müssen die Workflow-Ausgaben aktualisieren, um die Ausgabe des `cowpy`-Prozesses zu veröffentlichen.

##### 2.2.4.1. Aktualisiere den `publish:`-Abschnitt

Nimm im `workflow`-Block die folgende Codeänderung vor:

=== "Danach"

    ```groovy title="hello-containers.nf" linenums="34" hl_lines="6"
        publish:
        first_output = sayHello.out
        uppercased = convertToUpper.out
        collected = collectGreetings.out.outfile
        batch_report = collectGreetings.out.report
        cowpy_art = cowpy.out
    ```

=== "Vorher"

    ```groovy title="hello-containers.nf" linenums="34"
        publish:
        first_output = sayHello.out
        uppercased = convertToUpper.out
        collected = collectGreetings.out.outfile
        batch_report = collectGreetings.out.report
    ```

Der `cowpy`-Prozess produziert nur eine Ausgabe, also können wir auf die übliche Weise darauf verweisen, indem wir `.out` anhängen.

Aber lass uns zunächst die Aktualisierung der Workflow-Level-Ausgaben abschließen.

##### 2.2.4.2. Aktualisiere den `output`-Block

Wir müssen die finale `cowpy_art`-Ausgabe zum `output`-Block hinzufügen. Während wir dabei sind, bearbeiten wir auch die Veröffentlichungsziele, da unsere Pipeline jetzt vollständig ist und wir wissen, welche Ausgaben uns wirklich wichtig sind.

Nimm im `output`-Block die folgenden Codeänderungen vor:

=== "Danach"

    ```groovy title="hello-containers.nf" linenums="42" hl_lines="3 7 11 15 18-21"
    output {
        first_output {
            path 'hello_containers/intermediates'
            mode 'copy'
        }
        uppercased {
            path 'hello_containers/intermediates'
            mode 'copy'
        }
        collected {
            path 'hello_containers/intermediates'
            mode 'copy'
        }
        batch_report {
            path 'hello_containers'
            mode 'copy'
        }
        cowpy_art {
            path 'hello_containers'
            mode 'copy'
        }
    }
    ```

=== "Vorher"

    ```groovy title="hello-containers.nf" linenums="42" hl_lines="3 7 11 15"
    output {
        first_output {
            path 'hello_containers'
            mode 'copy'
        }
        uppercased {
            path 'hello_containers'
            mode 'copy'
        }
        collected {
            path 'hello_containers'
            mode 'copy'
        }
        batch_report {
            path 'hello_containers'
            mode 'copy'
        }
    }
    ```

Jetzt werden die veröffentlichten Ausgaben etwas organisierter sein.

#### 2.2.5. Führe den Workflow aus

Zur Zusammenfassung, das ist es, was wir anstreben:

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello_pipeline_complete.svg"
</figure>

Denkst du, es wird funktionieren?

Lass uns die vorherigen veröffentlichten Ausgaben löschen, um einen sauberen Ausgangspunkt zu haben, und den Workflow mit dem `-resume`-Flag ausführen.

```bash
rm -r hello_containers/
nextflow run hello-containers.nf -resume
```

??? failure "Befehlsausgabe (zur Klarheit bearbeitet)"

    ```console hl_lines="10 13 20-21 26-27"
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-containers.nf` [lonely_woese] DSL2 - revision: abf1dccf7f

    executor >  local (1)
    [c9/f5c686] sayHello (3)       [100%] 3 of 3, cached: 3 ✔
    [ef/3135a8] convertToUpper (3) [100%] 3 of 3, cached: 3 ✔
    [7f/f435e3] collectGreetings   [100%] 1 of 1, cached: 1 ✔
    [9b/02e776] cowpy              [  0%] 0 of 1 ✘
    ERROR ~ Error executing process > 'cowpy'

    Caused by:
      Process `cowpy` terminated with an error exit status (127)


    Command executed:

      cat COLLECTED-batch-output.txt | cowpy -c "turkey" > cowpy-COLLECTED-batch-output.txt

    Command exit status:
      127

    Command output:
      (empty)

    Command error:
      .command.sh: line 2: cowpy: command not found

    Work dir:
      /workspaces/training/hello-nextflow/work/9b/02e7761db848f82db3c3e59ff3a9b6

    Tip: when you have fixed the problem you can continue the execution adding the option `-resume` to the run command line

    -- Check '.nextflow.log' file for details
    ERROR ~ Cannot access first() element from an empty List

    -- Check '.nextflow.log' file for details
    ```

Oh nein, es gibt einen Fehler!
Der Fehlercode `error exit status (127)` bedeutet, dass die ausführbare Datei, die wir angefordert haben, nicht gefunden wurde.

Das macht Sinn, da wir das `cowpy`-Tool aufrufen, aber noch keinen Container angegeben haben (ups).

### 2.3. Verwende einen Container, um den `cowpy`-Prozess auszuführen

Wir müssen einen Container angeben und Nextflow mitteilen, dass er ihn für den `cowpy()`-Prozess verwenden soll.

#### 2.3.1. Gib einen Container für `cowpy` an

Wir können dasselbe Image verwenden, das wir direkt im ersten Abschnitt dieses Tutorials verwendet haben.

Bearbeite das `cowpy.nf`-Modul, um die `container`-Direktive zur Prozessdefinition wie folgt hinzuzufügen:

=== "Danach"

    ```groovy title="modules/cowpy.nf" linenums="4" hl_lines="3"
    process cowpy {

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'

        input:
        path input_file
        val character

        output:
        path "cowpy-${input_file}"

        script:
        """
        cat ${input_file} | cowpy -c "${character}" > cowpy-${input_file}
        """
    }
    ```

=== "Vorher"

    ```groovy title="modules/cowpy.nf" linenums="4"
    process cowpy {

        input:
        path input_file
        val character

        output:
        path "cowpy-${input_file}"

        script:
        """
        cat ${input_file} | cowpy -c "${character}" > cowpy-${input_file}
        """
    }
    ```

Dies teilt Nextflow mit, dass _wenn die Verwendung von Docker aktiviert ist_, es das hier angegebene Container-Image verwenden soll, um den Prozess auszuführen.

#### 2.3.2. Aktiviere die Verwendung von Docker über die `nextflow.config`-Datei

Beachte, dass wir sagten _'wenn die Verwendung von Docker aktiviert ist'_. Standardmäßig ist sie es nicht, also müssen wir Nextflow mitteilen, dass es Docker verwenden darf.
Zu diesem Zweck werden wir das Thema des nächsten und letzten Teils dieses Kurses (Teil 6), der die Konfiguration behandelt, leicht vorwegnehmen.

Eine der Hauptmöglichkeiten, die Nextflow für die Konfiguration der Workflow-Ausführung bietet, ist die Verwendung einer `nextflow.config`-Datei.
Wenn eine solche Datei im aktuellen Verzeichnis vorhanden ist, lädt Nextflow sie automatisch und wendet jede darin enthaltene Konfiguration an.

Wir haben eine `nextflow.config`-Datei mit einer einzigen Codezeile bereitgestellt, die Docker explizit deaktiviert: `docker.enabled = false`.

Lass uns das jetzt auf `true` ändern, um Docker zu aktivieren:

=== "Danach"

    ```console title="nextflow.config" linenums="1" hl_lines="1"
    docker.enabled = true
    ```

=== "Vorher"

    ```console title="nextflow.config" linenums="1" hl_lines="1"
    docker.enabled = false
    ```

!!! tip "Tipp"

    Es ist möglich, die Docker-Ausführung von der Befehlszeile aus zu aktivieren, auf Basis einzelner Ausführungen, mit dem Parameter `-with-docker <container>`.
    Allerdings erlaubt uns das nur, einen Container für den gesamten Workflow anzugeben, während der Ansatz, den wir dir gerade gezeigt haben, es uns ermöglicht, einen anderen Container pro Prozess anzugeben.
    Dies ist besser für Modularität, Code-Wartung und Reproduzierbarkeit.

#### 2.3.3. Führe den Workflow mit aktiviertem Docker aus

Führe den Workflow mit dem `-resume`-Flag aus:

```bash
nextflow run hello-containers.nf -resume
```

??? success "Befehlsausgabe"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-containers.nf` [drunk_perlman] DSL2 - revision: abf1dccf7f

    executor >  local (1)
    [c9/f5c686] sayHello (3)       [100%] 3 of 3, cached: 3 ✔
    [ef/3135a8] convertToUpper (3) [100%] 3 of 3, cached: 3 ✔
    [7f/f435e3] collectGreetings   [100%] 1 of 1, cached: 1 ✔
    [98/656c6c] cowpy              [100%] 1 of 1 ✔
    ```

Diesmal funktioniert es tatsächlich!
Wie üblich findest du die Workflow-Ausgaben im entsprechenden Ergebnisverzeichnis, obwohl sie diesmal etwas ordentlicher organisiert sind, mit nur dem Bericht und der finalen Ausgabe auf der obersten Ebene und allen Zwischendateien aus dem Weg in einem Unterverzeichnis.

??? abstract "Verzeichnisinhalt"

    ```console
    results/hello_containers/
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

Die finale ASCII-Art-Ausgabe befindet sich im Verzeichnis `results/hello_containers/` unter dem Namen `cowpy-COLLECTED-batch-output.txt`.

??? abstract "Dateiinhalt"

    ```console title="results/hello_containers/cowpy-COLLECTED-batch-output.txt"
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

Und da ist er, unser schöner Truthahn, der die Begrüßungen wie gewünscht sagt.

#### 2.3.4. Untersuche, wie Nextflow die containerisierte Aufgabe gestartet hat

Als finales Schlusswort zu diesem Abschnitt werfen wir einen Blick auf das Arbeitsunterverzeichnis für einen der `cowpy`-Prozessaufrufe, um etwas mehr Einblick zu bekommen, wie Nextflow unter der Haube mit Containern arbeitet.

Überprüfe die Ausgabe deines `nextflow run`-Befehls, um den Pfad zum Arbeitsunterverzeichnis für den `cowpy`-Prozess zu finden.
Wenn wir uns ansehen, was wir für die oben gezeigte Ausführung erhalten haben, beginnt die Konsolenprotokollzeile für den `cowpy`-Prozess mit `[98/656c6c]`.
Das entspricht dem folgenden gekürzten Verzeichnispfad: `work/98/656c6c`.

In diesem Verzeichnis findest du die `.command.run`-Datei, die alle Befehle enthält, die Nextflow in deinem Auftrag im Verlauf der Ausführung der Pipeline ausgeführt hat.

??? abstract "Dateiinhalt"

    ```console title="work/98/656c6c90cce1667c094d880f4b6dcc/.command.run"
    #!/bin/bash
    ### ---
    ### name: 'cowpy'
    ### container: 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
    ### outputs:
    ### - 'cowpy-COLLECTED-batch-output.txt'
    ### ...
    set -e
    set -u
    NXF_DEBUG=${NXF_DEBUG:=0}; [[ $NXF_DEBUG > 1 ]] && set -x
    NXF_ENTRY=${1:-nxf_main}


    nxf_sleep() {
      sleep $1 2>/dev/null || sleep 1;
    }

    nxf_date() {
        local ts=$(date +%s%3N);
        if [[ ${#ts} == 10 ]]; then echo ${ts}000
        elif [[ $ts == *%3N ]]; then echo ${ts/\%3N/000}
        elif [[ $ts == *3N ]]; then echo ${ts/3N/000}
        elif [[ ${#ts} == 13 ]]; then echo $ts
        else echo "Unexpected timestamp value: $ts"; exit 1
        fi
    }

    nxf_env() {
        echo '============= task environment ============='
        env | sort | sed "s/\(.*\)AWS\(.*\)=\(.\{6\}\).*/\1AWS\2=\3xxxxxxxxxxxxx/"
        echo '============= task output =================='
    }

    nxf_kill() {
        declare -a children
        while read P PP;do
            children[$PP]+=" $P"
        done < <(ps -e -o pid= -o ppid=)

        kill_all() {
            [[ $1 != $$ ]] && kill $1 2>/dev/null || true
            for i in ${children[$1]:=}; do kill_all $i; done
        }

        kill_all $1
    }

    nxf_mktemp() {
        local base=${1:-/tmp}
        mkdir -p "$base"
        if [[ $(uname) = Darwin ]]; then mktemp -d $base/nxf.XXXXXXXXXX
        else TMPDIR="$base" mktemp -d -t nxf.XXXXXXXXXX
        fi
    }

    nxf_fs_copy() {
      local source=$1
      local target=$2
      local basedir=$(dirname $1)
      mkdir -p $target/$basedir
      cp -fRL $source $target/$basedir
    }

    nxf_fs_move() {
      local source=$1
      local target=$2
      local basedir=$(dirname $1)
      mkdir -p $target/$basedir
      mv -f $source $target/$basedir
    }

    nxf_fs_rsync() {
      rsync -rRl $1 $2
    }

    nxf_fs_rclone() {
      rclone copyto $1 $2/$1
    }

    nxf_fs_fcp() {
      fcp $1 $2/$1
    }

    on_exit() {
        local last_err=$?
        local exit_status=${nxf_main_ret:=0}
        [[ ${exit_status} -eq 0 && ${nxf_unstage_ret:=0} -ne 0 ]] && exit_status=${nxf_unstage_ret:=0}
        [[ ${exit_status} -eq 0 && ${last_err} -ne 0 ]] && exit_status=${last_err}
        printf -- $exit_status > /workspaces/training/hello-nextflow/work/98/656c6c90cce1667c094d880f4b6dcc/.exitcode
        set +u
        docker rm $NXF_BOXID &>/dev/null || true
        exit $exit_status
    }

    on_term() {
        set +e
        docker stop $NXF_BOXID
    }

    nxf_launch() {
        docker run -i --cpu-shares 1024 -e "NXF_TASK_WORKDIR" -v /workspaces/training/hello-nextflow/work:/workspaces/training/hello-nextflow/work -w "$NXF_TASK_WORKDIR" --name $NXF_BOXID community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273 /bin/bash -ue /workspaces/training/hello-nextflow/work/98/656c6c90cce1667c094d880f4b6dcc/.command.sh
    }

    nxf_stage() {
        true
        # stage input files
        rm -f COLLECTED-batch-output.txt
        ln -s /workspaces/training/hello-nextflow/work/7f/f435e3f2cf95979b5f3d7647ae6696/COLLECTED-batch-output.txt COLLECTED-batch-output.txt
    }

    nxf_unstage_outputs() {
        true
    }

    nxf_unstage_controls() {
        true
    }

    nxf_unstage() {
        if [[ ${nxf_main_ret:=0} == 0 ]]; then
            (set -e -o pipefail; (nxf_unstage_outputs | tee -a .command.out) 3>&1 1>&2 2>&3 | tee -a .command.err)
            nxf_unstage_ret=$?
        fi
        nxf_unstage_controls
    }

    nxf_main() {
        trap on_exit EXIT
        trap on_term TERM INT USR2
        trap '' USR1

        [[ "${NXF_CHDIR:-}" ]] && cd "$NXF_CHDIR"
        export NXF_BOXID="nxf-$(dd bs=18 count=1 if=/dev/urandom 2>/dev/null | base64 | tr +/ 0A | tr -d '\r\n')"
        NXF_SCRATCH=''
        [[ $NXF_DEBUG > 0 ]] && nxf_env
        touch /workspaces/training/hello-nextflow/work/98/656c6c90cce1667c094d880f4b6dcc/.command.begin
        set +u
        set -u
        [[ $NXF_SCRATCH ]] && cd $NXF_SCRATCH
        export NXF_TASK_WORKDIR="$PWD"
        nxf_stage

        set +e
        (set -o pipefail; (nxf_launch | tee .command.out) 3>&1 1>&2 2>&3 | tee .command.err) &
        pid=$!
        wait $pid || nxf_main_ret=$?
        nxf_unstage
    }

    $NXF_ENTRY

    ```

Wenn du in dieser Datei nach `nxf_launch` suchst, solltest du etwas wie dies sehen:

```console
nxf_launch() {
    docker run -i --cpu-shares 1024 -e "NXF_TASK_WORKDIR" -v /workspaces/training/hello-nextflow/work:/workspaces/training/hello-nextflow/work -w "$NXF_TASK_WORKDIR" --name $NXF_BOXID community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273 /bin/bash -ue /workspaces/training/hello-nextflow/work/98/656c6c90cce1667c094d880f4b6dcc/.command.sh
}
```

Wie du sehen kannst, verwendet Nextflow den `docker run`-Befehl, um den Prozessaufruf zu starten.
Es mountet auch das entsprechende Arbeitsunterverzeichnis in den Container, setzt das Arbeitsverzeichnis innerhalb des Containers entsprechend und führt unser vorlagenbasiertes Bash-Skript in der `.command.sh`-Datei aus.

All die harte Arbeit, die wir im ersten Abschnitt manuell erledigen mussten? Nextflow erledigt sie für uns hinter den Kulissen!

```txt
 _______________________
< Hurray for robots...! >
 -----------------------
                                   ,-----.
                                   |     |
                                ,--|     |-.
                         __,----|  |     | |
                       ,;::     |  `_____' |
                       `._______|    i^i   |
                                `----| |---'| .
                           ,-------._| |== ||//
                           |       |_|P`.  /'/
                           `-------' 'Y Y/'/'
                                     .==\ /_\
   ^__^                             /   /'|  `i
   (oo)\_______                   /'   /  |   |
   (__)\       )\/\             /'    /   |   `i
       ||----w |           ___,;`----'.___L_,-'`\__
       ||     ||          i_____;----\.____i""\____\
```

### Fazit

Du weißt, wie man Container in Nextflow verwendet, um Prozesse auszuführen.

### Wie geht es weiter?

Mach eine Pause!

Wenn du bereit bist, gehe weiter zu [**Teil 6: Hello Config**](./06_hello_config.md), um zu lernen, wie du die Ausführung deiner Pipeline an deine Infrastruktur anpasst und die Konfiguration von Eingaben und Parametern verwaltest.

Es ist der allerletzte Teil, und dann bist du mit diesem Kurs fertig!

---

## Quiz

<quiz>
Was ist ein Container?
- [ ] Eine Art virtuelle Maschine
- [ ] Ein Dateikomprimierungsformat
- [x] Eine leichtgewichtige, eigenständige ausführbare Einheit, die alles enthält, was zum Ausführen einer Anwendung benötigt wird
- [ ] Ein Netzwerkprotokoll
</quiz>

<quiz>
Was ist der Unterschied zwischen einem Container-Image und einer Container-Instanz?
- [ ] Sie sind dasselbe
- [x] Ein Image ist eine Vorlage; eine Instanz ist ein laufender Container, der aus diesem Image erstellt wurde
- [ ] Eine Instanz ist eine Vorlage; ein Image ist ein laufender Container
- [ ] Images sind für Docker; Instanzen sind für Singularity
</quiz>

<quiz>
Was bewirkt das `-v`-Flag in einem `docker run`-Befehl?
- [ ] Aktiviert ausführliche Ausgabe
- [ ] Validiert den Container
- [x] Mountet ein Volume vom Host-System in den Container
- [ ] Gibt die Version des Containers an

Mehr erfahren: [1.3.4. Mounte Daten in den Container](#134-mount-data-into-the-container)
</quiz>

<quiz>
Warum musst du Volumes mounten, wenn du Container verwendest?
- [ ] Um die Container-Performance zu verbessern
- [ ] Um Speicherplatz zu sparen
- [x] Weil Container standardmäßig vom Host-Dateisystem isoliert sind
- [ ] Um Netzwerk zu aktivieren

Mehr erfahren: [1.3.4. Mounte Daten in den Container](#134-mount-data-into-the-container)
</quiz>

<quiz>
Wie gibst du einen Container für einen Nextflow-Prozess an?
- [ ] `docker 'container-uri'`
- [ ] `image 'container-uri'`
- [x] `container 'container-uri'`
- [ ] `use 'container-uri'`

Mehr erfahren: [2.3.1. Gib einen Container für cowpy an](#231-specify-a-container-for-cowpy)
</quiz>

<quiz>
Welche `nextflow.config`-Einstellung aktiviert Docker für deinen Workflow?
- [ ] `#!groovy process.docker = true`
- [x] `#!groovy docker.enabled = true`
- [ ] `#!groovy container.engine = 'docker'`
- [ ] `#!groovy docker.activate = true`

Mehr erfahren: [2.3.2. Aktiviere die Verwendung von Docker über die `nextflow.config`-Datei](#232-enable-use-of-docker-via-the-nextflowconfig-file)
</quiz>

<quiz>
Was übernimmt Nextflow automatisch, wenn ein Prozess in einem Container ausgeführt wird? (Wähle alle zutreffenden aus)
- [x] Pullen des Container-Images bei Bedarf
- [x] Mounten des Arbeitsverzeichnisses
- [x] Ausführen des Prozess-Skripts innerhalb des Containers
- [x] Aufräumen der Container-Instanz nach der Ausführung

Mehr erfahren: [2.3.4. Untersuche, wie Nextflow die containerisierte Aufgabe gestartet hat](#234-inspect-how-nextflow-launched-the-containerized-task)
</quiz>
