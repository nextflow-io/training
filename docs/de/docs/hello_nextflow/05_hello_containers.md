# Teil 5: Hello Containers

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } KI-gestützte Übersetzung - [mehr erfahren & Verbesserungen vorschlagen](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<!--
<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/5PyOWjKnNmg?si=QinuAnFwFj-Z8CrO&amp;list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

/// caption
:fontawesome-brands-youtube:{ .youtube } Sieh dir [die ganze Playlist](https://www.youtube.com/playlist?list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik) auf dem Nextflow YouTube-Kanal an.

:green_book: Das Videotranskript ist [hier](./transcripts/05_hello_containers.md) verfügbar.
///
-->

In den Teilen 1-4 dieses Trainings hast du gelernt, wie du die grundlegenden Bausteine von Nextflow verwendest, um einen einfachen Workflow zu erstellen, der Text verarbeiten, die Ausführung bei mehreren Eingaben parallelisieren und die Ergebnisse zur weiteren Verarbeitung sammeln kann.

Du warst jedoch auf grundlegende UNIX-Tools beschränkt, die in deiner Umgebung verfügbar sind.
Reale Aufgaben erfordern oft verschiedene Tools und Pakete, die standardmäßig nicht enthalten sind.
Normalerweise müsstest du diese Tools installieren, ihre Abhängigkeiten verwalten und eventuelle Konflikte lösen.

Das ist alles sehr mühsam und nervig, daher zeigen wir dir, wie du **Container** verwendest, um dieses Problem viel bequemer zu lösen.

Ein **Container** ist eine leichtgewichtige, eigenständige, ausführbare Softwareeinheit, die aus einem Container-**Image** erstellt wird und alles enthält, was zum Ausführen einer Anwendung benötigt wird, einschließlich Code, Systembibliotheken und Einstellungen.
Wie du dir vorstellen kannst, wird das sehr hilfreich sein, um deine Pipelines reproduzierbarer zu machen.

Beachte, dass wir dies mit [Docker](https://www.docker.com/get-started/) lehren, aber Nextflow unterstützt auch [mehrere andere Container-Technologien](https://www.nextflow.io/docs/latest/container.html#).

??? info "Wie du von diesem Abschnitt aus beginnst"

    Dieser Abschnitt des Kurses setzt voraus, dass du die Teile 1-4 des [Hello Nextflow](./index.md)-Kurses abgeschlossen hast und eine vollständig funktionierende Pipeline hast.

    Wenn du den Kurs von diesem Punkt aus beginnst, musst du das `modules`-Verzeichnis aus den Lösungen kopieren:

    ```bash
    cp -r solutions/4-hello-modules/modules .
    ```

---

## 0. Aufwärmen: `hello-containers.nf` ausführen

Wir werden das Workflow-Script `hello-containers.nf` als Ausgangspunkt verwenden.
Es entspricht dem Script, das durch Durcharbeiten von Teil 4 dieses Trainings entstanden ist, außer dass wir die Ausgabeziele geändert haben:

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

Um sicherzustellen, dass alles funktioniert, führe das Script einmal aus, bevor du Änderungen vornimmst:

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

Wie zuvor findest du die Ausgabedateien in dem im `output`-Block angegebenen Verzeichnis (`results/hello_containers/`).

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

## 1. Einen Container 'manuell' verwenden

Was wir tun möchten, ist einen Schritt zu unserem Workflow hinzuzufügen, der einen Container zur Ausführung verwendet.

Allerdings werden wir zuerst einige grundlegende Konzepte und Operationen durchgehen, um dein Verständnis davon zu festigen, was Container sind, bevor wir sie in Nextflow verwenden.

### 1.1. Das Container-Image herunterladen

Um einen Container zu verwenden, lädst du normalerweise ein Container-Image von einer Container-Registry herunter oder _pullst_ es, und führst dann das Container-Image aus, um eine Container-Instanz zu erstellen.

Die allgemeine Syntax ist wie folgt:

```bash title="Syntax"
docker pull '<container>'
```

Der `docker pull`-Teil ist die Anweisung an das Container-System, ein Container-Image aus einem Repository zu pullen.

Der `'<container>'`-Teil ist die URI-Adresse des Container-Images.

Als Beispiel pullen wir ein Container-Image, das [cowpy](https://github.com/jeffbuttars/cowpy) enthält, eine Python-Implementierung eines Tools namens `cowsay`, das ASCII-Kunst generiert, um beliebige Texteingaben auf lustige Weise anzuzeigen.

```txt title="Beispiel"
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
Wir haben den [Seqera Containers](https://seqera.io/containers/)-Service verwendet, um dieses Docker-Container-Image aus dem `cowpy`-Conda-Paket zu generieren: `'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'`.

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

### 1.2. Den Container verwenden, um `cowpy` als einmaligen Befehl auszuführen

Eine sehr häufige Art, wie Leute Container verwenden, ist sie direkt auszuführen, _d.h._ nicht-interaktiv.
Das ist großartig für einmalige Befehle.

Die allgemeine Syntax ist wie folgt:

```bash title="Syntax"
docker run --rm '<container>' [tool command]
```

Der `docker run --rm '<container>'`-Teil ist die Anweisung an das Container-System, eine Container-Instanz aus einem Container-Image zu starten und einen Befehl darin auszuführen.
Das `--rm`-Flag teilt dem System mit, die Container-Instanz nach Abschluss des Befehls herunterzufahren.

Die `[tool command]`-Syntax hängt vom verwendeten Tool und der Einrichtung des Containers ab.
Beginnen wir einfach mit `cowpy`.

Vollständig zusammengestellt sieht der Container-Ausführungsbefehl so aus; führe ihn aus.

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

### 1.3. Den Container verwenden, um `cowpy` interaktiv auszuführen

Du kannst einen Container auch interaktiv ausführen, was dir eine Shell-Eingabeaufforderung innerhalb des Containers gibt und es dir ermöglicht, mit dem Befehl zu experimentieren.

#### 1.3.1. Den Container starten

Um interaktiv auszuführen, fügen wir einfach `-it` zum `docker run`-Befehl hinzu.
Optional können wir die Shell angeben, die wir innerhalb des Containers verwenden möchten, indem wir z.B. `/bin/bash` an den Befehl anhängen.

```bash
docker run --rm -it 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273' /bin/bash
```

Beachte, dass sich deine Eingabeaufforderung zu etwas wie `(base) root@b645838b3314:/tmp#` ändert, was anzeigt, dass du dich jetzt innerhalb des Containers befindest.

Du kannst dies überprüfen, indem du `ls /` ausführst, um den Verzeichnisinhalt vom Wurzelverzeichnis des Dateisystems aufzulisten:

```bash
ls /
```

??? abstract "Befehlsausgabe"

    ```console
    bin  boot  dev  etc  home  lib  lib64  media  mnt  opt  proc  root  run  sbin  srv  sys  tmp  usr  var
    ```

Wir verwenden hier `ls` anstelle von `tree`, weil das `tree`-Tool in diesem Container nicht verfügbar ist.
Du kannst sehen, dass das Dateisystem innerhalb des Containers sich vom Dateisystem auf deinem Host-System unterscheidet.

Eine Einschränkung dessen, was wir gerade getan haben, ist, dass der Container standardmäßig vollständig vom Host-System isoliert ist.
Das bedeutet, dass der Container auf keine Dateien auf dem Host-System zugreifen kann, es sei denn, du erlaubst es ausdrücklich.

Wir zeigen dir gleich, wie das geht.

#### 1.3.2. Die gewünschten Tool-Befehle ausführen

Jetzt, da du dich innerhalb des Containers befindest, kannst du den `cowpy`-Befehl direkt ausführen und ihm einige Parameter geben.
Zum Beispiel sagt die Tool-Dokumentation, dass wir den Charakter ('Cowacter') mit `-c` ändern können.

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

Da du dich innerhalb des Containers befindest, kannst du den `cowpy`-Befehl beliebig oft ausführen und die Eingabeparameter variieren, ohne dich mit Docker-Befehlen befassen zu müssen.

!!! tip "Tipp"

    Verwende das '-c'-Flag, um einen anderen Charakter auszuwählen, einschließlich:
    `beavis`, `cheese`, `daemon`, `dragonandcow`, `ghostbusters`, `kitty`, `moose`, `milk`, `stegosaurus`, `turkey`, `turtle`, `tux`

Das ist cool. Noch cooler wäre es, wenn wir unsere `greetings.csv` als Eingabe verwenden könnten.
Aber da wir keinen Zugriff auf das Dateisystem haben, können wir das nicht.

Lass uns das beheben.

#### 1.3.3. Den Container verlassen

Um den Container zu verlassen, kannst du `exit` an der Eingabeaufforderung eingeben oder die Tastenkombination ++ctrl+d++ verwenden.

```bash
exit
```

Deine Eingabeaufforderung sollte jetzt wieder so aussehen wie vor dem Starten des Containers.

#### 1.3.4. Daten in den Container einbinden

Wie bereits erwähnt, ist der Container standardmäßig vom Host-System isoliert.

Um dem Container den Zugriff auf das Host-Dateisystem zu ermöglichen, kannst du ein **Volume** vom Host-System in den Container **einbinden** (mounten), indem du die folgende Syntax verwendest:

```bash title="Syntax"
-v <outside_path>:<inside_path>
```

In unserem Fall wird `<outside_path>` das aktuelle Arbeitsverzeichnis sein, also können wir einfach einen Punkt (`.`) verwenden, und `<inside_path>` ist nur ein Alias, den wir erfinden; nennen wir es `/my_project` (der innere Pfad muss absolut sein).

Um ein Volume einzubinden, ersetzen wir die Pfade und fügen das Volume-Mounting-Argument zum Docker-Run-Befehl wie folgt hinzu:

```bash
docker run --rm -it -v .:/my_project 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273' /bin/bash
```

Dies bindet das aktuelle Arbeitsverzeichnis als Volume ein, das unter `/my_project` innerhalb des Containers zugänglich sein wird.

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

Dies hat effektiv einen Tunnel durch die Containerwand etabliert, den du verwenden kannst, um auf diesen Teil deines Dateisystems zuzugreifen.

#### 1.3.5. Die eingebundenen Daten verwenden

Jetzt, da wir das Arbeitsverzeichnis in den Container eingebunden haben, können wir den `cowpy`-Befehl verwenden, um den Inhalt der `greetings.csv`-Datei anzuzeigen.

Dazu verwenden wir `cat /my_project/data/greetings.csv | `, um den Inhalt der CSV-Datei an den `cowpy`-Befehl zu übergeben.

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
              (      /  (_))^^)) )  )  ))^^^^^))^^^)__/     +^^
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

Dies erzeugt die gewünschte ASCII-Kunst eines Truthahns, der unsere Beispiel-Grüße aufsagt!
Nur dass hier der Truthahn die vollständigen Zeilen wiederholt anstatt nur die Grüße.
Wir wissen bereits, dass unser Nextflow-Workflow das besser machen wird!

Experimentiere gerne mit diesem Befehl.
Wenn du fertig bist, verlasse den Container wie zuvor:

```bash
exit
```

Du befindest dich wieder in deiner normalen Shell.

### Fazit

Du weißt, wie man Container pullt und entweder einmalig oder interaktiv ausführt. Du weißt auch, wie du deine Daten im Container zugänglich machst, was dich beliebige Tools mit echten Daten ausprobieren lässt, ohne Software auf deinem System installieren zu müssen.

### Wie geht es weiter?

Lerne, wie man Container für die Ausführung von Nextflow-Prozessen verwendet.

---

## 2. Container in Nextflow verwenden

Nextflow hat eingebaute Unterstützung für die Ausführung von Prozessen innerhalb von Containern, damit du Tools ausführen kannst, die nicht in deiner Rechenumgebung installiert sind.
Das bedeutet, dass du jedes beliebige Container-Image verwenden kannst, um deine Prozesse auszuführen, und Nextflow kümmert sich um das Pullen des Images, das Einbinden der Daten und das Ausführen des Prozesses darin.

Um dies zu demonstrieren, werden wir einen `cowpy`-Schritt zur Pipeline hinzufügen, die wir entwickelt haben, nach dem `collectGreetings`-Schritt.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-pipeline-cowpy.svg"
</figure>

### 2.1. Ein `cowpy`-Modul schreiben

Zuerst erstellen wir das `cowpy`-Prozess-Modul.

#### 2.1.1. Eine Datei-Vorlage für das neue Modul erstellen

Erstelle eine leere Datei für das Modul namens `cowpy.nf`.

```bash
touch modules/cowpy.nf
```

Das gibt uns einen Platz für den Prozess-Code.

#### 2.1.2. Den `cowpy`-Prozess-Code in die Modul-Datei kopieren

Wir können unseren `cowpy`-Prozess nach den anderen Prozessen modellieren, die wir zuvor geschrieben haben.

```groovy title="modules/cowpy.nf" linenums="1"
#!/usr/bin/env nextflow

// ASCII-Kunst mit cowpy generieren
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

Der Prozess erwartet eine `input_file`, die die Grüße enthält, sowie einen `character`-Wert.

Die Ausgabe wird eine neue Textdatei sein, die die vom `cowpy`-Tool generierte ASCII-Kunst enthält.

### 2.2. cowpy zum Workflow hinzufügen

Jetzt müssen wir das Modul importieren und den Prozess aufrufen.

#### 2.2.1. Den `cowpy`-Prozess in `hello-containers.nf` importieren

Füge die Import-Deklaration über dem Workflow-Block ein und fülle sie entsprechend aus.

=== "Nachher"

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

Jetzt ist das `cowpy`-Modul für die Verwendung im Workflow verfügbar.

#### 2.2.2. Einen Aufruf des `cowpy`-Prozesses im Workflow hinzufügen

Verbinden wir den `cowpy()`-Prozess mit der Ausgabe des `collectGreetings()`-Prozesses, der, wie du dich vielleicht erinnerst, zwei Ausgaben produziert:

- `collectGreetings.out.outfile` enthält die Ausgabedatei <--_was wir wollen_
- `collectGreetings.out.report` enthält die Berichtsdatei mit der Anzahl der Grüße pro Batch

Nimm im Workflow-Block folgende Code-Änderung vor:

=== "Nachher"

    ```groovy title="hello-containers.nf" linenums="19" hl_lines="12-13"
        main:
        // Einen Channel für Eingaben aus einer CSV-Datei erstellen
        greeting_ch = channel.fromPath(params.input)
                            .splitCsv()
                            .map { line -> line[0] }
        // Eine Begrüßung ausgeben
        sayHello(greeting_ch)
        // Die Begrüßung in Großbuchstaben umwandeln
        convertToUpper(sayHello.out)
        // Alle Begrüßungen in einer Datei sammeln
        collectGreetings(convertToUpper.out.collect(), params.batch)
        // ASCII-Kunst der Begrüßungen mit cowpy generieren
        cowpy(collectGreetings.out.outfile, params.character)
    ```

=== "Vorher"

    ```groovy title="hello-containers.nf" linenums="19"
        main:
        // Einen Channel für Eingaben aus einer CSV-Datei erstellen
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

Beachte, dass wir einen neuen CLI-Parameter, `params.character`, deklariert haben, um anzugeben, welcher Charakter die Grüße sagen soll.

#### 2.2.3. Den `character`-Parameter zum `params`-Block hinzufügen

Dies ist technisch optional, aber es ist die empfohlene Praxis und eine Gelegenheit, einen Standardwert für den Charakter festzulegen, während wir dabei sind.

=== "Nachher"

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

Jetzt können wir faul sein und das Tippen des character-Parameters in unseren Befehlszeilen überspringen.

#### 2.2.4. Die Workflow-Ausgaben aktualisieren

Wir müssen die Workflow-Ausgaben aktualisieren, um die Ausgabe des `cowpy`-Prozesses zu veröffentlichen.

##### 2.2.4.1. Den `publish:`-Abschnitt aktualisieren

Nimm im `workflow`-Block folgende Code-Änderung vor:

=== "Nachher"

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

Der `cowpy`-Prozess produziert nur eine Ausgabe, also können wir uns darauf auf die übliche Weise beziehen, indem wir `.out` anhängen.

Aber für jetzt beenden wir die Aktualisierung der Workflow-Level-Ausgaben.

##### 2.2.4.2. Den `output`-Block aktualisieren

Wir müssen die finale `cowpy_art`-Ausgabe zum `output`-Block hinzufügen. Während wir dabei sind, bearbeiten wir auch die Veröffentlichungsziele, da unsere Pipeline jetzt vollständig ist und wir wissen, welche Ausgaben uns wirklich wichtig sind.

Nimm im `output`-Block folgende Code-Änderungen vor:

=== "Nachher"

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

Jetzt werden die veröffentlichten Ausgaben etwas besser organisiert sein.

#### 2.2.5. Den Workflow ausführen

Nur zur Zusammenfassung, das ist, was wir anstreben:

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello_pipeline_complete.svg"
</figure>

Denkst du, es wird funktionieren?

Lass uns die vorherigen veröffentlichten Ausgaben löschen, um einen sauberen Stand zu haben, und den Workflow mit dem `-resume`-Flag ausführen.

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
Der Fehlercode `error exit status (127)` bedeutet, dass die angeforderte ausführbare Datei nicht gefunden wurde.

Das ergibt Sinn, da wir das `cowpy`-Tool aufrufen, aber wir haben noch keinen Container angegeben (ups).

### 2.3. Einen Container verwenden, um den `cowpy`-Prozess auszuführen

Wir müssen einen Container angeben und Nextflow anweisen, ihn für den `cowpy()`-Prozess zu verwenden.

#### 2.3.1. Einen Container für `cowpy` angeben

Wir können dasselbe Image verwenden, das wir direkt im ersten Abschnitt dieses Tutorials verwendet haben.

Bearbeite das `cowpy.nf`-Modul, um die `container`-Direktive zur Prozessdefinition wie folgt hinzuzufügen:

=== "Nachher"

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

#### 2.3.2. Die Verwendung von Docker über die `nextflow.config`-Datei aktivieren

Beachte, dass wir sagten _'wenn die Verwendung von Docker aktiviert ist'_. Standardmäßig ist sie es nicht, also müssen wir Nextflow mitteilen, dass es Docker verwenden darf.
Dazu werden wir das Thema des nächsten und letzten Teils dieses Kurses (Teil 6) leicht vorwegnehmen, der die Konfiguration behandelt.

Eine der Hauptmöglichkeiten, die Nextflow für die Konfiguration der Workflow-Ausführung bietet, ist die Verwendung einer `nextflow.config`-Datei.
Wenn eine solche Datei im aktuellen Verzeichnis vorhanden ist, lädt Nextflow sie automatisch und wendet die enthaltene Konfiguration an.

Wir haben eine `nextflow.config`-Datei mit einer einzigen Codezeile bereitgestellt, die Docker explizit deaktiviert: `docker.enabled = false`.

Jetzt ändern wir das auf `true`, um Docker zu aktivieren:

=== "Nachher"

    ```console title="nextflow.config" linenums="1" hl_lines="1"
    docker.enabled = true
    ```

=== "Vorher"

    ```console title="nextflow.config" linenums="1" hl_lines="1"
    docker.enabled = false
    ```

!!! tip "Tipp"

    Es ist möglich, die Docker-Ausführung über die Befehlszeile pro Ausführung mit dem Parameter `-with-docker <container>` zu aktivieren.
    Das erlaubt uns jedoch nur, einen Container für den gesamten Workflow anzugeben, während der gerade gezeigte Ansatz es uns ermöglicht, einen anderen Container pro Prozess anzugeben.
    Das ist besser für Modularität, Code-Wartung und Reproduzierbarkeit.

#### 2.3.3. Den Workflow mit aktiviertem Docker ausführen

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
Wie üblich findest du die Workflow-Ausgaben im entsprechenden Ergebnisverzeichnis, obwohl sie diesmal etwas übersichtlicher organisiert sind, mit nur dem Bericht und der finalen Ausgabe auf der obersten Ebene, und alle Zwischendateien in einem Unterverzeichnis verstaut.

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

Die finale ASCII-Kunst-Ausgabe befindet sich im Verzeichnis `results/hello_containers/`, unter dem Namen `cowpy-COLLECTED-batch-output.txt`.

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
              (      /  (_))^^)) )  )  ))^^^^^))^^^)__/     +^^
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

Und da ist er, unser wunderschöner Truthahn, der die Grüße wie gewünscht sagt.

#### 2.3.4. Untersuchen, wie Nextflow die containerisierte Aufgabe gestartet hat

Als letzten Abschluss dieses Abschnitts werfen wir einen Blick auf das Work-Unterverzeichnis für einen der `cowpy`-Prozessaufrufe, um etwas mehr Einblick zu bekommen, wie Nextflow unter der Haube mit Containern arbeitet.

Überprüfe die Ausgabe deines `nextflow run`-Befehls, um den Pfad zum Work-Unterverzeichnis für den `cowpy`-Prozess zu finden.
Wenn wir uns ansehen, was wir für die oben gezeigte Ausführung erhalten haben, beginnt die Konsolenprotokollzeile für den `cowpy`-Prozess mit `[98/656c6c]`.
Das entspricht dem folgenden abgekürzten Verzeichnispfad: `work/98/656c6c`.

In diesem Verzeichnis findest du die `.command.run`-Datei, die alle Befehle enthält, die Nextflow in deinem Namen während der Ausführung der Pipeline ausgeführt hat.

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

    ...

    nxf_launch() {
        docker run -i --cpu-shares 1024 -e "NXF_TASK_WORKDIR" -v /workspaces/training/hello-nextflow/work:/workspaces/training/hello-nextflow/work -w "$NXF_TASK_WORKDIR" --name $NXF_BOXID community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273 /bin/bash -ue /workspaces/training/hello-nextflow/work/98/656c6c90cce1667c094d880f4b6dcc/.command.sh
    }

    ...
    ```

Wenn du in dieser Datei nach `nxf_launch` suchst, solltest du so etwas sehen:

```console
nxf_launch() {
    docker run -i --cpu-shares 1024 -e "NXF_TASK_WORKDIR" -v /workspaces/training/hello-nextflow/work:/workspaces/training/hello-nextflow/work -w "$NXF_TASK_WORKDIR" --name $NXF_BOXID community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273 /bin/bash -ue /workspaces/training/hello-nextflow/work/98/656c6c90cce1667c094d880f4b6dcc/.command.sh
}
```

Wie du sehen kannst, verwendet Nextflow den `docker run`-Befehl, um den Prozessaufruf zu starten.
Es bindet auch das entsprechende Work-Unterverzeichnis in den Container ein, setzt das Arbeitsverzeichnis innerhalb des Containers entsprechend und führt unser templated Bash-Script in der `.command.sh`-Datei aus.

All die harte Arbeit, die wir im ersten Abschnitt manuell erledigen mussten? Nextflow erledigt das für uns hinter den Kulissen!

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
   ^__^                         /   /'|  `i
   (oo)\_______               /'   /  |   |
   (__)\       )\/\         /'    /   |   `i
       ||----w |           ___,;`----'.___L_,-'`\__
       ||     ||          i_____;----\.____i""\____\
```

### Fazit

Du weißt, wie man Container in Nextflow verwendet, um Prozesse auszuführen.

### Wie geht es weiter?

Mach eine Pause!

Wenn du bereit bist, gehe zu [**Teil 6: Hello Config**](./06_hello_config.md) weiter, um zu lernen, wie du die Ausführung deiner Pipeline an deine Infrastruktur anpasst sowie die Konfiguration von Eingaben und Parametern verwaltest.

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
Was macht das `-v`-Flag in einem `docker run`-Befehl?
- [ ] Aktiviert ausführliche Ausgabe
- [ ] Validiert den Container
- [x] Bindet ein Volume vom Host-System in den Container ein
- [ ] Gibt die Version des Containers an

Mehr erfahren: [1.3.4. Daten in den Container einbinden](#134-daten-in-den-container-einbinden)
</quiz>

<quiz>
Warum musst du Volumes einbinden, wenn du Container verwendest?
- [ ] Um die Container-Leistung zu verbessern
- [ ] Um Speicherplatz zu sparen
- [x] Weil Container standardmäßig vom Host-Dateisystem isoliert sind
- [ ] Um Netzwerke zu aktivieren

Mehr erfahren: [1.3.4. Daten in den Container einbinden](#134-daten-in-den-container-einbinden)
</quiz>

<quiz>
Wie gibst du einen Container für einen Nextflow-Prozess an?
- [ ] `docker 'container-uri'`
- [ ] `image 'container-uri'`
- [x] `container 'container-uri'`
- [ ] `use 'container-uri'`

Mehr erfahren: [2.3.1. Einen Container für cowpy angeben](#231-einen-container-für-cowpy-angeben)
</quiz>

<quiz>
Welche `nextflow.config`-Einstellung aktiviert Docker für deinen Workflow?
- [ ] `#!groovy process.docker = true`
- [x] `#!groovy docker.enabled = true`
- [ ] `#!groovy container.engine = 'docker'`
- [ ] `#!groovy docker.activate = true`

Mehr erfahren: [2.3.2. Die Verwendung von Docker über die `nextflow.config`-Datei aktivieren](#232-die-verwendung-von-docker-über-die-nextflowconfig-datei-aktivieren)
</quiz>

<quiz>
Was handhabt Nextflow automatisch, wenn ein Prozess in einem Container ausgeführt wird? (Wähle alle zutreffenden aus)
- [x] Das Container-Image bei Bedarf pullen
- [x] Das Work-Verzeichnis einbinden
- [x] Das Prozess-Script innerhalb des Containers ausführen
- [x] Die Container-Instanz nach der Ausführung bereinigen

Mehr erfahren: [2.3.4. Untersuchen, wie Nextflow die containerisierte Aufgabe gestartet hat](#234-untersuchen-wie-nextflow-die-containerisierte-aufgabe-gestartet-hat)
</quiz>
