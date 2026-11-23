# Teil 3: Workflow-Ausführungen verwalten

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } KI-gestützte Übersetzung - [mehr erfahren & Verbesserungen vorschlagen](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Wenn du Pipelines ausführst und erneut ausführst, sammeln sich Ausführungshistorie und alte `work/`-Verzeichnisse an.
In [Teil 1](./01_run_nextflow.md#23-use-resume-to-skip-completed-work) hast du bereits `-resume` verwendet, um bereits erledigte Arbeit zu überspringen.
Hier lernst du, wie du Berichte über eine Ausführung erstellst, die Historie vergangener Ausführungen mit [`nextflow log`](https://nextflow.io/docs/latest/reference/cli.html#log) einsehen und alte, nicht mehr benötigte work-Verzeichnisse mit [`nextflow clean`](https://nextflow.io/docs/latest/reference/cli.html#clean) löschen kannst.

---

## 1. Pipeline-Berichte erstellen

Nextflow kann verschiedene Arten von Berichten über eine Ausführung erstellen. Jeder wird mit einem eigenen `-with-*`-Flag hinzugefügt: ein Ausführungsbericht (`-with-report`), eine Ausführungs-Timeline (`-with-timeline`), eine Task-Trace-Datei (`-with-trace`) und ein Workflow-Diagramm (`-with-dag`).
Wir erstellen hier die ersten beiden; die übrigen findest du unter [Execution reports](https://nextflow.io/docs/latest/reports.html) in der Nextflow-Referenz.

### 1.1. Einen Ausführungsbericht erstellen

Füge `-with-report` zu einem beliebigen `nextflow run`-Befehl hinzu, um nach Abschluss der Pipeline einen HTML-Bericht zu erstellen:

```bash
nextflow run main.nf --input data/greetings.csv --character turkey -with-report
```

??? success "Befehlsausgabe"

    ```console hl_lines="6 7 8 9"
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [intergalactic_dalembert] revision: ce74f81996

    executor >  local (8)
    [34/23f10f] sayHello (2)       | 3 of 3 ✔
    [af/cab69d] convertToUpper (3) | 3 of 3 ✔
    [9e/d73afb] collectGreetings   | 1 of 1 ✔
    [3c/392db0] cowpy              | 1 of 1 ✔

    Outputs:

      /workspaces/training/nextflow-run/results

      first_output:
        - batch/intermediates/Bonjour-output.txt
        - batch/intermediates/Hello-output.txt
        - batch/intermediates/Hola-output.txt

      uppercased:
        - batch/intermediates/UPPER-Bonjour-output.txt
        - batch/intermediates/UPPER-Hello-output.txt
        - batch/intermediates/UPPER-Hola-output.txt

      collected: batch/intermediates/COLLECTED-batch-output.txt

      batch_report: batch/batch-report.txt

      cowpy_art: batch/cowpy-COLLECTED-batch-output.txt
    ```

Nextflow schreibt den Bericht in eine Datei namens `report-<timestamp>.html` im Arbeitsverzeichnis.
Öffne sie in einem Browser, um eine Ausführungszusammenfassung, eine Tabelle aller Tasks mit Status und Laufzeit sowie Diagramme zur Ressourcennutzung aufgeschlüsselt nach Prozess zu sehen.

Der Tab **Tasks** listet alle Tasks auf, die die Pipeline ausgeführt hat, mit Prozessname, Status und Ressourcennutzung:

![Aufgabentabelle im Ausführungsbericht](img/execution_report_tasks.png)

Der Bericht ist besonders nützlich, wenn eine Pipeline länger als erwartet dauert oder ein Task fehlschlägt: Die Task-Tabelle zeigt genau, wo Zeit verbraucht wurde und welche Tasks erfolgreich waren oder fehlgeschlagen sind.

### 1.2. Eine Ausführungs-Timeline erstellen

Füge `-with-timeline` zu einer Ausführung hinzu, um eine Gantt-Diagramm-ähnliche Ansicht zu erhalten, die zeigt, wann jeder Task ausgeführt wurde:

```bash
nextflow run main.nf --input data/greetings.csv --character turkey -with-timeline
```

??? success "Befehlsausgabe"

    ```console hl_lines="6 7 8 9"
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [jolly_noyce] revision: ce74f81996

    executor >  local (8)
    [ad/e92ef3] sayHello (3)       | 3 of 3 ✔
    [2a/df8a8d] convertToUpper (2) | 3 of 3 ✔
    [be/7fb72a] collectGreetings   | 1 of 1 ✔
    [63/dc9bd6] cowpy              | 1 of 1 ✔

    Outputs:

      /workspaces/training/nextflow-run/results

      first_output:
        - batch/intermediates/Bonjour-output.txt
        - batch/intermediates/Hello-output.txt
        - batch/intermediates/Hola-output.txt

      uppercased:
        - batch/intermediates/UPPER-Bonjour-output.txt
        - batch/intermediates/UPPER-Hello-output.txt
        - batch/intermediates/UPPER-Hola-output.txt

      collected: batch/intermediates/COLLECTED-batch-output.txt

      batch_report: batch/batch-report.txt

      cowpy_art: batch/cowpy-COLLECTED-batch-output.txt
    ```

Nextflow schreibt die Timeline in eine Datei namens `timeline-<timestamp>.html`.
Öffne sie in einem Browser, um für jeden Task einen Balken zu sehen, der nach Startzeitpunkt und Dauer positioniert und skaliert ist:

![Ausführungs-Timeline](img/execution_timeline.png)

Die Timeline macht die Fächer-auf-dann-Fächer-zu-Form aus [Teil 1](./01_run_nextflow.md#31-run-the-workflow) auf einen Blick sichtbar: Die drei `sayHello`-Tasks laufen parallel, dann die drei `convertToUpper`-Tasks, dann laufen `collectGreetings` und `cowpy` nacheinander, da jeder von allem vorherigen abhängt.

### Fazit

Du weißt, wie du mit `-with-report` einen HTML-Ausführungsbericht und mit `-with-timeline` eine Ausführungs-Timeline erstellst, und wo du die anderen von Nextflow unterstützten Berichtstypen findest.

### Wie geht es weiter?

Lerne, wie du die Historie vergangener Ausführungen einsehen kannst.

---

## 2. Das Log vergangener Ausführungen einsehen

Egal ob du eine Pipeline entwickelst oder produktiv betreibst – irgendwann musst du Informationen über vergangene Ausführungen nachschlagen.

### 2.1. Die Verlaufsdatei

Jedes Mal, wenn du einen Nextflow-Workflow startest, wird eine Zeile in eine Log-Datei namens `history` geschrieben, die sich in einem versteckten Verzeichnis namens `.nextflow` im aktuellen Arbeitsverzeichnis befindet.

??? abstract "Dateiinhalt"

    ```txt title=".nextflow/history" linenums="1"
    2026-09-13 01:47:58	2.8s	determined_ramanujan	OK	ce74f81996b8ce999853fdbb25ba7969	03078950-5011-414f-992b-8fe1bd219ac2	nextflow run main.nf --input data/greetings.csv --character turkey
    2026-09-13 01:48:03	3s	prickly_cuvier	OK	ce74f81996b8ce999853fdbb25ba7969	5a9e8bf7-a7db-4c2c-b1d8-3bd4c7266528	nextflow run main.nf --input data/greetings.csv --character tux
    2026-09-13 01:48:09	2.3s	elegant_panini	OK	ce74f81996b8ce999853fdbb25ba7969	517ccc29-9db7-4d8d-8ea0-c8db6e9ac653	nextflow run main.nf --input data/greetings.csv --character stegosaurus
    2026-09-13 01:48:14	1.5s	deadly_lamport	OK	ce74f81996b8ce999853fdbb25ba7969	517ccc29-9db7-4d8d-8ea0-c8db6e9ac653	nextflow run main.nf --input data/greetings.csv --character stegosaurus -resume
    ```

Jede Zeile enthält Zeitstempel, Dauer, Run-Name, Status, Revisions-ID, Session-ID und die vollständige Befehlszeile einer Ausführung, die aus diesem Verzeichnis gestartet wurde.

Schau dir die letzten beiden Zeilen an: Es sind zwei separate Aufrufe (einer normal, einer mit `-resume`) desselben Befehls, und sie teilen dieselbe Session-ID.
Die Session-ID ändert sich nur, wenn du eine wirklich neue Ausführung startest; bei Verwendung von `-resume` bleibt sie erhalten – so weiß Nextflow, welchen Cache es wiederverwenden soll.

### 2.2. `nextflow log` für eine übersichtlichere Ansicht verwenden

Die rohe Verlaufsdatei zu lesen funktioniert, aber `nextflow log` formatiert dieselben Informationen mit einer Kopfzeile:

```bash
nextflow log
```

??? success "Befehlsausgabe"

    ```console linenums="1"
    TIMESTAMP          	DURATION	RUN NAME            	STATUS	REVISION ID	SESSION ID                          	COMMAND
    2026-09-13 01:47:58	2.8s    	determined_ramanujan	OK    	ce74f81996 	03078950-5011-414f-992b-8fe1bd219ac2	nextflow run main.nf --input data/greetings.csv --character turkey
    2026-09-13 01:48:03	3s      	prickly_cuvier      	OK    	ce74f81996 	5a9e8bf7-a7db-4c2c-b1d8-3bd4c7266528	nextflow run main.nf --input data/greetings.csv --character tux
    2026-09-13 01:48:09	2.3s    	elegant_panini      	OK    	ce74f81996 	517ccc29-9db7-4d8d-8ea0-c8db6e9ac653	nextflow run main.nf --input data/greetings.csv --character stegosaurus
    2026-09-13 01:48:14	1.5s    	deadly_lamport      	OK    	ce74f81996 	517ccc29-9db7-4d8d-8ea0-c8db6e9ac653	nextflow run main.nf --input data/greetings.csv --character stegosaurus -resume
    ```

Nextflow gruppiert die Caching-Informationen, die es für `-resume` verwendet, unter `.nextflow/cache`, geordnet nach Session-ID.
Deshalb ist das Nachschlagen des richtigen Run-Namens oder der Session-ID hier der erste Schritt, wenn du eine vergangene Ausführung untersuchen oder bereinigen möchtest.

### Fazit

Du weißt, wo Nextflow die Historie vergangener Ausführungen aufzeichnet, und wie du sie mit `nextflow log` einsehen kannst.

### Wie geht es weiter?

Lerne, wie du alte, nicht mehr benötigte work-Verzeichnisse entfernen kannst.

---

## 3. Ältere work-Verzeichnisse löschen

Jede Ausführung hinterlässt ihre Task-Verzeichnisse unter `work/`, auch nachdem du die gewünschten Ausgaben nach `results/` kopiert hast.
Führst du während der Entwicklung genug Pipelines aus, summieren sich diese Unterverzeichnisse – deshalb bietet Nextflow `nextflow clean`, um nicht mehr benötigte zu entfernen.

### 3.1. Löschkriterien festlegen

`nextflow clean` unterstützt verschiedene Möglichkeiten, auszuwählen, was entfernt werden soll; die vollständige Liste findest du in der [Referenzdokumentation](https://www.nextflow.io/docs/latest/reference/cli.html#clean).
Hier löschst du alles von Ausführungen vor einer bestimmten Ausführung anhand ihres Run-Namens.

Suche mit `nextflow log` die aktuellste Ausführung, die du behalten möchtest; im [Beispiel aus 2.2](#22-use-nextflow-log-for-a-friendlier-view) ist das `elegant_panini`, die letzte normale Ausführung vor der `-resume`-Ausführung.
Der Run-Name ist der maschinell generierte zweiteilige String, der in der `Launching (...)`-Konsolenzeile oder in der Spalte `RUN NAME` von `nextflow log` angezeigt wird.

### 3.2. Einen Probelauf durchführen

Füge zuerst `-n` hinzu, um zu prüfen, was ein bestimmter Befehl löschen würde, ohne tatsächlich etwas zu löschen:

```bash
nextflow clean -before elegant_panini -n
```

??? success "Befehlsausgabe"

    ```console
    Would remove /workspaces/training/nextflow-run/work/e5/9ec3fe24deb5d1ffa478dff5e5e1d4
    Would remove /workspaces/training/nextflow-run/work/75/36b1ece86b213c331852d0a3dc659b
    Would remove /workspaces/training/nextflow-run/work/be/8157808f28699259b5ff37e4bf9f72
    Would remove /workspaces/training/nextflow-run/work/f3/d828b7ba6b315eb08e0f96168119f5
    Would remove /workspaces/training/nextflow-run/work/94/7a2d7ef67d81b6e66430f300816698
    Would remove /workspaces/training/nextflow-run/work/eb/e37c5c8638c45916de3d9794d6d22f
    Would remove /workspaces/training/nextflow-run/work/3f/f56a705a6370f6c615aff22d3f240a
    Would remove /workspaces/training/nextflow-run/work/25/4922fa35ba2c6087777ddce20e9f2f
    Would remove /workspaces/training/nextflow-run/work/a4/2de7bcc2f8fc974261f7280c1480f7
    Would remove /workspaces/training/nextflow-run/work/48/8c49bded446d4fcb78c8125318c5e3
    Would remove /workspaces/training/nextflow-run/work/8e/c9834a9bdd7c21de2991d959201a9d
    Would remove /workspaces/training/nextflow-run/work/d5/a9055dc8c817c6c1f436494685955c
    Would remove /workspaces/training/nextflow-run/work/27/9473df9e55e35c8869d9571b480926
    Would remove /workspaces/training/nextflow-run/work/d3/ebc066e1d844232b9f5ee0a7d87464
    Would remove /workspaces/training/nextflow-run/work/91/17ef1815ed06b8f177d0135d8ef0dc
    Would remove /workspaces/training/nextflow-run/work/12/0443f9fcff27a552d4bc73aaedf24a
    ```

Das sind 16 Task-Verzeichnisse: die 8 Tasks der `turkey`-Ausführung plus die 8 der `tux`-Ausführung – genau so viele, wie man für zwei vollständige Ausführungen dieser Vier-Prozess-Pipeline erwarten würde.
Die `elegant_panini`-Ausführung selbst und die gecachten Tasks, die die `-resume`-Ausführung davon wiederverwendet hat, bleiben unangetastet.

Deine Ausgabe wird andere Verzeichnisnamen auflisten, und wie viele Zeilen du erhältst, hängt davon ab, wie viele Ausführungen du durchgeführt hast. Wenn du keine Zeilen siehst, stimmt entweder der Run-Name nicht mit einem in deinem Log überein, oder es gibt nichts, das davor gelöscht werden könnte.

### 3.3. Mit dem Löschen fortfahren

Sobald der Probelauf korrekt aussieht, führe denselben Befehl erneut mit `-f` statt `-n` aus:

```bash
nextflow clean -before elegant_panini -f
```

??? success "Befehlsausgabe"

    ```console
    Removed /workspaces/training/nextflow-run/work/e5/9ec3fe24deb5d1ffa478dff5e5e1d4
    Removed /workspaces/training/nextflow-run/work/75/36b1ece86b213c331852d0a3dc659b
    Removed /workspaces/training/nextflow-run/work/be/8157808f28699259b5ff37e4bf9f72
    Removed /workspaces/training/nextflow-run/work/f3/d828b7ba6b315eb08e0f96168119f5
    Removed /workspaces/training/nextflow-run/work/94/7a2d7ef67d81b6e66430f300816698
    Removed /workspaces/training/nextflow-run/work/eb/e37c5c8638c45916de3d9794d6d22f
    Removed /workspaces/training/nextflow-run/work/3f/f56a705a6370f6c615aff22d3f240a
    Removed /workspaces/training/nextflow-run/work/25/4922fa35ba2c6087777ddce20e9f2f
    Removed /workspaces/training/nextflow-run/work/a4/2de7bcc2f8fc974261f7280c1480f7
    Removed /workspaces/training/nextflow-run/work/48/8c49bded446d4fcb78c8125318c5e3
    Removed /workspaces/training/nextflow-run/work/8e/c9834a9bdd7c21de2991d959201a9d
    Removed /workspaces/training/nextflow-run/work/d5/a9055dc8c817c6c1f436494685955c
    Removed /workspaces/training/nextflow-run/work/27/9473df9e55e35c8869d9571b480926
    Removed /workspaces/training/nextflow-run/work/d3/ebc066e1d844232b9f5ee0a7d87464
    Removed /workspaces/training/nextflow-run/work/91/17ef1815ed06b8f177d0135d8ef0dc
    Removed /workspaces/training/nextflow-run/work/12/0443f9fcff27a552d4bc73aaedf24a
    ```

`nextflow clean` leert die Task-Verzeichnisse, lässt aber die zweistelligen übergeordneten Verzeichnisse (wie `e5/`) bestehen.

!!! warning "Warnung"

    Das Löschen von work-Verzeichnissen vergangener Ausführungen entfernt sie aus dem Nextflow-Cache und löscht alle Ausgaben, die nur dort gespeichert sind.
    Dadurch verliert Nextflow die Fähigkeit, die Ausführung fortzusetzen, ohne die entsprechenden Prozesse erneut auszuführen. Bereinige daher nur Ausführungen, von denen du sicher bist, dass du sie nicht mehr fortsetzen musst.
    Das ist auch der Grund, warum es sich lohnt, alles Wichtige mit `mode 'copy'` nach `results/` zu veröffentlichen, anstatt sich auf das `work/`-Verzeichnis oder einen `symlink`-Publish-Modus zu verlassen.

### Fazit

Du weißt, wie du alte work-Verzeichnisse mit `nextflow clean` entfernen kannst, und warum du damit die Möglichkeit aufgibst, von diesen Ausführungen fortzusetzen.

### Wie geht es weiter?

Lerne in [Teil 4](./04_remote_repositories.md), wie du Pipelines direkt aus Remote-Repositories wie GitHub ausführen kannst.

---

## Zusammenfassung

In diesem Teil hast du gelernt:

- Einen HTML-Ausführungsbericht mit `-with-report` und eine Ausführungs-Timeline mit `-with-timeline` zu erstellen
- Die Historie vergangener Ausführungen mit `nextflow log` einzusehen
- Alte work-Verzeichnisse mit `nextflow clean` zu entfernen und den damit verbundenen Kompromiss beim Fortsetzen zu verstehen
