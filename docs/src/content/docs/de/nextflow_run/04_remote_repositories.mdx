# Teil 4: Remote-Pipelines ausführen

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } KI-gestützte Übersetzung - [mehr erfahren & Verbesserungen vorschlagen](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Bisher hast du Workflow-Skripte ausgeführt, die lokal gespeichert waren.
In der Praxis möchtest du oft Pipelines aus Remote-Repositories wie GitHub ausführen, ohne sie selbst herunterzuladen.

Nextflow macht das unkompliziert: Du kannst jede Pipeline direkt aus einer Git-Repository-URL ausführen.

---

## 1. Eine Pipeline von GitHub ausführen

Die grundlegende Syntax zum Ausführen einer Remote-Pipeline lautet `nextflow run <repository>`, wobei `<repository>` ein GitHub-Repository-Pfad wie `nextflow-io/hello`, eine vollständige URL oder ein Pfad zu GitLab, Bitbucket oder einem anderen Git-Hosting-Dienst sein kann.

### 1.1. Die Pipeline starten

Führe die offizielle Nextflow-Demo-Pipeline „hello" aus.
Das ist eine andere, viel einfachere Pipeline als die, die du in diesem Kurs verwendet hast: Sie ist älter als die „Hello"-Pipeline, die in diesem Training verwendet wird, und gibt einfach eine Begrüßung in einigen fest codierten Sprachen aus. Erwarte also weder die CSV-Eingabe noch die ASCII-Art, die du gewohnt bist.

```bash
nextflow run nextflow-io/hello
```

??? success "Befehlsausgabe"

    ```console
    N E X T F L O W   ~  version 26.04.4

    Pulling nextflow-io/hello ...
     downloaded from https://github.com/nextflow-io/hello.git
    Launching `https://github.com/nextflow-io/hello` [sleepy_swanson] revision: 3c2cdc9823 [master]

    executor >  local (4)
    [ba/08236d] sayHello (4) | 4 of 4 ✔
    Ciao world!

    Hello world!

    Bonjour world!

    Hola world!
    ```

### 1.2. Herausfinden, wo die Pipeline gecacht wird

Beim ersten Ausführen einer Remote-Pipeline lädt Nextflow sie herunter und speichert sie lokal im Cache.
Spätere Ausführungen verwenden die gecachte Version, es sei denn, du forderst explizit ein Update an.

Standardmäßig speichert Nextflow heruntergeladene Pipelines unter `$NXF_HOME/assets`.
Um herauszufinden, wo eine bestimmte Pipeline gespeichert wurde und welche Revisionen verfügbar sind, frage Nextflow direkt:

```bash
nextflow info nextflow-io/hello
```

??? success "Befehlsausgabe"

    ```console
     project name: nextflow-io/hello
     repository  : https://github.com/nextflow-io/hello
     local path  : /workspaces/.nextflow/assets/.repos/nextflow-io/hello
     main script : main.nf
     revisions   :
     > master (default)
       mybranch
       testing
       v1.1 [t]
       v1.2 [t]
       v1.3 [t]
    ```

    Nextflow markiert jede Revision, die du bereits lokal ausgecheckt hast, mit `>`; die übrigen sind verfügbar, wurden aber noch nicht in eine Arbeitskopie geholt.

Du kannst auch alle bisher heruntergeladenen Pipelines mit `nextflow list` auflisten:

```bash
nextflow list
```

??? success "Befehlsausgabe"

    ```console
    nextflow-io/hello
    ```

Der Kurs [Use nf-core](../nfcore_use/01_run_demo.md#12-retrieve-the-pipeline-code) behandelt diesen Caching-Mechanismus ausführlicher, einschließlich der Möglichkeit, den Quellcode einer heruntergeladenen Pipeline zu durchsuchen.

### Fazit

Du weißt jetzt, wie du eine Pipeline direkt aus einem GitHub-Repository ausführst, ohne sie selbst herunterzuladen, und wo du sie danach lokal findest.

### Wie geht es weiter?

Lerne, wie du eine bestimmte Version einer Remote-Pipeline für Reproduzierbarkeit festlegst.

---

## 2. Eine Version für Reproduzierbarkeit festlegen

Standardmäßig führt Nextflow die neueste Revision des Standard-Branches aus.
Du kannst eine bestimmte Version (Tag), einen Branch oder einen Commit mit dem Flag `-r` festlegen.

### 2.1. Eine bestimmte Revision festlegen

```bash
nextflow run nextflow-io/hello -r v1.3
```

??? success "Befehlsausgabe"

    ```console
    N E X T F L O W   ~  version 26.04.4

    Pulling nextflow-io/hello:v1.3 ...
     downloaded from https://github.com/nextflow-io/hello.git
    Launching `https://github.com/nextflow-io/hello` [sick_carson] revision: 2ce0b0e294 [v1.3]

    executor >  local (4)
    [61/e11f77] sayHello (4) | 4 of 4 ✔
    Ciao world!

    Bonjour world!

    Hello world!

    Hola world!
    ```

Nextflow lädt diese Revision beim ersten Aufruf herunter – daher die Zeilen `Pulling` und `downloaded from`. Bei einem erneuten Aufruf derselben Revision wird direkt zu `Launching` gesprungen.
Eine genaue Revision festzulegen ist entscheidend für Reproduzierbarkeit.
Es stellt sicher, dass du und deine Mitarbeiter\*innen exakt denselben Pipeline-Code ausführen, unabhängig davon, was sich seitdem im Repository geändert hat.

### 2.2. Revisionen gelten nur für den jeweiligen Aufruf

Eine Revision mit `-r` festzulegen wirkt sich nur auf den Aufruf aus, bei dem du sie angibst: Es ändert nicht, was ein späterer, einfacher `nextflow run`-Befehl verwendet.
Führe die Pipeline erneut ohne `-r` aus:

```bash
nextflow run nextflow-io/hello
```

??? success "Befehlsausgabe"

    ```console
    N E X T F L O W   ~  version 26.04.4

    Launching `https://github.com/nextflow-io/hello` [hungry_maxwell] revision: 3c2cdc9823 [master]

    executor >  local (4)
    [ba/08236d] sayHello (4) | 4 of 4 ✔
    Ciao world!

    Hello world!

    Bonjour world!

    Hola world!
    ```

Obwohl der vorherige Aufruf explizit `v1.3` festgelegt hat, kehrt dieser Aufruf direkt zum Standard-Branch (`master`) zurück.
Nextflow behält für jede verwendete Revision eine separate lokale Arbeitskopie, was die `>`-Markierungen in `nextflow info` zeigen – aber es merkt sich nie, welche du zuletzt ausgeführt hast.
Den Namen des Standard-Branches einer Pipeline findest du mit `nextflow info <pipeline>`; er ist mit `(default)` markiert.
Reproduzierbarkeit liegt ganz bei dir: Gib `-r` immer explizit an, wenn es darauf ankommt, anstatt davon auszugehen, dass eine in einem früheren Aufruf festgelegte Revision noch gilt.

### Fazit

Du weißt jetzt, wie du eine Remote-Pipeline auf eine bestimmte Version, einen Branch oder einen Commit festlegst, und dass diese Festlegung nur für den jeweiligen Aufruf gilt, nicht für spätere Ausführungen.

### Wie geht es weiter?

Du hast die Grundlagen zum Ausführen und Verwalten von Nextflow-Pipelines kennengelernt.
Unter [Kursübersicht](next_steps.md) erfährst du, wie es weitergeht.

---

## Zusammenfassung

In diesem Teil hast du gelernt:

- Eine Pipeline direkt aus einem GitHub-Repository auszuführen, ohne sie herunterzuladen
- Eine Remote-Pipeline auf eine bestimmte Revision festzulegen für reproduzierbare Ausführung, und zu verstehen, dass diese Festlegung nur für den jeweiligen Aufruf gilt
