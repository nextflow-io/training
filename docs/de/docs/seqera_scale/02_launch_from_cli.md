# Teil 2: Pipelines von der Kommandozeile starten

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } KI-gestützte Übersetzung - [mehr erfahren & Verbesserungen vorschlagen](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

In [Teil 1](./01_run_with_seqera.md) hast du nf-core/rnaseq über die Seqera-Weboberfläche gestartet.
Jetzt machen wir dasselbe von der Kommandozeile aus mit dem `tw` CLI und fügen deinem Workspace eine neue Pipeline hinzu.

---

## 1. Pipelines von der Kommandozeile starten

Klicke in der Run-Ansicht auf den Tab **Command line**.
Du siehst den genauen `nextflow run`-Befehl, den die Platform erstellt und in deinem Namen übermittelt hat – genau die Art von Befehl, den du im Use nf-core-Kurs manuell ausgeführt hast.

Die Platform ersetzt Nextflow nicht; sie orchestriert es.
Alles, was du über die Weboberfläche tun kannst, kannst du auch von einem Terminal aus mit dem `tw` CLI erledigen – dem Kommandozeilenwerkzeug für die Interaktion mit der Platform API.
Das ist nützlich, um Starts aus Skripten oder CI/CD-Pipelines zu automatisieren.

Wir machen das jetzt aus demselben Codespace, den du für die früheren Kurse verwendet hast.

### 1.1. Das tw CLI installieren

Führe die folgenden Befehle in deinem Codespace-Terminal aus, um das `tw`-Binary herunterzuladen und zu installieren:

```bash
curl -fsSL https://github.com/seqeralabs/tower-cli/releases/latest/download/tw-linux-x86_64 -o tw
chmod +x tw
sudo mv tw /usr/local/bin/
```

Überprüfe die Installation:

```bash
tw --version
```

??? success "Befehlsausgabe"

    ```console
    Tower CLI version 0.40.0 (build 26db579)
    ```

Das `tw` CLI ist installiert und bereit zur Konfiguration.

### 1.2. Einen Zugriffstoken erstellen

Das `tw` CLI authentifiziert sich bei Seqera mit einem persönlichen Zugriffstoken.

1. Klicke in der Seqera-Weboberfläche auf deinen Avatar in der oberen rechten Ecke und wähle **Your tokens**.
2. Klicke auf **Add token**, gib ihm einen Namen (z. B. `training`) und klicke auf **Add**.
3. Kopiere den Token-Wert – er wird nur einmal angezeigt.
   Wenn du ihn nicht sofort irgendwo speicherst, musst du einen neuen erstellen.

### 1.3. Das CLI konfigurieren

Der Einfachheit halber richten wir eine Konfigurationsdatei ein, die den
soeben erstellten Zugriffstoken und die Workspace-ID enthält.

Öffne die Datei `.seqera_config` in diesem Verzeichnis im Editor und setze die beiden Variablen:

- **`TOWER_ACCESS_TOKEN`**: der Token, den du in Abschnitt 1.2 erstellt hast
- **`TOWER_WORKSPACE_ID`**: die numerische ID deines Workspace (die Spalte `ID` in `tw workspaces list`, den du in Abschnitt 1.4 ausführst)

Sobald die Werte eingetragen sind, lade die Konfiguration:

```bash
source .seqera_config
```

Überprüfe die Verbindung:

```bash
tw info
```

??? success "Befehlsausgabe"

    ```console
        Details
    -------------------------+-----------------------------
     Tower API endpoint      | https://api.cloud.seqera.io
     Tower API version       | 1.150.0
     Tower version           | 26.1.0-cycle54
     CLI version             | 0.30.0 (fde9dec)
     CLI minimum API version | 1.148.0
     Authenticated user      | <your-name>

    System health status
    ---------------------------------------+----
     Remote API server connection check    | OK
     Tower API version check               | OK
     Authentication API credential's token | OK
    ```

Das `tw` CLI ist jetzt authentifiziert und mit deinem Seqera-Konto verbunden.
Führe `source .seqera_config` zu Beginn jeder Codespace-Sitzung aus, um die Konfiguration neu zu laden.

!!! tip "Tipp"

    Wenn dein Workspace keine primäre Compute-Umgebung festgelegt hat, kannst du `export TOWER_COMPUTE_ENV=<compute-env-name>` zu deiner Konfigurationsdatei hinzufügen, um eine Standardumgebung festzulegen.
    Jeder Konfigurationswert kann auf der Kommandozeile überschrieben werden, indem das Flag explizit übergeben wird (z. B. `--compute-env other-env`).
    Siehe die [tw CLI-Referenz](https://docs.seqera.io/platform/latest/cli/reference) für die vollständige Liste der Optionen und Umgebungsvariablen.

### 1.4. Deinen Workspace vom CLI aus erkunden

Liste die Workspaces auf, auf die du Zugriff hast:

```bash
tw workspaces list
```

??? success "Befehlsausgabe"

    ```console
    Available workspaces:
    ID              | Name            | Full Name                | Visibility
    --------------- | --------------- | ------------------------ | ----------
    <workspace-id>  | my-workspace    | my-org/my-workspace      | PRIVATE
    ```

Zeige die Runs in deinem Workspace an, einschließlich des nf-core/rnaseq-Runs, den du gerade gestartet hast:

```bash
tw runs list
```

??? success "Befehlsausgabe"

    ```console
    Pipeline runs for my-org/my-workspace workspace:
    ID        | Status   | Name          | Pipeline               | Run name
    --------- | -------- | ------------- | ---------------------- | --------
    <run-id>  | RUNNING  | ...           | nf-core/rnaseq         | happy_curie
    ```

Derselbe Run, den du in der Weboberfläche überwachst, ist hier sichtbar.

!!! note "Hinweis"

    Da `TOWER_WORKSPACE_ID` in `.seqera_config` gesetzt ist, kannst du `--workspace` bei allen `tw`-Befehlen weglassen.
    Ohne die Konfiguration würdest du es explizit übergeben:

    ```bash
    tw runs list --workspace <org>/<workspace>
    ```

Alles, was in der Weboberfläche sichtbar ist, ist auch vom CLI aus zugänglich.

### 1.5. nf-core/rnaseq vom CLI starten

Die Pipeline, die du in [Teil 1](./01_run_with_seqera.md) zu deinem Workspace hinzugefügt hast, ist im CLI über ihren Namen verfügbar.
Starte sie mit dem `test`-Profil:

```bash
tw launch nf-core-rnaseq --profile test
```

??? success "Befehlsausgabe"

    ```console
    Launching pipeline nf-core-rnaseq
    Run name: focused_einstein
    https://cloud.seqera.io/orgs/my-org/workspaces/my-workspace/watch/<run-id>
    ```

Öffne den Link in deinem Browser und bestätige, dass der Run im **Runs**-Panel erscheint.

Sobald du siehst, dass er läuft, hast du bestätigt, dass das CLI und die Weboberfläche zwei Ansichten desselben Workspace sind.

!!! note "Hinweis"

    Du kannst auch eine vollständige GitHub-URL direkt an `tw launch` übergeben, ohne die Pipeline vorher zu einem Workspace hinzuzufügen.
    Es ist jedoch generell besser, die Pipeline vor dem Start explizit hinzuzufügen: So wird die Pipeline-Konfiguration für zukünftige Runs gespeichert, sie ist über ihren Namen verfügbar und für alle Workspace-Mitglieder im Launchpad sichtbar.

    Es ist möglich, eine Pipeline direkt von der Kommandozeile aus mit `tw` zu einem Workspace hinzuzufügen.
    Der nächste Abschnitt zeigt, wie das mit der nf-core/demo-Pipeline geht.

### Fazit

Du weißt jetzt, wie du das `tw` CLI authentifizierst, deinen Workspace inspizierst und eine gespeicherte Pipeline vom Terminal aus startest.

### Wie geht es weiter?

Füge eine neue Pipeline von der Kommandozeile aus zu deinem Workspace hinzu und starte sie.

---

## 2. Eine neue Pipeline hinzufügen und ausführen

Jede Nextflow-Pipeline auf GitHub kann mit `tw pipelines add` zu deinem Workspace hinzugefügt werden, solange sie einen `main.nf`-Einstiegspunkt und eine `nextflow.config` im Stammverzeichnis hat.
nf-core/demo ist ein gutes Beispiel zum Üben: Du hast sie bereits im Use nf-core-Kurs ausgeführt, weißt also, was sie tut und was du erwarten kannst.

### 2.1. nf-core/demo zu deinem Workspace hinzufügen

Führe den folgenden Befehl aus, um die Pipeline in deinem Workspace zu registrieren:

```bash
tw pipelines add \
    --name nf-core-demo \
    https://github.com/nf-core/demo
```

??? success "Befehlsausgabe"

    ```console
    New pipeline 'nf-core-demo' added at my-org/my-workspace workspace
    ```

Die Pipeline ist jetzt registriert und erscheint im Launchpad.

### 2.2. Überprüfen, ob sie im Launchpad erscheint

Liste die Pipelines in deinem Workspace auf, um zu bestätigen, dass sie hinzugefügt wurde:

```bash
tw pipelines list
```

??? success "Befehlsausgabe"

    ```console
    Pipelines at my-org/my-workspace:
    ID   | Name            | Repository
    ---- | --------------- | -----------------------------------------
    ...  | nf-core-demo    | https://github.com/nf-core/demo
    ...  | nf-core-rnaseq  | ...
    ```

Öffne deinen Workspace im Browser und klicke auf **Launchpad**, um zu bestätigen, dass nf-core/demo jetzt neben nf-core/rnaseq erscheint.

!!! tip "Tipp"

    Du kannst Pipelines auch über die Weboberfläche hinzufügen: Klicke in der linken Seitenleiste auf **Launchpad**, dann auf **Add pipeline** und fülle das Formular entsprechend aus.

Klicke auf die Schaltfläche **Launch** beim nf-core/demo-Eintrag, um das Startformular zu öffnen.
Du wirst sehen, dass die Parameter `input` und `outdir` rot hervorgehoben sind – sie sind Pflichtfelder ohne Standardwerte, da `tw pipelines add` nur die Pipeline-Quelle registriert, ohne Parameter vorzukonfigurieren.
Die nächsten beiden Abschnitte zeigen, wie du diese Werte angibst: zuerst über das Webformular, dann von der Kommandozeile aus.

### 2.3. nf-core/demo über die Weboberfläche starten

Fülle bei geöffnetem Startformular die beiden Pflichtparameter aus.

Gib für `input` die URL des Test-Samplesheets aus dem nf-core/demo-Testprofil ein.
Du findest sie in `conf/test.config` im Pipeline-Repository, das du im Use nf-core-Kurs untersucht hast:

```
https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/samplesheet/samplesheet_test_illumina_amplicon.csv
```

Gib für `outdir` einen Cloud-Speicherpfad ein, in den die Pipeline ihre Ergebnisse schreiben kann.
Verwende den für deinen Workspace konfigurierten Bucket mit einem Unterverzeichnis, um Runs übersichtlich zu halten:

```
s3://my-bucket/demo-results
```

Sobald beide Felder ausgefüllt sind, klicke auf die blaue Schaltfläche **Launch**.

Der Run erscheint im **Runs**-Panel und sollte mit dem Testdatensatz in wenigen Minuten abgeschlossen sein.
Klicke auf den Run, um die Aufgabentabelle und eventuelle Ausführungsberichte zu erkunden.

### 2.4. nf-core/demo vom CLI starten

Im Gegensatz zu `nextflow run` akzeptiert der Befehl `tw launch` keine einzelnen Parameter-Flags wie `--input` oder `--outdir`.
Parameter müssen über eine Datei im YAML- oder JSON-Format angegeben werden, die mit `--params-file` übergeben wird.
Das fördert die Reproduzierbarkeit: Eine gespeicherte Parameterdatei dokumentiert genau, welche Werte für einen Run verwendet wurden, und macht es einfach, eine Run-Konfiguration zu wiederholen oder zu teilen.

Erstelle eine Parameterdatei in deinem Arbeitsverzeichnis:

```bash
touch params.yaml
```

Öffne sie im Editor und füge den Ausgabepfad hinzu:

```yaml title="params.yaml"
outdir: "s3://my-bucket/demo-results-cli"
```

Jetzt kannst du die Pipeline mit dem `test`-Profil (das das `input`-Samplesheet bereitstellt) und der Parameterdatei (die `outdir` bereitstellt) starten:

```bash
tw launch nf-core-demo --profile test --params-file params.yaml
```

??? success "Befehlsausgabe"

    ```console
    Launching pipeline nf-core-demo
    Run name: focused_feynman
    https://cloud.seqera.io/orgs/my-org/workspaces/my-workspace/watch/<run-id>
    ```

Öffne den Link, um zu bestätigen, dass der Run im **Runs**-Panel erscheint.

!!! tip "Tipp"

    Du kannst die Parameterdatei auch beim ersten Einrichtungsschritt angeben, wenn du einige Standardwerte festlegen möchtest, sowie einige zusätzliche Eigenschaften, die dem entsprechen, was wir zuvor über das Webformular gemacht haben:

    ```bash
    tw pipelines add \
      --name nf-core-demo-gg \
      --description "Demo pipeline with defaults for testing" \
      --params-file params.yaml \
      --profile test \
      https://github.com/nf-core/demo
    ```

### Fazit

Du weißt jetzt, wie du jede auf GitHub gehostete Nextflow-Pipeline zu deinem Workspace hinzufügst und startest – sowohl über die Weboberfläche durch manuelles Ausfüllen der Parameter als auch über das `tw` CLI durch die Kombination eines Profils mit einer Parameterdatei.

---

## Zusammenfassung

In diesem Teil hast du gelernt, wie du:

- Das `tw` CLI authentifizierst und eine gespeicherte Pipeline vom Terminal aus startest
- Eine neue Pipeline von GitHub über das CLI hinzufügst und überprüfst, ob sie im Launchpad erscheint
- Eine Pipeline über die Seqera-Weboberfläche startest, indem du die Pflichtparameter manuell ausfüllst
- Eine Pipeline über das CLI mit einem Nextflow-Profil und einer Parameterdatei startest
