# Teil 1: Weitere Container

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } KI-gestützte Übersetzung - [mehr erfahren & Verbesserungen vorschlagen](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

[TODO]

---

## 1. Wie man Container Images findet oder erstellt

Einige Softwareentwickler stellen Container Images für ihre Software bereit, die auf Container-Registries wie Docker Hub verfügbar sind, aber viele tun dies nicht.
In diesem optionalen Abschnitt zeigen wir dir zwei Möglichkeiten, wie du ein Container Image für Tools erhalten kannst, die du in deinen Nextflow Pipelines verwenden möchtest: die Verwendung von Seqera Containers und das selbstständige Erstellen des Container Images.

Du wirst ein Container Image für das `quote` pip-Paket erhalten/erstellen, das in der Übung am Ende dieses Abschnitts verwendet wird.

### 1.1. Ein Container Image von Seqera Containers erhalten

Seqera Containers ist ein kostenloser Service, der Container Images für pip- und conda-installierbare Tools (einschließlich bioconda) erstellt.
Navigiere zu [Seqera Containers](https://www.seqera.io/containers/) und suche nach dem `quote` pip-Paket.

![Seqera Containers](img/seqera-containers-1.png)

Klicke auf „+Add" und dann auf „Get Container", um ein Container Image für das `quote` pip-Paket anzufordern.

![Seqera Containers](img/seqera-containers-2.png)

Wenn dies das erste Mal ist, dass ein Community-Container für diese Version des Pakets erstellt wird, kann es einige Minuten dauern, bis er fertiggestellt ist.
Klicke, um die URI (z. B. `community.wave.seqera.io/library/pip_quote:ae07804021465ee9`) des für dich erstellten Container Images zu kopieren.

Du kannst nun das Container Image verwenden, um den `quote`-Befehl auszuführen und einen zufälligen Ausspruch von Grace Hopper zu erhalten.

```bash
docker run --rm community.wave.seqera.io/library/pip_quote:ae07804021465ee9 quote "Grace Hopper"
```

Ausgabe:

```console title="Ausgabe"
Humans are allergic to change. They love to say, 'We've always done it
this way.' I try to fight that. That's why I have a clock on my wall
that runs counter-clockwise.
```

### 1.2. Das Container Image selbst erstellen

Lass uns einige Build-Details von der Seqera Containers-Website verwenden, um das Container Image für das `quote` pip-Paket selbst zu erstellen.
Kehre zur Seqera Containers-Website zurück und klicke auf die Schaltfläche „Build Details".

Das erste Element, das wir uns ansehen, ist das `Dockerfile`, eine Art Skript-Datei, die alle Befehle enthält, die zum Erstellen des Container Images benötigt werden.
Wir haben dem untenstehenden Dockerfile einige erklärende Kommentare hinzugefügt, um dir zu helfen zu verstehen, was jeder Teil macht.

```Dockerfile title="Dockerfile"
# Starte vom micromamba-Basis-Docker-Image
FROM mambaorg/micromamba:1.5.10-noble
# Kopiere die Datei conda.yml in den Container
COPY --chown=$MAMBA_USER:$MAMBA_USER conda.yml /tmp/conda.yml
# Installiere verschiedene Utilities für Nextflow und die Pakete in der Datei conda.yml
RUN micromamba install -y -n base -f /tmp/conda.yml \
    && micromamba install -y -n base conda-forge::procps-ng \
    && micromamba env export --name base --explicit > environment.lock \
    && echo ">> CONDA_LOCK_START" \
    && cat environment.lock \
    && echo "<< CONDA_LOCK_END" \
    && micromamba clean -a -y
# Führe den Container als Root-Benutzer aus
USER root
# Setze die PATH-Umgebungsvariable, um das micromamba-Installationsverzeichnis einzuschließen
ENV PATH="$MAMBA_ROOT_PREFIX/bin:$PATH"
```

Das zweite Element, das wir uns ansehen, ist die Datei `conda.yml`, die die Liste der Pakete enthält, die im Container Image installiert werden müssen.

```conda.yml title="conda.yml"
channels:
- conda-forge
- bioconda
dependencies:
- pip
- pip:
  - quote==3.0.0 #
```

Kopiere den Inhalt dieser Dateien in die Stubs, die sich im Verzeichnis `containers/build` befinden, und führe dann den folgenden Befehl aus, um das Container Image selbst zu erstellen.

!!! note "Hinweis"

    Wir verwenden das Flag `-t quote:latest`, um das Container Image mit dem Namen `quote` und dem Tag `latest` zu versehen.
    Wir können diesen Tag verwenden, um uns auf das Container Image zu beziehen, wenn wir es auf diesem System ausführen.

```bash
docker build -t quote:latest containers/build
```

Nachdem es fertig erstellt wurde, kannst du das gerade erstellte Container Image ausführen.

```bash
docker run --rm quote:latest quote "Margaret Oakley Dayhoff"
```

### Zusammenfassung

Du hast zwei verschiedene Möglichkeiten kennengelernt, ein Container Image für ein Tool zu erhalten, das du in deinen Nextflow Pipelines verwenden möchtest: die Verwendung von Seqera Containers und das selbstständige Erstellen des Container Images.

### Wie geht es weiter?

Du hast alles, was du brauchst, um mit dem [nächsten Kapitel](./04_hello_genomics.md) dieser Trainingsreihe fortzufahren.
Du kannst auch mit einer optionalen Übung fortfahren, um Zitate über Computer-/Biologie-Pioniere mit dem `quote`-Container abzurufen und sie mit dem `cowsay`-Container auszugeben.

---

## 2. Lass die Kuh berühmte Wissenschaftler\*innen zitieren

Dieser Abschnitt enthält einige zusätzliche Übungen, um das bisher Gelernte zu üben.
Die Durchführung dieser Übungen ist _nicht erforderlich_, um spätere Teile des Trainings zu verstehen, bietet aber eine unterhaltsame Möglichkeit, dein Wissen zu festigen, indem du herausfindest, wie du die Kuh berühmte Wissenschaftler\*innen zitieren lassen kannst.

```console title="cowsay-output-Grace-Hopper.txt"
  _________________________________________________
 /                                                 \
| Humans are allergic to change. They love to       |
| say, 'We've always done it this way.' I try to fi |
| ght that. That's why I have a clock on my wall th |
| at runs counter-clockwise.                        |
| -Grace Hopper                                     |
 \                                                 /
  =================================================
                                                 \
                                                  \
                                                    ^__^
                                                    (oo)\_______
                                                    (__)\       )\/\
                                                        ||----w |
                                                        ||     ||
```

### 2.1. Modifiziere das `hello-containers.nf`-Skript, um einen getQuote-Prozess zu verwenden

Wir haben eine Liste von Computer- und Biologie-Pionieren in der Datei `containers/data/pioneers.csv`.
Auf hoher Ebene musst du für diese Übung Folgendes tun:

- Modifiziere die Standard-`params.input_file`, um auf die Datei `pioneers.csv` zu verweisen.
- Erstelle einen `getQuote`-Prozess, der den `quote`-Container verwendet, um ein Zitat für jede Eingabe abzurufen.
- Verbinde die Ausgabe des `getQuote`-Prozesses mit dem `cowsay`-Prozess, um das Zitat anzuzeigen.

Für das `quote`-Container Image kannst du entweder das verwenden, das du selbst in der vorherigen zusätzlichen Übung erstellt hast, oder das, das du von Seqera Containers erhalten hast.

!!! tip "Tipp"

    Eine gute Wahl für den `script`-Block deines getQuote-Prozesses könnte sein:
        ```groovy
        script:
            def safe_author = author.tokenize(' ').join('-')
            """
            quote "$author" > quote-${safe_author}.txt
            echo "-${author}" >> quote-${safe_author}.txt
            """
        ```

Du findest eine Lösung für diese Übung in `containers/solutions/hello-containers-4.1.nf`.

### 2.2. Modifiziere deine Nextflow-Pipeline, um sie im `quote`- und `sayHello`-Modus ausführen zu können.

Füge deiner Pipeline etwas Verzweigungslogik hinzu, damit sie Eingaben akzeptieren kann, die sowohl für `quote` als auch für `sayHello` vorgesehen sind.
Hier ist ein Beispiel, wie man eine `if`-Anweisung in einem Nextflow Workflow verwendet:

```groovy title="hello-containers.nf"
workflow {
    if (params.quote) {
        ...
    }
    else {
        ...
    }
    cowSay(text_ch)
}
```

!!! tip "Tipp"

    Du kannst `new_ch = processName.out` verwenden, um der Ausgabe eines Prozesses einen Namen zuzuweisen.

Du findest eine Lösung für diese Übung in `containers/solutions/hello-containers-4.2.nf`.

### Zusammenfassung

Du weißt, wie man Container in Nextflow verwendet, um Prozesse auszuführen, und wie man etwas Verzweigungslogik in deine Pipelines einbaut!

### Wie geht es weiter?

Feiere, mach eine Pause und trink etwas Wasser!

Wenn du bereit bist, fahre mit Teil 3 dieser Trainingsreihe fort, um zu lernen, wie du das bisher Gelernte auf einen realistischeren Datenanalyse-Anwendungsfall anwenden kannst.
