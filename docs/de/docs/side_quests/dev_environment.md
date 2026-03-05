# Entwicklungsumgebung

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } KI-gestützte Übersetzung - [mehr erfahren & Verbesserungen vorschlagen](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Moderne Integrierte Entwicklungsumgebungen (IDEs) können deine Nextflow-Entwicklungserfahrung dramatisch verändern. Dieser Side Quest konzentriert sich speziell darauf, VS Code und seine Nextflow-Erweiterung zu nutzen, um Code schneller zu schreiben, Fehler frühzeitig zu erkennen und effizient in komplexen Workflows zu navigieren.

!!! note "Dies ist kein traditionelles Tutorial"

    Anders als andere Trainingsmodule ist dieser Leitfaden als Sammlung von schnellen Hinweisen, Tipps und praktischen Beispielen organisiert, nicht als Schritt-für-Schritt-Tutorial. Jeder Abschnitt kann unabhängig basierend auf deinen Interessen und aktuellen Entwicklungsbedürfnissen erkundet werden. Spring gerne herum und konzentriere dich auf die Features, die für deine Workflow-Entwicklung am unmittelbarsten nützlich sein werden.

## Was du zuerst wissen solltest

Dieser Leitfaden setzt voraus, dass du den [Hello Nextflow](../hello_nextflow/)-Trainingskurs abgeschlossen hast und mit grundlegenden Nextflow-Konzepten vertraut bist, einschließlich:

- **Grundlegende Workflow-Struktur**: Verständnis von Processes, Workflows und wie sie zusammenarbeiten
- **Channel-Operationen**: Channels erstellen, Daten zwischen Processes weitergeben und grundlegende Operatoren verwenden
- **Module und Organisation**: Wiederverwendbare Module erstellen und include-Anweisungen verwenden
- **Grundlagen der Konfiguration**: `nextflow.config` für Parameter, Process-Direktiven und Profile verwenden

## Was du hier lernen wirst

Dieser Leitfaden konzentriert sich auf **IDE-Produktivitätsfeatures**, die dich zu einer effizienteren Nextflow-Entwickler\*in machen werden:

- **Erweitertes Syntax-Highlighting**: Verstehen, was VS Code dir über deine Codestruktur zeigt
- **Intelligente Auto-Vervollständigung**: Kontextbezogene Vorschläge für schnelleres Code-Schreiben nutzen
- **Fehlererkennung und Diagnose**: Syntaxfehler erkennen, bevor du deinen Workflow ausführst
- **Code-Navigation**: Schnell zwischen Processes, Modulen und Definitionen wechseln
- **Formatierung und Organisation**: Konsistenten, lesbaren Code-Stil beibehalten
- **KI-unterstützte Entwicklung** (optional): Moderne KI-Tools nutzen, die in deine IDE integriert sind

!!! info "Warum IDE-Features jetzt?"

    Du hast wahrscheinlich bereits VS Code während des [Hello Nextflow](../hello_nextflow/)-Kurses verwendet, aber wir haben uns auf das Erlernen der Nextflow-Grundlagen konzentriert anstatt auf IDE-Features. Jetzt, da du mit grundlegenden Nextflow-Konzepten wie Processes, Workflows, Channels und Modulen vertraut bist, bist du bereit, die ausgefeilten IDE-Features zu nutzen, die dich zu einer effizienteren Entwickler*in machen werden.

    Denke daran als "Level-up" deiner Entwicklungsumgebung - der gleiche Editor, den du verwendet hast, hat viel leistungsfähigere Funktionen, die wirklich wertvoll werden, sobald du verstehst, wobei sie dir helfen.

---

## 0. Setup und Aufwärmen

Lass uns einen Arbeitsbereich speziell für die Erkundung von IDE-Features einrichten:

```bash title="Navigiere zum IDE-Features-Verzeichnis"
cd side-quests/ide_features
```

Öffne dieses Verzeichnis in VS Code:

```bash title="Öffne VS Code im aktuellen Verzeichnis"
code .
```

Das `ide_features`-Verzeichnis enthält Beispiel-Workflows, die verschiedene IDE-Features demonstrieren:

```bash title="Zeige Verzeichnisstruktur"
tree .
```

```console title="Projektstruktur"
tree .
.
├── basic_workflow.nf
├── complex_workflow.nf
├── data
│   ├── sample_001.fastq.gz
│   ├── sample_002.fastq.gz
│   ├── sample_003.fastq.gz
│   ├── sample_004.fastq.gz
│   ├── sample_005.fastq.gz
│   └── sample_data.csv
├── modules
│   ├── fastqc.nf
│   ├── star.nf
│   └── utils.nf
└── nextflow.config

3 directories, 12 files
```

!!! note "Über die Beispieldateien"

    - `basic_workflow.nf` ist ein funktionierender einfacher Workflow, den du ausführen und ändern kannst
    - `complex_workflow.nf` ist nur zur Illustration gedacht, um Navigationsfunktionen zu demonstrieren - er läuft möglicherweise nicht erfolgreich, zeigt aber eine realistische Multi-Datei-Workflow-Struktur

### Tastaturkürzel

Einige der Features in diesem Leitfaden verwenden optionale Tastaturkürzel. Du greifst möglicherweise über GitHub Codespaces im Browser auf dieses Material zu, und in diesem Fall funktionieren die Tastaturkürzel manchmal nicht wie erwartet, weil sie für andere Dinge in deinem System verwendet werden.

Wenn du VS Code lokal ausführst, wie du es wahrscheinlich tun wirst, wenn du tatsächlich Workflows schreibst, funktionieren die Tastaturkürzel wie beschrieben.

Wenn du einen Mac verwendest, verwenden einige (nicht alle) Tastaturkürzel "cmd" anstelle von "ctrl", und wir werden dies im Text wie `Ctrl/Cmd` angeben.

### 0.1. Installation der Nextflow-Erweiterung

!!! note "Verwendest du bereits Devcontainer?"

    Wenn du in **GitHub Codespaces** arbeitest oder einen **lokalen Devcontainer** verwendest, ist die Nextflow-Erweiterung wahrscheinlich bereits für dich installiert und konfiguriert. Du kannst die manuellen Installationsschritte unten überspringen und direkt mit der Erkundung der Erweiterungsfeatures fortfahren.

Um die Erweiterung manuell zu installieren:

1. Öffne VS Code
2. Gehe zur Extensions-Ansicht, indem du auf das Extensions-Symbol links klickst: ![Extensions-Symbol](img/extensions_icon.png) (Tastaturkürzel `Ctrl/Cmd+Shift+X`, wenn du VSCode lokal ausführst)
3. Suche nach "Nextflow"
4. Installiere die offizielle Nextflow-Erweiterung

![Nextflow-Erweiterung installieren](img/install_extension.png)

### 0.2. Arbeitsbereich-Layout

Da du VS Code während Hello Nextflow verwendet hast, bist du bereits mit den Grundlagen vertraut. So organisierst du deinen Arbeitsbereich effizient für diese Session:

- **Editor-Bereich**: Zum Anzeigen und Bearbeiten von Dateien. Du kannst diesen in mehrere Bereiche aufteilen, um Dateien nebeneinander zu vergleichen.
- **Datei-Explorer** Klick (![Datei-Explorer-Symbol](img/files_icon.png)) (`Ctrl/Cmd+Shift+E`): Die lokalen Dateien und Ordner auf deinem System. Halte dies links offen, um zwischen Dateien zu navigieren
- **Integriertes Terminal** (`Ctrl+Shift+` Backtick für Windows und MacOS): Ein Terminal zur Interaktion mit dem Computer unten. Verwende dies, um Nextflow oder andere Befehle auszuführen.
- **Problems-Panel** (`Ctrl+Shift+M`): VS Code zeigt hier alle Fehler und Probleme an, die es erkennt. Dies ist nützlich, um Probleme auf einen Blick hervorzuheben.

Du kannst Panels herumziehen oder ausblenden (`Ctrl/Cmd+B` zum Umschalten der Seitenleiste), um dein Layout anzupassen, während wir die Beispiele durchgehen.

### Zusammenfassung

Du hast VS Code mit der Nextflow-Erweiterung eingerichtet und verstehst das Arbeitsbereich-Layout für effiziente Entwicklung.

### Was kommt als Nächstes?

Lerne, wie Syntax-Highlighting dir hilft, die Nextflow-Codestruktur auf einen Blick zu verstehen.

---

## 1. Syntax-Highlighting und Code-Struktur

Jetzt, da dein Arbeitsbereich eingerichtet ist, lass uns erkunden, wie das Syntax-Highlighting von VS Code dir hilft, Nextflow-Code effektiver zu lesen und zu schreiben.

### 1.1. Nextflow-Syntax-Elemente

Öffne `basic_workflow.nf`, um Syntax-Highlighting in Aktion zu sehen:

![Syntax Showcase](img/syntax_showcase.png)

Beachte, wie VS Code hervorhebt:

- **Schlüsselwörter** (`process`, `workflow`, `input`, `output`, `script`) in unterschiedlichen Farben
- **String-Literale** und **Parameter** mit unterschiedlicher Gestaltung
- **Kommentare** in einer gedämpften Farbe
- **Variablen** und **Funktionsaufrufe** mit angemessener Betonung
- **Codeblöcke** mit korrekten Einrückungslinien

!!! note "Theme-abhängige Farben"

    Die spezifischen Farben, die du siehst, hängen von deinem VS Code-Theme (Dark/Light-Modus), Farbeinstellungen und allen Anpassungen ab, die du vorgenommen hast. Das Wichtige ist, dass verschiedene Syntaxelemente visuell voneinander unterschieden werden, was die Codestruktur leichter verständlich macht, unabhängig von deinem gewählten Farbschema.

### 1.2. Verständnis der Code-Struktur

Das Syntax-Highlighting hilft dir, schnell zu identifizieren:

- **Process-Grenzen**: Klare Unterscheidung zwischen verschiedenen Processes
- **Input/Output-Blöcke**: Einfach zu erkennende Datenfluss-Definitionen
- **Script-Blöcke**: Die tatsächlich ausgeführten Befehle
- **Channel-Operationen**: Datentransformationsschritte
- **Konfigurationsdirektiven**: Process-spezifische Einstellungen

Diese visuelle Organisation wird unschätzbar wertvoll, wenn du mit komplexen Workflows arbeitest, die mehrere Processes und komplizierte Datenflüsse enthalten.

### Zusammenfassung

Du verstehst, wie das Syntax-Highlighting von VS Code dir hilft, die Nextflow-Codestruktur zu lesen und verschiedene Sprachelemente für schnellere Entwicklung zu identifizieren.

### Was kommt als Nächstes?

Lerne, wie intelligente Auto-Vervollständigung das Code-Schreiben mit kontextbezogenen Vorschlägen beschleunigt.

---

## 2. Intelligente Auto-Vervollständigung

Die Auto-Vervollständigungsfeatures von VS Code helfen dir, Code schneller und mit weniger Fehlern zu schreiben, indem sie basierend auf dem Kontext passende Optionen vorschlagen.

### 2.1. Kontextbezogene Vorschläge

Die Auto-Vervollständigungsoptionen variieren je nachdem, wo du dich in deinem Code befindest:

#### Channel-Operationen

Öffne `basic_workflow.nf` erneut und versuche, `channel.` im workflow-Block einzugeben:

![Channel-Auto-Vervollständigung](img/autocomplete_channel.png)

Du wirst Vorschläge sehen für:

- `fromPath()` - Channel aus Dateipfaden erstellen
- `fromFilePairs()` - Channel aus gepaarten Dateien erstellen
- `of()` - Channel aus Werten erstellen
- `fromSRA()` - Channel aus SRA-Accessions erstellen
- Und viele mehr...

Dies hilft dir, schnell die richtige Channel-Factory zu finden, ohne dir genaue Methodennamen merken zu müssen.

Du kannst auch die Operatoren entdecken, die auf Channels angewendet werden können. Gib zum Beispiel `FASTQC.out.html.` ein, um verfügbare Operationen zu sehen:

![Channel-Operationen-Auto-Vervollständigung](img/autocomplete_operators.png)

#### Process-Direktiven

Innerhalb eines Process-Script-Blocks, gib `task.` ein, um verfügbare Runtime-Eigenschaften zu sehen:

![Task-Eigenschaften-Auto-Vervollständigung](img/autocomplete_task.png)

#### Konfiguration

Öffne nextflow.config und gib `process.` irgendwo ein, um verfügbare Process-Direktiven zu sehen:

![Config-Auto-Vervollständigung](img/autocomplete_config.png)

Du wirst Vorschläge sehen für:

- `executor`
- `memory`
- `cpus`

Dies spart Zeit beim Konfigurieren von Processes und funktioniert über verschiedene Konfigurationsbereiche hinweg. Versuche zum Beispiel, `docker.` einzugeben, um Docker-spezifische Konfigurationsoptionen zu sehen.

### Zusammenfassung

Du kannst die intelligente Auto-Vervollständigung von VS Code verwenden, um verfügbare Channel-Operationen, Process-Direktiven und Konfigurationsoptionen zu entdecken, ohne Syntax auswendig zu lernen.

### Was kommt als Nächstes?

Lerne, wie Echtzeit-Fehlererkennung dir hilft, Probleme zu erkennen, bevor du deinen Workflow ausführst, einfach durch das Lesen des Codes.

## 3. Fehlererkennung und Diagnose

Die Echtzeit-Fehlererkennung von VS Code hilft dir, Probleme zu erkennen, bevor du deinen Workflow ausführst.

### 3.1. Syntaxfehler-Erkennung

Lass uns einen absichtlichen Fehler erstellen, um die Erkennung in Aktion zu sehen. Öffne `basic_workflow.nf` und ändere den Process-Namen von `FASTQC` zu `FASTQ` (oder einen anderen ungültigen Namen). VS Code wird den Fehler im workflow-Block sofort mit einer roten Wellenlinie hervorheben:

![Fehlerunterstreichung](img/error_underline.png)

### 3.2. Problems-Panel

Über die individuelle Fehlerhervorhebung hinaus bietet VS Code ein zentrales Problems-Panel, das alle Fehler, Warnungen und Informationsmeldungen über deinen Arbeitsbereich hinweg aggregiert. Öffne es mit `Ctrl/Cmd+Shift+M` und verwende das Filter-Symbol, um nur Fehler anzuzeigen, die für die aktuelle Datei relevant sind:

![Problems-Panel filtern](img/active_file.png)

Klicke auf ein Problem, um direkt zur problematischen Zeile zu springen

![Problems-Panel](img/problems_panel.png)

Behebe den Fehler, indem du den Process-Namen zurück zu `FASTQC` änderst.

### 3.3. Häufige Fehlermuster

Häufige Fehler in der Nextflow-Syntax umfassen:

- **Fehlende Klammern**: Nicht übereinstimmende `{` oder `}`
- **Unvollständige Blöcke**: Fehlende erforderliche Abschnitte in Processes
- **Ungültige Syntax**: Fehlerhafte Nextflow DSL
- **Tippfehler in Schlüsselwörtern**: Falsch geschriebene Process-Direktiven
- **Channel-Nichtübereinstimmungen**: Typ-Inkompatibilitäten

Der Nextflow Language Server hebt diese Probleme im Problems-Panel hervor. Du kannst diese frühzeitig überprüfen, um Syntaxfehler beim Ausführen einer Pipeline zu vermeiden.

### Zusammenfassung

Du kannst die Fehlererkennung und das Problems-Panel von VS Code verwenden, um Syntaxfehler und Probleme zu erkennen, bevor du deinen Workflow ausführst, was Zeit spart und Frustration verhindert.

### Was kommt als Nächstes?

Lerne, wie du effizient zwischen Processes, Modulen und Definitionen in komplexen Workflows navigierst.

---

## 4. Code-Navigation und Symbol-Verwaltung

Effiziente Navigation ist entscheidend, wenn du mit komplexen Workflows arbeitest, die sich über mehrere Dateien erstrecken. Um dies zu verstehen, ersetze die Process-Definition in `basic_workflow.nf` durch einen Import für das Modul, das wir dir bereitgestellt haben:

=== "Nachher"

    ```groovy title="basic_workflow.nf" linenums="3"
    include { FASTQC } from './modules/fastqc.nf'
    ```

=== "Vorher"

    ```groovy title="basic_workflow.nf" linenums="3"
    process FASTQC {
        tag "${sample_id}"
        publishDir "${params.output_dir}/fastqc", mode: 'copy'

        input:
        tuple val(sample_id), path(reads)

        output:
        tuple val(sample_id), path("*.html"), emit: html
        tuple val(sample_id), path("*.zip"), emit: zip

        script:
        def args = task.ext.args ?: ''
        """
        fastqc \\
            ${args} \\
            --threads ${task.cpus} \\
            ${reads}
        """
    }
    ```

### 4.1. Gehe zu Definition

Wenn du mit der Maus über einen Process-Namen wie `FASTQC` fährst, siehst du ein Popup mit dem Modul-Interface (Eingaben und Ausgaben):

![Gehe zu Definition](img/syntax.png)

Diese Funktion ist besonders wertvoll beim Erstellen von Workflows, da sie dir ermöglicht, das Modul-Interface zu verstehen, ohne die Moduldatei direkt zu öffnen.

Du kannst schnell zu jeder Process-, Modul- oder Variablendefinition navigieren, indem du **Ctrl/Cmd-Klick** verwendest. Fahre mit der Maus über den Link zur Moduldatei oben im Script und folge dem Link wie vorgeschlagen:

![Link folgen](img/follow_link.png)

Das Gleiche funktioniert für Process-Namen. Gehe zurück zu `basic_workflow.nf` und probiere dies beim `FASTQC`-Process-Namen im workflow-Block aus. Dies verlinkt dich direkt zum Process-Namen (der in diesem Beispiel mit der Moduldatei identisch ist, aber Teil einer viel größeren Datei sein könnte).

Um dorthin zurückzugehen, wo du warst, verwende **Alt+←** (oder **Ctrl+-** auf Mac). Dies ist eine leistungsstarke Möglichkeit, Code zu erkunden, ohne deinen Platz zu verlieren.

Lass uns nun die Navigation in einem komplexeren Workflow mit `complex_workflow.nf` erkunden (die zuvor erwähnte reine Illustrationsdatei). Dieser Workflow enthält mehrere Processes, die in separaten Moduldateien definiert sind, sowie einige Inline-Processes. Während komplexe Multi-Datei-Strukturen manuell schwer zu navigieren sein können, macht die Fähigkeit, zu Definitionen zu springen, die Erkundung viel handhabbarer.

1. Öffne `complex_workflow.nf`
2. Navigiere zu Modul-Definitionen
3. Verwende **Alt+←** (oder **Ctrl+-**), um zurückzunavigieren
4. Navigiere zum `FASTQC`-Process-Namen im workflow-Block. Dies verlinkt dich direkt zum Process-Namen (der in diesem Beispiel mit der Moduldatei identisch ist, aber Teil einer viel größeren Datei sein könnte).
5. Navigiere wieder zurück
6. Navigiere zum `TRIM_GALORE`-Process im workflow-Block. Dieser ist inline definiert, sodass er dich nicht zu einer separaten Datei bringt, aber er zeigt dir trotzdem die Process-Definition, und du kannst immer noch dorthin zurücknavigieren, wo du warst.

### 4.2. Symbol-Navigation

Während `complex_workflow.nf` noch geöffnet ist, kannst du eine Übersicht über alle Symbole in der Datei erhalten, indem du `@` in die Suchleiste oben in VSCode eingibst (das Tastaturkürzel ist `Ctrl/Cmd+Shift+O`, funktioniert aber möglicherweise nicht in Codespaces). Dies öffnet das Symbol-Navigationspanel, das alle Symbole in der aktuellen Datei auflistet:

![Symbol-Navigation](img/symbols.png)

Dies zeigt:

- Alle Process-Definitionen
- Workflow-Definitionen (in dieser Datei sind zwei Workflows definiert)
- Funktionsdefinitionen

Beginne mit dem Tippen, um Ergebnisse zu filtern.

### 4.3. Alle Referenzen finden

Zu verstehen, wo ein Process oder eine Variable in deiner Codebasis verwendet wird, kann sehr hilfreich sein. Wenn du beispielsweise alle Referenzen zum `FASTQC`-Process finden möchtest, beginne damit, zu seiner Definition zu navigieren. Du kannst dies tun, indem du `modules/fastqc.nf` direkt öffnest oder VS Codes Schnellnavigationsfunktion mit `Ctrl/Cmd-Klick` verwendest, wie wir es oben getan haben. Sobald du bei der Process-Definition bist, klicke mit der rechten Maustaste auf den `FASTQC`-Process-Namen und wähle "Find All References" aus dem Kontextmenü, um alle Instanzen zu sehen, wo er verwendet wird.

![Referenzen finden](img/references.png)

Diese Funktion zeigt alle Instanzen an, wo `FASTQC` innerhalb deines Arbeitsbereichs referenziert wird, einschließlich seiner Verwendung in den beiden unterschiedlichen Workflows. Diese Einsicht ist entscheidend, um die potenziellen Auswirkungen von Änderungen am `FASTQC`-Process zu beurteilen.

### 4.4. Outline-Panel

Das Outline-Panel, das sich in der Explorer-Seitenleiste befindet (klicke ![Explorer-Symbol](img/files_icon.png)), bietet eine bequeme Übersicht über alle Symbole in deiner aktuellen Datei. Diese Funktion ermöglicht es dir, schnell zu navigieren und die Struktur deines Codes zu verwalten, indem Funktionen, Variablen und andere Schlüsselelemente in einer hierarchischen Ansicht angezeigt werden.

![Outline-Panel](img/outline.png)

Verwende das Outline-Panel, um schnell zu verschiedenen Teilen deines Codes zu navigieren, ohne den Datei-Browser zu verwenden.

### 4.5. DAG-Visualisierung

Die Nextflow-Erweiterung von VS Code kann deinen Workflow als Directed Acyclic Graph (DAG) visualisieren. Dies hilft dir, den Datenfluss und die Abhängigkeiten zwischen Processes zu verstehen. Öffne `complex_workflow.nf` und klicke auf die "Preview DAG"-Schaltfläche über `workflow {` (der zweite `workflow`-Block in dieser Datei):

![DAG-Vorschau](img/dag_preview.png)

Dies ist nur der 'Entry'-Workflow, aber du kannst auch den DAG für die inneren Workflows in der Vorschau anzeigen, indem du auf die "Preview DAG"-Schaltfläche über dem workflow `RNASEQ_PIPELINE {` weiter oben klickst:

![DAG-Vorschau innerer Workflow](img/dag_preview_inner.png)

Für diesen Workflow kannst du die Knoten im DAG verwenden, um zu den entsprechenden Process-Definitionen im Code zu navigieren. Klicke auf einen Knoten, und er bringt dich zur relevanten Process-Definition im Editor. Besonders wenn ein Workflow zu einer großen Größe anwächst, kann dies dir wirklich helfen, im Code zu navigieren und zu verstehen, wie die Processes verbunden sind.

### Zusammenfassung

Du kannst komplexe Workflows effizient navigieren, indem du Gehe-zu-Definition, Symbolsuche, Referenzen finden und DAG-Visualisierung verwendest, um Codestruktur und Abhängigkeiten zu verstehen.

### Was kommt als Nächstes?

Lerne, wie du effektiv über mehrere miteinander verbundene Dateien in größeren Nextflow-Projekten hinweg arbeitest.

## 5. Arbeiten über mehrere Dateien hinweg

Echte Nextflow-Entwicklung beinhaltet die Arbeit mit mehreren miteinander verbundenen Dateien. Lass uns erkunden, wie VS Code dir hilft, komplexe Projekte effizient zu verwalten.

### 5.1. Schnelle Datei-Navigation

Während `complex_workflow.nf` geöffnet ist, wirst du bemerken, dass es mehrere Module importiert. Lass uns die schnelle Navigation zwischen ihnen üben.

Drücke **Ctrl+P** (oder **Cmd+P**) und beginne "fast" einzugeben:

VS Code zeigt dir passende Dateien an. Wähle `modules/fastqc.nf`, um sofort dorthin zu springen. Dies ist viel schneller als durch den Datei-Explorer zu klicken, wenn du ungefähr weißt, welche Datei du suchst.

Probiere dies mit anderen Mustern aus:

- Gib "star" ein, um die STAR-Alignment-Moduldatei zu finden (`star.nf`)
- Gib "utils" ein, um die Utility-Funktionsdatei zu finden (`utils.nf`)
- Gib "config" ein, um zu Konfigurationsdateien zu springen (`nextflow.config`)

### 5.2. Split-Editor für Multi-Datei-Entwicklung

Wenn du mit Modulen arbeitest, musst du oft sowohl den Haupt-Workflow als auch die Modul-Definitionen gleichzeitig sehen. Lass uns das einrichten:

1. Öffne `complex_workflow.nf`
2. Öffne `modules/fastqc.nf` in einem neuen Tab
3. Rechtsklicke auf den `modules/fastqc.nf`-Tab und wähle "Split Right"
4. Jetzt kannst du beide Dateien nebeneinander sehen

![Split-Editor](img/split_editor.png)

Dies ist unschätzbar wertvoll, wenn:

- Du Modul-Interfaces überprüfst, während du Workflow-Aufrufe schreibst, und die Vorschau nicht ausreicht
- Du ähnliche Processes über verschiedene Module hinweg vergleichst
- Du den Datenfluss zwischen Workflow und Modulen debuggst

### 5.3. Projektweite Suche

Manchmal musst du finden, wo bestimmte Muster in deinem gesamten Projekt verwendet werden. Drücke `Ctrl/Cmd+Shift+F`, um das Suchpanel zu öffnen.

Versuche, über den Arbeitsbereich hinweg nach `publishDir` zu suchen:

![Projektsuche](img/project_search.png)

Dies zeigt dir jede Datei, die Publish-Verzeichnisse verwendet, und hilft dir:

- Ausgabe-Organisationsmuster zu verstehen
- Beispiele für bestimmte Direktiven zu finden
- Konsistenz über Module hinweg sicherzustellen

### Zusammenfassung

Du kannst komplexe Multi-Datei-Projekte verwalten, indem du schnelle Datei-Navigation, Split-Editoren und projektweite Suche verwendest, um effizient über Workflows und Module hinweg zu arbeiten.

### Was kommt als Nächstes?

Lerne, wie Code-Formatierungs- und Wartungsfunktionen deine Workflows organisiert und lesbar halten.

---

## 6. Code-Formatierung und Wartung

Richtige Code-Formatierung ist nicht nur für die Ästhetik wichtig, sondern auch für die Verbesserung der Lesbarkeit, des Verständnisses und der Leichtigkeit der Aktualisierung komplexer Workflows.

### 6.1. Automatische Formatierung in Aktion

Öffne `basic_workflow.nf` und bringe die Formatierung absichtlich durcheinander:

- Entferne einige Einrückungen: Markiere das gesamte Dokument und drücke `shift+tab` viele Male, um so viele Einrückungen wie möglich zu entfernen.
- Füge zusätzliche Leerzeichen an zufälligen Stellen hinzu: Bei der `channel.fromPath`-Anweisung füge 30 Leerzeichen nach der `(` hinzu.
- Breche einige Zeilen unbeholfen: Füge eine neue Zeile zwischen dem `.view {`-Operator und dem `Processing sample:`-String hinzu, füge aber keine entsprechende neue Zeile vor der schließenden Klammer `}` hinzu.

Drücke jetzt `Shift+Alt+F` (oder `Shift+Option+F` auf MacOS) zum Auto-Formatieren:

VS Code:

- Behebt die Einrückung, um die Process-Struktur klar zu zeigen
- Richtet ähnliche Elemente konsistent aus
- Entfernt unnötige Leerzeichen
- Behält lesbare Zeilenumbrüche bei

Beachte, dass die automatische Formatierung möglicherweise nicht jedes Code-Style-Problem löst. Der Nextflow Language Server zielt darauf ab, deinen Code ordentlich zu halten, respektiert aber auch deine persönlichen Vorlieben in bestimmten Bereichen. Wenn du beispielsweise die Einrückung innerhalb des `script`-Blocks eines Processes entfernst, lässt der Formatierer sie so, da du diesen Stil möglicherweise absichtlich bevorzugst.

Derzeit gibt es keine strikte Style-Durchsetzung für Nextflow, daher bietet der Language Server etwas Flexibilität. Er wird jedoch konsequent Formatierungsregeln um Methoden- und Funktionsdefinitionen herum anwenden, um Klarheit zu erhalten.

### 6.2. Code-Organisationsfunktionen

#### Schnelles Kommentieren

Wähle einen Codeblock in deinem Workflow aus und drücke **Ctrl+/** (oder **Cmd+/**), um ihn auszukommentieren:

```groovy
// workflow {
//     ch_input = channel.fromPath(params.input)
//         .splitCsv(header: true)
//         .map { row -> [row.sample_id, file(row.fastq_path)] }
//
//     FASTQC(ch_input)
// }
```

Dies ist perfekt für:

- Vorübergehendes Deaktivieren von Teilen von Workflows während der Entwicklung
- Hinzufügen erklärender Kommentare zu komplexen Channel-Operationen
- Dokumentieren von Workflow-Abschnitten

Verwende **Ctrl+/** (oder **Cmd+/**) erneut, um den Code zu entkommentieren.

#### Code-Faltung für Übersicht

Beachte in `complex_workflow.nf` die kleinen Pfeile neben Process-Definitionen. Klicke darauf, um Processes zu falten (zusammenklappen):

![Code-Faltung](img/code_folding.png)

Dies gibt dir eine High-Level-Übersicht über deine Workflow-Struktur, ohne dich in Implementierungsdetails zu verlieren.

#### Klammer-Zuordnung

Platziere deinen Cursor neben einer beliebigen `{` oder `}`-Klammer, und VS Code hebt die passende Klammer hervor. Verwende **Ctrl+Shift+\\** (oder **Cmd+Shift+\\**), um zwischen passenden Klammern zu springen.

Dies ist entscheidend für:

- Verständnis von Process-Grenzen
- Finden fehlender oder zusätzlicher Klammern
- Navigation in verschachtelten Workflow-Strukturen

#### Mehrzeilige Auswahl und Bearbeitung

Für die gleichzeitige Bearbeitung mehrerer Zeilen bietet VS Code leistungsstarke Multi-Cursor-Fähigkeiten:

- **Mehrzeilige Auswahl**: Halte **Ctrl+Alt** (oder **Cmd+Option** für MacOS) und verwende Pfeiltasten, um mehrere Zeilen auszuwählen
- **Mehrzeilige Einrückung**: Wähle mehrere Zeilen aus und verwende **Tab** zum Einrücken oder **Shift+Tab** zum Ausrücken ganzer Blöcke

Dies ist besonders nützlich für:

- Konsistentes Einrücken ganzer Process-Blöcke
- Gleichzeitiges Hinzufügen von Kommentaren zu mehreren Zeilen
- Bearbeiten ähnlicher Parameter-Definitionen über mehrere Processes hinweg

### Zusammenfassung

Du kannst sauberen, lesbaren Code mit automatischer Formatierung, Kommentierungsfunktionen, Code-Faltung, Klammer-Zuordnung und mehrzeiliger Bearbeitung pflegen, um komplexe Workflows effizient zu organisieren.

### Was kommt als Nächstes?

Lerne, wie VS Code mit deinem breiteren Entwicklungs-Workflow über das bloße Bearbeiten von Code hinaus integriert.

---

## 7. Entwicklungs-Workflow-Integration

VS Code integriert sich gut mit deinem Entwicklungs-Workflow über das bloße Bearbeiten von Code hinaus.

### 7.1. Versionskontroll-Integration

!!! note "Codespaces und Git-Integration"

    Wenn du in **GitHub Codespaces** arbeitest, funktionieren einige Git-Integrationsfunktionen möglicherweise nicht wie erwartet, insbesondere Tastaturkürzel für Source Control. Du hast möglicherweise auch während des ersten Setups abgelehnt, das Verzeichnis als Git-Repository zu öffnen, was für Trainingszwecke in Ordnung ist.

Wenn dein Projekt ein Git-Repository ist (wie dieses), zeigt VS Code:

- Geänderte Dateien mit farbigen Indikatoren
- Git-Status in der Statusleiste
- Inline-Diff-Ansichten
- Commit- und Push-Funktionen

Öffne das Source Control-Panel mit der Source Control-Schaltfläche (![Source Control-Symbol](img/source_control_icon.png)) (`Ctrl+Shift+G` oder `Cmd+Shift+G`, wenn du lokal mit VSCode arbeitest), um Git-Änderungen zu sehen und Commits direkt im Editor zu stagen.

![Source Control-Panel](img/source_control.png)

### 7.2. Workflows ausführen und inspizieren

Lass uns einen Workflow ausführen und dann die Ergebnisse inspizieren. Im integrierten Terminal (`Ctrl+Shift+` Backtick in Windows und MacOS) führe den einfachen Workflow aus:

```bash title="Führe den einfachen Workflow aus"
nextflow run basic_workflow.nf --input data/sample_data.csv --output_dir results
```

Während der Workflow läuft, siehst du Echtzeit-Ausgaben im Terminal. Nach Abschluss kannst du VS Code verwenden, um Ergebnisse zu inspizieren, ohne deinen Editor zu verlassen:

1. **Navigiere zu work-Verzeichnissen**: Verwende den Datei-Explorer oder das Terminal, um `.nextflow/work` zu durchsuchen
2. **Öffne Log-Dateien**: Klicke auf Log-Dateipfade in der Terminal-Ausgabe, um sie direkt in VS Code zu öffnen
3. **Inspiziere Ausgaben**: Durchsuche veröffentlichte Ergebnisverzeichnisse im Datei-Explorer
4. **Zeige Ausführungsberichte an**: Öffne HTML-Berichte direkt in VS Code oder deinem Browser

Dies hält alles an einem Ort, anstatt zwischen mehreren Anwendungen zu wechseln.

### Zusammenfassung

Du kannst VS Code mit Versionskontrolle und Workflow-Ausführung integrieren, um deinen gesamten Entwicklungsprozess von einer einzigen Oberfläche aus zu verwalten.

### Was kommt als Nächstes?

Sieh, wie all diese IDE-Features in deinem täglichen Entwicklungs-Workflow zusammenarbeiten.

---

## 8. Zusammenfassung und Schnellnotizen

Hier sind einige Schnellnotizen zu jedem der oben besprochenen IDE-Features:

### 8.1. Starten eines neuen Features

1. **Schnelles Datei-Öffnen** (`Ctrl+P` oder `Cmd+P`), um relevante existierende Module zu finden
2. **Split-Editor**, um ähnliche Processes nebeneinander anzuzeigen
3. **Symbol-Navigation** (`Ctrl+Shift+O` oder `Cmd+Shift+O`), um die Dateistruktur zu verstehen
4. **Auto-Vervollständigung**, um neuen Code schnell zu schreiben

### 8.2. Debugging von Problemen

1. **Problems-Panel** (`Ctrl+Shift+M` oder `Cmd+Shift+M`), um alle Fehler auf einmal zu sehen
2. **Gehe zu Definition** (`Ctrl-Klick` oder `Cmd-Klick`), um Process-Interfaces zu verstehen
3. **Alle Referenzen finden**, um zu sehen, wie Processes verwendet werden
4. **Projektweite Suche**, um ähnliche Muster oder Probleme zu finden

### 8.3. Refactoring und Verbesserung

1. **Projektweite Suche** (`Ctrl+Shift+F` oder `Cmd+Shift+F`), um Muster zu finden
2. **Auto-Formatierung** (`Shift+Alt+F` oder `Shift+Option+F`), um Konsistenz zu wahren
3. **Code-Faltung**, um dich auf Struktur zu konzentrieren
4. **Git-Integration**, um Änderungen zu verfolgen

---

## Zusammenfassung

Du hast jetzt eine Schnelltour durch die IDE-Features von VS Code für Nextflow-Entwicklung erhalten. Diese Tools werden dich deutlich produktiver machen, indem sie:

- **Fehler reduzieren** durch Echtzeit-Syntaxprüfung
- **Entwicklung beschleunigen** mit intelligenter Auto-Vervollständigung
- **Navigation verbessern** in komplexen Multi-Datei-Workflows
- **Qualität erhalten** durch konsistente Formatierung
- **Verständnis verbessern** durch erweitertes Highlighting und Strukturvisualisierung

Wir erwarten nicht, dass du dich an alles erinnerst, aber jetzt weißt du, dass diese Features existieren, und du wirst sie finden können, wenn du sie brauchst. Während du weiterhin Nextflow-Workflows entwickelst, werden diese IDE-Features zur zweiten Natur, sodass du dich auf das Schreiben von hochwertigem Code konzentrieren kannst, anstatt mit Syntax und Struktur zu ringen.

### Was kommt als Nächstes?

Wende diese IDE-Fähigkeiten an, während du andere Trainingsmodule durcharbeitest, zum Beispiel:

- **[nf-test](nf-test.md)**: Erstelle umfassende Test-Suites für deine Workflows
- **[Hello nf-core](../../hello_nf-core/)**: Erstelle produktionsreife Pipelines mit Community-Standards

Die wahre Kraft dieser IDE-Features zeigt sich, wenn du an größeren, komplexeren Projekten arbeitest. Beginne damit, sie schrittweise in deinen Workflow zu integrieren - innerhalb weniger Sessions werden sie zur zweiten Natur werden und transformieren, wie du an Nextflow-Entwicklung herangehst.

Vom Erkennen von Fehlern, bevor sie dich verlangsamen, bis zum einfachen Navigieren in komplexen Codebasen - diese Tools werden dich zu einer selbstbewussteren und effizienteren Entwickler\*in machen.

Frohes Coden!
