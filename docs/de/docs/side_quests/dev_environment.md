# Entwicklungsumgebung

Moderne Integrierte Entwicklungsumgebungen (IDEs) können deine Nextflow-Entwicklung grundlegend verändern. Dieser Side Quest konzentriert sich speziell darauf, VS Code und seine Nextflow-Erweiterung zu nutzen, um Code schneller zu schreiben, Fehler frühzeitig zu erkennen und effizient durch komplexe Workflows zu navigieren.

!!! note "Dies ist kein traditionelles Tutorial"

    Anders als andere Trainingsmodule ist dieser Leitfaden als Sammlung von schnellen Hinweisen, Tipps und praktischen Beispielen organisiert, nicht als Schritt-für-Schritt-Anleitung. Jeder Abschnitt kann unabhängig erkundet werden, basierend auf deinen Interessen und aktuellen Entwicklungsbedürfnissen. Spring gerne herum und konzentriere dich auf die Features, die für deine Workflow-Entwicklung am nützlichsten sind.

## Was du zuerst wissen solltest

Dieser Leitfaden setzt voraus, dass du den Kurs [Hello Nextflow](../hello_nextflow/) abgeschlossen hast und mit grundlegenden Nextflow-Konzepten vertraut bist, einschließlich:

- **Grundlegende Workflow-Struktur**: Verständnis von Prozessen, Workflows und wie sie zusammenhängen
- **Kanal-Operationen**: Kanäle erstellen, Daten zwischen Prozessen übergeben und grundlegende Operatoren verwenden
- **Module und Organisation**: Wiederverwendbare Module erstellen und include-Anweisungen verwenden
- **Konfigurations-Grundlagen**: `nextflow.config` für Parameter, Prozess-Direktiven und Profile verwenden

## Was du hier lernst

Dieser Leitfaden konzentriert sich auf **IDE-Produktivitäts-Features**, die dich zu einem effizienteren Nextflow-Entwickler\*in machen:

- **Erweitertes Syntax-Highlighting**: Verstehen, was VS Code dir über deine Code-Struktur zeigt
- **Intelligente Auto-Vervollständigung**: Kontextbezogene Vorschläge für schnelleres Code-Schreiben nutzen
- **Fehlererkennung und Diagnose**: Syntaxfehler erkennen, bevor du deinen Workflow ausführst
- **Code-Navigation**: Schnell zwischen Prozessen, Modulen und Definitionen wechseln
- **Formatierung und Organisation**: Konsistenten, lesbaren Code-Stil beibehalten
- **KI-gestützte Entwicklung** (optional): Moderne KI-Tools nutzen, die in deine IDE integriert sind

!!! info "Warum jetzt IDE-Features?"

    Du hast wahrscheinlich bereits VS Code während des Kurses [Hello Nextflow](../hello_nextflow/) verwendet, aber wir haben uns auf das Lernen der Nextflow-Grundlagen konzentriert statt auf IDE-Features. Jetzt, wo du mit grundlegenden Nextflow-Konzepten wie Prozessen, Workflows, Kanälen und Modulen vertraut bist, bist du bereit, die ausgefeilten IDE-Features zu nutzen, die dich zu einem effizienteren Entwickler\*in machen.

    Betrachte dies als "Level-Up" deiner Entwicklungsumgebung - der Editor, den du bereits verwendet hast, hat viel leistungsfähigere Funktionen, die wirklich wertvoll werden, sobald du verstehst, wobei sie dir helfen.

---

## 0. Setup und Aufwärmen

Lass uns einen Workspace speziell für die Erkundung von IDE-Features einrichten:

```bash title="Navigate to the IDE features directory"
cd side-quests/ide_features
```

Öffne dieses Verzeichnis in VS Code:

```bash title="Open VS Code in current directory"
code .
```

Das `ide_features`-Verzeichnis enthält Beispiel-Workflows, die verschiedene IDE-Features demonstrieren:

```bash title="Show directory structure"
tree .
```

```console title="Project structure"
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

    - `basic_workflow.nf` ist ein funktionierender einfacher Workflow, den du ausführen und modifizieren kannst
    - `complex_workflow.nf` ist nur zur Illustration gedacht, um Navigations-Features zu demonstrieren - er läuft möglicherweise nicht erfolgreich, zeigt aber eine realistische Multi-File-Workflow-Struktur

### Tastenkombinationen

Einige Features in diesem Leitfaden verwenden optionale Tastenkombinationen. Du greifst möglicherweise über GitHub Codespaces im Browser auf dieses Material zu, und in diesem Fall funktionieren die Shortcuts manchmal nicht wie erwartet, weil sie für andere Dinge in deinem System verwendet werden.

Wenn du VS Code lokal ausführst, wie du es wahrscheinlich tun wirst, wenn du tatsächlich Workflows schreibst, funktionieren die Shortcuts wie beschrieben.

Wenn du einen Mac verwendest, nutzen einige (nicht alle) Tastenkombinationen "cmd" statt "ctrl", und wir werden dies im Text wie `Ctrl/Cmd` angeben.

### 0.1. Installation der Nextflow-Erweiterung

!!! note "Verwendest du bereits Devcontainer?"

    Wenn du in **GitHub Codespaces** arbeitest oder einen **lokalen Devcontainer** verwendest, ist die Nextflow-Erweiterung wahrscheinlich bereits installiert und für dich konfiguriert. Du kannst die manuellen Installationsschritte unten überspringen und direkt mit der Erkundung der Erweiterungs-Features fortfahren.

Um die Erweiterung manuell zu installieren:

1. Öffne VS Code
2. Gehe zur Extensions-Ansicht, indem du auf das Extensions-Symbol links klickst: ![Extensions-Symbol](img/extensions_icon.png) (Shortcut `Ctrl/Cmd+Shift+X`, wenn du VSCode lokal ausführst)
3. Suche nach "Nextflow"
4. Installiere die offizielle Nextflow-Erweiterung

![Nextflow-Erweiterung installieren](img/install_extension.png)

### 0.2. Workspace-Layout

Da du VS Code bereits während Hello Nextflow verwendet hast, kennst du bereits die Grundlagen. So organisierst du deinen Workspace effizient für diese Session:

- **Editor-Bereich**: Zum Anzeigen und Bearbeiten von Dateien. Du kannst diesen in mehrere Bereiche aufteilen, um Dateien nebeneinander zu vergleichen.
- **Datei-Explorer** klicke (![Datei-Explorer-Symbol](img/files_icon.png)) (`Ctrl/Cmd+Shift+E`): Die lokalen Dateien und Ordner auf deinem System. Halte dies links offen, um zwischen Dateien zu navigieren
- **Integriertes Terminal** (`Ctrl+Shift+` Backtick für Windows und MacOS): Ein Terminal zur Interaktion mit dem Computer unten. Verwende dies, um Nextflow oder andere Befehle auszuführen.
- **Problems-Panel** (`Ctrl+Shift+M`): VS Code zeigt hier alle erkannten Fehler und Probleme an. Dies ist nützlich, um Probleme auf einen Blick hervorzuheben.

Du kannst Panels herumziehen oder ausblenden (`Ctrl/Cmd+B` zum Umschalten der Seitenleiste), um dein Layout anzupassen, während wir die Beispiele durchgehen.

### Fazit

Du hast VS Code mit der Nextflow-Erweiterung eingerichtet und verstehst das Workspace-Layout für effiziente Entwicklung.

### Wie geht es weiter?

Lerne, wie Syntax-Highlighting dir hilft, die Nextflow-Code-Struktur auf einen Blick zu verstehen.

---

## 1. Syntax-Highlighting und Code-Struktur

Jetzt, wo dein Workspace eingerichtet ist, lass uns erkunden, wie VS Codes Syntax-Highlighting dir hilft, Nextflow-Code effektiver zu lesen und zu schreiben.

### 1.1. Nextflow-Syntax-Elemente

Öffne `basic_workflow.nf`, um Syntax-Highlighting in Aktion zu sehen:

![Syntax-Showcase](img/syntax_showcase.png)

Beachte, wie VS Code hervorhebt:

- **Schlüsselwörter** (`process`, `workflow`, `input`, `output`, `script`) in unterschiedlichen Farben
- **String-Literale** und **Parameter** mit unterschiedlichem Styling
- **Kommentare** in einer gedämpften Farbe
- **Variablen** und **Funktionsaufrufe** mit entsprechender Betonung
- **Code-Blöcke** mit korrekten Einrückungslinien

!!! note "Theme-abhängige Farben"

    Die spezifischen Farben, die du siehst, hängen von deinem VS Code-Theme (Dark/Light-Modus), Farbeinstellungen und allen Anpassungen ab, die du vorgenommen hast. Wichtig ist, dass verschiedene Syntax-Elemente visuell voneinander unterschieden werden, was die Code-Struktur unabhängig von deinem gewählten Farbschema leichter verständlich macht.

### 1.2. Code-Struktur verstehen

Das Syntax-Highlighting hilft dir, schnell zu identifizieren:

- **Prozess-Grenzen**: Klare Unterscheidung zwischen verschiedenen Prozessen
- **Input/Output-Blöcke**: Einfach zu erkennende Datenfluss-Definitionen
- **Script-Blöcke**: Die tatsächlich ausgeführten Befehle
- **Kanal-Operationen**: Datentransformations-Schritte
- **Konfigurations-Direktiven**: Prozess-spezifische Einstellungen

Diese visuelle Organisation wird unschätzbar wertvoll, wenn du mit komplexen Workflows arbeitest, die mehrere Prozesse und komplizierte Datenflüsse enthalten.

### Fazit

Du verstehst, wie VS Codes Syntax-Highlighting dir hilft, die Nextflow-Code-Struktur zu lesen und verschiedene Sprachelemente für schnellere Entwicklung zu identifizieren.

### Wie geht es weiter?

Lerne, wie intelligente Auto-Vervollständigung das Code-Schreiben mit kontextbezogenen Vorschlägen beschleunigt.

---

## 2. Intelligente Auto-Vervollständigung

VS Codes Auto-Vervollständigungs-Features helfen dir, Code schneller und mit weniger Fehlern zu schreiben, indem sie passende Optionen basierend auf dem Kontext vorschlagen.

### 2.1. Kontextbezogene Vorschläge

Die Auto-Vervollständigungs-Optionen variieren je nachdem, wo du dich in deinem Code befindest:

#### Kanal-Operationen

Öffne `basic_workflow.nf` erneut und versuche, `channel.` im Workflow-Block zu tippen:

![Kanal-Auto-Vervollständigung](img/autocomplete_channel.png)

Du siehst Vorschläge für:

- `fromPath()` - Kanal aus Dateipfaden erstellen
- `fromFilePairs()` - Kanal aus gepaarten Dateien erstellen
- `of()` - Kanal aus Werten erstellen
- `fromSRA()` - Kanal aus SRA-Accessions erstellen
- Und viele mehr...

Dies hilft dir, schnell die richtige Channel-Factory zu finden, ohne dir exakte Methodennamen merken zu müssen.

Du kannst auch die verfügbaren Operatoren entdecken, die auf Kanäle angewendet werden können. Tippe zum Beispiel `FASTQC.out.html.`, um verfügbare Operationen zu sehen:

![Kanal-Operationen-Auto-Vervollständigung](img/autocomplete_operators.png)

#### Prozess-Direktiven

Tippe innerhalb eines Prozess-Script-Blocks `task.`, um verfügbare Laufzeit-Eigenschaften zu sehen:

![Task-Eigenschaften-Auto-Vervollständigung](img/autocomplete_task.png)

#### Konfiguration

Öffne nextflow.config und tippe irgendwo `process.`, um verfügbare Prozess-Direktiven zu sehen:

![Config-Auto-Vervollständigung](img/autocomplete_config.png)

Du siehst Vorschläge für:

- `executor`
- `memory`
- `cpus`

Dies spart Zeit beim Konfigurieren von Prozessen und funktioniert über verschiedene Konfigurations-Scopes hinweg. Versuche zum Beispiel, `docker.` zu tippen, um Docker-spezifische Konfigurationsoptionen zu sehen.

### Fazit

Du kannst VS Codes intelligente Auto-Vervollständigung nutzen, um verfügbare Kanal-Operationen, Prozess-Direktiven und Konfigurationsoptionen zu entdecken, ohne Syntax auswendig zu lernen.

### Wie geht es weiter?

Lerne, wie Echtzeit-Fehlererkennung dir hilft, Probleme zu erkennen, bevor du deinen Workflow ausführst, einfach durch das Lesen des Codes.

## 3. Fehlererkennung und Diagnose

VS Codes Echtzeit-Fehlererkennung hilft dir, Probleme zu erkennen, bevor du deinen Workflow ausführst.

### 3.1. Syntaxfehler-Erkennung

Lass uns einen absichtlichen Fehler erstellen, um die Erkennung in Aktion zu sehen. Öffne `basic_workflow.nf` und ändere den Prozessnamen von `FASTQC` zu `FASTQ` (oder einem anderen ungültigen Namen). VS Code wird den Fehler im Workflow-Block sofort mit einer roten Wellenlinie hervorheben:

![Fehler-Unterstreichung](img/error_underline.png)

### 3.2. Problems-Panel

Über die individuelle Fehlerhervorhebung hinaus bietet VS Code ein zentrales Problems-Panel, das alle Fehler, Warnungen und Info-Meldungen in deinem Workspace aggregiert. Öffne es mit `Ctrl/Cmd+Shift+M` und verwende das Filter-Symbol, um nur Fehler anzuzeigen, die für die aktuelle Datei relevant sind:

![Problems-Panel filtern](img/active_file.png)

Klicke auf ein Problem, um direkt zur problematischen Zeile zu springen

![Problems-Panel](img/problems_panel.png)

Behebe den Fehler, indem du den Prozessnamen zurück zu `FASTQC` änderst.

### 3.3. Häufige Fehlermuster

Häufige Fehler in der Nextflow-Syntax umfassen:

- **Fehlende Klammern**: Nicht übereinstimmende `{` oder `}`
- **Unvollständige Blöcke**: Fehlende erforderliche Abschnitte in Prozessen
- **Ungültige Syntax**: Fehlerhafte Nextflow-DSL
- **Tippfehler in Schlüsselwörtern**: Falsch geschriebene Prozess-Direktiven
- **Kanal-Inkompatibilitäten**: Typ-Inkompatibilitäten

Der Nextflow-Language-Server hebt diese Probleme im Problems-Panel hervor. Du kannst diese frühzeitig überprüfen, um Syntaxfehler beim Ausführen einer Pipeline zu vermeiden.

### Fazit

Du kannst VS Codes Fehlererkennung und Problems-Panel nutzen, um Syntaxfehler und Probleme zu erkennen, bevor du deinen Workflow ausführst, was Zeit spart und Frustration verhindert.

### Wie geht es weiter?

Lerne, wie du effizient zwischen Prozessen, Modulen und Definitionen in komplexen Workflows navigierst.

---

## 4. Code-Navigation und Symbol-Verwaltung

Effiziente Navigation ist entscheidend, wenn du mit komplexen Workflows arbeitest, die sich über mehrere Dateien erstrecken. Um dies zu verstehen, ersetze die Prozess-Definition in `basic_workflow.nf` durch einen Import für das Modul, das wir dir bereitgestellt haben:

=== "Danach"

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

Wenn du mit der Maus über einen Prozessnamen wie `FASTQC` fährst, siehst du ein Popup mit dem Modul-Interface (Eingaben und Ausgaben):

![Gehe zu Definition](img/syntax.png)

Dieses Feature ist besonders wertvoll beim Erstellen von Workflows, da es dir ermöglicht, das Modul-Interface zu verstehen, ohne die Modul-Datei direkt zu öffnen.

Du kannst schnell zu jeder Prozess-, Modul- oder Variablen-Definition navigieren, indem du **Ctrl/Cmd-Klick** verwendest. Fahre mit der Maus über den Link zur Modul-Datei oben im Script und folge dem Link wie vorgeschlagen:

![Link folgen](img/follow_link.png)

Das Gleiche funktioniert für Prozessnamen. Gehe zurück zu `basic_workflow.nf` und probiere dies beim `FASTQC`-Prozessnamen im Workflow-Block aus. Dies verlinkt dich direkt zum Prozessnamen (der in diesem Beispiel mit der Modul-Datei identisch ist, aber auch mitten in einer viel größeren Datei sein könnte).

Um zurückzugehen, wo du warst, verwende **Alt+←** (oder **Ctrl+-** auf Mac). Dies ist eine leistungsstarke Möglichkeit, Code zu erkunden, ohne deine Position zu verlieren.

Lass uns nun die Navigation in einem komplexeren Workflow mit `complex_workflow.nf` erkunden (die zuvor erwähnte Nur-Illustrations-Datei). Dieser Workflow enthält mehrere Prozesse, die in separaten Modul-Dateien definiert sind, sowie einige Inline-Prozesse. Während komplexe Multi-File-Strukturen manuell schwierig zu navigieren sein können, macht die Fähigkeit, zu Definitionen zu springen, die Erkundung viel handhabbarer.

1. Öffne `complex_workflow.nf`
2. Navigiere zu Modul-Definitionen
3. Verwende **Alt+←** (oder **Ctrl+-**), um zurückzunavigieren
4. Navigiere zum `FASTQC`-Prozessnamen im Workflow-Block. Dies verlinkt dich direkt zum Prozessnamen (der in diesem Beispiel mit der Modul-Datei identisch ist, aber auch mitten in einer viel größeren Datei sein könnte).
5. Navigiere wieder zurück
6. Navigiere zum `TRIM_GALORE`-Prozess im Workflow-Block. Dieser ist inline definiert, also führt er dich nicht zu einer separaten Datei, aber er zeigt dir trotzdem die Prozess-Definition, und du kannst trotzdem zurücknavigieren, wo du warst.

### 4.2. Symbol-Navigation

Mit `complex_workflow.nf` noch geöffnet, kannst du eine Übersicht aller Symbole in der Datei erhalten, indem du `@` in die Suchleiste oben in VSCode tippst (die Tastenkombination ist `Ctrl/Cmd+Shift+O`, funktioniert aber möglicherweise nicht in Codespaces). Dies öffnet das Symbol-Navigations-Panel, das alle Symbole in der aktuellen Datei auflistet:

![Symbol-Navigation](img/symbols.png)

Dies zeigt:

- Alle Prozess-Definitionen
- Workflow-Definitionen (es sind zwei Workflows in dieser Datei definiert)
- Funktions-Definitionen

Beginne zu tippen, um Ergebnisse zu filtern.

### 4.3. Alle Referenzen finden

Zu verstehen, wo ein Prozess oder eine Variable in deiner Codebasis verwendet wird, kann sehr hilfreich sein. Wenn du zum Beispiel alle Referenzen zum `FASTQC`-Prozess finden möchtest, beginne damit, zu seiner Definition zu navigieren. Du kannst dies tun, indem du `modules/fastqc.nf` direkt öffnest oder VS Codes Schnellnavigations-Feature mit `Ctrl/Cmd-Klick` verwendest, wie wir es oben getan haben. Sobald du bei der Prozess-Definition bist, klicke mit der rechten Maustaste auf den `FASTQC`-Prozessnamen und wähle "Find All References" aus dem Kontextmenü, um alle Instanzen zu sehen, wo er verwendet wird.

![Referenzen finden](img/references.png)

Dieses Feature zeigt alle Instanzen, wo `FASTQC` in deinem Workspace referenziert wird, einschließlich seiner Verwendung in den zwei unterschiedlichen Workflows. Diese Einsicht ist entscheidend, um die potenzielle Auswirkung von Änderungen am `FASTQC`-Prozess zu bewerten.

### 4.4. Outline-Panel

Das Outline-Panel, das sich in der Explorer-Seitenleiste befindet (klicke ![Explorer-Symbol](img/files_icon.png)), bietet eine praktische Übersicht aller Symbole in deiner aktuellen Datei. Dieses Feature ermöglicht es dir, schnell zu navigieren und die Struktur deines Codes zu verwalten, indem es Funktionen, Variablen und andere Schlüsselelemente in einer hierarchischen Ansicht anzeigt.

![Outline-Panel](img/outline.png)

Verwende das Outline-Panel, um schnell zu verschiedenen Teilen deines Codes zu navigieren, ohne den Datei-Browser zu verwenden.

### 4.5. DAG-Visualisierung

VS Codes Nextflow-Erweiterung kann deinen Workflow als Directed Acyclic Graph (DAG) visualisieren. Dies hilft dir, den Datenfluss und die Abhängigkeiten zwischen Prozessen zu verstehen. Öffne `complex_workflow.nf` und klicke auf den "Preview DAG"-Button über `workflow {` (der zweite `workflow`-Block in dieser Datei):

![DAG-Vorschau](img/dag_preview.png)

Dies ist nur der 'Entry'-Workflow, aber du kannst auch den DAG für die inneren Workflows in der Vorschau anzeigen, indem du auf den "Preview DAG"-Button über dem Workflow `RNASEQ_PIPELINE {` weiter oben klickst:

![DAG-Vorschau innerer Workflow](img/dag_preview_inner.png)

Für diesen Workflow kannst du die Knoten im DAG verwenden, um zu den entsprechenden Prozess-Definitionen im Code zu navigieren. Klicke auf einen Knoten, und er führt dich zur relevanten Prozess-Definition im Editor. Besonders wenn ein Workflow zu einer großen Größe wächst, kann dies dir wirklich helfen, im Code zu navigieren und zu verstehen, wie die Prozesse verbunden sind.

### Fazit

Du kannst komplexe Workflows effizient navigieren, indem du Gehe-zu-Definition, Symbol-Suche, Referenzen-Suche und DAG-Visualisierung verwendest, um Code-Struktur und Abhängigkeiten zu verstehen.

### Wie geht es weiter?

Lerne, wie du effektiv über mehrere miteinander verbundene Dateien in größeren Nextflow-Projekten arbeitest.

## 5. Arbeiten über mehrere Dateien hinweg

Echte Nextflow-Entwicklung beinhaltet die Arbeit mit mehreren miteinander verbundenen Dateien. Lass uns erkunden, wie VS Code dir hilft, komplexe Projekte effizient zu verwalten.

### 5.1. Schnelle Datei-Navigation

Mit `complex_workflow.nf` geöffnet, wirst du bemerken, dass es mehrere Module importiert. Lass uns schnelle Navigation zwischen ihnen üben.

Drücke **Ctrl+P** (oder **Cmd+P**) und beginne "fast" zu tippen:

VS Code zeigt dir passende Dateien. Wähle `modules/fastqc.nf`, um sofort dorthin zu springen. Dies ist viel schneller als durch den Datei-Explorer zu klicken, wenn du ungefähr weißt, welche Datei du suchst.

Probiere dies mit anderen Mustern:

- Tippe "star", um die STAR-Alignment-Modul-Datei zu finden (`star.nf`)
- Tippe "utils", um die Utility-Funktionen-Datei zu finden (`utils.nf`)
- Tippe "config", um zu Konfigurationsdateien zu springen (`nextflow.config`)

### 5.2. Split-Editor für Multi-File-Entwicklung

Wenn du mit Modulen arbeitest, musst du oft sowohl den Haupt-Workflow als auch Modul-Definitionen gleichzeitig sehen. Lass uns dies einrichten:

1. Öffne `complex_workflow.nf`
2. Öffne `modules/fastqc.nf` in einem neuen Tab
3. Klicke mit der rechten Maustaste auf den `modules/fastqc.nf`-Tab und wähle "Split Right"
4. Jetzt kannst du beide Dateien nebeneinander sehen

![Split-Editor](img/split_editor.png)

Dies ist unschätzbar wertvoll, wenn:

- Du Modul-Interfaces überprüfst, während du Workflow-Aufrufe schreibst, und die Vorschau nicht ausreicht
- Du ähnliche Prozesse über verschiedene Module hinweg vergleichst
- Du den Datenfluss zwischen Workflow und Modulen debuggst

### 5.3. Projektweite Suche

Manchmal musst du finden, wo spezifische Muster in deinem gesamten Projekt verwendet werden. Drücke `Ctrl/Cmd+Shift+F`, um das Such-Panel zu öffnen.

Versuche, über den Workspace hinweg nach `publishDir` zu suchen:

![Projekt-Suche](img/project_search.png)

Dies zeigt dir jede Datei, die Publish-Directories verwendet, und hilft dir:

- Ausgabe-Organisations-Muster zu verstehen
- Beispiele für spezifische Direktiven zu finden
- Konsistenz über Module hinweg sicherzustellen

### Fazit

Du kannst komplexe Multi-File-Projekte verwalten, indem du schnelle Datei-Navigation, Split-Editoren und projektweite Suche verwendest, um effizient über Workflows und Module hinweg zu arbeiten.

### Wie geht es weiter?

Lerne, wie Code-Formatierungs- und Wartungs-Features deine Workflows organisiert und lesbar halten.

---

## 6. Code-Formatierung und Wartung

Korrekte Code-Formatierung ist nicht nur für die Ästhetik wichtig, sondern auch zur Verbesserung der Lesbarkeit, des Verständnisses und der Leichtigkeit, komplexe Workflows zu aktualisieren.

### 6.1. Automatische Formatierung in Aktion

Öffne `basic_workflow.nf` und bringe die Formatierung absichtlich durcheinander:

- Entferne einige Einrückungen: Markiere das gesamte Dokument und drücke `shift+tab` sehr oft, um so viele Einrückungen wie möglich zu entfernen.
- Füge zusätzliche Leerzeichen an zufälligen Stellen hinzu: Füge in der `channel.fromPath`-Anweisung 30 Leerzeichen nach der `(` hinzu.
- Breche einige Zeilen ungeschickt: Füge eine neue Zeile zwischen dem `.view {`-Operator und dem `Processing sample:`-String hinzu, aber füge keine entsprechende neue Zeile vor der schließenden Klammer `}` hinzu.

Drücke jetzt `Shift+Alt+F` (oder `Shift+Option+F` auf MacOS), um automatisch zu formatieren:

VS Code behebt sofort:

- Einrückung, um die Prozess-Struktur klar zu zeigen
- Richtet ähnliche Elemente konsistent aus
- Entfernt unnötige Leerzeichen
- Behält lesbare Zeilenumbrüche bei

Beachte, dass die automatische Formatierung möglicherweise nicht jedes Code-Stil-Problem löst. Der Nextflow-Language-Server zielt darauf ab, deinen Code ordentlich zu halten, respektiert aber auch deine persönlichen Präferenzen in bestimmten Bereichen. Wenn du zum Beispiel die Einrückung innerhalb des `script`-Blocks eines Prozesses entfernst, lässt der Formatter sie so, wie sie ist, da du diesen Stil möglicherweise absichtlich bevorzugst.

Derzeit gibt es keine strikte Stil-Durchsetzung für Nextflow, daher bietet der Language-Server etwas Flexibilität. Er wird jedoch konsistent Formatierungsregeln um Methoden- und Funktions-Definitionen herum anwenden, um Klarheit zu erhalten.

### 6.2. Code-Organisations-Features

#### Schnelles Kommentieren

Wähle einen Code-Block in deinem Workflow aus und drücke **Ctrl+/** (oder **Cmd+/**), um ihn auszukommentieren:

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

- Temporäres Deaktivieren von Teilen von Workflows während der Entwicklung
- Hinzufügen erklärender Kommentare zu komplexen Kanal-Operationen
- Dokumentieren von Workflow-Abschnitten

Verwende **Ctrl+/** (oder **Cmd+/**) erneut, um den Code zu entkommentieren.

#### Code-Faltung für Übersicht

In `complex_workflow.nf` beachte die kleinen Pfeile neben Prozess-Definitionen. Klicke darauf, um Prozesse zu falten (zusammenzuklappen):

![Code-Faltung](img/code_folding.png)

Dies gibt dir eine High-Level-Übersicht deiner Workflow-Struktur, ohne dich in Implementierungsdetails zu verlieren.

#### Klammer-Matching

Platziere deinen Cursor neben einer `{`- oder `}`-Klammer und VS Code hebt die passende Klammer hervor. Verwende **Ctrl+Shift+\\** (oder **Cmd+Shift+\\**), um zwischen passenden Klammern zu springen.

Dies ist entscheidend für:

- Verstehen von Prozess-Grenzen
- Finden fehlender oder zusätzlicher Klammern
- Navigieren verschachtelter Workflow-Strukturen

#### Mehrzeilige Auswahl und Bearbeitung

Für die Bearbeitung mehrerer Zeilen gleichzeitig bietet VS Code leistungsstarke Multi-Cursor-Funktionen:

- **Mehrzeilige Auswahl**: Halte **Ctrl+Alt** (oder **Cmd+Option** für MacOS) und verwende Pfeiltasten, um mehrere Zeilen auszuwählen
- **Mehrzeiliges Einrücken**: Wähle mehrere Zeilen aus und verwende **Tab** zum Einrücken oder **Shift+Tab** zum Ausrücken ganzer Blöcke

Dies ist besonders nützlich für:

- Konsistentes Einrücken ganzer Prozess-Blöcke
- Hinzufügen von Kommentaren zu mehreren Zeilen auf einmal
- Bearbeiten ähnlicher Parameter-Definitionen über mehrere Prozesse hinweg

### Fazit

Du kannst sauberen, lesbaren Code pflegen, indem du automatische Formatierung, Kommentar-Features, Code-Faltung, Klammer-Matching und mehrzeilige Bearbeitung verwendest, um komplexe Workflows effizient zu organisieren.

### Wie geht es weiter?

Lerne, wie VS Code sich in deinen breiteren Entwicklungs-Workflow integriert, über das bloße Bearbeiten von Code hinaus.

---

## 7. Entwicklungs-Workflow-Integration

VS Code integriert sich gut in deinen Entwicklungs-Workflow über das bloße Bearbeiten von Code hinaus.

### 7.1. Versionskontroll-Integration

!!! note "Codespaces und Git-Integration"

    Wenn du in **GitHub Codespaces** arbeitest, funktionieren einige Git-Integrations-Features möglicherweise nicht wie erwartet, insbesondere Tastenkombinationen für Source Control. Du hast möglicherweise auch abgelehnt, das Verzeichnis während des initialen Setups als Git-Repository zu öffnen, was für Trainingszwecke in Ordnung ist.

Wenn dein Projekt ein Git-Repository ist (wie dieses), zeigt VS Code:

- Geänderte Dateien mit farbigen Indikatoren
- Git-Status in der Statusleiste
- Inline-Diff-Ansichten
- Commit- und Push-Funktionen

Öffne das Source Control-Panel mit dem Source Control-Button (![Source Control-Symbol](img/source_control_icon.png)) (`Ctrl+Shift+G` oder `Cmd+Shift+G`, wenn du mit VSCode lokal arbeitest), um Git-Änderungen zu sehen und Commits direkt im Editor zu stagen.

![Source Control-Panel](img/source_control.png)

### 7.2. Workflows ausführen und inspizieren

Lass uns einen Workflow ausführen und dann die Ergebnisse inspizieren. Im integrierten Terminal (`Ctrl+Shift+` Backtick in Windows und MacOS), führe den einfachen Workflow aus:

```bash title="Run the basic workflow"
nextflow run basic_workflow.nf --input data/sample_data.csv --output_dir results
```

Während der Workflow läuft, siehst du Echtzeit-Ausgabe im Terminal. Nach Abschluss kannst du VS Code verwenden, um Ergebnisse zu inspizieren, ohne deinen Editor zu verlassen:

1. **Navigiere zu Work-Verzeichnissen**: Verwende den Datei-Explorer oder das Terminal, um `.nextflow/work` zu durchsuchen
2. **Öffne Log-Dateien**: Klicke auf Log-Dateipfade in der Terminal-Ausgabe, um sie direkt in VS Code zu öffnen
3. **Inspiziere Ausgaben**: Durchsuche veröffentlichte Ergebnis-Verzeichnisse im Datei-Explorer
4. **Zeige Ausführungs-Reports an**: Öffne HTML-Reports direkt in VS Code oder deinem Browser

Dies hält alles an einem Ort, anstatt zwischen mehreren Anwendungen zu wechseln.

### Fazit

Du kannst VS Code mit Versionskontrolle und Workflow-Ausführung integrieren, um deinen gesamten Entwicklungsprozess von einer einzigen Oberfläche aus zu verwalten.

### Wie geht es weiter?

Sieh, wie all diese IDE-Features in deinem täglichen Entwicklungs-Workflow zusammenarbeiten.

---

## 8. Zusammenfassung und schnelle Notizen

Hier sind einige schnelle Notizen zu jedem der oben besprochenen IDE-Features:

### 8.1. Ein neues Feature starten

1. **Schnelles Datei-Öffnen** (`Ctrl+P` oder `Cmd+P`), um relevante existierende Module zu finden
2. **Split-Editor**, um ähnliche Prozesse nebeneinander anzuzeigen
3. **Symbol-Navigation** (`Ctrl+Shift+O` oder `Cmd+Shift+O`), um die Datei-Struktur zu verstehen
4. **Auto-Vervollständigung**, um neuen Code schnell zu schreiben

### 8.2. Probleme debuggen

1. **Problems-Panel** (`Ctrl+Shift+M` oder `Cmd+Shift+M`), um alle Fehler auf einmal zu sehen
2. **Gehe zu Definition** (`Ctrl-Klick` oder `Cmd-Klick`), um Prozess-Interfaces zu verstehen
3. **Alle Referenzen finden**, um zu sehen, wie Prozesse verwendet werden
4. **Projektweite Suche**, um ähnliche Muster oder Probleme zu finden

### 8.3. Refactoring und Verbesserung

1. **Projektweite Suche** (`Ctrl+Shift+F` oder `Cmd+Shift+F`), um Muster zu finden
2. **Auto-Formatierung** (`Shift+Alt+F` oder `Shift+Option+F`), um Konsistenz zu erhalten
3. **Code-Faltung**, um sich auf die Struktur zu konzentrieren
4. **Git-Integration**, um Änderungen zu verfolgen

---

## Zusammenfassung

Du hast jetzt eine Schnelltour durch VS Codes IDE-Features für Nextflow-Entwicklung gemacht. Diese Tools werden dich deutlich produktiver machen, indem sie:

- **Fehler reduzieren** durch Echtzeit-Syntax-Überprüfung
- **Entwicklung beschleunigen** mit intelligenter Auto-Vervollständigung
- **Navigation verbessern** in komplexen Multi-File-Workflows
- **Qualität erhalten** durch konsistente Formatierung
- **Verständnis verbessern** durch erweitertes Highlighting und Struktur-Visualisierung

Wir erwarten nicht, dass du dich an alles erinnerst, aber jetzt weißt du, dass diese Features existieren, und du wirst sie finden können, wenn du sie brauchst. Während du weiterhin Nextflow-Workflows entwickelst, werden diese IDE-Features zur zweiten Natur, sodass du dich auf das Schreiben von hochwertigem Code konzentrieren kannst, anstatt mit Syntax und Struktur zu kämpfen.

### Wie geht es weiter?

Wende diese IDE-Fähigkeiten an, während du andere Trainingsmodule durcharbeitest, zum Beispiel:

- **[nf-test](nf-test.md)**: Erstelle umfassende Test-Suites für deine Workflows
- **[Hello nf-core](../../hello_nf-core/)**: Baue produktionsreife Pipelines mit Community-Standards

Die wahre Kraft dieser IDE-Features zeigt sich, wenn du an größeren, komplexeren Projekten arbeitest. Beginne damit, sie schrittweise in deinen Workflow zu integrieren - innerhalb weniger Sessions werden sie zur zweiten Natur und transformieren, wie du Nextflow-Entwicklung angehst.

Von der Fehlererkennung, bevor sie dich verlangsamen, bis zur Navigation durch komplexe Codebasen mit Leichtigkeit - diese Tools werden dich zu einem selbstbewussteren und effizienteren Entwickler\*in machen.

Viel Spaß beim Coden!
