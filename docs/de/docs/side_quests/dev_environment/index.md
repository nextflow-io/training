# Entwicklungsumgebung

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } KI-gestützte Übersetzung - [mehr erfahren & Verbesserungen vorschlagen](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Moderne Integrated Development Environments (IDEs) können deine Nextflow-Entwicklung grundlegend verändern. Dieser Side Quest konzentriert sich darauf, wie du VS Code und seine Nextflow-Erweiterung nutzen kannst, um schneller Code zu schreiben, Fehler frühzeitig zu erkennen und komplexe Workflows effizient zu navigieren.

!!! note "Das ist kein klassisches Tutorial"

    Anders als andere Trainingsmodule ist dieser Leitfaden als Sammlung von kurzen Hinweisen, Tipps und praktischen Beispielen aufgebaut – kein Schritt-für-Schritt-Tutorial. Jeder Abschnitt kann unabhängig erkundet werden, je nach deinen Interessen und aktuellen Entwicklungsbedürfnissen. Spring ruhig hin und her und konzentriere dich auf die Funktionen, die für deine Workflow-Entwicklung am nützlichsten sind.

## Was du vorher wissen solltest

Dieser Leitfaden setzt voraus, dass du den [Hello Nextflow](../hello_nextflow/)-Kurs abgeschlossen hast und mit den grundlegenden Nextflow-Konzepten vertraut bist, darunter:

- **Grundlegende Workflow-Struktur**: Prozesse, Workflows und wie sie miteinander verbunden sind
- **Kanal-Operationen**: Kanäle erstellen, Daten zwischen Prozessen weitergeben und grundlegende Operatoren verwenden
- **Module und Organisation**: Wiederverwendbare Module erstellen und `include`-Anweisungen verwenden
- **Konfigurationsgrundlagen**: `nextflow.config` für Parameter, Prozess-Direktiven und Profile verwenden

## Was du hier lernst

Dieser Leitfaden konzentriert sich auf **IDE-Produktivitätsfunktionen**, die dich zu einer effizienteren Nextflow-Entwickler\*in machen:

- **Erweitertes Syntax-Highlighting**: Verstehen, was VS Code dir über deine Code-Struktur zeigt
- **Intelligente Auto-Vervollständigung**: Kontextbezogene Vorschläge für schnelleres Schreiben nutzen
- **Fehlererkennung und Diagnose**: Syntaxfehler erkennen, bevor du deinen Workflow ausführst
- **Code-Navigation**: Schnell zwischen Prozessen, Modulen und Definitionen wechseln
- **Formatierung und Organisation**: Konsistenten, lesbaren Code-Stil beibehalten
- **KI-gestützte Entwicklung** (optional): Moderne KI-Tools nutzen, die in deine IDE integriert sind

!!! info "Warum jetzt IDE-Funktionen?"

    Du hast VS Code wahrscheinlich schon während des [Hello Nextflow](../hello_nextflow/)-Kurses verwendet, aber wir haben uns dort auf die Nextflow-Grundlagen konzentriert statt auf IDE-Funktionen. Jetzt, wo du mit grundlegenden Nextflow-Konzepten wie Prozessen, Workflows, Kanälen und Modulen vertraut bist, kannst du die leistungsstarken IDE-Funktionen nutzen, die dich effizienter machen.

    Stell dir das als „Level-Up" für deine Entwicklungsumgebung vor – der Editor, den du schon verwendest, hat viel mächtigere Fähigkeiten, die wirklich wertvoll werden, sobald du verstehst, wobei sie dir helfen.

---

## 0. Einrichtung und Aufwärmen

Richten wir einen Arbeitsbereich speziell zum Erkunden der IDE-Funktionen ein:

```bash title="Navigate to the IDE features directory"
cd side-quests/ide_features
```

Öffne dieses Verzeichnis in VS Code:

```bash title="Open VS Code in current directory"
code .
```

Das Verzeichnis `ide_features` enthält Beispiel-Workflows, die verschiedene IDE-Funktionen demonstrieren:

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

    - `basic_workflow.nf` ist ein funktionierender einfacher Workflow, den du ausführen und anpassen kannst
    - `complex_workflow.nf` dient nur zur Veranschaulichung von Navigationsfunktionen – er läuft möglicherweise nicht erfolgreich, zeigt aber eine realistische Workflow-Struktur mit mehreren Dateien

### Tastaturkürzel

Einige der Funktionen in diesem Leitfaden verwenden optionale Tastaturkürzel. Wenn du über GitHub Codespaces im Browser arbeitest, funktionieren manche Kürzel möglicherweise nicht wie erwartet, weil sie dort für andere Dinge verwendet werden.

Wenn du VS Code lokal ausführst – was wahrscheinlich der Fall sein wird, wenn du tatsächlich Workflows schreibst – funktionieren die Kürzel wie beschrieben.

Auf einem Mac verwenden manche (nicht alle) Tastaturkürzel „cmd" statt „ctrl". Wir weisen darauf im Text hin, z. B. mit `Ctrl/Cmd`.

### 0.1. Die Nextflow-Erweiterung installieren

!!! note "Verwendest du bereits Devcontainer?"

    Wenn du in **GitHub Codespaces** arbeitest oder einen **lokalen Devcontainer** verwendest, ist die Nextflow-Erweiterung wahrscheinlich bereits installiert und konfiguriert. Du kannst die manuellen Installationsschritte überspringen und direkt mit den Erweiterungsfunktionen loslegen.

So installierst du die Erweiterung manuell:

1. Öffne VS Code
2. Öffne die Erweiterungsansicht, indem du auf das Erweiterungssymbol links klickst: ![Erweiterungssymbol](img/extensions_icon.png) (Kürzel `Ctrl/Cmd+Shift+X`, wenn du VS Code lokal ausführst)
3. Suche nach „Nextflow"
4. Installiere die offizielle Nextflow-Erweiterung

![Nextflow-Erweiterung installieren](img/install_extension.png)

### 0.2. Arbeitsbereich-Layout

Da du VS Code bereits während Hello Nextflow verwendet hast, kennst du die Grundlagen. So organisierst du deinen Arbeitsbereich effizient für diese Sitzung:

- **Editor-Bereich**: Zum Anzeigen und Bearbeiten von Dateien. Du kannst ihn in mehrere Bereiche aufteilen, um Dateien nebeneinander zu vergleichen.
- **Datei-Explorer** (![Datei-Explorer-Symbol](img/files_icon.png)) (`Ctrl/Cmd+Shift+E`): Die lokalen Dateien und Ordner auf deinem System. Halte ihn links offen, um zwischen Dateien zu navigieren.
- **Integriertes Terminal** (`Ctrl+Shift+` Backtick für Windows und MacOS): Ein Terminal am unteren Rand für die Interaktion mit dem Computer. Verwende es, um Nextflow oder andere Befehle auszuführen.
- **Probleme-Panel** (`Ctrl+Shift+M`): VS Code zeigt hier alle erkannten Fehler und Probleme an. Nützlich, um Probleme auf einen Blick zu sehen.

Du kannst Panels verschieben oder ausblenden (`Ctrl/Cmd+B` zum Ein-/Ausblenden der Seitenleiste), um dein Layout anzupassen, während wir die Beispiele durcharbeiten.

### Fazit

Du hast VS Code mit der Nextflow-Erweiterung eingerichtet und verstehst das Arbeitsbereich-Layout für effiziente Entwicklung.

### Wie geht es weiter?

Lerne, wie Syntax-Highlighting dir hilft, die Nextflow-Code-Struktur auf einen Blick zu verstehen.

---

## 1. Syntax-Highlighting und Code-Struktur

Jetzt, wo dein Arbeitsbereich eingerichtet ist, erkunden wir, wie VS Code's Syntax-Highlighting dir hilft, Nextflow-Code effektiver zu lesen und zu schreiben.

### 1.1. Nextflow-Syntaxelemente

Öffne `basic_workflow.nf`, um Syntax-Highlighting in Aktion zu sehen:

![Syntax-Übersicht](img/syntax_showcase.png)

Beachte, wie VS Code hervorhebt:

- **Schlüsselwörter** (`process`, `workflow`, `input`, `output`, `script`) in verschiedenen Farben
- **String-Literale** und **Parameter** mit unterschiedlicher Formatierung
- **Kommentare** in einer gedämpften Farbe
- **Variablen** und **Funktionsaufrufe** mit passender Betonung
- **Code-Blöcke** mit korrekten Einrückungslinien

!!! note "Theme-abhängige Farben"

    Die genauen Farben hängen von deinem VS Code-Theme (Hell-/Dunkelmodus), Farbeinstellungen und Anpassungen ab. Wichtig ist, dass verschiedene Syntaxelemente visuell voneinander unterschieden werden, sodass die Code-Struktur unabhängig vom gewählten Farbschema leichter zu verstehen ist.

### 1.2. Code-Struktur verstehen

Das Syntax-Highlighting hilft dir, schnell zu erkennen:

- **Prozess-Grenzen**: Klare Unterscheidung zwischen verschiedenen Prozessen
- **Eingabe-/Ausgabe-Blöcke**: Datenflussdefinitionen leicht erkennbar
- **Script-Blöcke**: Die tatsächlich ausgeführten Befehle
- **Kanal-Operationen**: Datentransformationsschritte
- **Konfigurations-Direktiven**: Prozessspezifische Einstellungen

Diese visuelle Organisation wird unverzichtbar, wenn du mit komplexen Workflows arbeitest, die mehrere Prozesse und komplizierte Datenflüsse enthalten.

### Fazit

Du verstehst, wie VS Code's Syntax-Highlighting dir hilft, die Nextflow-Code-Struktur zu lesen und verschiedene Sprachelemente für schnellere Entwicklung zu erkennen.

### Wie geht es weiter?

Lerne, wie intelligente Auto-Vervollständigung das Code-Schreiben mit kontextbezogenen Vorschlägen beschleunigt.

---

## 2. Intelligente Auto-Vervollständigung

VS Code's Auto-Vervollständigung hilft dir, Code schneller und mit weniger Fehlern zu schreiben, indem passende Optionen basierend auf dem Kontext vorgeschlagen werden.

### 2.1. Kontextbezogene Vorschläge

Die Auto-Vervollständigungsoptionen variieren je nachdem, wo du dich im Code befindest:

#### Kanal-Operationen

Öffne `basic_workflow.nf` erneut und tippe `channel.` im Workflow-Block:

![Kanal-Auto-Vervollständigung](img/autocomplete_channel.png)

Du siehst Vorschläge für:

- `fromPath()` – Kanal aus Dateipfaden erstellen
- `fromFilePairs()` – Kanal aus Dateipaaren erstellen
- `of()` – Kanal aus Werten erstellen
- `fromSRA()` – Kanal aus SRA-Accessions erstellen
- Und viele mehr...

Das hilft dir, schnell die richtige Kanal-Factory zu finden, ohne dir genaue Methodennamen merken zu müssen.

Du kannst auch die verfügbaren Operatoren für Kanäle entdecken. Tippe zum Beispiel `FASTQC.out.html.`, um verfügbare Operationen zu sehen:

![Kanal-Operationen Auto-Vervollständigung](img/autocomplete_operators.png)

#### Prozess-Direktiven

Tippe innerhalb eines Prozess-Script-Blocks `task.`, um verfügbare Laufzeiteigenschaften zu sehen:

![Task-Eigenschaften Auto-Vervollständigung](img/autocomplete_task.png)

#### Konfiguration

Öffne `nextflow.config` und tippe `process.` an einer beliebigen Stelle, um verfügbare Prozess-Direktiven zu sehen:

![Konfigurations-Auto-Vervollständigung](img/autocomplete_config.png)

Du siehst Vorschläge für:

- `executor`
- `memory`
- `cpus`

Das spart Zeit beim Konfigurieren von Prozessen und funktioniert über verschiedene Konfigurationsbereiche hinweg. Tippe zum Beispiel `docker.`, um Docker-spezifische Konfigurationsoptionen zu sehen.

### Fazit

Du kannst VS Code's intelligente Auto-Vervollständigung nutzen, um verfügbare Kanal-Operationen, Prozess-Direktiven und Konfigurationsoptionen zu entdecken, ohne Syntax auswendig lernen zu müssen.

### Wie geht es weiter?

Lerne, wie Echtzeit-Fehlererkennung dir hilft, Probleme zu erkennen, bevor du deinen Workflow ausführst – einfach durch Lesen des Codes.

## 3. Fehlererkennung und Diagnose

VS Code's Echtzeit-Fehlererkennung hilft dir, Probleme zu erkennen, bevor du deinen Workflow ausführst.

### 3.1. Syntaxfehler-Erkennung

Erstellen wir einen absichtlichen Fehler, um die Erkennung in Aktion zu sehen. Öffne `basic_workflow.nf` und ändere den Prozessnamen von `FASTQC` zu `FASTQ` (oder einem anderen ungültigen Namen). VS Code hebt den Fehler im Workflow-Block sofort mit einer roten Wellenlinie hervor:

![Fehler-Unterstreichung](img/error_underline.png)

### 3.2. Probleme-Panel

Über die individuelle Fehlerhervorhebung hinaus bietet VS Code ein zentrales Probleme-Panel, das alle Fehler, Warnungen und Informationsmeldungen in deinem Arbeitsbereich zusammenfasst. Öffne es mit `Ctrl/Cmd+Shift+M` und verwende das Filter-Symbol, um nur Fehler der aktuellen Datei anzuzeigen:

![Probleme-Panel filtern](img/active_file.png)

Klicke auf ein Problem, um direkt zur betroffenen Zeile zu springen:

![Probleme-Panel](img/problems_panel.png)

Behebe den Fehler, indem du den Prozessnamen wieder auf `FASTQC` änderst.

### 3.3. Häufige Fehlermuster

Häufige Fehler in der Nextflow-Syntax sind:

- **Fehlende Klammern**: Nicht übereinstimmende `{` oder `}`
- **Unvollständige Blöcke**: Fehlende Pflichtabschnitte in Prozessen
- **Ungültige Syntax**: Fehlerhaftes Nextflow DSL
- **Tippfehler in Schlüsselwörtern**: Falsch geschriebene Prozess-Direktiven
- **Kanal-Inkompatibilitäten**: Typinkompatibilitäten

Der Nextflow Language Server hebt diese Probleme im Probleme-Panel hervor. Du kannst sie frühzeitig prüfen, um Syntaxfehler beim Ausführen einer Pipeline zu vermeiden.

### Fazit

Du kannst VS Code's Fehlererkennung und das Probleme-Panel nutzen, um Syntaxfehler und Probleme zu erkennen, bevor du deinen Workflow ausführst – das spart Zeit und verhindert Frust.

### Wie geht es weiter?

Lerne, wie du effizient zwischen Prozessen, Modulen und Definitionen in komplexen Workflows navigierst.

---

## 4. Code-Navigation und Symbol-Verwaltung

Effiziente Navigation ist entscheidend, wenn du mit komplexen Workflows arbeitest, die sich über mehrere Dateien erstrecken. Um das zu verstehen, ersetze die Prozessdefinition in `basic_workflow.nf` durch einen Import des bereitgestellten Moduls:

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

### 4.1. Zur Definition springen

Wenn du mit der Maus über einen Prozessnamen wie `FASTQC` fährst, siehst du ein Popup mit der Modul-Schnittstelle (Eingaben und Ausgaben):

![Zur Definition springen](img/syntax.png)

Diese Funktion ist besonders wertvoll beim Schreiben von Workflows, da du die Modul-Schnittstelle verstehen kannst, ohne die Moduldatei direkt öffnen zu müssen.

Du kannst schnell zu jeder Prozess-, Modul- oder Variablendefinition navigieren, indem du **Ctrl/Cmd-klickst**. Fahre mit der Maus über den Link zur Moduldatei am Anfang des Skripts und folge dem Link wie vorgeschlagen:

![Link folgen](img/follow_link.png)

Das Gleiche funktioniert für Prozessnamen. Gehe zurück zu `basic_workflow.nf` und probiere es beim `FASTQC`-Prozessnamen im Workflow-Block aus. Das führt dich direkt zum Prozessnamen (der in diesem Beispiel identisch mit der Moduldatei ist, aber auch mitten in einer viel größeren Datei liegen könnte).

Um zurückzukehren, verwende **Alt+←** (oder **Ctrl+-** auf dem Mac). Das ist eine leistungsstarke Möglichkeit, Code zu erkunden, ohne deinen Platz zu verlieren.

Erkunden wir nun die Navigation in einem komplexeren Workflow mit `complex_workflow.nf` (die zuvor erwähnte Datei nur zur Veranschaulichung). Dieser Workflow enthält mehrere Prozesse, die in separaten Moduldateien definiert sind, sowie einige inline definierte. Während komplexe Mehrfachdatei-Strukturen manuell schwer zu navigieren sein können, macht die Möglichkeit, zu Definitionen zu springen, die Erkundung viel überschaubarer.

1. Öffne `complex_workflow.nf`
2. Navigiere zu Moduldefinitionen
3. Verwende **Alt+←** (oder **Ctrl+-**), um zurückzunavigieren
4. Navigiere zum `FASTQC`-Prozessnamen im Workflow-Block. Das führt dich direkt zum Prozessnamen (der in diesem Beispiel identisch mit der Moduldatei ist, aber auch mitten in einer viel größeren Datei liegen könnte).
5. Navigiere wieder zurück
6. Navigiere zum `TRIM_GALORE`-Prozess im Workflow-Block. Dieser ist inline definiert, führt dich also nicht zu einer separaten Datei, zeigt dir aber trotzdem die Prozessdefinition – und du kannst trotzdem zurücknavigieren.

### 4.2. Symbol-Navigation

Mit noch geöffnetem `complex_workflow.nf` kannst du eine Übersicht aller Symbole in der Datei erhalten, indem du `@` in die Suchleiste oben in VS Code eingibst (das Tastaturkürzel ist `Ctrl/Cmd+Shift+O`, funktioniert aber möglicherweise nicht in Codespaces). Das öffnet das Symbol-Navigationspanel, das alle Symbole in der aktuellen Datei auflistet:

![Symbol-Navigation](img/symbols.png)

Das zeigt:

- Alle Prozessdefinitionen
- Workflow-Definitionen (in dieser Datei sind zwei Workflows definiert)
- Funktionsdefinitionen

Tippe, um die Ergebnisse zu filtern.

### 4.3. Alle Referenzen finden

Es kann sehr hilfreich sein zu verstehen, wo ein Prozess oder eine Variable in deiner Codebasis verwendet wird. Wenn du zum Beispiel alle Referenzen auf den `FASTQC`-Prozess finden möchtest, navigiere zunächst zu seiner Definition. Du kannst das tun, indem du `modules/fastqc.nf` direkt öffnest oder VS Code's Schnellnavigation mit `Ctrl/Cmd-Klick` verwendest, wie oben beschrieben. Sobald du bei der Prozessdefinition bist, klicke mit der rechten Maustaste auf den `FASTQC`-Prozessnamen und wähle „Find All References" aus dem Kontextmenü, um alle Stellen zu sehen, wo er verwendet wird.

![Referenzen finden](img/references.png)

Diese Funktion zeigt alle Stellen, wo `FASTQC` in deinem Arbeitsbereich referenziert wird, einschließlich seiner Verwendung in den zwei verschiedenen Workflows. Das ist entscheidend, um die möglichen Auswirkungen von Änderungen am `FASTQC`-Prozess einzuschätzen.

### 4.4. Outline-Panel

Das Outline-Panel in der Explorer-Seitenleiste (klicke auf ![Explorer-Symbol](img/files_icon.png)) bietet eine praktische Übersicht aller Symbole in der aktuellen Datei. Es ermöglicht dir, schnell durch die Struktur deines Codes zu navigieren, indem Funktionen, Variablen und andere wichtige Elemente in einer hierarchischen Ansicht angezeigt werden.

![Outline-Panel](img/outline.png)

Verwende das Outline-Panel, um schnell zu verschiedenen Teilen deines Codes zu navigieren, ohne den Datei-Browser zu verwenden.

### 4.5. DAG-Visualisierung

VS Code's Nextflow-Erweiterung kann deinen Workflow als Directed Acyclic Graph (DAG) visualisieren. Das hilft dir, den Datenfluss und die Abhängigkeiten zwischen Prozessen zu verstehen. Öffne `complex_workflow.nf` und klicke auf die Schaltfläche „Preview DAG" über `workflow {` (dem zweiten `workflow`-Block in dieser Datei):

![DAG-Vorschau](img/dag_preview.png)

Das ist nur der „Entry"-Workflow, aber du kannst auch den DAG für die inneren Workflows anzeigen, indem du auf die Schaltfläche „Preview DAG" über dem Workflow `RNASEQ_PIPELINE {` weiter oben klickst:

![DAG-Vorschau innerer Workflow](img/dag_preview_inner.png)

Bei diesem Workflow kannst du die Knoten im DAG verwenden, um zu den entsprechenden Prozessdefinitionen im Code zu navigieren. Klicke auf einen Knoten, und er führt dich zur relevanten Prozessdefinition im Editor. Besonders wenn ein Workflow sehr groß wird, hilft das wirklich dabei, im Code zu navigieren und zu verstehen, wie die Prozesse miteinander verbunden sind.

### Fazit

Du kannst komplexe Workflows effizient navigieren, indem du „Zur Definition springen", Symbol-Suche, „Alle Referenzen finden" und DAG-Visualisierung verwendest, um Code-Struktur und Abhängigkeiten zu verstehen.

### Wie geht es weiter?

Lerne, wie du effektiv mit mehreren miteinander verbundenen Dateien in größeren Nextflow-Projekten arbeitest.

## 5. Arbeiten mit mehreren Dateien

Echte Nextflow-Entwicklung bedeutet, mit mehreren miteinander verbundenen Dateien zu arbeiten. Erkunden wir, wie VS Code dir hilft, komplexe Projekte effizient zu verwalten.

### 5.1. Schnelle Datei-Navigation

Mit geöffnetem `complex_workflow.nf` siehst du, dass es mehrere Module importiert. Üben wir die schnelle Navigation zwischen ihnen.

Drücke **Ctrl+P** (oder **Cmd+P**) und tippe „fast":

VS Code zeigt dir passende Dateien. Wähle `modules/fastqc.nf`, um sofort dorthin zu springen. Das ist viel schneller als durch den Datei-Explorer zu klicken, wenn du ungefähr weißt, welche Datei du suchst.

Probiere das mit anderen Mustern:

- Tippe „star", um die STAR-Alignment-Moduldatei zu finden (`star.nf`)
- Tippe „utils", um die Datei mit Hilfsfunktionen zu finden (`utils.nf`)
- Tippe „config", um zu Konfigurationsdateien zu springen (`nextflow.config`)

### 5.2. Geteilter Editor für Mehrfachdatei-Entwicklung

Wenn du mit Modulen arbeitest, musst du oft sowohl den Haupt-Workflow als auch Moduldefinitionen gleichzeitig sehen. So richtest du das ein:

1. Öffne `complex_workflow.nf`
2. Öffne `modules/fastqc.nf` in einem neuen Tab
3. Klicke mit der rechten Maustaste auf den Tab `modules/fastqc.nf` und wähle „Split Right"
4. Jetzt siehst du beide Dateien nebeneinander

![Geteilter Editor](img/split_editor.png)

Das ist unverzichtbar, wenn du:

- Modul-Schnittstellen prüfst, während du Workflow-Aufrufe schreibst, und die Vorschau nicht ausreicht
- Ähnliche Prozesse über verschiedene Module vergleichst
- Den Datenfluss zwischen Workflow und Modulen debuggst

### 5.3. Projektweite Suche

Manchmal musst du herausfinden, wo bestimmte Muster in deinem gesamten Projekt verwendet werden. Drücke `Ctrl/Cmd+Shift+F`, um das Suchpanel zu öffnen.

Suche nach `publishDir` im gesamten Arbeitsbereich:

![Projektsuche](img/project_search.png)

Das zeigt dir jede Datei, die Ausgabeverzeichnisse verwendet, und hilft dir:

- Ausgabe-Organisationsmuster zu verstehen
- Beispiele für bestimmte Direktiven zu finden
- Konsistenz über Module hinweg sicherzustellen

### Fazit

Du kannst komplexe Mehrfachdatei-Projekte mit schneller Datei-Navigation, geteilten Editoren und projektweiter Suche verwalten, um effizient über Workflows und Module hinweg zu arbeiten.

### Wie geht es weiter?

Lerne, wie Code-Formatierung und Wartungsfunktionen deine Workflows organisiert und lesbar halten.

---

## 6. Code-Formatierung und Wartung

Korrekte Code-Formatierung ist nicht nur für die Ästhetik wichtig, sondern verbessert auch die Lesbarkeit, das Verständnis und die einfache Aktualisierung komplexer Workflows.

### 6.1. Automatische Formatierung in Aktion

Öffne `basic_workflow.nf` und mache die Formatierung absichtlich kaputt:

- Entferne einige Einrückungen: Markiere das gesamte Dokument und drücke `Shift+Tab` viele Male, um so viele Einrückungen wie möglich zu entfernen.
- Füge an zufälligen Stellen zusätzliche Leerzeichen ein: Füge in der `channel.fromPath`-Anweisung 30 Leerzeichen nach der `(` ein.
- Breche einige Zeilen ungeschickt um: Füge eine neue Zeile zwischen dem `.view {`-Operator und dem `Processing sample:`-String ein, aber füge keine entsprechende neue Zeile vor der schließenden Klammer `}` hinzu.

Drücke jetzt `Shift+Alt+F` (oder `Shift+Option+F` auf MacOS) zur automatischen Formatierung:

VS Code:

- Korrigiert die Einrückung, um die Prozessstruktur klar darzustellen
- Richtet ähnliche Elemente konsistent aus
- Entfernt unnötige Leerzeichen
- Behält lesbare Zeilenumbrüche bei

Beachte, dass die automatische Formatierung nicht alle Code-Stil-Probleme löst. Der Nextflow Language Server versucht, deinen Code ordentlich zu halten, respektiert aber in bestimmten Bereichen auch deine persönlichen Vorlieben. Wenn du zum Beispiel die Einrückung innerhalb des `script`-Blocks eines Prozesses entfernst, lässt der Formatter das so – vielleicht bevorzugst du diesen Stil absichtlich.

Derzeit gibt es keine strenge Stil-Durchsetzung für Nextflow, daher bietet der Language Server etwas Flexibilität. Er wendet jedoch konsistent Formatierungsregeln rund um Methoden- und Funktionsdefinitionen an, um die Übersichtlichkeit zu wahren.

### 6.2. Code-Organisationsfunktionen

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

Das ist perfekt für:

- Vorübergehendes Deaktivieren von Workflow-Teilen während der Entwicklung
- Hinzufügen erklärender Kommentare zu komplexen Kanal-Operationen
- Dokumentieren von Workflow-Abschnitten

Drücke **Ctrl+/** (oder **Cmd+/**) erneut, um den Code wieder einzukommentieren.

#### Code-Faltung für den Überblick

In `complex_workflow.nf` siehst du kleine Pfeile neben Prozessdefinitionen. Klicke darauf, um Prozesse zu falten (einzuklappen):

![Code-Faltung](img/code_folding.png)

Das gibt dir einen Überblick über deine Workflow-Struktur, ohne dich in Implementierungsdetails zu verlieren.

#### Klammer-Matching

Platziere deinen Cursor neben einer `{`- oder `}`-Klammer, und VS Code hebt die passende Klammer hervor. Verwende **Ctrl+Shift+\\** (oder **Cmd+Shift+\\**), um zwischen passenden Klammern zu springen.

Das ist entscheidend für:

- Prozess-Grenzen zu verstehen
- Fehlende oder überschüssige Klammern zu finden
- Verschachtelte Workflow-Strukturen zu navigieren

#### Mehrzeilige Auswahl und Bearbeitung

Für die gleichzeitige Bearbeitung mehrerer Zeilen bietet VS Code leistungsstarke Mehrcursor-Funktionen:

- **Mehrzeilige Auswahl**: Halte **Ctrl+Alt** (oder **Cmd+Option** auf MacOS) gedrückt und verwende Pfeiltasten, um mehrere Zeilen auszuwählen
- **Mehrzeiliges Einrücken**: Wähle mehrere Zeilen aus und verwende **Tab** zum Einrücken oder **Shift+Tab** zum Ausrücken ganzer Blöcke

Das ist besonders nützlich für:

- Konsistentes Einrücken ganzer Prozessblöcke
- Gleichzeitiges Hinzufügen von Kommentaren zu mehreren Zeilen
- Bearbeiten ähnlicher Parameter-Definitionen über mehrere Prozesse hinweg

### Fazit

Du kannst sauberen, lesbaren Code mit automatischer Formatierung, Kommentarfunktionen, Code-Faltung, Klammer-Matching und Mehrzeilenbearbeitung pflegen, um komplexe Workflows effizient zu organisieren.

### Wie geht es weiter?

Lerne, wie VS Code sich in deinen breiteren Entwicklungs-Workflow integriert – über das reine Code-Bearbeiten hinaus.

---

## 7. Integration in den Entwicklungs-Workflow

VS Code lässt sich gut in deinen Entwicklungs-Workflow integrieren – über das reine Code-Bearbeiten hinaus.

### 7.1. Versionskontrolle

!!! note "Codespaces und Git-Integration"

    Wenn du in **GitHub Codespaces** arbeitest, funktionieren manche Git-Integrationsfunktionen möglicherweise nicht wie erwartet, insbesondere Tastaturkürzel für die Quellcodeverwaltung. Du hast das Verzeichnis beim ersten Einrichten möglicherweise auch nicht als Git-Repository geöffnet – das ist für Trainingszwecke in Ordnung.

Wenn dein Projekt ein Git-Repository ist (wie dieses), zeigt VS Code:

- Geänderte Dateien mit farbigen Indikatoren
- Git-Status in der Statusleiste
- Inline-Diff-Ansichten
- Commit- und Push-Funktionen

Öffne das Quellcodeverwaltungs-Panel über die Schaltfläche für Quellcodeverwaltung (![Quellcodeverwaltungs-Symbol](img/source_control_icon.png)) (`Ctrl+Shift+G` oder `Cmd+Shift+G`, wenn du VS Code lokal verwendest), um Git-Änderungen zu sehen und Commits direkt im Editor zu erstellen.

![Quellcodeverwaltungs-Panel](img/source_control.png)

### 7.2. Workflows ausführen und Ergebnisse prüfen

Führen wir einen Workflow aus und prüfen dann die Ergebnisse. Führe im integrierten Terminal (`Ctrl+Shift+` Backtick in Windows und MacOS) den einfachen Workflow aus:

```bash title="Run the basic workflow"
nextflow run basic_workflow.nf --input data/sample_data.csv --output_dir results
```

Während der Workflow läuft, siehst du die Ausgabe in Echtzeit im Terminal. Nach Abschluss kannst du VS Code verwenden, um Ergebnisse zu prüfen, ohne den Editor zu verlassen:

1. **Zu Work-Verzeichnissen navigieren**: Verwende den Datei-Explorer oder das Terminal, um `.nextflow/work` zu durchsuchen
2. **Log-Dateien öffnen**: Klicke auf Log-Dateipfade in der Terminal-Ausgabe, um sie direkt in VS Code zu öffnen
3. **Ausgaben prüfen**: Durchsuche veröffentlichte Ergebnisverzeichnisse im Datei-Explorer
4. **Ausführungsberichte anzeigen**: Öffne HTML-Berichte direkt in VS Code oder deinem Browser

Das hält alles an einem Ort, anstatt zwischen mehreren Anwendungen zu wechseln.

### Fazit

Du kannst VS Code mit Versionskontrolle und Workflow-Ausführung integrieren, um deinen gesamten Entwicklungsprozess über eine einzige Oberfläche zu verwalten.

### Wie geht es weiter?

Sieh, wie all diese IDE-Funktionen in deinem täglichen Entwicklungs-Workflow zusammenwirken.

---

## 8. Zusammenfassung und kurze Hinweise

Hier sind einige kurze Hinweise zu den oben besprochenen IDE-Funktionen:

### 8.1. Ein neues Feature beginnen

1. **Schnelles Öffnen von Dateien** (`Ctrl+P` oder `Cmd+P`), um relevante vorhandene Module zu finden
2. **Geteilter Editor**, um ähnliche Prozesse nebeneinander zu sehen
3. **Symbol-Navigation** (`Ctrl+Shift+O` oder `Cmd+Shift+O`), um die Dateistruktur zu verstehen
4. **Auto-Vervollständigung**, um neuen Code schnell zu schreiben

### 8.2. Probleme debuggen

1. **Probleme-Panel** (`Ctrl+Shift+M` oder `Cmd+Shift+M`), um alle Fehler auf einmal zu sehen
2. **Zur Definition springen** (`Ctrl-Klick` oder `Cmd-Klick`), um Prozess-Schnittstellen zu verstehen
3. **Alle Referenzen finden**, um zu sehen, wie Prozesse verwendet werden
4. **Projektweite Suche**, um ähnliche Muster oder Probleme zu finden

### 8.3. Refactoring und Verbesserung

1. **Projektweite Suche** (`Ctrl+Shift+F` oder `Cmd+Shift+F`), um Muster zu finden
2. **Automatische Formatierung** (`Shift+Alt+F` oder `Shift+Option+F`), um Konsistenz zu wahren
3. **Code-Faltung**, um sich auf die Struktur zu konzentrieren
4. **Git-Integration**, um Änderungen zu verfolgen

---

## Zusammenfassung

Du hast jetzt einen Schnelldurchlauf durch VS Code's IDE-Funktionen für die Nextflow-Entwicklung gemacht. Diese Tools machen dich deutlich produktiver:

- **Fehler reduzieren** durch Echtzeit-Syntaxprüfung
- **Entwicklung beschleunigen** mit intelligenter Auto-Vervollständigung
- **Navigation verbessern** in komplexen Mehrfachdatei-Workflows
- **Qualität sichern** durch konsistente Formatierung
- **Verständnis fördern** durch erweitertes Highlighting und Strukturvisualisierung

Wir erwarten nicht, dass du dir alles merkst, aber jetzt weißt du, dass diese Funktionen existieren, und kannst sie finden, wenn du sie brauchst. Während du weiter Nextflow-Workflows entwickelst, werden diese IDE-Funktionen zur zweiten Natur – und du kannst dich auf das Schreiben von hochwertigem Code konzentrieren, anstatt mit Syntax und Struktur zu kämpfen.

### Wie geht es weiter?

Wende diese IDE-Kenntnisse an, während du andere Trainingsmodule durcharbeitest, zum Beispiel:

- **[nf-test](nf-test.md)**: Erstelle umfassende Test-Suites für deine Workflows
- **[Hello nf-core](../../hello_nf-core/)**: Baue produktionsreife Pipelines mit Community-Standards

Die wahre Stärke dieser IDE-Funktionen zeigt sich, wenn du an größeren, komplexeren Projekten arbeitest. Integriere sie schrittweise in deinen Workflow – nach wenigen Sitzungen werden sie zur zweiten Natur und verändern, wie du an die Nextflow-Entwicklung herangehst.

Vom Erkennen von Fehlern, bevor sie dich aufhalten, bis hin zur mühelosen Navigation in komplexen Codebasen – diese Tools machen dich zu einer selbstsichereren und effizienteren Entwickler\*in.

Viel Spaß beim Coden!
