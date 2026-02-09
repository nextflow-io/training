# Teil 3: Ein nf-core-Modul verwenden

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } KI-gestützte Übersetzung - [mehr erfahren & Verbesserungen vorschlagen](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

In diesem dritten Teil des Hello nf-core Trainingskurses zeigen wir dir, wie du ein vorhandenes nf-core-Modul findest, installierst und in deiner Pipeline verwendest.

Einer der großen Vorteile der Arbeit mit nf-core ist die Möglichkeit, vorgefertigte, getestete Module aus dem [nf-core/modules](https://github.com/nf-core/modules) Repository zu nutzen.
Anstatt jeden Prozess von Grund auf neu zu schreiben, kannst du von der Community gepflegte Module installieren und verwenden, die Best Practices folgen.

Um zu zeigen, wie das funktioniert, ersetzen wir das benutzerdefinierte `collectGreetings`-Modul durch das `cat/cat`-Modul aus nf-core/modules in der `core-hello`-Pipeline.

??? info "Wie du von diesem Abschnitt aus beginnst"

    Dieser Abschnitt des Kurses setzt voraus, dass du [Teil 2: Hello für nf-core umschreiben](./02_rewrite_hello.md) abgeschlossen hast und eine funktionierende `core-hello`-Pipeline besitzt.

    Falls du Teil 2 nicht abgeschlossen hast oder für diesen Teil neu beginnen möchtest, kannst du die `core-hello-part2`-Lösung als Ausgangspunkt verwenden.
    Führe diesen Befehl im `hello-nf-core/`-Verzeichnis aus:

    ```bash
    cp -r solutions/core-hello-part2 core-hello
    cd core-hello
    ```

    Damit erhältst du eine voll funktionsfähige nf-core-Pipeline, die bereit ist für das Hinzufügen von Modulen.
    Du kannst testen, ob sie erfolgreich läuft, indem du folgenden Befehl ausführst:

    ```bash
    nextflow run . --outdir core-hello-results -profile test,docker --validate_params false
    ```

---

## 1. Ein passendes nf-core-Modul finden und installieren

Zuerst lernen wir, wie man ein vorhandenes nf-core-Modul findet und in unsere Pipeline installiert.

Wir wollen den `collectGreetings`-Prozess ersetzen, der den Unix-`cat`-Befehl verwendet, um mehrere Begrüßungsdateien zu einer zusammenzufügen.
Das Zusammenfügen von Dateien ist eine sehr häufige Operation, daher ist es naheliegend, dass es bereits ein Modul in nf-core für diesen Zweck geben könnte.

Lass uns eintauchen.

### 1.1. Verfügbare Module auf der nf-core-Website durchsuchen

Das nf-core-Projekt pflegt einen zentralen Katalog von Modulen unter [https://nf-co.re/modules](https://nf-co.re/modules).

Navigiere in deinem Webbrowser zur Modulseite und verwende die Suchleiste, um nach 'concatenate' zu suchen.

![Modul-Suchergebnisse](./img/module-search-results.png)

Wie du sehen kannst, gibt es einige Ergebnisse, viele davon Module, die für das Zusammenfügen sehr spezifischer Dateitypen entwickelt wurden.
Unter ihnen solltest du eines namens `cat_cat` sehen, das universell einsetzbar ist.

!!! note "Namenskonvention für Module"

    Der Unterstrich (`_`) wird als Platzhalter für das Schrägstrich-Zeichen (`/`) in Modulnamen verwendet.

    nf-core-Module folgen der Namenskonvention `software/command`, wenn ein Tool mehrere Befehle bereitstellt, wie `samtools/view` (samtools-Paket, view-Befehl) oder `gatk/haplotypecaller` (GATK-Paket, HaplotypeCaller-Befehl).
    Für Tools, die nur einen Hauptbefehl bereitstellen, verwenden Module eine einzelne Ebene wie `fastqc` oder `multiqc`.

Klicke auf das `cat_cat`-Modulfeld, um die Moduldokumentation anzuzeigen.

Die Modulseite zeigt:

- Eine kurze Beschreibung: "A module for concatenation of gzipped or uncompressed files"
- Installationsbefehl: `nf-core modules install cat/cat`
- Eingabe- und Ausgabekanalstruktur
- Verfügbare Parameter

### 1.2. Verfügbare Module über die Kommandozeile auflisten

Alternativ kannst du auch direkt über die Kommandozeile mit nf-core-Tools nach Modulen suchen.

```bash
nf-core modules list remote
```

Dies zeigt eine Liste aller verfügbaren Module im nf-core/modules-Repository an, ist aber etwas weniger praktisch, wenn du den Namen des gesuchten Moduls noch nicht kennst.
Wenn du ihn jedoch kennst, kannst du die Liste an `grep` weiterleiten, um bestimmte Module zu finden:

```bash
nf-core modules list remote | grep 'cat/cat'
```

??? success "Befehlsausgabe"

    ```console
    │ cat/cat
    ```

Beachte jedoch, dass der `grep`-Ansatz nur Ergebnisse mit dem Suchbegriff im Namen findet, was für `cat_cat` nicht funktionieren würde.

### 1.3. Detaillierte Informationen über das Modul abrufen

Um detaillierte Informationen über ein bestimmtes Modul über die Kommandozeile zu sehen, verwende den `info`-Befehl:

```bash
nf-core modules info cat/cat
```

Dies zeigt die Dokumentation über das Modul an, einschließlich seiner Eingaben, Ausgaben und grundlegenden Nutzungsinformationen.

??? success "Befehlsausgabe"

    ```console

                                              ,--./,-.
              ___     __   __   __   ___     /,-._.--~\
        |\ | |__  __ /  ` /  \ |__) |__         }  {
        | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                              `._,._,'

        nf-core/tools version 3.4.1 - https://nf-co.re


    ╭─ Module: cat/cat  ─────────────────────────────────────────────────╮
    │ 🌐 Repository: https://github.com/nf-core/modules.git              │
    │ 🔧 Tools: cat                                                      │
    │ 📖 Description: A module for concatenation of gzipped or           │
    │ uncompressed files                                                 │
    ╰────────────────────────────────────────────────────────────────────╯
                      ╷                                          ╷
    📥 Inputs        │Description                               │Pattern
    ╺━━━━━━━━━━━━━━━━━┿━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┿━━━━━━━╸
    input[0]         │                                          │
    ╶─────────────────┼──────────────────────────────────────────┼───────╴
      meta  (map)     │Groovy Map containing sample information  │
                      │e.g. [ id:'test', single_end:false ]      │
    ╶─────────────────┼──────────────────────────────────────────┼───────╴
      files_in  (file)│List of compressed / uncompressed files   │      *
                      ╵                                          ╵
                          ╷                                 ╷
    📥 Outputs           │Description                      │     Pattern
    ╺━━━━━━━━━━━━━━━━━━━━━┿━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┿━━━━━━━━━━━━╸
    file_out             │                                 │
    ╶─────────────────────┼─────────────────────────────────┼────────────╴
      meta  (map)         │Groovy Map containing sample     │
                          │information                      │
    ╶─────────────────────┼─────────────────────────────────┼────────────╴
      ${prefix}  (file)   │Concatenated file. Will be       │ ${file_out}
                          │gzipped if file_out ends with    │
                          │".gz"                            │
    ╶─────────────────────┼─────────────────────────────────┼────────────╴
    versions             │                                 │
    ╶─────────────────────┼─────────────────────────────────┼────────────╴
      versions.yml  (file)│File containing software versions│versions.yml
                          ╵                                 ╵

    💻  Installation command: nf-core modules install cat/cat

    ```

Dies sind genau die gleichen Informationen, die du auch auf der Website finden kannst.

### 1.4. Das cat/cat-Modul installieren

Jetzt, da wir das gewünschte Modul gefunden haben, müssen wir es zum Quellcode unserer Pipeline hinzufügen.

Die gute Nachricht ist, dass das nf-core-Projekt einige Tools enthält, die diesen Teil einfach machen.
Konkret ermöglicht der Befehl `nf-core modules install`, das Abrufen des Codes und das Verfügbarmachen für dein Projekt in einem einzigen Schritt zu automatisieren.

Navigiere zu deinem Pipeline-Verzeichnis und führe den Installationsbefehl aus:

```bash
cd core-hello
nf-core modules install cat/cat
```

Das Tool fordert dich möglicherweise zuerst auf, einen Repository-Typ anzugeben.
(Falls nicht, springe zu "Schließlich wird das Tool mit der Installation des Moduls fortfahren.")

??? success "Befehlsausgabe"

    ```console

                                          ,--./,-.
          ___     __   __   __   ___     /,-._.--~\
    |\ | |__  __ /  ` /  \ |__) |__         }  {
    | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                          `._,._,'

    nf-core/tools version 3.4.1 - https://nf-co.re


    WARNING  'repository_type' not defined in .nf-core.yml
    ? Is this repository a pipeline or a modules repository? (Use arrow keys)
    » Pipeline
      Modules repository
    ```

Falls ja, drücke Enter, um die Standardantwort (`Pipeline`) zu akzeptieren und fortzufahren.

Das Tool bietet dann an, die Konfiguration deines Projekts zu ändern, um diese Aufforderung in Zukunft zu vermeiden.

??? success "Befehlsausgabe"

    ```console
        INFO     To avoid this prompt in the future, add the 'repository_type' key to your .nf-core.yml file.
        ? Would you like me to add this config now? [y/n] (y):
    ```

Warum nicht von diesem praktischen Tool profitieren!
Drücke Enter, um die Standardantwort (ja) zu akzeptieren.

Schließlich wird das Tool mit der Installation des Moduls fortfahren.

??? success "Befehlsausgabe"

    ```console
    INFO Config added to '.nf-core.yml'
    INFO Reinstalling modules found in 'modules.json' but missing from directory:
    INFO Installing 'cat/cat'
    INFO Use the following statement to include this module:

        include { CAT_CAT } from '../modules/nf-core/cat/cat/main'
    ```

Der Befehl erledigt automatisch:

- Herunterladen der Moduldateien nach `modules/nf-core/cat/cat/`
- Aktualisieren von `modules.json`, um das installierte Modul zu verfolgen
- Bereitstellen der korrekten `include`-Anweisung zur Verwendung im Workflow

!!! tip

    Stelle immer sicher, dass dein aktuelles Arbeitsverzeichnis das Stammverzeichnis deines Pipeline-Projekts ist, bevor du den Modulinstallationsbefehl ausführst.

Lass uns überprüfen, ob das Modul korrekt installiert wurde:

```bash
tree -L 4 modules
```

??? abstract "Verzeichnisinhalt"

    ```console
    modules
    ├── local
    │   ├── collectGreetings.nf
    │   ├── convertToUpper.nf
    │   ├── cowpy.nf
    │   └── sayHello.nf
    └── nf-core
        └── cat
            └── cat
                ├── environment.yml
                ├── main.nf
                ├── meta.yml
                └── tests

    5 directories, 7 files
    ```

Du kannst die Installation auch überprüfen, indem du das nf-core-Dienstprogramm bittest, lokal installierte Module aufzulisten:

```bash
nf-core modules list local
```

??? success "Befehlsausgabe"

    ```console
    INFO     Repository type: pipeline
    INFO     Modules installed in '.':

    ┏━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━┓
    ┃ Module Name ┃ Repository      ┃ Version SHA ┃ Message                                ┃ Date       ┃
    ┡━━━━━━━━━━━━━╇━━━━━━━━━━━━━━━━━╇━━━━━━━━━━━━━╇━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╇━━━━━━━━━━━━┩
    │ cat/cat     │ nf-core/modules │ 41dfa3f     │ update meta.yml of all modules (#8747) │ 2025-07-07 │
    └─────────────┴─────────────────┴─────────────┴────────────────────────────────────────┴────────────┘
    ```

Dies bestätigt, dass das `cat/cat`-Modul jetzt Teil des Quellcodes deines Projekts ist.

Um das neue Modul tatsächlich zu verwenden, müssen wir es jedoch in unsere Pipeline importieren.

### 1.5. Die Modulimporte aktualisieren

Lass uns die `include`-Anweisung für das `collectGreetings`-Modul durch die für `CAT_CAT` im Importbereich des `workflows/hello.nf`-Workflows ersetzen.

Zur Erinnerung: Das Modulinstallationstool hat uns die exakte Anweisung gegeben:

```groovy title="Import-Anweisung vom Installationsbefehl erzeugt"
include { CAT_CAT } from '../modules/nf-core/cat/cat/main'`
```

Beachte, dass die nf-core-Konvention darin besteht, beim Importieren von Modulen Großbuchstaben für Modulnamen zu verwenden.

Öffne [core-hello/workflows/hello.nf](core-hello/workflows/hello.nf) und nimm folgende Ersetzung vor:

=== "Danach"

    ```groovy title="core-hello/workflows/hello.nf" linenums="1" hl_lines="10"
    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    include { paramsSummaryMap       } from 'plugin/nf-schema'
    include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
    include { sayHello               } from '../modules/local/sayHello.nf'
    include { convertToUpper         } from '../modules/local/convertToUpper.nf'
    include { CAT_CAT                } from '../modules/nf-core/cat/cat/main'
    include { cowpy                  } from '../modules/local/cowpy.nf'
    ```

=== "Vorher"

    ```groovy title="core-hello/workflows/hello.nf" linenums="1" hl_lines="10"
    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    include { paramsSummaryMap       } from 'plugin/nf-schema'
    include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
    include { sayHello               } from '../modules/local/sayHello.nf'
    include { convertToUpper         } from '../modules/local/convertToUpper.nf'
    include { collectGreetings       } from '../modules/local/collectGreetings.nf'
    include { cowpy                  } from '../modules/local/cowpy.nf'
    ```

Beachte, wie sich der Pfad für das nf-core-Modul von den lokalen Modulen unterscheidet:

- **nf-core-Modul**: `'../modules/nf-core/cat/cat/main'` (verweist auf `main.nf`)
- **Lokales Modul**: `'../modules/local/collectGreetings.nf'` (Einzeldatei-Verweis)

Das Modul ist jetzt für den Workflow verfügbar, also müssen wir nur noch den Aufruf von `collectGreetings` durch `CAT_CAT` ersetzen. Richtig?

Nicht so schnell.

An diesem Punkt könntest du versucht sein, direkt mit der Codebearbeitung zu beginnen, aber es lohnt sich, einen Moment innezuhalten und sorgfältig zu prüfen, was das neue Modul erwartet und was es produziert.

Wir werden das als separaten Abschnitt behandeln, da es einen neuen Mechanismus beinhaltet, den wir noch nicht behandelt haben: Metadaten-Maps.

!!! note

    Du kannst optional die Datei `collectGreetings.nf` löschen:

    ```bash
    rm modules/local/collectGreetings.nf
    ```

    Möglicherweise möchtest du sie jedoch als Referenz behalten, um die Unterschiede zwischen lokalen und nf-core-Modulen zu verstehen.

### Fazit

Du weißt, wie du ein nf-core-Modul findest und es für dein Projekt verfügbar machst.

### Wie geht es weiter?

Bewerte, was ein neues Modul benötigt, und identifiziere wichtige Änderungen, die erforderlich sind, um es in eine Pipeline zu integrieren.

---

## 2. Die Anforderungen des neuen Moduls bewerten

Konkret müssen wir die **Schnittstelle** des Moduls untersuchen, d.h. seine Eingabe- und Ausgabedefinitionen, und sie mit der Schnittstelle des Moduls vergleichen, das wir ersetzen möchten.
Dies ermöglicht uns festzustellen, ob wir das neue Modul einfach als direkten Ersatz behandeln können oder ob wir einige Anpassungen an der Verkabelung vornehmen müssen.

Idealerweise solltest du dies tun, _bevor_ du das Modul überhaupt installierst, aber besser spät als nie.
(Übrigens gibt es einen `uninstall`-Befehl, um Module loszuwerden, die du nicht mehr möchtest.)

!!! note

    Der CAT_CAT-Prozess enthält eine ziemlich clevere Handhabung verschiedener Kompressionstypen, Dateierweiterungen usw., die für das, was wir dir hier zeigen möchten, nicht streng relevant sind. Daher ignorieren wir das meiste davon und konzentrieren uns nur auf die wichtigen Teile.

### 2.1. Die Schnittstellen der beiden Module vergleichen

Zur Erinnerung: So sieht die Schnittstelle unseres `collectGreetings`-Moduls aus:

```groovy title="modules/local/collectGreetings.nf (Auszug)" linenums="1" hl_lines="6-7 10"
process collectGreetings {

    publishDir 'results', mode: 'copy'

    input:
        path input_files
        val batch_name

    output:
        path "COLLECTED-${batch_name}-output.txt" , emit: outfile
```

Das `collectGreetings`-Modul nimmt zwei Eingaben entgegen:

- `input_files` enthält eine oder mehrere zu verarbeitende Eingabedateien;
- `batch_name` ist ein Wert, den wir verwenden, um der Ausgabedatei einen laufspezifischen Namen zuzuweisen, was eine Form von Metadaten ist.

Nach Abschluss gibt `collectGreetings` einen einzelnen Dateipfad aus, der mit dem Tag `outfile` ausgegeben wird.

Im Vergleich dazu ist die Schnittstelle des `cat/cat`-Moduls komplexer:

```groovy title="modules/nf-core/cat/cat/main.nf (Auszug)" linenums="1" hl_lines="11 14"
process CAT_CAT {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pigz:2.3.4' :
        'biocontainers/pigz:2.3.4' }"

    input:
    tuple val(meta), path(files_in)

    output:
    tuple val(meta), path("${prefix}"), emit: file_out
    path "versions.yml"               , emit: versions
```

Das CAT_CAT-Modul nimmt eine einzelne Eingabe entgegen, aber diese Eingabe ist ein Tupel, das zwei Dinge enthält:

- `meta` ist eine Struktur, die Metadaten enthält, genannt Metamap;
- `files_in` enthält eine oder mehrere zu verarbeitende Eingabedateien, entsprechend `input_files` von `collectGreetings`.

Nach Abschluss liefert CAT_CAT seine Ausgaben in zwei Teilen:

- Ein weiteres Tupel, das die Metamap und die zusammengefügte Ausgabedatei enthält, ausgegeben mit dem Tag `file_out`;
- Eine `versions.yml`-Datei, die Informationen über die verwendete Softwareversion erfasst, ausgegeben mit dem Tag `versions`.

Beachte auch, dass die Ausgabedatei standardmäßig basierend auf einem Identifikator benannt wird, der Teil der Metadaten ist (Code hier nicht gezeigt).

Das mag viel erscheinen, wenn man nur den Code betrachtet, daher hier ein Diagramm, das dir hilft zu visualisieren, wie alles zusammenpasst.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nf-core/img/module_comparison.svg"
</figure>

Du kannst sehen, dass die beiden Module ähnliche Eingabeanforderungen in Bezug auf den Inhalt haben (eine Menge von Eingabedateien plus einige Metadaten), aber sehr unterschiedliche Erwartungen daran, wie dieser Inhalt verpackt ist.
Wenn wir die Versionsdatei vorerst ignorieren, ist ihre Hauptausgabe ebenfalls gleichwertig (eine zusammengefügte Datei), außer dass CAT_CAT auch die Metamap zusammen mit der Ausgabedatei ausgibt.

Die Verpackungsunterschiede werden ziemlich einfach zu handhaben sein, wie du gleich sehen wirst.
Um jedoch den Metamap-Teil zu verstehen, müssen wir dir etwas zusätzlichen Kontext geben.

### 2.2. Metamaps verstehen

Wir haben dir gerade gesagt, dass das CAT_CAT-Modul eine Metadaten-Map als Teil seines Eingabe-Tupels erwartet.
Lass uns ein paar Minuten nehmen, um genauer zu betrachten, was das ist.

Die **Metadaten-Map**, oft kurz als **Metamap** bezeichnet, ist eine Groovy-Map, die Informationen über Dateneinheiten enthält.
Im Kontext von Nextflow-Pipelines können Dateneinheiten alles sein, was du möchtest: einzelne Proben, Probenstapel oder ganze Datensätze.

Per Konvention wird eine nf-core-Metamap `meta` genannt und enthält das erforderliche Feld `id`, das zur Benennung von Ausgaben und zur Verfolgung von Dateneinheiten verwendet wird.

Zum Beispiel könnte eine typische Metadaten-Map so aussehen:

```groovy title="Beispiel einer Metamap auf Probenebene"
[id: 'sample1', single_end: false, strandedness: 'forward']
```

Oder in einem Fall, in dem die Metadaten auf Stapelebene angehängt sind:

```groovy title="Beispiel einer Metamap auf Stapelebene"
[id: 'batch1', date: '25.10.01']
```

Lass uns dies nun in den Kontext des `CAT_CAT`-Prozesses setzen, der erwartet, dass die Eingabedateien in ein Tupel mit einer Metamap verpackt werden, und die Metamap auch als Teil des Ausgabe-Tupels ausgibt.

```groovy title="modules/nf-core/cat/cat/main.nf (Auszug)" linenums="1" hl_lines="2 5"
input:
tuple val(meta), path(files_in)

output:
tuple val(meta), path("${prefix}"), emit: file_out
```

Dadurch durchläuft jede Dateneinheit die Pipeline mit den relevanten Metadaten.
Nachfolgende Prozesse können dann auch problemlos auf diese Metadaten zugreifen.

Erinnerst du dich, wie wir dir gesagt haben, dass die von `CAT_CAT` ausgegebene Datei basierend auf einem Identifikator benannt wird, der Teil der Metadaten ist?
Dies ist der relevante Code:

```groovy title="modules/nf-core/cat/cat/main.nf (Auszug)" linenums="35"
prefix   = task.ext.prefix ?: "${meta.id}${getFileSuffix(file_list[0])}"
```

Dies bedeutet ungefähr Folgendes: Wenn ein `prefix` über das externe Aufgabenparametersystem (`task.ext`) bereitgestellt wird, verwende diesen zur Benennung der Ausgabedatei; andernfalls erstelle einen mit `${meta.id}`, was dem `id`-Feld in der Metamap entspricht.

Du kannst dir den Eingabekanal vorstellen, der in dieses Modul kommt, mit Inhalten wie diesem:

```groovy title="Beispiel für Eingabekanalinhalte"
ch_input = [[[id: 'batch1', date: '25.10.01'], ['file1A.txt', 'file1B.txt']],
            [[id: 'batch2', date: '25.10.26'], ['file2A.txt', 'file2B.txt']],
            [[id: 'batch3', date: '25.11.14'], ['file3A.txt', 'file3B.txt']]]
```

Dann kommen die Ausgabekanalinhalte so heraus:

```groovy title="Beispiel für Ausgabekanalinhalte"
ch_input = [[[id: 'batch1', date: '25.10.01'], 'batch1.txt'],
            [[id: 'batch2', date: '25.10.26'], 'batch2.txt'],
            [[id: 'batch3', date: '25.11.14'], 'batch3.txt']]
```

Wie bereits erwähnt, ist das `tuple val(meta), path(files_in)`-Eingabesetup ein Standardmuster, das in allen nf-core-Modulen verwendet wird.

Hoffentlich kannst du langsam sehen, wie nützlich dies sein kann.
Es ermöglicht dir nicht nur, Ausgaben basierend auf Metadaten zu benennen, sondern du kannst auch Dinge tun wie unterschiedliche Parameterwerte anwenden, und in Kombination mit bestimmten Operatoren kannst du sogar Daten gruppieren, sortieren oder filtern, während sie durch die Pipeline fließen.

!!! note "Mehr über Metadaten erfahren"

    Für eine umfassende Einführung in die Arbeit mit Metadaten in Nextflow-Workflows, einschließlich wie man Metadaten aus Samplesheets liest und sie zur Anpassung der Verarbeitung verwendet, siehe die [Metadaten in Workflows](../side_quests/metadata) Side Quest.

### 2.3. Vorzunehmende Änderungen zusammenfassen

Basierend auf dem, was wir überprüft haben, sind dies die wichtigsten Änderungen, die wir an unserer Pipeline vornehmen müssen, um das `cat/cat`-Modul zu nutzen:

- Eine Metamap erstellen, die den Stapelnamen enthält;
- Die Metamap in ein Tupel mit der Menge der zusammenzufügenden Eingabedateien verpacken (die aus `convertToUpper` kommen);
- Den Aufruf von `collectGreetings()` zu `CAT_CAT` ändern;
- Die Ausgabedatei aus dem vom `CAT_CAT`-Prozess erzeugten Tupel extrahieren, bevor sie an `cowpy` übergeben wird.

Das sollte reichen! Jetzt, da wir einen Plan haben, sind wir bereit loszulegen.

### Fazit

Du weißt, wie du die Eingabe- und Ausgabeschnittstelle eines neuen Moduls bewertest, um seine Anforderungen zu identifizieren, und du hast gelernt, wie Metamaps von nf-core-Pipelines verwendet werden, um Metadaten eng mit den Daten zu verknüpfen, während sie durch eine Pipeline fließen.

### Wie geht es weiter?

Integriere das neue Modul in einen Workflow.

---

## 3. CAT_CAT in den `hello.nf`-Workflow integrieren

Jetzt, da du alles über Metamaps weißt (oder zumindest genug für die Zwecke dieses Kurses), ist es an der Zeit, die oben skizzierten Änderungen tatsächlich umzusetzen.

Der Klarheit halber werden wir dies aufschlüsseln und jeden Schritt separat behandeln.

!!! note

    Alle unten gezeigten Änderungen werden an der Workflow-Logik im `main`-Block in der Workflow-Datei `core-hello/workflows/hello.nf` vorgenommen.

### 3.1. Eine Metadaten-Map erstellen

Zuerst müssen wir eine Metadaten-Map für `CAT_CAT` erstellen, wobei wir bedenken, dass nf-core-Module mindestens ein `id`-Feld in der Metamap benötigen.

Da wir keine anderen Metadaten benötigen, können wir es einfach halten und so etwas verwenden:

```groovy title="Syntaxbeispiel"
def cat_meta = [id: 'test']
```

Außer dass wir den `id`-Wert nicht fest codieren wollen; wir wollen den Wert des `params.batch`-Parameters verwenden.
Also wird der Code zu:

```groovy title="Syntaxbeispiel"
def cat_meta = [id: params.batch]
```

Ja, es ist buchstäblich so einfach, eine grundlegende Metamap zu erstellen.

Lass uns diese Zeilen nach dem `convertToUpper`-Aufruf hinzufügen und den `collectGreetings`-Aufruf entfernen:

=== "Danach"

    ```groovy title="core-hello/workflows/hello.nf" linenums="26" hl_lines="7-8"
        // Eine Begrüßung ausgeben
        sayHello(ch_samplesheet)

        // Die Begrüßung in Großbuchstaben umwandeln
        convertToUpper(sayHello.out)

        // Metadaten-Map mit Stapelname als ID erstellen
        def cat_meta = [ id: params.batch ]

        // ASCII-Art der Begrüßungen mit cowpy generieren
        cowpy(collectGreetings.out.outfile, params.character)
    ```

=== "Vorher"

    ```groovy title="core-hello/workflows/hello.nf" linenums="26" hl_lines="7-8"
        // Eine Begrüßung ausgeben
        sayHello(ch_samplesheet)

        // Die Begrüßung in Großbuchstaben umwandeln
        convertToUpper(sayHello.out)

        // Alle Begrüßungen in einer Datei sammeln
        collectGreetings(convertToUpper.out.collect(), params.batch)

        // ASCII-Art der Begrüßungen mit cowpy generieren
        cowpy(collectGreetings.out.outfile, params.character)
    ```

Dies erstellt eine einfache Metadaten-Map, bei der die `id` auf unseren Stapelnamen gesetzt ist (der `test` sein wird, wenn das Testprofil verwendet wird).

### 3.2. Einen Kanal mit Metadaten-Tupeln erstellen

Als Nächstes transformiere den Kanal von Dateien in einen Kanal von Tupeln, die Metadaten und Dateien enthalten:

=== "Danach"

    ```groovy title="core-hello/workflows/hello.nf" linenums="26" hl_lines="10-11"
        // Eine Begrüßung ausgeben
        sayHello(ch_samplesheet)

        // Die Begrüßung in Großbuchstaben umwandeln
        convertToUpper(sayHello.out)

        // Metadaten-Map mit Stapelname als ID erstellen
        def cat_meta = [ id: params.batch ]

        // Einen Kanal mit Metadaten und Dateien im Tupelformat erstellen
        ch_for_cat = convertToUpper.out.collect().map { files -> tuple(cat_meta, files) }

        // ASCII-Art der Begrüßungen mit cowpy generieren
        cowpy(collectGreetings.out.outfile, params.character)
    ```

=== "Vorher"

    ```groovy title="core-hello/workflows/hello.nf" linenums="26"
        // Eine Begrüßung ausgeben
        sayHello(ch_samplesheet)

        // Die Begrüßung in Großbuchstaben umwandeln
        convertToUpper(sayHello.out)

        // Metadaten-Map mit Stapelname als ID erstellen
        def cat_meta = [ id: params.batch ]

        // ASCII-Art der Begrüßungen mit cowpy generieren
        cowpy(collectGreetings.out.outfile, params.character)
    ```

Die Zeile, die wir hinzugefügt haben, erreicht zwei Dinge:

- `.collect()` sammelt alle Dateien aus der `convertToUpper`-Ausgabe in einer einzigen Liste
- `.map { files -> tuple(cat_meta, files) }` erstellt ein Tupel von `[metadata, files]` im Format, das `CAT_CAT` erwartet

Das ist alles, was wir tun müssen, um das Eingabe-Tupel für `CAT_CAT` einzurichten.

### 3.3. Das CAT_CAT-Modul aufrufen

Rufe nun `CAT_CAT` auf dem neu erstellten Kanal auf:

=== "Danach"

    ```groovy title="core-hello/workflows/hello.nf" linenums="26" hl_lines="13-14"
        // Eine Begrüßung ausgeben
        sayHello(ch_samplesheet)

        // Die Begrüßung in Großbuchstaben umwandeln
        convertToUpper(sayHello.out)

        // Metadaten-Map mit Stapelname als ID erstellen
        def cat_meta = [ id: params.batch ]

        // Einen Kanal mit Metadaten und Dateien im Tupelformat erstellen
        ch_for_cat = convertToUpper.out.collect().map { files -> tuple(cat_meta, files) }

        // Dateien mit dem nf-core cat/cat-Modul zusammenfügen
        CAT_CAT(ch_for_cat)

        // ASCII-Art der Begrüßungen mit cowpy generieren
        cowpy(collectGreetings.out.outfile, params.character)
    ```

=== "Vorher"

    ```groovy title="core-hello/workflows/hello.nf" linenums="26"
        // Eine Begrüßung ausgeben
        sayHello(ch_samplesheet)

        // Die Begrüßung in Großbuchstaben umwandeln
        convertToUpper(sayHello.out)

        // Metadaten-Map mit Stapelname als ID erstellen
        def cat_meta = [ id: params.batch ]

        // Einen Kanal mit Metadaten und Dateien im Tupelformat erstellen
        ch_for_cat = convertToUpper.out.collect().map { files -> tuple(cat_meta, files) }

        // ASCII-Art der Begrüßungen mit cowpy generieren
        cowpy(collectGreetings.out.outfile, params.character)
    ```

Dies vervollständigt den kniffligsten Teil dieser Ersetzung, aber wir sind noch nicht ganz fertig: Wir müssen noch aktualisieren, wie wir die zusammengefügte Ausgabe an den `cowpy`-Prozess übergeben.

### 3.4. Die Ausgabedatei aus dem Tupel für `cowpy` extrahieren

Zuvor hat der `collectGreetings`-Prozess einfach eine Datei erzeugt, die wir direkt an `cowpy` übergeben konnten.
Der `CAT_CAT`-Prozess erzeugt jedoch ein Tupel, das zusätzlich zur Ausgabedatei die Metamap enthält.

Da `cowpy` noch keine Metadaten-Tupel akzeptiert (wir werden dies im nächsten Teil des Kurses beheben), müssen wir die Ausgabedatei aus dem von `CAT_CAT` erzeugten Tupel extrahieren, bevor wir sie an `cowpy` übergeben:

=== "Danach"

    ```groovy title="core-hello/workflows/hello.nf" linenums="26" hl_lines="16-17 20"
        // Eine Begrüßung ausgeben
        sayHello(ch_samplesheet)

        // Die Begrüßung in Großbuchstaben umwandeln
        convertToUpper(sayHello.out)

        // Metadaten-Map mit Stapelname als ID erstellen
        def cat_meta = [ id: params.batch ]

        // Einen Kanal mit Metadaten und Dateien im Tupelformat erstellen
        ch_for_cat = convertToUpper.out.collect().map { files -> tuple(cat_meta, files) }

        // Die Begrüßungen zusammenfügen
        CAT_CAT(ch_for_cat)

        // Die Datei aus dem Tupel extrahieren, da cowpy noch keine Metadaten verwendet
        ch_for_cowpy = CAT_CAT.out.file_out.map{ meta, file -> file }

        // ASCII-Art der Begrüßungen mit cowpy generieren
        cowpy(ch_for_cowpy, params.character)
    ```

=== "Vorher"

    ```groovy title="core-hello/workflows/hello.nf" linenums="26" hl_lines="17"
        // Eine Begrüßung ausgeben
        sayHello(ch_samplesheet)

        // Die Begrüßung in Großbuchstaben umwandeln
        convertToUpper(sayHello.out)

        // Metadaten-Map mit Stapelname als ID erstellen
        def cat_meta = [ id: params.batch ]

        // Einen Kanal mit Metadaten und Dateien im Tupelformat erstellen
        ch_for_cat = convertToUpper.out.collect().map { files -> tuple(cat_meta, files) }

        // Die Begrüßungen zusammenfügen
        CAT_CAT(ch_for_cat)

        // ASCII-Art der Begrüßungen mit cowpy generieren
        cowpy(collectGreetings.out.outfile, params.character)
    ```

Die `.map{ meta, file -> file }`-Operation extrahiert die Datei aus dem `[metadata, file]`-Tupel, das von `CAT_CAT` erzeugt wurde, in einen neuen Kanal, `ch_for_cowpy`.

Dann ist es nur noch eine Frage der Übergabe von `ch_for_cowpy` an `cowpy` anstelle von `collectGreetings.out.outfile` in dieser letzten Zeile.

!!! note

    Im nächsten Teil des Kurses werden wir `cowpy` aktualisieren, um direkt mit Metadaten-Tupeln zu arbeiten, sodass dieser Extraktionsschritt nicht mehr notwendig sein wird.

### 3.5. Den Workflow testen

Lass uns testen, ob der Workflow mit dem neu integrierten `cat/cat`-Modul funktioniert:

```bash
nextflow run . --outdir core-hello-results -profile test,docker --validate_params false
```

Dies sollte ziemlich schnell laufen.

??? success "Befehlsausgabe"

    ```console
    N E X T F L O W ~ version 25.04.3

        Launching `./main.nf` [evil_pike] DSL2 - revision: b9e9b3b8de

        Input/output options
          input                     : /workspaces/training/hello-nf-core/core-hello/assets/greetings.csv
          outdir                    : core-hello-results

        Institutional config options
          config_profile_name       : Test profile
          config_profile_description: Minimal test dataset to check pipeline function

        Generic options
          validate_params           : false
          trace_report_suffix       : 2025-10-30_18-50-58

        Core Nextflow options
          runName                   : evil_pike
          containerEngine           : docker
          launchDir                 : /workspaces/training/hello-nf-core/core-hello
          workDir                   : /workspaces/training/hello-nf-core/core-hello/work
          projectDir                : /workspaces/training/hello-nf-core/core-hello
          userName                  : root
          profile                   : test,docker
          configFiles               : /workspaces/training/hello-nf-core/core-hello/nextflow.config

        !! Only displaying parameters that differ from the pipeline defaults !!
        ------------------------------------------------------
        executor >  local (8)
        [b3/f005fd] CORE_HELLO:HELLO:sayHello (3)       [100%] 3 of 3 ✔
        [08/f923d0] CORE_HELLO:HELLO:convertToUpper (3) [100%] 3 of 3 ✔
        [34/3729a9] CORE_HELLO:HELLO:CAT_CAT (test)     [100%] 1 of 1 ✔
        [24/df918a] CORE_HELLO:HELLO:cowpy              [100%] 1 of 1 ✔
        -[core/hello] Pipeline completed successfully-
    ```

Beachte, dass `CAT_CAT` jetzt in der Prozessausführungsliste anstelle von `collectGreetings` erscheint.

Und das war's! Wir verwenden jetzt ein robustes, von der Community kuratiertes Modul anstelle von benutzerdefiniertem Prototyp-Code für diesen Schritt in der Pipeline.

### Fazit

Du weißt jetzt, wie du:

- nf-core-Module findest und installierst
- Die Anforderungen eines nf-core-Moduls bewertest
- Eine einfache Metadaten-Map zur Verwendung mit einem nf-core-Modul erstellst
- Ein nf-core-Modul in deinen Workflow integrierst

### Wie geht es weiter?

Lerne, deine lokalen Module anzupassen, um nf-core-Konventionen zu folgen.
Wir zeigen dir auch, wie du neue nf-core-Module aus einer Vorlage mit den nf-core-Tools erstellst.
