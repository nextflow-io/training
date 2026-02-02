# Teil 3: Ein nf-core-Modul verwenden

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } KI-gestützte Übersetzung - [mehr erfahren & Verbesserungen vorschlagen](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

In diesem dritten Teil des Hello nf-core Trainingskurses zeigen wir dir, wie du ein existierendes nf-core-Modul findest, installierst und in deiner Pipeline verwendest.

Einer der großen Vorteile der Arbeit mit nf-core ist die Möglichkeit, vorgefertigte, getestete Module aus dem [nf-core/modules](https://github.com/nf-core/modules) Repository zu nutzen.
Anstatt jeden Prozess von Grund auf neu zu schreiben, kannst du von der Community gepflegte Module installieren und verwenden, die Best Practices folgen.

Um zu demonstrieren, wie das funktioniert, werden wir das eigene `collectGreetings`-Modul durch das `cat/cat`-Modul aus nf-core/modules in der `core-hello`-Pipeline ersetzen.

??? info "Wie du von diesem Abschnitt aus beginnst"

    Dieser Abschnitt des Kurses setzt voraus, dass du [Teil 2: Hello für nf-core umschreiben](./02_rewrite_hello.md) abgeschlossen hast und eine funktionierende `core-hello`-Pipeline besitzt.

    Falls du Teil 2 nicht abgeschlossen hast oder für diesen Teil neu beginnen möchtest, kannst du die `core-hello-part2`-Lösung als Ausgangspunkt verwenden.
    Führe diesen Befehl innerhalb des `hello-nf-core/`-Verzeichnisses aus:

    ```bash
    cp -r solutions/core-hello-part2 core-hello
    cd core-hello
    ```

    Das gibt dir eine voll funktionsfähige nf-core-Pipeline, die bereit ist für das Hinzufügen von Modulen.
    Du kannst testen, ob sie erfolgreich läuft, indem du folgenden Befehl ausführst:

    ```bash
    nextflow run . --outdir core-hello-results -profile test,docker --validate_params false
    ```

---

## 1. Ein passendes nf-core-Modul finden und installieren

Zuerst lernen wir, wie man ein existierendes nf-core-Modul findet und in unsere Pipeline installiert.

Wir werden versuchen, den `collectGreetings`-Prozess zu ersetzen, der den Unix-`cat`-Befehl verwendet, um mehrere Grußdateien in eine einzige zu verketten.
Das Verketten von Dateien ist eine sehr häufige Operation, daher liegt es nahe, dass es bereits ein Modul in nf-core gibt, das für diesen Zweck entwickelt wurde.

Lass uns eintauchen.

### 1.1. Verfügbare Module auf der nf-core-Website durchsuchen

Das nf-core-Projekt pflegt einen zentralen Katalog von Modulen unter [https://nf-co.re/modules](https://nf-co.re/modules).

Navigiere zur Modulseite in deinem Webbrowser und verwende die Suchleiste, um nach 'concatenate' zu suchen.

![Modul-Suchergebnisse](./img/module-search-results.png)

Wie du sehen kannst, gibt es einige Ergebnisse, viele davon Module, die für das Verketten sehr spezifischer Dateitypen entwickelt wurden.
Unter ihnen solltest du eines namens `cat_cat` sehen, das universell einsetzbar ist.

!!! note "Namenskonvention für Module"

    Der Unterstrich (`_`) wird als Platzhalter für den Schrägstrich (`/`) in Modulnamen verwendet.

    nf-core-Module folgen der Namenskonvention `software/command`, wenn ein Tool mehrere Befehle bereitstellt, wie `samtools/view` (samtools-Paket, view-Befehl) oder `gatk/haplotypecaller` (GATK-Paket, HaplotypeCaller-Befehl).
    Für Tools, die nur einen Hauptbefehl bereitstellen, verwenden Module eine einzelne Ebene wie `fastqc` oder `multiqc`.

Klicke auf das `cat_cat`-Modulfeld, um die Moduldokumentation anzuzeigen.

Die Modulseite zeigt:

- Eine kurze Beschreibung: "A module for concatenation of gzipped or uncompressed files"
- Installationsbefehl: `nf-core modules install cat/cat`
- Struktur der Eingabe- und Ausgabekanäle
- Verfügbare Parameter

### 1.2. Verfügbare Module von der Kommandozeile auflisten

Alternativ kannst du auch direkt von der Kommandozeile aus mit nf-core-Tools nach Modulen suchen.

```bash
nf-core modules list remote
```

Dies zeigt eine Liste aller verfügbaren Module im nf-core/modules Repository an, ist aber etwas weniger praktisch, wenn du den Namen des Moduls, nach dem du suchst, nicht bereits kennst.
Wenn du ihn jedoch kennst, kannst du die Liste mit `grep` filtern, um spezifische Module zu finden:

```bash
nf-core modules list remote | grep 'cat/cat'
```

??? success "Befehlsausgabe"

    ```console
    │ cat/cat
    ```

Beachte nur, dass der `grep`-Ansatz nur Ergebnisse mit dem Suchbegriff im Namen herausfiltern wird, was bei `cat_cat` nicht funktionieren würde.

### 1.3. Detaillierte Informationen über das Modul erhalten

Um detaillierte Informationen über ein bestimmtes Modul von der Kommandozeile aus zu sehen, verwende den `info`-Befehl:

```bash
nf-core modules info cat/cat
```

Dies zeigt die Dokumentation über das Modul an, einschließlich seiner Eingaben, Ausgaben und grundlegender Nutzungsinformationen.

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

Das ist genau die gleiche Information, die du auch auf der Website findest.

### 1.4. Das cat/cat-Modul installieren

Jetzt, da wir das gewünschte Modul gefunden haben, müssen wir es zum Quellcode unserer Pipeline hinzufügen.

Die gute Nachricht ist, dass das nf-core-Projekt Werkzeuge enthält, die diesen Teil einfach machen.
Speziell der `nf-core modules install`-Befehl ermöglicht es, das Abrufen des Codes und das Verfügbarmachen für dein Projekt in einem einzigen Schritt zu automatisieren.

Navigiere zu deinem Pipeline-Verzeichnis und führe den Installationsbefehl aus:

```bash
cd core-hello
nf-core modules install cat/cat
```

Das Tool wird dich möglicherweise zuerst auffordern, einen Repository-Typ anzugeben.
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

Das Tool wird dann anbieten, die Konfiguration deines Projekts zu ändern, um diese Aufforderung in Zukunft zu vermeiden.

??? success "Befehlsausgabe"

    ```console
        INFO     To avoid this prompt in the future, add the 'repository_type' key to your .nf-core.yml file.
        ? Would you like me to add this config now? [y/n] (y):
    ```

Warum nicht diese praktische Automatisierung nutzen!
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
- Aktualisierung der `modules.json`, um das installierte Modul zu verfolgen
- Bereitstellung der korrekten `include`-Anweisung zur Verwendung im Workflow

!!! tip

    Stelle immer sicher, dass dein aktuelles Arbeitsverzeichnis das Wurzelverzeichnis deines Pipeline-Projekts ist, bevor du den Modulinstallationsbefehl ausführst.

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

Dies bestätigt, dass das `cat/cat`-Modul nun Teil des Quellcodes deines Projekts ist.

Um das neue Modul jedoch tatsächlich zu verwenden, müssen wir es in unsere Pipeline importieren.

### 1.5. Die Modulimporte aktualisieren

Lass uns die `include`-Anweisung für das `collectGreetings`-Modul durch die für `CAT_CAT` im Importbereich des `workflows/hello.nf` Workflows ersetzen.

Zur Erinnerung, das Modulinstallationstool hat uns die exakte Anweisung gegeben, die wir verwenden sollen:

```groovy title="Import-Anweisung, die vom Installationsbefehl erzeugt wurde"
include { CAT_CAT } from '../modules/nf-core/cat/cat/main'`
```

Beachte, dass die nf-core-Konvention darin besteht, Großbuchstaben für Modulnamen beim Importieren zu verwenden.

Öffne [core-hello/workflows/hello.nf](core-hello/workflows/hello.nf) und nimm folgende Ersetzung vor:

=== "Nachher"

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
- **Lokales Modul**: `'../modules/local/collectGreetings.nf'` (einzelne Dateireferenz)

Das Modul ist jetzt für den Workflow verfügbar, also müssen wir nur noch den Aufruf von `collectGreetings` durch `CAT_CAT` ersetzen. Richtig?

Nicht so schnell.

An diesem Punkt könntest du versucht sein, direkt einzusteigen und Code zu bearbeiten, aber es lohnt sich, einen Moment innezuhalten und sorgfältig zu prüfen, was das neue Modul erwartet und was es produziert.

Wir werden das als separaten Abschnitt behandeln, weil es einen neuen Mechanismus beinhaltet, den wir noch nicht behandelt haben: Metadaten-Maps.

!!! note

    Du kannst optional die Datei `collectGreetings.nf` löschen:

    ```bash
    rm modules/local/collectGreetings.nf
    ```

    Du möchtest sie jedoch vielleicht als Referenz behalten, um die Unterschiede zwischen lokalen und nf-core-Modulen zu verstehen.

### Zusammenfassung

Du weißt jetzt, wie du ein nf-core-Modul findest und es für dein Projekt verfügbar machst.

### Wie geht's weiter?

Bewerte, was ein neues Modul benötigt, und identifiziere wichtige Änderungen, die nötig sind, um es in eine Pipeline zu integrieren.

---

## 2. Die Anforderungen des neuen Moduls bewerten

Konkret müssen wir die **Schnittstelle** des Moduls untersuchen, d.h. seine Eingabe- und Ausgabedefinitionen, und sie mit der Schnittstelle des Moduls vergleichen, das wir ersetzen möchten.
Dadurch können wir feststellen, ob wir das neue Modul einfach als direkten Ersatz behandeln können oder ob wir einige Anpassungen in der Verkabelung vornehmen müssen.

Idealerweise solltest du das tun, _bevor_ du das Modul überhaupt installierst, aber hey, besser spät als nie.
(Übrigens gibt es einen `uninstall`-Befehl, um Module loszuwerden, die du nicht mehr möchtest.)

!!! note

    Der CAT_CAT-Prozess enthält eine ziemlich clevere Handhabung verschiedener Komprimierungstypen, Dateierweiterungen usw., die für das, was wir dir hier zeigen wollen, nicht streng relevant sind, daher werden wir das meiste davon ignorieren und uns nur auf die wichtigen Teile konzentrieren.

### 2.1. Die Schnittstellen der beiden Module vergleichen

Zur Erinnerung, so sieht die Schnittstelle zu unserem `collectGreetings`-Modul aus:

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

- `input_files` enthält eine oder mehrere Eingabedateien zur Verarbeitung;
- `batch_name` ist ein Wert, den wir verwenden, um der Ausgabedatei einen ausführungsspezifischen Namen zuzuweisen, was eine Form von Metadaten darstellt.

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
- `files_in` enthält eine oder mehrere Eingabedateien zur Verarbeitung, entsprechend `collectGreetings`' `input_files`.

Nach Abschluss liefert CAT_CAT seine Ausgaben in zwei Teilen:

- Ein weiteres Tupel, das die Metamap und die verkettete Ausgabedatei enthält, ausgegeben mit dem Tag `file_out`;
- Eine `versions.yml`-Datei, die Informationen über die verwendete Softwareversion erfasst, ausgegeben mit dem Tag `versions`.

Beachte auch, dass die Ausgabedatei standardmäßig basierend auf einem Identifikator benannt wird, der Teil der Metadaten ist (Code hier nicht gezeigt).

Das mag wie viel erscheinen, wenn man nur den Code betrachtet, daher hier ein Diagramm, das dir hilft zu visualisieren, wie alles zusammenpasst.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nf-core/img/module_comparison.svg"
</figure>

Du kannst sehen, dass die beiden Module ähnliche Eingabeanforderungen in Bezug auf den Inhalt haben (eine Menge von Eingabedateien plus einige Metadaten), aber sehr unterschiedliche Erwartungen dafür, wie dieser Inhalt verpackt ist.
Wenn wir die Versionsdatei vorerst ignorieren, ist auch ihre Hauptausgabe gleichwertig (eine verkettete Datei), außer dass CAT_CAT auch die Metamap zusammen mit der Ausgabedatei ausgibt.

Die Verpackungsunterschiede werden relativ einfach zu handhaben sein, wie du gleich sehen wirst.
Um den Metamap-Teil zu verstehen, müssen wir dir jedoch zusätzlichen Kontext geben.

### 2.2. Metamaps verstehen

Wir haben dir gerade gesagt, dass das CAT_CAT-Modul eine Metadaten-Map als Teil seines Eingabe-Tupels erwartet.
Lass uns ein paar Minuten nehmen, um genauer zu betrachten, was das ist.

Die **Metadaten-Map**, oft kurz als **Metamap** bezeichnet, ist eine Map im Groovy-Stil, die Informationen über Dateneinheiten enthält.
Im Kontext von Nextflow-Pipelines können Dateneinheiten alles sein, was du möchtest: einzelne Proben, Chargen von Proben oder ganze Datensätze.

Per Konvention wird eine nf-core-Metamap `meta` genannt und enthält das erforderliche Feld `id`, das zum Benennen von Ausgaben und Verfolgen von Dateneinheiten verwendet wird.

Zum Beispiel könnte eine typische Metadaten-Map so aussehen:

```groovy title="Beispiel einer Metamap auf Probenebene"
[id: 'sample1', single_end: false, strandedness: 'forward']
```

Oder in einem Fall, in dem die Metadaten auf Chargenebene angehängt sind:

```groovy title="Beispiel einer Metamap auf Chargenebene"
[id: 'batch1', date: '25.10.01']
```

Lass uns das nun in den Kontext des `CAT_CAT`-Prozesses setzen, der erwartet, dass die Eingabedateien in ein Tupel mit einer Metamap verpackt sind, und die Metamap auch als Teil des Ausgabe-Tupels ausgibt.

```groovy title="modules/nf-core/cat/cat/main.nf (Auszug)" linenums="1" hl_lines="2 5"
input:
tuple val(meta), path(files_in)

output:
tuple val(meta), path("${prefix}"), emit: file_out
```

Dadurch durchläuft jede Dateneinheit die Pipeline mit den relevanten Metadaten versehen.
Nachfolgende Prozesse können dann auch problemlos auf diese Metadaten zugreifen.

Erinnerst du dich daran, dass wir dir gesagt haben, dass die von `CAT_CAT` ausgegebene Datei basierend auf einem Identifikator benannt wird, der Teil der Metadaten ist?
Das ist der relevante Code:

```groovy title="modules/nf-core/cat/cat/main.nf (Auszug)" linenums="35"
prefix   = task.ext.prefix ?: "${meta.id}${getFileSuffix(file_list[0])}"
```

Das bedeutet ungefähr Folgendes: Wenn ein `prefix` über das externe Task-Parametersystem (`task.ext`) bereitgestellt wird, verwende diesen, um die Ausgabedatei zu benennen; andernfalls erstelle einen mit `${meta.id}`, was dem `id`-Feld in der Metamap entspricht.

Du kannst dir den eingehenden Eingabekanal in dieses Modul mit Inhalten wie diesem vorstellen:

```groovy title="Beispiel für Eingabekanal-Inhalte"
ch_input = [[[id: 'batch1', date: '25.10.01'], ['file1A.txt', 'file1B.txt']],
            [[id: 'batch2', date: '25.10.26'], ['file2A.txt', 'file2B.txt']],
            [[id: 'batch3', date: '25.11.14'], ['file3A.txt', 'file3B.txt']]]
```

Dann kommen die Ausgabekanal-Inhalte so heraus:

```groovy title="Beispiel für Ausgabekanal-Inhalte"
ch_input = [[[id: 'batch1', date: '25.10.01'], 'batch1.txt'],
            [[id: 'batch2', date: '25.10.26'], 'batch2.txt'],
            [[id: 'batch3', date: '25.11.14'], 'batch3.txt']]
```

Wie bereits erwähnt, ist das `tuple val(meta), path(files_in)` Eingabe-Setup ein Standardmuster, das in allen nf-core-Modulen verwendet wird.

Hoffentlich kannst du allmählich sehen, wie nützlich das sein kann.
Es ermöglicht dir nicht nur, Ausgaben basierend auf Metadaten zu benennen, sondern du kannst auch Dinge tun wie verschiedene Parameterwerte anzuwenden, und in Kombination mit bestimmten Operatoren kannst du sogar Daten gruppieren, sortieren oder herausfiltern, während sie durch die Pipeline fließen.

!!! note "Mehr über Metadaten erfahren"

    Für eine umfassende Einführung in die Arbeit mit Metadaten in Nextflow-Workflows, einschließlich wie man Metadaten aus Samplesheets liest und sie zur Anpassung der Verarbeitung verwendet, siehe die Side Quest [Metadaten in Workflows](../side_quests/metadata).

### 2.3. Vorzunehmende Änderungen zusammenfassen

Basierend auf dem, was wir überprüft haben, sind dies die wichtigsten Änderungen, die wir an unserer Pipeline vornehmen müssen, um das `cat/cat`-Modul zu nutzen:

- Eine Metamap erstellen, die den Chargennamen enthält;
- Die Metamap in ein Tupel mit der Menge der zu verkettenden Eingabedateien verpacken (kommend aus `convertToUpper`);
- Den Aufruf von `collectGreetings()` zu `CAT_CAT` ändern;
- Die Ausgabedatei aus dem vom `CAT_CAT`-Prozess produzierten Tupel extrahieren, bevor sie an `cowpy` übergeben wird.

Das sollte reichen! Jetzt, da wir einen Plan haben, sind wir bereit loszulegen.

### Zusammenfassung

Du weißt jetzt, wie du die Eingabe- und Ausgabeschnittstelle eines neuen Moduls bewertest, um seine Anforderungen zu identifizieren, und du hast gelernt, wie Metamaps von nf-core-Pipelines verwendet werden, um Metadaten eng mit den Daten verbunden zu halten, während sie durch eine Pipeline fließen.

### Wie geht's weiter?

Das neue Modul in einen Workflow integrieren.

---

## 3. CAT_CAT in den `hello.nf` Workflow integrieren

Jetzt, da du alles über Metamaps weißt (oder zumindest genug für die Zwecke dieses Kurses), ist es an der Zeit, die oben skizzierten Änderungen tatsächlich zu implementieren.

Der Klarheit halber werden wir dies aufschlüsseln und jeden Schritt separat behandeln.

!!! note

    Alle unten gezeigten Änderungen werden an der Workflow-Logik im `main`-Block in der `core-hello/workflows/hello.nf` Workflow-Datei vorgenommen.

### 3.1. Eine Metadaten-Map erstellen

Zuerst müssen wir eine Metadaten-Map für `CAT_CAT` erstellen, wobei wir bedenken, dass nf-core-Module mindestens ein `id`-Feld in der Metamap benötigen.

Da wir keine anderen Metadaten benötigen, können wir es einfach halten und so etwas verwenden:

```groovy title="Syntaxbeispiel"
def cat_meta = [id: 'test']
```

Außer dass wir den `id`-Wert nicht fest codieren wollen; wir wollen den Wert des `params.batch`-Parameters verwenden.
Der Code wird also:

```groovy title="Syntaxbeispiel"
def cat_meta = [id: params.batch]
```

Ja, es ist buchstäblich so einfach, eine grundlegende Metamap zu erstellen.

Lass uns diese Zeilen nach dem `convertToUpper`-Aufruf hinzufügen und den `collectGreetings`-Aufruf entfernen:

=== "Nachher"

    ```groovy title="core-hello/workflows/hello.nf" linenums="26" hl_lines="7-8"
        // Eine Begrüßung ausgeben
        sayHello(ch_samplesheet)

        // Die Begrüßung in Großbuchstaben umwandeln
        convertToUpper(sayHello.out)

        // Metadaten-Map mit Chargennamen als ID erstellen
        def cat_meta = [ id: params.batch ]

        // ASCII-Kunst der Begrüßungen mit cowpy generieren
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

        // ASCII-Kunst der Begrüßungen mit cowpy generieren
        cowpy(collectGreetings.out.outfile, params.character)
    ```

Dies erstellt eine einfache Metadaten-Map, bei der die `id` auf unseren Chargennamen gesetzt ist (was bei Verwendung des Testprofils `test` sein wird).

### 3.2. Einen Kanal mit Metadaten-Tupeln erstellen

Als Nächstes transformieren wir den Kanal von Dateien in einen Kanal von Tupeln, die Metadaten und Dateien enthalten:

=== "Nachher"

    ```groovy title="core-hello/workflows/hello.nf" linenums="26" hl_lines="10-11"
        // Eine Begrüßung ausgeben
        sayHello(ch_samplesheet)

        // Die Begrüßung in Großbuchstaben umwandeln
        convertToUpper(sayHello.out)

        // Metadaten-Map mit Chargennamen als ID erstellen
        def cat_meta = [ id: params.batch ]

        // Kanal mit Metadaten und Dateien im Tupel-Format erstellen
        ch_for_cat = convertToUpper.out.collect().map { files -> tuple(cat_meta, files) }

        // ASCII-Kunst der Begrüßungen mit cowpy generieren
        cowpy(collectGreetings.out.outfile, params.character)
    ```

=== "Vorher"

    ```groovy title="core-hello/workflows/hello.nf" linenums="26"
        // Eine Begrüßung ausgeben
        sayHello(ch_samplesheet)

        // Die Begrüßung in Großbuchstaben umwandeln
        convertToUpper(sayHello.out)

        // Metadaten-Map mit Chargennamen als ID erstellen
        def cat_meta = [ id: params.batch ]

        // ASCII-Kunst der Begrüßungen mit cowpy generieren
        cowpy(collectGreetings.out.outfile, params.character)
    ```

Die Zeile, die wir hinzugefügt haben, erreicht zwei Dinge:

- `.collect()` sammelt alle Dateien aus der `convertToUpper`-Ausgabe in eine einzelne Liste
- `.map { files -> tuple(cat_meta, files) }` erstellt ein Tupel `[metadata, files]` im Format, das `CAT_CAT` erwartet

Das ist alles, was wir tun müssen, um das Eingabe-Tupel für `CAT_CAT` einzurichten.

### 3.3. Das CAT_CAT-Modul aufrufen

Rufe nun `CAT_CAT` auf dem neu erstellten Kanal auf:

=== "Nachher"

    ```groovy title="core-hello/workflows/hello.nf" linenums="26" hl_lines="13-14"
        // Eine Begrüßung ausgeben
        sayHello(ch_samplesheet)

        // Die Begrüßung in Großbuchstaben umwandeln
        convertToUpper(sayHello.out)

        // Metadaten-Map mit Chargennamen als ID erstellen
        def cat_meta = [ id: params.batch ]

        // Kanal mit Metadaten und Dateien im Tupel-Format erstellen
        ch_for_cat = convertToUpper.out.collect().map { files -> tuple(cat_meta, files) }

        // Dateien mit dem nf-core cat/cat-Modul verketten
        CAT_CAT(ch_for_cat)

        // ASCII-Kunst der Begrüßungen mit cowpy generieren
        cowpy(collectGreetings.out.outfile, params.character)
    ```

=== "Vorher"

    ```groovy title="core-hello/workflows/hello.nf" linenums="26"
        // Eine Begrüßung ausgeben
        sayHello(ch_samplesheet)

        // Die Begrüßung in Großbuchstaben umwandeln
        convertToUpper(sayHello.out)

        // Metadaten-Map mit Chargennamen als ID erstellen
        def cat_meta = [ id: params.batch ]

        // Kanal mit Metadaten und Dateien im Tupel-Format erstellen
        ch_for_cat = convertToUpper.out.collect().map { files -> tuple(cat_meta, files) }

        // ASCII-Kunst der Begrüßungen mit cowpy generieren
        cowpy(collectGreetings.out.outfile, params.character)
    ```

Dies vervollständigt den kniffligsten Teil dieser Ersetzung, aber wir sind noch nicht ganz fertig: wir müssen noch aktualisieren, wie wir die verkettete Ausgabe an den `cowpy`-Prozess übergeben.

### 3.4. Die Ausgabedatei aus dem Tupel für `cowpy` extrahieren

Zuvor hat der `collectGreetings`-Prozess einfach eine Datei produziert, die wir direkt an `cowpy` übergeben konnten.
Der `CAT_CAT`-Prozess produziert jedoch ein Tupel, das zusätzlich zur Ausgabedatei die Metamap enthält.

Da `cowpy` noch keine Metadaten-Tupel akzeptiert (das werden wir im nächsten Teil des Kurses beheben), müssen wir die Ausgabedatei aus dem von `CAT_CAT` produzierten Tupel extrahieren, bevor wir sie an `cowpy` übergeben:

=== "Nachher"

    ```groovy title="core-hello/workflows/hello.nf" linenums="26" hl_lines="16-17 20"
        // Eine Begrüßung ausgeben
        sayHello(ch_samplesheet)

        // Die Begrüßung in Großbuchstaben umwandeln
        convertToUpper(sayHello.out)

        // Metadaten-Map mit Chargennamen als ID erstellen
        def cat_meta = [ id: params.batch ]

        // Kanal mit Metadaten und Dateien im Tupel-Format erstellen
        ch_for_cat = convertToUpper.out.collect().map { files -> tuple(cat_meta, files) }

        // Die Begrüßungen verketten
        CAT_CAT(ch_for_cat)

        // Datei aus dem Tupel extrahieren, da cowpy noch keine Metadaten verwendet
        ch_for_cowpy = CAT_CAT.out.file_out.map{ meta, file -> file }

        // ASCII-Kunst der Begrüßungen mit cowpy generieren
        cowpy(ch_for_cowpy, params.character)
    ```

=== "Vorher"

    ```groovy title="core-hello/workflows/hello.nf" linenums="26" hl_lines="17"
        // Eine Begrüßung ausgeben
        sayHello(ch_samplesheet)

        // Die Begrüßung in Großbuchstaben umwandeln
        convertToUpper(sayHello.out)

        // Metadaten-Map mit Chargennamen als ID erstellen
        def cat_meta = [ id: params.batch ]

        // Kanal mit Metadaten und Dateien im Tupel-Format erstellen
        ch_for_cat = convertToUpper.out.collect().map { files -> tuple(cat_meta, files) }

        // Die Begrüßungen verketten
        CAT_CAT(ch_for_cat)

        // ASCII-Kunst der Begrüßungen mit cowpy generieren
        cowpy(collectGreetings.out.outfile, params.character)
    ```

Die `.map{ meta, file -> file }`-Operation extrahiert die Datei aus dem `[metadata, file]`-Tupel, das von `CAT_CAT` produziert wurde, in einen neuen Kanal, `ch_for_cowpy`.

Dann ist es nur noch eine Sache, `ch_for_cowpy` anstelle von `collectGreetings.out.outfile` in dieser letzten Zeile an `cowpy` zu übergeben.

!!! note

    Im nächsten Teil des Kurses werden wir `cowpy` aktualisieren, damit es direkt mit Metadaten-Tupeln funktioniert, sodass dieser Extraktionsschritt nicht mehr notwendig sein wird.

### 3.5. Den Workflow testen

Lass uns testen, ob der Workflow mit dem neu integrierten `cat/cat`-Modul funktioniert:

```bash
nextflow run . --outdir core-hello-results -profile test,docker --validate_params false
```

Dies sollte relativ schnell laufen.

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

Beachte, dass `CAT_CAT` jetzt in der Liste der Prozessausführungen anstelle von `collectGreetings` erscheint.

Und das war's! Wir verwenden jetzt ein robustes, von der Community kuratiertes Modul anstelle von eigenem Prototyp-Code für diesen Schritt in der Pipeline.

### Zusammenfassung

Du weißt jetzt, wie du:

- nf-core-Module findest und installierst
- Die Anforderungen eines nf-core-Moduls bewertest
- Eine einfache Metadaten-Map zur Verwendung mit einem nf-core-Modul erstellst
- Ein nf-core-Modul in deinen Workflow integrierst

### Wie geht's weiter?

Lerne, deine lokalen Module anzupassen, um nf-core-Konventionen zu folgen.
Wir zeigen dir auch, wie du neue nf-core-Module aus einer Vorlage mit den nf-core-Werkzeugen erstellst.
