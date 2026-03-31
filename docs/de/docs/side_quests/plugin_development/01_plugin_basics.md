# Teil 1: Plugin-Grundlagen

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } KI-gestützte Übersetzung - [mehr erfahren & Verbesserungen vorschlagen](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

In diesem Abschnitt lernst du, wie Plugins Nextflow erweitern, und probierst drei verschiedene Plugins aus, um sie in Aktion zu sehen.

---

## 1. Wie Plugins funktionieren

Plugins erweitern Nextflow durch verschiedene Arten von Erweiterungen:

| Erweiterungstyp    | Was es tut                                              | Beispiel                     |
| ------------------ | ------------------------------------------------------- | ---------------------------- |
| Functions          | Fügt benutzerdefinierte Funktionen hinzu, die aus Workflows aufgerufen werden können | `samplesheetToList()`        |
| Workflow monitors  | Reagiert auf Ereignisse wie den Abschluss von Aufgaben  | Benutzerdefiniertes Logging, Slack-Benachrichtigungen |
| Executors          | Fügt Backends für die Aufgabenausführung hinzu          | AWS Batch, Kubernetes        |
| Filesystems        | Fügt Speicher-Backends hinzu                            | S3, Azure Blob               |

Functions und Workflow-Monitore (im Nextflow-API als „trace observers" bezeichnet) sind die häufigsten Typen für Plugin-Entwickler\*innen.
Executors und Filesystems werden typischerweise von Plattformanbietern erstellt.

Die nächsten Übungen zeigen dir Function-Plugins und ein Observer-Plugin, damit du beide Typen in Aktion sehen kannst.

---

## 2. Function-Plugins verwenden

Function-Plugins stellen aufrufbare Funktionen bereit, die du in deine Workflows importierst.
Du probierst zwei aus: nf-hello (ein einfaches Beispiel) und nf-schema (ein weit verbreitetes Plugin aus der Praxis).
Beide Übungen ändern dieselbe `hello.nf`-Pipeline, damit du sehen kannst, wie Plugins einen bestehenden Workflow erweitern.

### 2.1. nf-hello: handgeschriebenen Code ersetzen

Das [nf-hello](https://github.com/nextflow-io/nf-hello)-Plugin stellt eine `randomString`-Funktion bereit, die zufällige Zeichenketten generiert.
Die Pipeline definiert bereits eine eigene inline Version dieser Funktion, die du durch die Version aus dem Plugin ersetzen wirst.

#### 2.1.1. Den Ausgangspunkt ansehen

Schau dir die Pipeline an:

```bash
cat hello.nf
```

```groovy title="Output"
#!/usr/bin/env nextflow

params.input = 'greetings.csv'

/**
 * Einen zufälligen alphanumerischen String generieren
 */
def randomString(int length) {
    def chars = ('a'..'z') + ('A'..'Z') + ('0'..'9')
    def random = new Random()
    return (1..length).collect { chars[random.nextInt(chars.size())] }.join()
}

process SAY_HELLO {
    input:
        val greeting
    output:
        stdout
    script:
    """
    echo '$greeting'
    """
}

workflow {
    greeting_ch = channel.fromPath(params.input)
        .splitCsv(header: true)
        .map { row -> "${row.greeting}_${randomString(8)}" }
    SAY_HELLO(greeting_ch)
    SAY_HELLO.out.view { result -> "Output: ${result.trim()}" }
}
```

Die Pipeline definiert ihre eigene `randomString`-Funktion inline und verwendet sie, um jeder Begrüßung eine zufällige ID anzuhängen.

Ausführen:

```bash
nextflow run hello.nf
```

```console title="Output"
Output: Hello_aBcDeFgH
Output: Bonjour_xYzWvUtS
Output: Holà_qRsPdMnK
Output: Ciao_jLhGfEcB
Output: Hallo_tNwOiAuR
```

Die Reihenfolge der Ausgabe und die zufälligen Zeichenketten werden bei dir anders sein. Wenn du das Skript erneut ausführst, erhältst du andere zufällige Begrüßungen.

#### 2.1.2. Das Plugin konfigurieren

Ersetze die inline Funktion durch eine aus dem Plugin. Füge folgendes zu deiner `nextflow.config` hinzu:

```groovy title="nextflow.config"
// Konfiguration für Plugin-Entwicklungsübungen
plugins {
    id 'nf-hello@0.5.0'
}
```

Plugins werden in `nextflow.config` mit dem `plugins {}`-Block deklariert.
Nextflow lädt sie automatisch aus der [Nextflow Plugin Registry](https://registry.nextflow.io/) herunter, einem zentralen Repository für Community- und offizielle Plugins.

#### 2.1.3. Die Plugin-Funktion verwenden

Ersetze die inline `randomString`-Funktion durch die Plugin-Version:

=== "Danach"

    ```groovy title="hello.nf" hl_lines="3"
    #!/usr/bin/env nextflow

    include { randomString } from 'plugin/nf-hello'

    params.input = 'greetings.csv'

    process SAY_HELLO {
        input:
            val greeting
        output:
            stdout
        script:
        """
        echo '$greeting'
        """
    }

    workflow {
        greeting_ch = channel.fromPath(params.input)
            .splitCsv(header: true)
            .map { row -> "${row.greeting}_${randomString(8)}" }
        SAY_HELLO(greeting_ch)
        SAY_HELLO.out.view { result -> "Output: ${result.trim()}" }
    }
    ```

=== "Vorher"

    ```groovy title="hello.nf" hl_lines="5-12"
    #!/usr/bin/env nextflow

    params.input = 'greetings.csv'

    /**
     * Einen zufälligen alphanumerischen String generieren
     */
    def randomString(int length) {
        def chars = ('a'..'z') + ('A'..'Z') + ('0'..'9')
        def random = new Random()
        return (1..length).collect { chars[random.nextInt(chars.size())] }.join()
    }

    process SAY_HELLO {
        input:
            val greeting
        output:
            stdout
        script:
        """
        echo '$greeting'
        """
    }

    workflow {
        greeting_ch = channel.fromPath(params.input)
            .splitCsv(header: true)
            .map { row -> "${row.greeting}_${randomString(8)}" }
        SAY_HELLO(greeting_ch)
        SAY_HELLO.out.view { result -> "Output: ${result.trim()}" }
    }
    ```

Die `include`-Anweisung importiert `randomString` aus einer Bibliothek, die erprobt, getestet und von einer größeren Gruppe von Mitwirkenden gepflegt wird, die Fehler erkennen und beheben können.
Anstatt dass jede Pipeline ihre eigene Kopie der Funktion pflegt, erhält jede Pipeline, die das Plugin verwendet, dieselbe geprüfte Implementierung.
Das reduziert doppelten Code und den damit verbundenen Wartungsaufwand.
Die Syntax `#!groovy include { function } from 'plugin/plugin-id'` ist dasselbe `include`, das für Nextflow-Module verwendet wird, mit dem Präfix `plugin/`.
Den [Quellcode von `randomString`](https://github.com/nextflow-io/nf-hello/blob/e67bddebfa589c7ae51f41bf780c92068dc09e93/plugins/nf-hello/src/main/nextflow/hello/HelloExtension.groovy#L110) findest du im nf-hello-Repository auf GitHub.

#### 2.1.4. Ausführen

```bash
nextflow run hello.nf
```

```console title="Output"
Pipeline is starting! 🚀
Output: Hello_yqvtclcc
Output: Bonjour_vwwpyzcs
Output: Holà_wrghmgab
Output: Ciao_noniajuy
Output: Hallo_tvrtuxtp
Pipeline complete! 👋
```

(Deine zufälligen Zeichenketten werden anders sein.)

Die Ausgabe hat weiterhin zufällige Suffixe, aber jetzt kommt `randomString` aus dem nf-hello-Plugin statt aus inline Code.
Die Meldungen „Pipeline is starting!" und „Pipeline complete!" sind neu.
Sie stammen von der Observer-Komponente des Plugins, die du in Teil 5 genauer kennenlernen wirst.

Nextflow lädt Plugins beim ersten Einsatz automatisch herunter. Jede Pipeline, die `nf-hello@0.5.0` deklariert, erhält so dieselbe getestete `randomString`-Funktion, ohne Code zwischen Projekten kopieren zu müssen.

Du hast jetzt die drei Schritte zur Verwendung eines Function-Plugins gesehen: es in `nextflow.config` deklarieren, die Funktion mit `include` importieren und sie im Workflow aufrufen.
Die nächste Übung wendet dieselben Schritte auf ein Plugin aus der Praxis an.

### 2.2. nf-schema: validiertes CSV-Parsing

Das [nf-schema](https://github.com/nextflow-io/nf-schema)-Plugin ist eines der am häufigsten verwendeten Nextflow-Plugins.
Es stellt `samplesheetToList` bereit, eine Funktion, die CSV/TSV-Dateien anhand eines JSON-Schemas parst, das die erwarteten Spalten und Typen definiert.

Die Pipeline liest `greetings.csv` derzeit mit `splitCsv` und einem manuellen `map`, aber nf-schema kann das durch validiertes, schemabasiertes Parsing ersetzen.
Eine JSON-Schema-Datei (`greetings_schema.json`) ist bereits im Übungsverzeichnis vorhanden.

??? info "Was ist ein Schema?"

    Ein Schema ist eine formale Beschreibung, wie gültige Daten aussehen sollen.
    Es definiert Dinge wie welche Spalten erwartet werden, welchen Typ jeder Wert haben soll (String, Zahl usw.) und welche Felder erforderlich sind.

    Stell es dir als Vertrag vor: Wenn die Eingabedaten nicht zum Schema passen, kann das Tool das Problem frühzeitig erkennen, anstatt es später in der Pipeline zu verwirrenden Fehlern kommen zu lassen.

#### 2.2.1. Das Schema ansehen

```bash
cat greetings_schema.json
```

```json title="Output"
{
  "$schema": "https://json-schema.org/draft/2020-12/schema",
  "type": "array",
  "items": {
    "type": "object",
    "properties": {
      "greeting": {
        "type": "string",
        "description": "The greeting text"
      },
      "language": {
        "type": "string",
        "description": "The language of the greeting"
      }
    },
    "required": ["greeting"]
  }
}
```

Das Schema definiert zwei Spalten (`greeting` und `language`) und markiert `greeting` als erforderlich.
Wenn jemand eine CSV ohne die Spalte `greeting` übergibt, erkennt nf-schema den Fehler, bevor die Pipeline ausgeführt wird.

#### 2.2.2. nf-schema zur Konfiguration hinzufügen

Aktualisiere `nextflow.config`, um beide Plugins einzuschließen:

=== "Danach"

    ```groovy title="nextflow.config" hl_lines="3"
    plugins {
        id 'nf-hello@0.5.0'
        id 'nf-schema@2.6.1'
    }
    ```

=== "Vorher"

    ```groovy title="nextflow.config"
    plugins {
        id 'nf-hello@0.5.0'
    }
    ```

#### 2.2.3. hello.nf aktualisieren, um samplesheetToList zu verwenden

Ersetze die `splitCsv`-Eingabe durch `samplesheetToList`:

=== "Danach"

    ```groovy title="hello.nf" hl_lines="4 20 21 22"
    #!/usr/bin/env nextflow

    include { randomString } from 'plugin/nf-hello'
    include { samplesheetToList } from 'plugin/nf-schema'

    params.input = 'greetings.csv'

    process SAY_HELLO {
        input:
            val greeting
        output:
            stdout
        script:
        """
        echo '$greeting'
        """
    }

    workflow {
        def samplesheet_list = samplesheetToList(params.input, 'greetings_schema.json')
        greeting_ch = Channel.fromList(samplesheet_list)
            .map { row -> "${row[0]}_${randomString(8)}" }
        SAY_HELLO(greeting_ch)
        SAY_HELLO.out.view { result -> "Output: ${result.trim()}" }
    }
    ```

=== "Vorher"

    ```groovy title="hello.nf" hl_lines="19 20 21"
    #!/usr/bin/env nextflow

    include { randomString } from 'plugin/nf-hello'

    params.input = 'greetings.csv'

    process SAY_HELLO {
        input:
            val greeting
        output:
            stdout
        script:
        """
        echo '$greeting'
        """
    }

    workflow {
        greeting_ch = channel.fromPath(params.input)
            .splitCsv(header: true)
            .map { row -> "${row.greeting}_${randomString(8)}" }
        SAY_HELLO(greeting_ch)
        SAY_HELLO.out.view { result -> "Output: ${result.trim()}" }
    }
    ```

Der benutzerdefinierte `splitCsv`- und `map`-Parsing-Code wird durch `samplesheetToList` ersetzt, eine erprobte und getestete Funktion, die das Samplesheet außerdem gegen das Schema validiert, bevor die Pipeline ausgeführt wird.
Das reduziert den Wartungsaufwand für handgeschriebene Parsing-Logik und verbessert gleichzeitig die Erfahrung für Pipeline-Nutzer\*innen, die klare Fehlermeldungen erhalten, wenn ihre Eingabe nicht dem erwarteten Format entspricht.
Jede Zeile wird zu einer Liste von Werten in Spaltenreihenfolge, also ist `row[0]` die Begrüßung und `row[1]` die Sprache.

#### 2.2.4. Ausführen

```bash
nextflow run hello.nf
```

```console title="Output"
Pipeline is starting! 🚀
Output: Hello_diozjdwm
Output: Bonjour_speathmm
Output: Holà_dllxnzap
Output: Ciao_wzueddzc
Output: Hallo_hsxwrjbh
Pipeline complete! 👋
```

(Deine zufälligen Zeichenketten werden anders sein.)

Die Ausgabe ist dieselbe, aber jetzt validiert das Schema die CSV-Struktur, bevor die Pipeline ausgeführt wird.
In echten Pipelines mit komplexen Samplesheets und vielen Spalten verhindert diese Art der Validierung Fehler, die manuelles `splitCsv` + `map` übersehen würde.

#### 2.2.5. Validierung in Aktion sehen

Um zu sehen, was die Schema-Validierung erkennt, füge Fehler in `greetings.csv` ein.

Benenne die erforderliche Spalte `greeting` in `message` um:

```csv title="greetings.csv" hl_lines="1"
message,language
Hello,English
Bonjour,French
Holà,Spanish
Ciao,Italian
Hallo,German
```

Führe die Pipeline aus:

```bash
nextflow run hello.nf
```

```console title="Output"
ERROR ~ Validation of samplesheet failed!

The following errors have been detected in greetings.csv:

-> Entry 1: Missing required field(s): greeting
-> Entry 2: Missing required field(s): greeting
-> Entry 3: Missing required field(s): greeting
-> Entry 4: Missing required field(s): greeting
-> Entry 5: Missing required field(s): greeting
```

Die Pipeline verweigert die Ausführung, weil das Schema eine Spalte `greeting` erfordert und keine findet.

Stelle jetzt die erforderliche Spalte wieder her, aber benenne die optionale Spalte `language` in `lang` um:

```csv title="greetings.csv" hl_lines="1"
greeting,lang
Hello,English
Bonjour,French
Holà,Spanish
Ciao,Italian
Hallo,German
```

```bash
nextflow run hello.nf
```

Diesmal wird die Pipeline ausgeführt, gibt aber eine Warnung aus:

```console title="Output (partial)"
WARN: Found the following unidentified headers in greetings.csv:
	- lang
```

Erforderliche Spalten verursachen harte Fehler; optionale Spalten verursachen Warnungen.
Das ist die Art von frühzeitigem Feedback, das in echten Pipelines mit Dutzenden von Spalten Debugging-Zeit spart.

#### 2.2.6. Validierungsverhalten konfigurieren

Die Warnung über `lang` ist nützlich, aber du kannst ihren Schweregrad über die Konfiguration steuern.
Plugins können eigene Konfigurationsbereiche definieren, die ihr Verhalten steuern.
Das nf-schema-Plugin enthält den Konfigurationsbereich `validation`; durch Ändern der Einstellungen hier kannst du das Verhalten von nf-schema anpassen.

Füge einen `validation`-Block zu `nextflow.config` hinzu, damit nicht erkannte Spaltenüberschriften einen Fehler statt einer Warnung verursachen:

=== "Danach"

    ```groovy title="nextflow.config" hl_lines="6-10"
    plugins {
        id 'nf-hello@0.5.0'
        id 'nf-schema@2.6.1'
    }

    validation {
        logging {
            unrecognisedHeaders = "error"
        }
    }
    ```

=== "Vorher"

    ```groovy title="nextflow.config"
    plugins {
        id 'nf-hello@0.5.0'
        id 'nf-schema@2.6.1'
    }
    ```

Führe die Pipeline erneut aus, während die Spalte `lang` noch vorhanden ist:

```bash
nextflow run hello.nf
```

```console title="Output (partial)"
Found the following unidentified headers in greetings.csv:
	- lang
 -- Check script 'hello.nf' at line: 20 or see '.nextflow.log' file for more details
```

Die Pipeline schlägt jetzt fehl, anstatt nur zu warnen.
Der Pipeline-Code hat sich nicht geändert; nur die Konfiguration.

Stelle `greetings.csv` in den ursprünglichen Zustand zurück und entferne den `validation`-Block, bevor du weitermachst:

```csv title="greetings.csv"
greeting,language
Hello,English
Bonjour,French
Holà,Spanish
Ciao,Italian
Hallo,German
```

```groovy title="nextflow.config"
plugins {
    id 'nf-hello@0.5.0'
    id 'nf-schema@2.6.1'
}
```

Sowohl nf-hello als auch nf-schema sind Function-Plugins: Sie stellen Funktionen bereit, die du mit `include` importierst und in deinem Workflow-Code aufrufst.
Die nächste Übung zeigt einen anderen Plugin-Typ, der ganz ohne `include`-Anweisungen auskommt.

---

## 3. Ein Observer-Plugin verwenden: nf-co2footprint

Nicht alle Plugins stellen Funktionen zum Importieren bereit.
Das [nf-co2footprint](https://github.com/nextflow-io/nf-co2footprint)-Plugin verwendet einen **trace observer**, um die Ressourcennutzung deiner Pipeline zu überwachen und ihren CO₂-Fußabdruck zu schätzen.
Du musst keinen Pipeline-Code ändern; füge es einfach zur Konfiguration hinzu.

### 3.1. nf-co2footprint zur Konfiguration hinzufügen

Aktualisiere `nextflow.config`:

=== "Danach"

    ```groovy title="nextflow.config" hl_lines="4"
    plugins {
        id 'nf-hello@0.5.0'
        id 'nf-schema@2.6.1'
        id 'nf-co2footprint@1.2.0'
    }
    ```

=== "Vorher"

    ```groovy title="nextflow.config"
    plugins {
        id 'nf-hello@0.5.0'
        id 'nf-schema@2.6.1'
    }
    ```

### 3.2. Die Pipeline ausführen

```bash
nextflow run hello.nf
```

Das Plugin gibt während der Ausführung mehrere INFO- und WARN-Meldungen aus.
Diese sind für ein kleines Beispiel auf einem lokalen Rechner normal:

```console title="Output (partial)"
nf-co2footprint plugin  ~  version 1.2.0
WARN - [nf-co2footprint] Target zone null not found. Attempting to retrieve carbon intensity for fallback zone GLOBAL.
INFO - [nf-co2footprint] Using fallback carbon intensity from GLOBAL from CI table: 480.0 gCO₂eq/kWh.
WARN - [nf-co2footprint] Executor 'null' not mapped.
WARN - [nf-co2footprint] Fallback to: `machineType = null`, `pue = 1.0`. ...
...
WARN - [nf-co2footprint] No CPU model detected. Using default CPU power draw value (11.41 W).
WARN - [nf-co2footprint] 🔁 Requested memory is null for task 2. Using maximum consumed memory/`peak_rss` (0 GB) for CO₂e footprint computation.
```

Die Warnungen über Zone, Executor, CPU-Modell und Arbeitsspeicher erscheinen, weil das Plugin die vollständigen Hardware-Details einer lokalen Trainingsumgebung nicht erkennen kann.
In einer Produktionsumgebung (z. B. einem HPC-Cluster oder in der Cloud) wären diese Werte verfügbar und die Schätzungen genauer.

Am Ende sollte eine Zeile wie diese erscheinen:

```console title="Output (partial)"
🌱 The workflow run used 126.76 uWh of electricity, resulting in the release of 60.84 ug of CO₂ equivalents into the atmosphere.
```

(Deine Zahlen werden anders sein.)

### 3.3. Den Bericht ansehen

Das Plugin erzeugt Ausgabedateien in deinem Arbeitsverzeichnis:

```bash
ls co2footprint_*
```

```console title="Output"
co2footprint_report_<timestamp>.html
co2footprint_summary_<timestamp>.txt
co2footprint_trace_<timestamp>.txt
```

Schau dir die Zusammenfassung an:

```bash
cat co2footprint_summary_*.txt
```

```console title="Output"
Total CO₂e footprint measures of this workflow run (including cached tasks):
  CO₂e emissions: 60.84 ug
  Energy consumption: 126.76 uWh
  CO₂e emissions (market): -

Which equals:
  - 3.48E-7 km travelled by car
  - It takes one tree 0.17s to sequester the equivalent amount of CO₂ from the atmosphere
  - 1.22E-7 % of a flight from Paris to London
```

(Deine Zahlen werden anders sein.)

Der erste Abschnitt zeigt die rohen Energie- und Emissionswerte.
Der Abschnitt „Which equals" setzt diese Zahlen in Perspektive, indem er sie in vertraute Äquivalente umrechnet.
Die Zusammenfassung enthält außerdem einen Abschnitt mit den Konfigurationsoptionen des Plugins und einen Verweis auf das [Green Algorithms](https://doi.org/10.1002/advs.202100707)-Forschungspapier, auf dem die Berechnungsmethode basiert.

### 3.4. Das Plugin konfigurieren

Die Warnung „Target zone null" aus Abschnitt 3.2 erschien, weil für das Plugin kein Standort konfiguriert war.
Das nf-co2footprint-Plugin definiert einen Konfigurationsbereich `co2footprint`, in dem du deinen geografischen Standort festlegen kannst.

Füge einen `co2footprint`-Block zu `nextflow.config` hinzu:

=== "Danach"

    ```groovy title="nextflow.config" hl_lines="7-9"
    plugins {
        id 'nf-hello@0.5.0'
        id 'nf-schema@2.6.1'
        id 'nf-co2footprint@1.2.0'
    }

    co2footprint {
        location = 'GB'
    }
    ```

=== "Vorher"

    ```groovy title="nextflow.config"
    plugins {
        id 'nf-hello@0.5.0'
        id 'nf-schema@2.6.1'
        id 'nf-co2footprint@1.2.0'
    }
    ```

!!! tip "Tipp"

    Verwende deinen eigenen Ländercode, wenn du möchtest (z. B. `'US'`, `'DE'`, `'FR'`).

Führe die Pipeline aus:

```bash
nextflow run hello.nf
```

```console title="Output (partial)"
INFO - [nf-co2footprint] Using fallback carbon intensity from GB from CI table: 163.92 gCO₂eq/kWh.
```

Die Zonenwarnung ist verschwunden.
Das Plugin verwendet jetzt die GB-spezifische CO₂-Intensität (163.92 gCO₂eq/kWh) statt des globalen Fallback-Werts (480.0 gCO₂eq/kWh).

!!! note "Hinweis"

    Möglicherweise siehst du auch eine Meldung `WARN: Unrecognized config option 'co2footprint.location'`.
    Diese ist kosmetischer Natur und kann bedenkenlos ignoriert werden; das Plugin liest den Wert trotzdem korrekt.

In Teil 6 wirst du einen Konfigurationsbereich für dein eigenes Plugin erstellen.

Dieses Plugin funktioniert vollständig über den Observer-Mechanismus: Es hängt sich in Workflow-Lifecycle-Ereignisse ein, um Ressourcenmetriken zu sammeln und seinen Bericht zu generieren, wenn die Pipeline abgeschlossen ist.

Du hast jetzt Function-Plugins (importiert mit `include`) und ein Observer-Plugin (allein durch die Konfiguration aktiviert) ausprobiert.
Das sind die zwei häufigsten Erweiterungstypen, aber wie die Tabelle in Abschnitt 1 zeigt, können Plugins auch Executors und Filesystems hinzufügen.

---

## 4. Plugins entdecken

Die [Nextflow Plugin Registry](https://registry.nextflow.io/) ist die zentrale Anlaufstelle, um verfügbare Plugins zu finden.

![Die nf-hello-Plugin-Seite auf registry.nextflow.io](img/plugin-registry-nf-hello.png)

Jede Plugin-Seite zeigt die Beschreibung, verfügbare Versionen, Installationsanweisungen und Links zur Dokumentation.

---

## 5. Für die Plugin-Entwicklung vorbereiten

Die folgenden Abschnitte (Teile 2–6) verwenden eine separate Pipeline-Datei, `greet.nf`, die nf-schema, aber nicht nf-hello oder nf-co2footprint benötigt.

Aktualisiere `nextflow.config`, um nur nf-schema zu behalten:

```groovy title="nextflow.config"
// Konfiguration für Plugin-Entwicklungsübungen
plugins {
    id 'nf-schema@2.6.1'
}
```

Entferne die co2footprint-Ausgabedateien:

```bash
rm -f co2footprint_*
```

Die Datei `hello.nf` behält deine Arbeit aus Teil 1 als Referenz; von nun an arbeitest du mit `greet.nf`.

---

## Fazit

Du hast drei verschiedene Plugins verwendet:

- **nf-hello**: Ein Function-Plugin, das `randomString` bereitstellt, importiert mit `include`
- **nf-schema**: Ein Function-Plugin, das `samplesheetToList` für schemavalidiertes CSV-Parsing bereitstellt
- **nf-co2footprint**: Ein Observer-Plugin, das die Ressourcennutzung automatisch überwacht, ohne `include`

Wichtige Muster:

- Plugins werden in `nextflow.config` mit `#!groovy plugins { id 'plugin-name@version' }` deklariert
- Function-Plugins erfordern `#!groovy include { function } from 'plugin/plugin-id'`
- Observer-Plugins funktionieren automatisch, sobald sie in der Konfiguration deklariert sind
- Plugins können Konfigurationsbereiche definieren (z. B. `#!groovy validation {}`, `#!groovy co2footprint {}`), um ihr Verhalten anzupassen
- Die [Nextflow Plugin Registry](https://registry.nextflow.io/) listet verfügbare Plugins auf

---

## Wie geht es weiter?

Die folgenden Abschnitte zeigen dir, wie du dein eigenes Plugin erstellst.
Wenn du kein Interesse an der Plugin-Entwicklung hast, kannst du hier aufhören oder direkt zur [Zusammenfassung](summary.md) springen.

[Weiter zu Teil 2 :material-arrow-right:](02_create_project.md){ .md-button .md-button--primary }
