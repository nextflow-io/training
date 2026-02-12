# Teil 2: Verarbeitung einzelner Proben

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } KI-gestützte Übersetzung - [mehr erfahren & Verbesserungen vorschlagen](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

In Teil 1 hast du die Befehle für {TOOL_A} und {TOOL_B} manuell in ihren jeweiligen Containern getestet.
Jetzt werden wir diese Befehle in einen Nextflow-Workflow einbinden.

## Aufgabe

In diesem Teil des Kurses entwickeln wir einen Workflow, der Folgendes tut:

1. {PROCESS_1_DESCRIPTION}
2. {PROCESS_2_DESCRIPTION}

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/{DOMAIN_DIR}/img/{DIAGRAM_1}.svg"
</figure>

Dies entspricht den Schritten aus Teil 1, in denen du diese Befehle manuell in ihren Containern ausgeführt hast.

Als Ausgangspunkt stellen wir dir eine Workflow-Datei `{DOMAIN_DIR}.nf` zur Verfügung, die die Hauptbestandteile des Workflows skizziert, sowie zwei Moduldateien {TOOL_A_MODULE}.nf und {TOOL_B_MODULE}.nf, die die Struktur der Module umreißen.
Diese Dateien sind nicht funktionsfähig; ihr Zweck besteht lediglich darin, als Gerüst zu dienen, das du mit den interessanten Teilen des Codes ausfüllen kannst.

## Lektionsplan

Um den Entwicklungsprozess lehrreicher zu gestalten, haben wir ihn in {N} Schritte unterteilt:

1. **Schreibe einen einstufigen Workflow, der {TOOL_A_ACTION} ausführt.**
   Dies umfasst das Erstellen eines Moduls, das Importieren und das Aufrufen im Workflow.
2. **Füge einen zweiten Prozess hinzu, um {TOOL_B_ACTION} auszuführen.**
   Dies führt die Verkettung von Prozessausgaben mit Eingaben ein und behandelt zusätzliche Dateien.
3. **Passe den Workflow an, um ihn auf mehreren Proben auszuführen.**
   Dies behandelt parallele Ausführung und führt Tupel ein, um zusammengehörige Dateien zusammenzuhalten.
4. **Ermögliche dem Workflow, ein Samplesheet mit mehreren Eingabedateien zu akzeptieren.**
   Dies demonstriert ein gängiges Muster zur Bereitstellung von Eingaben in großen Mengen.

Jeder Schritt konzentriert sich auf einen bestimmten Aspekt der Workflow-Entwicklung.

---

## 1. Schreibe einen einstufigen Workflow, der {TOOL_A_ACTION} ausführt

Dieser erste Schritt konzentriert sich auf die Grundlagen: Laden von {PRIMARY_INPUT_TYPE} und {TOOL_A_OUTPUT_DESCRIPTION}.

Erinnere dich an den Befehl `{TOOL_A_COMMAND_NAME}` aus [Teil 1](01_method.md):

```bash
{TOOL_A_COMMAND_SUMMARY}
```

Der Befehl nimmt {INPUT_DESCRIPTION} und erzeugt {OUTPUT_DESCRIPTION}.
Die Container-URI war `{TOOL_A_CONTAINER_URI}`.

Wir werden diese Informationen in drei Schritten in Nextflow einbinden:

1. Eingabe einrichten
2. Prozess schreiben und im Workflow aufrufen
3. Ausgabeverarbeitung konfigurieren

### 1.1. Eingabe einrichten

Wir müssen einen Eingabeparameter deklarieren, ein Testprofil erstellen, um einen praktischen Standardwert bereitzustellen, und einen Eingabekanal erstellen.

#### 1.1.1. Eingabeparameter-Deklaration hinzufügen

Deklariere in der Haupt-Workflow-Datei `{DOMAIN_DIR}.nf` unter dem Abschnitt `Pipeline parameters` einen CLI-Parameter namens `{PRIMARY_PARAM_NAME}`.

=== "Danach"

    ```groovy title="{DOMAIN_DIR}.nf" linenums="5" hl_lines="4-7"
    /*
     * Pipeline parameters
     */
    params {
        // Primary input
        {PRIMARY_PARAM_NAME}: Path
    }
    ```

=== "Vorher"

    ```groovy title="{DOMAIN_DIR}.nf" linenums="5"
    /*
     * Pipeline parameters
     */

    // Primary input
    ```

Das richtet den CLI-Parameter ein, aber wir wollen nicht jedes Mal den Dateipfad eingeben, wenn wir den Workflow während der Entwicklung ausführen.
Es gibt mehrere Möglichkeiten, einen Standardwert bereitzustellen; hier verwenden wir ein Testprofil.

#### 1.1.2. Testprofil mit Standardwert in `nextflow.config` erstellen

Ein Testprofil bietet praktische Standardwerte zum Ausprobieren eines Workflows, ohne Eingaben auf der Kommandozeile angeben zu müssen.
Dies ist eine gängige Konvention im Nextflow-Ökosystem (siehe [Hello Config](../../hello_nextflow/06_hello_config.md) für weitere Details).

Füge einen `profiles`-Block zu `nextflow.config` hinzu mit einem `test`-Profil, das den Parameter `{PRIMARY_PARAM_NAME}` auf eine der Test-{PRIMARY_INPUT_TYPE}s setzt.

=== "Danach"

    ```groovy title="nextflow.config" linenums="1" hl_lines="3-7"
    docker.enabled = true

    profiles {
        test {
            params.{PRIMARY_PARAM_NAME} = "${projectDir}/data/{TEST_INPUT_PATH}"
        }
    }
    ```

=== "Vorher"

    ```groovy title="nextflow.config" linenums="1"
    docker.enabled = true
    ```

Hier verwenden wir `${projectDir}`, eine eingebaute Nextflow-Variable, die auf das Verzeichnis zeigt, in dem sich das Workflow-Skript befindet.
Dies erleichtert das Referenzieren von Datendateien und anderen Ressourcen, ohne absolute Pfade fest zu codieren.

#### 1.1.3. Eingabekanal einrichten

{INPUT_CHANNEL_INSTRUCTIONS}

### 1.2. Das {TOOL_A_NAME}-Modul schreiben

{MODULE_INSTRUCTIONS_COVERING: container, input, output, script}

### 1.3. Modul im Workflow importieren und aufrufen

{IMPORT_AND_CALL_INSTRUCTIONS}

### 1.4. Workflow ausführen

An diesem Punkt haben wir einen einstufigen Workflow, der vollständig funktionsfähig sein sollte.

Wir können ihn mit `-profile test` ausführen, um den im Testprofil eingerichteten Standardwert zu verwenden und zu vermeiden, den Pfad auf der Kommandozeile eingeben zu müssen.

```bash
nextflow run {DOMAIN_DIR}.nf -profile test
```

??? success "Befehlsausgabe"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `{DOMAIN_DIR}.nf` [{RUN_NAME}] DSL2 - revision: {HASH}

    executor >  local (1)
    [{HASH}] {PROCESS_A_NAME} (1) | 1 of 1 ✔
    ```

{VERIFY_OUTPUT_INSTRUCTIONS}

### Fazit

Du weißt jetzt, wie du ein Modul mit einem Prozess erstellst, es in einen Workflow importierst, es mit einem Eingabekanal aufrufst und die Ergebnisse veröffentlichst.

### Wie geht es weiter?

Füge einen zweiten Prozess hinzu, um zusätzliche Analyseschritte zu verketten.

---

## 2. Füge einen zweiten Prozess hinzu, um {TOOL_B_ACTION} auszuführen

{DESCRIBE_WHAT_THIS_STEP_ADDS}

Erinnere dich an den Befehl `{TOOL_B_COMMAND_NAME}` aus [Teil 1](01_method.md):

```bash
{TOOL_B_COMMAND_SUMMARY}
```

### 2.1. Das {TOOL_B_NAME}-Modul schreiben

{MODULE_INSTRUCTIONS}

### 2.2. Modul im Workflow importieren und aufrufen

{IMPORT_AND_CALL_INSTRUCTIONS_CHAINING_OUTPUTS}

### 2.3. Workflow ausführen

```bash
nextflow run {DOMAIN_DIR}.nf -profile test
```

??? success "Befehlsausgabe"

    ```console
    {EXPECTED_OUTPUT}
    ```

{VERIFY_OUTPUT_INSTRUCTIONS}

### Fazit

Du weißt jetzt, wie du Prozessausgaben mit Eingaben verkettest und zusätzliche Dateien im Workflow verarbeitest.

### Wie geht es weiter?

Skaliere den Workflow, um mehrere Proben parallel zu verarbeiten.

---

## 3. Passe den Workflow an, um ihn auf mehreren Proben auszuführen

Bisher verarbeitet der Workflow eine einzelne Probe.
Um mehrere Proben zu verarbeiten, müssen wir ändern, wie Eingaben bereitgestellt werden, und Nextflows Dataflow-Paradigma nutzen, um die Ausführung zu parallelisieren.

### 3.1. {BATCH_INPUT_INSTRUCTIONS}

{INSTRUCTIONS_FOR_SWITCHING_TO_SAMPLESHEET_OR_BATCH_INPUT}

### 3.2. Workflow auf mehreren Proben ausführen

```bash
nextflow run {DOMAIN_DIR}.nf -profile test
```

??? success "Befehlsausgabe"

    ```console
    {EXPECTED_OUTPUT_MULTIPLE_SAMPLES}
    ```

### Fazit

Du weißt jetzt, wie du Nextflows Dataflow-Paradigma nutzt, um die Verarbeitung einzelner Proben über mehrere Eingabeproben hinweg zu parallelisieren.

### Wie geht es weiter?

In [Teil 3](03_multi_sample.md) wirst du eine probenübergreifende Aggregation hinzufügen, um Ergebnisse über alle Proben hinweg zu kombinieren.
