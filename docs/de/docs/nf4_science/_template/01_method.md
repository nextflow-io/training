# Teil 1: Methodenübersicht und manuelles Testen

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } KI-gestützte Übersetzung - [mehr erfahren & Verbesserungen vorschlagen](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

{BRIEF_METHOD_DESCRIPTION}

![Pipeline-Übersicht](img/{PIPELINE_DIAGRAM}.png)

{OPTIONAL_EXTENDED_METHOD_DESCRIPTION}

### Methoden

{DESCRIBE_THE_TWO_APPROACHES: single-sample and multi-sample/aggregation}

Bevor wir mit dem Schreiben von Workflow-Code beginnen, probieren wir die Befehle manuell mit einigen Testdaten aus.

### Datensatz

Wir stellen die folgenden Daten und zugehörigen Ressourcen bereit:

- **{PRIMARY_INPUT_DESCRIPTION}**
- **{SAMPLE_DESCRIPTION}**
- **{ADDITIONAL_RESOURCES_DESCRIPTION}**

### Software

Die wichtigsten Tools sind [{TOOL_A}]({TOOL_A_URL}) und [{TOOL_B}]({TOOL_B_URL}).

Diese Tools sind in der GitHub Codespaces-Umgebung nicht installiert, daher verwenden wir sie über Container (siehe [Hello Containers](../../hello_nextflow/05_hello_containers.md)).

!!! note "Hinweis"

    Stelle sicher, dass du dich im Verzeichnis `nf4-science/{DOMAIN_DIR}` befindest, sodass der letzte Teil des Pfads, der angezeigt wird, wenn du `pwd` eingibst, `{DOMAIN_DIR}` ist.

---

## 1. {SINGLE_SAMPLE_SECTION_TITLE}

{BRIEF_DESCRIPTION_OF_SINGLE_SAMPLE_APPROACH}

In diesem Abschnitt testen wir die Befehle, aus denen der Ansatz zur Verarbeitung einzelner Proben besteht.
Dies sind die Befehle, die wir in Teil 2 dieses Kurses in einen Nextflow-Workflow einbinden werden.

1. {STEP_1_SUMMARY}
2. {STEP_2_SUMMARY}

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/{DOMAIN_DIR}/img/{DIAGRAM_1}.svg"
</figure>

Wir beginnen damit, die Befehle an nur einer Probe zu testen.

### 1.1. {FIRST_TOOL_STEP_TITLE}

{BRIEF_DESCRIPTION_OF_STEP}

#### 1.1.1. Container herunterladen

Führe den Befehl `docker pull` aus, um das Container-Image herunterzuladen:

```bash
docker pull {TOOL_A_CONTAINER_URI}
```

??? success "Befehlsausgabe"

    ```console
    {EXPECTED_PULL_OUTPUT}
    ```

#### 1.1.2. Container interaktiv starten

Starte den Container und mounte das Verzeichnis `data`, damit die Tools auf die Eingabedateien zugreifen können:

```bash
docker run -it -v ./data:/data {TOOL_A_CONTAINER_URI}
```

Deine Eingabeaufforderung ändert sich, um anzuzeigen, dass du dich im Container befindest.

#### 1.1.3. Befehl ausführen

```bash
{TOOL_A_COMMAND}
```

??? success "Befehlsausgabe"

    ```console
    {EXPECTED_COMMAND_OUTPUT}
    ```

{DESCRIBE_EXPECTED_OUTPUT_FILES}

#### 1.1.4. Container verlassen

Um den Container zu verlassen, gib `exit` ein.

```bash
exit
```

Deine Eingabeaufforderung sollte wieder normal sein.

### 1.2. {SECOND_TOOL_STEP_TITLE}

{BRIEF_DESCRIPTION_OF_STEP}

#### 1.2.1. Container herunterladen

```bash
docker pull {TOOL_B_CONTAINER_URI}
```

??? success "Befehlsausgabe"

    ```console
    {EXPECTED_PULL_OUTPUT}
    ```

#### 1.2.2. Container interaktiv starten

```bash
docker run -it -v ./data:/data {TOOL_B_CONTAINER_URI}
```

#### 1.2.3. Befehl ausführen

```bash
{TOOL_B_COMMAND}
```

??? success "Befehlsausgabe"

    ```console
    {EXPECTED_COMMAND_OUTPUT}
    ```

{DESCRIBE_EXPECTED_OUTPUT_FILES}

#### 1.2.4. Container verlassen

Um den Container zu verlassen, gib `exit` ein.

```bash
exit
```

Deine Eingabeaufforderung sollte wieder normal sein.
Damit ist der Test zur Verarbeitung einzelner Proben abgeschlossen.

---

## 2. {MULTI_SAMPLE_SECTION_TITLE}

{BRIEF_DESCRIPTION_OF_MULTI_SAMPLE_APPROACH}

{EXPLAIN_WHY_MULTI_SAMPLE_IS_NEEDED}

In diesem Abschnitt testen wir die zusätzlichen Befehle, die für die Verarbeitung mehrerer Proben benötigt werden.
Dies sind die Befehle, die wir in Teil 3 dieses Kurses in einen Nextflow-Workflow einbinden werden.

1. {AGGREGATION_STEP_1_SUMMARY}
2. {AGGREGATION_STEP_2_SUMMARY}

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/{DOMAIN_DIR}/img/{DIAGRAM_2}.svg"
</figure>

### 2.1. {AGGREGATION_STEP_TITLE}

{STEP_INSTRUCTIONS_FOLLOWING_SAME_PATTERN_AS_ABOVE}

---

### Fazit

Du weißt jetzt, wie du die Befehle von {TOOL_A} und {TOOL_B} in ihren jeweiligen Containern testest, einschließlich {MULTI_SAMPLE_SUMMARY}.

### Wie geht es weiter?

Lerne, wie du dieselben Befehle in Workflows einbindest, die Container zur Ausführung der Arbeit verwenden.
