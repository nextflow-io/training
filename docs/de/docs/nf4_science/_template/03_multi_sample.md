# Teil 3: Aggregation mehrerer Proben

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } KI-gestützte Übersetzung - [mehr erfahren & Verbesserungen vorschlagen](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

In Teil 2 hast du eine Pipeline zur Verarbeitung einzelner Proben erstellt, die jede Probe unabhängig verarbeitet hat.
Jetzt erweitern wir sie, um die {AGGREGATION_METHOD} mehrerer Proben zu implementieren, wie in [Teil 1](01_method.md) behandelt.

## Aufgabe

In diesem Teil des Kurses erweitern wir den Workflow, um Folgendes zu tun:

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/{DOMAIN_DIR}/img/{DIAGRAM_2}.svg"
</figure>

1. {PER_SAMPLE_STEP_1}
2. {PER_SAMPLE_STEP_2}
3. {AGGREGATION_STEP_1}
4. {AGGREGATION_STEP_2}

Dieser Teil baut direkt auf dem Workflow aus Teil 2 auf.

??? info "Wie du von diesem Abschnitt aus beginnst"

    Dieser Abschnitt des Kurses setzt voraus, dass du [Teil 2: Verarbeitung einzelner Proben](./02_single_sample.md) abgeschlossen hast und eine funktionierende `{DOMAIN_DIR}.nf`-Pipeline besitzt.

    Falls du Teil 2 nicht abgeschlossen hast oder für diesen Teil neu beginnen möchtest, kannst du die Lösung von Teil 2 als Ausgangspunkt verwenden.
    Führe diese Befehle im Verzeichnis `nf4-science/{DOMAIN_DIR}/` aus:

    ```bash
    cp solutions/part2/{DOMAIN_DIR}-2.nf {DOMAIN_DIR}.nf
    cp solutions/part2/nextflow.config .
    cp solutions/part2/modules/* modules/
    ```

    Damit erhältst du einen vollständigen Workflow zur Verarbeitung einzelner Proben.
    Du kannst testen, ob er erfolgreich läuft, indem du folgenden Befehl ausführst:

    ```bash
    nextflow run {DOMAIN_DIR}.nf -profile test
    ```

## Lektionsplan

Wir haben dies in zwei Schritte unterteilt:

1. **{MODIFICATION_STEP_SUMMARY}.**
   Dies umfasst die Aktualisierung von Prozessbefehlen und Ausgaben.
2. **{AGGREGATION_STEP_SUMMARY}.**
   Dies führt den `collect()`-Operator {AND_OTHER_CONCEPTS} ein.

!!! note "Hinweis"

     Stelle sicher, dass du im richtigen Arbeitsverzeichnis bist:
     `cd /workspaces/training/nf4-science/{DOMAIN_DIR}`

---

## 1. {MODIFICATION_STEP_TITLE}

{DESCRIPTION_OF_MODIFICATION_TO_EXISTING_PROCESSES}

Erinnere dich an den modifizierten Befehl aus [Teil 1](01_method.md):

```bash
{MODIFIED_COMMAND}
```

{EXPLAIN_DIFFERENCES_FROM_PART_2}

### 1.1. {MODIFICATION_SUBSTEP}

{INSTRUCTIONS_WITH_BEFORE_AFTER_TABS}

=== "Danach"

    ```groovy title="modules/{MODULE_FILE}.nf" linenums="{N}" hl_lines="{LINES}"
    {UPDATED_CODE}
    ```

=== "Vorher"

    ```groovy title="modules/{MODULE_FILE}.nf" linenums="{N}" hl_lines="{LINES}"
    {ORIGINAL_CODE}
    ```

### 1.2. Führe den Workflow aus, um die Änderung zu überprüfen

```bash
nextflow run {DOMAIN_DIR}.nf -profile test
```

??? success "Befehlsausgabe"

    ```console
    {EXPECTED_OUTPUT}
    ```

{VERIFY_OUTPUT_INSTRUCTIONS}

### Fazit

Du weißt, wie du Prozessbefehle und Ausgaben änderst, um das Verhalten des Workflows anzupassen.

### Wie geht es weiter?

Füge den Schritt zur Aggregation mehrerer Proben hinzu.

---

## 2. {AGGREGATION_STEP_TITLE}

{DESCRIPTION_OF_AGGREGATION}

### 2.1. Schreibe das Aggregationsmodul

{MODULE_INSTRUCTIONS}

### 2.2. Sammle die Ausgaben der einzelnen Proben und übergebe sie an den Aggregationsprozess

{INSTRUCTIONS_USING_COLLECT_OPERATOR}

### 2.3. Führe den vollständigen Workflow aus

```bash
nextflow run {DOMAIN_DIR}.nf -profile test
```

??? success "Befehlsausgabe"

    ```console
    {EXPECTED_OUTPUT}
    ```

{VERIFY_OUTPUT_INSTRUCTIONS}

---

### Fazit

Du hast eine vollständige Pipeline, die Proben einzeln verarbeitet und die Ergebnisse über alle Proben hinweg aggregiert.
Du weißt, wie du Kanaloperatoren wie `collect()` verwendest, um Ausgaben einzelner Proben für die Analyse mehrerer Proben zu aggregieren.

### Wie geht es weiter?

Herzlichen Glückwunsch zum Abschluss dieses Kurses! Gehe zur [Kurszusammenfassung](next_steps.md), um zu wiederholen, was du gelernt hast, und weitere Schritte zu erkunden.
