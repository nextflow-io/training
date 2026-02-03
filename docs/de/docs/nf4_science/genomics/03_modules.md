# Teil 3: Code in Module verschieben

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } KI-gestützte Übersetzung - [mehr erfahren & Verbesserungen vorschlagen](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Im ersten Teil dieses Kurses hast du eine Variant-Calling-Pipeline erstellt, die vollständig linear war und die Daten jeder Probe unabhängig von den anderen verarbeitet hat.

Im zweiten Teil haben wir dir gezeigt, wie du Channels und Channel-Operatoren verwendest, um Joint-Variant-Calling mit GATK zu implementieren und dabei auf der Pipeline aus Teil 1 aufzubauen.

In diesem Teil zeigen wir dir, wie du den Code in diesem Workflow in Module umwandelst. Um diesem Teil des Trainings zu folgen, solltest du Teil 1 und Teil 2 sowie [Hello Modules](../../../hello_nextflow/hello_modules.md) abgeschlossen haben, das die Grundlagen von Modulen behandelt.

---

## 0. Aufwärmen

Als wir begonnen haben, unseren Workflow zu entwickeln, haben wir alles in eine einzige Code-Datei gepackt.
Jetzt ist es an der Zeit, unseren Code zu **modularisieren**, _d.h._ die Prozess-Definitionen in Module zu extrahieren.

Wir beginnen mit demselben Workflow wie in Teil 2, den wir für dich in der Datei `genomics-3.nf` bereitgestellt haben.

!!! note "Hinweis"

     Stelle sicher, dass du im richtigen Arbeitsverzeichnis bist:
     `cd /workspaces/training/nf4-science/genomics`

Führe den Workflow aus, um den Ausgangspunkt zu überprüfen:

```bash
nextflow run genomics-3.nf -resume
```

```console title="Ausgabe"
 N E X T F L O W   ~  version 25.10.2

Launching `genomics-3.nf` [serene_borg] DSL2 - revision: 0cbebb67a1

executor >  local (7)
[6f/83ee72] SAMTOOLS_INDEX (3)       | 3 of 3 ✔
[53/b9d342] GATK_HAPLOTYPECALLER (1) | 3 of 3 ✔
[0c/fa6d15] GATK_JOINTGENOTYPING     | 1 of 1 ✔
```

Es gibt nun ein `work`-Verzeichnis und ein `results_genomics`-Verzeichnis in deinem Projektverzeichnis.

### Zusammenfassung

Du bist bereit, mit der Modularisierung deines Workflows zu beginnen.

### Wie geht es weiter?

Verschiebe die Prozesse des Genomics-Workflows in Module.

---

## 1. Prozesse in Module verschieben

Wie du in [Hello Modules](../../../hello_nextflow/hello_modules.md) gelernt hast, kannst du ein Modul einfach erstellen, indem du die Prozess-Definition in eine eigene Datei kopierst, in einem beliebigen Verzeichnis, und du kannst diese Datei beliebig benennen.

Aus Gründen, die später klar werden (insbesondere wenn wir zum Testen kommen), werden wir in diesem Training der Konvention folgen, die Datei `main.nf` zu nennen und sie in eine Verzeichnisstruktur zu legen, die nach dem Toolkit und dem Befehl benannt ist.

### 1.1. Erstelle ein Modul für den `SAMTOOLS_INDEX`-Prozess

Im Fall des `SAMTOOLS_INDEX`-Prozesses ist 'samtools' das Toolkit und 'index' der Befehl. Wir erstellen also eine Verzeichnisstruktur `modules/samtools/index` und legen die `SAMTOOLS_INDEX`-Prozess-Definition in die `main.nf`-Datei in diesem Verzeichnis.

```bash
mkdir -p modules/samtools/index
touch modules/samtools/index/main.nf
```

Öffne die `main.nf`-Datei und kopiere die `SAMTOOLS_INDEX`-Prozess-Definition hinein.

```groovy title="modules/samtools/index/main.nf" linenums="1"
/*
 * BAM-Index-Datei generieren
 */
process SAMTOOLS_INDEX {

    container 'community.wave.seqera.io/library/samtools:1.20--b5dfbd93de237464'

    input:
    path input_bam

    output:
    tuple path(input_bam), path("${input_bam}.bai")

    script:
    """
    samtools index '$input_bam'
    """
}
```

Entferne dann die `SAMTOOLS_INDEX`-Prozess-Definition aus `genomics-3.nf` und füge eine Import-Deklaration für das Modul vor der nächsten Prozess-Definition hinzu, so wie hier:

=== "Nachher"

    ```groovy title="genomics-3.nf" linenums="1" hl_lines="1 2"
    // Module einbinden
    include { SAMTOOLS_INDEX } from './modules/samtools/index/main.nf'

    /*
     * Varianten mit GATK HaplotypeCaller aufrufen
     */
    process GATK_HAPLOTYPECALLER {
    ```

=== "Vorher"

    ```groovy title="genomics-3.nf" linenums="1" hl_lines="1"
    /*
     * Varianten mit GATK HaplotypeCaller aufrufen
     */
    process GATK_HAPLOTYPECALLER {
    ```

Du kannst den Workflow jetzt erneut ausführen, und er sollte immer noch auf die gleiche Weise funktionieren. Wenn du das Flag `-resume` verwendest, sollten nicht einmal neue Aufgaben ausgeführt werden müssen:

```bash
nextflow run genomics-3.nf -resume
```

??? success "Befehlsausgabe"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `genomics-3.nf` [sleepy_snyder] DSL2 - revision: aa68d06c43

    [0f/71b55e] SAMTOOLS_INDEX (1)       | 3 of 3, cached: 3 ✔
    [f1/18971b] GATK_HAPLOTYPECALLER (3) | 3 of 3, cached: 3 ✔
    [0c/fa6d15] GATK_JOINTGENOTYPING     | 1 of 1, cached: 1 ✔
    ```

### 1.2. Erstelle Module für die Prozesse `GATK_HAPLOTYPECALLER` und `GATK_JOINTGENOTYPING`

Wiederhole die gleichen Schritte für die verbleibenden Prozesse.
Für jeden Prozess:

1. Erstelle die Verzeichnisstruktur (`modules/gatk/haplotypecaller/` und `modules/gatk/jointgenotyping/`)
2. Erstelle eine `main.nf`-Datei mit der Prozess-Definition
3. Entferne die Prozess-Definition aus `genomics-3.nf`
4. Füge eine Import-Deklaration für das Modul hinzu

Überprüfe nach Abschluss, ob deine Modulverzeichnisstruktur korrekt ist, indem du Folgendes ausführst:

```bash
tree modules/
```

??? abstract "Verzeichnisinhalte"

    ```console
    modules/
    ├── gatk
    │   ├── haplotypecaller
    │   │   └── main.nf
    │   └── jointgenotyping
    │       └── main.nf
    └── samtools
        └── index
            └── main.nf

    5 directories, 3 files
    ```

Du solltest auch etwas wie dies in der Haupt-Workflow-Datei haben, nach dem Parameterabschnitt:

```
include { SAMTOOLS_INDEX } from './modules/samtools/index/main.nf'
include { GATK_HAPLOTYPECALLER } from './modules/gatk/haplotypecaller/main.nf'
include { GATK_JOINTGENOTYPING } from './modules/gatk/jointgenotyping/main.nf'

workflow {
```

### Zusammenfassung

Du hast das Modularisieren eines Workflows geübt, mit dem Genomics-Workflow als Beispiel.

### Wie geht es weiter?

Teste den modularisierten Workflow.

---

## 2. Den modularisierten Workflow testen

Führe den modularisierten Workflow aus, um zu überprüfen, dass alles noch funktioniert.

```bash
nextflow run genomics-3.nf -resume
```

```console title="Ausgabe"
 N E X T F L O W   ~  version 25.10.2

Launching `genomics-3.nf` [astonishing_venter] DSL2 - revision: ca27264c13

[6f/83ee72] SAMTOOLS_INDEX (3)       | 3 of 3, cached: 3 ✔
[53/b9d342] GATK_HAPLOTYPECALLER (3) | 3 of 3, cached: 3 ✔
[0c/fa6d15] GATK_JOINTGENOTYPING     | 1 of 1, cached: 1 ✔
```

Alles funktioniert noch, einschließlich der Wiederaufnahmefähigkeit der Pipeline.
Die Ergebnisse werden weiterhin im Verzeichnis `results_genomics` veröffentlicht.

```console title="Verzeichnisinhalte"
results_genomics/
├── family_trio.joint.vcf
├── family_trio.joint.vcf.idx
├── gvcf
│   ├── reads_father.bam.g.vcf
│   ├── reads_father.bam.g.vcf.idx
│   ├── reads_mother.bam.g.vcf
│   ├── reads_mother.bam.g.vcf.idx
│   ├── reads_son.bam.g.vcf
│   └── reads_son.bam.g.vcf.idx
└── indexed_bam
    ├── reads_father.bam
    ├── reads_father.bam.bai
    ├── reads_mother.bam
    ├── reads_mother.bam.bai
    ├── reads_son.bam
    └── reads_son.bam.bai
```

### Zusammenfassung

Du hast einen Workflow modularisiert und überprüft, dass er immer noch genauso funktioniert wie zuvor.

### Wie geht es weiter?

Fasse zusammen, was du gelernt hast, und schau voraus auf das Testen.

---

## 3. Zusammenfassung

Du hast den Workflow modularisiert, und nichts hat sich an der Funktionsweise der Pipeline geändert.
Das ist beabsichtigt: Du hast den Code umstrukturiert, ohne seine Funktion zu beeinflussen.

Die Module enthalten nur die Prozess-Logik, was sie sauber und wiederverwendbar macht.
Das Haupt-Script kontrolliert, was wo veröffentlicht wird, während die Module auf ihre Berechnungsaufgabe fokussiert bleiben.

Du hast eine Grundlage für Dinge geschaffen, die deinen Code einfacher wartbar machen.
Du kannst jetzt zum Beispiel Tests zu deiner Pipeline hinzufügen, indem du das nf-test-Framework verwendest.
Das werden wir uns im nächsten Teil dieses Kurses ansehen.
