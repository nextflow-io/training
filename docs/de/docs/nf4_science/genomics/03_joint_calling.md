# Teil 3: Joint Calling für eine Kohorte

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } KI-gestützte Übersetzung - [mehr erfahren & Verbesserungen vorschlagen](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

In Teil 2 hast du eine Pipeline für die Variantenanalyse pro Probe erstellt, die die Daten jeder Probe unabhängig verarbeitet hat.
Jetzt erweitern wir sie, um Joint Variant Calling zu implementieren, wie in [Teil 1](01_method.md) behandelt.

## Aufgabe

In diesem Teil des Kurses erweitern wir den Workflow, um Folgendes zu tun:

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/genomics/img/hello-gatk-2.svg"
</figure>

1. Eine Indexdatei für jede BAM-Eingabedatei mit Samtools erstellen
2. GATK HaplotypeCaller auf jeder BAM-Eingabedatei ausführen, um ein GVCF mit genomischen Variantenaufrufen pro Probe zu erzeugen
3. Alle GVCFs sammeln und in einem GenomicsDB-Datenspeicher kombinieren
4. Joint Genotyping auf dem kombinierten GVCF-Datenspeicher ausführen, um ein VCF auf Kohortenebene zu erzeugen

Dieser Teil baut direkt auf dem Workflow aus Teil 2 auf.

??? info "Wie du von diesem Abschnitt aus beginnst"

    Dieser Abschnitt des Kurses setzt voraus, dass du [Teil 2: Variantenanalyse pro Probe](./02_per_sample_variant_calling.md) abgeschlossen hast und eine funktionierende `genomics.nf`-Pipeline besitzt.

    Falls du Teil 2 nicht abgeschlossen hast oder für diesen Teil neu beginnen möchtest, kannst du die Lösung von Teil 2 als Ausgangspunkt verwenden.
    Führe diese Befehle im Verzeichnis `nf4-science/genomics/` aus:

    ```bash
    cp solutions/part2/genomics-2.nf genomics.nf
    cp solutions/part2/nextflow.config .
    cp solutions/part2/modules/* modules/
    ```

    Damit erhältst du einen vollständigen Workflow für die Variantenanalyse pro Probe.
    Du kannst testen, ob er erfolgreich läuft, indem du folgenden Befehl ausführst:

    ```bash
    nextflow run genomics.nf -profile test
    ```

## Lektionsplan

Wir haben dies in zwei Schritte unterteilt:

1. **Den Schritt für die Variantenanalyse pro Probe so ändern, dass ein GVCF erzeugt wird.**
   Dies umfasst die Aktualisierung von Prozessbefehlen und Ausgaben.
2. **Einen Joint-Genotyping-Schritt hinzufügen, der die GVCFs pro Probe kombiniert und genotypisiert.**
   Dies führt den `collect()`-Operator, Groovy-Closures für die Befehlszeilenkonstruktion und Prozesse mit mehreren Befehlen ein.

!!! note "Hinweis"

     Stelle sicher, dass du im richtigen Arbeitsverzeichnis bist:
     `cd /workspaces/training/nf4-science/genomics`

---

## 1. Den Schritt für die Variantenanalyse pro Probe ändern, um ein GVCF zu erzeugen

Die Pipeline aus Teil 2 erzeugt VCF-Dateien, aber Joint Calling erfordert GVCF-Dateien.
Wir müssen den GVCF-Variantenaufrufmodus aktivieren und die Ausgabedateierweiterung aktualisieren.

Erinnere dich an den GVCF-Variantenaufrufbefehl aus [Teil 1](01_method.md):

```bash
gatk HaplotypeCaller \
        -R /data/ref/ref.fasta \
        -I /data/bam/reads_mother.bam \
        -O reads_mother.g.vcf \
        -L /data/ref/intervals.bed \
        -ERC GVCF
```

Im Vergleich zum Basis-HaplotypeCaller-Befehl, den wir in Teil 2 eingebunden haben, sind die Unterschiede der Parameter `-ERC GVCF` und die Ausgabeerweiterung `.g.vcf`.

### 1.1. HaplotypeCaller anweisen, ein GVCF auszugeben und die Ausgabeerweiterung aktualisieren

Öffne die Moduldatei `modules/gatk_haplotypecaller.nf`, um zwei Änderungen vorzunehmen:

- Füge den Parameter `-ERC GVCF` zum GATK HaplotypeCaller-Befehl hinzu;
- Aktualisiere den Ausgabedateipfad, um die entsprechende Erweiterung `.g.vcf` zu verwenden, gemäß GATK-Konvention.

Stelle sicher, dass du am Ende der vorherigen Zeile einen Backslash (`\`) hinzufügst, wenn du `-ERC GVCF` ergänzt.

=== "Danach"

    ```groovy title="modules/gatk_haplotypecaller.nf" linenums="22" hl_lines="5 7"
        """
        gatk HaplotypeCaller \
            -R ${ref_fasta} \
            -I ${input_bam} \
            -O ${input_bam}.g.vcf \
            -L ${interval_list} \
            -ERC GVCF
        """
    ```

=== "Vorher"

    ```groovy title="modules/gatk_haplotypecaller.nf" linenums="22" hl_lines="5"
        """
        gatk HaplotypeCaller \
            -R ${ref_fasta} \
            -I ${input_bam} \
            -O ${input_bam}.vcf \
            -L ${interval_list}
        """
    ```

Wir müssen auch den Output-Block aktualisieren, damit er zur neuen Dateierweiterung passt.
Da wir die Befehlsausgabe von `.vcf` auf `.g.vcf` geändert haben, muss der `output:`-Block des Prozesses dieselbe Änderung widerspiegeln.

### 1.2. Die Ausgabedateierweiterung im Prozess-Output-Block aktualisieren

=== "Danach"

    ```groovy title="modules/gatk_haplotypecaller.nf" linenums="17" hl_lines="2 3"
        output:
        path "${input_bam}.g.vcf"     , emit: vcf
        path "${input_bam}.g.vcf.idx" , emit: idx
    ```

=== "Vorher"

    ```groovy title="modules/gatk_haplotypecaller.nf" linenums="17" hl_lines="2 3"
        output:
        path "${input_bam}.vcf"     , emit: vcf
        path "${input_bam}.vcf.idx" , emit: idx
    ```

Wir müssen auch die Publish- und Output-Konfiguration des Workflows aktualisieren, um die neuen GVCF-Ausgaben widerzuspiegeln.

### 1.3. Die Publish-Ziele für die neuen GVCF-Ausgaben aktualisieren

Da wir jetzt GVCFs statt VCFs erzeugen, sollten wir den `publish:`-Abschnitt des Workflows aktualisieren, um aussagekräftigere Namen zu verwenden.
Wir organisieren die GVCF-Dateien auch in einem eigenen Unterverzeichnis für mehr Klarheit.

=== "Danach"

    ```groovy title="genomics.nf" linenums="47" hl_lines="3 4"
        publish:
        indexed_bam = SAMTOOLS_INDEX.out
        gvcf = GATK_HAPLOTYPECALLER.out.vcf
        gvcf_idx = GATK_HAPLOTYPECALLER.out.idx
    ```

=== "Vorher"

    ```groovy title="genomics.nf" linenums="47"
        publish:
        indexed_bam = SAMTOOLS_INDEX.out
        vcf = GATK_HAPLOTYPECALLER.out.vcf
        vcf_idx = GATK_HAPLOTYPECALLER.out.idx
    ```

Aktualisiere nun den Output-Block entsprechend.

### 1.4. Den Output-Block für die neue Verzeichnisstruktur aktualisieren

Wir müssen auch den `output`-Block aktualisieren, um die GVCF-Dateien in einem `gvcf`-Unterverzeichnis abzulegen.

=== "Danach"

    ```groovy title="genomics.nf" linenums="53" hl_lines="3 5 6 8 9"
    output {
        indexed_bam {
            path 'indexed_bam'
        }
        gvcf {
            path 'gvcf'
        }
        gvcf_idx {
            path 'gvcf'
        }
    }
    ```

=== "Vorher"

    ```groovy title="genomics.nf" linenums="53"
    output {
        indexed_bam {
            path 'bam'
        }
        vcf {
            path 'vcf'
        }
        vcf_idx {
            path 'vcf'
        }
    }
    ```

Mit dem aktualisierten Modul, den Publish-Zielen und dem Output-Block können wir die Änderungen testen.

### 1.5. Die Pipeline ausführen

Führe den Workflow aus, um zu überprüfen, ob die Änderungen funktionieren.

```bash
nextflow run genomics.nf
```

??? success "Befehlsausgabe"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics.nf` [nostalgic_franklin] DSL2 - revision: f2c0a93c6a

    executor >  local (6)
    [cc/fbc705] SAMTOOLS_INDEX (3)       | 3 of 3 ✔
    [27/0d7eb9] GATK_HAPLOTYPECALLER (2) | 3 of 3 ✔
    ```

Die Nextflow-Ausgabe sieht genauso aus wie zuvor, aber die `.g.vcf`-Dateien und ihre Indexdateien sind jetzt in Unterverzeichnissen organisiert.

??? abstract "Verzeichnisinhalt (Symlinks gekürzt)"

    ```console
    results/
    ├── gvcf/
    │   ├── reads_father.bam.g.vcf -> */27/0d7eb9*/reads_father.bam.g.vcf
    │   ├── reads_father.bam.g.vcf.idx -> */27/0d7eb9*/reads_father.bam.g.vcf.idx
    │   ├── reads_mother.bam.g.vcf -> */e4/4ed55e*/reads_mother.bam.g.vcf
    │   ├── reads_mother.bam.g.vcf.idx -> */e4/4ed55e*/reads_mother.bam.g.vcf.idx
    │   ├── reads_son.bam.g.vcf -> */08/e95962*/reads_son.bam.g.vcf
    │   └── reads_son.bam.g.vcf.idx -> */08/e95962*/reads_son.bam.g.vcf.idx
    └── indexed_bam/
        ├── reads_father.bam -> */9a/c7a873*/reads_father.bam
        ├── reads_father.bam.bai -> */9a/c7a873*/reads_father.bam.bai
        ├── reads_mother.bam -> */f1/8d8486*/reads_mother.bam
        ├── reads_mother.bam.bai -> */f1/8d8486*/reads_mother.bam.bai
        ├── reads_son.bam -> */cc/fbc705*/reads_son.bam
        └── reads_son.bam.bai -> */cc/fbc705*/reads_son.bam.bai
    ```

Wenn du eine der GVCF-Dateien öffnest und durchscrollst, kannst du überprüfen, dass GATK HaplotypeCaller wie gewünscht GVCF-Dateien erzeugt hat.

### Fazit

Wenn du den Ausgabedateinamen eines Tool-Befehls änderst, müssen der `output:`-Block des Prozesses und die Publish/Output-Konfiguration entsprechend aktualisiert werden.

### Wie geht es weiter?

Lerne, wie du den Inhalt eines Kanals sammelst und als einzelne Eingabe an den nächsten Prozess weitergibst.

---

## 2. Einen Joint-Genotyping-Schritt hinzufügen

Wir müssen jetzt die GVCFs pro Probe sammeln, sie in einem GenomicsDB-Datenspeicher kombinieren und Joint Genotyping ausführen, um ein VCF auf Kohortenebene zu erzeugen.
Wie in [Teil 1](01_method.md) behandelt, ist dies eine Operation mit zwei Tools: GenomicsDBImport kombiniert die GVCFs, dann erzeugt GenotypeGVCFs die finalen Variantenaufrufe.
Wir binden beide Tools in einem einzigen Prozess namens `GATK_JOINTGENOTYPING` ein.

Erinnere dich an die beiden Befehle aus [Teil 1](01_method.md):

```bash
gatk GenomicsDBImport \
    -V reads_mother.g.vcf \
    -V reads_father.g.vcf \
    -V reads_son.g.vcf \
    -L /data/ref/intervals.bed \
    --genomicsdb-workspace-path family_trio_gdb
```

```bash
gatk GenotypeGVCFs \
    -R /data/ref/ref.fasta \
    -V gendb://family_trio_gdb \
    -O family_trio.vcf
```

Der erste Befehl nimmt die GVCFs pro Probe und eine Intervalldatei und erzeugt einen GenomicsDB-Datenspeicher.
Der zweite nimmt diesen Datenspeicher, ein Referenzgenom und erzeugt das finale VCF auf Kohortenebene.
Die Container-URI ist dieselbe wie für HaplotypeCaller: `community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867`.

### 2.1. Die Eingaben einrichten

Der Joint-Genotyping-Prozess benötigt zwei Arten von Eingaben, die wir noch nicht haben: einen beliebigen Kohortennamen und die gesammelten GVCF-Ausgaben aller Proben zusammen gebündelt.

#### 2.1.1. Einen `cohort_name`-Parameter hinzufügen

Wir müssen einen beliebigen Namen für die Kohorte angeben.
Später in der Trainingsreihe lernst du, wie du Probenmetadaten für solche Dinge verwendest, aber vorerst deklarieren wir einfach einen CLI-Parameter mit `params` und geben ihm einen Standardwert zur Vereinfachung.

=== "Danach"

    ```groovy title="genomics.nf" linenums="14" hl_lines="3-4"
        intervals: Path = "${projectDir}/data/ref/intervals.bed"

        // Basisname für die finale Ausgabedatei
        cohort_name: String = "family_trio"
    }
    ```

=== "Vorher"

    ```groovy title="genomics.nf" linenums="14"
        intervals: Path = "${projectDir}/data/ref/intervals.bed"
    }
    ```

#### 2.1.2. Die HaplotypeCaller-Ausgaben über Proben hinweg sammeln

Wenn wir den Ausgabekanal von `GATK_HAPLOTYPECALLER` direkt in den neuen Prozess einspeisen würden, würde Nextflow den Prozess für jedes Proben-GVCF separat aufrufen.
Wir möchten alle drei GVCFs (und ihre Indexdateien) bündeln, sodass Nextflow sie alle zusammen an einen einzigen Prozessaufruf übergibt.

Das können wir mit dem `collect()`-Kanal-Operator tun.
Füge die folgenden Zeilen zum `workflow`-Body hinzu, direkt nach dem Aufruf von GATK_HAPLOTYPECALLER:

=== "Danach"

    ```groovy title="genomics.nf" hl_lines="4-6"
            intervals_file
        )

        // Variantenaufruf-Ausgaben über Proben hinweg sammeln
        all_gvcfs_ch = GATK_HAPLOTYPECALLER.out.vcf.collect()
        all_idxs_ch = GATK_HAPLOTYPECALLER.out.idx.collect()
    ```

=== "Vorher"

    ```groovy title="genomics.nf"
            intervals_file
        )
    ```

Das aufgeschlüsselt:

1. Wir nehmen den Ausgabekanal von `GATK_HAPLOTYPECALLER` mit der `.out`-Eigenschaft.
2. Da wir die Ausgaben in Abschnitt 1 mit `emit:` benannt haben, können wir die GVCFs mit `.vcf` und die Indexdateien mit `.idx` auswählen. Ohne benannte Ausgaben müssten wir `.out[0]` und `.out[1]` verwenden.
3. Der `collect()`-Operator bündelt alle Dateien in einem einzigen Element, sodass `all_gvcfs_ch` alle drei GVCFs zusammen enthält und `all_idxs_ch` alle drei Indexdateien zusammen enthält.

Wir können die GVCFs und ihre Indexdateien separat sammeln (anstatt sie in Tupeln zusammenzuhalten), weil Nextflow alle Eingabedateien zusammen für die Ausführung bereitstellt, sodass die Indexdateien neben den GVCFs vorhanden sein werden.

!!! tip "Tipp"

    Du kannst den `view()`-Operator verwenden, um den Inhalt von Kanälen vor und nach der Anwendung von Kanal-Operatoren zu inspizieren.

### 2.2. Den Joint-Genotyping-Prozess schreiben und im Workflow aufrufen

Nach demselben Muster, das wir in Teil 2 verwendet haben, schreiben wir die Prozessdefinition in eine Moduldatei, importieren sie in den Workflow und rufen sie mit den gerade vorbereiteten Eingaben auf.

#### 2.2.1. Einen String konstruieren, um jedem GVCF ein `-V`-Argument zu geben

Bevor wir beginnen, die Prozessdefinition auszufüllen, gibt es eine Sache zu klären.
Der GenomicsDBImport-Befehl erwartet für jede GVCF-Datei ein separates `-V`-Argument, so:

```bash
gatk GenomicsDBImport \
    -V reads_mother.bam.g.vcf \
    -V reads_father.bam.g.vcf \
    -V reads_son.bam.g.vcf \
    ...
```

Wenn wir `-V ${all_gvcfs_ch}` schreiben würden, würde Nextflow einfach die Dateinamen verketten und dieser Teil des Befehls würde so aussehen:

```groovy
-V reads_mother.bam.g.vcf reads_father.bam.g.vcf reads_son.bam.g.vcf
```

Aber wir brauchen den String so:

```groovy
-V reads_mother.bam.g.vcf -V reads_father.bam.g.vcf -V reads_son.bam.g.vcf
```

Wichtig ist, dass wir diesen String dynamisch aus den Dateien im gesammelten Kanal konstruieren müssen.
Nextflow (über Groovy) bietet eine prägnante Möglichkeit, dies zu tun:

```groovy
def gvcfs_line = all_gvcfs.collect { gvcf -> "-V ${gvcf}" }.join(' ')
```

Das aufgeschlüsselt:

1. `all_gvcfs.collect { gvcf -> "-V ${gvcf}" }` iteriert über jeden Dateipfad und stellt `-V ` voran, was `["-V A.g.vcf", "-V B.g.vcf", "-V C.g.vcf"]` ergibt.
2. `.join(' ')` verkettet sie mit Leerzeichen: `"-V A.g.vcf -V B.g.vcf -V C.g.vcf"`.
3. Das Ergebnis wird einer lokalen Variable `gvcfs_line` (definiert mit `def`) zugewiesen, die wir in die Befehlsvorlage interpolieren können.

Diese Zeile kommt in den `script:`-Block des Prozesses, vor die Befehlsvorlage.
Du kannst beliebigen Groovy-Code zwischen `script:` und dem öffnenden `"""` der Befehlsvorlage platzieren.

Dann kannst du auf diesen ganzen String als `gvcfs_line` im `script:`-Block des Prozesses verweisen.

#### 2.2.2. Das Modul für den Joint-Genotyping-Prozess ausfüllen

Jetzt können wir uns daran machen, den vollständigen Prozess zu schreiben.

Öffne `modules/gatk_jointgenotyping.nf` und untersuche die Gliederung der Prozessdefinition.

Fülle die Prozessdefinition mit den oben bereitgestellten Informationen aus und überprüfe dann deine Arbeit anhand der Lösung im Tab "Danach" unten.

=== "Vorher"

    ```groovy title="modules/gatk_jointgenotyping.nf" linenums="1"
    #!/usr/bin/env nextflow

    /*
     * GVCFs in GenomicsDB-Datenspeicher kombinieren und Joint Genotyping ausführen, um Aufrufe auf Kohortenebene zu erzeugen
     */
    process GATK_JOINTGENOTYPING {

        container

        input:

        output:

        script:
        """

        """
    }
    ```

=== "Danach"

    ```groovy title="modules/gatk_jointgenotyping.nf" linenums="1" hl_lines="8 11-17 20-21 24-25 29-33"
    #!/usr/bin/env nextflow

    /*
     * GVCFs in GenomicsDB-Datenspeicher kombinieren und Joint Genotyping ausführen, um Aufrufe auf Kohortenebene zu erzeugen
     */
    process GATK_JOINTGENOTYPING {

        container "community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867"

        input:
        path all_gvcfs
        path all_idxs
        path interval_list
        val cohort_name
        path ref_fasta
        path ref_index
        path ref_dict

        output:
        path "${cohort_name}.joint.vcf"     , emit: vcf
        path "${cohort_name}.joint.vcf.idx" , emit: idx

        script:
        def gvcfs_line = all_gvcfs.collect { gvcf -> "-V ${gvcf}" }.join(' ')
        """
        gatk GenomicsDBImport \
            ${gvcfs_line} \
            -L ${interval_list} \
            --genomicsdb-workspace-path ${cohort_name}_gdb

        gatk GenotypeGVCFs \
            -R ${ref_fasta} \
            -V gendb://${cohort_name}_gdb \
            -L ${interval_list} \
            -O ${cohort_name}.joint.vcf
        """
    }
    ```

Hier gibt es mehrere Dinge, die erwähnenswert sind.

Wie zuvor sind mehrere Eingaben aufgelistet, obwohl die Befehle nicht direkt auf sie verweisen: `all_idxs`, `ref_index` und `ref_dict`.
Ihre Auflistung stellt sicher, dass Nextflow diese Dateien im Arbeitsverzeichnis neben den Dateien bereitstellt, die in den Befehlen erscheinen, was GATK basierend auf Namenskonventionen erwartet.

Die Variable `gvcfs_line` verwendet die oben beschriebene Groovy-Closure, um die `-V`-Argumente für GenomicsDBImport zu konstruieren.

Dieser Prozess führt zwei Befehle seriell aus, genau wie du es im Terminal tun würdest.
GenomicsDBImport kombiniert die GVCFs pro Probe in einem Datenspeicher, dann liest GenotypeGVCFs diesen Datenspeicher und erzeugt das finale VCF auf Kohortenebene.
Der GenomicsDB-Datenspeicher (`${cohort_name}_gdb`) ist ein Zwischenartefakt, das nur innerhalb des Prozesses verwendet wird; er erscheint nicht im Output-Block.

Sobald du dies abgeschlossen hast, ist der Prozess einsatzbereit.
Um ihn im Workflow zu verwenden, musst du das Modul importieren und einen Prozessaufruf hinzufügen.

#### 2.2.3. Das Modul importieren

Füge die Import-Anweisung zu `genomics.nf` hinzu, unterhalb der bestehenden Import-Anweisungen:

=== "Danach"

    ```groovy title="genomics.nf" linenums="21" hl_lines="3"
    include { SAMTOOLS_INDEX } from './modules/samtools_index.nf'
    include { GATK_HAPLOTYPECALLER } from './modules/gatk_haplotypecaller.nf'
    include { GATK_JOINTGENOTYPING } from './modules/gatk_jointgenotyping.nf'
    ```

=== "Vorher"

    ```groovy title="genomics.nf" linenums="21"
    include { SAMTOOLS_INDEX } from './modules/samtools_index.nf'
    include { GATK_HAPLOTYPECALLER } from './modules/gatk_haplotypecaller.nf'
    ```

Der Prozess ist jetzt im Workflow-Scope verfügbar.

#### 2.2.4. Den Prozessaufruf hinzufügen

Füge den Aufruf von `GATK_JOINTGENOTYPING` im Workflow-Body hinzu, nach den `collect()`-Zeilen:

=== "Danach"

    ```groovy title="genomics.nf" hl_lines="3-12"
        all_idxs_ch = GATK_HAPLOTYPECALLER.out.idx.collect()

        // GVCFs in einem GenomicsDB-Datenspeicher kombinieren und Joint Genotyping anwenden
        GATK_JOINTGENOTYPING(
            all_gvcfs_ch,
            all_idxs_ch,
            intervals_file,
            params.cohort_name,
            ref_file,
            ref_index_file,
            ref_dict_file
        )
    ```

=== "Vorher"

    ```groovy title="genomics.nf"
        all_idxs_ch = GATK_HAPLOTYPECALLER.out.idx.collect()
    ```

Der Prozess ist jetzt vollständig verdrahtet.
Als Nächstes konfigurieren wir, wie die Ausgaben veröffentlicht werden.

### 2.3. Die Ausgabebehandlung konfigurieren

Wir müssen die Joint-VCF-Ausgaben veröffentlichen.
Füge Publish-Ziele und Output-Block-Einträge für die Joint-Genotyping-Ergebnisse hinzu.

#### 2.3.1. Publish-Ziele für das Joint-VCF hinzufügen

Füge das Joint-VCF und seinen Index zum `publish:`-Abschnitt des Workflows hinzu:

=== "Danach"

    ```groovy title="genomics.nf" hl_lines="5-6"
        publish:
        indexed_bam = SAMTOOLS_INDEX.out
        gvcf = GATK_HAPLOTYPECALLER.out.vcf
        gvcf_idx = GATK_HAPLOTYPECALLER.out.idx
        joint_vcf = GATK_JOINTGENOTYPING.out.vcf
        joint_vcf_idx = GATK_JOINTGENOTYPING.out.idx
    ```

=== "Vorher"

    ```groovy title="genomics.nf"
        publish:
        indexed_bam = SAMTOOLS_INDEX.out
        gvcf = GATK_HAPLOTYPECALLER.out.vcf
        gvcf_idx = GATK_HAPLOTYPECALLER.out.idx
    ```

Aktualisiere nun den Output-Block entsprechend.

#### 2.3.2. Output-Block-Einträge für das Joint-VCF hinzufügen

Füge Einträge für die Joint-VCF-Dateien hinzu.
Wir legen sie im Stammverzeichnis des Ergebnisverzeichnisses ab, da dies die finale Ausgabe ist.

=== "Danach"

    ```groovy title="genomics.nf" hl_lines="11-16"
    output {
        indexed_bam {
            path 'indexed_bam'
        }
        gvcf {
            path 'gvcf'
        }
        gvcf_idx {
            path 'gvcf'
        }
        joint_vcf {
            path '.'
        }
        joint_vcf_idx {
            path '.'
        }
    }
    ```

=== "Vorher"

    ```groovy title="genomics.nf"
    output {
        indexed_bam {
            path 'indexed_bam'
        }
        gvcf {
            path 'gvcf'
        }
        gvcf_idx {
            path 'gvcf'
        }
    }
    ```

Mit dem Prozess, den Publish-Zielen und dem Output-Block an Ort und Stelle können wir den vollständigen Workflow testen.

### 2.4. Den Workflow ausführen

Führe den Workflow aus, um zu überprüfen, dass alles funktioniert.

```bash
nextflow run genomics.nf -resume
```

??? success "Befehlsausgabe"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics.nf` [crazy_marconi] DSL2 - revision: 5da9afc841

    executor >  local (1)
    [9a/c7a873] SAMTOOLS_INDEX (2)       | 3 of 3, cached: 3 ✔
    [e4/4ed55e] GATK_HAPLOTYPECALLER (2) | 3 of 3, cached: 3 ✔
    [a6/7cc8ed] GATK_JOINTGENOTYPING     | 1 of 1 ✔
    ```

Die ersten beiden Schritte sind aus dem vorherigen Lauf gecacht, und der neue `GATK_JOINTGENOTYPING`-Schritt läuft einmal auf den gesammelten Eingaben aller drei Proben.
Die finale Ausgabedatei, `family_trio.joint.vcf` (und ihr Index), befinden sich im Ergebnisverzeichnis.

??? abstract "Verzeichnisinhalt (Symlinks gekürzt)"

    ```console
    results/
    ├── family_trio.joint.vcf -> */a6/7cc8ed*/family_trio.joint.vcf
    ├── family_trio.joint.vcf.idx -> */a6/7cc8ed*/family_trio.joint.vcf.idx
    ├── gvcf/
    │   ├── reads_father.bam.g.vcf -> */27/0d7eb9*/reads_father.bam.g.vcf
    │   ├── reads_father.bam.g.vcf.idx -> */27/0d7eb9*/reads_father.bam.g.vcf.idx
    │   ├── reads_mother.bam.g.vcf -> */e4/4ed55e*/reads_mother.bam.g.vcf
    │   ├── reads_mother.bam.g.vcf.idx -> */e4/4ed55e*/reads_mother.bam.g.vcf.idx
    │   ├── reads_son.bam.g.vcf -> */08/e95962*/reads_son.bam.g.vcf
    │   └── reads_son.bam.g.vcf.idx -> */08/e95962*/reads_son.bam.g.vcf.idx
    └── indexed_bam/
        ├── reads_father.bam -> */9a/c7a873*/reads_father.bam
        ├── reads_father.bam.bai -> */9a/c7a873*/reads_father.bam.bai
        ├── reads_mother.bam -> */f1/8d8486*/reads_mother.bam
        ├── reads_mother.bam.bai -> */f1/8d8486*/reads_mother.bam.bai
        ├── reads_son.bam -> */cc/fbc705*/reads_son.bam
        └── reads_son.bam.bai -> */cc/fbc705*/reads_son.bam.bai
    ```

Wenn du die Joint-VCF-Datei öffnest, kannst du überprüfen, dass der Workflow die erwarteten Variantenaufrufe erzeugt hat.

```console title="family_trio.joint.vcf" linenums="40"
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	reads_father	reads_mother	reads_son
20_10037292_10066351	3480	.	C	CT	1625.89	.	AC=5;AF=0.833;AN=6;BaseQRankSum=0.220;DP=85;ExcessHet=0.0000;FS=2.476;MLEAC=5;MLEAF=0.833;MQ=60.00;MQRankSum=0.00;QD=21.68;ReadPosRankSum=-1.147e+00;SOR=0.487	GT:AD:DP:GQ:PL	0/1:15,16:31:99:367,0,375	1/1:0,18:18:54:517,54,0	1/1:0,26:26:78:756,78,0
20_10037292_10066351	3520	.	AT	A	1678.89	.	AC=5;AF=0.833;AN=6;BaseQRankSum=1.03;DP=80;ExcessHet=0.0000;FS=2.290;MLEAC=5;MLEAF=0.833;MQ=60.00;MQRankSum=0.00;QD=22.39;ReadPosRankSum=0.701;SOR=0.730	GT:AD:DP:GQ:PL	0/1:18,13:31:99:296,0,424	1/1:0,18:18:54:623,54,0	1/1:0,26:26:78:774,78,0
20_10037292_10066351	3529	.	T	A	154.29	.	AC=1;AF=0.167;AN=6;BaseQRankSum=-5.440e-01;DP=104;ExcessHet=0.0000;FS=1.871;MLEAC=1;MLEAF=0.167;MQ=60.00;MQRankSum=0.00;QD=7.71;ReadPosRankSum=-1.158e+00;SOR=1.034	GT:AD:DP:GQ:PL	0/0:44,0:44:99:0,112,1347	0/1:12,8:20:99:163,0,328	0/0:39,0:39:99:0,105,1194
```

Du hast jetzt einen automatisierten, vollständig reproduzierbaren Joint-Variant-Calling-Workflow!

!!! note "Hinweis"

    Bedenke, dass die Datendateien, die wir dir gegeben haben, nur einen winzigen Teil von Chromosom 20 abdecken.
    Die tatsächliche Größe eines Variantenaufruf-Sets würde in Millionen von Varianten gezählt werden.
    Deshalb verwenden wir nur winzige Teilmengen von Daten für Trainingszwecke!

### Fazit

Du weißt jetzt, wie du Ausgaben aus einem Kanal sammelst und sie als einzelne Eingabe an einen anderen Prozess bündelst.
Du weißt auch, wie du eine Befehlszeile mit Groovy-Closures konstruierst und wie du mehrere Befehle in einem einzigen Prozess ausführst.

### Wie geht es weiter?

Klopf dir selbst auf die Schulter! Du hast den Nextflow-for-Genomics-Kurs abgeschlossen.

Gehe weiter zur finalen [Kurszusammenfassung](./next_steps.md), um zu überprüfen, was du gelernt hast, und herauszufinden, was als Nächstes kommt.
