# Kurszusammenfassung

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } KI-gestützte Übersetzung - [mehr erfahren & Verbesserungen vorschlagen](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Herzlichen Glückwunsch zum Abschluss des Nextflow for Genomics Trainingskurses! 🎉

## Deine Lernreise

Du hast damit begonnen, Variant-Calling-Tools manuell im Terminal auszuführen, um die Methodik zu verstehen.
Dann hast du eine Nextflow-Pipeline für eine einzelne Probe erstellt, um den Prozess zu automatisieren, sie für die parallele Verarbeitung mehrerer Proben skaliert und mit Channel-Operatoren ein gemeinsames Genotyping für mehrere Proben hinzugefügt.

### Was du erstellt hast

- Eine Variant-Calling-Pipeline, die BAM-Dateien als Eingabe nimmt und gemeinsam aufgerufene VCFs als Ausgabe erzeugt.
- Drei Prozesse (`SAMTOOLS_INDEX`, `GATK_HAPLOTYPECALLER` und `GATK_JOINTGENOTYPING`), die in separaten Moduldateien gespeichert sind.
- Die Pipeline parallelisiert automatisch die Verarbeitung der Eingabeproben durch Nextflows Dataflow-Paradigma.
- Die Ergebnisse werden in einem Verzeichnis namens `results/` veröffentlicht.

### Erworbene Fähigkeiten

Durch diesen praxisorientierten Kurs hast du gelernt, wie du:

- Einen linearen Workflow schreibst, um Variant Calling auf eine einzelne Probe anzuwenden
- Zusatzdateien wie Index-Dateien und Referenzgenom-Ressourcen angemessen handhabst
- Nextflows Dataflow-Paradigma nutzt, um Variant Calling pro Probe zu parallelisieren
- Gemeinsames Genotyping für mehrere Proben mit relevanten Channel-Operatoren implementierst

Du bist jetzt bereit, Nextflow auf Genomik-Analyse-Workflows in deiner eigenen Arbeit anzuwenden.

## Nächste Schritte zum Ausbau deiner Fähigkeiten

Hier sind unsere Top-Empfehlungen für die nächsten Schritte:

- Wende Nextflow auf andere wissenschaftliche Analyse-Anwendungsfälle an mit [Nextflow for Science](../index.md)
- Starte mit nf-core durch [Hello nf-core](../../hello_nf-core/index.md)
- Erkunde fortgeschrittenere Nextflow-Features mit den [Side Quests](../../side_quests/index.md)

Abschließend empfehlen wir dir, einen Blick auf [**Seqera Platform**](https://seqera.io/) zu werfen, eine cloudbasierte Plattform, die von den Entwickler\*innen von Nextflow erstellt wurde und es noch einfacher macht, deine Workflows zu starten und zu verwalten sowie deine Daten zu managen und Analysen interaktiv in jeder Umgebung auszuführen.

## Hilfe erhalten

Hilfsressourcen und Community-Support findest du auf der [Hilfe-Seite](../../help.md).

## Feedback-Umfrage

Bevor du weitermachst, nimm dir bitte eine Minute Zeit, um die Kursumfrage auszufüllen! Dein Feedback hilft uns, unsere Trainingsmaterialien für alle zu verbessern.

[Zur Umfrage :material-arrow-right:](survey.md){ .md-button .md-button--primary }
