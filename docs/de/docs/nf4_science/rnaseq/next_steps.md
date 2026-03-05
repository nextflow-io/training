# Kurszusammenfassung

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } KI-gestützte Übersetzung - [mehr erfahren & Verbesserungen vorschlagen](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Herzlichen Glückwunsch zum Abschluss des Nextflow for RNAseq Trainingskurses!

## Deine Reise

Du hast damit begonnen, RNAseq-Verarbeitungstools manuell im Terminal auszuführen, um die Methodik zu verstehen.
Dann hast du eine Nextflow-Pipeline für eine einzelne Probe erstellt, um den Prozess zu automatisieren, sie auf die parallele Verarbeitung mehrerer Proben skaliert und erweitert, um Paired-End-Daten zu verarbeiten und QC-Berichte über alle Proben hinweg zu aggregieren.

### Was du erstellt hast

- Eine RNAseq-Verarbeitungspipeline, die FASTQ-Dateien als Eingabe nimmt und getrimmte Reads, Alignments und aggregierte QC-Berichte als Ausgabe produziert.
- Prozesse für Trimming (Trim Galore), Alignment (HISAT2), Qualitätskontrolle (FastQC) und Berichtsaggregation (MultiQC), die in separaten Moduldateien gespeichert sind.
- Die Pipeline parallelisiert automatisch die Verarbeitung der Eingabeproben mithilfe von Nextflows Dataflow-Paradigma.
- Die finale Pipeline verarbeitet Paired-End-Sequenzierungsdaten.

### Erworbene Fähigkeiten

Durch diesen praxisorientierten Kurs hast du gelernt, wie du:

- Einen linearen Workflow schreibst, um grundlegende RNAseq-Verarbeitungs- und QC-Methoden anzuwenden
- Domänenspezifische Dateien wie FASTQ und Referenzgenom-Ressourcen angemessen handhabst
- Single-End- und Paired-End-Sequenzierungsdaten verarbeitest
- Nextflows Dataflow-Paradigma nutzt, um die RNAseq-Verarbeitung pro Probe zu parallelisieren
- QC-Berichte über mehrere Schritte und Proben hinweg mithilfe relevanter Kanal-Operatoren aggregierst

Du bist jetzt bereit, Nextflow auf RNAseq-Analyse-Workflows in deiner eigenen Arbeit anzuwenden.

## Nächste Schritte zum Ausbau deiner Kenntnisse

Hier sind unsere wichtigsten Empfehlungen, was du als Nächstes tun kannst:

- Wende Nextflow auf andere wissenschaftliche Analyse-Anwendungsfälle an mit [Nextflow for Science](../index.md)
- Starte mit nf-core durch [Hello nf-core](../../hello_nf-core/index.md)
- Erkunde fortgeschrittenere Nextflow-Features mit den [Side Quests](../../side_quests/index.md)

Abschließend empfehlen wir dir, einen Blick auf die [**Seqera Platform**](https://seqera.io/) zu werfen, eine cloudbasierte Plattform, die von den Entwickler\*innen von Nextflow erstellt wurde und es noch einfacher macht, deine Workflows zu starten und zu verwalten sowie deine Daten zu verwalten und Analysen interaktiv in jeder Umgebung auszuführen.

## Hilfe erhalten

Für Hilferessourcen und Community-Support siehe die [Hilfe-Seite](../../help.md).

## Feedback-Umfrage

Bevor du weitermachst, nimm dir bitte eine Minute Zeit, um die Kursumfrage auszufüllen! Dein Feedback hilft uns, unsere Trainingsmaterialien für alle zu verbessern.

[Zur Umfrage :material-arrow-right:](survey.md){ .md-button .md-button--primary }
