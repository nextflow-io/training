# Kurszusammenfassung

Herzlichen Glückwunsch zum Abschluss des Nextflow For Genomics Trainingskurses! 🎉

## Deine Lernreise

Du hast damit begonnen, Variant-Calling-Tools manuell im Terminal auszuführen, um die Methodik zu verstehen.
Dann hast du eine Nextflow-Pipeline für eine einzelne Probe erstellt, um den Prozess zu automatisieren, sie auf die parallele Verarbeitung mehrerer Proben skaliert und mit Channel-Operatoren Multi-Sample-Joint-Genotyping hinzugefügt.

### Was du erstellt hast

- Eine Variant-Calling-Pipeline, die BAM-Dateien als Eingabe nimmt und gemeinsam gecallte VCF-Dateien als Ausgabe produziert.
- Drei Prozesse (`SAMTOOLS_INDEX`, `GATK_HAPLOTYPECALLER` und `GATK_JOINTGENOTYPING`), die in separaten Moduldateien gespeichert sind.
- Die Pipeline skaliert automatisch auf eine beliebige Anzahl von Eingabeproben unter Verwendung des Dataflow-Paradigmas von Nextflow.
- Die Ergebnisse werden in einem Verzeichnis namens `results/` veröffentlicht.

### Erworbene Fähigkeiten

Durch diesen praxisorientierten Kurs hast du gelernt, wie du:

- Einen linearen Workflow schreibst, um Variant Calling auf eine einzelne Probe anzuwenden
- Zusatzdateien wie Index-Dateien und Referenzgenom-Ressourcen angemessen handhabst
- Das Dataflow-Paradigma von Nextflow nutzt, um probenspezifisches Variant Calling zu parallelisieren
- Multi-Sample-Joint-Calling mit relevanten Channel-Operatoren implementierst

Du bist jetzt bereit, Nextflow auf Genomanalyse-Workflows in deiner eigenen Arbeit anzuwenden.

## Nächste Schritte, um deine Fähigkeiten auszubauen

Hier sind unsere wichtigsten Empfehlungen, was du als Nächstes tun kannst:

- Wende Nextflow auf andere wissenschaftliche Anwendungsfälle an mit [Nextflow for Science](../index.md)
- Starte mit nf-core mit [Hello nf-core](../../hello_nf-core/index.md)
- Erkunde fortgeschrittenere Nextflow-Features mit den [Side Quests](../../side_quests/index.md)

Abschließend empfehlen wir dir, einen Blick auf [**Seqera Platform**](https://seqera.io/) zu werfen – eine cloudbasierte Plattform, die vom Nextflow-Entwicklungsteam erstellt wurde und es dir noch einfacher macht, deine Workflows zu starten und zu verwalten, sowie deine Daten zu managen und Analysen interaktiv in jeder Umgebung auszuführen.

## Hilfe erhalten

Für Hilferessourcen und Community-Support siehe die [Hilfe-Seite](../../help.md).

## Feedback-Umfrage

Bevor du weitermachst, nimm dir bitte eine Minute Zeit, um die Kursumfrage auszufüllen! Dein Feedback hilft uns, unsere Trainingsmaterialien für alle zu verbessern.

[Zur Umfrage :material-arrow-right:](survey.md){ .md-button .md-button--primary }
