---
title: Nextflow für Genomik
hide:
  - toc
---

# Nextflow für Genomik

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } KI-gestützte Übersetzung - [mehr erfahren & Verbesserungen vorschlagen](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Dieser Trainingskurs richtet sich an Forschende im Bereich Genomik und verwandten Feldern, die daran interessiert sind, Datenanalyse-Pipelines zu entwickeln oder anzupassen.
Er baut auf dem [Hello Nextflow](../../hello_nextflow/) Einsteiger-Training auf und zeigt, wie du Nextflow im spezifischen Kontext der Genomik-Domäne nutzen kannst.

Konkret demonstriert dieser Kurs, wie du eine einfache Variant-Calling-Pipeline mit [GATK](https://gatk.broadinstitute.org/) (Genome Analysis Toolkit) implementierst, einem weit verbreiteten Softwarepaket zur Analyse von High-Throughput-Sequenzierungsdaten.

Los geht's! Klicke auf die Schaltfläche „Open in GitHub Codespaces" unten, um die Trainingsumgebung zu starten (am besten in einem separaten Tab), und lies dann weiter, während sie lädt.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

## Lernziele

Durch die Bearbeitung dieses Kurses lernst du, wie du grundlegende Nextflow-Konzepte und -Werkzeuge auf einen typischen Genomik-Anwendungsfall anwendest.

Am Ende dieses Workshops kannst du:

- Einen linearen Workflow zu schreiben, um Variant Calling auf eine einzelne Probe anzuwenden
- Zusatzdateien wie Index-Dateien und Referenzgenom-Ressourcen angemessen zu handhaben
- Das Dataflow-Paradigma von Nextflow zu nutzen, um Variant Calling pro Probe zu parallelisieren
- Multi-Sample Variant Calling mit relevanten Channel-Operatoren zu implementieren
- Tests pro Schritt und Ende-zu-Ende-Pipeline-Tests zu implementieren, die genomik-spezifische Besonderheiten angemessen behandeln

<!-- TODO for future expansion: add metadata/samplesheet handling -->

## Voraussetzungen

Der Kurs setzt ein Mindestmaß an Vertrautheit mit Folgendem voraus:

- Werkzeuge und Dateiformate, die in diesem wissenschaftlichen Bereich häufig verwendet werden
- Erfahrung mit der Kommandozeile
- Grundlegende Nextflow-Konzepte und -Werkzeuge, die im [Hello Nextflow](../../hello_nextflow/) Einsteiger-Training behandelt werden.

Für technische Anforderungen und Umgebungs-Setup siehe den Mini-Kurs [Environment Setup](../../envsetup/).
