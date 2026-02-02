---
title: Nextflow für RNAseq
hide:
  - toc
---

# Nextflow für RNAseq

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } KI-gestützte Übersetzung - [mehr erfahren & Verbesserungen vorschlagen](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Dieser Trainingskurs richtet sich an Forschende im Bereich Transkriptomik und verwandter Fachgebiete, die daran interessiert sind, Datenanalyse-Pipelines zu entwickeln oder anzupassen.
Er baut auf dem [Hello Nextflow](../../hello_nextflow/) Einsteiger-Training auf und demonstriert, wie man Nextflow im spezifischen Kontext der Bulk-RNAseq-Analyse einsetzt.

Konkret zeigt dieser Kurs, wie man eine einfache Bulk-RNAseq-Verarbeitungs-Pipeline implementiert, um Adaptersequenzen zu trimmen, die Reads an eine Genomreferenz zu alignieren und Qualitätskontrollen (QC) in mehreren Schritten durchzuführen.

Lass uns anfangen! Klicke auf den Button "Open in GitHub Codespaces" unten, um die Trainingsumgebung zu starten (vorzugsweise in einem separaten Tab), und lies dann weiter, während sie lädt.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

## Lernziele

Durch die Bearbeitung dieses Kurses lernst du, wie man grundlegende Nextflow-Konzepte und Tools auf einen typischen RNAseq-Anwendungsfall anwendet.

Am Ende dieses Workshops kannst du:

- Einen linearen Workflow zu schreiben, um grundlegende RNAseq-Verarbeitungs- und QC-Methoden anzuwenden
- Domänenspezifische Dateien wie FASTQ und Genomreferenzen angemessen zu handhaben
- Single-End- und Paired-End-Sequenzierungsdaten zu verarbeiten
- Nextflows Dataflow-Paradigma zu nutzen, um die RNAseq-Verarbeitung pro Probe zu parallelisieren
- QC-Reports über mehrere Schritte und Proben hinweg mit relevanten Channel-Operatoren zu aggregieren

<!-- TODO
- Configure pipeline execution and manage and optimize resource allocations
- Implement per-step and end-to-end pipeline tests that handle RNAseq-specific idiosyncrasies appropriately
-->
<!-- TODO for future expansion: add metadata/samplesheet handling -->

## Voraussetzungen

Der Kurs setzt eine minimale Vertrautheit mit Folgendem voraus:

- Tools und Dateiformate, die in diesem wissenschaftlichen Bereich häufig verwendet werden
- Erfahrung mit der Kommandozeile
- Grundlegende Nextflow-Konzepte und Tools, die im [Hello Nextflow](../../hello_nextflow/) Einsteiger-Training behandelt werden.

Für technische Anforderungen und Umgebungseinrichtung siehe den Mini-Kurs [Environment Setup](../../envsetup/).
