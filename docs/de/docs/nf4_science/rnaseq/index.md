---
title: RNAseq
hide:
  - toc
page_type: index_page
index_type: course
additional_information:
  technical_requirements: true
  learning_objectives:
    - Einen linearen Workflow zu schreiben, um grundlegende RNAseq-Verarbeitungs- und QC-Methoden anzuwenden
    - Domänenspezifische Dateien wie FASTQ und Genomreferenzen angemessen zu handhaben
    - Single-End- und Paired-End-Sequenzierungsdaten zu verarbeiten
    - Nextflows Dataflow-Paradigma zu nutzen, um die RNAseq-Verarbeitung pro Probe zu parallelisieren
    - QC-Reports über mehrere Schritte und Proben hinweg mit relevanten Channel-Operatoren zu aggregieren
  audience_prerequisites:
    - "**Zielgruppe:** Dieser Kurs richtet sich an Forschende im Bereich Transkriptomik und verwandter Fachgebiete, die Datenanalyse-Pipelines entwickeln oder anpassen möchten."
    - "**Kenntnisse:** Grundlegende Vertrautheit mit der Kommandozeile, einfachen Scripting-Konzepten und gängigen RNAseq-Dateiformaten wird vorausgesetzt."
    - "**Voraussetzungen:** Grundlegende Nextflow-Konzepte und Tools, die im [Hello Nextflow](../../hello_nextflow/) Einsteiger-Training behandelt werden."
---

# Nextflow für RNAseq

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } KI-gestützte Übersetzung - [mehr erfahren & Verbesserungen vorschlagen](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

**Ein praxisorientierter Kurs, der Nextflow auf einen realen Transkriptomik-Anwendungsfall anwendet: Bulk-RNAseq-Verarbeitung mit Trim Galore, HISAT2 und FastQC.**

Dieser Kurs baut auf dem [Hello Nextflow](../../hello_nextflow/) Einsteiger-Training auf und zeigt, wie man Nextflow im spezifischen Kontext der Bulk-RNAseq-Analyse einsetzt.
Du implementierst eine Verarbeitungs-Pipeline, die Adaptersequenzen trimmt, Reads an eine Genomreferenz aligniert und Qualitätskontrollen (QC) in mehreren Schritten durchführt.

<!-- additional_information -->

## Kursübersicht

Dieser Kurs ist praxisorientiert, mit zielgerichteten Übungen, die so strukturiert sind, dass Informationen schrittweise eingeführt werden.

Du beginnst damit, die Verarbeitungstools manuell im Terminal auszuführen, um die Methodik zu verstehen. Dann baust du schrittweise eine Nextflow-Pipeline auf, die die Analyse automatisiert und skaliert.

### Lektionsplan

Wir haben den Kurs in drei Teile unterteilt, die sich jeweils auf spezifische Aspekte der Anwendung von Nextflow auf einen RNAseq-Anwendungsfall konzentrieren.

| Kurskapitel                                                           | Zusammenfassung                                                                                                        | Geschätzte Dauer |
| --------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------- | ---------------- |
| [Teil 1: Methodenübersicht](./01_method.md)                           | Die RNAseq-Verarbeitungsmethodik verstehen und die Tools manuell ausführen                                             | 30 Min.          |
| [Teil 2: Einzelproben-Implementierung](./02_single-sample.md)         | Eine Pipeline erstellen, die eine einzelne Probe trimmt, aligniert und QC durchführt, dann auf mehrere Proben skaliert | 60 Min.          |
| [Teil 3: Mehrproben-Paired-End-Implementierung](./03_multi-sample.md) | Die Pipeline erweitern, um Paired-End-Daten zu verarbeiten und QC-Reports über Proben hinweg zu aggregieren            | 45 Min.          |

Am Ende dieses Kurses kannst du grundlegende Nextflow-Konzepte und Tools auf einen typischen RNAseq-Anwendungsfall anwenden.

Bereit, den Kurs zu beginnen?

[Los geht's :material-arrow-right:](00_orientation.md){ .md-button .md-button--primary }

<!-- Clearfix for float -->
<div style="content: ''; clear: both; display: table;"></div>
