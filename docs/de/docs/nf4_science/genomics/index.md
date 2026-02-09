---
title: Nextflow für Genomik
hide:
  - toc
page_type: index_page
index_type: course
additional_information:
  technical_requirements: true
  learning_objectives:
    - Einen linearen Workflow zu schreiben, um Variant Calling auf eine einzelne Probe anzuwenden
    - Zusatzdateien wie Index-Dateien und Referenzgenom-Ressourcen angemessen zu handhaben
    - Das Dataflow-Paradigm von Nextflow zu nutzen, um Variant Calling pro Probe zu parallelisieren
    - Multi-Sample Joint Calling mit relevanten Channel-Operatoren zu implementieren
  audience_prerequisites:
    - "**Zielgruppe:** Dieser Kurs richtet sich an Forschende im Bereich Genomik und verwandten Feldern, die Datenanalyse-Pipelines entwickeln oder anpassen möchten."
    - "**Kenntnisse:** Grundlegende Vertrautheit mit der Kommandozeile, grundlegenden Scripting-Konzepten und gängigen Dateiformaten der Genomik wird vorausgesetzt."
    - "**Voraussetzungen:** Grundlegende Nextflow-Konzepte und -Werkzeuge, die im Kurs [Hello Nextflow](../../hello_nextflow/) behandelt werden."
---

# Nextflow für Genomik

**Ein praxisorientierter Kurs, der Nextflow auf einen realen Genomik-Anwendungsfall anwendet: Variant Calling mit GATK.**

Dieser Kurs baut auf dem [Hello Nextflow](../../hello_nextflow/) Einsteiger-Training auf und zeigt, wie du Nextflow im spezifischen Kontext der Genomik-Domäne nutzen kannst.
Du wirst eine Variant-Calling-Pipeline mit [GATK](https://gatk.broadinstitute.org/) (Genome Analysis Toolkit) implementieren, einem weit verbreiteten Softwarepaket zur Analyse von High-Throughput-Sequenzierungsdaten.

<!-- additional_information -->

## Kursübersicht

Dieser Kurs ist praxisorientiert, mit zielgerichteten Übungen, die darauf ausgelegt sind, Informationen schrittweise einzuführen.

Du beginnst damit, die Variant-Calling-Tools manuell im Terminal auszuführen, um die Methodik zu verstehen, und baust dann schrittweise eine Nextflow-Pipeline auf, die die Analyse automatisiert und skaliert.

### Lektionsplan

Wir haben dies in drei Teile unterteilt, die sich jeweils auf spezifische Aspekte der Anwendung von Nextflow auf einen Genomik-Anwendungsfall konzentrieren.

| Kurskapitel                                                             | Zusammenfassung                                                                                            | Geschätzte Dauer |
| ----------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------- | ---------------- |
| [Teil 1: Methodenübersicht](./01_method.md)                             | Verständnis der Variant-Calling-Methodik und manuelle Ausführung der Tools                                 | 30 Min.          |
| [Teil 2: Variant Calling pro Probe](./02_per_sample_variant_calling.md) | Aufbau einer Pipeline, die BAM-Dateien indexiert und Varianten aufruft, dann Skalierung auf mehrere Proben | 60 Min.          |
| [Teil 3: Joint Calling auf einer Kohorte](./03_joint_calling.md)        | Hinzufügen von Multi-Sample Joint Genotyping mit Channel-Operatoren zur Aggregation von Ausgaben pro Probe | 45 Min.          |

Am Ende dieses Kurses kannst du grundlegende Nextflow-Konzepte und -Werkzeuge auf einen typischen Genomik-Anwendungsfall anwenden.

Bereit, den Kurs zu starten?

[Los geht's :material-arrow-right:](00_orientation.md){ .md-button .md-button--primary }

<!-- Clearfix for float -->
<div style="content: ''; clear: both; display: table;"></div>
