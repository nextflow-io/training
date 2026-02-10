---
title: Nextflow for Genomics
hide:
  - toc
page_type: index_page
index_type: course
additional_information:
  technical_requirements: true
  learning_objectives:
    - Write a linear workflow to apply variant calling to a single sample
    - Handle accessory files such as index files and reference genome resources appropriately
    - Leverage Nextflow's dataflow paradigm to parallelize per-sample variant calling
    - Implement multi-sample joint calling using relevant channel operators
  audience_prerequisites:
    - "**Audience:** This course is designed for researchers in genomics and related fields who want to develop or customize data analysis pipelines."
    - "**Skills:** Some familiarity with the command line, basic scripting concepts, and common genomics file formats is assumed."
    - "**Prerequisites:** Foundational Nextflow concepts and tooling covered in [Hello Nextflow](../../hello_nextflow/)."
---

---

title: Nextflow für Genomik
hide:

- toc
  page_type: index_page
  index_type: course
  additional_information:
  technical_requirements: true
  learning_objectives: - Schreibe einen linearen Workflow, um Variantenaufruf auf eine einzelne Probe anzuwenden - Verwalte Hilfsdateien wie Indexdateien und Referenzgenom-Ressourcen angemessen - Nutze Nextflows Dataflow-Paradigma, um den Variantenaufruf pro Probe zu parallelisieren - Implementiere Multi-Sample Joint Calling mit relevanten Channel-Operatoren
  audience_prerequisites: - "**Zielgruppe:** Dieser Kurs richtet sich an Forscher\*innen in der Genomik und verwandten Bereichen, die Datenanalyse-Pipelines entwickeln oder anpassen möchten." - "**Fähigkeiten:** Grundkenntnisse der Kommandozeile, grundlegende Skripting-Konzepte und gängige Genomik-Dateiformate werden vorausgesetzt." - "**Voraussetzungen:** Grundlegende Nextflow-Konzepte und Werkzeuge, die in [Hello Nextflow](../../hello_nextflow/) behandelt werden."

---

# Nextflow für Genomik

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } KI-gestützte Übersetzung - [mehr erfahren & Verbesserungen vorschlagen](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

**Ein praxisorientierter Kurs, der Nextflow auf einen realen Genomik-Anwendungsfall anwendet: Variantenaufruf mit GATK.**

Dieser Kurs baut auf dem [Hello Nextflow](../../hello_nextflow/) Einsteiger-Training auf und zeigt, wie du Nextflow im spezifischen Kontext der Genomik einsetzen kannst.
Du wirst eine Variantenaufruf-Pipeline mit [GATK](https://gatk.broadinstitute.org/) (Genome Analysis Toolkit) implementieren, einem weit verbreiteten Softwarepaket zur Analyse von Hochdurchsatz-Sequenzierungsdaten.

<!-- additional_information -->

## Kursübersicht

Dieser Kurs ist praxisorientiert, mit zielgerichteten Übungen, die so strukturiert sind, dass Informationen schrittweise eingeführt werden.

Du beginnst damit, die Variantenaufruf-Tools manuell im Terminal auszuführen, um die Methodik zu verstehen. Anschließend baust du schrittweise eine Nextflow-Pipeline auf, die die Analyse automatisiert und skaliert.

### Lektionsplan

Wir haben dies in drei Teile unterteilt, die sich jeweils auf spezifische Aspekte der Anwendung von Nextflow auf einen Genomik-Anwendungsfall konzentrieren.

| Kurskapitel                                                             | Zusammenfassung                                                                                             | Geschätzte Dauer |
| ----------------------------------------------------------------------- | ----------------------------------------------------------------------------------------------------------- | ---------------- |
| [Teil 1: Methodenübersicht](./01_method.md)                             | Die Variantenaufruf-Methodik verstehen und die Tools manuell ausführen                                      | 30 Min.          |
| [Teil 2: Variantenaufruf pro Probe](./02_per_sample_variant_calling.md) | Eine Pipeline erstellen, die BAM-Dateien indiziert und Varianten aufruft, dann auf mehrere Proben skalieren | 60 Min.          |
| [Teil 3: Joint Calling auf einer Kohorte](./03_joint_calling.md)        | Multi-Sample Joint Genotyping hinzufügen mit Channel-Operatoren, um Ausgaben pro Probe zu aggregieren       | 45 Min.          |

Am Ende dieses Kurses kannst du grundlegende Nextflow-Konzepte und Werkzeuge auf einen typischen Genomik-Anwendungsfall anwenden.

Bereit, den Kurs zu beginnen?

[Los geht's :material-arrow-right:](00_orientation.md){ .md-button .md-button--primary }

<!-- Clearfix for float -->
<div style="content: ''; clear: both; display: table;"></div>
