---
title: Nextflow for {DOMAIN}
hide:
  - toc
page_type: index_page
index_type: course
additional_information:
  technical_requirements: true
  learning_objectives:
    - Write a linear workflow to apply {METHOD} to a single sample
    - Handle accessory files such as {ACCESSORY_FILES} appropriately
    - Leverage Nextflow's dataflow paradigm to parallelize per-sample processing
    - Implement multi-sample aggregation using relevant channel operators
  audience_prerequisites:
    - "**Audience:** This course is designed for researchers in {DOMAIN} and related fields who want to develop or customize data analysis pipelines."
    - "**Skills:** Some familiarity with the command line, basic scripting concepts, and common {DOMAIN} file formats is assumed."
    - "**Prerequisites:** Foundational Nextflow concepts and tooling covered in [Hello Nextflow](../../hello_nextflow/)."
---

---
title: Nextflow für {DOMAIN}
hide:
  - toc
page_type: index_page
index_type: course
additional_information:
  technical_requirements: true
  learning_objectives:
    - Schreibe einen linearen Workflow, um {METHOD} auf eine einzelne Probe anzuwenden
    - Verarbeite zusätzliche Dateien wie {ACCESSORY_FILES} angemessen
    - Nutze Nextflows Dataflow-Paradigma, um die Verarbeitung pro Probe zu parallelisieren
    - Implementiere Multi-Sample-Aggregation mit relevanten Channel-Operatoren
  audience_prerequisites:
    - "**Zielgruppe:** Dieser Kurs richtet sich an Forscher\*innen im Bereich {DOMAIN} und verwandten Feldern, die Datenanalyse-Pipelines entwickeln oder anpassen möchten."
    - "**Fähigkeiten:** Grundkenntnisse der Kommandozeile, grundlegende Skripting-Konzepte und gängige {DOMAIN}-Dateiformate werden vorausgesetzt."
    - "**Voraussetzungen:** Grundlegende Nextflow-Konzepte und -Werkzeuge, die in [Hello Nextflow](../../hello_nextflow/) behandelt werden."
---

# Nextflow für {DOMAIN}

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } KI-gestützte Übersetzung - [mehr erfahren & Verbesserungen vorschlagen](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

**Ein praxisorientierter Kurs, der Nextflow auf einen realen {DOMAIN}-Anwendungsfall anwendet: {METHOD_SHORT_DESCRIPTION}.**

Dieser Kurs baut auf dem [Hello Nextflow](../../hello_nextflow/)-Einsteiger-Training auf und zeigt, wie du Nextflow im spezifischen Kontext der {DOMAIN}-Domäne verwendest.
Du wirst eine {METHOD}-Pipeline mit [{TOOL_A}]({TOOL_A_URL}) und [{TOOL_B}]({TOOL_B_URL}) implementieren.

<!-- additional_information -->

## Kursübersicht

Dieser Kurs ist praxisorientiert, mit zielgerichteten Übungen, die so strukturiert sind, dass Informationen schrittweise eingeführt werden.

Du beginnst damit, die Analyse-Tools manuell im Terminal auszuführen, um die Methodik zu verstehen. Anschließend baust du schrittweise eine Nextflow-Pipeline auf, die die Analyse automatisiert und skaliert.

### Lektionsplan

Wir haben den Kurs in drei Teile unterteilt, die sich jeweils auf spezifische Aspekte der Anwendung von Nextflow auf einen {DOMAIN}-Anwendungsfall konzentrieren.

| Kurskapitel                                               | Zusammenfassung                                                                                                                | Geschätzte Dauer |
| --------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------ | ---------------- |
| [Teil 1: Methodenübersicht](./01_method.md)               | Verstehen der {METHOD}-Methodik und manuelles Ausführen der Tools                                                              | 30 Min.          |
| [Teil 2: Einzelproben-Verarbeitung](./02_single_sample.md) | Erstellen einer Pipeline, die {PART2_SUMMARY}, und anschließende Skalierung auf mehrere Proben                                | 60 Min.          |
| [Teil 3: Multi-Sample-Aggregation](./03_multi_sample.md)  | Hinzufügen von Multi-Sample-{AGGREGATION_SUMMARY} mithilfe von Channel-Operatoren zur Aggregation von Proben-spezifischen Ausgaben | 45 Min.          |

Nach Abschluss dieses Kurses kannst du grundlegende Nextflow-Konzepte und -Werkzeuge auf einen typischen {DOMAIN}-Anwendungsfall anwenden.

Bereit, den Kurs zu beginnen?

[Los geht's :material-arrow-right:](00_orientation.md){ .md-button .md-button--primary }

<!-- Clearfix for float -->
<div style="content: ''; clear: both; display: table;"></div>
