---
title: nf-core verwenden
hide:
  - toc
page_type: index_page
index_type: course
additional_information:
  technical_requirements: true
  learning_objectives:
    - nf-core Community-Pipelines finden, abrufen und ausführen
    - Pipeline-Ausführung mit Parametern und Konfigurationsdateien anpassen
    - Verstehen, wie nf-core-Pipelines Parameter und Eingabedaten validieren
    - Eine produktionsreife Pipeline (nf-core/rnaseq) ausführen und ihre Standard-Ressourcenzuweisungen überschreiben
  audience_prerequisites:
    - "**Zielgruppe:** Dieser Kurs richtet sich an Lernende, die bereits wissen, wie man lokale Nextflow-Pipelines ausführt, neu bei nf-core sind und bestehende Community-Pipelines nutzen möchten."
    - "**Kenntnisse:** Grundlegende Vertrautheit mit der Kommandozeile, einfachen Scripting-Konzepten und gängigen Dateiformaten wird vorausgesetzt."
    - "**Kurse:** [Nextflow Run](../nextflow_run/index.md) muss abgeschlossen sein oder du solltest anderweitig mit dem Ausführen einer lokalen Pipeline mit `nextflow run` vertraut sein."
    - "**Fachgebiet:** Die Übungen verwenden Bioinformatik-Pipelines, aber kein wissenschaftliches Vorwissen ist erforderlich."
---

# nf-core verwenden

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } KI-gestützte Übersetzung - [mehr erfahren & Verbesserungen vorschlagen](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

**„nf-core verwenden" ist eine praxisorientierte Einführung in das Finden, Ausführen und Konfigurieren von nf-core Community-Pipelines.**

Anhand praktischer Beispiele und geführter Übungen lernst du, nf-core-Pipelines zu finden und abzurufen, sie mit ihren integrierten Test-Profilen auszuführen und ihre Ausführung über Parameter und Konfigurationsdateien anzupassen.

Du wirst die Fähigkeiten und das Vertrauen mitnehmen, um nf-core-Pipelines für deine eigenen Analysen einzusetzen.

<!-- additional_information -->

## Kursübersicht

Dieser Kurs ist praxisorientiert, mit zielgerichteten Übungen, die Informationen schrittweise einführen.

Du beginnst mit `nf-core/demo`, einer minimalen Pipeline, die vom nf-core-Projekt zu Trainingszwecken gepflegt wird, und wendest das Gelernte dann auf `nf-core/rnaseq` an – eine weit verbreitete Produktions-Pipeline für die Bulk-RNA-Sequenzierungsanalyse.

Dieser Kurs konzentriert sich auf das Ausführen von Pipelines.
Wenn du eine Einführung in die Entwicklung nf-core-kompatibler Pipelines suchst, schau dir [Build with nf-core](../nfcore_build/index.md) an.

### Lehrplan

| Kurskapitel                                                                    | Zusammenfassung                                                                                  | Geschätzte Dauer |
| ------------------------------------------------------------------------------ | ------------------------------------------------------------------------------------------------ | ---------------- |
| [Teil 1: Eine Demo-Pipeline ausführen](./01_run_demo.md)                       | Eine nf-core-Pipeline finden und abrufen und sie mit ihrem Test-Profil ausführen                 | 20 Min.          |
| [Teil 2: Pipeline-Ausführung konfigurieren](./02_configure_execution.md)       | Parameter setzen, Validierung verstehen und Ressourcenzuweisung sowie Tool-Argumente anpassen    | 20 Min.          |
| [Teil 3: Eine Produktions-Pipeline ausführen](./03_run_production_pipeline.md) | nf-core/rnaseq herunterladen und ausführen und ihre Standard-Ressourcenzuweisungen überschreiben | 20 Min.          |

Am Ende dieses Kurses kannst du die Vielzahl an Community-Pipelines des nf-core-Projekts für dich nutzen.

Bereit, den Kurs zu starten?

[Jetzt loslegen :material-arrow-right:](00_orientation.md){ .md-button .md-button--primary }
