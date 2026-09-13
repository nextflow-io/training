---
title: Nextflow Run
hide:
  - toc
page_type: index_page
index_type: course
additional_information:
  technical_requirements: true
  learning_objectives:
    - Nextflow-Pipelines von der Befehlszeile starten und verwalten
    - Verstehen, wie Kanäle und Operatoren effiziente Workflows mit mehreren Eingaben und Schritten ermöglichen
    - Container verwenden, um Software-Abhängigkeiten zu verwalten und Reproduzierbarkeit sicherzustellen
    - Pipeline-Ausführung und Ausgaben konfigurieren
    - Ausführungsberichte erstellen, den Verlauf vergangener Ausführungen einsehen und alte Work-Verzeichnisse bereinigen
    - Pipelines direkt aus Remote-Repositories wie GitHub ausführen
  audience_prerequisites:
    - "**Zielgruppe:** Dieser Kurs richtet sich an Lernende, die völlig neu bei Nextflow sind und bestehende Pipelines ausführen möchten."
    - "**Kenntnisse:** Grundlegende Vertrautheit mit der Befehlszeile, grundlegenden Scripting-Konzepten und gängigen Dateiformaten wird vorausgesetzt."
    - "**Fachgebiet:** Die Übungen sind alle fachunabhängig, sodass keine wissenschaftlichen Vorkenntnisse erforderlich sind."
---

# Nextflow Run

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } KI-gestützte Übersetzung - [mehr erfahren & Verbesserungen vorschlagen](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

**Nextflow Run ist eine praktische Einführung in die Ausführung reproduzierbarer und skalierbarer Datenanalyse-Workflows.**

In einer Reihe zielgerichteter Übungen lernst du die Grundlagen zum Starten und Verwalten von Nextflow-Pipelines, verstehst, wie Kanäle und Operatoren die parallele Verarbeitung mehrerer Eingaben ermöglichen, und verwendest Container zur Verwaltung von Software-Abhängigkeiten.

Nach diesem Kurs hast du die Fähigkeiten und das Vertrauen, um Workflows mit Nextflow auszuführen.

<!-- additional_information -->

## Kursübersicht

Dieser Kurs ist praxisorientiert, mit zielgerichteten Übungen, die Informationen schrittweise einführen.

Du wirst mehrere Versionen einer Nextflow-Pipeline ausführen, die Texteingaben verarbeitet. Du beginnst mit einer einfachen Version, die aus einem einzelnen Schritt besteht, und gelangst schließlich zu einer mehrstufigen Version, die eine CSV-Datei mit Eingaben nimmt, einige Transformationsschritte ausführt und eine einzelne Textdatei ausgibt, die ASCII-Art enthält, die von einem containerisierten Tool erzeugt wird.

Dieser Kurs konzentriert sich auf das Ausführen von Pipelines (benannt nach dem Kernbefehl `nextflow run`).
Wenn du eine Einführung in die Entwicklung von Nextflow-Pipelines suchst, siehe [Hello Nextflow](../hello_nextflow/index.md).

!!! note "Hinweis"

    Suchst du die frühere Version dieses Kurses? Sie wurde durch die Version auf dieser Seite ersetzt, ist aber noch im [Release 3.6.1](https://training.nextflow.io/3.6.1/nextflow_run/) der Trainingsseite verfügbar.

### Lehrplan

| Kurskapitel                                                          | Zusammenfassung                                                                                                     | Geschätzte Dauer |
| -------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------- | ---------------- |
| [Teil 1: Nextflow ausführen](./01_run_nextflow.md)                   | Nextflow-Pipelines starten und verwalten sowie grundlegende Workflow-Mechanismen verstehen                          | 25 Min           |
| [Teil 2: Pipeline konfigurieren](./02_configure_pipeline.md)         | Pipeline-Ausführung und Ausgaben mit `nextflow.config` konfigurieren                                                | 20 Min           |
| [Teil 3: Workflow-Ausführungen verwalten](./03_manage_executions.md) | Ausführungsberichte erstellen, den Verlauf vergangener Ausführungen einsehen und alte Work-Verzeichnisse bereinigen | 10 Min           |
| [Teil 4: Remote-Pipelines ausführen](./04_remote_repositories.md)    | Eine Pipeline direkt von GitHub ausführen und auf eine bestimmte Revision festlegen                                 | 10 Min           |

Nach Abschluss dieses Kurses bist du gut vorbereitet, um reproduzierbare Workflows für deine eigenen Projekte auszuführen.

Bereit, den Kurs zu beginnen?

[Jetzt lernen :material-arrow-right:](00_orientation.md){ .md-button .md-button--primary }
