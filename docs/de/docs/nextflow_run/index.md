---
title: Nextflow Run
hide:
  - toc
page_type: index_page
index_type: course
additional_information:
  technical_requirements: true
  learning_objectives:
    - Nextflow-Workflows starten und ihre Ausführung verwalten
    - Ausgaben (Ergebnisse) und Protokolldateien finden und interpretieren
    - Kernkomponenten von Nextflow in einem einfachen mehrstufigen Workflow erkennen
    - Pipeline-Ausführung für gängige Rechenplattformen einschließlich HPC und Cloud konfigurieren
    - Best Practices für Reproduzierbarkeit, Portabilität und Code-Wiederverwendung zusammenfassen, die Pipelines FAIR machen, einschließlich Code-Modularität und Software-Container
  audience_prerequisites:
    - "**Zielgruppe:** Dieser Kurs richtet sich an Lernende, die völlig neu bei Nextflow sind und bestehende Pipelines ausführen möchten."
    - "**Kenntnisse:** Grundlegende Vertrautheit mit der Befehlszeile, grundlegenden Scripting-Konzepten und gängigen Dateiformaten wird vorausgesetzt."
    - "**Fachgebiet:** Die Übungen sind alle fachunabhängig, sodass keine wissenschaftlichen Vorkenntnisse erforderlich sind."
---

# Nextflow Run

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } KI-gestützte Übersetzung - [mehr erfahren & Verbesserungen vorschlagen](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

**Nextflow Run ist eine praktische Einführung in die Ausführung reproduzierbarer und skalierbarer Datenanalyse-Workflows.**

Durch praktische Beispiele und geführte Übungen lernst du die Grundlagen der Verwendung von Nextflow, einschließlich der Ausführung von Pipelines, der Verwaltung von Dateien und Software-Abhängigkeiten, der mühelosen Parallelisierung der Ausführung und der Ausführung von Workflows in verschiedenen Rechenumgebungen.

Nach diesem Kurs kannst du Workflows mit Nextflow ausführen.

<!-- additional_information -->

## Kursübersicht

### Was du tun wirst

Dieser Kurs ist praxisorientiert, mit zielgerichteten Übungen, die Informationen schrittweise einführen.

Du wirst mehrere Versionen einer Nextflow-Pipeline ausführen, die Texteingaben verarbeitet.
Du beginnst mit einer einfachen Version, die aus einem einzelnen Schritt besteht, und gelangst schließlich zu einer mehrstufigen Version, die eine CSV-Datei mit tabellarischen Texteingaben nimmt, einige Transformationsschritte ausführt und eine einzelne Textdatei ausgibt, die ein ASCII-Bild eines Charakters enthält, der den transformierten Text sagt.

Dieser Kurs konzentriert sich auf das Ausführen von Pipelines (benannt nach dem Kernbefehl `nextflow run`).
Wenn du eine Einführung in die Entwicklung von Nextflow-Pipelines suchst, siehe [Hello Nextflow](../hello_nextflow/index.md).

### Lehrplan

Wir haben dies in drei Teile unterteilt, die sich jeweils auf bestimmte Aspekte des Ausführens und Verwaltens von Pipelines konzentrieren, die in Nextflow geschrieben sind.

| Kurskapitel                                           | Zusammenfassung                                                                                                 | Geschätzte Dauer |
| ----------------------------------------------------- | --------------------------------------------------------------------------------------------------------------- | ---------------- |
| [Teil 1: Grundlegende Operationen](./01_basics.md)    | Starten und Verwalten der Ausführung eines einfachen Workflows                                                  | 30 Min           |
| [Teil 2: Echte Pipelines ausführen](./02_pipeline.md) | Komplexe Eingaben verarbeiten, mehrstufige Workflows ausführen, Container verwenden und mühelos parallelisieren | 60 Min           |
| [Teil 3: Konfiguration](./03_config.md)               | Pipeline-Verhalten anpassen und Nutzung in verschiedenen Rechenumgebungen optimieren                            | 60 Min           |

Nach Abschluss dieses Kurses bist du bereit, reproduzierbare Workflows für deine eigenen Projekte auszuführen.

Bereit, den Kurs zu beginnen?

[Jetzt lernen :material-arrow-right:](00_orientation.md){ .md-button .md-button--primary }
