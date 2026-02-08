---
title: Hello Nextflow
hide:
  - toc
page_type: index_page
index_type: course
additional_information:
  technical_requirements: true
  learning_objectives:
    - Starten und Verwalten der Ausführung von Nextflow-Workflows
    - Finden und Interpretieren von Ausgaben (Ergebnissen) und Log-Dateien, die von Nextflow generiert werden
    - Beheben grundlegender Probleme
    - Erstellen eines einfachen mehrstufigen Workflows aus Nextflow-Kernkomponenten
    - Unterscheiden zwischen wesentlichen Arten von channel factories und Operatoren und deren effektive Nutzung in einem einfachen Workflow
    - Konfigurieren der Pipeline-Ausführung für gängige Rechenplattformen einschließlich HPC und Cloud
    - Anwenden von Best Practices für Reproduzierbarkeit, Portabilität und Code-Wiederverwendung, die Pipelines FAIR machen, einschließlich Code-Modularität und Software-Container
  audience_prerequisites:
    - "**Zielgruppe:** Dieser Kurs ist für Lernende konzipiert, die komplett neu bei Nextflow sind und eigene Pipelines entwickeln möchten."
    - "**Fähigkeiten:** Grundlegende Vertrautheit mit der Befehlszeile, grundlegenden Scripting-Konzepten und gängigen Dateiformaten wird vorausgesetzt."
    - "**Fachgebiet:** Die Übungen sind alle fachgebietsunabhängig, daher ist kein wissenschaftliches Vorwissen erforderlich."
  videos_playlist: https://youtube.com/playlist?list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&si=eF7cLR62goy-lc6n
---

# Hello Nextflow

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } KI-gestützte Übersetzung - [mehr erfahren & Verbesserungen vorschlagen](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

**Hello Nextflow ist eine praxisorientierte Einführung in den Aufbau reproduzierbarer und skalierbarer Datenanalyse-Workflows.**

In praktischen Übungen lernst du die Pipeline-Entwicklung mit Nextflow: Prozesse definieren, Pipelines verbinden, Dateien und Dependencies verwalten, Ausführung parallelisieren und Workflows in verschiedenen Umgebungen starten.

Nach diesem Kurs kannst du eigene Workflows mit Nextflow entwickeln und ausführen.

<!-- additional_information -->

## Kursübersicht

Dieser Kurs ist praxisorientiert konzipiert, mit zielgerichteten Übungen, die Informationen schrittweise einführen.

Du wirst eine einfache Nextflow-Pipeline entwickeln, die einige Texteingaben nimmt, einige Transformationsschritte ausführt und eine einzelne Textdatei mit einem ASCII-Bild einer Figur ausgibt, die den transformierten Text sagt.

### Lektionsplan

Um dich nicht mit Konzepten und Code zu überfordern, haben wir dies in sechs Teile aufgeteilt, die sich jeweils auf bestimmte Aspekte der Pipeline-Entwicklung mit Nextflow konzentrieren.

| Kurskapitel                                          | Zusammenfassung                                                                                                          | Geschätzte Dauer |
| ---------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------ | ---------------- |
| [Teil 1: Hello World](./01_hello_world.md)           | Grundlegende Komponenten und Prinzipien beim Zusammenstellen und Ausführen eines Nextflow-Workflows                      | 30 Min.          |
| [Teil 2: Hello Channels](./02_hello_channels.md)     | Verwendung von channels und Operatoren zur Verarbeitung von Eingaben und müheloser Parallelisierung                      | 45 Min.          |
| [Teil 3: Hello Workflow](./03_hello_workflow.md)     | Verwendung von channels zum Verketten mehrerer Schritte und Handhabung des Datentransfers zwischen Schritten             | 60 Min.          |
| [Teil 4: Hello Modules](./04_hello_modules.md)       | Anwendung von Code-Modularitätsprinzipien zur Erhöhung der Wiederverwendbarkeit und Verringerung des Wartungsaufwands    | 20 Min.          |
| [Teil 5: Hello Containers](./05_hello_containers.md) | Verwendung von Containern als Mechanismus zur Verwaltung von Software-Abhängigkeiten und Erhöhung der Reproduzierbarkeit | 60 Min.          |
| [Teil 6: Hello Config](./06_hello_config.md)         | Anpassung des Pipeline-Verhaltens und Optimierung der Nutzung in verschiedenen Rechenumgebungen                          | 60 Min.          |

Nach Abschluss dieses Kurses bist du bereit, reproduzierbare Workflows für deine eigenen Projekte zu entwickeln.

Bereit, den Kurs zu starten?

[Los geht's :material-arrow-right:](00_orientation.md){ .md-button .md-button--primary }

<!-- Clearfix for float -->
<div style="content: ''; clear: both; display: table;"></div>
