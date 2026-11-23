---
title: Skalieren mit Seqera
hide:
  - toc
page_type: index_page
index_type: course
additional_information:
  technical_requirements: true
  learning_objectives:
    - Registriere dich bei Seqera Platform und erkunde die Community Showcase
    - Füge eine Pipeline zu einem Workspace hinzu und starte sie über die Web-Oberfläche
    - Authentifiziere dich und starte Pipelines über die Befehlszeile mit der `tw` CLI
    - Registriere eine auf GitHub gehostete Pipeline und starte sie auf beide Arten
  audience_prerequisites:
    - "**Zielgruppe:** Dieser Kurs richtet sich an Lernende, die Nextflow-Pipelines mit Seqera Platform in großem Maßstab ausführen möchten."
    - "**Kenntnisse:** Grundlegende Erfahrung mit dem Ausführen von nf-core-Pipelines über die Befehlszeile wird vorausgesetzt."
    - "**Kurse:** Du musst [Nextflow Run](../nextflow_run/index.md) und [Use nf-core](../nfcore_use/index.md) abgeschlossen haben oder anderweitig mit dem Ausführen lokaler und `nf-core/rnaseq`-Pipelines vertraut sein."
---

# Skalieren mit Seqera

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } KI-gestützte Übersetzung - [mehr erfahren & Verbesserungen vorschlagen](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

**Skalieren mit Seqera ist eine praxisorientierte Einführung in das Starten und Überwachen von Nextflow-Pipelines mit Seqera Platform.**

Anhand praktischer Beispiele richtest du den Zugang zu Seqera Platform ein, startest eine produktionsreife Pipeline sowohl über die Web-Oberfläche als auch über die Befehlszeile und fügst eine neue Pipeline zu deinem Workspace hinzu.

Du wirst die Fähigkeiten und das Vertrauen mitnehmen, um deine eigenen Pipelines auf Seqera Platform auszuführen und zu überwachen.

<!-- additional_information -->

## Kursübersicht

Dieser Kurs ist praxisorientiert und baut auf den Pipelines auf, die du bereits in [Use nf-core](../nfcore_use/index.md) ausgeführt hast.

Du beginnst damit, dich bei Seqera Platform zu registrieren und `nf-core/rnaseq`, eine produktionsreife Pipeline, über die Web-Oberfläche zu starten.
Dann wechselst du zum `tw`-Befehlszeilentool, um dasselbe über ein Terminal zu tun, und registrierst abschließend eine neue Pipeline, `nf-core/demo`, und startest sie auf beide Arten.

### Lehrplan

| Kurskapitel                                                                  | Zusammenfassung                                                                                                    | Geschätzte Dauer |
| ---------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------ | ---------------- |
| [Teil 1: Pipelines über die Web-Oberfläche starten](./01_run_with_seqera.md) | Seqera Platform-Zugang einrichten und eine produktionsreife Pipeline über die Web-Oberfläche starten               | 20 Min.          |
| [Teil 2: Pipelines über die Befehlszeile starten](./02_launch_from_cli.md)   | Die `tw` CLI authentifizieren, eine gespeicherte Pipeline starten und eine neue Pipeline über die CLI registrieren | 25 Min.          |

Am Ende dieses Kurses kannst du Nextflow-Pipelines auf Seqera Platform sicher starten und überwachen – egal ob du lieber die Web-Oberfläche oder die Befehlszeile verwendest.

Bereit, den Kurs zu starten?

[Jetzt loslegen :material-arrow-right:](00_orientation.md){ .md-button .md-button--primary }
