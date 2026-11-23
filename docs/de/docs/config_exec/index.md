---
title: Configure Execution
hide:
  - toc
page_type: index_page
index_type: course
additional_information:
  technical_requirements: true
  learning_objectives:
    - Zwischen Docker und Conda als Software-Packaging-Technologie wechseln
    - Eine Ausführungsplattform auswählen und verstehen, wie Nextflow die Aufgabenausführung daran anpasst
    - Rechenressourcen steuern und fehlgeschlagene Aufgaben automatisch wiederholen
    - Profile definieren und kombinieren, um zwischen vordefinierten Konfigurationen zu wechseln
  audience_prerequisites:
    - "**Zielgruppe:** Dieser Kurs richtet sich an Lernende, die bereits wissen, wie man lokale Nextflow-Pipelines startet, und die Ausführung tiefer konfigurieren möchten."
    - "**Kenntnisse:** Grundlegende Kenntnisse der Kommandozeile werden vorausgesetzt."
    - "**Kurse:** Du musst [Nextflow Run](../nextflow_run/index.md) abgeschlossen haben oder anderweitig mit dem Ausführen einer lokalen Pipeline mit `nextflow run` vertraut sein."
---

# Configure Execution

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } KI-gestützte Übersetzung - [mehr erfahren & Verbesserungen vorschlagen](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

**Configure Execution ist eine praxisorientierte Einführung in die Anpassung der Nextflow-Pipeline-Ausführung an verschiedene Rechenumgebungen.**

In zielorientierten Übungen lernst du, wie du die Software-Packaging-Technologie wechselst, eine Ausführungsplattform auswählst, Rechenressourcen und Wiederholungsversuche steuerst und Konfigurationen in umschaltbare Profile bündelst.

Du wirst die Fähigkeiten und das Vertrauen mitnehmen, um die Ausführung von Nextflow-Pipelines wie ein Profi zu konfigurieren.

<!-- additional_information -->

## Kursübersicht

Dieser Kurs ist praxisorientiert und baut auf den Kenntnissen auf, die in [Nextflow Run](../nextflow_run/index.md) vermittelt werden.

Du nimmst dieselbe mehrstufige Pipeline aus diesem Kurs und passt ihre Konfiguration schrittweise an verschiedene Rechenumgebungen an. Anschließend bündelst du alles in Profile, zwischen denen du zur Laufzeit wechseln kannst.

### Lehrplan

| Kurskapitel                                                                    | Zusammenfassung                                                                          | Geschätzte Dauer |
| ------------------------------------------------------------------------------ | ---------------------------------------------------------------------------------------- | ---------------- |
| [Teil 1: An deine Rechenumgebung anpassen](./01_packaging_and_execution.md)    | Software-Packaging-Technologie wechseln und eine Ausführungsplattform auswählen          | 20 Min.          |
| [Teil 2: Rechenressourcen und Fehler verwalten](./02_resources_and_retries.md) | Ressourcenzuweisung steuern und fehlgeschlagene Aufgaben automatisch wiederholen         | 15 Min.          |
| [Teil 3: Profile zum Wechseln von Konfigurationen nutzen](./03_profiles.md)    | Profile definieren und kombinieren sowie die vollständig aufgelöste Konfiguration prüfen | 15 Min.          |

Am Ende dieses Kurses kannst du Nextflow-Pipelines sicher für verschiedene Rechenumgebungen konfigurieren und mit minimalem Aufwand zwischen ihnen wechseln.

Bereit, den Kurs zu starten?

[Jetzt lernen :material-arrow-right:](00_orientation.md){ .md-button .md-button--primary }
