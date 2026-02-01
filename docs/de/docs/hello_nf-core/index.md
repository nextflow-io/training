---
title: Hello nf-core
hide:
  - toc
page_type: index_page
index_type: course
additional_information:
  technical_requirements: true
  learning_objectives:
    - nf-core Pipelines abrufen, ausführen und verwalten
    - Die Code-Struktur und Projektorganisation von nf-core Pipelines beschreiben
    - Eine grundlegende nf-core kompatible Pipeline aus einem Template erstellen
    - Einen einfachen Nextflow Workflow auf nf-core Standards upgraden
    - nf-core Module zu einer nf-core kompatiblen Pipeline hinzufügen
    - Eigene Module zu nf-core beisteuern
    - Eingaben und Parameter mit nf-core Tooling validieren
  audience_prerequisites:
    - "**Zielgruppe:** Dieser Kurs richtet sich an Lernende, die bereits mit den Grundlagen von Nextflow vertraut sind und lernen möchten, nf-core Ressourcen und Best Practices zu nutzen."
    - "**Fähigkeiten:** Vertrautheit mit der Kommandozeile, grundlegende Scripting-Konzepte und gängige Dateiformate werden vorausgesetzt."
    - "**Kurse:** Der Kurs [Hello Nextflow](../hello_nextflow/index.md) oder gleichwertige Kenntnisse müssen abgeschlossen sein."
    - "**Fachgebiet:** Die Übungen sind alle fachgebietsneutral, daher ist kein spezifisches wissenschaftliches Vorwissen erforderlich."
---

# Hello nf-core

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } KI-gestützte Übersetzung - [mehr erfahren & Verbesserungen vorschlagen](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

**Hello nf-core ist eine praktische Einführung in die Nutzung von nf-core Ressourcen und Best Practices.**

![nf-core logo](./img/nf-core-logo.png)

Durch praktische Beispiele und angeleitete Übungen lernst du, nf-core kompatible Module und Pipelines zu nutzen und zu entwickeln sowie nf-core Tooling effektiv einzusetzen.

Nach diesem Kurs kannst du Pipelines nach nf-core Best Practices entwickeln.

<!-- additional_information -->

## Kursübersicht

Dieser Kurs ist darauf ausgelegt, praktisch zu sein, mit zielorientierten Übungen, die so strukturiert sind, dass Informationen schrittweise eingeführt werden.

Du wirst in [**nf-core**](https://nf-co.re/) eingeführt, eine Community-Initiative zur Entwicklung und Pflege einer kuratierten Sammlung wissenschaftlicher Pipelines, die mit Nextflow erstellt wurden, sowie relevantes Tooling und Richtlinien, die offene Entwicklung, Testing und Peer Review fördern ([Nat Biotechnol 38, 276–278 (2020)](https://www.nature.com/articles/s41587-020-0439-x), [Genome Biol 26, 228 (2025)](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-025-03673-9)).

nf-core Pipelines sind modular, skalierbar und portabel. So lassen sie sich leicht anpassen und mit eigenen Daten ausführen. Die Best-Practice-Richtlinien sorgen für robuste, gut dokumentierte und validierte Pipelines. Das verbessert die Reproduzierbarkeit wissenschaftlicher Analysen.

Wir werden in diesem Kurs nicht alles abdecken, was es über nf-core Pipelines zu wissen gibt, denn nf-core umfasst viele Features und Konventionen, die die Community über Jahre entwickelt hat.
Stattdessen konzentrieren wir uns auf die wesentlichen Konzepte, die dir helfen werden, loszulegen und zu verstehen, wie nf-core funktioniert.

### Lektionsplan

Wir haben dies in fünf Teile unterteilt, die sich jeweils auf bestimmte Aspekte der Nutzung von nf-core Ressourcen konzentrieren.

| Kurskapitel                                                       | Zusammenfassung                                                                                                                                                                   | Geschätzte Dauer |
| ----------------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ---------------- |
| [Teil 1: Eine Demo-Pipeline ausführen](./01_run_demo.md)          | Führe eine bestehende nf-core Pipeline aus und untersuche ihre Code-Struktur, um ein Gefühl dafür zu bekommen, was diese Pipelines von einfachen Nextflow Workflows unterscheidet | 30 Min.          |
| [Teil 2: Hello für nf-core umschreiben](./02_rewrite_hello.md)    | Passe einen bestehenden Workflow an das nf-core Template-Gerüst an, ausgehend vom einfachen Workflow aus dem Kurs [Hello Nextflow](../hello_nextflow/index.md)                    | 60 Min.          |
| [Teil 3: Ein nf-core Modul verwenden](./03_use_module.md)         | Erkunde die Community-Modulbibliothek und lerne, vorgefertigte, getestete Module zu integrieren, die gängige Bioinformatik-Tools wrappen                                          | 30 Min.          |
| [Teil 4: Ein nf-core Modul erstellen](./04_make_module.md)        | Erstelle dein eigenes nf-core-style Modul unter Verwendung der spezifischen Struktur, Namenskonventionen und Metadaten-Anforderungen, die von nf-core etabliert wurden            | 30 Min.          |
| [Teil 5: Eingabevalidierung hinzufügen](./05_input_validation.md) | Implementiere Eingabevalidierung sowohl für Kommandozeilen-Parameter als auch für Eingabedateien mit nf-schema                                                                    | 30 Min.          |

Am Ende dieses Kurses kannst du den enormen Reichtum an Ressourcen nutzen, die das nf-core Projekt bietet.

Bereit, den Kurs zu beginnen?

[Mit dem Lernen beginnen :material-arrow-right:](00_orientation.md){ .md-button .md-button--primary }
