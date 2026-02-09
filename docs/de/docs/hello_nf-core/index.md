---
title: Hello nf-core
hide:
  - toc
page_type: index_page
index_type: course
additional_information:
  technical_requirements: true
  learning_objectives:
    - nf-core-Pipelines abrufen, starten und deren Ausführung verwalten
    - Die Code-Struktur und Projektorganisation von nf-core-Pipelines beschreiben
    - Eine einfache nf-core-kompatible Pipeline aus einem Template erstellen
    - Einen einfachen Nextflow-Workflow auf nf-core-Standards upgraden
    - nf-core-Module zu einer nf-core-kompatiblen Pipeline hinzufügen
    - Eigene Module zu nf-core beitragen
    - Eingaben und Parameter mit nf-core-Tools validieren
  audience_prerequisites:
    - "**Zielgruppe:** Dieser Kurs richtet sich an Lernende, die bereits mit den Grundlagen von Nextflow vertraut sind und lernen möchten, nf-core-Ressourcen und Best Practices zu nutzen."
    - "**Fähigkeiten:** Vertrautheit mit der Kommandozeile, grundlegende Scripting-Konzepte und gängige Dateiformate werden vorausgesetzt."
    - "**Kurse:** Du musst den Kurs [Hello Nextflow](../hello_nextflow/index.md) oder einen gleichwertigen Kurs abgeschlossen haben."
    - "**Fachbereich:** Alle Übungen sind fachbereichsunabhängig, daher sind keine wissenschaftlichen Vorkenntnisse erforderlich."
---

# Hello nf-core

**Hello nf-core ist eine praktische Einführung in die Nutzung von nf-core-Ressourcen und Best Practices.**

![nf-core logo](./img/nf-core-logo.png#only-light)
![nf-core logo](./img/nf-core-logo-darkbg.png#only-dark)

Anhand praktischer Beispiele und angeleiteter Übungen lernst du, nf-core-kompatible Module und Pipelines zu nutzen und zu entwickeln sowie nf-core-Tools effektiv einzusetzen.

Du wirst die Fähigkeiten und das Selbstvertrauen erlangen, um Pipelines nach nf-core Best Practices zu entwickeln.

<!-- additional_information -->

## Kursübersicht

Dieser Kurs ist praxisorientiert gestaltet, mit zielgerichteten Übungen, die Informationen schrittweise einführen.

Du wirst in [**nf-core**](https://nf-co.re/) eingeführt, eine Community-Initiative zur Entwicklung und Pflege einer kuratierten Sammlung wissenschaftlicher Pipelines, die mit Nextflow erstellt wurden, sowie relevanter Tools und Richtlinien, die offene Entwicklung, Testing und Peer Review fördern ([Nat Biotechnol 38, 276–278 (2020)](https://www.nature.com/articles/s41587-020-0439-x), [Genome Biol 26, 228 (2025)](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-025-03673-9)).

Die von der nf-core-Community entwickelten Pipelines sind modular, skalierbar und portabel gestaltet, sodass Forscher\*innen sie einfach anpassen und mit ihren eigenen Daten und Rechenressourcen ausführen können.
Die vom Projekt durchgesetzten Best-Practice-Richtlinien stellen außerdem sicher, dass die Pipelines robust, gut dokumentiert und gegen reale Datensätze validiert sind.
Dies trägt dazu bei, die Zuverlässigkeit und Reproduzierbarkeit wissenschaftlicher Analysen zu erhöhen und ermöglicht es Forscher\*innen letztendlich, ihre wissenschaftlichen Entdeckungen zu beschleunigen.

Wir werden in diesem Kurs nicht alles abdecken, was es über nf-core-Pipelines zu wissen gibt, denn nf-core umfasst viele Features und Konventionen, die von der Community über Jahre entwickelt wurden.
Stattdessen konzentrieren wir uns auf die wesentlichen Konzepte, die dir den Einstieg erleichtern und dir helfen zu verstehen, wie nf-core funktioniert.

### Lektionsplan

Wir haben dies in fünf Teile unterteilt, die sich jeweils auf spezifische Aspekte der Nutzung von nf-core-Ressourcen konzentrieren.

| Kurskapitel                                                | Zusammenfassung                                                                                                                                                                                  | Geschätzte Dauer |
| ---------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ | ---------------- |
| [Teil 1: Eine Demo-Pipeline ausführen](./01_run_demo.md)   | Führe eine bestehende nf-core-Pipeline aus und untersuche ihre Code-Struktur, um ein Gefühl dafür zu bekommen, was diese Pipelines von einfachen Nextflow-Workflows unterscheidet               | 30 Min.          |
| [Teil 2: Hello für nf-core umschreiben](./02_rewrite_hello.md) | Passe einen bestehenden Workflow an das nf-core-Template-Gerüst an, ausgehend vom einfachen Workflow aus dem Kurs [Hello Nextflow](../hello_nextflow/index.md)                                  | 60 Min.          |
| [Teil 3: Ein nf-core-Modul verwenden](./03_use_module.md)  | Erkunde die Community-Modul-Bibliothek und lerne, vorgefertigte, getestete Module zu integrieren, die gängige Bioinformatik-Tools umschließen                                                    | 30 Min.          |
| [Teil 4: Ein nf-core-Modul erstellen](./04_make_module.md) | Erstelle dein eigenes nf-core-Modul unter Verwendung der spezifischen Struktur, Namenskonventionen und Metadaten-Anforderungen, die von nf-core festgelegt wurden                               | 30 Min.          |
| [Teil 5: Eingabevalidierung hinzufügen](./05_input_validation.md) | Implementiere Eingabevalidierung sowohl für Kommandozeilen-Parameter als auch für Eingabedateien mit nf-schema                                                                                  | 30 Min.          |

Am Ende dieses Kurses wirst du in der Lage sein, den enormen Reichtum an Ressourcen zu nutzen, die das nf-core-Projekt bietet.

Bereit, den Kurs zu beginnen?

[Jetzt loslegen :material-arrow-right:](00_orientation.md){ .md-button .md-button--primary }
