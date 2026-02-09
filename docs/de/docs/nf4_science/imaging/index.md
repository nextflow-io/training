---
title: Nextflow run für Bildgebung
hide:
  - toc
---

# Nextflow run für Bildgebung

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } KI-gestützte Übersetzung - [mehr erfahren & Verbesserungen vorschlagen](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Dieser Trainingskurs richtet sich an Forschende im Bereich Bildgebung und räumliche Biologie, die daran interessiert sind, Datenanalyse-Pipelines auszuführen und anzupassen.
Er vermittelt grundlegende Nextflow-Konzepte für das Ausführen, Organisieren und Konfigurieren von Workflows anhand von [nf-core/molkart](https://nf-co.re/molkart), einer Pipeline zur Verarbeitung von Molecular Cartography Spatial-Transkriptomik-Daten.
Die hier erlernten Fähigkeiten sind auf jede Nextflow- oder nf-core-Pipeline übertragbar.

Los geht's! Klicke unten auf die Schaltfläche "Open in GitHub Codespaces", um die Trainingsumgebung zu starten (vorzugsweise in einem separaten Tab), und lies dann weiter, während sie lädt.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

## Lernziele

Durch die Bearbeitung dieses Kurses lernst du, wie du grundlegende Nextflow-Konzepte und -Tools auf die Ausführung von Bildgebungsanalyse-Pipelines anwendest.

Am Ende dieses Workshops kannst du:

- Einen Nextflow-Workflow lokal zu starten und die Ausführung zu überwachen
- Von Nextflow generierte Ausgaben (Ergebnisse) und Log-Dateien zu finden und zu interpretieren
- Eine nf-core-Pipeline mit Testdaten und benutzerdefinierten Eingaben auszuführen
- Die Pipeline-Ausführung mithilfe von Profilen und Parameterdateien zu konfigurieren
- Eingaben mithilfe von Samplesheets und Befehlszeilenparametern zu verwalten

## Zielgruppe & Voraussetzungen

Dieser Kurs setzt ein Mindestmaß an Vertrautheit mit Folgendem voraus:

- Erfahrung mit der Befehlszeile
- Grundlegende Vertrautheit mit Bildgebungsdateiformaten (TIFF-Bilder, tabellarische Daten)

Für technische Anforderungen und die Einrichtung der Umgebung siehe den Mini-Kurs [Umgebungseinrichtung](../../envsetup/).
