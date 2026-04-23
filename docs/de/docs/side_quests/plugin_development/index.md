---
title: Plugin-Entwicklung
hide:
  - toc
---

# Plugin-Entwicklung

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } KI-gestützte Übersetzung - [mehr erfahren & Verbesserungen vorschlagen](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Das Plugin-System von Nextflow ermöglicht es dir, die Sprache um eigene Funktionen, Monitoring-Hooks, Ausführungs-Backends und mehr zu erweitern.
Plugins ermöglichen es der Community, Nextflow um neue Features zu ergänzen, ohne den Kern zu verändern – ideal, um wiederverwendbare Funktionalität über Pipelines hinweg zu teilen.

In diesem Training lernst du, wie du vorhandene Plugins verwendest und optional eigene erstellst.

## Zielgruppe & Voraussetzungen

Teil 1 behandelt die Verwendung vorhandener Plugins und ist für alle Nextflow-Nutzer\*innen relevant.
Die Teile 2–6 befassen sich mit der Entwicklung eigener Plugins und beinhalten Groovy-Code sowie Build-Tools.
Vorkenntnisse in Java oder Groovy sind nicht erforderlich.

**Voraussetzungen**

- Ein GitHub-Konto ODER eine lokale Installation, wie [hier](../../envsetup/02_local) beschrieben.
- Abgeschlossener [Hello Nextflow](../../hello_nextflow/index.md)-Kurs oder gleichwertige Kenntnisse.
- Java 21 oder höher (in der Trainingsumgebung enthalten; nur für die Teile 2–6 benötigt).

**Arbeitsverzeichnis:** `side-quests/plugin_development`

## Lernziele

Nach Abschluss dieses Trainings kannst du:

**Plugins verwenden (Teil 1):**

- Vorhandene Plugins in deinen Workflows installieren und konfigurieren
- Plugin-Funktionen importieren und verwenden

**Plugins entwickeln (Teile 2–6):**

- Ein neues Plugin-Projekt mit dem integrierten Projektgenerator von Nextflow erstellen
- Eigene Funktionen implementieren, die aus Workflows aufgerufen werden können
- Dein Plugin lokal bauen, testen und installieren
- Workflow-Ereignisse überwachen (z. B. Aufgabenabschluss, Pipeline-Start/-Ende) für eigenes Logging oder Benachrichtigungen
- Konfigurationsoptionen hinzufügen, um Plugins anpassbar zu machen
- Dein Plugin verteilen

## Kursplan

#### Teil 1: Plugin-Grundlagen

Verwende vorhandene Plugins in einem Nextflow-Workflow und konfiguriere ihr Verhalten.

#### Teil 2: Ein Plugin-Projekt erstellen

Generiere ein neues Plugin-Projekt und untersuche seine Struktur.

#### Teil 3: Eigene Funktionen

Implementiere eigene Funktionen, baue dein Plugin und führe es in einem Workflow aus.

#### Teil 4: Testen

Schreibe und führe Unit-Tests mit dem Spock-Framework aus.

#### Teil 5: Workflow-Monitoring

Reagiere auf Ereignisse wie den Abschluss von Aufgaben, um einen Aufgabenzähler zu erstellen.

#### Teil 6: Konfiguration & Verteilung

Lese Einstellungen aus `nextflow.config`, um dein Plugin anpassbar zu machen, und lerne, wie du es teilst.

Bereit, den Kurs zu starten?

[Jetzt loslegen :material-arrow-right:](01_plugin_basics.md){ .md-button .md-button--primary }
