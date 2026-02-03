---
title: Startseite
description: Willkommen im Nextflow Community-Trainingsportal!
hide:
  - toc
  - footer
---

# Nextflow Training

<div class="grid cards" markdown>

-   :material-book-open-variant:{ .lg .middle } __Selbstlernkurse__

    ---

    **Willkommen im Nextflow Community-Trainingsportal!**

    Die unten aufgeführten Trainingskurse sind als Selbstlernressource konzipiert.
    Du kannst sie jederzeit eigenständig durcharbeiten, entweder in der webbasierten Umgebung, die wir über Github Codespaces bereitstellen, oder in deiner eigenen Umgebung.

    [Kurse erkunden :material-arrow-right:](#katalog-der-nextflow-trainingskurse){ .md-button .md-button--primary .mt-1 }

-   :material-information-outline:{ .lg .middle } __Zusätzliche Informationen__

    ---

    ??? warning "Versionskompatibilität"

        <!-- Any update to this content needs to be copied to the local installation page -->
        **Ab Januar 2026 erfordern alle unsere Nextflow-Trainingskurse Nextflow Version 25.10.2 oder höher mit aktivierter strikter Syntax, sofern nicht anders angegeben.**

        Weitere Informationen zu Versionsanforderungen und strikter Syntax findest du im [Nextflow-Dokumentations-Migrationshandbuch](https://nextflow.io/docs/latest/strict-syntax.html).

        Ältere Versionen des Trainingsmaterials, die früherer Syntax entsprechen, sind über den Versionsauswähler in der Menüleiste dieser Webseite verfügbar.

    ??? terminal "Umgebungsoptionen"

        Wir bieten eine webbasierte Trainingsumgebung an, in der alles vorinstalliert ist, was du für das Training benötigst. Diese ist über Github Codespaces verfügbar (erfordert ein kostenloses GitHub-Konto).

        [![In GitHub Codespaces öffnen](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

        Wenn dies nicht deinen Anforderungen entspricht, siehe die anderen [Umgebungsoptionen](./envsetup/index.md).

    ??? learning "Trainingsveranstaltungen"

        Wenn du lieber an einem Nextflow-Training im Rahmen einer strukturierten Veranstaltung teilnehmen möchtest, gibt es viele Möglichkeiten dazu. Wir empfehlen, die folgenden Optionen zu prüfen:

        - **[Training Weeks]()** werden vierteljährlich vom Community-Team organisiert
        - **[Seqera Events](https://seqera.io/events/)** beinhalten Präsenztrainings, die von Seqera organisiert werden (suche nach 'Seqera Sessions' und 'Nextflow Summit')
        - **[Nextflow Ambassadors]()** organisieren Veranstaltungen für ihre lokale Community
        - **[nf-core events](https://nf-co.re/events)** beinhalten Community-Hackathons

    ??? people "Informationen für Trainer"

        Wenn du als Dozent eigene Trainings durchführst, kannst du unsere Materialien gerne direkt aus dem Trainingsportal verwenden, solange du die entsprechende Quellenangabe machst. Siehe 'Credits und Beiträge' unten für Details.

        Außerdem würden wir gerne von dir erfahren, wie wir deine Trainingsarbeit besser unterstützen können! Bitte kontaktiere uns unter [community@seqera.io](mailto:community@seqera.io) oder im Community-Forum (siehe Seite [Hilfe](help.md)).

    ??? licensing "Open-Source-Lizenz und Beitragsrichtlinien"

        [![Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International (CC BY-NC-SA 4.0)](assets/img/cc_by-nc-sa.svg){ align=right }](https://creativecommons.org/licenses/by-nc-sa/4.0/)

        Dieses Trainingsmaterial wird von [Seqera](https://seqera.io) entwickelt und gepflegt und unter einer Open-Source-Lizenz ([CC BY-NC-SA](https://creativecommons.org/licenses/by-nc-sa/4.0/)) zum Nutzen der Community veröffentlicht. Wenn du dieses Material auf eine Weise nutzen möchtest, die außerhalb des Geltungsbereichs der Lizenz liegt (beachte die Einschränkungen bei kommerzieller Nutzung und Weiterverbreitung), kontaktiere uns bitte unter [community@seqera.io](mailto:community@seqera.io), um deine Anfrage zu besprechen.

        Wir freuen uns über Verbesserungen, Korrekturen und Fehlermeldungen aus der Community. Jede Seite hat ein :material-file-edit-outline: Symbol oben rechts, das zum Code-Repository verlinkt, wo du Probleme melden oder Änderungen am Trainingsmaterial per Pull Request vorschlagen kannst. Siehe die \`README.md\` im Repository für weitere Details.

</div>

!!! note "KI-gestützte Übersetzung"

    Diese Übersetzung wurde mit künstlicher Intelligenz erstellt und von menschlichen Übersetzern überprüft.
    Wir freuen uns über Feedback und Verbesserungsvorschläge.
    Weitere Informationen findest du in unserer [Übersetzungsanleitung](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md).

## Katalog der Nextflow-Trainingskurse

<div class="grid cards" markdown>

-   :material-walk:{ .lg .middle } __Einführungskurse__

    ---

    ### :material-compass:{.nextflow-primary} Nextflow für Einsteiger {.mt-1}

    Fachbereichsunabhängige Kurse für alle, die komplett neu bei Nextflow sind. Jeder Kurs besteht aus mehreren Modulen, die Fähigkeiten schrittweise aufbauen.

    ??? courses "**Hello Nextflow:** Lerne, eigene Pipelines zu entwickeln"

        Dieser Kurs behandelt die Kernkomponenten der Nextflow-Sprache ausführlich genug, um einfache, aber voll funktionsfähige Pipelines zu entwickeln, plus Schlüsselelemente des Pipeline-Designs, der Entwicklung und der Konfigurationspraktiken.

        [Hello Nextflow Training starten :material-arrow-right:](hello_nextflow/index.md){ .md-button .md-button--secondary }

    ??? courses "**Nextflow Run:** Lerne, bestehende Pipelines auszuführen"

        Eine kompakte Einführung in die Ausführung und Konfiguration von Nextflow-Pipelines, basierend auf dem Hello Nextflow-Kurs für Entwickler*innen, aber mit weniger Fokus auf Code. Behandelt Ausführung, Ausgaben, grundlegende Codestruktur und Konfiguration für verschiedene Rechenumgebungen.

        [Nextflow Run Training starten :material-arrow-right:](nextflow_run/index.md){ .md-button .md-button--secondary }

    ---

    ### :material-microscope:{.nextflow-primary} Nextflow für die Wissenschaft {.mt-1}

    Lerne, die in 'Hello Nextflow' vorgestellten Konzepte und Komponenten auf spezifische wissenschaftliche Anwendungsfälle anzuwenden.

    ??? courses "**Nextflow für Genomik** (Variantenanalyse)"

        Für Forschende, die lernen möchten, eigene Genomik-Pipelines zu entwickeln. Der Kurs verwendet einen Anwendungsfall zur Variantenanalyse, um zu demonstrieren, wie man eine einfache, aber funktionale Genomik-Pipeline entwickelt.

        [Nextflow für Genomik Training starten :material-arrow-right:](nf4_science/genomics/){ .md-button .md-button--secondary }

    ??? courses "**Nextflow für RNAseq** (Bulk-RNAseq)"

        Für Forschende, die lernen möchten, eigene RNAseq-Pipelines zu entwickeln. Der Kurs verwendet einen Bulk-RNAseq-Verarbeitungs-Anwendungsfall, um zu demonstrieren, wie man eine einfache, aber funktionale RNAseq-Pipeline entwickelt.

        [Nextflow für RNAseq Training starten :material-arrow-right:](nf4_science/rnaseq/){ .md-button .md-button--secondary }

    ??? courses "**Nextflow für Bildgebung** (Spatial Omics)"

        Für Forschende in der Bildgebung und Spatial Omics, die lernen möchten, Analyse-Pipelines auszuführen und anzupassen. Der Kurs verwendet die nf-core/molkart-Pipeline, um eine biologisch relevante Pipeline bereitzustellen und zu demonstrieren, wie man Nextflow-Pipelines und Workflows ausführt, konfiguriert und Eingaben verwaltet.

        [Nextflow für Bildgebung Training starten :material-arrow-right:](nf4_science/imaging/){ .md-button .md-button--secondary }

-   :material-run:{ .lg .middle } __Fortgeschrittenenkurse__

    ---

    ### :material-bridge:{.nextflow-primary} Von Nextflow zu nf-core {.mt-1}

    Lerne, Code und Best Practices aus dem [nf-core](https://nf-co.re/) Community-Projekt zu nutzen.

    Diese Kurse helfen dir, von den Nextflow-Grundlagen zu den nf-core Best Practices zu gelangen.
    Verstehe, wie und warum die nf-core-Community Pipelines entwickelt, und wie du beitragen und diese Techniken wiederverwenden kannst.

    ??? courses "**Hello nf-core:** Erste Schritte mit nf-core"

        Für Entwickler*innen, die lernen möchten, [nf-core](https://nf-co.re/)-konforme Pipelines auszuführen und zu entwickeln. Der Kurs behandelt die Struktur von nf-core-Pipelines ausführlich genug, um einfache, aber voll funktionsfähige Pipelines zu entwickeln, die dem nf-core-Template und den Entwicklungs-Best-Practices folgen, sowie bestehende nf-core-Module zu verwenden.

        [Hello nf-core Training starten :material-arrow-right:](hello_nf-core/index.md){ .md-button .md-button--secondary }

    ---

    ### :material-rocket-launch:{.nextflow-primary} Fortgeschrittenes Nextflow-Training {.mt-1}

    Lerne fortgeschrittene Konzepte und Mechanismen für die Entwicklung und Bereitstellung von Nextflow-Pipelines zur Lösung realer Anwendungsfälle.

    ??? courses "**Side Quests:** Vertiefungen zu eigenständigen Themen"

        Eigenständige Minikurse für Nextflow-Entwickler*innen, die ihr Repertoire erweitern und/oder ihre Fähigkeiten zu bestimmten Themen vertiefen möchten. Sie werden linear präsentiert, können aber in beliebiger Reihenfolge absolviert werden (siehe Abhängigkeiten in jeder Minikurs-Übersicht).

        [Side Quests durchstöbern :material-arrow-right:](side_quests/){ .md-button .md-button--secondary }

    ??? courses "**Training Collections:** Empfohlene Lernpfade durch die Side Quests"

        Training Collections kombinieren mehrere Side Quests, um eine umfassende Lernerfahrung zu einem bestimmten Thema oder Anwendungsfall zu bieten.

        [Training Collections durchstöbern :material-arrow-right:](training_collections/){ .md-button .md-button--secondary }

</div>

!!! info "Auf der Suche nach archivierten Trainingsmaterialien?"

    Ältere Trainingsmaterialien (Fundamentals Training, Advanced Training und andere experimentelle Kurse) wurden aus dem Trainingsportal entfernt, da sie nicht mit der strikten Nextflow 3.0-Syntax kompatibel sind.
    Wenn du Zugriff auf diese Materialien benötigst, sind sie in der [Git-Historie](https://github.com/nextflow-io/training) vor Januar 2026 verfügbar.

---

<div markdown class="homepage_logos">

![Seqera](assets/img/seqera_logo.png#only-light)

![Seqera](assets/img/seqera_logo_dark.png#only-dark)

</div>
