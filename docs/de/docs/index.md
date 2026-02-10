---
title: Startseite
description: Willkommen im Nextflow Community-Trainingsportal!
hide:
  - toc
  - footer
---

# Nextflow Training

<div class="grid cards" markdown>

-   :material-book-open-variant:{ .lg .middle } __Self-Service-Kurse__

    ---

    **Willkommen im Nextflow Community-Trainingsportal!**

    Die unten aufgeführten Trainingskurse sind als Self-Service-Ressource konzipiert.
    Du kannst sie jederzeit selbstständig durcharbeiten – entweder in der webbasierten Umgebung, die wir über Github Codespaces bereitstellen, oder in deiner eigenen Umgebung.

    [Kurse erkunden :material-arrow-right:](#catalog-of-nextflow-training-courses){ .md-button .md-button--primary .mt-1 }

-   :material-information-outline:{ .lg .middle } __Zusätzliche Informationen__

    ---

    ??? warning "Versionskompatibilität"

        <!-- Any update to this content needs to be copied to the local installation page -->
        **Seit Januar 2026 erfordern alle unsere Nextflow-Trainingskurse Nextflow Version 25.10.2 oder höher mit aktivierter strikter Syntax, sofern nicht anders angegeben.**

        Weitere Informationen zu Versionsanforderungen und strikter Syntax findest du im [Nextflow Docs Migration Guide](https://nextflow.io/docs/latest/strict-syntax.html).

        Ältere Versionen des Trainingsmaterials, die der früheren Syntax entsprechen, sind über den Versionsauswähler in der Menüleiste dieser Webseite verfügbar.

    ??? terminal "Umgebungsoptionen"

        Wir stellen eine webbasierte Trainingsumgebung bereit, in der alles vorinstalliert ist, was du für das Training benötigst. Diese ist über Github Codespaces verfügbar (erfordert ein kostenloses GitHub-Konto).

        [![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

        Falls dies nicht deinen Anforderungen entspricht, sieh dir bitte die anderen [Umgebungsoptionen](./envsetup/index.md) an.

    ??? learning "Trainingsveranstaltungen"

        Wenn du Nextflow-Training lieber im Rahmen einer strukturierten Veranstaltung absolvieren möchtest, gibt es viele Möglichkeiten dazu. Wir empfehlen die folgenden Optionen:

        - **[Training Weeks]()**, die vierteljährlich vom Community-Team organisiert werden
        - **[Seqera Events](https://seqera.io/events/)** umfassen Präsenz-Trainingsveranstaltungen, die von Seqera organisiert werden (suche nach 'Seqera Sessions' und 'Nextflow Summit')
        - **[Nextflow Ambassadors]()** organisieren Veranstaltungen für ihre lokale Community
        - **[nf-core events](https://nf-co.re/events)** umfassen Community-Hackathons

    ??? people "Informationen für Trainer\*innen"

        Wenn du als Trainer\*in eigene Schulungen durchführst, kannst du unsere Materialien gerne direkt vom Trainingsportal verwenden, solange du die entsprechende Quellenangabe machst. Details findest du unten unter 'Credits und Beiträge'.

        Außerdem würden wir gerne von dir hören, wie wir deine Trainingsaktivitäten besser unterstützen können! Bitte kontaktiere uns unter [community@seqera.io](mailto:community@seqera.io) oder im Community-Forum (siehe [Hilfe](help.md)-Seite).

    ??? licensing "Open-Source-Lizenz und Beitragsrichtlinien"

        [![Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International (CC BY-NC-SA 4.0)](assets/img/cc_by-nc-sa.svg){ align=right }](https://creativecommons.org/licenses/by-nc-sa/4.0/)

        Dieses Trainingsmaterial wird von [Seqera](https://seqera.io) entwickelt und gepflegt und unter einer Open-Source-Lizenz ([CC BY-NC-SA](https://creativecommons.org/licenses/by-nc-sa/4.0/)) zum Nutzen der Community veröffentlicht. Wenn du dieses Material auf eine Weise verwenden möchtest, die außerhalb des Lizenzumfangs liegt (beachte die Einschränkungen bezüglich kommerzieller Nutzung und Weiterverbreitung), kontaktiere uns bitte unter [community@seqera.io](mailto:community@seqera.io), um deine Anfrage zu besprechen.

        Wir freuen uns über Verbesserungen, Korrekturen und Fehlerberichte aus der Community. Jede Seite hat ein :material-file-edit-outline:-Symbol oben rechts, das zum Code-Repository verlinkt, wo du Probleme melden oder Änderungen am Trainingsquellmaterial über einen Pull Request vorschlagen kannst. Weitere Details findest du in der `README.md` im Repository.

</div>

!!! note "KI-gestützte Übersetzung"

    Diese Übersetzung wurde mit künstlicher Intelligenz erstellt und von menschlichen Übersetzern überprüft.
    Wir freuen uns über Feedback und Verbesserungsvorschläge.
    Weitere Informationen findest du in unserer [Übersetzungsanleitung](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md).

## Katalog der Nextflow-Trainingskurse

<div class="grid cards" markdown>

-   :material-walk:{ .lg .middle } __Einsteiger-Track__

    ---

    ### :material-compass:{.nextflow-primary} Nextflow für Neueinsteiger\*innen {.mt-1}

    Domänenunabhängige Kurse für alle, die völlig neu bei Nextflow sind. Jeder Kurs besteht aus einer Reihe von Trainingsmodulen, die Lernenden helfen sollen, ihre Fähigkeiten schrittweise aufzubauen.

    ??? courses "**Hello Nextflow:** Lerne, eigene Pipelines zu entwickeln"

        Dieser Kurs behandelt die Kernkomponenten der Nextflow-Sprache ausführlich genug, um die Entwicklung einfacher, aber voll funktionsfähiger Pipelines zu ermöglichen, sowie wichtige Elemente des Pipeline-Designs, der Entwicklung und der Konfigurationspraktiken.

        [Hello Nextflow Training starten :material-arrow-right:](hello_nextflow/index.md){ .md-button .md-button--secondary }

    ??? courses "**Nextflow Run:** Lerne, bestehende Pipelines auszuführen"

        Eine kompakte Einführung in das Ausführen und Konfigurieren von Nextflow-Pipelines, basierend auf dem Hello Nextflow Entwickler\*innen-Kurs, aber mit weniger Fokus auf Code. Behandelt Ausführung, Ausgaben, grundlegende Code-Struktur und Konfiguration für verschiedene Rechenumgebungen.

        [Nextflow Run Training starten :material-arrow-right:](nextflow_run/index.md){ .md-button .md-button--secondary }

    ---

    ### :material-microscope:{.nextflow-primary} Nextflow für die Wissenschaft {.mt-1}

    Lerne, die in 'Hello Nextflow' vorgestellten Konzepte und Komponenten auf spezifische wissenschaftliche Anwendungsfälle anzuwenden.

    ??? courses "**Nextflow für Genomik** (Variant Calling)"

        Für Forscher\*innen, die lernen möchten, wie sie eigene Genomik-Pipelines entwickeln. Der Kurs verwendet einen Variant-Calling-Anwendungsfall, um zu demonstrieren, wie man eine einfache, aber funktionale Genomik-Pipeline entwickelt.

        [Nextflow für Genomik Training starten :material-arrow-right:](nf4_science/genomics/){ .md-button .md-button--secondary }

    ??? courses "**Nextflow für RNAseq** (Bulk RNAseq)"

        Für Forscher\*innen, die lernen möchten, wie sie eigene RNAseq-Pipelines entwickeln. Der Kurs verwendet einen Bulk-RNAseq-Verarbeitungsfall, um zu demonstrieren, wie man eine einfache, aber funktionale RNAseq-Pipeline entwickelt.

        [Nextflow für RNAseq Training starten :material-arrow-right:](nf4_science/rnaseq/){ .md-button .md-button--secondary }

    ??? courses "**Nextflow für Bildgebung** (Spatial Omics)"

        Für Forscher\*innen in der Bildgebung und Spatial Omics, die lernen möchten, wie sie Analyse-Pipelines ausführen und anpassen. Der Kurs verwendet die nf-core/molkart-Pipeline, um anhand einer biologisch relevanten Pipeline zu demonstrieren, wie man Nextflow-Pipelines ausführt, konfiguriert und Eingaben verwaltet.

        [Nextflow für Bildgebung Training starten :material-arrow-right:](nf4_science/imaging/){ .md-button .md-button--secondary }

-   :material-run:{ .lg .middle } __Fortgeschrittenen-Track__

    ---

    ### :material-bridge:{.nextflow-primary} Von Nextflow zu nf-core {.mt-1}

    Lerne, Code und Best Practices aus dem [nf-core](https://nf-co.re/)-Community-Projekt zu nutzen.

    Diese Kurse helfen dir, von Nextflow-Grundlagen zu nf-core-Best-Practices zu gelangen.
    Verstehe, wie und warum die nf-core-Community Pipelines erstellt und wie du beitragen und diese Techniken wiederverwenden kannst.

    ??? courses "**Hello nf-core:** Erste Schritte mit nf-core"

        Für Entwickler\*innen, die lernen möchten, [nf-core](https://nf-co.re/)-konforme Pipelines auszuführen und zu entwickeln. Der Kurs behandelt die Struktur von nf-core-Pipelines ausführlich genug, um die Entwicklung einfacher, aber voll funktionsfähiger Pipelines zu ermöglichen, die dem nf-core-Template und den Entwicklungs-Best-Practices folgen, sowie die Verwendung bestehender nf-core-Module.

        [Hello nf-core Training starten :material-arrow-right:](hello_nf-core/index.md){ .md-button .md-button--secondary }

    ---

    ### :material-rocket-launch:{.nextflow-primary} Fortgeschrittenes Nextflow Training {.mt-1}

    Lerne fortgeschrittene Konzepte und Mechanismen für die Entwicklung und Bereitstellung von Nextflow-Pipelines zur Bewältigung realer Anwendungsfälle.

    ??? courses "**Side Quests:** Vertiefungen in eigenständige Themen"

        Eigenständige Mini-Kurse für Nextflow-Entwickler\*innen, die ihre Reichweite erweitern und/oder ihre Fähigkeiten zu bestimmten Themen vertiefen möchten. Sie sind linear präsentiert, können aber in beliebiger Reihenfolge absolviert werden (siehe Abhängigkeiten in jeder Mini-Kurs-Übersicht).

        [Side Quests durchstöbern :material-arrow-right:](side_quests/){ .md-button .md-button--secondary }

    ??? courses "**Training Collections:** Empfohlene Lernpfade durch die Side Quests"

        Training Collections kombinieren mehrere Side Quests, um eine umfassende Lernerfahrung zu einem bestimmten Thema oder Anwendungsfall zu bieten.

        [Training Collections durchstöbern :material-arrow-right:](training_collections/){ .md-button .md-button--secondary }

</div>

!!! info "Suchst du nach archivierten Trainingsmaterialien?"

    Ältere Trainingsmaterialien (Fundamentals Training, Advanced Training und andere experimentelle Kurse) wurden aus dem Trainingsportal entfernt, da sie mit der strikten Syntax von Nextflow 3.0 nicht kompatibel sind.
    Falls du Zugriff auf diese Materialien benötigst, sind sie in der [Git-Historie](https://github.com/nextflow-io/training) vor Januar 2026 verfügbar.

---

<div markdown class="homepage_logos">

![Seqera](assets/img/seqera_logo.png#only-light)

![Seqera](assets/img/seqera_logo_dark.png#only-dark)

</div>
