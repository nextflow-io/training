---
title: Home
description: Willkommen im Nextflow Community Training Portal!
hide:
  - toc
  - footer
---

# Nextflow Training

<div class="grid cards" markdown>

-   :material-book-open-variant:{ .lg .middle } __Self-service-Kurse__

    ---

    **Willkommen im Nextflow Community Training Portal!**

    Die unten aufgeführten Trainingskurse sind als Self-service-Ressource konzipiert.
    Du kannst sie jederzeit in deinem eigenen Tempo durcharbeiten – entweder in der webbasierten Umgebung, die wir über GitHub Codespaces bereitstellen, oder in deiner eigenen Umgebung.

    [Kurse entdecken :material-arrow-right:](#catalog-of-nextflow-training-courses){ .md-button .md-button--primary .mt-1 }

-   :material-information-outline:{ .lg .middle } __Weitere Informationen__

    ---

    ??? warning "Versionskompatibilität"

        <!-- Any update to this content needs to be copied to the local installation page -->
        **Ab Januar 2026 erfordern alle unsere Nextflow-Trainingskurse Nextflow Version 25.10.2 oder höher mit aktivierter strikter Syntax, sofern nicht anders angegeben.**

        Weitere Informationen zu Versionsanforderungen und strikter Syntax findest du im [Nextflow-Docs-Migrationsleitfaden](https://nextflow.io/docs/latest/strict-syntax.html).

        Ältere Versionen des Trainingsmaterials, die der vorherigen Syntax entsprechen, sind über die Versionsauswahl in der Menüleiste dieser Webseite verfügbar.

    ??? terminal "Umgebungsoptionen"

        Wir stellen eine webbasierte Trainingsumgebung bereit, in der alles vorinstalliert ist, was du für das Training benötigst. Sie ist über GitHub Codespaces verfügbar (erfordert ein kostenloses GitHub-Konto).

        [![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

        Falls das nicht deinen Anforderungen entspricht, sieh dir die anderen [Umgebungsoptionen](./envsetup/index.md) an.

    ??? learning "Trainingsveranstaltungen"

        Wenn du Nextflow-Training lieber im Rahmen einer strukturierten Veranstaltung absolvieren möchtest, gibt es viele Möglichkeiten dazu. Wir empfehlen folgende Optionen:

        - **[Training Weeks]()** – vierteljährlich vom Community-Team organisiert
        - **[Seqera Events](https://seqera.io/events/)** – Präsenz-Trainingsveranstaltungen von Seqera (suche nach „Seqera Sessions" und „Nextflow Summit")
        - **[Nextflow Ambassadors]()** – organisieren Veranstaltungen für ihre lokale Community
        - **[nf-core events](https://nf-co.re/events)** – Community-Hackathons

    ??? people "Informationen für Trainer\*innen"

        Wenn du als Trainer\*in eigene Schulungen durchführst, kannst du unser Material direkt vom Training Portal verwenden, solange du die entsprechenden Quellen angibst. Weitere Details findest du unter „Credits und Beiträge" weiter unten.

        Außerdem würden wir gerne von dir hören, wie wir deine Trainingsarbeit besser unterstützen können! Kontaktiere uns unter [community@seqera.io](mailto:community@seqera.io) oder im Community-Forum (siehe Seite [Hilfe](help.md)).

    ??? licensing "Open-Source-Lizenz und Beitragsrichtlinie"

        [![Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International (CC BY-NC-SA 4.0)](assets/img/cc_by-nc-sa.svg){ align=right }](https://creativecommons.org/licenses/by-nc-sa/4.0/)

        Dieses Trainingsmaterial wird von [Seqera](https://seqera.io) entwickelt und gepflegt und unter einer Open-Source-Lizenz ([CC BY-NC-SA](https://creativecommons.org/licenses/by-nc-sa/4.0/)) zum Wohl der Community veröffentlicht. Wenn du dieses Material auf eine Weise nutzen möchtest, die außerhalb des Lizenzrahmens liegt (beachte die Einschränkungen zur kommerziellen Nutzung und Weiterverbreitung), kontaktiere uns bitte unter [community@seqera.io](mailto:community@seqera.io), um dein Anliegen zu besprechen.

        Wir freuen uns über Verbesserungen, Korrekturen und Fehlerberichte aus der Community. Jede Seite hat ein :material-file-edit-outline:-Symbol oben rechts, das zum Code-Repository verlinkt. Dort kannst du Probleme melden oder Änderungen am Trainingsmaterial per Pull Request vorschlagen. Weitere Details findest du in der `README.md` im Repository.

</div>

!!! note "KI-gestützte Übersetzung"

    Diese Übersetzung wurde mit künstlicher Intelligenz erstellt und von menschlichen Übersetzern überprüft.
    Wir freuen uns über Feedback und Verbesserungsvorschläge.
    Weitere Informationen findest du in unserer [Übersetzungsanleitung](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md).

## Catalog of Nextflow training courses

<div class="grid cards" markdown>

-   :material-walk:{ .lg .middle } __Einsteiger-Track__

    ---

    ### :material-compass:{.nextflow-primary} Nextflow für Einsteiger\*innen {.mt-1}

    Domänenunabhängige Kurse für alle, die noch keine Erfahrung mit Nextflow haben. Jeder Kurs besteht aus einer Reihe von Trainingsmodulen, die darauf ausgelegt sind, Lernenden schrittweise neue Fähigkeiten zu vermitteln.

    ??? courses "**Hello Nextflow:** Lerne, eigene Pipelines zu entwickeln"

        Dieser Kurs behandelt die Kernkomponenten der Nextflow-Sprache in ausreichendem Detail, um einfache, aber voll funktionsfähige Pipelines zu entwickeln – plus wichtige Aspekte des Pipeline-Designs, der Entwicklung und der Konfiguration.

        [Hello Nextflow Training starten :material-arrow-right:](hello_nextflow/index.md){ .md-button .md-button--secondary }

    ??? courses "**Nextflow Run:** Lerne, bestehende Pipelines auszuführen"

        Eine kompakte Einführung in das Ausführen und Konfigurieren von Nextflow-Pipelines, basierend auf dem Hello Nextflow-Entwicklerkurs, aber mit weniger Fokus auf Code. Behandelt werden Ausführung, Ausgaben, grundlegende Code-Struktur und Konfiguration für verschiedene Rechenumgebungen.

        [Nextflow Run Training starten :material-arrow-right:](nextflow_run/index.md){ .md-button .md-button--secondary }

    ---

    ### :material-microscope:{.nextflow-primary} Nextflow für die Wissenschaft {.mt-1}

    Lerne, die in „Hello Nextflow" vorgestellten Konzepte und Komponenten auf konkrete wissenschaftliche Anwendungsfälle anzuwenden.

    ??? courses "**Nextflow for Genomics** (Variantenanalyse)"

        Für Forscher\*innen, die lernen möchten, eigene Genomik-Pipelines zu entwickeln. Der Kurs verwendet einen Variantenanalyse-Anwendungsfall, um zu zeigen, wie eine einfache, aber funktionsfähige Genomik-Pipeline entwickelt wird.

        [Nextflow for Genomics Training starten :material-arrow-right:](nf4_science/genomics/index.md){ .md-button .md-button--secondary }

    ??? courses "**Nextflow for RNAseq** (Bulk RNAseq)"

        Für Forscher\*innen, die lernen möchten, eigene RNAseq-Pipelines zu entwickeln. Der Kurs verwendet einen Bulk-RNAseq-Verarbeitungsanwendungsfall, um zu zeigen, wie eine einfache, aber funktionsfähige RNAseq-Pipeline entwickelt wird.

        [Nextflow for RNAseq Training starten :material-arrow-right:](nf4_science/rnaseq/index.md){ .md-button .md-button--secondary }

    ??? courses "**Nextflow for Imaging** (Spatial Omics)"

        Für Forscher\*innen in den Bereichen Bildgebung und Spatial Omics, die lernen möchten, Analyse-Pipelines auszuführen und anzupassen. Der Kurs verwendet die nf-core/molkart-Pipeline, um zu zeigen, wie Nextflow-Pipelines ausgeführt, konfiguriert und mit Eingabedaten versorgt werden.

        [Nextflow for Imaging Training starten :material-arrow-right:](nf4_science/imaging/index.md){ .md-button .md-button--secondary }

-   :material-run:{ .lg .middle } __Fortgeschrittenen-Track__

    ---

    ### :material-bridge:{.nextflow-primary} Von Nextflow zu nf-core {.mt-1}

    Lerne, Code und Best Practices aus dem [nf-core](https://nf-co.re/)-Community-Projekt zu nutzen.

    Diese Kurse führen dich von den Nextflow-Grundlagen zu den nf-core-Best-Practices.
    Verstehe, wie und warum die nf-core-Community Pipelines entwickelt, und lerne, wie du diese Techniken einsetzen und dazu beitragen kannst.

    ??? courses "**Hello nf-core:** Erste Schritte mit nf-core"

        Für Entwickler\*innen, die lernen möchten, [nf-core](https://nf-co.re/)-konforme Pipelines auszuführen und zu entwickeln. Der Kurs behandelt die Struktur von nf-core-Pipelines in ausreichendem Detail, um einfache, aber voll funktionsfähige Pipelines zu entwickeln, die dem nf-core-Template und den Entwicklungs-Best-Practices folgen, sowie bestehende nf-core-Module zu verwenden.

        [Hello nf-core Training starten :material-arrow-right:](hello_nf-core/index.md){ .md-button .md-button--secondary }

    ---

    ### :material-rocket-launch:{.nextflow-primary} Fortgeschrittenes Nextflow Training {.mt-1}

    Lerne fortgeschrittene Konzepte und Mechanismen für die Entwicklung und den Einsatz von Nextflow-Pipelines für reale Anwendungsfälle.

    ??? courses "**Side Quests:** Tiefe Einblicke in einzelne Themen"

        Eigenständige Mini-Kurse für Nextflow-Entwickler\*innen, die ihr Repertoire erweitern und/oder ihre Kenntnisse zu bestimmten Themen vertiefen möchten. Sie sind linear aufgebaut, können aber in beliebiger Reihenfolge absolviert werden (siehe Abhängigkeiten in der jeweiligen Mini-Kurs-Übersicht).

        [Side Quests durchsuchen :material-arrow-right:](side_quests/index.md){ .md-button .md-button--secondary }

    ??? courses "**Training Collections:** Empfohlene Lernpfade durch die Side Quests"

        Training Collections kombinieren mehrere Side Quests, um ein umfassendes Lernerlebnis rund um ein bestimmtes Thema oder einen Anwendungsfall zu bieten.

        [Training Collections durchsuchen :material-arrow-right:](training_collections/index.md){ .md-button .md-button--secondary }

</div>

!!! info "Suchst du nach archivierten Trainingsmaterialien?"

    Ältere Trainingsmaterialien (Fundamentals Training, Advanced Training und andere experimentelle Kurse) wurden aus dem Training Portal entfernt, da sie mit der strikten Syntax von Nextflow 3.0 nicht kompatibel sind.
    Falls du Zugang zu diesen Materialien benötigst, sind sie in der [Git-Historie](https://github.com/nextflow-io/training) vor Januar 2026 verfügbar.

---

<div markdown class="homepage_logos">

![Seqera](assets/img/seqera_logo.png#only-light)

![Seqera](assets/img/seqera_logo_dark.png#only-dark)

</div>
