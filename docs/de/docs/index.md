---
title: Home
description: Willkommen im Nextflow Community Training Portal!
hide:
  - toc
  - footer
---

# Nextflow Training

<div class="grid cards" markdown>

-   :material-book-open-variant:{ .lg .middle } __Selbstlernkurse__

    ---

    **Willkommen im Nextflow Community Training Portal!**

    Arbeite die Kurse unten in deinem eigenen Tempo durch – in unserer webbasierten Umgebung oder in deiner eigenen.
    Jeder Kurs ist praxisorientiert, mit zielgerichteten Übungen, die du eigenständig absolvieren kannst.

    [Kurse entdecken :material-arrow-down:](#catalog-of-nextflow-training-courses){ .md-button .md-button--primary .mt-1 }

-   :material-account-group-outline:{ .lg .middle } __Training Events__

    ---

    **Suchst du etwas über Selbstlernkurse hinaus?**

    Finde strukturierte Training Events, Anleitungen für eigene Trainings sowie unsere Open-Source-Lizenz und Beitragsrichtlinien.

    [Training Events ansehen :material-arrow-right:](training_events.md){ .md-button .md-button--secondary .mt-1 }

</div>

!!! note "KI-gestützte Übersetzung"

    Diese Übersetzung wurde mit künstlicher Intelligenz erstellt und von menschlichen Übersetzern überprüft.
    Wir freuen uns über Feedback und Verbesserungsvorschläge.
    Weitere Informationen findest du in unserer [Übersetzungsanleitung](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md).

## Catalog of Nextflow training courses

<div class="grid cards" markdown>

-   :material-account:{ .lg .middle } __Für Anwender\*innen__

    ---

    ### :material-play-circle:{.nextflow-primary} Pipelines ausführen {.mt-1}

    Lerne, bestehende Pipelines auszuführen, ohne Code schreiben zu müssen.

    ??? courses "**Nextflow Run:** Pipelines mit Nextflow ausführen"

        Eine kompakte Einführung in das Ausführen von Nextflow-Pipelines – ohne Programmierkenntnisse. Themen: Pipelines starten, Ergebnisse abrufen, Container verwenden und die Ausführung grundlegend konfigurieren.

        [Zum Training :material-arrow-right:](nextflow_run/index.md){ .md-button .md-button--secondary }

    ??? courses "**Use nf-core:** Community-Pipelines finden und ausführen"

        Eine kompakte Einführung in das Finden, Ausführen und Konfigurieren von Pipelines aus dem nf-core Community-Projekt – beginnend mit einer minimalen Demo-Pipeline bis hin zu einer produktionsreifen Analyse-Pipeline.

        [Zum Training :material-arrow-right:](nfcore_use/index.md){ .md-button .md-button--secondary }

    ??? courses "**Scale with Seqera:** Pipelines skaliert starten und überwachen"

        Eine praxisorientierte Einführung in das Starten und Überwachen von Nextflow-Pipelines mit der Seqera Platform – sowohl über die Weboberfläche als auch über die Kommandozeile.

        [Zum Training :material-arrow-right:](seqera_scale/index.md){ .md-button .md-button--secondary }

    ---

    ### :material-tune:{.nextflow-primary} Ausführung verwalten {.mt-1}

    Lerne, Pipeline-Ausführungen effektiv zu verwalten.

    ??? courses "**Configure Execution:** Ressourcen, Wiederholungen und Ausführungsprofile konfigurieren"

        Eine praxisorientierte Einführung in die Konfiguration von Nextflow-Pipelines: Anpassung an verschiedene Rechenumgebungen, Steuerung von Ressourcenzuweisungen und Wiederholungen sowie Wechsel zwischen vordefinierten Konfigurationsprofilen.

        [Zum Training :material-arrow-right:](config_exec/index.md){ .md-button .md-button--secondary }

    !!! info compact "Weitere Themen folgen"

        Performance-Tuning, HPC/Cloud-Ausführung und mehr sind für diesen Bereich geplant.
        Stimme ab, was als nächstes behandelt werden soll, in unserer [kurzen Umfrage](https://seqera.typeform.com/to/JCs91e8v).

-   :material-code-tags:{ .lg .middle } __Für Entwickler\*innen__

    ---

    ### :material-wrench:{.nextflow-primary} Pipelines schreiben {.mt-1}

    Lerne, eigene Nextflow-Pipelines zu entwickeln.

    ??? courses "**Hello Nextflow:** Eigene Pipelines von Grund auf entwickeln"

        Dieser Kurs behandelt die Kernkomponenten der Nextflow-Sprache in ausreichender Tiefe, um einfache, aber vollständig funktionsfähige Pipelines zu entwickeln – plus wichtige Aspekte des Pipeline-Designs, der Entwicklung und der Konfiguration.

        [Zum Training :material-arrow-right:](hello_nextflow/index.md){ .md-button .md-button--secondary }

    ??? courses "**Build with nf-core:** nf-core-Tools und -Regeln verwenden"

        Für Nextflow-Entwickler\*innen, die lernen möchten, [nf-core](https://nf-co.re/)-konforme Pipelines zu entwickeln.
        Der Kurs behandelt die Struktur von nf-core-Pipelines in ausreichender Tiefe, um einfache, aber vollständig funktionsfähige Pipelines zu entwickeln, die das nf-core-Template und Best Practices nutzen – einschließlich der Verwendung bestehender nf-core-Module.

        [Zum Training :material-arrow-right:](nfcore_build/index.md){ .md-button .md-button--secondary }

    ??? catalog "**Side Quests:** Fortgeschrittene Nextflow-Themen vertiefen"

        Eine Sammlung eigenständiger Mini-Kurse für Nextflow-Entwickler\*innen, die ihr Wissen erweitern und/oder ihre Kenntnisse zu bestimmten Themen vertiefen möchten.
        Sie sind linear aufgebaut, können aber in beliebiger Reihenfolge absolviert werden (siehe Abhängigkeiten in der jeweiligen Mini-Kurs-Übersicht).

        [Side Quests durchsuchen :material-arrow-right:](side_quests/index.md){ .md-button .md-button--secondary }

    ---

    ### :material-microscope:{.nextflow-primary} Nextflow für die Wissenschaft {.mt-1}

    Lerne, Nextflow-Pipelines für spezifische wissenschaftliche Anwendungen zu entwickeln.

    ??? courses "**Genomics:** Eine Variant-Calling-Pipeline entwickeln"

        Ein Kurs für Forscher\*innen, die lernen möchten, eigene Genomik-Pipelines zu entwickeln – anhand eines Variant-Calling-Anwendungsfalls werden wesentliche Nextflow-Entwicklungsmuster demonstriert.

        [Zum Training :material-arrow-right:](nf4_science/genomics/index.md){ .md-button .md-button--secondary }

    ??? courses "**RNAseq:** Eine Bulk-RNAseq-Verarbeitungs-Pipeline entwickeln"

        Ein Kurs für Forscher\*innen, die lernen möchten, eigene RNAseq-Pipelines zu entwickeln – anhand eines Bulk-RNAseq-Verarbeitungsanwendungsfalls werden wesentliche Nextflow-Entwicklungsmuster demonstriert.

        [Zum Training :material-arrow-right:](nf4_science/rnaseq/index.md){ .md-button .md-button--secondary }

    ??? courses "**Bioimaging:** Imaging-Pipelines ausführen und konfigurieren"

        Ein Kurs für Forscher\*innen, die lernen möchten, Bioimaging-Pipelines auszuführen und zu konfigurieren – anhand von nf-core/molkart werden wesentliche Nextflow-Nutzungsmuster demonstriert.

        [Zum Training :material-arrow-right:](nf4_science/imaging/index.md){ .md-button .md-button--secondary }

</div>

## Setup & Hilfe

<div class="grid cards mb-4" markdown>

-   :material-cog-outline:{ .lg .middle } __Training-Umgebung__

    ---

    Optionen zum Einrichten deiner Umgebung für die Nextflow-Trainings.

    [Training-Umgebungen ansehen :material-arrow-right:](envsetup/index.md){ .md-button .md-button--secondary }

-   :material-tag-outline:{ .lg .middle } __Nextflow-Versionen__

    ---

    Die Entwicklung der Nextflow-Syntaxversionen verstehen und verwalten.

    [Versionsanforderungen prüfen :material-arrow-right:](info/nxf_versions.md){ .md-button .md-button--secondary }

-   :material-file-code-outline:{ .lg .middle } __Die Hello-Pipeline__

    ---

    Überblick darüber, was die Hello-Pipeline macht und wie sie aufgebaut ist.

    [Überblick lesen :material-arrow-right:](info/hello_pipeline.md){ .md-button .md-button--secondary }

-   :material-lifebuoy:{ .lg .middle } __Hilfe erhalten__

    ---

    Hilfreiche Ressourcen, wenn du ein Problem mit dem Nextflow-Training hast.

    [Hilfe finden :material-arrow-right:](help.md){ .md-button .md-button--secondary }

</div>
