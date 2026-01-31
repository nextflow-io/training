# Manuelle Installation

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } KI-gestützte Übersetzung - [mehr erfahren & Verbesserungen vorschlagen](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Es ist möglich, alles, was du zum Ausführen des Trainings benötigst, manuell in deiner eigenen lokalen Umgebung zu installieren.

Hier haben wir dokumentiert, wie das auf Standard-POSIX-kompatiblen Systemen geht (unter der Annahme eines persönlichen Rechners wie eines Laptops).
Beachte, dass einige Details je nach deinem spezifischen System unterschiedlich sein können.

!!! tip "Tipp"

    Bevor du fortfährst, hast du den [Devcontainers-Ansatz](03_devcontainer.md) in Betracht gezogen?
    Er stellt alle notwendigen Tools und Abhängigkeiten bereit, ohne dass eine manuelle Installation erforderlich ist.

## Allgemeine Softwareanforderungen

Nextflow kann auf jedem POSIX-kompatiblen System (Linux, macOS, Windows Subsystem for Linux usw.) mit installiertem Java verwendet werden.
Unsere Trainingskurse haben einige zusätzliche Anforderungen.

Insgesamt musst du die folgende Software installiert haben:

- Bash oder eine gleichwertige Shell
- [Java 11 (oder höher, bis 21)](https://www.oracle.com/technetwork/java/javase/downloads/index.html)
- [Git](https://git-scm.com/)
- [Docker](https://docs.docker.com/get-docker/)
- [Conda](https://conda.io/) 4.5 (oder höher)
- [VSCode](https://code.visualstudio.com) mit der [Nextflow-Erweiterung](https://www.nextflow.io/docs/latest/developer-env.html#devenv-nextflow)

Die VSCode-Anwendung ist technisch optional, aber wir empfehlen dringend, sie sowohl für das Durcharbeiten der Kurse als auch für deine Nextflow-Entwicklungsarbeit im Allgemeinen zu verwenden.

Das Nextflow-Dokumentationshandbuch enthält Anweisungen zur Installation dieser Abhängigkeiten unter [Umgebungseinrichtung](https://www.nextflow.io/docs/latest/developer-env.html).

## Nextflow und nf-core tools

Du musst Nextflow selbst sowie die nf-core tools installieren, wie in den unten verlinkten Artikeln beschrieben:

- [Nextflow-Installation](https://www.nextflow.io/docs/latest/install.html)
- [nf-core tools](https://nf-co.re/docs/nf-core-tools/installation)

Wir empfehlen die Self-Install-Option für Nextflow und die PyPI-Option für nf-core tools.

!!! warning "Versionskompatibilität"

    <!-- Any update to this content needs to be copied to the home page -->
    **Ab Januar 2026 erfordern alle unsere Nextflow-Trainingskurse Nextflow Version 25.10.2 oder höher mit aktivierter strikter v2-Syntax, sofern nicht anders angegeben.**

    Weitere Informationen zu Versionsanforderungen und strikter v2-Syntax findest du im [Nextflow-Versionen](../info/nxf_versions.md)-Leitfaden.

    Ältere Versionen des Trainingsmaterials, die früherer Syntax entsprechen, sind über den Versionsauswähler in der Menüleiste dieser Webseite verfügbar.

## Trainingsmaterialien

Der einfachste Weg, die Trainingsmaterialien herunterzuladen, ist das Klonen des gesamten Repositorys mit diesem Befehl:

```bash
git clone https://github.com/nextflow-io/training.git
```

Jeder Kurs hat sein eigenes Verzeichnis.
Um einen Kurs durchzuarbeiten, öffne ein Terminalfenster (idealerweise innerhalb der VSCode-Anwendung) und wechsle mit `cd` in das entsprechende Verzeichnis.

Du kannst dann den Kursanweisungen auf der Website folgen.
