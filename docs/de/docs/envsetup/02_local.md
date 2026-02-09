# Manuelle Installation

Es ist möglich, alles, was du für das Training benötigst, manuell in deiner eigenen lokalen Umgebung zu installieren.

Hier haben wir dokumentiert, wie das auf Standard-POSIX-kompatiblen Systemen funktioniert (ausgehend von einem persönlichen Rechner wie einem Laptop).
Beachte, dass einige Details je nach deinem spezifischen System unterschiedlich sein können.

!!! tip

    Bevor du fortfährst: Hast du die [Devcontainer-Methode](03_devcontainer.md) in Betracht gezogen?
    Sie stellt alle notwendigen Tools und Abhängigkeiten bereit, ohne dass eine manuelle Installation erforderlich ist.

## Allgemeine Softwareanforderungen

Nextflow kann auf jedem POSIX-kompatiblen System (Linux, macOS, Windows Subsystem for Linux usw.) mit installiertem Java verwendet werden.
Unsere Trainingskurse haben einige zusätzliche Anforderungen.

Insgesamt musst du folgende Software installiert haben:

- Bash oder eine gleichwertige Shell
- [Java 11 (oder neuer, bis Version 21)](https://www.oracle.com/technetwork/java/javase/downloads/index.html)
- [Git](https://git-scm.com/)
- [Docker](https://docs.docker.com/get-docker/)
- [Conda](https://conda.io/) 4.5 (oder neuer)
- [VSCode](https://code.visualstudio.com) mit der [Nextflow-Erweiterung](https://www.nextflow.io/docs/latest/developer-env.html#devenv-nextflow)

Die VSCode-Anwendung ist technisch gesehen optional, aber wir empfehlen dir dringend, sie sowohl für die Kurse als auch für deine allgemeine Nextflow-Entwicklungsarbeit zu verwenden.

Das Nextflow-Dokumentationshandbuch bietet Anleitungen zur Installation dieser Abhängigkeiten unter [Environment setup](https://www.nextflow.io/docs/latest/developer-env.html).

## Nextflow und nf-core Tools

Du musst Nextflow selbst sowie die nf-core Tools installieren, wie in den unten verlinkten Artikeln beschrieben:

- [Nextflow installation](https://www.nextflow.io/docs/latest/install.html)
- [nf-core tools](https://nf-co.re/docs/nf-core-tools/installation)

Wir empfehlen die Selbstinstallationsoption für Nextflow und die PyPI-Option für nf-core Tools.

!!! warning "Versionskompatibilität"

    <!-- Any update to this content needs to be copied to the home page -->
    **Ab Januar 2026 erfordern alle unsere Nextflow-Trainingskurse Nextflow Version 25.10.2 oder neuer mit aktivierter strikter v2-Syntax, sofern nicht anders angegeben.**

    Weitere Informationen zu Versionsanforderungen und strikter v2-Syntax findest du im Leitfaden [Nextflow versions](../info/nxf_versions.md).

    Ältere Versionen des Trainingsmaterials, die früheren Syntaxversionen entsprechen, sind über den Versionsauswähler in der Menüleiste dieser Webseite verfügbar.

## Trainingsmaterialien

Der einfachste Weg, die Trainingsmaterialien herunterzuladen, ist das Klonen des gesamten Repositorys mit diesem Befehl:

```bash
git clone https://github.com/nextflow-io/training.git
```

Jeder Kurs hat sein eigenes Verzeichnis.
Um einen Kurs durchzuarbeiten, öffne ein Terminalfenster (idealerweise innerhalb der VSCode-Anwendung) und wechsle mit `cd` in das entsprechende Verzeichnis.

Du kannst dann den Kursanweisungen folgen, die auf der Website bereitgestellt werden.
