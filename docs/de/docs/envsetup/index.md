# Umgebungsoptionen

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } KI-gestützte Übersetzung - [mehr erfahren & Verbesserungen vorschlagen](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Wir möchten eine konsistente und gründlich getestete Umgebung bereitstellen, die es dir ermöglicht, dich auf das Erlernen von Nextflow zu konzentrieren, ohne Zeit und Mühe für die Verwaltung von Software aufwenden zu müssen.
Zu diesem Zweck haben wir eine containerisierte Umgebung entwickelt, die alle notwendige Software, Code-Dateien und Beispieldaten enthält, um alle unsere Kurse durchzuarbeiten.

Diese containerisierte Umgebung kann sofort auf Github Codespaces oder lokal in VS Code mit der Devcontainers-Erweiterung ausgeführt werden.

<div class="grid cards" markdown>

-   :material-cloud-outline:{ .lg .middle } __Github Codespaces__

    ---

    GitHub Codespaces ist ein webbasierter Dienst, der es uns ermöglicht, eine vorgefertigte Umgebung für das Training bereitzustellen, mit allen Tools und Daten inklusive, unterstützt durch virtuelle Maschinen in der Cloud. Der Zugang ist kostenlos für alle mit einem Github-Konto.

    [Github Codespaces verwenden:material-arrow-right:](01_setup.md){ .md-button .md-button--primary .mt-1 }

-   :material-laptop:{ .lg .middle } __Lokale Devcontainers__

    ---

    VS Code mit Devcontainers bietet eine lokal ausgeführte containerisierte Entwicklungsumgebung mit allen vorkonfigurierten Training-Tools. Es bietet die gleiche vorgefertigte Umgebung wie Codespaces, läuft aber vollständig auf deiner lokalen Hardware.

    [Devcontainers lokal verwenden :material-arrow-right:](03_devcontainer.md){ .md-button .md-button--primary .mt-1 }

</div>

## Anleitung für manuelle Installation

Wenn keine der oben genannten Optionen deinen Anforderungen entspricht, kannst du diese Umgebung auf deinem eigenen lokalen System replizieren, indem du die Software-Abhängigkeiten manuell installierst und das Training-Repository klonst.

[Manuelle Installation :material-arrow-right:](02_local.md){ .md-button .md-button--primary .mt-1 }

---

!!! info "Einstellung von Gitpod"

    Nextflow Training nutzte bis Februar 2025 [Gitpod](https://gitpod.io).
    Die Entwickler\*innen von Gitpod haben jedoch beschlossen, die kostenlose Funktionalität zugunsten des [Gitpod Flex](https://www.gitpod.io/blog/introducing-gitpod-flex)-Systems einzustellen.
    Aus diesem Grund sind wir zu GitHub Codespaces gewechselt, das ebenfalls eine Entwicklungsumgebung mit einem Klick und ohne vorherige Einrichtung bietet.

    Je nachdem, wann du dich bei Gitpod angemeldet hast und wann genau der Dienst eingestellt wird, kannst du das Training möglicherweise noch in der alten Cloud-IDE starten, obwohl wir keinen zuverlässigen Zugang in Zukunft garantieren können:
    [In Gitpod öffnen](https://gitpod.io/#https://github.com/nextflow-io/training).
