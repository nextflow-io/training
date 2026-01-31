# Lokale Devcontainers

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } KI-gestützte Übersetzung - [mehr erfahren & Verbesserungen vorschlagen](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Wenn du eine lokale Docker-Installation hast oder bereit bist, eine zu installieren, ist der einfachste Weg, lokal mit diesen Materialien zu arbeiten, die Devcontainer-Funktion von Visual Studio Code zu verwenden. Dieser Ansatz stellt alle notwendigen Tools und Abhängigkeiten bereit, ohne dass eine manuelle Installation erforderlich ist.

## Anforderungen

Um das lokale Devcontainer-Setup zu verwenden, benötigst du:

- [Visual Studio Code](https://code.visualstudio.com/)
- Eine lokale Docker-Installation, zum Beispiel:
  - [Docker Desktop](https://docs.docker.com/get-docker/) (für Windows/macOS)
  - [Docker Engine](https://docs.docker.com/engine/install/) (für Linux)
  - [Colima](https://github.com/abiosoft/colima) (Alternative für macOS)
- [Docker Buildx](https://docs.docker.com/build/concepts/overview/#install-buildx) (in Docker Desktop enthalten, kann aber bei anderen Docker-Setups eine separate Installation erfordern)
- [Dev Containers-Erweiterung](https://marketplace.visualstudio.com/items?itemName=ms-vscode-remote.remote-containers) für VS Code

Deine Docker-Installation muss laufen, bevor du versuchst, den Devcontainer zu öffnen.

Um zu überprüfen, ob Docker buildx verfügbar ist, führe aus:

```bash
docker buildx version
```

Wenn dieser Befehl fehlschlägt, musst du die buildx-Erweiterung installieren, bevor du fortfährst.

## Einrichtungsanleitung

Folge diesen Schritten, um deine lokale Umgebung mit VS Code Devcontainers einzurichten:

### Installiere die "Dev Containers"-Erweiterung in VS Code

- Öffne VS Code
- Gehe zu Extensions (Strg+Umschalt+X oder Cmd+Umschalt+X auf macOS)
- Suche nach "Dev Containers"
- Klicke auf "Install"

![Installation der Dev Containers-Erweiterung in VS Code](img/install_extension.png)

### Klone das Repository:

```bash
git clone https://github.com/nextflow-io/training.git
cd training
```

### Öffne das Repository in VS Code:

- Starte VS Code
- Wähle **Datei -> Ordner öffnen** aus dem Menü
- Navigiere zum gerade geklonten Training-Repository-Ordner und wähle ihn aus
- Klicke auf **Öffnen**

### In Container erneut öffnen

Wenn VS Code dich auffordert, "In Container erneut öffnen", klicke darauf. Alternativ:

- Drücke F1 (oder Strg+Umschalt+P / Cmd+Umschalt+P auf macOS)
- Tippe "Dev Containers: Reopen in Container"
- **Wichtig**: Wenn du aufgefordert wirst, eine Konfiguration auszuwählen, wähle die **local-dev** Devcontainer-Konfiguration

![Aufforderung zum erneuten Öffnen im Container](img/reopen_prompt.png)

![Auswahl der lokalen Konfiguration](img/select_local_config.png)

Warte, bis der Container erstellt ist. Dies kann beim ersten Mal einige Minuten dauern, da alle notwendigen Komponenten heruntergeladen und eingerichtet werden.

Sobald der Container erstellt und ausgeführt wird, hast du eine vollständig konfigurierte Umgebung mit allen installierten notwendigen Tools, einschließlich:

- Java
- Nextflow
- Docker
- Git
- Und alle anderen für das Training erforderlichen Abhängigkeiten

![VS Code mit laufendem Devcontainer](img/running_container.png)

## Vorteile der Verwendung von Devcontainers

Die Verwendung des Devcontainer-Ansatzes bietet mehrere Vorteile:

- **Konsistenz**: Stellt eine konsistente Entwicklungsumgebung über verschiedene Maschinen hinweg sicher
- **Einfachheit**: Alle Abhängigkeiten sind vorinstalliert und konfiguriert
- **Isolation**: Die Entwicklungsumgebung ist von deinem lokalen System isoliert
- **Reproduzierbarkeit**: Jeder, der den Devcontainer verwendet, erhält dasselbe Setup
- **Keine manuelle Installation**: Keine Notwendigkeit, Java, Nextflow und andere Tools manuell zu installieren

## Überprüfung deiner Umgebung

Sobald dein Devcontainer läuft, kannst du überprüfen, ob alles korrekt eingerichtet ist, indem du ausführst:

```bash
nextflow info
```

Dies sollte die Nextflow-Version und Laufzeitinformationen anzeigen und bestätigen, dass deine Umgebung korrekt konfiguriert ist.

## Fehlerbehebung

Wenn du Probleme mit dem Devcontainer-Setup hast:

1. Stelle sicher, dass deine Docker-Installation (Docker Desktop, Colima, Docker Engine usw.) läuft, bevor du den Devcontainer öffnest
2. Überprüfe, dass du die **local-dev**-Konfiguration ausgewählt hast, wenn du dazu aufgefordert wurdest
3. Überprüfe, ob Docker buildx installiert ist und funktioniert, indem du `docker buildx version` ausführst
4. Wenn der Container nicht erstellt werden kann, versuche ihn neu zu erstellen, indem du den Befehl "Dev Containers: Rebuild Container" ausführst
5. Bei anhaltenden Problemen siehe den [VS Code Dev Containers-Fehlerbehebungsleitfaden](https://code.visualstudio.com/docs/devcontainers/troubleshooting)
