# Lokale Devcontainer

Wenn du eine lokale Docker-Installation hast oder bereit bist, eine zu installieren, ist die einfachste Möglichkeit, lokal mit diesen Materialien zu arbeiten, die Devcontainer-Funktion von Visual Studio Code zu verwenden. Dieser Ansatz stellt alle notwendigen Tools und Abhängigkeiten bereit, ohne dass eine manuelle Installation erforderlich ist.

## Voraussetzungen

Um das lokale Devcontainer-Setup zu verwenden, benötigst du:

- [Visual Studio Code](https://code.visualstudio.com/)
- Eine lokale Docker-Installation, zum Beispiel:
  - [Docker Desktop](https://docs.docker.com/get-docker/) (für Windows/macOS)
  - [Docker Engine](https://docs.docker.com/engine/install/) (für Linux)
  - [Colima](https://github.com/abiosoft/colima) (Alternative für macOS)
- [Docker Buildx](https://docs.docker.com/build/concepts/overview/#install-buildx) (in Docker Desktop enthalten, bei anderen Docker-Setups möglicherweise separate Installation erforderlich)
- [Dev Containers Extension](https://marketplace.visualstudio.com/items?itemName=ms-vscode-remote.remote-containers) für VS Code

Deine Docker-Installation muss laufen, bevor du versuchst, den Devcontainer zu öffnen.

Um zu überprüfen, ob Docker buildx verfügbar ist, führe aus:

```bash
docker buildx version
```

Wenn dieser Befehl fehlschlägt, musst du die buildx-Erweiterung installieren, bevor du fortfährst.

## Einrichtungsanleitung

Folge diesen Schritten, um deine lokale Umgebung mit VS Code Devcontainern einzurichten:

### Installiere die "Dev Containers" Extension in VS Code

- Öffne VS Code
- Gehe zu Extensions (Strg+Umschalt+X oder Cmd+Umschalt+X auf macOS)
- Suche nach "Dev Containers"
- Klicke auf "Install"

![Installation der Dev Containers Extension in VS Code](img/install_extension.png)

### Klone das Repository:

```bash
git clone https://github.com/nextflow-io/training.git
cd training
```

### Öffne das Repository in VS Code:

- Starte VS Code
- Wähle **File -> Open Folder** aus dem Menü
- Navigiere zum gerade geklonten training-Repository-Ordner und wähle ihn aus
- Klicke auf **Open**

### Im Container erneut öffnen

Wenn VS Code dich auffordert, "Reopen in Container" auszuwählen, klicke darauf. Alternativ:

- Drücke F1 (oder Strg+Umschalt+P / Cmd+Umschalt+P auf macOS)
- Tippe "Dev Containers: Reopen in Container"
- **Wichtig**: Wenn du aufgefordert wirst, eine Konfiguration auszuwählen, wähle die **local-dev** Devcontainer-Konfiguration

![Reopen in Container Aufforderung](img/reopen_prompt.png)

![Auswahl der lokalen Konfiguration](img/select_local_config.png)

Warte, bis der Container gebaut ist. Dies kann beim ersten Mal einige Minuten dauern, da alle notwendigen Komponenten heruntergeladen und eingerichtet werden.

Sobald der Container gebaut ist und läuft, hast du eine vollständig konfigurierte Umgebung mit allen notwendigen installierten Tools, einschließlich:

- Java
- Nextflow
- Docker
- Git
- Und allen anderen für das Training erforderlichen Abhängigkeiten

![VS Code mit laufendem Devcontainer](img/running_container.png)

## Vorteile der Verwendung von Devcontainern

Die Verwendung des Devcontainer-Ansatzes bietet mehrere Vorteile:

- **Konsistenz**: Gewährleistet eine konsistente Entwicklungsumgebung auf verschiedenen Rechnern
- **Einfachheit**: Alle Abhängigkeiten sind vorinstalliert und konfiguriert
- **Isolation**: Die Entwicklungsumgebung ist von deinem lokalen System isoliert
- **Reproduzierbarkeit**: Jede\*r, der/die den Devcontainer verwendet, erhält das gleiche Setup
- **Keine manuelle Installation**: Keine Notwendigkeit, Java, Nextflow und andere Tools manuell zu installieren

## Überprüfung deiner Umgebung

Sobald dein Devcontainer läuft, kannst du überprüfen, ob alles korrekt eingerichtet ist, indem du ausführst:

```bash
nextflow info
```

Dies sollte die Nextflow-Version und Laufzeitinformationen anzeigen und bestätigen, dass deine Umgebung ordnungsgemäß konfiguriert ist.

## Fehlerbehebung

Wenn du Probleme mit dem Devcontainer-Setup hast:

1. Stelle sicher, dass deine Docker-Installation (Docker Desktop, Colima, Docker Engine usw.) läuft, bevor du den Devcontainer öffnest
2. Überprüfe, dass du die **local-dev** Konfiguration ausgewählt hast, als du dazu aufgefordert wurdest
3. Verifiziere, dass Docker buildx installiert ist und funktioniert, indem du `docker buildx version` ausführst
4. Wenn der Container nicht gebaut werden kann, versuche ihn neu zu bauen, indem du den Befehl "Dev Containers: Rebuild Container" ausführst
5. Bei anhaltenden Problemen siehe den [VS Code Dev Containers Fehlerbehebungsleitfaden](https://code.visualstudio.com/docs/devcontainers/troubleshooting)
