# Devcontainers locals

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducció assistida per IA - [més informació i suggeriments](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Si teniu una instal·lació local de Docker o esteu disposats a instal·lar-ne una, la manera més fàcil de treballar localment amb aquests materials és utilitzar la funcionalitat de devcontainer de Visual Studio Code. Aquest enfocament proporciona totes les eines i dependències necessàries sense requerir instal·lació manual.

## Requisits

Per utilitzar la configuració de devcontainer local, necessitareu:

- [Visual Studio Code](https://code.visualstudio.com/)
- Una instal·lació local de Docker, per exemple:
  - [Docker Desktop](https://docs.docker.com/get-docker/) (per a Windows/macOS)
  - [Docker Engine](https://docs.docker.com/engine/install/) (per a Linux)
  - [Colima](https://github.com/abiosoft/colima) (alternativa per a macOS)
- [Docker Buildx](https://docs.docker.com/build/concepts/overview/#install-buildx) (inclòs a Docker Desktop, però pot necessitar instal·lació separada amb altres configuracions de Docker)
- [Extensió Dev Containers](https://marketplace.visualstudio.com/items?itemName=ms-vscode-remote.remote-containers) per a VS Code

La vostra instal·lació de Docker ha d'estar en execució abans d'intentar obrir el devcontainer.

Per verificar que Docker buildx està disponible, executeu:

```bash
docker buildx version
```

Si aquesta comanda falla, haureu d'instal·lar l'extensió buildx abans de continuar.

## Instruccions de configuració

Seguiu aquests passos per configurar el vostre entorn local utilitzant devcontainers de VS Code:

### Instal·leu l'extensió "Dev Containers" a VS Code

- Obriu VS Code
- Aneu a Extensions (Ctrl+Shift+X o Cmd+Shift+X a macOS)
- Cerqueu "Dev Containers"
- Feu clic a "Install"

![Instal·lant l'extensió Dev Containers a VS Code](img/install_extension.png)

### Cloneu el repositori:

```bash
git clone https://github.com/nextflow-io/training.git
cd training
```

### Obriu el repositori a VS Code:

- Inicieu VS Code
- Seleccioneu **File -> Open Folder** del menú
- Navegueu fins a la carpeta del repositori de formació que acabeu de clonar i seleccioneu-la
- Feu clic a **Open**

### Reobriu en contenidor

Si VS Code us demana "Reopen in Container", feu-hi clic. Alternativament:

- Premeu F1 (o Ctrl+Shift+P / Cmd+Shift+P a macOS)
- Escriviu "Dev Containers: Reopen in Container"
- **Important**: Quan se us demani seleccionar una configuració, trieu la configuració de devcontainer **local-dev**

![Avís de Reopen in Container](img/reopen_prompt.png)

![Seleccionant la configuració local](img/select_local_config.png)

Espereu que es construeixi el contenidor. Això pot trigar uns minuts la primera vegada, ja que descarrega i configura tots els components necessaris.

Un cop el contenidor estigui construït i en execució, tindreu un entorn completament configurat amb totes les eines necessàries instal·lades, incloent:

- Java
- Nextflow
- Docker
- Git
- I totes les altres dependències necessàries per a la formació

![VS Code amb devcontainer en execució](img/running_container.png)

## Avantatges d'utilitzar devcontainers

L'ús de l'enfocament de devcontainer ofereix diversos avantatges:

- **Consistència**: Garanteix un entorn de desenvolupament consistent entre diferents màquines
- **Simplicitat**: Totes les dependències estan preinstal·lades i configurades
- **Aïllament**: L'entorn de desenvolupament està aïllat del vostre sistema local
- **Reproductibilitat**: Tothom que utilitzi el devcontainer obté la mateixa configuració
- **Sense instal·lació manual**: No cal instal·lar manualment Java, Nextflow i altres eines

## Comprovació del vostre entorn

Un cop el vostre devcontainer estigui en execució, podeu verificar que tot està configurat correctament executant:

```bash
nextflow info
```

Això hauria de mostrar la versió de Nextflow i la informació d'execució, confirmant que el vostre entorn està configurat correctament.

## Resolució de problemes

Si trobeu problemes amb la configuració del devcontainer:

1. Assegureu-vos que la vostra instal·lació de Docker (Docker Desktop, Colima, Docker Engine, etc.) està en execució abans d'obrir el devcontainer
2. Comproveu que heu seleccionat la configuració **local-dev** quan se us demani
3. Verifiqueu que Docker buildx està instal·lat i funciona executant `docker buildx version`
4. Si el contenidor no es construeix, proveu de reconstruir-lo executant la comanda "Dev Containers: Rebuild Container"
5. Per a problemes persistents, consulteu la [guia de resolució de problemes de VS Code Dev Containers](https://code.visualstudio.com/docs/devcontainers/troubleshooting)
