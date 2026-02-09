# Devcontainer Locali

Se avete un'installazione locale di Docker o siete disposti a installarne una, il modo più semplice per lavorare localmente con questi materiali è utilizzare la funzionalità devcontainer di Visual Studio Code. Questo approccio fornisce tutti gli strumenti e le dipendenze necessarie senza richiedere un'installazione manuale.

## Requisiti

Per utilizzare la configurazione devcontainer locale, avrete bisogno di:

- [Visual Studio Code](https://code.visualstudio.com/)
- Un'installazione locale di Docker, per esempio:
  - [Docker Desktop](https://docs.docker.com/get-docker/) (per Windows/macOS)
  - [Docker Engine](https://docs.docker.com/engine/install/) (per Linux)
  - [Colima](https://github.com/abiosoft/colima) (alternativa per macOS)
- [Docker Buildx](https://docs.docker.com/build/concepts/overview/#install-buildx) (incluso in Docker Desktop, ma potrebbe richiedere un'installazione separata con altre configurazioni Docker)
- [Estensione Dev Containers](https://marketplace.visualstudio.com/items?itemName=ms-vscode-remote.remote-containers) per VS Code

La vostra installazione Docker deve essere in esecuzione prima di tentare di aprire il devcontainer.

Per verificare che Docker buildx sia disponibile, eseguite:

```bash
docker buildx version
```

Se questo comando fallisce, dovrete installare l'estensione buildx prima di procedere.

## Istruzioni per la Configurazione

Seguite questi passaggi per configurare il vostro ambiente locale utilizzando i devcontainer di VS Code:

### Installare l'estensione "Dev Containers" in VS Code

- Aprite VS Code
- Andate su Estensioni (Ctrl+Shift+X o Cmd+Shift+X su macOS)
- Cercate "Dev Containers"
- Cliccate su "Install"

![Installazione dell'estensione Dev Containers in VS Code](img/install_extension.png)

### Clonare il repository:

```bash
git clone https://github.com/nextflow-io/training.git
cd training
```

### Aprire il repository in VS Code:

- Avviate VS Code
- Selezionate **File -> Open Folder** dal menu
- Navigate fino alla cartella del repository training che avete appena clonato e selezionatela
- Cliccate su **Open**

### Riaprire nel Container

Se VS Code vi chiede di "Reopen in Container", cliccate su di esso. In alternativa:

- Premete F1 (o Ctrl+Shift+P / Cmd+Shift+P su macOS)
- Digitate "Dev Containers: Reopen in Container"
- **Importante**: Quando vi viene chiesto di selezionare una configurazione, scegliete la configurazione devcontainer **local-dev**

![Prompt Reopen in Container](img/reopen_prompt.png)

![Selezione della configurazione locale](img/select_local_config.png)

Attendete che il container venga costruito. La prima volta potrebbe richiedere alcuni minuti poiché scarica e configura tutti i componenti necessari.

Una volta che il container è costruito e in esecuzione, avrete un ambiente completamente configurato con tutti gli strumenti necessari installati, tra cui:

- Java
- Nextflow
- Docker
- Git
- E tutte le altre dipendenze richieste per la formazione

![VS Code con devcontainer in esecuzione](img/running_container.png)

## Vantaggi dell'Utilizzo dei Devcontainer

L'utilizzo dell'approccio devcontainer offre diversi vantaggi:

- **Coerenza**: Garantisce un ambiente di sviluppo coerente su macchine diverse
- **Semplicità**: Tutte le dipendenze sono preinstallate e configurate
- **Isolamento**: L'ambiente di sviluppo è isolato dal vostro sistema locale
- **Riproducibilità**: Tutti coloro che utilizzano il devcontainer ottengono la stessa configurazione
- **Nessuna installazione manuale**: Non è necessario installare manualmente Java, Nextflow e altri strumenti

## Verifica del Vostro Ambiente

Una volta che il vostro devcontainer è in esecuzione, potete verificare che tutto sia configurato correttamente eseguendo:

```bash
nextflow info
```

Questo dovrebbe visualizzare la versione di Nextflow e le informazioni di runtime, confermando che il vostro ambiente è configurato correttamente.

## Risoluzione dei Problemi

Se riscontrate problemi con la configurazione del devcontainer:

1. Assicuratevi che la vostra installazione Docker (Docker Desktop, Colima, Docker Engine, ecc.) sia in esecuzione prima di aprire il devcontainer
2. Verificate di aver selezionato la configurazione **local-dev** quando richiesto
3. Verificate che Docker buildx sia installato e funzionante eseguendo `docker buildx version`
4. Se il container non riesce a costruirsi, provate a ricostruirlo eseguendo il comando "Dev Containers: Rebuild Container"
5. Per problemi persistenti, fate riferimento alla [guida alla risoluzione dei problemi di VS Code Dev Containers](https://code.visualstudio.com/docs/devcontainers/troubleshooting)
