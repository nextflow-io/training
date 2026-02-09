# Primi passi

## Avviare un ambiente di formazione

Per utilizzare l'ambiente preconfigurato che forniamo su GitHub Codespaces, cliccate sul pulsante "Open in GitHub Codespaces" qui sotto. Per altre opzioni, consultate [Opzioni di ambiente](../envsetup/index.md).

Vi consigliamo di aprire l'ambiente di formazione in una nuova scheda o finestra del browser (utilizzate il clic destro, ctrl-clic o cmd-clic a seconda del vostro dispositivo) in modo da poter continuare a leggere mentre l'ambiente si carica.
Dovrete tenere queste istruzioni aperte in parallelo per seguire il corso.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

### Nozioni di base sull'ambiente

Questo ambiente di formazione contiene tutto il software, il codice e i dati necessari per seguire il corso di formazione, quindi non dovete installare nulla voi stessi.

Il codespace è configurato con un'interfaccia VSCode, che include un esploratore del filesystem, un editor di codice e una shell del terminale.
Tutte le istruzioni fornite durante il corso (ad es. 'aprite il file', 'modificate il codice' o 'eseguite questo comando') si riferiscono a queste tre parti dell'interfaccia VSCode salvo diversa indicazione.

Se state seguendo questo corso da soli, vi preghiamo di familiarizzare con le [nozioni di base sull'ambiente](../envsetup/01_setup.md) per ulteriori dettagli.

### Requisiti di versione

Questa formazione è progettata per Nextflow 25.10.2 o successivo **con il parser di sintassi v2 ABILITATO**.
Se state utilizzando un ambiente locale o personalizzato, assicuratevi di utilizzare le impostazioni corrette come documentato [qui](../info/nxf_versions.md).

## Preparatevi a lavorare

Una volta che il vostro codespace è in esecuzione, ci sono due cose che dovete fare prima di immergervi nella formazione: impostare la directory di lavoro per questo corso specifico e dare un'occhiata ai materiali forniti.

### Impostare la directory di lavoro

Per impostazione predefinita, il codespace si apre con la directory di lavoro impostata alla radice di tutti i corsi di formazione, ma per questo corso lavoreremo nella directory `nextflow-run/`.

Cambiate directory ora eseguendo questo comando nel terminale:

```bash
cd nextflow-run/
```

Potete impostare VSCode per concentrarsi su questa directory, in modo che solo i file pertinenti vengano visualizzati nella barra laterale dell'esploratore file:

```bash
code .
```

!!! tip "Suggerimento"

    Se per qualsiasi motivo vi spostate da questa directory (ad es. il vostro codespace va in sospensione), potete sempre utilizzare il percorso completo per ritornarvi, supponendo che stiate eseguendo questo nell'ambiente di formazione Github Codespaces:

    ```bash
    cd /workspaces/training/nextflow-run
    ```

Ora diamo un'occhiata ai contenuti.

### Esplorare i materiali forniti

Potete esplorare i contenuti di questa directory utilizzando l'esploratore file sul lato sinistro dell'area di lavoro di formazione.
In alternativa, potete utilizzare il comando `tree`.

Durante il corso, utilizziamo l'output di `tree` per rappresentare la struttura e i contenuti della directory in una forma leggibile, a volte con modifiche minori per chiarezza.

Qui generiamo un indice dei contenuti fino al secondo livello:

```bash
tree . -L 2
```

??? abstract "Directory contents"

    ```console
    .
    ├── 1-hello.nf
    ├── 2a-inputs.nf
    ├── 2b-multistep.nf
    ├── 2c-modules.nf
    ├── 2d-container.nf
    ├── 3-main.nf
    ├── data
    │   └── greetings.csv
    ├── modules
    │   ├── collectGreetings.nf
    │   ├── convertToUpper.nf
    │   ├── cowpy.nf
    │   └── sayHello.nf
    ├── nextflow.config
    ├── solutions
    │   ├── 3-main.nf
    │   ├── modules
    │   └── nextflow.config
    ├── test-params.json
    └── test-params.yaml
    ```

Cliccate sulla casella colorata per espandere la sezione e visualizzarne i contenuti.
Utilizziamo sezioni comprimibili come questa per visualizzare l'output previsto dei comandi così come i contenuti di directory e file in modo conciso.

- **I file `.nf`** sono script di flusso di lavoro numerati in base alla parte del corso in cui vengono utilizzati.

- **Il file `nextflow.config`** è un file di configurazione che imposta proprietà minime dell'ambiente.
  Potete ignorarlo per ora.

- **Il file `greetings.csv`** sotto `data/` contiene i dati di input che utilizzeremo nella maggior parte del corso. È descritto nella Parte 2 (Eseguire pipeline), quando lo introduciamo per la prima volta.

- **I file `test-params.*`** sono file di configurazione che utilizzeremo nella Parte 3 (Configurazione). Potete ignorarli per ora.

- **La directory `solutions`** contiene lo stato finale del flusso di lavoro e dei suoi file accessori (config e moduli) che risultano dal completamento del corso.
  Sono destinati ad essere utilizzati come riferimento per verificare il vostro lavoro e risolvere eventuali problemi.

## Lista di controllo della preparazione

Pensate di essere pronti per iniziare?

- [ ] Comprendo l'obiettivo di questo corso e i suoi prerequisiti
- [ ] Il mio ambiente è attivo e funzionante
- [ ] Ho impostato la mia directory di lavoro in modo appropriato

Se potete spuntare tutte le caselle, siete pronti per partire.

**Per continuare alla [Parte 1: Eseguire operazioni di base](./01_basics.md), cliccate sulla freccia nell'angolo in basso a destra di questa pagina.**
