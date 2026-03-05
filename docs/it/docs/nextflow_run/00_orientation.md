# Per iniziare

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduzione assistita da IA - [scopri di più e suggerisci miglioramenti](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

## Avvia un ambiente di formazione

Per utilizzare l'ambiente pre-costruito che forniamo su GitHub Codespaces, clicca sul pulsante "Open in GitHub Codespaces" qui sotto. Per altre opzioni, consulta [Opzioni di ambiente](../envsetup/index.md).

Ti consigliamo di aprire l'ambiente di formazione in una nuova scheda o finestra del browser (usa clic destro, ctrl-clic o cmd-clic a seconda del tuo dispositivo) in modo da poter continuare a leggere mentre l'ambiente si carica.
Dovrai tenere queste istruzioni aperte in parallelo per seguire il corso.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

### Nozioni base sull'ambiente

Questo ambiente di formazione contiene tutto il software, il codice e i dati necessari per seguire il corso, quindi non devi installare nulla.

Il codespace è configurato con un'interfaccia VSCode, che include un esploratore del filesystem, un editor di codice e un terminale shell.
Tutte le istruzioni fornite durante il corso (es. 'apri il file', 'modifica il codice' o 'esegui questo comando') si riferiscono a queste tre parti dell'interfaccia VSCode salvo diversa indicazione.

Se stai seguendo questo corso da solo, ti preghiamo di familiarizzare con le [nozioni base sull'ambiente](../envsetup/01_setup.md) per ulteriori dettagli.

### Requisiti di versione

Questa formazione è progettata per Nextflow 25.10.2 o successivo **con il parser della sintassi v2 ABILITATO**.
Se stai usando un ambiente locale o personalizzato, assicurati di utilizzare le impostazioni corrette come documentato [qui](../info/nxf_versions.md).

## Preparati al lavoro

Una volta che il tuo codespace è in esecuzione, ci sono due cose da fare prima di immergerti nella formazione: impostare la directory di lavoro per questo corso specifico e dare un'occhiata ai materiali forniti.

### Imposta la directory di lavoro

Per impostazione predefinita, il codespace si apre con la directory di lavoro impostata alla radice di tutti i corsi di formazione, ma per questo corso lavoreremo nella directory `nextflow-run/`.

Cambia directory ora eseguendo questo comando nel terminale:

```bash
cd nextflow-run/
```

Puoi impostare VSCode per concentrarsi su questa directory, in modo che solo i file rilevanti vengano mostrati nella barra laterale dell'esploratore file:

```bash
code .
```

!!! tip "Suggerimento"

    Se per qualsiasi motivo esci da questa directory (es. il tuo codespace va in sospensione), puoi sempre usare il percorso completo per tornare, assumendo che tu stia eseguendo questo all'interno dell'ambiente di formazione di Github Codespaces:

    ```bash
    cd /workspaces/training/nextflow-run
    ```

Ora diamo un'occhiata ai contenuti.

### Esplora i materiali forniti

Puoi esplorare i contenuti di questa directory usando l'esploratore file sul lato sinistro dello spazio di lavoro di formazione.
In alternativa, puoi usare il comando `tree`.

Durante il corso, usiamo l'output di `tree` per rappresentare la struttura e i contenuti delle directory in forma leggibile, a volte con modifiche minori per chiarezza.

Qui generiamo un indice dei contenuti fino al secondo livello:

```bash
tree . -L 2
```

??? abstract "Contenuti della directory"

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

Clicca sulla casella colorata per espandere la sezione e visualizzarne i contenuti.
Usiamo sezioni espandibili come questa per mostrare in modo conciso l'output dei comandi previsto e i contenuti di directory e file.

- **I file `.nf`** sono script di workflow numerati in base alla parte del corso in cui vengono utilizzati.

- **Il file `nextflow.config`** è un file di configurazione che imposta proprietà minime dell'ambiente.
  Puoi ignorarlo per ora.

- **Il file `greetings.csv`** sotto `data/` contiene i dati di input che useremo nella maggior parte del corso. Viene descritto nella Parte 2 (Eseguire pipeline), quando lo introduciamo per la prima volta.

- **I file `test-params.*`** sono file di configurazione che useremo nella Parte 3 (Configurazione). Puoi ignorarli per ora.

- **La directory `solutions`** contiene lo stato finale del workflow e dei suoi file accessori (config e moduli) che risultano dal completamento del corso.
  Sono pensati per essere usati come riferimento per controllare il tuo lavoro e risolvere eventuali problemi.

## Checklist di preparazione

Pensi di essere pronto a iniziare?

- [ ] Comprendo l'obiettivo di questo corso e i suoi prerequisiti
- [ ] Il mio ambiente è attivo e funzionante
- [ ] Ho impostato la mia directory di lavoro in modo appropriato

Se puoi spuntare tutte le caselle, sei pronto per partire.

**Per continuare alla [Parte 1: Operazioni base](./01_basics.md), clicca sulla freccia nell'angolo in basso a destra di questa pagina.**
