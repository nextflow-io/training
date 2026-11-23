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

Questo corso richiede Nextflow 25.10.2 o successivo, con il parser della sintassi v2 abilitato (impostazione predefinita dalla versione 25.10+).
Se stai usando un ambiente locale o personalizzato, assicurati di utilizzare le impostazioni corrette come documentato [qui](../info/nxf_versions.md).

## Preparati al lavoro

Una volta che il tuo codespace è in esecuzione, ci sono due cose da fare prima di immergerti nella formazione: impostare la directory di lavoro e dare un'occhiata ai materiali forniti.

### Imposta la directory di lavoro

Per impostazione predefinita, il codespace si apre alla radice di tutti i corsi di formazione.
Per questo corso, cambia nella directory `nextflow-run/`:

```bash
cd nextflow-run/
```

Poi imposta VSCode per concentrarsi su questa directory, in modo che solo i file rilevanti vengano mostrati nella barra laterale dell'esploratore file:

```bash
code .
```

!!! tip "Suggerimento"

    Se per qualsiasi motivo esci da questa directory (es. il tuo codespace va in sospensione), puoi sempre usare il percorso completo per tornare, assumendo che tu stia eseguendo questo all'interno dell'ambiente di formazione di Github Codespaces:

    ```bash
    cd /workspaces/training/nextflow-run
    ```

### Esplora i materiali forniti

Puoi esplorare i materiali del corso usando l'esploratore file sul lato sinistro, oppure con il comando `tree`.
Esegui il seguente comando dal terminale per vedere la struttura completa:

```bash
tree . -L 2
```

??? abstract "Contenuto della directory"

    ```console
    .
    ├── 1-hello.nf
    ├── 2-inputs.nf
    ├── data
    │   ├── greetings-extended.csv
    │   └── greetings.csv
    ├── main.nf
    ├── modules
    │   ├── collectGreetings.nf
    │   ├── convertToUpper.nf
    │   ├── cowpy.nf
    │   └── sayHello.nf
    └── nextflow.config
    ```

I **file `.nf`** sono script di workflow di complessità crescente, utilizzati in quest'ordine durante il corso.

La **directory `data/`** contiene i file CSV di input che useremo a partire dalla sezione 2.

La **directory `modules/`** contiene le definizioni dei processi utilizzate da `main.nf`.

Il **file `nextflow.config`** è un file di configurazione che imposta proprietà minime dell'ambiente. Puoi ignorarlo per ora; lo esamineremo nella sezione 4.

## Checklist di preparazione

Pensi di essere pronto a iniziare?

- [ ] Comprendo l'obiettivo di questo corso e i suoi prerequisiti
- [ ] Il mio ambiente è attivo e funzionante
- [ ] Ho impostato la mia directory di lavoro in modo appropriato

Se puoi spuntare tutte le caselle, sei pronto per partire.

**Per continuare alla [Parte 1: Eseguire Nextflow](./01_run_nextflow.md), clicca sulla freccia nell'angolo in basso a destra di questa pagina.**
