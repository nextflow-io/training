# Per iniziare

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduzione assistita da IA - [scopri di più e suggerisci miglioramenti](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>


## Avviare un ambiente di formazione

Per utilizzare l'ambiente preconfigurato che forniamo su GitHub Codespaces, clicca sul pulsante "Open in GitHub Codespaces" qui sotto. Per altre opzioni, consulta [Opzioni di ambiente](../envsetup/index.md).

Consigliamo di aprire l'ambiente di formazione in una nuova scheda o finestra del browser (usa il tasto destro del mouse, ctrl-click o cmd-click a seconda del dispositivo) in modo da poter continuare a leggere mentre l'ambiente si carica.
Sarà necessario tenere queste istruzioni aperte in parallelo per seguire il corso.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

### Nozioni di base sull'ambiente

Questo ambiente di formazione contiene tutto il software, il codice e i dati necessari per seguire il corso, quindi non è necessario installare nulla.

Il codespace è configurato con un'interfaccia VSCode, che include un esplora file, un editor di codice e un terminale.
Tutte le istruzioni fornite durante il corso (ad es. 'apri il file', 'modifica il codice' o 'esegui questo comando') si riferiscono a queste tre parti dell'interfaccia VSCode, salvo diversa indicazione.

Se stai seguendo questo corso da solo, ti invitiamo a familiarizzare con le [nozioni di base sull'ambiente](../envsetup/01_setup.md) per ulteriori dettagli.

### Requisiti di versione

Questo corso richiede Nextflow 25.10.2 o versioni successive, con il parser della sintassi v2 abilitato (impostazione predefinita dalla versione 25.10+).
Se stai utilizzando un ambiente locale o personalizzato, assicurati di utilizzare le impostazioni corrette come documentato [qui](../info/nxf_versions.md).

## Prepararsi al lavoro

Una volta avviato il codespace, ci sono due cose da fare prima di iniziare: impostare la directory di lavoro e dare un'occhiata ai materiali forniti.

### Impostare la directory di lavoro

Per impostazione predefinita, il codespace si apre alla radice di tutti i corsi di formazione.
Per questo corso, spostatevi nella directory `config-exec/`:

```bash
cd config-exec/
```

Poi impostate VSCode per concentrarsi su questa directory, in modo che nella barra laterale dell'esplora file appaiano solo i file pertinenti:

```bash
code .
```

!!! tip "Suggerimento"

    Se per qualsiasi motivo vi spostate da questa directory (ad es. il codespace va in sospensione), potete sempre usare il percorso completo per tornarci, supponendo che stiate eseguendo il corso nell'ambiente di formazione GitHub Codespaces:

    ```bash
    cd /workspaces/training/config-exec
    ```

### Esplorare i materiali forniti

Potete esplorare i materiali del corso usando l'esplora file sulla sinistra, oppure con il comando `tree`.
Eseguite il seguente comando dal terminale per vedere la struttura completa:

```bash
tree . -L 2
```

??? abstract "Contenuto della directory"

    ```console
    .
    ├── data
    │   └── greetings.csv
    ├── main.nf
    ├── modules
    │   ├── collectGreetings.nf
    │   ├── convertToUpper.nf
    │   ├── cowpy.nf
    │   └── sayHello.nf
    └── nextflow.config
    ```

I file **`main.nf`** e **`modules/`** sono la stessa pipeline multi-step di [Nextflow Run](../nextflow_run/index.md), e il file **`nextflow.config`** è la stessa configurazione che avete già visto lì.
Estenderete entrambi nel corso di questi esercizi.

La directory **`data/`** contiene il file di input CSV da cui la pipeline legge i dati.

## Lista di controllo per la preparazione

Pensate di essere pronti a iniziare?

- [ ] Comprendo l'obiettivo di questo corso e i suoi prerequisiti
- [ ] Il mio ambiente è attivo e funzionante
- [ ] Ho impostato correttamente la mia directory di lavoro

Se riuscite a spuntare tutte le caselle, siete pronti per partire.

**Per continuare alla [Parte 1: Adattarsi al proprio ambiente di calcolo](./01_packaging_and_execution.md), cliccate sulla freccia nell'angolo in basso a destra di questa pagina.**
