# Per iniziare

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduzione assistita da IA - [scopri di più e suggerisci miglioramenti](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/PIjOdFaYwWA?si=gZxlXgkVxuLEzOsC" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

/// caption
:fontawesome-brands-youtube:{ .youtube } Guarda [l'intera playlist](https://youtube.com/playlist?list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&si=eF7cLR62goy-lc6n) sul canale YouTube di Nextflow.

:green_book: La trascrizione del video è disponibile [qui](./transcripts/00_orientation.md).
///

!!! tip

    I video di YouTube hanno alcune funzionalità speciali!

    - :fontawesome-solid-closed-captioning: Sottotitoli di alta qualità (curati manualmente). Li attivi con l'icona :material-subtitles:
    - :material-bookmark: Capitoli video nella timeline che corrispondono alle intestazioni della pagina.

## Avviare un ambiente di formazione

Per utilizzare l'ambiente preconfigurato che forniamo su GitHub Codespaces, cliccate sul pulsante "Open in GitHub Codespaces" qui sotto. Per altre opzioni, consultate [Opzioni per l'ambiente](../envsetup/index.md).

Vi consigliamo di aprire l'ambiente di formazione in una nuova scheda o finestra del browser (usate il clic destro, ctrl-clic o cmd-clic a seconda del vostro dispositivo) in modo da poter continuare a leggere mentre l'ambiente si carica.
Dovrete tenere queste istruzioni aperte in parallelo per lavorare attraverso il corso.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

### Nozioni di base sull'ambiente

Questo ambiente di formazione contiene tutto il software, il codice e i dati necessari per lavorare attraverso il corso di formazione, quindi non è necessario installare nulla.

Il codespace è configurato con un'interfaccia VSCode, che include un esplora file, un editor di codice e una shell del terminale.
Tutte le istruzioni fornite durante il corso (ad es. 'aprite il file', 'modificate il codice' o 'eseguite questo comando') si riferiscono a queste tre parti dell'interfaccia VScode, salvo diversa indicazione.

Se state seguendo questo corso autonomamente, familiarizzatevi con le [nozioni di base sull'ambiente](../envsetup/01_setup.md) per ulteriori dettagli.

### Requisiti di versione

Questa formazione è progettata per Nextflow 25.10.2 o successivo **con il parser di sintassi v2 ABILITATO**.
Se state utilizzando un ambiente locale o personalizzato, assicuratevi di utilizzare le impostazioni corrette come documentato [qui](../info/nxf_versions.md).

## Prepararsi a lavorare

Una volta che il codespace è in esecuzione, ci sono due cose da fare prima di immergersi nella formazione: impostare la directory di lavoro per questo corso specifico e dare un'occhiata ai materiali forniti.

### Impostare la directory di lavoro

Per impostazione predefinita, il codespace si apre con la directory di lavoro impostata alla radice di tutti i corsi di formazione, ma per questo corso lavoreremo nella directory `hello-nextflow/`.

Cambiate directory ora eseguendo questo comando nel terminale:

```bash
cd hello-nextflow/
```

Potete impostare VSCode per focalizzarsi su questa directory, in modo che solo i file rilevanti vengano mostrati nella barra laterale dell'esplora file:

```bash
code .
```

!!! tip

    Se per qualsiasi motivo uscite da questa directory (ad es. il codespace va in sleep), potete sempre usare il percorso completo per tornarci, assumendo che stiate eseguendo all'interno dell'ambiente di formazione Github Codespaces:

    ```bash
    cd /workspaces/training/hello-nextflow
    ```

Ora diamo un'occhiata ai contenuti.

### Esplorare i materiali forniti

Potete esplorare i contenuti di questa directory utilizzando l'esplora file sul lato sinistro dell'area di lavoro di formazione.
In alternativa, potete usare il comando `tree`.

Durante il corso, utilizziamo l'output di `tree` per rappresentare la struttura della directory e i contenuti in una forma leggibile, a volte con piccole modifiche per chiarezza.

Qui generiamo una tabella dei contenuti fino al secondo livello:

```bash
tree . -L 2
```

??? abstract "Contenuti della directory"

    ```console
    .
    ├── data
    │   └── greetings.csv
    ├── hello-channels.nf
    ├── hello-config.nf
    ├── hello-containers.nf
    ├── hello-modules.nf
    ├── hello-workflow.nf
    ├── hello-world.nf
    ├── nextflow.config
    ├── solutions
    │   ├── 1-hello-world
    │   ├── 2-hello-channels
    │   ├── 3-hello-workflow
    │   ├── 4-hello-modules
    │   ├── 5-hello-containers
    │   └── 6-hello-config
    ├── test-params.json
    └── test-params.yaml
    ```

Cliccate sul riquadro colorato per espandere la sezione e visualizzarne i contenuti.
Utilizziamo sezioni espandibili come questa per includere l'output dei comandi atteso in modo conciso.

- **I file `.nf`** sono script di workflow denominati in base alla parte del corso in cui vengono utilizzati.

- **Il file `nextflow.config`** è un file di configurazione che imposta proprietà minime dell'ambiente.
  Potete ignorarlo per ora.

- **Il file `greetings.csv`** sotto `data/` contiene i dati di input che useremo nella maggior parte del corso. È descritto nella Parte 2 (Channels), quando lo introduciamo per la prima volta.

- **I file `test-params.*`** sono file di configurazione che useremo nella Parte 6 (Configuration). Potete ignorarli per ora.

- **La directory `solutions`** contiene gli script di workflow completati che risultano da ogni passaggio del corso.
  Sono pensati per essere usati come riferimento per verificare il vostro lavoro e risolvere eventuali problemi.

## Checklist di preparazione

Pensate di essere pronti a iniziare?

- [ ] Comprendo l'obiettivo di questo corso e i suoi prerequisiti
- [ ] Il mio ambiente è attivo e funzionante
- [ ] Ho impostato la mia directory di lavoro in modo appropriato

Se potete spuntare tutte le caselle, siete pronti a partire.

**Per continuare alla [Parte 1: Hello World](./01_hello_world.md), cliccate sulla freccia nell'angolo in basso a destra di questa pagina.**
