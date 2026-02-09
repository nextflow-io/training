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

Questa formazione è progettata per **Nextflow 25.10.2** o successivo **con il parser di sintassi v2 DISABILITATO**.

#### Se state utilizzando il nostro ambiente di formazione:

DOVETE eseguire il seguente comando prima di procedere:

```bash
export NXF_SYNTAX_PARSER=v1
```

#### Se state utilizzando un ambiente locale o personalizzato:

Assicuratevi di utilizzare le impostazioni corrette come documentato [qui](../info/nxf_versions.md).

La formazione richiede inoltre **nf-core tools 3.4.1**.
Se utilizzate una versione diversa degli strumenti nf-core, potreste avere difficoltà a seguire.

Potete verificare quale versione è installata nel vostro ambiente utilizzando il comando `nf-core --version`.

## Prepararsi a lavorare

Una volta che il vostro codespace è in esecuzione, ci sono due cose che dovete fare prima di immergervi nella formazione: impostare la directory di lavoro per questo corso specifico e dare un'occhiata ai materiali forniti.

### Impostare la directory di lavoro

Per impostazione predefinita, il codespace si apre con la directory di lavoro impostata alla radice di tutti i corsi di formazione, ma per questo corso lavoreremo nella directory `hello-nf-core/`.

Cambiate directory ora eseguendo questo comando nel terminale:

```bash
cd hello-nf-core/
```

!!! tip "Suggerimento"

    Se per qualsiasi motivo vi spostate da questa directory (ad es. il vostro codespace va in sospensione), potete sempre utilizzare il percorso completo per ritornarvi, supponendo che stiate eseguendo questo nell'ambiente di formazione Github Codespaces:

    ```bash
    cd /workspaces/training/hello-nf-core
    ```

Ora diamo un'occhiata ai contenuti di questa directory.

### Esplorare i materiali forniti

Potete esplorare i contenuti di questa directory utilizzando l'esploratore di file sul lato sinistro dell'area di lavoro di formazione.
In alternativa, potete utilizzare il comando `tree`.

Durante il corso, utilizziamo l'output di `tree` per rappresentare la struttura e i contenuti della directory in forma leggibile, a volte con modifiche minori per chiarezza.

Qui generiamo un indice dei contenuti fino al secondo livello:

```bash
tree . -L 2
```

??? abstract "Directory contents"

    ```console
    .
    ├── greetings.csv
    ├── original-hello
    │   ├── hello.nf
    │   ├── modules
    │   └── nextflow.config
    └── solutions
        ├── composable-hello
        ├── core-hello-part2
        ├── core-hello-part3
        ├── core-hello-part4
        ├── core-hello-part5
        └── core-hello-start
    ```

Cliccate sulla casella colorata per espandere la sezione e visualizzarne i contenuti.
Utilizziamo sezioni comprimibili come questa per includere l'output previsto dei comandi in modo conciso.

- **Il file `greetings.csv`** è un CSV contenente alcuni dati colonnari minimi che utilizziamo per scopi di test.

- **La directory `original-hello`** contiene una copia del codice sorgente prodotto lavorando attraverso la serie completa di formazione Hello Nextflow (con Docker abilitato).

- **La directory `solutions`** contiene gli script del flusso di lavoro completati che risultano da ogni passaggio del corso.
  Sono destinati ad essere utilizzati come riferimento per verificare il vostro lavoro e risolvere eventuali problemi.

## Lista di controllo della preparazione

Pensate di essere pronti per iniziare?

- [ ] Comprendo l'obiettivo di questo corso e i suoi prerequisiti
- [ ] Il mio ambiente è attivo e funzionante
- [ ] Mi sono assicurato che il parser di sintassi sia impostato su **v1**
- [ ] Ho impostato la mia directory di lavoro in modo appropriato

Se potete spuntare tutte le caselle, siete pronti per partire.

**Per continuare alla Parte 1, cliccate sulla freccia nell'angolo in basso a destra di questa pagina.**
