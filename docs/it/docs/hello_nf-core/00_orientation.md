# Iniziare

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduzione assistita da IA - [scopri di più e suggerisci miglioramenti](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

## Avviare un ambiente di formazione

Per utilizzare l'ambiente preconfigurato che forniamo su GitHub Codespaces, cliccate sul pulsante "Open in GitHub Codespaces" qui sotto. Per altre opzioni, consultate [Opzioni di ambiente](../envsetup/index.md).

Consigliamo di aprire l'ambiente di formazione in una nuova scheda o finestra del browser (utilizzi il clic destro, ctrl-clic o cmd-clic a seconda della sua dotazione) in modo da poter continuare a leggere mentre l'ambiente si carica.
Dovrà tenere queste istruzioni aperte in parallelo per seguire il corso.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

### Nozioni di base sull'ambiente

Questo ambiente di formazione contiene tutto il software, il codice e i dati necessari per seguire il corso di formazione, quindi non dovrà installare nulla da sola/o.

Il codespace è configurato con un'interfaccia VSCode, che include un esploratore del filesystem, un editor di codice e una shell di terminale.
Tutte le istruzioni fornite durante il corso (ad es. 'aprite il file', 'modificate il codice' o 'eseguite questo comando') si riferiscono a queste tre parti dell'interfaccia VSCode, salvo diversa indicazione.

Se sta seguendo questo corso autonomamente, Vi preghiamo di familiarizzare con le [nozioni di base sull'ambiente](../envsetup/01_setup.md) per ulteriori dettagli.

### Requisiti di versione

Questa formazione è progettata per **Nextflow 25.10.2** o successivo **con il parser di sintassi v2 DISABILITATO**.

#### Se sta utilizzando il nostro ambiente di formazione:

DEVE eseguire il seguente comando prima di procedere oltre:

```bash
export NXF_SYNTAX_PARSER=v1
```

#### Se state utilizzando un ambiente locale o personalizzato:

Assicuratevi di utilizzare le impostazioni corrette come documentato [qui](../info/nxf_versions.md).

La formazione richiede inoltre **nf-core tools 3.4.1**.
Se utilizzate una versione diversa degli strumenti nf-core, potreste avere difficoltà a seguire.

Potete verificare quale versione è installata nel vostro ambiente utilizzando il comando `nf-core --version`.

## Prepararsi a lavorare

Una volta che il vostro codespace è in esecuzione, ci sono due cose da fare prima di immergersi nella formazione: impostare la directory di lavoro per questo corso specifico e dare un'occhiata ai materiali forniti.

### Impostare la directory di lavoro

Per impostazione predefinita, il codespace si apre con la directory di lavoro impostata alla radice di tutti i corsi di formazione, ma per questo corso lavoreremo nella directory `hello-nf-core/`.

Cambiate directory ora eseguendo questo comando nel terminale:

```bash
cd hello-nf-core/
```

!!! tip "Suggerimento"

    Se per qualsiasi motivo vi spostate fuori da questa directory (ad es. il vostro codespace va in sospensione), potete sempre utilizzare il percorso completo per ritornarvi, supponendo che stiate eseguendo questo all'interno dell'ambiente di formazione Github Codespaces:

    ```bash
    cd /workspaces/training/hello-nf-core
    ```

Ora diamo un'occhiata al contenuto di questa directory.

### Esplorare i materiali forniti

Potete esplorare il contenuto di questa directory utilizzando l'esploratore di file sul lato sinistro dello spazio di lavoro di formazione.
In alternativa, potete utilizzare il comando `tree`.

Durante il corso, utilizziamo l'output di `tree` per rappresentare la struttura e il contenuto delle directory in forma leggibile, a volte con lievi modifiche per chiarezza.

Qui generiamo un indice dei contenuti fino al secondo livello:

```bash
tree . -L 2
```

??? abstract "Contenuto della directory"

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

Clicchi sulla casella colorata per espandere la sezione e visualizzarne il contenuto.
Utilizziamo sezioni espandibili come questa per includere l'output previsto dei comandi in modo conciso.

- **Il file `greetings.csv`** è un CSV contenente alcuni dati colonnari minimi che utilizziamo a scopo di test.

- **La directory `original-hello`** contiene una copia del codice sorgente prodotto lavorando attraverso la serie completa di formazione Hello Nextflow (con Docker abilitato).

- **La directory `solutions`** contiene gli script del flusso di lavoro completati che risultano da ogni fase del corso.
  Sono destinati a essere utilizzati come riferimento per verificare il suo lavoro e risolvere eventuali problemi.

## Lista di controllo della preparazione

Pensa di essere pronto/a per iniziare?

- [ ] Comprendo l'obiettivo di questo corso e i suoi prerequisiti
- [ ] Il mio ambiente è attivo e funzionante
- [ ] Mi sono assicurato/a che il parser di sintassi sia impostato su **v1**
- [ ] Ho impostato la mia directory di lavoro in modo appropriato

Se potete spuntare tutte le caselle, siete pronti per iniziare.

**Per continuare con la Parte 1, cliccate sulla freccia nell'angolo in basso a destra di questa pagina.**
