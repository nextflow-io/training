# Per iniziare

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduzione assistita da IA - [scopri di più e suggerisci miglioramenti](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

## Avviare un ambiente di formazione

Per utilizzare l'ambiente preconfigurato che forniamo su GitHub Codespaces, clicca sul pulsante "Open in GitHub Codespaces" qui sotto. Per altre opzioni, consulta [Opzioni di ambiente](../envsetup/index.md).

Consigliamo di aprire l'ambiente di formazione in una nuova scheda o finestra del browser (usa il tasto destro del mouse, ctrl-click o cmd-click a seconda del dispositivo) in modo da poter continuare a leggere mentre l'ambiente si carica.
Sarà necessario tenere queste istruzioni aperte in parallelo per seguire il corso.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

### Nozioni di base sull'ambiente

Questo ambiente di formazione contiene tutto il software, il codice e i dati necessari per seguire il corso, quindi non è necessario installare nulla.

Il codespace è configurato con un'interfaccia VSCode, che include un file explorer, un editor di codice e un terminale.
Tutte le istruzioni fornite durante il corso (ad es. "apri il file", "modifica il codice" o "esegui questo comando") si riferiscono a queste tre parti dell'interfaccia VSCode, salvo diversa indicazione.

Se stai seguendo questo corso da solo, ti invitiamo a familiarizzare con le [nozioni di base sull'ambiente](../envsetup/01_setup.md) per ulteriori dettagli.

### Requisiti di versione

Questa formazione funziona con Nextflow 25.10.2 o versioni successive **con il parser della sintassi v2**, che è quello predefinito a partire da Nextflow 26.04.
Nel nostro ambiente di formazione non è necessario fare nulla: viene eseguito Nextflow 26.04.4 con il parser v2. Se stai utilizzando un ambiente locale o personalizzato, consulta le [note sulla versione](../info/nxf_versions.md).

!!! warning "nf-core/demo richiede Nextflow 25.10.4 o versioni successive"

    La pipeline `nf-core/demo` utilizzata nella Parte 1 impone la propria versione minima di Nextflow (`>=25.10.4`), che è più restrittiva rispetto al requisito minimo generale della formazione di 25.10.2.
    Il nostro ambiente di formazione soddisfa già questo requisito; se stai utilizzando un ambiente locale o personalizzato, assicurati di avere Nextflow 25.10.4 o versioni successive.

Questa formazione richiede inoltre **nf-core tools 4.0.2**.
Se utilizzi una versione diversa degli strumenti nf-core, potresti avere difficoltà a seguire il corso.

Puoi verificare quale versione è installata nel tuo ambiente con il comando `nf-core --version`.

!!! warning "Compatibilità con il parser v2"

    Molte pipeline nf-core non supportano ancora il parser della sintassi v2.
    Se esegui una pipeline nf-core diversa da quelle utilizzate in questo corso e riscontri errori, potrebbe essere necessario passare al parser v1 impostando `export NXF_SYNTAX_PARSER=v1`.
    Consulta le [note sulla versione](../info/nxf_versions.md) per i dettagli.

## Prepararsi al lavoro

Una volta avviato il codespace, ci sono due cose da fare prima di immergersi nella formazione: impostare la directory di lavoro per questo corso specifico e dare un'occhiata ai materiali forniti.

### Impostare la directory di lavoro

Per impostazione predefinita, il codespace si apre con la directory di lavoro impostata alla radice di tutti i corsi di formazione, ma per questo corso lavoreremo nella directory `nfcore-use/`.

Cambia directory ora eseguendo questo comando nel terminale:

```bash
cd nfcore-use/
```

!!! tip "Suggerimento"

    Se per qualsiasi motivo esci da questa directory (ad es. il codespace va in sospensione), puoi sempre usare il percorso completo per tornarci, supponendo che tu stia lavorando nell'ambiente di formazione su Github Codespaces:

    ```bash
    cd /workspaces/training/nfcore-use
    ```

Successivamente, esplora il contenuto di questa directory.

### Esplorare i materiali forniti

Puoi esplorare il contenuto di questa directory utilizzando il file explorer sul lato sinistro dell'area di lavoro.
In alternativa, puoi usare il comando `tree`.

```bash
tree .
```

??? abstract "Contenuto della directory"

    ```console
    .
    ├── custom.config
    ├── laptop.config
    ├── malformed_samplesheet.csv
    └── my_params.yml
    ```

- **Il file `laptop.config`** è un file di configurazione che utilizzeremo nella sezione 4 per limitare l'utilizzo delle risorse quando si esegue una pipeline su scala produttiva in locale.
  Puoi ignorarlo fino ad allora.
- **I file `my_params.yml`, `malformed_samplesheet.csv` e `custom.config`** sono utilizzati nella Parte 2, per illustrare come impostare i parametri da un file, la validazione dell'input e le configurazioni a livello di processo.
  Puoi ignorarli fino ad allora.

## Lista di controllo per la preparazione

Pensi di essere pronto a iniziare?

- [ ] Comprendo l'obiettivo di questo corso e i suoi prerequisiti
- [ ] Il mio ambiente è attivo e funzionante
- [ ] Sto usando nf-core tools 4.0.2 (verifica con `nf-core --version`)
- [ ] Ho impostato correttamente la mia directory di lavoro

Se puoi spuntare tutte le caselle, sei pronto per partire.

**Per continuare alla Parte 1, clicca sulla freccia nell'angolo in basso a destra di questa pagina.**
