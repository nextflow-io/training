# Primi passi

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduzione assistita da IA - [scopri di più e suggerisci miglioramenti](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

## Avviare un ambiente di formazione

Per utilizzare l'ambiente preconfigurato che forniamo su GitHub Codespaces, cliccate sul pulsante "Open in GitHub Codespaces" qui sotto. Per altre opzioni, consultate [Opzioni di ambiente](../../envsetup/index.md).

Vi consigliamo di aprire l'ambiente di formazione in una nuova scheda o finestra del browser (utilizzate il clic destro, ctrl-clic o cmd-clic a seconda del vostro dispositivo) in modo da poter continuare a leggere mentre l'ambiente si carica.
Dovrete tenere queste istruzioni aperte in parallelo per seguire il corso.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

### Nozioni di base sull'ambiente

Questo ambiente di formazione contiene tutto il software, il codice e i dati necessari per seguire il corso di formazione, quindi non è necessario installare nulla autonomamente.

Il codespace è configurato con un'interfaccia VSCode, che include un esploratore del filesystem, un editor di codice e una shell del terminale.
Tutte le istruzioni fornite durante il corso (ad es. 'aprite il file', 'modificate il codice' o 'eseguite questo comando') si riferiscono a queste tre parti dell'interfaccia VSCode salvo diversa indicazione.

Se state seguendo questo corso autonomamente, vi preghiamo di familiarizzare con le [nozioni di base sull'ambiente](../../envsetup/01_setup.md) per ulteriori dettagli.

### Requisiti di versione

Questa formazione è progettata per Nextflow 25.10.2 o successivo **con il parser di sintassi v2 ABILITATO**.
Se state utilizzando un ambiente locale o personalizzato, assicuratevi di utilizzare le impostazioni corrette come documentato [qui](../../info/nxf_versions.md).

## Preparatevi a lavorare

Una volta che il vostro codespace è in esecuzione, ci sono due cose che dovete fare prima di immergervi nella formazione: impostare la directory di lavoro per questo corso specifico ed esaminare i materiali forniti.

### Impostare la directory di lavoro

Per impostazione predefinita, il codespace si apre con la directory di lavoro impostata alla radice di tutti i corsi di formazione, ma per questo corso lavoreremo nella directory `nf4-science/rnaseq/`.

Cambiate directory ora eseguendo questo comando nel terminale:

```bash
cd nf4-science/rnaseq/
```

Potete impostare VSCode per concentrarsi su questa directory, in modo che solo i file pertinenti vengano visualizzati nella barra laterale dell'esploratore di file:

```bash
code .
```

!!! tip "Suggerimento"

    Se per qualsiasi motivo vi spostate fuori da questa directory (ad es. il vostro codespace va in sospensione), potete sempre utilizzare il percorso completo per ritornarvi, supponendo che stiate eseguendo questo nell'ambiente di formazione GitHub Codespaces:

    ```bash
    cd /workspaces/training/nf4-science/rnaseq
    ```

Ora diamo un'occhiata ai contenuti.

### Esplorare i materiali forniti

Potete esplorare i contenuti di questa directory utilizzando l'esploratore di file sul lato sinistro dell'area di lavoro di formazione.
In alternativa, potete utilizzare il comando `tree`.

Nel corso della formazione, utilizziamo l'output di `tree` per rappresentare la struttura e i contenuti delle directory in una forma leggibile, talvolta con modifiche minori per chiarezza.

Qui generiamo un sommario dei contenuti fino al terzo livello:

```bash
tree . -L 3
```

??? abstract "Contenuto della directory"

    ```console
    .
    ├── data
    │   ├── genome.fa
    │   ├── paired-end.csv
    │   ├── reads
    │   │   ├── ENCSR000COQ1_1.fastq.gz
    │   │   ├── ENCSR000COQ1_2.fastq.gz
    │   │   ├── ENCSR000COQ2_1.fastq.gz
    │   │   ├── ENCSR000COQ2_2.fastq.gz
    │   │   ├── ENCSR000COR1_1.fastq.gz
    │   │   ├── ENCSR000COR1_2.fastq.gz
    │   │   ├── ENCSR000COR2_1.fastq.gz
    │   │   ├── ENCSR000COR2_2.fastq.gz
    │   │   ├── ENCSR000CPO1_1.fastq.gz
    │   │   ├── ENCSR000CPO1_2.fastq.gz
    │   │   ├── ENCSR000CPO2_1.fastq.gz
    │   │   └── ENCSR000CPO2_2.fastq.gz
    │   └── single-end.csv
    ├── nextflow.config
    ├── rnaseq.nf
    └── solutions
        ├── modules
        │   ├── fastqc.nf
        │   ├── fastqc_pe.nf
        │   ├── hisat2_align.nf
        │   ├── hisat2_align_pe.nf
        │   ├── multiqc.nf
        │   ├── trim_galore.nf
        │   └── trim_galore_pe.nf
        ├── rnaseq-2.1.nf
        ├── rnaseq-2.2.nf
        ├── rnaseq-2.3.nf
        ├── rnaseq-3.1.nf
        ├── rnaseq-3.2.nf
        └── rnaseq_pe-3.3.nf
    ```

Cliccate sulla casella colorata per espandere la sezione e visualizzarne i contenuti.
Utilizziamo sezioni comprimibili come questa per visualizzare l'output previsto dei comandi, nonché i contenuti di directory e file in modo conciso.

- **Il file `rnaseq.nf`** è una struttura per uno script del flusso di lavoro che costruirete man mano che procedete nel corso.

- **La directory `modules`** contiene strutture per i moduli di processo che compilerete durante il corso.

- **Il file `nextflow.config`** è un file di configurazione che imposta proprietà minime dell'ambiente.
  Potete ignorarlo per ora.

- **La directory `data`** contiene i dati di input e le risorse correlate, descritte più avanti nel corso.

- **La directory `solutions`** contiene gli script del flusso di lavoro completi e i moduli risultanti da ogni fase del corso.
  Sono destinati a essere utilizzati come riferimento per verificare il vostro lavoro e risolvere eventuali problemi.
  La soluzione della Parte 2 può essere utilizzata come punto di partenza per la Parte 3.

## Lista di controllo della preparazione

Pensate di essere pronti per iniziare?

- [ ] Comprendo l'obiettivo di questo corso e i suoi prerequisiti
- [ ] Il mio ambiente è attivo e funzionante
- [ ] Ho impostato la mia directory di lavoro in modo appropriato

Se potete spuntare tutte le caselle, siete pronti per iniziare.

**Per continuare alla [Parte 1: Panoramica del metodo](./01_method.md), cliccate sulla freccia nell'angolo in basso a destra di questa pagina.**
