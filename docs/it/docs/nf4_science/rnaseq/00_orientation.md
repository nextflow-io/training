# Orientamento

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduzione assistita da IA - [scopri di più e suggerisci miglioramenti](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

L'ambiente di formazione contiene tutto il software, il codice e i dati necessari per seguire questo corso di formazione, quindi non è necessario installare nulla autonomamente.
Tuttavia, è necessario un account (gratuito) per effettuare l'accesso, e si dovrebbe dedicare qualche minuto per familiarizzare con l'interfaccia.

Se non lo ha ancora fatto, completi il mini-corso [Configurazione dell'Ambiente](../../envsetup/) prima di procedere.

## Materiali forniti

Nel corso di questa formazione, lavoreremo nella directory `nf4-science/rnaseq/`, nella quale è necessario spostarsi all'apertura dell'ambiente di formazione.
Questa directory contiene tutti i file di codice, i dati di test e i file accessori necessari.

Si senta libero di esplorare i contenuti di questa directory; il modo più semplice per farlo è utilizzare l'esploratore di file sul lato sinistro dell'ambiente di formazione nell'interfaccia VSCode.
In alternativa, potete utilizzare il comando `tree`.
Nel corso della formazione, utilizziamo l'output di `tree` per rappresentare la struttura e i contenuti delle directory in una forma leggibile, talvolta con modifiche minori per chiarezza.

Qui generiamo un sommario dei contenuti fino al secondo livello:

```bash
tree . -L 3
```

??? success "Contenuti della directory"

    ```console
    rnaseq
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

!!!note

    Non si preoccupi se questo sembra molto; esamineremo le parti pertinenti in ogni fase del corso.
    Questo è solo per fornire una panoramica generale.

**Ecco un riepilogo di ciò che dovrebbe sapere per iniziare:**

- **Il file `rnaseq.nf`** è la struttura dello script del workflow che svilupperemo.

- **Il file `nextflow.config`** è un file di configurazione che imposta proprietà minime dell'ambiente. Può ignorarlo per ora.

- **La directory `data`** contiene i dati di input e le risorse correlate:

  - _Un genoma di riferimento_ chiamato `genome.fa` costituito da una piccola regione del cromosoma umano 20 (da hg19/b37).
  - _Dati RNAseq_ che sono stati ridotti a una piccola regione per mantenere ridotte le dimensioni dei file, nella directory `reads/`.
  - _File CSV_ che elencano gli ID e i percorsi dei file di dati di esempio, per l'elaborazione in batch.

- **La directory `solutions`** contiene gli script del workflow completi e i moduli risultanti da ogni fase del corso.
  Sono destinati a essere utilizzati come riferimento per verificare il proprio lavoro e risolvere eventuali problemi.
  Il numero nel nome del file corrisponde alla fase della parte pertinente del corso.

!!!tip

    Se per qualsiasi motivo vi spostate fuori da questa directory, potete sempre eseguire questo comando per ritornarvi:

    ```bash
    cd /workspaces/training/nf4-science/rnaseq
    ```

Ora, per iniziare il corso, cliccate sulla freccia nell'angolo in basso a destra di questa pagina.
