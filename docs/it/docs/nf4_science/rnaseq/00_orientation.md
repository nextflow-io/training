# Orientamento

L'ambiente di formazione contiene tutto il software, il codice e i dati necessari per seguire questo corso di formazione, quindi non è necessario installare nulla.
Tuttavia, è necessario un account (gratuito) per effettuare il login, e dovreste prendervi qualche minuto per familiarizzare con l'interfaccia.

Se non l'avete ancora fatto, completate il mini-corso [Configurazione dell'Ambiente](../../envsetup/) prima di procedere.

## Materiali forniti

Durante questo corso di formazione, lavoreremo nella directory `nf4-science/rnaseq/`, nella quale dovrete spostarvi quando aprite lo spazio di lavoro di formazione.
Questa directory contiene tutti i file di codice, i dati di test e i file accessori di cui avrete bisogno.

Sentitevi liberi di esplorare il contenuto di questa directory; il modo più semplice per farlo è utilizzare l'esploratore di file sul lato sinistro dello spazio di lavoro di formazione nell'interfaccia VSCode.
In alternativa, potete utilizzare il comando `tree`.
Durante il corso, utilizziamo l'output di `tree` per rappresentare la struttura e il contenuto delle directory in una forma leggibile, a volte con piccole modifiche per maggiore chiarezza.

Qui generiamo un indice dei contenuti fino al secondo livello:

```bash
tree . -L 3
```

??? success "Directory contents"

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

    Non preoccupatevi se questo sembra molto; esamineremo le parti rilevanti ad ogni passo del corso.
    Questo è solo per darvi una panoramica.

**Ecco un riepilogo di ciò che dovreste sapere per iniziare:**

- **Il file `rnaseq.nf`** è la struttura dello script del flusso di lavoro che svilupperemo.

- **Il file `nextflow.config`** è un file di configurazione che imposta proprietà minime dell'ambiente. Potete ignorarlo per ora.

- **La directory `data`** contiene i dati di input e le risorse correlate:

  - _Un genoma di riferimento_ chiamato `genome.fa` costituito da una piccola regione del cromosoma umano 20 (da hg19/b37).
  - _Dati RNAseq_ che sono stati ridotti a una piccola regione per mantenere le dimensioni dei file contenute, nella directory `reads/`.
  - _File CSV_ che elencano gli ID e i percorsi dei file di dati di esempio, per l'elaborazione in batch.

- **La directory `solutions`** contiene gli script del flusso di lavoro completi e i moduli che risultano da ogni passo del corso.
  Sono destinati ad essere utilizzati come riferimento per verificare il vostro lavoro e risolvere eventuali problemi.
  Il numero nel nome del file corrisponde al passo della parte rilevante del corso.

!!!tip

    Se per qualsiasi motivo vi spostate da questa directory, potete sempre eseguire questo comando per ritornarvi:

    ```bash
    cd /workspaces/training/nf4-science/rnaseq
    ```

Ora, per iniziare il corso, cliccate sulla freccia nell'angolo in basso a destra di questa pagina.
