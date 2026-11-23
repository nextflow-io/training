# Parte 3: Eseguire una pipeline di produzione

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduzione assistita da IA - [scopri di più e suggerisci miglioramenti](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Nella [Parte 2](./02_configure_execution.md), hai imparato a impostare i parametri e personalizzare la configurazione per nf-core/demo.
Ora applichiamo ciò che abbiamo imparato a una vera pipeline di produzione, nf-core/rnaseq.

---

## 1. Scaricare ed eseguire nf-core/rnaseq

Finora abbiamo utilizzato `nf-core/demo`, che è una pipeline minimale progettata per la formazione.
Ora scarichiamo una vera pipeline di produzione e la eseguiamo con il suo profilo di test.

La pipeline `nf-core/rnaseq` esegue i passaggi principali dell'analisi bulk RNA sequencing: controllo qualità, trimming degli adattatori, allineamento delle reads e quantificazione a livello genico.
È probabilmente la pipeline nf-core più utilizzata fino ad oggi.

### 1.1. Scaricare la pipeline

Eseguiamo il seguente comando per scaricarla.

```bash
nextflow pull nf-core/rnaseq
```

??? success "Output del comando"

    ```console
    Checking nf-core/rnaseq ...
     downloaded from https://github.com/nf-core/rnaseq.git - revision: 1f03b53ef7 [master]
    ```

La pipeline è ora memorizzata nella cache locale ed è pronta per essere eseguita.

### 1.2. Eseguire il profilo di test

Eseguiamola con il profilo di test e Docker:

```bash
nextflow run nf-core/rnaseq -profile test,docker --outdir rnaseq-results
```

??? failure "Output del comando"

    ```console
     N E X T F L O W   ~  version 26.04.4

    Launching `https://github.com/nf-core/rnaseq` [suspicious_dijkstra] revision: 1f03b53ef7 [master]

    ------------------------------------------------------
                                            ,--./,-.
            ___     __   __   __   ___     /,-._.--~'
      |\ | |__  __ /  ` /  \ |__) |__         }  {
      | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                            `._,._,'
      nf-core/rnaseq 3.26.0
    ------------------------------------------------------
    ...
    [-        ] NFCORE_RNASEQ:RNASEQ:FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS:FQ_LINT  -
    Plus 47 more processes waiting for tasks…

    Execution cancelled -- Finishing pending tasks before exit
    -[nf-core/rnaseq] Pipeline completed with errors-
    ERROR ~ Error executing process > 'NFCORE_RNASEQ:RNASEQ:FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS:FQ_LINT (WT_REP2)'

    Caused by:
      Process requirement exceeds available memory -- req: 12 GB; avail: 7.7 GB


    Command executed:

      fq lint \
          --disable-validator P001 \
          SRR6357072_1.fastq.gz SRR6357072_2.fastq.gz > WT_REP2.fq_lint.txt

    Command exit status:
      -

    Command output:
      (empty)

    Work dir:
      /workspaces/training/nfcore-use/work/xx/xxxxxxxxxxxxxxxxxxxxxx

    Container:
      quay.io/biocontainers/fq:0.12.0--h9ee0642_0

    Tip: you can replicate the issue by changing to the process work dir and entering the command `bash .command.run`

     -- Check '.nextflow.log' file for details
    ERROR ~ Pipeline failed. Please refer to troubleshooting docs for common issues: https://nf-co.re/docs/running/troubleshooting
    ```

La riga chiave in quell'errore è:

```console
Process requirement exceeds available memory -- req: 12 GB; avail: 7.7 GB
```

La macchina Codespaces predefinita ha 8 GB di RAM, che è anche il valore predefinito tipico per Docker Desktop.
La pipeline sta richiedendo 12 GB per il processo `FQ_LINT` — più di quanto la macchina possa fornire.

Quei 12 GB provengono dall'etichetta di risorse `process_low` definita in `conf/base.config`:

```groovy title="conf/base.config"
withLabel:process_low {
    cpus   = { 2     * task.attempt }
    memory = { 12.GB * task.attempt }
    time   = { 4.h   * task.attempt }
}
```

Un'opzione sarebbe quella di utilizzare un tipo di macchina più grande, ma a scopo di test vogliamo poter eseguire su qualsiasi hardware disponibile.
L'approccio migliore è sovrascrivere le risorse predefinite in un file di configurazione personalizzato.

### 1.3. Rieseguire con una configurazione personalizzata

Vi forniamo un file di configurazione personalizzato che sovrascrive le risorse predefinite basate sulle etichette.

??? full-code "laptop.config"

    ```groovy title="laptop.config"
    process {
        withLabel: 'process_low' {
            cpus   = 2
            memory = 6.GB
        }
        withLabel: 'process_medium' {
            cpus   = 4
            memory = 6.GB
        }
        withLabel: 'process_high' {
            cpus   = 6
            memory = 6.GB
        }
        withLabel: 'process_high_memory' {
            memory = 6.GB
        }
    }
    ```

La [Parte 2](./02_configure_execution.md) ha introdotto `withName:` per selezionare un singolo processo per nome.
Qui usiamo `withLabel:` per selezionare tutti i processi che condividono un'etichetta contemporaneamente.

Questo file è già presente nella directory di lavoro.
Passiamolo con `-c` per applicare le sovrascritture:

```bash
nextflow run nf-core/rnaseq -profile test,docker -c laptop.config --outdir rnaseq-results
```

??? success "Output del comando (pipeline in avvio)"

    ```console
     N E X T F L O W   ~  version 26.04.4

    Launching `https://github.com/nf-core/rnaseq` [romantic_faraday] revision: 1f03b53ef7 [master]

    ------------------------------------------------------
                                            ,--./,-.
            ___     __   __   __   ___     /,-._.--~'
      |\ | |__  __ /  ` /  \ |__) |__         }  {
      | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                            `._,._,'
      nf-core/rnaseq 3.26.0
    ------------------------------------------------------
    ...
    executor >  local (7)
    [xx/xxxxxx] NFCORE_RNASEQ:RNASEQ:FQ_LINT (RAP1_IAA_30M_REP1)   | 3 of 5
    [xx/xxxxxx] NFCORE_RNASEQ:RNASEQ:FASTQC (RAP1_IAA_30M_REP1)    | 2 of 5
    ...
    ```

La pipeline è ora in esecuzione e potete osservare le attività completarsi una per una.
Su questo dataset di test minimale si completerà in 15–20 minuti, eseguendo oltre 200 attività in totale.

I veri esperimenti di RNA-seq tipicamente coinvolgono decine di campioni e richiedono ore o giorni di esecuzione.
Nextflow supporta gli scheduler HPC (SLURM, PBS, LSF) e le piattaforme cloud (AWS, Google Cloud, Azure), che possono ridurre drasticamente il tempo di esecuzione distribuendo il lavoro su molti nodi.
La configurazione di questi ambienti, tuttavia, aggiunge una complessità significativa.

La piattaforma Seqera (sviluppata dai creatori di Nextflow) fornisce un'interfaccia web per lanciare pipeline Nextflow su infrastrutture HPC o cloud (proprie o gestite per voi), con capacità di gestione del calcolo e dei dati che semplificano il processo di esecuzione delle pipeline su larga scala.

!!! tip "Suggerimento"

    I ricercatori accademici possono accedere a Seqera Platform gratuitamente tramite il [programma accademico Seqera](https://seqera.io/academic-program/).

### Takeaway

Avete scaricato `nf-core/rnaseq`, visto come funzionano le etichette di risorse nf-core e imparato a sovrascriverle con un file di configurazione personalizzato.
Ancora più importante, avete visto perché l'esecuzione locale è un punto di partenza piuttosto che una destinazione per analisi su scala reale.

### Cosa c'è dopo?

Avete coperto i fondamentali dell'esecuzione delle pipeline nf-core.
Consultate i [Passi successivi](next_steps.md) per sapere dove andare da qui.

---

## Riepilogo

In questa parte avete imparato a:

- Scaricare ed eseguire una pipeline su scala di produzione (nf-core/rnaseq) e sovrascrivere le sue etichette di risorse predefinite
