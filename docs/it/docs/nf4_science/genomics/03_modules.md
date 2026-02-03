# Parte 3: Spostare il codice in moduli

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduzione assistita da IA - [scopri di più e suggerisci miglioramenti](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Nella prima parte di questo corso, avete costruito una pipeline di variant calling completamente lineare che processava i dati di ciascun campione indipendentemente dagli altri.

Nella seconda parte, Le abbiamo mostrato come utilizzare i canali e gli operatori di canale per implementare il joint variant calling con GATK, costruendo sulla pipeline della Parte 1.

In questa parte, Le mostreremo come convertire il codice in quel workflow in moduli. Per seguire questa parte della formazione, dovrebbe aver completato la Parte 1 e la Parte 2, così come [Hello Modules](../../../hello_nextflow/hello_modules.md), che copre le basi dei moduli.

---

## 0. Preparazione

Quando abbiamo iniziato a sviluppare il nostro workflow, abbiamo inserito tutto in un unico file di codice.
Ora è il momento di affrontare la **modularizzazione** del nostro codice, _cioè_ estrarre le definizioni dei processi in moduli.

Inizieremo con lo stesso workflow della Parte 2, che Le abbiamo fornito nel file `genomics-3.nf`.

!!! note "Nota"

     Assicuratevi di essere nella directory di lavoro corretta:
     `cd /workspaces/training/nf4-science/genomics`

Esegua il workflow per verificare il punto di partenza:

```bash
nextflow run genomics-3.nf -resume
```

```console title="Output"
 N E X T F L O W   ~  version 25.10.2

Launching `genomics-3.nf` [serene_borg] DSL2 - revision: 0cbebb67a1

executor >  local (7)
[6f/83ee72] SAMTOOLS_INDEX (3)       | 3 of 3 ✔
[53/b9d342] GATK_HAPLOTYPECALLER (1) | 3 of 3 ✔
[0c/fa6d15] GATK_JOINTGENOTYPING     | 1 of 1 ✔
```

Ora ci sarà una directory `work` e una directory `results_genomics` all'interno della vostra directory di progetto.

### Takeaway

È pronto per iniziare a modularizzare il vostro workflow.

### Prossimi passi

Spostare i processi del workflow Genomics in moduli.

---

## 1. Spostare i processi in moduli

Come avete appreso in [Hello Modules](../../../hello_nextflow/hello_modules.md), potete creare un modulo semplicemente copiando la definizione del processo in un proprio file, in qualsiasi directory, e potete nominare quel file come desiderate.

Per motivi che diventeranno chiari più avanti (in particolare quando arriveremo ai test), in questa formazione seguiremo la convenzione di nominare il file `main.nf`, e posizionarlo in una struttura di directory denominata secondo il toolkit e il comando.

### 1.1. Creare un modulo per il processo `SAMTOOLS_INDEX`

Nel caso del processo `SAMTOOLS_INDEX`, 'samtools' è il toolkit e 'index' è il comando. Quindi, creeremo una struttura di directory `modules/samtools/index` e inseriremo la definizione del processo `SAMTOOLS_INDEX` nel file `main.nf` all'interno di quella directory.

```bash
mkdir -p modules/samtools/index
touch modules/samtools/index/main.nf
```

Aprite il file `main.nf` e copi la definizione del processo `SAMTOOLS_INDEX` al suo interno.

```groovy title="modules/samtools/index/main.nf" linenums="1"
/*
 * Genera il file indice BAM
 */
process SAMTOOLS_INDEX {

    container 'community.wave.seqera.io/library/samtools:1.20--b5dfbd93de237464'

    input:
    path input_bam

    output:
    tuple path(input_bam), path("${input_bam}.bai")

    script:
    """
    samtools index '$input_bam'
    """
}
```

Poi, rimuova la definizione del processo `SAMTOOLS_INDEX` da `genomics-3.nf`, e aggiunga una dichiarazione di importazione per il modulo prima della definizione del processo successivo, in questo modo:

=== "Dopo"

    ```groovy title="genomics-3.nf" linenums="1" hl_lines="1 2"
    // Includi moduli
    include { SAMTOOLS_INDEX } from './modules/samtools/index/main.nf'

    /*
     * Chiama varianti con GATK HaplotypeCaller
     */
    process GATK_HAPLOTYPECALLER {
    ```

=== "Prima"

    ```groovy title="genomics-3.nf" linenums="1" hl_lines="1"
    /*
     * Chiama varianti con GATK HaplotypeCaller
     */
    process GATK_HAPLOTYPECALLER {
    ```

Ora potete eseguire nuovamente il workflow, e dovrebbe continuare a funzionare nello stesso modo di prima. Se fornisce il flag `-resume`, non dovrebbe nemmeno essere necessario eseguire nuove attività:

```bash
nextflow run genomics-3.nf -resume
```

??? success "Output del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `genomics-3.nf` [sleepy_snyder] DSL2 - revision: aa68d06c43

    [0f/71b55e] SAMTOOLS_INDEX (1)       | 3 of 3, cached: 3 ✔
    [f1/18971b] GATK_HAPLOTYPECALLER (3) | 3 of 3, cached: 3 ✔
    [0c/fa6d15] GATK_JOINTGENOTYPING     | 1 of 1, cached: 1 ✔
    ```

### 1.2. Creare moduli per i processi `GATK_HAPLOTYPECALLER` e `GATK_JOINTGENOTYPING`

Ripeta gli stessi passaggi per i processi rimanenti.
Per ciascun processo:

1. Crei la struttura di directory (`modules/gatk/haplotypecaller/` e `modules/gatk/jointgenotyping/`)
2. Crei un file `main.nf` contenente la definizione del processo
3. Rimuova la definizione del processo da `genomics-3.nf`
4. Aggiunga una dichiarazione di importazione per il modulo

Una volta terminato, verifichi che la struttura della directory dei moduli sia corretta eseguendo:

```bash
tree modules/
```

??? abstract "Contenuto della directory"

    ```console
    modules/
    ├── gatk
    │   ├── haplotypecaller
    │   │   └── main.nf
    │   └── jointgenotyping
    │       └── main.nf
    └── samtools
        └── index
            └── main.nf

    5 directories, 3 files
    ```

Dovrebbe anche avere qualcosa di simile nel file del workflow principale, dopo la sezione dei parametri:

```
include { SAMTOOLS_INDEX } from './modules/samtools/index/main.nf'
include { GATK_HAPLOTYPECALLER } from './modules/gatk/haplotypecaller/main.nf'
include { GATK_JOINTGENOTYPING } from './modules/gatk/jointgenotyping/main.nf'

workflow {
```

### Takeaway

Ha praticato la modularizzazione di un workflow, con il workflow genomics come esempio.

### Prossimi passi

Testare il workflow modularizzato.

---

## 2. Testare il workflow modularizzato

Esegua il workflow modularizzato per verificare che tutto funzioni ancora.

```bash
nextflow run genomics-3.nf -resume
```

```console title="Output"
 N E X T F L O W   ~  version 25.10.2

Launching `genomics-3.nf` [astonishing_venter] DSL2 - revision: ca27264c13

[6f/83ee72] SAMTOOLS_INDEX (3)       | 3 of 3, cached: 3 ✔
[53/b9d342] GATK_HAPLOTYPECALLER (3) | 3 of 3, cached: 3 ✔
[0c/fa6d15] GATK_JOINTGENOTYPING     | 1 of 1, cached: 1 ✔
```

Tutto funziona ancora, inclusa la capacità di riprendere la pipeline.
I risultati continuano ad essere pubblicati nella directory `results_genomics`.

```console title="Contenuto della directory"
results_genomics/
├── family_trio.joint.vcf
├── family_trio.joint.vcf.idx
├── gvcf
│   ├── reads_father.bam.g.vcf
│   ├── reads_father.bam.g.vcf.idx
│   ├── reads_mother.bam.g.vcf
│   ├── reads_mother.bam.g.vcf.idx
│   ├── reads_son.bam.g.vcf
│   └── reads_son.bam.g.vcf.idx
└── indexed_bam
    ├── reads_father.bam
    ├── reads_father.bam.bai
    ├── reads_mother.bam
    ├── reads_mother.bam.bai
    ├── reads_son.bam
    └── reads_son.bam.bai
```

### Takeaway

Ha modularizzato un workflow e verificato che funzioni ancora nello stesso modo di prima.

### Prossimi passi

Rivedere ciò che avete appreso e guardare avanti ai test.

---

## 3. Riepilogo

Ha modularizzato il workflow, e nulla è cambiato nel modo in cui la pipeline funziona.
Questo è intenzionale: ha ristrutturato il codice senza impattare la sua funzione.

I moduli contengono solo la logica del processo, rendendoli puliti e riutilizzabili.
Lo script principale controlla cosa viene pubblicato e dove, mentre i moduli rimangono focalizzati sulla loro attività computazionale.

Ha posto le fondamenta per cose che renderanno il vostro codice più facile da mantenere.
Ad esempio, ora potete aggiungere test alla vostra pipeline utilizzando il framework nf-test.
Questo è ciò che esamineremo nella prossima parte di questo corso.
