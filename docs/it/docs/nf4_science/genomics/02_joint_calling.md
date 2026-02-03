# Parte 2: Chiamata congiunta su una coorte

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduzione assistita da IA - [scopri di più e suggerisci miglioramenti](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Nella prima parte di questo corso, avete costruito una pipeline di chiamata delle varianti completamente lineare che processava i dati di ciascun campione indipendentemente dagli altri.
Tuttavia, in un caso d'uso genomico reale, tipicamente sarà necessario esaminare le chiamate di varianti di più campioni insieme.

In questa seconda parte, Le mostriamo come utilizzare i canali e gli operatori di canale per implementare la chiamata congiunta delle varianti con GATK, basandosi sulla pipeline della Parte 1.

### Panoramica del metodo

Il metodo di chiamata delle varianti GATK che abbiamo utilizzato nella prima parte di questo corso generava semplicemente chiamate di varianti per campione.
Questo va bene se si desidera solo esaminare le varianti di ciascun campione in isolamento, ma ciò fornisce informazioni limitate.
Spesso è più interessante esaminare come le chiamate di varianti differiscono tra più campioni, e per fare ciò, GATK offre un metodo alternativo chiamato chiamata congiunta delle varianti, che dimostriamo qui.

La chiamata congiunta delle varianti comporta la generazione di un tipo speciale di output di varianti chiamato GVCF (per Genomic VCF) per ciascun campione, quindi la combinazione dei dati GVCF di tutti i campioni e infine, l'esecuzione di un'analisi statistica di 'genotipizzazione congiunta'.

![Analisi congiunta](img/joint-calling.png)

Ciò che è speciale in un GVCF di un campione è che contiene record che riassumono le statistiche dei dati di sequenza su tutte le posizioni nell'area mirata del genoma, non solo le posizioni in cui il programma ha trovato evidenza di variazione.
Questo è fondamentale per il calcolo della genotipizzazione congiunta ([ulteriori informazioni](https://gatk.broadinstitute.org/hc/en-us/articles/360035890431-The-logic-of-joint-calling-for-germline-short-variants)).

Il GVCF è prodotto da GATK HaplotypeCaller, lo stesso strumento che abbiamo utilizzato nella Parte 1, con un parametro aggiuntivo (`-ERC GVCF`).
La combinazione dei GVCF viene eseguita con GATK GenomicsDBImport, che combina le chiamate per campione in un archivio dati (analogo a un database), quindi l'analisi vera e propria di 'genotipizzazione congiunta' viene eseguita con GATK GenotypeGVCFs.

### Workflow

Quindi per ricapitolare, in questa parte del corso, svilupperemo un workflow che fa quanto segue:

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/genomics/img/hello-gatk-2.svg"
</figure>

1. Generare un file di indice per ciascun file BAM di input utilizzando Samtools
2. Eseguire GATK HaplotypeCaller su ciascun file BAM di input per generare un GVCF delle chiamate di varianti genomiche per campione
3. Raccogliere tutti i GVCF e combinarli in un archivio dati GenomicsDB
4. Eseguire la genotipizzazione congiunta sull'archivio dati GVCF combinato per produrre un VCF a livello di coorte

Applicheremo questo allo stesso dataset della Parte 1.

---

## 0. Riscaldamento: Eseguire Samtools e GATK direttamente

Proprio come in precedenza, vogliamo provare i comandi manualmente prima di tentare di incorporarli in un workflow.

!!! note

     Assicuratevi di essere nella directory di lavoro corretta:
     `cd /workspaces/training/nf4-science/genomics`

### 0.1. Indicizzare un file BAM di input con Samtools

Questo primo passo è lo stesso della Parte 1, quindi dovrebbe risultare molto familiare, ma questa volta dobbiamo farlo per tutti e tre i campioni.

!!! note

    Tecnicamente abbiamo già generato file di indice per i tre campioni attraverso la nostra pipeline, quindi potremmo andare a recuperarli dalla directory dei risultati. Tuttavia, è più pulito semplicemente rifarlo manualmente, e ci vorrà solo un minuto.

#### 0.1.1. Avviare il container Samtools in modo interattivo

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/samtools:1.20--b5dfbd93de237464
```

<!--
??? success "Output del comando"

    ```console

    ```
-->

#### 0.1.2. Eseguire il comando di indicizzazione per i tre campioni

```bash
samtools index /data/bam/reads_mother.bam
samtools index /data/bam/reads_father.bam
samtools index /data/bam/reads_son.bam
```

Proprio come in precedenza, questo dovrebbe produrre i file di indice nella stessa directory dei file BAM corrispondenti.

??? abstract "Contenuto della directory"

    ```console
    data/bam/
    ├── reads_father.bam
    ├── reads_father.bam.bai
    ├── reads_mother.bam
    ├── reads_mother.bam.bai
    ├── reads_son.bam
    └── reads_son.bam.bai
    ```

Ora che abbiamo file di indice per tutti e tre i campioni, possiamo procedere alla generazione dei GVCF per ciascuno di essi.

#### 0.1.3. Uscire dal container Samtools

```bash
exit
```

### 0.2. Chiamare le varianti con GATK HaplotypeCaller in modalità GVCF

Questo secondo passo è molto simile a ciò che abbiamo fatto nella Parte 1: Hello Genomics, ma ora eseguiremo GATK in 'modalità GVCF'.

#### 0.2.1. Avviare il container GATK in modo interattivo

```bash
docker run -it -v ./data:/data community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867
```

<!--
??? success "Output del comando"

    ```console

    ```
-->

#### 0.2.2. Eseguire il comando di chiamata delle varianti con l'opzione GVCF

Per produrre un VCF genomico (GVCF), aggiungiamo l'opzione `-ERC GVCF` al comando di base, che attiva la modalità GVCF di HaplotypeCaller.

Modifichiamo anche l'estensione del file per il file di output da `.vcf` a `.g.vcf`.
Tecnicamente questo non è un requisito, ma è una convenzione fortemente raccomandata.

```bash
gatk HaplotypeCaller \
        -R /data/ref/ref.fasta \
        -I /data/bam/reads_mother.bam \
        -O reads_mother.g.vcf \
        -L /data/ref/intervals.bed \
        -ERC GVCF
```

<!--
??? success "Output del comando"

    ```console

    ```
-->

Questo crea il file di output GVCF `reads_mother.g.vcf` nella directory di lavoro corrente nel container.

Se lo visualizzate con `cat` per vederne il contenuto, vedrà che è molto più lungo del VCF equivalente che abbiamo generato nella Parte 1. Non potete nemmeno scorrere fino all'inizio del file, e la maggior parte delle righe appare abbastanza diversa da ciò che abbiamo visto nel VCF della Parte 1.

```console title="Output" linenums="1674"
20_10037292_10066351    14714   .       T       <NON_REF>       .       .       END=14718       GT:DP:GQ:MIN_DP:PL       0/0:37:99:37:0,99,1192
20_10037292_10066351    14719   .       T       <NON_REF>       .       .       END=14719       GT:DP:GQ:MIN_DP:PL       0/0:36:82:36:0,82,1087
20_10037292_10066351    14720   .       T       <NON_REF>       .       .       END=14737       GT:DP:GQ:MIN_DP:PL       0/0:42:99:37:0,100,1160
```

Questi rappresentano regioni non varianti dove il chiamante di varianti non ha trovato evidenza di variazione, quindi ha catturato alcune statistiche che descrivono il suo livello di confidenza nell'assenza di variazione. Questo rende possibile distinguere tra due casi molto diversi: (1) ci sono dati di buona qualità che mostrano che il campione è omozigote-riferimento, e (2) non ci sono abbastanza dati di buona qualità disponibili per fare una determinazione in un modo o nell'altro.

In un GVCF, ci sono tipicamente molte di queste righe non varianti, con un numero minore di record di varianti sparse tra di esse. Provi a eseguire `head -176` sul GVCF per caricare solo le prime 176 righe del file per trovare una chiamata di variante effettiva.

```console title="Output" linenums="174"
20_10037292_10066351    3479    .       T       <NON_REF>       .       .       END=3479        GT:DP:GQ:MIN_DP:PL       0/0:34:36:34:0,36,906
20_10037292_10066351    3480    .       C       CT,<NON_REF>    503.03  .       DP=23;ExcessHet=0.0000;MLEAC=2,0;MLEAF=1.00,0.00;RAW_MQandDP=82800,23    GT:AD:DP:GQ:PL:SB       1/1:0,18,0:18:54:517,54,0,517,54,517:0,0,7,11
20_10037292_10066351    3481    .       T       <NON_REF>       .       .       END=3481        GT:DP:GQ:MIN_DP:PL       0/0:21:51:21:0,51,765
```

La seconda riga mostra il primo record di variante nel file, che corrisponde alla prima variante nel file VCF che abbiamo esaminato nella Parte 1.

Proprio come il VCF originale, anche il file GVCF di output è accompagnato da un file di indice, chiamato `reads_mother.g.vcf.idx`.

#### 0.2.3. Ripetere il processo sugli altri due campioni

Per testare il passo di genotipizzazione congiunta, abbiamo bisogno di GVCF per tutti e tre i campioni, quindi generiamoli manualmente ora.

```bash
gatk HaplotypeCaller \
        -R /data/ref/ref.fasta \
        -I /data/bam/reads_father.bam \
        -O reads_father.g.vcf \
        -L /data/ref/intervals.bed \
        -ERC GVCF
```

<!--
??? success "Output del comando"

    ```console

    ```
-->

```bash
gatk HaplotypeCaller \
        -R /data/ref/ref.fasta \
        -I /data/bam/reads_son.bam \
        -O reads_son.g.vcf \
        -L /data/ref/intervals.bed \
        -ERC GVCF
```

<!--
??? success "Output del comando"

    ```console

    ```
-->

Una volta completato, dovrebbe avere tre file che terminano con `.g.vcf` nella vostra directory corrente (uno per campione) e i rispettivi file di indice che terminano con `.g.vcf.idx`.

### 0.3. Eseguire la genotipizzazione congiunta

Ora che abbiamo tutti i GVCF, possiamo finalmente provare l'approccio di genotipizzazione congiunta per generare chiamate di varianti per una coorte di campioni.
Come promemoria, è un metodo in due fasi che consiste nel combinare i dati di tutti i GVCF in un archivio dati, quindi eseguire l'analisi di genotipizzazione congiunta vera e propria per generare il VCF finale delle varianti chiamate congiuntamente.

#### 0.3.1. Combinare tutti i GVCF per campione

Questa prima fase utilizza un altro strumento GATK, chiamato GenomicsDBImport, per combinare i dati di tutti i GVCF in un archivio dati GenomicsDB.

```bash
gatk GenomicsDBImport \
    -V reads_mother.g.vcf \
    -V reads_father.g.vcf \
    -V reads_son.g.vcf \
    -L /data/ref/intervals.bed \
    --genomicsdb-workspace-path family_trio_gdb
```

<!--
??? success "Output del comando"

    ```console

    ```
-->

L'output di questa fase è effettivamente una directory contenente un insieme di ulteriori directory nidificate che contengono i dati di varianti combinati sotto forma di più file diversi.
Può esplorarlo ma vedrà rapidamente che questo formato di archivio dati non è destinato a essere letto direttamente dagli esseri umani.

!!! note

    GATK include strumenti che rendono possibile ispezionare ed estrarre dati di chiamate di varianti dall'archivio dati secondo necessità.

#### 0.3.2. Eseguire l'analisi di genotipizzazione congiunta vera e propria

Questa seconda fase utilizza un altro strumento GATK, chiamato GenotypeGVCFs, per ricalcolare le statistiche delle varianti e i genotipi individuali alla luce dei dati disponibili in tutti i campioni della coorte.

```bash
gatk GenotypeGVCFs \
    -R /data/ref/ref.fasta \
    -V gendb://family_trio_gdb \
    -O family_trio.vcf
```

<!--
??? success "Output del comando"

    ```console

    ```
-->

Questo crea il file di output VCF `family_trio.vcf` nella directory di lavoro corrente nel container.
È un altro file ragionevolmente piccolo quindi potete visualizzare questo file con `cat` per vederne il contenuto, e scorrere verso l'alto per trovare le prime righe di varianti.

```console title="family_trio.vcf" linenums="40"
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  reads_father    reads_mother    reads_son
20_10037292_10066351    3480    .       C       CT      1625.89 .       AC=5;AF=0.833;AN=6;BaseQRankSum=0.220;DP=85;ExcessHet=0.0000;FS=2.476;MLEAC=5;MLEAF=0.833;MQ=60.00;MQRankSum=0.00;QD=21.68;ReadPosRankSum=-1.147e+00;SOR=0.487    GT:AD:DP:GQ:PL  0/1:15,16:31:99:367,0,375       1/1:0,18:18:54:517,54,0 1/1:0,26:26:78:756,78,0
20_10037292_10066351    3520    .       AT      A       1678.89 .       AC=5;AF=0.833;AN=6;BaseQRankSum=1.03;DP=80;ExcessHet=0.0000;FS=2.290;MLEAC=5;MLEAF=0.833;MQ=60.00;MQRankSum=0.00;QD=22.39;ReadPosRankSum=0.701;SOR=0.730 GT:AD:DP:GQ:PL   0/1:18,13:31:99:296,0,424       1/1:0,18:18:54:623,54,0 1/1:0,26:26:78:774,78,0
20_10037292_10066351    3529    .       T       A       154.29  .       AC=1;AF=0.167;AN=6;BaseQRankSum=-5.440e-01;DP=104;ExcessHet=0.0000;FS=1.871;MLEAC=1;MLEAF=0.167;MQ=60.00;MQRankSum=0.00;QD=7.71;ReadPosRankSum=-1.158e+00;SOR=1.034       GT:AD:DP:GQ:PL  0/0:44,0:44:99:0,112,1347       0/1:12,8:20:99:163,0,328        0/0:39,0:39:99:0,105,1194
```

Questo assomiglia di più al VCF originale che abbiamo generato nella Parte 1, tranne che questa volta abbiamo informazioni a livello di genotipo per tutti e tre i campioni.
Le ultime tre colonne nel file sono i blocchi di genotipo per i campioni, elencati in ordine alfabetico.

Se guardiamo ai genotipi chiamati per il nostro trio familiare di test per la primissima variante, vediamo che il padre è eterozigote-variante (`0/1`), e la madre e il figlio sono entrambi omozigoti-variante (`1/1`).

Questa è in definitiva l'informazione che stiamo cercando di estrarre dal dataset! Quindi andiamo a incorporare tutto questo in un workflow Nextflow in modo da poterlo fare su larga scala.

#### 0.3.3. Uscire dal container GATK

```bash
exit
```

### Takeaway

Sa come eseguire i singoli comandi coinvolti nella chiamata congiunta delle varianti nel terminale per verificare che producano le informazioni desiderate.

### Quali sono i prossimi passi?

Incorporare questi comandi in una pipeline vera e propria.

---

## 1. Modificare il passo di chiamata delle varianti per campione per produrre un GVCF

La buona notizia è che non dobbiamo ricominciare da zero, poiché abbiamo già scritto un workflow che fa parte di questo lavoro nella Parte 1.
Tuttavia, quella pipeline produce file VCF, mentre ora vogliamo file GVCF per fare la genotipizzazione congiunta.
Quindi dobbiamo iniziare attivando la modalità di chiamata delle varianti GVCF e aggiornando l'estensione del file di output.

!!! note

    Per convenienza, lavoreremo con una copia fresca del workflow GATK come si presenta alla fine della Parte 1, ma con un nome diverso: `genomics-2.nf`.

### 1.1. Dire a HaplotypeCaller di emettere un GVCF e aggiornare l'estensione di output

Apriamo il file `genomics-2.nf` nell'editor di codice.
Dovrebbe apparire molto familiare, ma si senta libero di eseguirlo se vuole accertarsi che funzioni come previsto.

Inizieremo facendo due modifiche:

- Aggiungere il parametro `-ERC GVCF` al comando GATK HaplotypeCaller;
- Aggiornare il percorso del file di output per utilizzare l'estensione corrispondente `.g.vcf`, secondo la convenzione GATK.

Assicuratevi di aggiungere una barra rovesciata (`\`) alla fine della riga precedente quando aggiunge `-ERC GVCF`.

=== "Dopo"

    ```groovy title="genomics-2.nf" linenums="56" hl_lines="4 6"
        """
        gatk HaplotypeCaller \
            -R ${ref_fasta} \
            -I ${input_bam} \
            -O ${input_bam}.g.vcf \
            -L ${interval_list} \
            -ERC GVCF
        """
    ```

=== "Prima"

    ```groovy title="genomics-2.nf" linenums="56" hl_lines="4"
        """
        gatk HaplotypeCaller \
            -R ${ref_fasta} \
            -I ${input_bam} \
            -O ${input_bam}.vcf \
            -L ${interval_list}
        """
    ```

E questo è tutto ciò che serve per passare HaplotypeCaller alla generazione di GVCF invece di VCF, giusto?

### 1.2. Eseguire la pipeline per verificare che si possano generare GVCF

Il comando di esecuzione Nextflow è lo stesso di prima, salvo per il nome del file del workflow stesso.
Assicuratevi di aggiornarlo appropriatamente.

```bash
nextflow run genomics-2.nf
```

??? failure "Output del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics-2.nf` [crazy_venter] DSL2 - revision: a2d6f6f09f

    executor >  local (6)
    [f1/8d8486] SAMTOOLS_INDEX (1)       | 3 of 3 ✔
    [72/3249ca] GATK_HAPLOTYPECALLER (3) | 0 of 3
    ERROR ~ Error executing process > 'GATK_HAPLOTYPECALLER (2)'

    Caused by:
      Missing output file(s) `reads_son.bam.vcf` expected by process `GATK_HAPLOTYPECALLER (2)`

    Command executed:

      gatk HaplotypeCaller         -R ref.fasta         -I reads_son.bam         -O reads_son.bam.g.vcf         -L intervals.bed         -ERC GVCF
    ```

E l'output è... tutto rosso! Oh no.

Il comando che è stato eseguito è corretto, quindi avevamo ragione sul fatto che quello fosse sufficiente per cambiare il comportamento dello strumento GATK.
Ma guardi quella riga sul file di output mancante. Nota qualcosa?

Esatto, abbiamo dimenticato di dire a Nextflow di aspettarsi un nuovo nome di file. Ops.

### 1.3. Aggiornare l'estensione del file di output anche nel blocco degli output del processo

Perché non è sufficiente cambiare solo l'estensione del file nel comando dello strumento stesso, bisogna anche dire a Nextflow che il nome del file di output previsto è cambiato.

=== "Dopo"

    ```groovy title="genomics-2.nf" linenums="50" hl_lines="2 3"
        output:
        path "${input_bam}.g.vcf"     , emit: vcf
        path "${input_bam}.g.vcf.idx" , emit: idx
    ```

=== "Prima"

    ```groovy title="genomics-2.nf" linenums="50" hl_lines="2 3"
        output:
        path "${input_bam}.vcf"     , emit: vcf
        path "${input_bam}.vcf.idx" , emit: idx
    ```

### 1.4. Aggiornare i target di pubblicazione per i nuovi output GVCF

Poiché ora stiamo producendo GVCF invece di VCF, dovremmo aggiornare la sezione `publish:` del workflow per utilizzare nomi più descrittivi.
Organizzeremo anche i file GVCF nella loro sottodirectory per chiarezza.

=== "Dopo"

    ```groovy title="genomics-2.nf" linenums="88" hl_lines="3 4"
        publish:
        indexed_bam = SAMTOOLS_INDEX.out
        gvcf = GATK_HAPLOTYPECALLER.out.vcf
        gvcf_idx = GATK_HAPLOTYPECALLER.out.idx
    ```

=== "Prima"

    ```groovy title="genomics-2.nf" linenums="88"
        publish:
        indexed_bam = SAMTOOLS_INDEX.out
        vcf = GATK_HAPLOTYPECALLER.out.vcf
        vcf_idx = GATK_HAPLOTYPECALLER.out.idx
    ```

### 1.5. Aggiornare il blocco output per la nuova struttura di directory

Dobbiamo anche aggiornare il blocco `output` per mettere i file GVCF in una sottodirectory `gvcf`.

=== "Dopo"

    ```groovy title="genomics-2.nf" linenums="94" hl_lines="3 5 6 8 9"
    output {
        indexed_bam {
            path 'indexed_bam'
        }
        gvcf {
            path 'gvcf'
        }
        gvcf_idx {
            path 'gvcf'
        }
    }
    ```

=== "Prima"

    ```groovy title="genomics-2.nf" linenums="94"
    output {
        indexed_bam {
            path '.'
        }
        vcf {
            path '.'
        }
        vcf_idx {
            path '.'
        }
    }
    ```

### 1.6. Eseguire nuovamente la pipeline

Eseguiamola con `-resume` questa volta.

```bash
nextflow run genomics-2.nf -resume
```

??? success "Output del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics-2.nf` [nostalgic_franklin] DSL2 - revision: f2c0a93c6a

    executor >  local (3)
    [cc/fbc705] SAMTOOLS_INDEX (3)       | 3 of 3, cached: 3 ✔
    [27/0d7eb9] GATK_HAPLOTYPECALLER (2) | 3 of 3 ✔
    ```

Questa volta funziona.

L'output di Nextflow stesso non appare diverso (rispetto a un'esecuzione riuscita in modalità VCF normale), ma ora possiamo trovare i file `.g.vcf` e i rispettivi file di indice, per tutti e tre i campioni, organizzati in sottodirectory.

??? abstract "Contenuto della directory (symlink abbreviati)"

    ```console
    results_genomics/
    ├── gvcf/
    │   ├── reads_father.bam.g.vcf -> */27/0d7eb9*/reads_father.bam.g.vcf
    │   ├── reads_father.bam.g.vcf.idx -> */27/0d7eb9*/reads_father.bam.g.vcf.idx
    │   ├── reads_mother.bam.g.vcf -> */e4/4ed55e*/reads_mother.bam.g.vcf
    │   ├── reads_mother.bam.g.vcf.idx -> */e4/4ed55e*/reads_mother.bam.g.vcf.idx
    │   ├── reads_son.bam.g.vcf -> */08/e95962*/reads_son.bam.g.vcf
    │   └── reads_son.bam.g.vcf.idx -> */08/e95962*/reads_son.bam.g.vcf.idx
    └── indexed_bam/
        ├── reads_father.bam -> */9a/c7a873*/reads_father.bam
        ├── reads_father.bam.bai -> */9a/c7a873*/reads_father.bam.bai
        ├── reads_mother.bam -> */f1/8d8486*/reads_mother.bam
        ├── reads_mother.bam.bai -> */f1/8d8486*/reads_mother.bam.bai
        ├── reads_son.bam -> */cc/fbc705*/reads_son.bam
        └── reads_son.bam.bai -> */cc/fbc705*/reads_son.bam.bai
    ```

Se apre uno dei file GVCF e lo scorre, potete verificare che GATK HaplotypeCaller abbia prodotto file GVCF come richiesto.

### Takeaway

Ok, questo è stato minimo in termini di apprendimento di Nextflow...
Ma è stata una bella opportunità per ribadire l'importanza del blocco di output del processo!

### Quali sono i prossimi passi?

Imparare a raccogliere il contenuto di un canale e trasmetterlo al processo successivo come singolo input.

---

## 2. Raccogliere e combinare i dati GVCF su tutti i campioni

Ora dobbiamo combinare i dati di tutti i GVCF per campione in una forma che supporti l'analisi di genotipizzazione congiunta che vogliamo fare.

### 2.1. Definire il processo che combinerà i GVCF

Come promemoria di ciò che abbiamo fatto in precedenza nella sezione di riscaldamento, combinare i GVCF è un compito per lo strumento GATK GenomicsDBImport, che produrrà un archivio dati nel cosiddetto formato GenomicsDB.

Scriviamo un nuovo processo per definire come funzionerà, basato sul comando che abbiamo usato in precedenza nella sezione di riscaldamento.

```groovy title="genomics-2.nf" linenums="66"
/*
 * Combinare i GVCF in un archivio dati GenomicsDB
 */
process GATK_GENOMICSDB {

    container "community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867"

    input:
    path all_gvcfs
    path all_idxs
    path interval_list
    val cohort_name

    output:
    path "${cohort_name}_gdb"

    script:
    """
    gatk GenomicsDBImport \
        -V ${all_gvcfs} \
        -L ${interval_list} \
        --genomicsdb-workspace-path ${cohort_name}_gdb
    """
}
```

Cosa ne pensa, sembra ragionevole?

Colleghiamolo e vediamo cosa succede.

### 2.2. Aggiungere un parametro `cohort_name` con un valore predefinito

Dobbiamo fornire un nome arbitrario per la coorte.
Più avanti nella serie di formazione imparerà come utilizzare i metadati dei campioni per questo tipo di cosa, ma per ora dichiariamo semplicemente un parametro CLI usando `params` e diamogli un valore predefinito per convenienza.

```groovy title="genomics-2.nf" linenums="16"
    // Nome base per il file di output finale
    cohort_name: String = "family_trio"
```

### 2.3. Raccogliere gli output di GATK_HAPLOTYPECALLER attraverso i campioni

Se dovessimo semplicemente collegare il canale di output dal processo `GATK_HAPLOTYPECALLER` così com'è, Nextflow chiamerebbe il processo su ciascun GVCF del campione separatamente.
Tuttavia, vogliamo raggruppare tutti e tre i GVCF (e i loro file di indice) in modo tale che Nextflow li consegni tutti insieme a una singola chiamata di processo.

Buone notizie: possiamo farlo usando l'operatore di canale `collect()`. Aggiungiamo le seguenti righe al corpo del `workflow`, subito dopo la chiamata a GATK_HAPLOTYPECALLER:

```groovy title="genomics-2.nf" linenums="118"
// Raccogliere gli output di chiamata delle varianti attraverso i campioni
all_gvcfs_ch = GATK_HAPLOTYPECALLER.out.vcf.collect()
all_idxs_ch = GATK_HAPLOTYPECALLER.out.idx.collect()
```

Sembra un po' complicato? Scomponiamolo e traduciamolo in linguaggio semplice.

1. Stiamo prendendo il canale di output dal processo `GATK_HAPLOTYPECALLER`, a cui si fa riferimento usando la proprietà `.out`.
2. Ogni 'elemento' che esce dal canale è una coppia di file: il GVCF e il suo file di indice, in quest'ordine perché è l'ordine in cui sono elencati nel blocco di output del processo. Convenientemente, poiché nell'ultima sessione abbiamo nominato gli output di questo processo (usando `emit:`), possiamo selezionare i GVCF da un lato aggiungendo `.vcf` e i file di indice dall'altro aggiungendo `.idx` dopo la proprietà `.out`. Se non avessimo nominato quegli output, avremmo dovuto fare riferimento ad essi rispettivamente con `.out[0]` e `.out[1]`.
3. Aggiungiamo l'operatore di canale `collect()` per raggruppare tutti i file GVCF insieme in un singolo elemento in un nuovo canale chiamato `all_gvcfs_ch`, e facciamo lo stesso con i file di indice per formare il nuovo canale chiamato `all_idxs_ch`.

!!! tip

    Se avete difficoltà a visualizzare esattamente cosa sta succedendo qui, ricordate che potete utilizzare l'operatore `view()` per ispezionare il contenuto dei canali prima e dopo l'applicazione degli operatori di canale.

I canali risultanti `all_gvcfs_ch` e `all_idxs_ch` sono ciò che collegheremo al processo `GATK_GENOMICSDB` che abbiamo appena scritto.

!!! note

    Nel caso si stesse chiedendo, raccogliamo i GVCF e i loro file di indice separatamente perché il comando GATK GenomicsDBImport vuole vedere solo i percorsi dei file GVCF. Fortunatamente, poiché Nextflow organizzerà tutti i file insieme per l'esecuzione, non dobbiamo preoccuparci dell'ordine dei file come abbiamo fatto per i BAM e il loro indice nella Parte 1.

### 2.4. Aggiungere una chiamata al blocco workflow per eseguire GATK_GENOMICSDB

Abbiamo un processo e abbiamo canali di input. Dobbiamo solo aggiungere la chiamata al processo.

```groovy title="genomics-2.nf" linenums="122"
    // Combinare i GVCF in un archivio dati GenomicsDB
    GATK_GENOMICSDB(
        all_gvcfs_ch,
        all_idxs_ch,
        intervals_file,
        params.cohort_name
    )
```

Ok, tutto è collegato.

### 2.5. Eseguire il workflow

Vediamo se funziona.

```bash
nextflow run genomics-2.nf -resume
```

??? failure "Output del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics-2.nf` [disturbed_bell] DSL2 - revision: 57942246cc

    executor >  local (1)
    [f1/8d8486] SAMTOOLS_INDEX (1)       | 3 of 3, cached: 3 ✔
    [e4/4ed55e] GATK_HAPLOTYPECALLER (2) | 3 of 3, cached: 3 ✔
    [51/d350ea] GATK_GENOMICSDB          | 0 of 1
    ERROR ~ Error executing process > 'GATK_GENOMICSDB'

    Caused by:
      Process `GATK_GENOMICSDB` terminated with an error exit status (1)

    Command executed:

      gatk GenomicsDBImport         -V reads_son.bam.g.vcf reads_father.bam.g.vcf reads_mother.bam.g.vcf         -L intervals.bed         --genomicsdb-workspace-path family_trio_gdb
    ```

Viene eseguito abbastanza rapidamente, poiché stiamo eseguendo con `-resume`, ma fallisce!

Ah. Dal lato positivo, vediamo che Nextflow ha rilevato il processo `GATK_GENOMICSDB`, e in particolare lo ha chiamato solo una volta.
Ciò suggerisce che l'approccio `collect()` ha funzionato, fino a un certo punto.
Ma, ed è un grande ma, la chiamata al processo è fallita.

Quando approfondiamo l'output della console sopra, possiamo vedere che il comando eseguito non è corretto.

Riesce a individuare l'errore?
Guardi questo bit: `-V reads_father.bam.g.vcf reads_son.bam.g.vcf reads_mother.bam.g.vcf`

Abbiamo dato a `gatk GenomicsDBImport` più file GVCF per un singolo argomento `-V`, ma lo strumento si aspetta un argomento `-V` separato per ciascun file GVCF.

Come promemoria, questo era il comando che abbiamo eseguito nel container:

```bash
gatk GenomicsDBImport \
    -V reads_mother.g.vcf \
    -V reads_father.g.vcf \
    -V reads_son.g.vcf \
    -L /data/ref/intervals.bed \
    --genomicsdb-workspace-path family_trio_gdb
```

Quindi significa che dobbiamo in qualche modo trasformare il nostro pacchetto di file GVCF in una stringa di comando formattata correttamente.

### 2.6. Costruire una riga di comando con un argomento `-V` separato per ciascun GVCF di input

È qui che Nextflow essendo basato su Groovy torna utile, perché ci permetterà di utilizzare alcune manipolazioni di stringhe abbastanza semplici per costruire la stringa di comando necessaria.

Specificamente, usando questa sintassi: `all_gvcfs.collect { gvcf -> "-V ${gvcf}" }.join(' ')`

Ancora una volta, scomponiamolo nei suoi componenti.

1. Prima, prendiamo il contenuto del canale di input `all_gvcfs` e applichiamo `.collect()` su di esso (proprio come prima).
2. Ciò ci consente di passare ciascun percorso di file GVCF individuale nel pacchetto alla **closure**, `{ gvcf -> "-V ${gvcf}" }`, dove `gvcf` si riferisce a quel percorso di file GVCF.
   La closure è una mini-funzione che usiamo per anteporre `-V ` al percorso del file, nella forma di `"-V ${gvcf}"`.
3. Quindi usiamo `.join(' ')` per concatenare tutte e tre le stringhe con un singolo spazio come separatore.

Con un esempio concreto, appare così:

1. Abbiamo tre file:

   `[A.ext, B.ext, C.ext]`

2. La closure modifica ciascuno per creare le stringhe:

   `"-V A.ext", "-V B.ext", "-V C.ext"`

3. L'operazione `.join(' ')` genera la stringa finale:

   `"-V A.ext -V B.ext -V C.ext"`

Una volta che abbiamo quella stringa, possiamo assegnarla a una variabile locale, `gvcfs_line`, definita con la parola chiave `def`:

`def gvcfs_line = all_gvcfs.collect { gvcf -> "-V ${gvcf}" }.join(' ')`

Ok, quindi abbiamo la nostra cosa di manipolazione delle stringhe. Dove la mettiamo?

Vogliamo che questo vada da qualche parte all'interno della definizione del processo, perché vogliamo farlo _dopo_ aver canalizzato i percorsi dei file GVCF nel processo.
Questo perché Nextflow deve vederli come percorsi di file per organizzare i file stessi correttamente per l'esecuzione.

Ma _dove_ nel processo possiamo aggiungere questo?

Fatto divertente: potete aggiungere codice arbitrario dopo `script:` e prima del `"""` !

Ottimo, aggiungiamo la nostra riga di manipolazione delle stringhe lì allora, e aggiorniamo il comando `gatk GenomicsDBImport` per utilizzare la stringa concatenata che produce.

=== "Dopo"

    ```groovy title="genomics-2.nf" linenums="87"  hl_lines="2 5"
        script:
        def gvcfs_line = all_gvcfs.collect { gvcf -> "-V ${gvcf}" }.join(' ')
        """
        gatk GenomicsDBImport \
            ${gvcfs_line} \
            -L ${interval_list} \
            --genomicsdb-workspace-path ${cohort_name}_gdb
        """
    ```

=== "Prima"

    ```groovy title="genomics-2.nf" linenums="87"  hl_lines="4"
        script:
        """
        gatk GenomicsDBImport \
            -V ${all_gvcfs} \
            -L ${interval_list} \
            --genomicsdb-workspace-path ${cohort_name}_gdb
        """
    ```

Questo dovrebbe essere tutto ciò che è necessario per fornire correttamente gli input a `gatk GenomicsDBImport`.

!!! tip

    Quando aggiorna il comando `gatk GenomicsDBImport`, assicuratevi di rimuovere il prefisso `-V ` quando sostituisce con la variabile `${gvcfs_line}`.

### 2.7. Eseguire il workflow per verificare che generi l'output GenomicsDB come previsto

Va bene, vediamo se questo ha risolto il problema.

```bash
nextflow run genomics-2.nf -resume
```

??? success "Output del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics-2.nf` [peaceful_gates] DSL2 - revision: ca0bf847ed

    executor >  local (1)
    [cc/fbc705] SAMTOOLS_INDEX (3)       | 3 of 3, cached: 3 ✔
    [27/0d7eb9] GATK_HAPLOTYPECALLER (2) | 3 of 3, cached: 3 ✔
    [76/d13861] GATK_GENOMICSDB          | 1 of 1 ✔
    ```

Aha! Sembra funzionare ora.

I primi due passi sono stati saltati con successo, e il terzo passo ha funzionato perfettamente questa volta.
L'archivio dati GenomicsDB viene creato nella directory di lavoro ma non pubblicato nei risultati, poiché è solo un formato intermedio che useremo per la genotipizzazione congiunta.

Tra l'altro, non abbiamo dovuto fare nulla di speciale per gestire l'output che è una directory invece di un singolo file.

### Takeaway

Ora sa come raccogliere output da un canale e raggrupparli come singolo input a un altro processo.
Sa anche come costruire una riga di comando per fornire input a un determinato strumento con la sintassi appropriata.

### Quali sono i prossimi passi?

Imparare come aggiungere un secondo comando allo stesso processo.

---

## 3. Eseguire il passo di genotipizzazione congiunta come parte dello stesso processo

Ora che abbiamo le chiamate di varianti genomiche combinate, possiamo eseguire lo strumento di genotipizzazione congiunta, che produrrà l'output finale che ci interessa veramente: il VCF delle chiamate di varianti a livello di coorte.

Per ragioni logistiche, decidiamo di includere la genotipizzazione congiunta all'interno dello stesso processo.

### 3.1. Rinominare il processo da GATK_GENOMICSDB a GATK_JOINTGENOTYPING

Poiché il processo eseguirà più di uno strumento, cambiamo il suo nome per fare riferimento all'operazione complessiva piuttosto che a un singolo nome di strumento.

=== "Dopo"

    ```groovy title="genomics-2.nf"
    /*
     * Combinare i GVCF in un archivio dati GenomicsDB ed eseguire la genotipizzazione congiunta per produrre chiamate a livello di coorte
     */
    process GATK_JOINTGENOTYPING {
    ```

=== "Prima"

    ```groovy title="genomics-2.nf"
    /*
     * Combinare i GVCF in un archivio dati GenomicsDB
     */
    process GATK_GENOMICSDB {
    ```

Ricordi di mantenere i nomi dei vostri processi il più descrittivi possibile, per massimizzare la leggibilità per i vostri colleghi — e il vostro futuro sé!

### 3.2. Aggiungere il comando di genotipizzazione congiunta al processo GATK_JOINTGENOTYPING

Semplicemente aggiunga il secondo comando dopo il primo all'interno della sezione script.

=== "Dopo"

    ```groovy title="genomics-2.nf" linenums="89"  hl_lines="6-10"
        """
        gatk GenomicsDBImport \
            ${gvcfs_line} \
            -L ${interval_list} \
            --genomicsdb-workspace-path ${cohort_name}_gdb

        gatk GenotypeGVCFs \
            -R ${ref_fasta} \
            -V gendb://${cohort_name}_gdb \
            -L ${interval_list} \
            -O ${cohort_name}.joint.vcf
        """
    ```

=== "Prima"

    ```groovy title="genomics-2.nf" linenums="89"
        """
        gatk GenomicsDBImport \
            ${gvcfs_line} \
            -L ${interval_list} \
            --genomicsdb-workspace-path ${cohort_name}_gdb
        """
    ```

I due comandi verranno eseguiti in serie, nello stesso modo in cui lo farebbero se dovessimo eseguirli manualmente nel terminale.

### 3.3. Aggiungere i file del genoma di riferimento alle definizioni di input del processo GATK_JOINTGENOTYPING

Il secondo comando richiede i file del genoma di riferimento, quindi dobbiamo aggiungerli agli input del processo.

=== "Dopo"

    ```groovy title="genomics-2.nf" linenums="78"  hl_lines="6-8"
        input:
        path all_gvcfs
        path all_idxs
        path interval_list
        val cohort_name
        path ref_fasta
        path ref_index
        path ref_dict
    ```

=== "Prima"

    ```groovy title="genomics-2.nf" linenums="78"
        input:
        path all_gvcfs
        path all_idxs
        path interval_list
        val cohort_name
    ```

Può sembrare fastidioso digitarli, ma ricordi, li digita una volta, e poi potete eseguire il workflow un milione di volte. Ne vale la pena?

### 3.4. Aggiornare la definizione di output del processo per emettere il VCF delle chiamate di varianti a livello di coorte

Non ci interessa veramente salvare l'archivio dati GenomicsDB, che è solo un formato intermedio che esiste solo per ragioni logistiche, quindi possiamo semplicemente rimuoverlo dal blocco di output se vogliamo.

L'output a cui siamo effettivamente interessati è il VCF prodotto dal comando di genotipizzazione congiunta.

=== "Dopo"

    ```groovy title="genomics-2.nf" linenums="87" hl_lines="2 3"
        output:
        path "${cohort_name}.joint.vcf"     , emit: vcf
        path "${cohort_name}.joint.vcf.idx" , emit: idx
    ```

=== "Prima"

    ```groovy title="genomics-2.nf" linenums="87"
        output:
        path "${cohort_name}_gdb"
    ```

Abbiamo quasi finito!

### 3.5. Aggiornare la chiamata al processo da GATK_GENOMICSDB a GATK_JOINTGENOTYPING

Non dimentichiamo di rinominare la chiamata al processo nel corpo del workflow da GATK_GENOMICSDB a GATK_JOINTGENOTYPING. E mentre ci siamo, dovremmo anche aggiungere i file del genoma di riferimento come input, poiché dobbiamo fornirli allo strumento di genotipizzazione congiunta.

=== "Dopo"

    ```groovy title="genomics-2.nf" linenums="126"
    // Combinare i GVCF in un archivio dati GenomicsDB e applicare la genotipizzazione congiunta
    GATK_JOINTGENOTYPING(
        all_gvcfs_ch,
        all_idxs_ch,
        intervals_file,
        params.cohort_name,
        ref_file,
        ref_index_file,
        ref_dict_file
    )
    ```

=== "Prima"

    ```groovy title="genomics-2.nf" linenums="126"
    // Combinare i GVCF in un archivio dati GenomicsDB
    GATK_GENOMICSDB(
        all_gvcfs_ch,
        all_idxs_ch,
        intervals_file,
        params.cohort_name
    )
    ```

Ora il processo è completamente collegato.

### 3.6. Aggiungere il VCF congiunto alla sezione publish

Dobbiamo pubblicare gli output del VCF congiunto dal nuovo processo.
Aggiunga queste righe alla sezione `publish:` del workflow:

```groovy title="genomics-2.nf" linenums="145"
    joint_vcf = GATK_JOINTGENOTYPING.out.vcf
    joint_vcf_idx = GATK_JOINTGENOTYPING.out.idx
```

### 3.7. Aggiungere i target del VCF congiunto al blocco output

Infine, aggiunga target di output per i file VCF congiunti.
Li metteremo alla radice della directory dei risultati poiché questo è l'output finale.

```groovy title="genomics-2.nf" linenums="157"
    joint_vcf {
        path '.'
    }
    joint_vcf_idx {
        path '.'
    }
```

Ora tutto dovrebbe essere completamente collegato.

### 3.8. Eseguire il workflow

Finalmente, possiamo eseguire il workflow modificato...

```bash
nextflow run genomics-2.nf -resume
```

??? success "Output del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics-2.nf` [crazy_marconi] DSL2 - revision: 5da9afc841

    executor >  local (1)
    [9a/c7a873] SAMTOOLS_INDEX (2)       | 3 of 3, cached: 3 ✔
    [e4/4ed55e] GATK_HAPLOTYPECALLER (2) | 3 of 3, cached: 3 ✔
    [a6/7cc8ed] GATK_JOINTGENOTYPING     | 1 of 1 ✔
    ```

E funziona!

Troverà il file di output finale, `family_trio.joint.vcf` (e il suo indice di file), nella directory dei risultati.

??? abstract "Contenuto della directory (symlink abbreviati)"

    ```console
    results_genomics/
    ├── family_trio.joint.vcf -> */a6/7cc8ed*/family_trio.joint.vcf
    ├── family_trio.joint.vcf.idx -> */a6/7cc8ed*/family_trio.joint.vcf.idx
    ├── gvcf/
    │   ├── reads_father.bam.g.vcf -> */27/0d7eb9*/reads_father.bam.g.vcf
    │   ├── reads_father.bam.g.vcf.idx -> */27/0d7eb9*/reads_father.bam.g.vcf.idx
    │   ├── reads_mother.bam.g.vcf -> */e4/4ed55e*/reads_mother.bam.g.vcf
    │   ├── reads_mother.bam.g.vcf.idx -> */e4/4ed55e*/reads_mother.bam.g.vcf.idx
    │   ├── reads_son.bam.g.vcf -> */08/e95962*/reads_son.bam.g.vcf
    │   └── reads_son.bam.g.vcf.idx -> */08/e95962*/reads_son.bam.g.vcf.idx
    └── indexed_bam/
        ├── reads_father.bam -> */9a/c7a873*/reads_father.bam
        ├── reads_father.bam.bai -> */9a/c7a873*/reads_father.bam.bai
        ├── reads_mother.bam -> */f1/8d8486*/reads_mother.bam
        ├── reads_mother.bam.bai -> */f1/8d8486*/reads_mother.bam.bai
        ├── reads_son.bam -> */cc/fbc705*/reads_son.bam
        └── reads_son.bam.bai -> */cc/fbc705*/reads_son.bam.bai
    ```

Se è del tipo scettico, potete fare clic sul file VCF congiunto per aprirlo e verificare che il workflow abbia generato le stesse chiamate di varianti che avete ottenuto eseguendo gli strumenti manualmente all'inizio di questa sezione.

```console title="family_trio.joint.vcf" linenums="40"
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	reads_father	reads_mother	reads_son
20_10037292_10066351	3480	.	C	CT	1625.89	.	AC=5;AF=0.833;AN=6;BaseQRankSum=0.220;DP=85;ExcessHet=0.0000;FS=2.476;MLEAC=5;MLEAF=0.833;MQ=60.00;MQRankSum=0.00;QD=21.68;ReadPosRankSum=-1.147e+00;SOR=0.487	GT:AD:DP:GQ:PL	0/1:15,16:31:99:367,0,375	1/1:0,18:18:54:517,54,0	1/1:0,26:26:78:756,78,0
20_10037292_10066351	3520	.	AT	A	1678.89	.	AC=5;AF=0.833;AN=6;BaseQRankSum=1.03;DP=80;ExcessHet=0.0000;FS=2.290;MLEAC=5;MLEAF=0.833;MQ=60.00;MQRankSum=0.00;QD=22.39;ReadPosRankSum=0.701;SOR=0.730	GT:AD:DP:GQ:PL	0/1:18,13:31:99:296,0,424	1/1:0,18:18:54:623,54,0	1/1:0,26:26:78:774,78,0
20_10037292_10066351	3529	.	T	A	154.29	.	AC=1;AF=0.167;AN=6;BaseQRankSum=-5.440e-01;DP=104;ExcessHet=0.0000;FS=1.871;MLEAC=1;MLEAF=0.167;MQ=60.00;MQRankSum=0.00;QD=7.71;ReadPosRankSum=-1.158e+00;SOR=1.034	GT:AD:DP:GQ:PL	0/0:44,0:44:99:0,112,1347	0/1:12,8:20:99:163,0,328	0/0:39,0:39:99:0,105,1194
```

Ora ha un workflow di chiamata congiunta delle varianti automatizzato, completamente riproducibile!

!!! note

    Tenga presente che i file di dati che Le abbiamo fornito coprono solo una piccola porzione del cromosoma 20.
    La dimensione reale di un set di chiamate di varianti sarebbe contata in milioni di varianti.
    Ecco perché usiamo solo piccoli sottoinsiemi di dati per scopi formativi!

### Takeaway

Sa come utilizzare alcuni operatori comuni così come le closure di Groovy per controllare il flusso di dati nel vostro workflow.

### Quali sono i prossimi passi?

Celebrate il vostro successo e prendetevi una meritata pausa.

Nella prossima parte di questo corso, imparerà come modularizzare il vostro workflow estraendo le definizioni dei processi in moduli riutilizzabili.
