# Parte 3: Joint calling su una coorte

Nella Parte 2 avete costruito una pipeline di variant calling per campione che processava i dati di ciascun campione in modo indipendente.
Ora la estenderemo per implementare il joint variant calling, come trattato nella [Parte 1](01_method.md).

## Compito

In questa parte del corso estenderemo il flusso di lavoro per fare quanto segue:

<figure class="excalidraw">
--8<-- "docs/en/docs/nf4_science/genomics/img/hello-gatk-2.svg"
</figure>

1. Generare un file indice per ciascun file BAM di input utilizzando Samtools
2. Eseguire GATK HaplotypeCaller su ciascun file BAM di input per generare un GVCF di chiamate di varianti genomiche per campione
3. Raccogliere tutti i GVCF e combinarli in un data store GenomicsDB
4. Eseguire il joint genotyping sui dati GVCF combinati per produrre un VCF a livello di coorte

Questa parte si basa direttamente sul flusso di lavoro prodotto nella Parte 2.

??? info "Come iniziare da questa sezione"

    Questa sezione del corso presuppone che abbiate completato la [Parte 2: Variant calling per campione](./02_per_sample_variant_calling.md) e abbiate una pipeline `genomics.nf` funzionante.

    Se non avete completato la Parte 2 o volete ripartire da zero per questa parte, potete utilizzare la soluzione della Parte 2 come punto di partenza.
    Eseguite questi comandi dall'interno della directory `nf4-science/genomics/`:

    ```bash
    cp solutions/part2/genomics-2.nf genomics.nf
    cp solutions/part2/nextflow.config .
    cp solutions/part2/modules/* modules/
    ```

    Questo vi fornisce un flusso di lavoro completo di variant calling per campione.
    Potete testare che funzioni correttamente eseguendo il seguente comando:

    ```bash
    nextflow run genomics.nf -profile test
    ```

## Piano della lezione

Abbiamo suddiviso questo in due passaggi:

1. **Modificare lo step di variant calling per campione per produrre un GVCF.**
   Questo copre l'aggiornamento dei comandi e degli output del processo.
2. **Aggiungere uno step di joint genotyping che combini e genotipi i GVCF per campione.**
   Questo introduce l'operatore `collect()`, le closure Groovy per la costruzione della riga di comando e i processi multi-comando.

!!! note

     Assicuratevi di essere nella directory di lavoro corretta:
     `cd /workspaces/training/nf4-science/genomics`

---

## 1. Modificare lo step di variant calling per campione per produrre un GVCF

La pipeline della Parte 2 produce file VCF, ma il joint calling richiede file GVCF.
Dobbiamo attivare la modalità di variant calling GVCF e aggiornare l'estensione del file di output.

Richiamiamo il comando di variant calling GVCF dalla [Parte 1](01_method.md):

```bash
gatk HaplotypeCaller \
        -R /data/ref/ref.fasta \
        -I /data/bam/reads_mother.bam \
        -O reads_mother.g.vcf \
        -L /data/ref/intervals.bed \
        -ERC GVCF
```

Rispetto al comando base HaplotypeCaller che abbiamo incapsulato nella Parte 2, le differenze sono il parametro `-ERC GVCF` e l'estensione di output `.g.vcf`.

### 1.1. Dire a HaplotypeCaller di emettere un GVCF e aggiornare l'estensione di output

Aprite il file del modulo `modules/gatk_haplotypecaller.nf` per apportare due modifiche:

- Aggiungere il parametro `-ERC GVCF` al comando GATK HaplotypeCaller;
- Aggiornare il percorso del file di output per utilizzare l'estensione `.g.vcf` corrispondente, come da convenzione GATK.

Assicuratevi di aggiungere un backslash (`\`) alla fine della riga precedente quando aggiungete `-ERC GVCF`.

=== "Dopo"

    ```groovy title="modules/gatk_haplotypecaller.nf" linenums="22" hl_lines="5 7"
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

    ```groovy title="modules/gatk_haplotypecaller.nf" linenums="22" hl_lines="5"
        """
        gatk HaplotypeCaller \
            -R ${ref_fasta} \
            -I ${input_bam} \
            -O ${input_bam}.vcf \
            -L ${interval_list}
        """
    ```

Dobbiamo anche aggiornare il blocco output per corrispondere alla nuova estensione del file.
Poiché abbiamo cambiato l'output del comando da `.vcf` a `.g.vcf`, il blocco `output:` del processo deve riflettere la stessa modifica.

### 1.2. Aggiornare l'estensione del file di output nel blocco outputs del processo

=== "Dopo"

    ```groovy title="modules/gatk_haplotypecaller.nf" linenums="17" hl_lines="2 3"
        output:
        path "${input_bam}.g.vcf"     , emit: vcf
        path "${input_bam}.g.vcf.idx" , emit: idx
    ```

=== "Prima"

    ```groovy title="modules/gatk_haplotypecaller.nf" linenums="17" hl_lines="2 3"
        output:
        path "${input_bam}.vcf"     , emit: vcf
        path "${input_bam}.vcf.idx" , emit: idx
    ```

Dobbiamo anche aggiornare la configurazione di publish e output del flusso di lavoro per riflettere i nuovi output GVCF.

### 1.3. Aggiornare i target di publish per i nuovi output GVCF

Poiché ora stiamo producendo GVCF invece di VCF, dovremmo aggiornare la sezione `publish:` del flusso di lavoro per utilizzare nomi più descrittivi.
Organizzeremo anche i file GVCF nella loro subdirectory per maggiore chiarezza.

=== "Dopo"

    ```groovy title="genomics.nf" linenums="47" hl_lines="3 4"
        publish:
        indexed_bam = SAMTOOLS_INDEX.out
        gvcf = GATK_HAPLOTYPECALLER.out.vcf
        gvcf_idx = GATK_HAPLOTYPECALLER.out.idx
    ```

=== "Prima"

    ```groovy title="genomics.nf" linenums="47"
        publish:
        indexed_bam = SAMTOOLS_INDEX.out
        vcf = GATK_HAPLOTYPECALLER.out.vcf
        vcf_idx = GATK_HAPLOTYPECALLER.out.idx
    ```

Ora aggiorniamo il blocco output per corrispondere.

### 1.4. Aggiornare il blocco output per la nuova struttura di directory

Dobbiamo anche aggiornare il blocco `output` per inserire i file GVCF in una subdirectory `gvcf`.

=== "Dopo"

    ```groovy title="genomics.nf" linenums="53" hl_lines="3 5 6 8 9"
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

    ```groovy title="genomics.nf" linenums="53"
    output {
        indexed_bam {
            path 'bam'
        }
        vcf {
            path 'vcf'
        }
        vcf_idx {
            path 'vcf'
        }
    }
    ```

Con il modulo, i target di publish e il blocco output tutti aggiornati, possiamo testare le modifiche.

### 1.5. Eseguire la pipeline

Eseguiamo il flusso di lavoro per verificare che le modifiche funzionino.

```bash
nextflow run genomics.nf
```

??? success "Output del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics.nf` [nostalgic_franklin] DSL2 - revision: f2c0a93c6a

    executor >  local (6)
    [cc/fbc705] SAMTOOLS_INDEX (3)       | 3 of 3 ✔
    [27/0d7eb9] GATK_HAPLOTYPECALLER (2) | 3 of 3 ✔
    ```

L'output di Nextflow appare uguale a prima, ma i file `.g.vcf` e i loro file indice sono ora organizzati in subdirectory.

??? abstract "Directory contents (symlink abbreviati)"

    ```console
    results/
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

Se aprite uno dei file GVCF e scorrete al suo interno, potete verificare che GATK HaplotypeCaller abbia prodotto file GVCF come richiesto.

### Takeaway

Quando modificate il nome del file di output di un comando di uno strumento, il blocco `output:` del processo e la configurazione publish/output devono essere aggiornati di conseguenza.

### Cosa c'è dopo?

Impariamo a raccogliere i contenuti di un canale e passarli al processo successivo come singolo input.

---

## 2. Aggiungere uno step di joint genotyping

Ora dobbiamo raccogliere i GVCF per campione, combinarli in un data store GenomicsDB ed eseguire il joint genotyping per produrre un VCF a livello di coorte.
Come trattato nella [Parte 1](01_method.md), questa è un'operazione che utilizza due strumenti: GenomicsDBImport combina i GVCF, poi GenotypeGVCFs produce le chiamate di varianti finali.
Incapsuleremo entrambi gli strumenti in un singolo processo chiamato `GATK_JOINTGENOTYPING`.

Richiamiamo i due comandi dalla [Parte 1](01_method.md):

```bash
gatk GenomicsDBImport \
    -V reads_mother.g.vcf \
    -V reads_father.g.vcf \
    -V reads_son.g.vcf \
    -L /data/ref/intervals.bed \
    --genomicsdb-workspace-path family_trio_gdb
```

```bash
gatk GenotypeGVCFs \
    -R /data/ref/ref.fasta \
    -V gendb://family_trio_gdb \
    -O family_trio.vcf
```

Il primo comando prende i GVCF per campione e un file di intervalli, e produce un data store GenomicsDB.
Il secondo prende quel data store, un genoma di riferimento, e produce il VCF finale a livello di coorte.
L'URI del container è lo stesso di HaplotypeCaller: `community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867`.

### 2.1. Configurare gli input

Il processo di joint genotyping necessita di due tipi di input che non abbiamo ancora: un nome arbitrario per la coorte e gli output GVCF raccolti da tutti i campioni raggruppati insieme.

#### 2.1.1. Aggiungere un parametro `cohort_name`

Dobbiamo fornire un nome arbitrario per la coorte.
Più avanti nella serie di formazione imparerete come utilizzare i metadati dei campioni per questo tipo di cose, ma per ora dichiariamo semplicemente un parametro CLI utilizzando `params` e gli diamo un valore predefinito per comodità.

=== "Dopo"

    ```groovy title="genomics.nf" linenums="14" hl_lines="3-4"
        intervals: Path = "${projectDir}/data/ref/intervals.bed"

        // Base name for final output file
        cohort_name: String = "family_trio"
    }
    ```

=== "Prima"

    ```groovy title="genomics.nf" linenums="14"
        intervals: Path = "${projectDir}/data/ref/intervals.bed"
    }
    ```

#### 2.1.2. Raccogliere gli output di HaplotypeCaller attraverso i campioni

Se collegassimo direttamente il canale di output da `GATK_HAPLOTYPECALLER` al nuovo processo, Nextflow chiamerebbe il processo su ciascun GVCF del campione separatamente.
Vogliamo raggruppare tutti e tre i GVCF (e i loro file indice) in modo che Nextflow li passi tutti insieme a una singola chiamata del processo.

Possiamo farlo utilizzando l'operatore di canale `collect()`.
Aggiungete le seguenti righe al corpo del `workflow`, subito dopo la chiamata a GATK_HAPLOTYPECALLER:

=== "Dopo"

    ```groovy title="genomics.nf" hl_lines="4-6"
            intervals_file
        )

        // Collect variant calling outputs across samples
        all_gvcfs_ch = GATK_HAPLOTYPECALLER.out.vcf.collect()
        all_idxs_ch = GATK_HAPLOTYPECALLER.out.idx.collect()
    ```

=== "Prima"

    ```groovy title="genomics.nf"
            intervals_file
        )
    ```

Analizziamo questo:

1. Prendiamo il canale di output da `GATK_HAPLOTYPECALLER` utilizzando la proprietà `.out`.
2. Poiché abbiamo nominato gli output usando `emit:` nella sezione 1, possiamo selezionare i GVCF con `.vcf` e i file indice con `.idx`. Senza output nominati, dovremmo usare `.out[0]` e `.out[1]`.
3. L'operatore `collect()` raggruppa tutti i file in un singolo elemento, quindi `all_gvcfs_ch` contiene tutti e tre i GVCF insieme, e `all_idxs_ch` contiene tutti e tre i file indice insieme.

Possiamo raccogliere i GVCF e i loro file indice separatamente (invece di mantenerli insieme in tuple) perché Nextflow posizionerà tutti i file di input insieme per l'esecuzione, quindi i file indice saranno presenti accanto ai GVCF.

!!! tip

    Potete utilizzare l'operatore `view()` per ispezionare i contenuti dei canali prima e dopo l'applicazione degli operatori di canale.

### 2.2. Scrivere il processo di joint genotyping e chiamarlo nel flusso di lavoro

Seguendo lo stesso schema che abbiamo utilizzato nella Parte 2, scriveremo la definizione del processo in un file modulo, lo importeremo nel flusso di lavoro e lo chiameremo sugli input che abbiamo appena preparato.

#### 2.2.1. Costruire una stringa per dare a ciascun GVCF un argomento `-V`

Prima di iniziare a compilare la definizione del processo, c'è una cosa da capire.
Il comando GenomicsDBImport si aspetta un argomento `-V` separato per ciascun file GVCF, così:

```bash
gatk GenomicsDBImport \
    -V reads_mother.bam.g.vcf \
    -V reads_father.bam.g.vcf \
    -V reads_son.bam.g.vcf \
    ...
```

Se scrivessimo `-V ${all_gvcfs_ch}`, Nextflow concatenerebbe semplicemente i nomi dei file e quella parte del comando apparirebbe così:

```groovy
-V reads_mother.bam.g.vcf reads_father.bam.g.vcf reads_son.bam.g.vcf
```

Ma dobbiamo che la stringa appaia così:

```groovy
-V reads_mother.bam.g.vcf -V reads_father.bam.g.vcf -V reads_son.bam.g.vcf
```

Cosa importante, dobbiamo costruire questa stringa dinamicamente da qualunque file sia nel canale raccolto.
Nextflow (tramite Groovy) fornisce un modo conciso per farlo:

```groovy
def gvcfs_line = all_gvcfs.collect { gvcf -> "-V ${gvcf}" }.join(' ')
```

Analizziamo questo:

1. `all_gvcfs.collect { gvcf -> "-V ${gvcf}" }` itera su ciascun percorso di file e antepone `-V ` ad esso, producendo `["-V A.g.vcf", "-V B.g.vcf", "-V C.g.vcf"]`.
2. `.join(' ')` li concatena con spazi: `"-V A.g.vcf -V B.g.vcf -V C.g.vcf"`.
3. Il risultato viene assegnato a una variabile locale `gvcfs_line` (definita con `def`), che possiamo interpolare nel template del comando.

Questa riga va all'interno del blocco `script:` del processo, prima del template del comando.
Potete inserire codice Groovy arbitrario tra `script:` e le `"""` di apertura del template del comando.

Poi sarete in grado di riferirvi a quell'intera stringa come `gvcfs_line` nel blocco `script:` del processo.

#### 2.2.2. Compilare il modulo per il processo di joint genotyping

Ora possiamo affrontare la scrittura del processo completo.

Aprite `modules/gatk_jointgenotyping.nf` ed esaminate la struttura della definizione del processo.

Procedete e compilate la definizione del processo utilizzando le informazioni fornite sopra, quindi verificate il vostro lavoro confrontandolo con la soluzione nella scheda "Dopo" qui sotto.

=== "Prima"

    ```groovy title="modules/gatk_jointgenotyping.nf" linenums="1"
    #!/usr/bin/env nextflow

    /*
     * Combina i GVCF in un datastore GenomicsDB ed esegui il joint genotyping per produrre chiamate a livello di coorte
     */
    process GATK_JOINTGENOTYPING {

        container

        input:

        output:

        script:
        """

        """
    }
    ```

=== "Dopo"

    ```groovy title="modules/gatk_jointgenotyping.nf" linenums="1" hl_lines="8 11-17 20-21 24-25 29-33"
    #!/usr/bin/env nextflow

    /*
     * Combina i GVCF in un datastore GenomicsDB ed esegui il joint genotyping per produrre chiamate a livello di coorte
     */
    process GATK_JOINTGENOTYPING {

        container "community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867"

        input:
        path all_gvcfs
        path all_idxs
        path interval_list
        val cohort_name
        path ref_fasta
        path ref_index
        path ref_dict

        output:
        path "${cohort_name}.joint.vcf"     , emit: vcf
        path "${cohort_name}.joint.vcf.idx" , emit: idx

        script:
        def gvcfs_line = all_gvcfs.collect { gvcf -> "-V ${gvcf}" }.join(' ')
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
    }
    ```

Ci sono diverse cose che vale la pena evidenziare qui.

Come in precedenza, diversi input sono elencati anche se i comandi non vi fanno riferimento direttamente: `all_idxs`, `ref_index` e `ref_dict`.
Elencarli assicura che Nextflow posizioni questi file nella directory di lavoro accanto ai file che appaiono nei comandi, che GATK si aspetta di trovare in base alle convenzioni di denominazione.

La variabile `gvcfs_line` utilizza la closure Groovy descritta sopra per costruire gli argomenti `-V` per GenomicsDBImport.

Questo processo esegue due comandi in serie, proprio come fareste nel terminale.
GenomicsDBImport combina i GVCF per campione in un data store, poi GenotypeGVCFs legge quel data store e produce il VCF finale a livello di coorte.
Il data store GenomicsDB (`${cohort_name}_gdb`) è un artefatto intermedio utilizzato solo all'interno del processo; non appare nel blocco output.

Una volta completato questo, il processo è pronto all'uso.
Per utilizzarlo nel flusso di lavoro, dovrete importare il modulo e aggiungere una chiamata al processo.

#### 2.2.3. Importare il modulo

Aggiungete l'istruzione import a `genomics.nf`, sotto le istruzioni import esistenti:

=== "Dopo"

    ```groovy title="genomics.nf" linenums="21" hl_lines="3"
    include { SAMTOOLS_INDEX } from './modules/samtools_index.nf'
    include { GATK_HAPLOTYPECALLER } from './modules/gatk_haplotypecaller.nf'
    include { GATK_JOINTGENOTYPING } from './modules/gatk_jointgenotyping.nf'
    ```

=== "Prima"

    ```groovy title="genomics.nf" linenums="21"
    include { SAMTOOLS_INDEX } from './modules/samtools_index.nf'
    include { GATK_HAPLOTYPECALLER } from './modules/gatk_haplotypecaller.nf'
    ```

Il processo è ora disponibile nello scope del flusso di lavoro.

#### 2.2.4. Aggiungere la chiamata al processo

Aggiungete la chiamata a `GATK_JOINTGENOTYPING` nel corpo del flusso di lavoro, dopo le righe `collect()`:

=== "Dopo"

    ```groovy title="genomics.nf" hl_lines="3-12"
        all_idxs_ch = GATK_HAPLOTYPECALLER.out.idx.collect()

        // Combine GVCFs into a GenomicsDB data store and apply joint genotyping
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

    ```groovy title="genomics.nf"
        all_idxs_ch = GATK_HAPLOTYPECALLER.out.idx.collect()
    ```

Il processo è ora completamente collegato.
Successivamente, configuriamo come vengono pubblicati gli output.

### 2.3. Configurare la gestione degli output

Dobbiamo pubblicare gli output del joint VCF.
Aggiungete target di publish ed entry del blocco output per i risultati del joint genotyping.

#### 2.3.1. Aggiungere target di publish per il joint VCF

Aggiungete il joint VCF e il suo indice alla sezione `publish:` del flusso di lavoro:

=== "Dopo"

    ```groovy title="genomics.nf" hl_lines="5-6"
        publish:
        indexed_bam = SAMTOOLS_INDEX.out
        gvcf = GATK_HAPLOTYPECALLER.out.vcf
        gvcf_idx = GATK_HAPLOTYPECALLER.out.idx
        joint_vcf = GATK_JOINTGENOTYPING.out.vcf
        joint_vcf_idx = GATK_JOINTGENOTYPING.out.idx
    ```

=== "Prima"

    ```groovy title="genomics.nf"
        publish:
        indexed_bam = SAMTOOLS_INDEX.out
        gvcf = GATK_HAPLOTYPECALLER.out.vcf
        gvcf_idx = GATK_HAPLOTYPECALLER.out.idx
    ```

Ora aggiorniamo il blocco output per corrispondere.

#### 2.3.2. Aggiungere entry del blocco output per il joint VCF

Aggiungete entry per i file joint VCF.
Li metteremo alla radice della directory results poiché questo è l'output finale.

=== "Dopo"

    ```groovy title="genomics.nf" hl_lines="11-16"
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
        joint_vcf {
            path '.'
        }
        joint_vcf_idx {
            path '.'
        }
    }
    ```

=== "Prima"

    ```groovy title="genomics.nf"
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

Con il processo, i target di publish e il blocco output tutti al loro posto, possiamo testare il flusso di lavoro completo.

### 2.4. Eseguire il flusso di lavoro

Eseguiamo il flusso di lavoro per verificare che tutto funzioni.

```bash
nextflow run genomics.nf -resume
```

??? success "Output del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    ┃ Launching `genomics.nf` [crazy_marconi] DSL2 - revision: 5da9afc841

    executor >  local (1)
    [9a/c7a873] SAMTOOLS_INDEX (2)       | 3 of 3, cached: 3 ✔
    [e4/4ed55e] GATK_HAPLOTYPECALLER (2) | 3 of 3, cached: 3 ✔
    [a6/7cc8ed] GATK_JOINTGENOTYPING     | 1 of 1 ✔
    ```

I primi due step sono in cache dall'esecuzione precedente, e il nuovo step `GATK_JOINTGENOTYPING` viene eseguito una volta sugli input raccolti da tutti e tre i campioni.
Il file di output finale, `family_trio.joint.vcf` (e il suo indice), sono nella directory results.

??? abstract "Directory contents (symlink abbreviati)"

    ```console
    results/
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

Se aprite il file joint VCF, potete verificare che il flusso di lavoro abbia prodotto le chiamate di varianti attese.

```console title="family_trio.joint.vcf" linenums="40"
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	reads_father	reads_mother	reads_son
20_10037292_10066351	3480	.	C	CT	1625.89	.	AC=5;AF=0.833;AN=6;BaseQRankSum=0.220;DP=85;ExcessHet=0.0000;FS=2.476;MLEAC=5;MLEAF=0.833;MQ=60.00;MQRankSum=0.00;QD=21.68;ReadPosRankSum=-1.147e+00;SOR=0.487	GT:AD:DP:GQ:PL	0/1:15,16:31:99:367,0,375	1/1:0,18:18:54:517,54,0	1/1:0,26:26:78:756,78,0
20_10037292_10066351	3520	.	AT	A	1678.89	.	AC=5;AF=0.833;AN=6;BaseQRankSum=1.03;DP=80;ExcessHet=0.0000;FS=2.290;MLEAC=5;MLEAF=0.833;MQ=60.00;MQRankSum=0.00;QD=22.39;ReadPosRankSum=0.701;SOR=0.730	GT:AD:DP:GQ:PL	0/1:18,13:31:99:296,0,424	1/1:0,18:18:54:623,54,0	1/1:0,26:26:78:774,78,0
20_10037292_10066351	3529	.	T	A	154.29	.	AC=1;AF=0.167;AN=6;BaseQRankSum=-5.440e-01;DP=104;ExcessHet=0.0000;FS=1.871;MLEAC=1;MLEAF=0.167;MQ=60.00;MQRankSum=0.00;QD=7.71;ReadPosRankSum=-1.158e+00;SOR=1.034	GT:AD:DP:GQ:PL	0/0:44,0:44:99:0,112,1347	0/1:12,8:20:99:163,0,328	0/0:39,0:39:99:0,105,1194
```

Ora avete un flusso di lavoro di joint variant calling automatizzato e completamente riproducibile!

!!! note

    Tenete presente che i file di dati che vi abbiamo fornito coprono solo una piccola porzione del cromosoma 20.
    La dimensione reale di un callset di varianti sarebbe contata in milioni di varianti.
    Ecco perché utilizziamo solo piccoli sottoinsiemi di dati per scopi di formazione!

### Takeaway

Sapete come raccogliere output da un canale e raggrupparli come singolo input per un altro processo.
Sapete anche come costruire una riga di comando utilizzando closure Groovy e come eseguire comandi multipli in un singolo processo.

### Cosa c'è dopo?

Datevi una bella pacca sulla spalla! Avete completato il corso Nextflow for Genomics.

Passate al [riepilogo finale del corso](./next_steps.md) per rivedere ciò che avete imparato e scoprire cosa viene dopo.
