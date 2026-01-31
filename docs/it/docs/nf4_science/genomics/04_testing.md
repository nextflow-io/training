# Parte 4: Aggiunta di test

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduzione assistita da IA - [scopri di più e suggerisci miglioramenti](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Nella prima parte di questo corso, avete costruito una pipeline di variant calling completamente lineare che processava i dati di ciascun campione indipendentemente dagli altri.

Nella seconda parte, Le abbiamo mostrato come utilizzare i canali e gli operatori di canale per implementare il joint variant calling con GATK.

Nella terza parte, abbiamo modularizzato la pipeline.

In questa parte della formazione, Le mostreremo come utilizzare [**nf-test**](https://www.nf-test.com/), un framework di testing che si integra bene con Nextflow e rende semplice aggiungere test sia a livello di modulo che a livello di workflow alla vostra pipeline. Per seguire questa parte della formazione, dovrebbe aver completato la Parte 1, la Parte 2 e la Parte 3, nonché la [side quest su nf-test](../../side_quests/nf-test.md), che copre le basi di nf-test e il motivo per cui il testing è importante.

---

## 0. Riscaldamento

!!! note "Nota"

    Assicuratevi di essere nella directory di lavoro corretta:
    `cd /workspaces/training/nf4-science/genomics`

Se ha completato le parti precedenti di questo corso di formazione, dovrebbe avere una versione funzionante della pipeline genomics con l'appropriata struttura di directory dei moduli.

??? abstract "Contenuti della directory"

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
    ```

Questa directory dei moduli può essere trovata nella directory `solutions` se ne avete bisogno.

Inizieremo con lo stesso workflow della Parte 3, che Le abbiamo fornito nel file `genomics-4.nf`. Esattamente come per la [side quest su nf-test](../../side_quests/nf-test.md), aggiungeremo alcuni tipi diversi di test ai tre processi in questa pipeline, nonché un test a livello di workflow.

### 0.1. Verificare che il workflow funzioni

Prima di iniziare ad aggiungere test, assicuratevi che il workflow funzioni come previsto.

```bash
nextflow run genomics-4.nf -resume
```

Questo dovrebbe sembrare molto familiare ormai se avete seguito questo corso di formazione dall'inizio.

??? success "Output del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `genomics-4.nf` [gloomy_poincare] DSL2 - revision: 43203316e0

    executor >  local (7)
    [18/89dfa4] SAMTOOLS_INDEX (1)       | 3 of 3 ✔
    [30/b2522b] GATK_HAPLOTYPECALLER (2) | 3 of 3 ✔
    [a8/d2c189] GATK_JOINTGENOTYPING     | 1 of 1 ✔
    ```

Come in precedenza, ci sarà ora una directory `work` e una directory `results_genomics` all'interno della vostra directory di progetto. Utilizzeremo effettivamente questi risultati più avanti nei nostri test. Ma da ora in poi utilizzeremo il pacchetto `nf-test` per testare la pipeline.

### 0.2. Inizializzare `nf-test`

Come per la [side quest su nf-test](../../side_quests/nf-test.md), dobbiamo inizializzare il pacchetto `nf-test`.

```bash
nf-test init
```

??? success "Output del comando"

    ```console
    🚀 nf-test 0.9.3
    https://www.nf-test.com
    (c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr

    Project configured. Configuration is stored in nf-test.config
    ```

??? abstract "Contenuti di nf-test.config"

    ```groovy title="nf-test.config"
    config {

        testsDir "tests"
        workDir ".nf-test"
        configFile "tests/nextflow.config"
        profile ""

    }
    ```

Crea anche una directory `tests` contenente uno stub del file di configurazione.

### Takeaway

Ora siamo pronti per iniziare a scrivere test per la nostra pipeline genomics.

### Quali sono i prossimi passi?

Scrivere test di base che valutano se le chiamate ai processi sono riuscite e hanno prodotto gli output corretti.

---

## 1. Testare un processo per successo e corrispondenza degli output

Inizieremo testando il processo `SAMTOOLS_INDEX`, che crea file indice per i file BAM per consentire un accesso casuale efficiente. Questo è un buon primo caso di test perché:

1. Ha un singolo input ben definito (un file BAM)
2. Produce un output prevedibile (un file indice BAI)
3. L'output dovrebbe essere identico per input identici

### 1.1. Generare uno stub del file di test

Prima, generiamo uno stub del file di test:

```bash
nf-test generate process modules/samtools/index/main.nf
```

??? success "Output del comando"

    ```console
    Load source file '/workspaces/training/nf4-science/genomics/modules/samtools/index/main.nf'
    Wrote process test file '/workspaces/training/nf4-science/genomics/tests/modules/samtools/index/main.nf.test'

    SUCCESS: Generated 1 test files.
    ```

Questo crea un file nella stessa directory di `main.nf`.
Può navigare alla directory nel file explorer e aprire il file, che dovrebbe contenere il seguente codice:

```groovy title="tests/modules/samtools/index/main.nf.test" linenums="1"
nextflow_process {

    name "Test Process SAMTOOLS_INDEX"
    script "modules/samtools/index/main.nf"
    process "SAMTOOLS_INDEX"

    test("Should run without failures") {

        when {
            params {
                // define parameters here. Example:
                // outdir = "tests/results"
            }
            process {
                """
                // define inputs of the process here. Example:
                // input[0] = file("test-file.txt")
                """
            }
        }

        then {
            assert process.success
            assert snapshot(process.out).match()
        }

    }

}
```

Le asserzioni iniziali dovrebbero essere familiari dalla [side quest su nf-test](../../side_quests/nf-test.md):

- `assert process.success` afferma che ci aspettiamo che il processo venga eseguito con successo e si completi senza errori.
- `snapshot(process.out).match()` afferma che ci aspettiamo che il risultato dell'esecuzione sia identico al risultato ottenuto in un'esecuzione precedente (se applicabile).
  Discuteremo questo in dettaglio più avanti.

Utilizzando questo come punto di partenza, dobbiamo aggiungere gli input di test corretti per il processo samtools index, e qualsiasi parametro se applicabile.

### 1.2. Spostare il file di test e aggiornare il percorso dello script

Prima di metterci al lavoro per compilare il test, dobbiamo spostare il file nella sua posizione definitiva. Parte del motivo per cui abbiamo aggiunto una directory per ciascun modulo è che ora possiamo fornire i test in una directory `tests` co-localizzata con il file `main.nf` di ciascun modulo. Crei quella directory e sposti il file di test lì.

```bash
mkdir -p modules/samtools/index/tests
mv tests/modules/samtools/index/main.nf.test modules/samtools/index/tests/
```

Ora possiamo semplificare la sezione `script` del file di test con un percorso relativo:

=== "Dopo"

    ```groovy title="modules/samtools/index/tests/main.nf.test" linenums="3" hl_lines="2"
    name "Test Process SAMTOOLS_INDEX"
    script "../main.nf"
    process "SAMTOOLS_INDEX"
    ```

=== "Prima"

    ```groovy title="modules/samtools/index/tests/main.nf.test" linenums="3" hl_lines="2"
    name "Test Process SAMTOOLS_INDEX"
    script "modules/samtools/index/main.nf"
    process "SAMTOOLS_INDEX"
    ```

Questo indica al test dove trovare il file `main.nf` del modulo, senza dover specificare il percorso completo.

### 1.3. Fornire input di test per SAMTOOLS_INDEX

Il file stub include un segnaposto che dobbiamo sostituire con un input di test effettivo, appropriato all'input di `samtools index`. L'input appropriato è un file BAM, che abbiamo disponibile nella directory `data/bam`.

=== "Dopo"

    ```groovy title="modules/samtools/index/tests/main.nf.test" linenums="14"
    process {
        """
        input[0] = file("${projectDir}/data/bam/reads_son.bam")
        """
    }
    ```

=== "Prima"

    ```groovy title="modules/samtools/index/tests/main.nf.test" linenums="14"
    process {
        """
        // define inputs of the process here. Example:
        // input[0] = file("test-file.txt")
        """
    }
    ```

### 1.4. Denominare il test in base alla funzionalità

Come abbiamo appreso in precedenza, è buona pratica rinominare il test con qualcosa che abbia senso nel contesto del test.

=== "Dopo"

    ```groovy title="modules/samtools/index/tests/main.nf.test" linenums="7"
    test("Should index reads_son.bam correctly") {
    ```

    Questo accetta una stringa arbitraria, quindi potremmo inserire qualsiasi cosa vogliamo.
    Qui scegliamo di fare riferimento al nome del file e al suo formato.

=== "Prima"

    ```groovy title="modules/samtools/index/tests/main.nf.test" linenums="7"
    test("Should run without failures") {
    ```

### 1.5. Eseguire il test ed esaminare l'output

Esegua il test:

```bash
nf-test test modules/samtools/index/tests/main.nf.test
```

??? success "Output del comando"

    ```console
    🚀 nf-test 0.9.3
    https://www.nf-test.com
    (c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


    Test Process SAMTOOLS_INDEX

      Test [625e39ee] 'Should index reads_son.bam correctly' PASSED (7.717s)
      Snapshots:
        1 created [Should index reads_son.bam correctly]


    Snapshot Summary:
      1 created

    SUCCESS: Executed 1 tests in 7.727s
    ```

Come abbiamo appreso in precedenza, questo ha verificato l'asserzione di base sul successo del processo e ha creato un file snapshot basato sull'output del processo. Possiamo vedere i contenuti del file snapshot nel file `tests/modules/samtools/index/tests/main.nf.test.snap`:

```json title="modules/samtools/index/tests/main.nf.test.snap" linenums="1"
{
  "Should index reads_son.bam correctly": {
    "content": [
      {
        "0": [
          [
            "reads_son.bam:md5,af5956d9388ba017944bef276b71d809",
            "reads_son.bam.bai:md5,a2ca7b84998218ee77eff14af8eb8ca2"
          ]
        ]
      }
    ],
    "meta": {
      "nf-test": "0.9.3",
      "nextflow": "25.10.2"
    },
    "timestamp": "2026-01-27T15:09:48.394063389"
  }
}
```

Possiamo anche eseguire nuovamente il test e vedere che passa, perché l'output è identico allo snapshot:

```bash
nf-test test modules/samtools/index/tests/main.nf.test
```

??? success "Output del comando"

    ```console
    🚀 nf-test 0.9.3
    https://www.nf-test.com
    (c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


    Test Process SAMTOOLS_INDEX

      Test [625e39ee] 'Should index reads_son.bam correctly' PASSED (7.938s)


    SUCCESS: Executed 1 tests in 7.987s
    ```

### 1.6. Aggiungere più test a `SAMTOOLS_INDEX`

A volte è utile testare una gamma di file di input diversi per assicurarsi di testare una varietà di potenziali problemi. Aggiunga test per i file BAM della madre e del padre nel trio dei nostri dati di test.

```groovy
    test("Should index reads_mother.bam correctly") {

        when {
            params {
                // define parameters here
            }
            process {
                """
                input[0] = file("${projectDir}/data/bam/reads_mother.bam")
                """
            }
        }

        then {
            assert process.success
            assert snapshot(process.out).match()
        }

    }

    test("Should index reads_father.bam correctly") {

        when {
            params {
                // define parameters here
            }
            process {
                """
                input[0] = file("${projectDir}/data/bam/reads_father.bam")
                """
            }
        }

        then {
            assert process.success
            assert snapshot(process.out).match()
        }

    }
```

Quindi potete eseguire nuovamente il test:

```bash
nf-test test modules/samtools/index/tests/main.nf.test
```

??? success "Output del comando"

    ```console
    🚀 nf-test 0.9.3
    https://www.nf-test.com
    (c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


    Test Process SAMTOOLS_INDEX

      Test [625e39ee] 'Should index reads_son.bam correctly' PASSED (7.185s)
      Test [a8b28f36] 'Should index reads_mother.bam correctly' PASSED (6.576s)
      Test [c15852a1] 'Should index reads_father.bam correctly' PASSED (6.31s)
      Snapshots:
        2 created [Should index reads_father.bam correctly, Should index reads_mother.bam correctly]


    Snapshot Summary:
      2 created

    SUCCESS: Executed 3 tests in 20.117s
    ```

Noti l'avviso, che si riferisce all'effetto del parametro `--update-snapshot`.

!!! note "Nota"

    Qui stiamo utilizzando dati di test che abbiamo usato in precedenza per dimostrare gli output scientifici della pipeline.
    Se avessimo pianificato di utilizzare questi test in un ambiente di produzione, avremmo generato input più piccoli per scopi di testing.

    In generale è importante mantenere i test unitari il più leggeri possibile utilizzando i pezzi di dati più piccoli necessari e sufficienti per valutare la funzionalità del processo, altrimenti il tempo di esecuzione totale può accumularsi in modo significativo.
    Una suite di test che richiede troppo tempo per essere eseguita regolarmente è una suite di test che probabilmente verrà saltata nell'interesse della rapidità.

### Takeaway

Ha scritto il vostro primo test di modulo per un processo genomics, verificando che `SAMTOOLS_INDEX` crei correttamente file indice per diversi file BAM. La suite di test garantisce che:

1. Il processo venga eseguito con successo
2. I file indice vengano creati
3. Gli output siano coerenti tra le esecuzioni
4. Il processo funzioni per tutti i file BAM dei campioni

### Quali sono i prossimi passi?

Imparare a scrivere test per altri processi nel nostro workflow genomics, utilizzando il metodo setup per gestire processi concatenati. Valuteremo anche se gli output, in particolare i nostri file VCF, contengono le chiamate di varianti attese.

---

## 2. Aggiungere test a un processo concatenato e testare i contenuti

Per testare `GATK_HAPLOTYPECALLER`, dobbiamo fornire al processo l'output di `SAMTOOLS_INDEX` come input. Potremmo farlo eseguendo `SAMTOOLS_INDEX`, recuperando i suoi output e memorizzandoli con i dati di test per il workflow. Questo è in realtà l'approccio raccomandato per una pipeline rifinita, ma nf-test fornisce un approccio alternativo, utilizzando il metodo `setup`.

Con il metodo setup, possiamo attivare il processo `SAMTOOLS_INDEX` come parte della configurazione del test, e quindi utilizzare il suo output come input per `GATK_HAPLOTYPECALLER`. Questo ha un costo: dovremo eseguire il processo `SAMTOOLS_INDEX` ogni volta che eseguiamo il test per `GATK_HAPLOTYPECALLER`. Tuttavia, forse stiamo ancora sviluppando il workflow e non vogliamo pre-generare dati di test che potremmo dover cambiare successivamente. Il processo `SAMTOOLS_INDEX` è anche molto veloce, quindi forse i benefici di pre-generare e memorizzare i suoi output sono trascurabili. Ecco come funziona il metodo setup.

### 2.1. Generare e posizionare il file di test

Come in precedenza, prima generiamo lo stub del file:

```bash
nf-test generate process modules/gatk/haplotypecaller/main.nf
```

??? success "Output del comando"

    ```console
    Load source file '/workspaces/training/nf4-science/genomics/modules/gatk/haplotypecaller/main.nf'
    Wrote process test file '/workspaces/training/nf4-science/genomics/tests/modules/gatk/haplotypecaller/main.nf.test'

    SUCCESS: Generated 1 test files.
    ```

Questo produce il seguente stub di test:

```groovy title="tests/modules/gatk/haplotypecaller/main.nf.test" linenums="1"
nextflow_process {

    name "Test Process GATK_HAPLOTYPECALLER"
    script "modules/gatk/haplotypecaller/main.nf"
    process "GATK_HAPLOTYPECALLER"

    test("Should run without failures") {

        when {
            params {
                // define parameters here. Example:
                // outdir = "tests/results"
            }
            process {
                """
                // define inputs of the process here. Example:
                // input[0] = file("test-file.txt")
                """
            }
        }

        then {
            assert process.success
            assert snapshot(process.out).match()
        }

    }

}
```

### 2.2. Spostare il file di test e aggiornare il percorso dello script

Creiamo una directory per il file di test co-localizzata con il file `main.nf` del modulo:

```bash
mkdir -p modules/gatk/haplotypecaller/tests
```

E spostiamo il file stub di test lì:

```bash
mv tests/modules/gatk/haplotypecaller/main.nf.test modules/gatk/haplotypecaller/tests/
```

Infine, non dimentichi di aggiornare il percorso dello script:

=== "Dopo"

    ```groovy title="modules/gatk/haplotypecaller/tests/main.nf.test" linenums="3" hl_lines="2"
        name "Test Process GATK_HAPLOTYPECALLER"
        script "../main.nf"
        process "GATK_HAPLOTYPECALLER"
    ```

=== "Prima"

    ```groovy title="modules/gatk/haplotypecaller/tests/main.nf.test" linenums="3" hl_lines="2"
        name "Test Process GATK_HAPLOTYPECALLER"
        script "modules/gatk/haplotypecaller/main.nf"
        process "GATK_HAPLOTYPECALLER"
    ```

### 2.3. Fornire input utilizzando il metodo setup

Inseriamo un blocco `setup` prima del blocco `when`, dove possiamo attivare un'esecuzione del processo `SAMTOOLS_INDEX` su uno dei nostri file di input originali. Inoltre, ricordate come prima di cambiare il nome del test in qualcosa di significativo.

=== "Dopo"

    ```groovy title="modules/gatk/haplotypecaller/tests/main.nf.test" linenums="7" hl_lines="1-12"
        test("Should call son's haplotype correctly") {

            setup {
                run("SAMTOOLS_INDEX") {
                    script "../../../samtools/index/main.nf"
                    process {
                        """
                        input[0] =  file("${projectDir}/data/bam/reads_son.bam")
                        """
                    }
                }
            }

            when {
    ```

=== "Prima"

    ```groovy title="modules/gatk/haplotypecaller/tests/main.nf.test" linenums="7"  hl_lines="1"
    test("Should run without failures") {

        when {
    ```

Quindi possiamo fare riferimento all'output di quel processo nel blocco `when` dove specifichiamo gli input del test:

```groovy title="modules/gatk/haplotypecaller/tests/main.nf.test" linenums="20"
        when {
            params {
                // define parameters here
            }
            process {
                """
                input[0] = SAMTOOLS_INDEX.out
                input[1] = file("${projectDir}/data/ref/ref.fasta")
                input[2] = file("${projectDir}/data/ref/ref.fasta.fai")
                input[3] = file("${projectDir}/data/ref/ref.dict")
                input[4] = file("${projectDir}/data/ref/intervals.bed")
                """
            }
        }
```

Effettuate quella modifica ed eseguite nuovamente il test:

```bash
nf-test test modules/gatk/haplotypecaller/tests/main.nf.test
```

??? success "Output del comando"

    ```console
    🚀 nf-test 0.9.3
    https://www.nf-test.com
    (c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


    Test Process GATK_HAPLOTYPECALLER

      Test [c5156c2b] 'Should call son's haplotype correctly' PASSED (40.53s)
      Snapshots:
        1 created [Should call son's haplotype correctly]


    Snapshot Summary:
      1 created

    SUCCESS: Executed 1 tests in 40.555s
    ```

Produce anche un file snapshot come prima.

### 2.4. Eseguire nuovamente e osservare il fallimento

Interessante, se esegue esattamente lo stesso comando di nuovo, questa volta il test fallirà.

```bash
nf-test test modules/gatk/haplotypecaller/tests/main.nf.test
```

??? failure "Output del comando"

    ```console
    🚀 nf-test 0.9.3
    https://www.nf-test.com
    (c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


    Test Process GATK_HAPLOTYPECALLER

      Test [c5156c2b] 'Should call son's haplotype correctly' FAILED (40.123s)

      java.lang.RuntimeException: Different Snapshot:
      [                                                                                           [
          {                                                                                           {
              "0": [                                                                                      "0": [
                  "reads_son.bam.g.vcf:md5,069316cdd4328542ffc6ae247b1dac39"                         |                 "reads_son.bam.g.vcf:md5,005f1a13ee39f11b0fc9bea094850eac"
              ],                                                                                          ],
              "1": [                                                                                      "1": [
                  "reads_son.bam.g.vcf.idx:md5,dc36c18f2afdc546f41e68b2687e9334"                     |                 "reads_son.bam.g.vcf.idx:md5,dbad4b76a4b90c158ffc9c9740764242"
              ],                                                                                          ],
              "idx": [                                                                                    "idx": [
                  "reads_son.bam.g.vcf.idx:md5,dc36c18f2afdc546f41e68b2687e9334"                     |                 "reads_son.bam.g.vcf.idx:md5,dbad4b76a4b90c158ffc9c9740764242"
              ],                                                                                          ],
              "vcf": [                                                                                    "vcf": [
                  "reads_son.bam.g.vcf:md5,069316cdd4328542ffc6ae247b1dac39"                         |                 "reads_son.bam.g.vcf:md5,005f1a13ee39f11b0fc9bea094850eac"
              ]                                                                                           ]
          }                                                                                           }
      ]                                                                                           ]

      Nextflow stdout:

      Nextflow stderr:


        Obsolete snapshots can only be checked if all tests of a file are executed successful.


    FAILURE: Executed 1 tests in 40.156s (1 failed)
    ```

Il messaggio di errore Le indica che c'erano differenze tra gli snapshot per le due esecuzioni; in particolare, i valori md5sum sono diversi per i file VCF.

Perché? Per farla breve, lo strumento HaplotypeCaller include un timestamp nell'intestazione del VCF che è diverso ogni volta (per definizione).
Di conseguenza, non possiamo semplicemente aspettarci che i file abbiano md5sum identici anche se hanno contenuti identici in termini delle chiamate di varianti stesse.

Come affrontiamo questo?

### 2.5. Utilizzare un metodo di asserzione del contenuto per verificare una variante specifica

Un modo per risolvere il problema è utilizzare un [diverso tipo di asserzione](https://nf-co.re/docs/contributing/tutorials/nf-test_assertions).
In questo caso, verificheremo un contenuto specifico invece di asserire l'identità.
Più esattamente, faremo leggere allo strumento le righe del file VCF e verificheremo l'esistenza di righe specifiche.

In pratica, sostituiamo la seconda asserzione nel blocco `then` come segue:

=== "Dopo"

    ```groovy title="modules/gatk/haplotypecaller/tests/main.nf.test" linenums="35" hl_lines="3 4"
            then {
                assert process.success
                assert path(process.out[0][0]).readLines().contains('#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	reads_son')
                assert path(process.out[0][0]).readLines().contains('20_10037292_10066351	3277	.	G	<NON_REF>	.	.	END=3282	GT:DP:GQ:MIN_DP:PL	0/0:25:72:24:0,72,719')
            }
    ```

=== "Prima"

    ```groovy title="modules/gatk/haplotypecaller/tests/main.nf.test" linenums="35" hl_lines="3"
    then {
        assert process.success
        assert snapshot(process.out).match()
    }
    ```

Qui stiamo leggendo l'intero contenuto del file di output VCF e cercando una corrispondenza di contenuto, il che va bene su un piccolo file di test, ma non vorrebbe farlo su un file più grande.
Potrebbe invece scegliere di leggere righe specifiche.

Questo approccio richiede di scegliere più attentamente cosa vogliamo usare come 'segnale' da testare.
Il lato positivo è che può essere usato per testare con grande precisione se uno strumento di analisi può identificare in modo coerente caratteristiche 'difficili' (come varianti rare) mentre subisce ulteriore sviluppo.

### 2.6. Eseguire nuovamente e osservare il successo

Una volta che abbiamo modificato il test in questo modo, possiamo eseguire il test più volte, e passerà costantemente.

```bash
nf-test test modules/gatk/haplotypecaller/tests/main.nf.test
```

??? success "Output del comando"

    ```console
    🚀 nf-test 0.9.3
    https://www.nf-test.com
    (c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


    Test Process GATK_HAPLOTYPECALLER

      Test [c5156c2b] 'Should call son's haplotype correctly' PASSED (40.53s)


    SUCCESS: Executed 1 tests in 40.555s
    ```

### 2.7. Aggiungere più test

Aggiunga test simili per i campioni della madre e del padre:

```groovy title="modules/gatk/haplotypecaller/tests/main.nf.test" linenums="43"
    test("Should call mother's haplotype correctly") {

        setup {
            run("SAMTOOLS_INDEX") {
                script "../../../samtools/index/main.nf"
                process {
                    """
                    input[0] =  file("${projectDir}/data/bam/reads_mother.bam")
                    """
                }
            }
        }

        when {
            params {
                // define parameters here
            }
            process {
                """
                input[0] = SAMTOOLS_INDEX.out
                input[1] = file("${projectDir}/data/ref/ref.fasta")
                input[2] = file("${projectDir}/data/ref/ref.fasta.fai")
                input[3] = file("${projectDir}/data/ref/ref.dict")
                input[4] = file("${projectDir}/data/ref/intervals.bed")
                """
            }
        }

        then {
            assert process.success
            assert path(process.out[0][0]).readLines().contains('#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	reads_mother')
            assert path(process.out[0][0]).readLines().contains('20_10037292_10066351	3277	.	G	<NON_REF>	.	.	END=3278	GT:DP:GQ:MIN_DP:PL	0/0:38:99:37:0,102,1530')
        }
    }

    test("Should call father's haplotype correctly") {

        setup {
            run("SAMTOOLS_INDEX") {
                script "../../../samtools/index/main.nf"
                process {
                    """
                    input[0] =  file("${projectDir}/data/bam/reads_father.bam")
                    """
                }
            }
        }

        when {
            params {
                // define parameters here
            }
            process {
                """
                input[0] = SAMTOOLS_INDEX.out
                input[1] = file("${projectDir}/data/ref/ref.fasta")
                input[2] = file("${projectDir}/data/ref/ref.fasta.fai")
                input[3] = file("${projectDir}/data/ref/ref.dict")
                input[4] = file("${projectDir}/data/ref/intervals.bed")
                """
            }
        }

        then {
            assert process.success
            assert path(process.out[0][0]).readLines().contains('#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	reads_father')
            assert path(process.out[0][0]).readLines().contains('20_10037292_10066351	3277	.	G	<NON_REF>	.	.	END=3281	GT:DP:GQ:MIN_DP:PL	0/0:44:99:42:0,120,1800')
        }
    }
```

### 2.8. Eseguire il comando di test

```bash
nf-test test modules/gatk/haplotypecaller/tests/main.nf.test
```

??? success "Output del comando"

    ```console
    🚀 nf-test 0.9.3
    https://www.nf-test.com
    (c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


    Test Process GATK_HAPLOTYPECALLER

      Test [c5156c2b] 'Should call son's haplotype correctly' PASSED (40.53s)
      Test [10de94a8] 'Should call mother's haplotype correctly' PASSED (41.47s)
      Test [c0386fc7] 'Should call father's haplotype correctly' PASSED (45.556s)


    SUCCESS: Executed 3 tests in 127.586s
    ```

Questo completa il piano di test di base per questo secondo passaggio nella pipeline. Avanti verso il terzo e ultimo test a livello di modulo!

### Takeaway

Ha imparato a:

1. Testare processi che dipendono da output di altri processi
2. Verificare varianti genomiche specifiche nei file di output VCF
3. Gestire output non deterministici verificando contenuti specifici
4. Testare il variant calling su più campioni

### Quali sono i prossimi passi?

Imparare a scrivere test che utilizzano dati di test pre-generati per il passaggio di joint genotyping.

---

## 3. Utilizzare dati di test pre-generati

Per il passaggio di joint genotyping, utilizzeremo un approccio diverso - utilizzando dati di test pre-generati. Questo è spesso preferibile per:

1. Processi complessi con più dipendenze
2. Processi che richiedono molto tempo per essere eseguiti
3. Processi che fanno parte di una pipeline stabile e di produzione

### 3.1. Generare dati di test

Ispezioni i risultati che abbiamo generato all'inizio di questa sezione:

```bash
tree results_genomics/
```

```console title="Contenuti della directory dei risultati"
results_genomics/
├── family_trio.joint.vcf
├── family_trio.joint.vcf.idx
├── gvcf
│   ├── reads_father.bam.g.vcf -> /workspaces/training/nf4-science/genomics/work/30/b2522b83c63baff8c3cf75704512a2/reads_father.bam.g.vcf
│   ├── reads_father.bam.g.vcf.idx -> /workspaces/training/nf4-science/genomics/work/30/b2522b83c63baff8c3cf75704512a2/reads_father.bam.g.vcf.idx
│   ├── reads_mother.bam.g.vcf -> /workspaces/training/nf4-science/genomics/work/f6/be2efa58e625d08cf8d0da1d0e9f09/reads_mother.bam.g.vcf
│   ├── reads_mother.bam.g.vcf.idx -> /workspaces/training/nf4-science/genomics/work/f6/be2efa58e625d08cf8d0da1d0e9f09/reads_mother.bam.g.vcf.idx
│   ├── reads_son.bam.g.vcf -> /workspaces/training/nf4-science/genomics/work/fe/2f22d56aa16ed45f8bc419312894f6/reads_son.bam.g.vcf
│   └── reads_son.bam.g.vcf.idx -> /workspaces/training/nf4-science/genomics/work/fe/2f22d56aa16ed45f8bc419312894f6/reads_son.bam.g.vcf.idx
└── indexed_bam
    ├── reads_father.bam -> /workspaces/training/nf4-science/genomics/work/42/a3bf19dbfaf1f3672b16a5d5e6a8be/reads_father.bam
    ├── reads_father.bam.bai -> /workspaces/training/nf4-science/genomics/work/cf/289c2d264f496d60a69e3e9ba6463e/reads_father.bam.bai
    ├── reads_mother.bam -> /workspaces/training/nf4-science/genomics/work/af/f31a6ade82cc0cf853c4f61c8bc473/reads_mother.bam
    ├── reads_mother.bam.bai -> /workspaces/training/nf4-science/genomics/work/18/89dfa40a3def17e45421e54431a126/reads_mother.bam.bai
    ├── reads_son.bam -> /workspaces/training/nf4-science/genomics/work/9f/9615dd553d6f13d8bec4f006ac395f/reads_son.bam
    └── reads_son.bam.bai -> /workspaces/training/nf4-science/genomics/work/4d/cb384a97db5687cc9daab002017c7c/reads_son.bam.bai

2 directories, 14 files
```

Il passaggio di joint genotyping necessita dei file VCF prodotti dai passaggi di haplotype caller come input, insieme agli indici. Quindi copiamo i risultati che abbiamo nella directory di test del modulo `jointgenotyping`.

```bash
mkdir -p modules/gatk/jointgenotyping/tests/inputs/
cp results_genomics/gvcf/*.g.vcf results_genomics/gvcf/*.g.vcf.idx modules/gatk/jointgenotyping/tests/inputs/
```

Ora possiamo utilizzare questi file come input per il test che stiamo per scrivere per il passaggio di joint genotyping.

### 3.2. Generare lo stub del file di test

Come in precedenza, prima generiamo lo stub del file:

```bash
nf-test generate process modules/gatk/jointgenotyping/main.nf
```

??? success "Output del comando"

    ```console
    Load source file '/workspaces/training/nf4-science/genomics/modules/gatk/jointgenotyping/main.nf'
    Wrote process test file '/workspaces/training/nf4-science/genomics/tests/modules/gatk/jointgenotyping/main.nf.test'

    SUCCESS: Generated 1 test files.
    ```

Questo produce il seguente stub di test:

```groovy title="tests/modules/gatk/jointgenotyping/main.nf.test" linenums="1"
nextflow_process {

    name "Test Process GATK_JOINTGENOTYPING"
    script "modules/gatk/jointgenotyping/main.nf"
    process "GATK_JOINTGENOTYPING"

    test("Should run without failures") {

        when {
            params {
                // define parameters here. Example:
                // outdir = "tests/results"
            }
            process {
                """
                // define inputs of the process here. Example:
                // input[0] = file("test-file.txt")
                """
            }
        }

        then {
            assert process.success
            assert snapshot(process.out).match()
        }

    }

}
```

### 3.3. Spostare il file di test e aggiornare il percorso dello script

Questa volta abbiamo già una directory per i test co-localizzata con il file `main.nf` del modulo, quindi possiamo spostare il file stub di test lì:

```bash
mv tests/modules/gatk/jointgenotyping/main.nf.test modules/gatk/jointgenotyping/tests/
```

E non dimentichi di aggiornare il percorso dello script:

=== "Dopo"

    ```groovy title="modules/gatk/jointgenotyping/tests/main.nf.test" linenums="3" hl_lines="2"
    name "Test Process GATK_JOINTGENOTYPING"
    script "../main.nf"
    process "GATK_JOINTGENOTYPING"
    ```

=== "Prima"

    ```groovy title="modules/gatk/jointgenotyping/tests/main.nf.test" linenums="3" hl_lines="2"
    name "Test Process GATK_JOINTGENOTYPING"
    script "modules/gatk/jointgenotyping/main.nf"
    process "GATK_JOINTGENOTYPING"
    ```

### 3.4. Fornire input

Compili gli input in base alle definizioni degli input del processo e rinomini il test di conseguenza:

```groovy title="modules/gatk/jointgenotyping/tests/main.nf.test" linenums="7"
    test("Should call trio's joint genotype correctly") {

        when {
            params {
                // define parameters here
            }
            process {
                """
                input[0] = [
                    file("${projectDir}/modules/gatk/jointgenotyping/tests/inputs/reads_father.bam.g.vcf"),
                    file("${projectDir}/modules/gatk/jointgenotyping/tests/inputs/reads_mother.bam.g.vcf"),
                    file("${projectDir}/modules/gatk/jointgenotyping/tests/inputs/reads_son.bam.g.vcf")
                ]
                input[1] = [
                    file("${projectDir}/modules/gatk/jointgenotyping/tests/inputs/reads_father.bam.g.vcf.idx"),
                    file("${projectDir}/modules/gatk/jointgenotyping/tests/inputs/reads_mother.bam.g.vcf.idx"),
                    file("${projectDir}/modules/gatk/jointgenotyping/tests/inputs/reads_son.bam.g.vcf.idx")
                ]
                input[2] = file("${projectDir}/data/ref/intervals.bed")
                input[3] = "family_trio"
                input[4] = file("${projectDir}/data/ref/ref.fasta")
                input[5] = file("${projectDir}/data/ref/ref.fasta.fai")
                input[6] = file("${projectDir}/data/ref/ref.dict")
                """
            }
        }
```

### 3.5. Utilizzare asserzioni di contenuto

L'output del passaggio di joint genotyping è un altro file VCF, quindi utilizzeremo nuovamente un'asserzione di contenuto.

```groovy title="modules/gatk/jointgenotyping/tests/main.nf.test" linenums="25"
    then {
        assert process.success
        assert path(process.out[0][0]).readLines().contains('#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	reads_father	reads_mother	reads_son')
        assert path(process.out[0][0]).readLines().contains('20_10037292_10066351	3480	.	C	CT	1625.89	.	AC=5;AF=0.833;AN=6;BaseQRankSum=0.220;DP=85;ExcessHet=0.0000;FS=2.476;MLEAC=5;MLEAF=0.833;MQ=60.00;MQRankSum=0.00;QD=21.68;ReadPosRankSum=-1.147e+00;SOR=0.487	GT:AD:DP:GQ:PL	0/1:15,16:31:99:367,0,375	1/1:0,18:18:54:517,54,0	1/1:0,26:26:78:756,78,0')
    }
```

Verificando il contenuto di una variante specifica nel file di output, questo test verifica che:

1. Il processo di joint genotyping venga eseguito con successo
2. L'output VCF contenga tutti e tre i campioni nell'ordine corretto
3. Una variante specifica sia chiamata correttamente con:
   - Genotipi accurati per ciascun campione (0/1 per il padre, 1/1 per madre e figlio)
   - Profondità di lettura e qualità dei genotipi corrette
   - Statistiche a livello di popolazione come la frequenza allelica (AF=0.833)

Non abbiamo fatto uno snapshot dell'intero file, ma verificando una variante specifica, possiamo essere sicuri che il processo di joint genotyping funzioni come previsto.

### 3.6. Eseguire il test

```bash
nf-test test modules/gatk/jointgenotyping/tests/main.nf.test
```

??? success "Output del comando"

    ```console
    🚀 nf-test 0.9.3
    https://www.nf-test.com
    (c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


    Test Process GATK_JOINTGENOTYPING

      Test [ac2067de] 'Should call trio's joint genotype correctly' PASSED (53.827s)


    SUCCESS: Executed 1 tests in 53.837s
    ```

Il test passa, verificando che il nostro processo di joint genotyping correttamente:

1. Combina i VCF dei singoli campioni
2. Esegue il variant calling congiunto
3. Produce un VCF multi-campione con chiamate genotipiche coerenti tra le esecuzioni

### Takeaway

Sa come:

- Utilizzare risultati generati in precedenza come input per i test
- Scrivere test utilizzando dati di test pre-generati

### Quali sono i prossimi passi?

Aggiungere un test a livello di workflow per verificare che l'intera pipeline di variant calling funzioni end-to-end.

---

## 4. Aggiungere un test a livello di workflow

Ora testeremo la pipeline completa di variant calling, dai file BAM ai genotipi congiunti. Questo verifica che:

1. Tutti i processi funzionino correttamente insieme
2. I dati fluiscano correttamente tra i passaggi
3. Le chiamate di varianti finali siano coerenti

### 4.1. Generare il test del workflow

Generi un file di test per la pipeline completa:

```bash
nf-test generate pipeline genomics-4.nf
```

??? success "Output del comando"

    ```console
    Load source file '/workspaces/training/nf4-science/genomics/genomics-4.nf'
    Wrote pipeline test file '/workspaces/training/nf4-science/genomics/tests/genomics-4.nf.test'

    SUCCESS: Generated 1 test files.
    ```

Questo crea uno stub di test di base:

```groovy title="tests/genomics-4.nf.test" linenums="1"
nextflow_pipeline {

    name "Test Workflow genomics-4.nf"
    script "genomics-4.nf"

    test("Should run without failures") {

        when {
            params {
                // define parameters here. Example:
                // outdir = "tests/results"
            }
        }

        then {
            assert workflow.success
        }

    }

}
```

Corregga solo il nome in qualcosa di significativo (vedrà perché questo è utile a breve).

=== "Dopo"

    ```groovy title="tests/genomics-4.nf.test" linenums="6" hl_lines="1"
        test("Should run the pipeline without failures") {
    ```

=== "Prima"

    ```groovy title="tests/genomics-4.nf.test" linenums="6" hl_lines="1"
        test("Should run without failures") {
    ```

!!! note "Nota"

    In questo caso il file di test può rimanere dove `nf-test` l'ha creato.

### 4.2. Specificare i parametri di input

Dobbiamo ancora specificare gli input, il che viene fatto in modo leggermente diverso a livello di workflow rispetto ai test a livello di modulo.
Ci sono diversi modi per farlo, incluso specificando un profilo.
Tuttavia, un modo più semplice è configurare un blocco `params {}` nel file `nextflow.config` che `nf-test init` ha originariamente creato nella directory `tests`.

```groovy title="tests/nextflow.config" linenums="1"
/*
========================================================================================
    Nextflow config file for running tests
========================================================================================
*/

// Output directory for workflow outputs
outputDir = 'results_genomics'

/*
 * Pipeline parameters
 */

params {
    // Input primario (file di file di input, uno per riga)
    reads_bam = "${projectDir}/data/sample_bams.txt"

    // Accessory files
    reference = "${projectDir}/data/ref/ref.fasta"
    reference_index = "${projectDir}/data/ref/ref.fasta.fai"
    reference_dict = "${projectDir}/data/ref/ref.dict"
    intervals = "${projectDir}/data/ref/intervals.bed"

    // Base name for final output file
    cohort_name = "family_trio"
}
```

Quando eseguiamo il test, `nf-test` preleverà questo file di configurazione e caricherà gli input di conseguenza.

### 4.3. Eseguire il test del workflow

```bash
nf-test test tests/genomics-4.nf.test
```

??? success "Output del comando"

    ```console
    🚀 nf-test 0.9.3
    https://www.nf-test.com
    (c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


    Test Workflow genomics-4.nf

      Test [1b4c6936] 'Should run the pipeline without failures' PASSED (171.019s)


    SUCCESS: Executed 1 tests in 171.056s
    ```

Il test passa, confermando che la nostra pipeline completa di variant calling:

1. Processa con successo tutti i campioni
2. Concatena correttamente tutti i passaggi

### 4.4. Eseguire TUTTI i test

nf-test ha un altro asso nella manica. Possiamo eseguire tutti i test contemporaneamente! Modifichi il file `nf-test.config` in modo che nf-test cerchi in ogni directory i file nf-test. Può farlo modificando il parametro `testsDir`:

=== "Dopo"

    ```groovy title="nf-test.config" linenums="1" hl_lines="3"
    config {

        testsDir "."
        workDir ".nf-test"
        configFile "tests/nextflow.config"
        profile ""

    }
    ```

=== "Prima"

    ```groovy title="nf-test.config" linenums="1" hl_lines="3"
    config {

        testsDir "tests"
        workDir ".nf-test"
        configFile "tests/nextflow.config"
        profile ""

    }
    ```

Ora, possiamo semplicemente eseguire nf-test e eseguirà _ogni singolo test_ nel nostro repository:

```bash
nf-test test
```

??? success "Output del comando"

    ```console
    🚀 nf-test 0.9.3
    https://www.nf-test.com
    (c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


    Test Process GATK_HAPLOTYPECALLER

      Test [c5156c2b] 'Should call son's haplotype correctly' PASSED (39.947s)
      Test [10de94a8] 'Should call mother's haplotype correctly' PASSED (43.17s)
      Test [c0386fc7] 'Should call father's haplotype correctly' PASSED (44.244s)

    Test Process GATK_JOINTGENOTYPING

      Test [ac2067de] 'Should call trio's joint genotype correctly' PASSED (61.129s)

    Test Process SAMTOOLS_INDEX

      Test [625e39ee] 'Should index reads_son.bam correctly' PASSED (8.671s)
      Test [a8b28f36] 'Should index reads_mother.bam correctly' PASSED (8.518s)
      Test [c15852a1] 'Should index reads_father.bam correctly' PASSED (5.378s)

    Test Workflow genomics-4.nf

      Test [1b4c6936] 'Should run the pipeline without failures' PASSED (169.714s)


    SUCCESS: Executed 8 tests in 380.801s
    ```

8 test in 1 comando! Abbiamo speso molto tempo configurando molti test, ma quando è arrivato il momento di eseguirli è stato molto rapido e facile. Può vedere quanto sia utile quando si mantiene una grande pipeline, che potrebbe includere centinaia di elementi diversi. Spendiamo tempo scrivendo i test una volta così possiamo risparmiare tempo eseguendoli molte volte.

Inoltre, possiamo automatizzare questo! Immagini che i test vengano eseguiti ogni volta che voi o un collega cerca di aggiungere nuovo codice. Questo è il modo in cui garantiamo che le nostre pipeline mantengano uno standard elevato.

## Takeaway

Ora sa come scrivere ed eseguire diversi tipi di test per la vostra pipeline genomics utilizzando nf-test. Questo framework di testing aiuta a garantire che il vostro workflow di variant calling produca risultati coerenti e affidabili in diversi ambienti e mentre apporta modifiche al codice.

Ha imparato a testare componenti critici come:

- Il processo `SAMTOOLS_INDEX` che prepara i file BAM per il variant calling
- Il processo `GATK_HAPLOTYPECALLER` che identifica varianti nei singoli campioni
- Il processo `GATK_JOINTGENOTYPING` che combina le chiamate di varianti attraverso una coorte

Ha anche implementato diverse strategie di testing specifiche per i dati genomici:

- Verificare che i file VCF contengano le chiamate di varianti attese nonostante elementi non deterministici come i timestamp
- Testare con un dataset di trio familiare per garantire la corretta identificazione delle varianti tra campioni correlati
- Verificare coordinate genomiche specifiche e informazioni sulle varianti nei vostri file di output

Queste competenze di testing sono essenziali per sviluppare pipeline bioinformatiche robuste che possono processare in modo affidabile dati genomici e produrre chiamate di varianti accurate. Man mano che continua a lavorare con Nextflow per l'analisi genomica, questa base di testing La aiuterà a mantenere codice di alta qualità che produce risultati scientifici affidabili.
