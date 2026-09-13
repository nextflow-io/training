# Parte 3: Usare i profili per cambiare configurazione

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduzione assistita da IA - [scopri di più e suggerisci miglioramenti](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Nelle [Parte 1](./01_packaging_and_execution.md) e [Parte 2](./02_resources_and_retries.md), avete accumulato alcune opzioni di configurazione: packaging del software, piattaforma di esecuzione e allocazione delle risorse.
In pratica, spesso vorrete passare da un intero insieme di queste opzioni a un altro a seconda di dove state eseguendo il flusso di lavoro, ad esempio un laptop per lo sviluppo e un cluster HPC per la produzione.

Nextflow vi permette di definire un numero qualsiasi di [profili](https://nextflow.io/docs/latest/config.html#profiles) che descrivono configurazioni diverse, e di selezionarne uno (o più) al momento dell'esecuzione con un singolo flag.

Ne avete già usato uno: il profilo `test` di [Nextflow Run](../nextflow_run/index.md) sovrascrive i parametri di input con un insieme piccolo e ben definito.
Ora creerete i vostri profili di infrastruttura e li combinerete con esso.

---

## 1. Creare profili per ambienti diversi

### 1.1. Configurare i profili

Aggiungete due profili a `nextflow.config`: uno per eseguire su un normale laptop con Docker, e uno per un cluster HPC universitario con uno scheduler Slurm e Conda.

=== "Dopo"

    ```groovy title="nextflow.config" linenums="35" hl_lines="10-19"
    profiles {
        test {
            params.input = 'data/greetings.csv'
            params.batch = 'test'
            params.character = 'tux'
        }
        conda {
            docker.enabled = false
            conda.enabled = true
        }
        my_laptop {
            process.executor = 'local'
            docker.enabled = true
        }
        univ_hpc {
            process.executor = 'slurm'
            conda.enabled = true
            process.resourceLimits = [
                memory: 750.GB,
                cpus: 200,
                time: 30.d
            ]
        }
    }
    ```

=== "Prima"

    ```groovy title="nextflow.config" linenums="35"
    profiles {
        test {
            params.input = 'data/greetings.csv'
            params.batch = 'test'
            params.character = 'tux'
        }
        conda {
            docker.enabled = false
            conda.enabled = true
        }
    }
    ```

Il profilo `univ_hpc` imposta anche i limiti delle risorse, poiché ciò è tipicamente richiesto su infrastrutture HPC condivise.

### 1.2. Eseguire il flusso di lavoro con un profilo

Selezionate un profilo al momento dell'esecuzione con `-profile`.

```bash
nextflow run main.nf -profile my_laptop
```

??? success "Output del comando"

    ```console
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [nice_kimura] revision: c3c85dec78

    executor >  local (8)
    [df/b86b2f] sayHello (2)       | 3 of 3 ✔
    [b2/e2cb9c] convertToUpper (3) | 3 of 3 ✔
    [af/2ecff7] collectGreetings   | 1 of 1 ✔
    [c6/713c96] cowpy              | 1 of 1 ✔

    Outputs:

      /workspaces/training/execution-config/results

      first_output:
        - full_pipeline/intermediates/Bonjour-output.txt
        - full_pipeline/intermediates/Hello-output.txt
        - full_pipeline/intermediates/Hola-output.txt

      uppercased:
        - full_pipeline/intermediates/UPPER-Bonjour-output.txt
        - full_pipeline/intermediates/UPPER-Hello-output.txt
        - full_pipeline/intermediates/UPPER-Hola-output.txt

      collected: full_pipeline/intermediates/COLLECTED-batch-output.txt

      batch_report: full_pipeline/batch-report.txt

      cowpy_art: full_pipeline/cowpy-COLLECTED-batch-output.txt
    ```

!!! warning "Avviso"

    Il profilo `univ_hpc` non funzionerà nell'ambiente di formazione, poiché non è disponibile uno scheduler Slurm.

Se trovate altre impostazioni che appartengono sempre insieme, aggiungetele al profilo corrispondente.
Potete anche creare profili aggiuntivi per raggruppare qualsiasi altra combinazione di cui avete bisogno.

### 1.3. Eseguire con più profili

I profili non si escludono a vicenda.
Potete attivarne più di uno contemporaneamente con `-profile <profilo1>,<profilo2>`.
Combinate `my_laptop` con il profilo `test` che già conoscete da Nextflow Run.

```bash
nextflow run main.nf -profile my_laptop,test
```

??? success "Output del comando"

    ```console
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [romantic_volhard] revision: c3c85dec78

    executor >  local (8)
    [22/81de4f] sayHello (3)       | 3 of 3 ✔
    [93/e16b04] convertToUpper (1) | 3 of 3 ✔
    [ca/3e6a91] collectGreetings   | 1 of 1 ✔
    [16/0cc7e3] cowpy              | 1 of 1 ✔

    Outputs:

      /workspaces/training/execution-config/results

      first_output:
        - full_pipeline/intermediates/Hola-output.txt
        - full_pipeline/intermediates/Hello-output.txt
        - full_pipeline/intermediates/Bonjour-output.txt

      uppercased:
        - full_pipeline/intermediates/UPPER-Hello-output.txt
        - full_pipeline/intermediates/UPPER-Bonjour-output.txt
        - full_pipeline/intermediates/UPPER-Hola-output.txt

      collected: full_pipeline/intermediates/COLLECTED-test-output.txt

      batch_report: full_pipeline/test-report.txt

      cowpy_art: full_pipeline/cowpy-COLLECTED-test-output.txt
    ```

I nomi dei singoli file riflettono correttamente `batch = 'test'` dal profilo `test` (`COLLECTED-test-output.txt`, e così via).

Se combinate profili che impostano la stessa opzione, Nextflow risolve il conflitto usando l'ultimo valore letto, ovvero quello che compare più avanti nel file.
Se le impostazioni in conflitto provengono da sorgenti di configurazione completamente diverse, si applica il [ordine di precedenza](https://www.nextflow.io/docs/latest/config.html) standard.

### Takeaway

Sapete come definire profili che raggruppano configurazioni specifiche per l'infrastruttura, selezionarne uno al momento dell'esecuzione con `-profile`, combinare più profili in una singola esecuzione, e come Nextflow risolve i conflitti quando più di un profilo imposta la stessa opzione.

### Cosa c'è dopo?

Imparate come ispezionare la configurazione completamente risolta prima di eseguire qualsiasi cosa.

---

## 2. Ispezionare la configurazione risolta

Avete già usato `nextflow config -profile test` in [Nextflow Run](../nextflow_run/02_configure_pipeline.md) per verificare a cosa si risolve un singolo profilo.
Quel comando diventa particolarmente utile quando si combinano più profili: come avete appena visto, quando due profili impostano la stessa opzione, può essere difficile determinare manualmente quale valore prevalga effettivamente.
Il comando `nextflow config` risolve tutto questo per voi, senza eseguire la pipeline.

### 2.1. Risolvere la configurazione predefinita

```bash
nextflow config
```

??? success "Output del comando"

    ```groovy
    params {
       input = 'data/greetings.csv'
       batch = 'batch'
       character = 'turkey'
    }

    docker {
       enabled = true
    }

    process {
       memory = '1 GB'
       withName:cowpy {
          conda = 'conda-forge::cowpy==1.1.5'
          memory = '2 GB'
          cpus = 2
       }
    }
    ```

Questa è esattamente la configurazione che verrebbe applicata se eseguiste la pipeline senza flag aggiuntivi.

### 2.2. Risolvere la configurazione con i profili attivati

Aggiungete gli stessi profili che usereste per un'esecuzione reale.

```bash
nextflow config -profile my_laptop,test
```

??? success "Output del comando"

    ```groovy
    params {
       input = 'data/greetings.csv'
       batch = 'test'
       character = 'tux'
    }

    docker {
       enabled = true
    }

    process {
       memory = '1 GB'
       withName:cowpy {
          conda = 'conda-forge::cowpy==1.1.5'
          memory = '2 GB'
          cpus = 2
       }
       executor = 'local'
    }
    ```

Confrontando i due output si conferma cosa è cambiato: `params.batch`, `params.character` e `process.executor` riflettono tutti i profili `my_laptop,test`.
Questo diventa particolarmente prezioso per le pipeline con molti livelli di configurazione, dove determinare manualmente le impostazioni risolte sarebbe tedioso e soggetto a errori.

### Takeaway

Sapete come usare `nextflow config` per ispezionare la configurazione completamente risolta per qualsiasi combinazione di profili, prima di eseguire qualsiasi cosa.

### Cosa c'è dopo?

Avete coperto gli elementi essenziali della configurazione delle pipeline Nextflow.
Consultate il [Riepilogo del corso](next_steps.md) per sapere come proseguire da qui.

---

## Riepilogo

In questa parte avete imparato a:

- Definire profili che raggruppano configurazioni specifiche per l'infrastruttura
- Combinare più profili in una singola esecuzione e capire come vengono risolti i conflitti tra di essi
- Usare `nextflow config` per ispezionare la configurazione completamente risolta
