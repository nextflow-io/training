# Test con nf-test

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduzione assistita da IA - [scopri di più e suggerisci miglioramenti](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Essere in grado di verificare sistematicamente che ogni parte del flusso di lavoro faccia ciò che dovrebbe fare è fondamentale per la riproducibilità e la manutenzione a lungo termine, e può essere di grande aiuto durante il processo di sviluppo.

Prendiamoci un momento per parlare di perché i test sono così importanti. Se state sviluppando un flusso di lavoro, una delle prime cose che farete è prendere alcuni dati di test che sapete essere validi e che dovrebbero produrre un risultato. Aggiungete il primo processo alla pipeline e lo collegate ai vostri input per farlo funzionare. Poi, per verificare che tutto funzioni, lo eseguite sui dati di test. Supponendo che funzioni, passate al processo successivo ed eseguite di nuovo i dati di test. Ripetete questo processo finché non avete una pipeline con cui siete soddisfatti.

Poi, magari aggiungete un semplice parametro vero o falso come `--skip_process`. Ora dovete eseguire la pipeline due volte, una per ogni parametro, per assicurarvi che funzioni come previsto. Ma aspettate, come facciamo a verificare se `--skip_process` salta effettivamente il processo? Dobbiamo esaminare gli output o controllare i file di log! Questo è fastidioso ed è soggetto a errori.

Man mano che sviluppate la vostra pipeline, diventerà rapidamente così complessa che testare manualmente ogni iterazione è lento e soggetto a errori. Inoltre, se trovate un errore, sarà molto difficile individuare esattamente dove nella pipeline si trova. È qui che entrano in gioco i test.

I test vi permettono di verificare sistematicamente che ogni parte della pipeline funzioni come previsto. I vantaggi per uno sviluppatore di test ben scritti sono enormi:

- **Fiducia**: Poiché i test coprono l'intera pipeline, potete essere certi che modificare qualcosa non influenzi nient'altro
- **Affidabilità**: Quando più sviluppatori lavorano sulla pipeline, sanno che gli altri sviluppatori non hanno compromesso la pipeline e ogni componente.
- **Trasparenza**: I test mostrano dove una pipeline sta fallendo e rendono più facile individuare il problema. Funzionano anche come forma di documentazione, mostrando come eseguire un processo o un flusso di lavoro.
- **Velocità**: Poiché i test sono automatizzati, possono essere eseguiti molto rapidamente e ripetutamente. Potete iterare rapidamente con meno paura di introdurre nuovi bug.

Ci sono molti tipi diversi di test che possiamo scrivere:

1. **Test a livello di modulo**: Per i singoli processi
2. **Test a livello di workflow**: Per un singolo flusso di lavoro
3. **Test a livello di pipeline**: Per la pipeline nel suo insieme
4. **Test di performance**: Per la velocità e l'efficienza della pipeline
5. **Stress test**: Per valutare le prestazioni della pipeline in condizioni estreme e determinarne i limiti

Testare i singoli processi è analogo agli unit test in altri linguaggi. Testare il flusso di lavoro o l'intera pipeline è analogo a quelli che vengono chiamati integration test in altri linguaggi, dove testiamo le interazioni tra i componenti.

[**nf-test**](https://www.nf-test.com/) è uno strumento che vi permette di scrivere test a livello di modulo, workflow e pipeline. In breve, vi permette di verificare sistematicamente che ogni singola parte della pipeline funzioni come previsto, _in isolamento_.

### Obiettivi di apprendimento

In questa side quest, imparerete a usare nf-test per scrivere un test a livello di workflow per la pipeline, nonché test a livello di modulo per i tre processi che essa richiama.

Al termine di questa side quest, sarete in grado di usare efficacemente le seguenti tecniche:

- Inizializzare nf-test nel vostro progetto
- Generare test a livello di modulo e di workflow
- Aggiungere tipi comuni di asserzioni
- Capire quando usare gli snapshot rispetto alle asserzioni sul contenuto
- Eseguire test per un intero progetto

Queste competenze vi aiuteranno a implementare una strategia di test completa nei vostri progetti di pipeline, rendendoli più robusti e manutenibili.

### Prerequisiti

Prima di affrontare questa side quest, dovreste:

- Aver completato il tutorial [Hello Nextflow](../hello_nextflow/README.md) o un corso equivalente per principianti.
- Essere a proprio agio con i concetti e i meccanismi di base di Nextflow (processi, canali, operatori, gestione di file e metadati)

---

## 0. Iniziamo

#### Aprite il codespace di formazione

Se non lo avete ancora fatto, assicuratevi di aprire l'ambiente di formazione come descritto nella sezione [Configurazione dell'ambiente](../envsetup/index.md).

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### Spostatevi nella directory del progetto

Spostiamoci nella directory dove si trovano i file per questo tutorial.

```bash
cd side-quests/nf-test
```

Potete impostare VSCode per concentrarsi su questa directory:

```bash
code .
```

#### Esaminate i materiali

Troverete un file del flusso di lavoro principale e un file CSV chiamato `greetings.csv` che contiene l'input della pipeline.

```console title="Directory contents"
.
├── greetings.csv
└── main.nf
```

Per una descrizione dettagliata dei file, consultate il [riscaldamento di Hello Nextflow](../hello_nextflow/00_orientation.md).

Il flusso di lavoro che testeremo è un sottoinsieme del flusso di lavoro Hello costruito in [Hello Workflow](../hello_nextflow/03_hello_workflow.md).

??? example "Cosa fa il flusso di lavoro Hello Nextflow?"

    Se non avete seguito la formazione [Hello Nextflow](../hello_nextflow/index.md), ecco una rapida panoramica di ciò che fa questo semplice flusso di lavoro.

    Il flusso di lavoro prende un file CSV contenente saluti, esegue quattro passaggi di trasformazione consecutivi su di essi e produce un singolo file di testo contenente un'immagine ASCII di un personaggio divertente che dice i saluti.

    I quattro passaggi sono implementati come processi Nextflow (`sayHello`, `convertToUpper`, `collectGreetings` e `cowpy`) memorizzati in file di modulo separati.

    1. **`sayHello`:** Scrive ogni saluto nel proprio file di output (es. "Hello-output.txt")
    2. **`convertToUpper`:** Converte ogni saluto in maiuscolo (es. "HELLO")
    3. **`collectGreetings`:** Raccoglie tutti i saluti in maiuscolo in un singolo file batch
    4. **`cowpy`:** Genera arte ASCII usando lo strumento `cowpy`

    I risultati vengono pubblicati in una directory chiamata `results/`, e l'output finale della pipeline (quando eseguita con i parametri predefiniti) è un file di testo normale contenente arte ASCII di un personaggio che dice i saluti in maiuscolo.

    In questa side quest, utilizziamo una forma intermedia del flusso di lavoro Hello che contiene solo i primi due processi. <!-- TODO: change this to use the full finished workflow as suggested in https://github.com/nextflow-io/training/issues/735 -->

Il sottoinsieme con cui lavoreremo è composto da due processi: `sayHello` e `convertToUpper`.
Potete vedere il codice completo del flusso di lavoro qui sotto.

??? example "Codice del workflow"

    ```groovy title="main.nf"
    /*
    * Parametri della pipeline
    */
    params.input_file = "greetings.csv"

    /*
    * Usa echo per stampare 'Hello World!' sullo standard output
    */
    process sayHello {

        publishDir 'results', mode: 'copy'

        input:
            val greeting

        output:
            path "${greeting}-output.txt"

        script:
        """
        echo '$greeting' > '$greeting-output.txt'
        """
    }

    /*
    * Usa un'utilità di sostituzione testo per convertire il saluto in maiuscolo
    */
    process convertToUpper {

        publishDir 'results', mode: 'copy'

        input:
            path input_file

        output:
            path "UPPER-${input_file}"

        script:
        """
        cat '$input_file' | tr '[a-z]' '[A-Z]' > UPPER-${input_file}
        """
    }

    workflow {

        // crea un canale per gli input da un file CSV
        greeting_ch = channel.fromPath(params.input_file).splitCsv().flatten()

        // emette un saluto
        sayHello(greeting_ch)

        // converte il saluto in maiuscolo
        convertToUpper(sayHello.out)
    }
    ```

#### Eseguite il flusso di lavoro

Eseguiamo il flusso di lavoro per assicurarci che funzioni come previsto.

```bash
nextflow run main.nf
```

```console title="Result of running the workflow"
 N E X T F L O W   ~  version 24.10.2

Launching `main.nf` [soggy_linnaeus] DSL2 - revision: bbf79d5c31

executor >  local (6)
[f7/c3be66] sayHello (3)       | 3 of 3 ✔
[cd/e15303] convertToUpper (3) | 3 of 3 ✔
```

CONGRATULAZIONI! Avete appena eseguito un test!

"Aspettate, cosa? Ho appena eseguito il flusso di lavoro e ha funzionato! Come può essere un test?"

Ottima domanda!

Analizziamo quello che è appena successo.

Avete eseguito il flusso di lavoro con i parametri predefiniti, avete confermato che ha funzionato e siete soddisfatti dei risultati. Questa è l'essenza del testing. Se avete seguito il corso di formazione Hello Nextflow, avrete notato che abbiamo sempre iniziato ogni sezione eseguendo il flusso di lavoro che stavamo usando come punto di partenza, per confermare che tutto fosse configurato correttamente.

Il testing del software fa essenzialmente questo processo per noi.

#### Esaminate il compito

La vostra sfida è aggiungere test standardizzati a questo flusso di lavoro usando nf-test, in modo da rendere facile verificare che ogni parte continui a funzionare come previsto nel caso vengano apportate ulteriori modifiche.

<!-- TODO: give a bit more details, similar to how it's done in the Metadata side quest -->

#### Lista di controllo per la preparazione

Pensate di essere pronti a iniziare?

- [ ] Capisco l'obiettivo di questo corso e i suoi prerequisiti
- [ ] Il mio codespace è attivo e funzionante
- [ ] Ho impostato correttamente la mia directory di lavoro
- [ ] Ho eseguito il flusso di lavoro con successo
- [ ] Capisco il compito

Se potete spuntare tutte le caselle, siete pronti per partire.

---

## 1. Inizializzare `nf-test`

Il pacchetto `nf-test` fornisce un comando di inizializzazione che configura alcune cose per permetterci di iniziare a sviluppare test per il nostro progetto.

```bash
nf-test init
```

Questo dovrebbe produrre il seguente output:

```bash
🚀 nf-test 0.9.3
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr

Project configured. Configuration is stored in nf-test.config
```

Crea anche una directory `tests` contenente uno stub del file di configurazione.

### 1.1. Generare un nf-test

`nf-test` viene fornito con un insieme di strumenti per costruire file nf-test, risparmiandoci la maggior parte del lavoro. Questi si trovano sotto il sottocomando `generate`. Generiamo un test per la pipeline:

```bash
nf-test generate pipeline main.nf
```

```console title="Output"
> nf-test generate pipeline main.nf

Load source file '/workspaces/training/side-quests/nf-test/main.nf'
Wrote pipeline test file '/workspaces/training/side-quests/nf-test/tests/main.nf.test

SUCCESS: Generated 1 test files.
```

Questo creerà un file `main.nf.test` all'interno della directory `tests`. Questo è il nostro file di test a livello di pipeline. Se eseguite `tree tests/` dovreste vedere qualcosa di simile:

```console title="Test directory contents"
tests/
├── main.nf.test
└── nextflow.config
```

Il file `main.nf.test` è il nostro file di test a livello di pipeline. Apriamolo e diamo un'occhiata al contenuto.

```groovy title="tests/main.nf.test"
nextflow_pipeline {

    name "Test Workflow main.nf"
    script "main.nf"

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

Prendiamoci un momento per capire la struttura del file di test.

Il blocco `nextflow_pipeline` è il punto di ingresso per tutti i test a livello di pipeline. Contiene quanto segue:

- `name`: Il nome del test.
- `script`: Il percorso allo script della pipeline.

Il blocco `test` è il test vero e proprio. Contiene quanto segue:

- `when`: Le condizioni in cui il test deve essere eseguito. Questo include i parametri che verranno usati per eseguire la pipeline.
- `then`: Le asserzioni che devono essere fatte. Questo include i risultati attesi della pipeline.

In italiano semplice, la logica del test si legge come segue:
"**Quando** questi _parametri_ vengono forniti a questa _pipeline_, **allora** ci aspettiamo di vedere questi risultati."

Questo non è un test funzionale; dimostreremo come trasformarlo in uno nella prossima sezione.

### Una nota sui nomi dei test

Nell'esempio precedente, abbiamo usato il nome predefinito "Should run without failures" che è appropriato per un test di base che verifica solo se la pipeline viene eseguita con successo. Tuttavia, man mano che aggiungiamo casi di test più specifici, dovremmo usare nomi più descrittivi che indichino cosa stiamo effettivamente testando. Per esempio:

- "Should convert input to uppercase" - quando si testa una funzionalità specifica
- "Should handle empty input gracefully" - quando si testano casi limite
- "Should respect max memory parameter" - quando si testano i vincoli sulle risorse
- "Should create expected output files" - quando si testa la generazione di file

I buoni nomi dei test dovrebbero:

1. Iniziare con "Should" per chiarire qual è il comportamento atteso
2. Descrivere la funzionalità specifica o lo scenario che viene testato
3. Essere sufficientemente chiari che, se il test fallisce, si sappia quale funzionalità è compromessa

Man mano che aggiungiamo più asserzioni e casi di test specifici in seguito, useremo questi nomi più descrittivi per chiarire cosa verifica ogni test.

### 1.2. Eseguire il test

Eseguiamo il test per vedere cosa succede.

```bash
nf-test test tests/main.nf.test
```

```console title="nf-test pipeline fail"
> nf-test test tests/main.nf.test

🚀 nf-test 0.9.3
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Workflow main.nf

  Test [693ba951] 'Should run without failures' FAILED (4.652s)

  Assertion failed:

  assert workflow.success
         |        |
         workflow false

  Nextflow stdout:

  ERROR ~ No such file or directory: /workspaces/training/side-quests/nf-test/.nf-test/tests/693ba951a20fec36a5a9292ed1cc8a9f/greetings.csv

   -- Check '/workspaces/training/side-quests/nf-test/.nf-test/tests/693ba951a20fec36a5a9292ed1cc8a9f/meta/nextflow.log' file for details
  Nextflow stderr:

FAILURE: Executed 1 tests in 4.679s (1 failed)
```

Il test fallisce! Cosa è successo?

1. nf-test ha provato a eseguire la pipeline così com'è, usando le impostazioni nel blocco `when`:

```groovy title="tests/main.nf.test"
when {
    params {
        // define parameters here. Example:
        // outdir = "tests/results"
    }
}
```

2. nf-test ha verificato lo stato della pipeline e lo ha confrontato con il blocco `when`:

```groovy title="tests/main.nf.test"
then {
    assert workflow.success
}
```

Notate come nf-test ha segnalato che la pipeline è fallita e ha fornito il messaggio di errore di Nextflow:

```console title="Error"
ERROR ~ No such file or directory: /workspaces/training/side-quests/nf-test/.nf-test/tests/693ba951a20fec36a5a9292ed1cc8a9f/greetings.csv
```

Qual era il problema? Ricordate che la pipeline ha un file greetings.csv nella directory del progetto. Quando nf-test esegue la pipeline, cercherà questo file, ma non riesce a trovarlo. Il file è lì, cosa sta succedendo? Beh, se guardiamo il percorso possiamo vedere che il test si sta svolgendo nel percorso `./nf-test/tests/longHashString/`. Proprio come Nextflow, nf-test crea una nuova directory per ogni test per mantenere tutto isolato. Il file di dati non si trova lì, quindi dobbiamo correggere il percorso al file nel test originale.

Torniamo al file di test e cambiamo il percorso al file nel blocco `when`.

Potreste chiedervi come faremo a puntare alla radice della pipeline nel test. Poiché questa è una situazione comune, nf-test ha una serie di variabili globali che possiamo usare per semplificarci la vita. Potete trovare l'elenco completo [qui](https://www.nf-test.com/docs/testcases/global_variables/) ma nel frattempo useremo la variabile `projectDir`, che indica la radice del progetto della pipeline.

=== "Dopo"

    ```groovy title="tests/main.nf.test" linenums="1" hl_lines="3"
    when {
        params {
            input_file = "${projectDir}/greetings.csv"
        }
    }
    ```

=== "Prima"

    ```groovy title="tests/main.nf.test" linenums="1" hl_lines="3 4"
    when {
        params {
            // define parameters here. Example:
            // outdir = "tests/results"
        }
    }
    ```

Eseguiamo di nuovo il test per vedere se funziona.

```bash title="nf-test pipeline pass"
nf-test test tests/main.nf.test
```

```console title="Pipeline passes"
> nf-test test tests/main.nf.test

🚀 nf-test 0.9.3
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Workflow main.nf

  Test [1d4aaf12] 'Should run without failures' PASSED (1.619s)


SUCCESS: Executed 1 tests in 1.626s
```

Successo! La pipeline viene eseguita con successo e il test passa. Eseguitelo tutte le volte che volete e otterrete sempre lo stesso risultato!

Per impostazione predefinita, l'output di Nextflow è nascosto, ma per convincervi che nf-test sta effettivamente eseguendo il flusso di lavoro, potete usare il flag `--verbose`:

```bash
nf-test test tests/main.nf.test --verbose
```

```console title="Pipeline runs all processes"
> nf-test test tests/main.nf.test

🚀 nf-test 0.9.3
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Workflow main.nf

  Test [693ba951] 'Should run without failures'
    > Nextflow 24.10.4 is available - Please consider updating your version to it
    > N E X T F L O W  ~  version 24.10.0
    > Launching `/workspaces/training/side-quests/nf-test/main.nf` [zen_ampere] DSL2 - revision: bbf79d5c31
    > [2b/61e453] Submitted process > sayHello (2)
    > [31/4e1606] Submitted process > sayHello (1)
    > [bb/5209ee] Submitted process > sayHello (3)
    > [83/83db6f] Submitted process > convertToUpper (2)
    > [9b/3428b1] Submitted process > convertToUpper (1)
    > [ca/0ba51b] Submitted process > convertToUpper (3)
    PASSED (5.206s)


SUCCESS: Executed 1 tests in 5.239s
```

### 1.3. Aggiungere asserzioni

Un semplice controllo è assicurarsi che la nostra pipeline stia eseguendo tutti i processi che ci aspettiamo e non ne stia saltando nessuno silenziosamente. Ricordate che la nostra pipeline esegue 6 processi, uno chiamato `sayHello` e uno chiamato `convertToUpper` per ciascuno dei 3 saluti.

Aggiungiamo un'asserzione al nostro test per verificare che la pipeline esegua il numero atteso di processi. Aggiorneremo anche il nome del test per riflettere meglio ciò che stiamo testando.

=== "Dopo"

    ```groovy title="tests/main.nf.test" linenums="1" hl_lines="1 11"
        test("Should run successfully with correct number of processes") {

            when {
                params {
                    input_file = "${projectDir}/greetings.csv"
                }
            }

            then {
                assert workflow.success
                assert workflow.trace.tasks().size() == 6
            }

        }
    ```

=== "Prima"

    ```groovy title="tests/main.nf.test" linenums="1" hl_lines="1"
        test("Should run without failures") {

            when {
                params {
                    input_file = "${projectDir}/greetings.csv"
                }
            }

            then {
                assert workflow.success
            }

        }
    ```

Il nome del test ora riflette meglio ciò che stiamo effettivamente verificando: non solo che la pipeline viene eseguita senza fallire, ma che esegue il numero atteso di processi.

Eseguiamo di nuovo il test per vedere se funziona.

```bash title="nf-test pipeline pass"
nf-test test tests/main.nf.test
```

```console title="Pipeline passes with assertions"
🚀 nf-test 0.9.3
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Workflow main.nf

  Test [1d4aaf12] 'Should run successfully with correct number of processes' PASSED (1.567s)


SUCCESS: Executed 1 tests in 1.588s
```

Successo! La pipeline viene eseguita con successo e il test passa. Ora abbiamo iniziato a testare i dettagli della pipeline, oltre allo stato generale.

### 1.4. Testare l'output

Aggiungiamo un'asserzione al nostro test per verificare che il file di output sia stato creato. Lo aggiungeremo come test separato, con un nome informativo, per rendere i risultati più facili da interpretare.

=== "Dopo"

    ```groovy title="tests/main.nf.test" linenums="1" hl_lines="14-33"
        test("Should run successfully with correct number of processes") {

            when {
                params {
                    input_file = "${projectDir}/greetings.csv"
                }
            }

            then {
                assert workflow.success
                assert workflow.trace.tasks().size() == 6
            }

        }

        test("Should produce correct output files") {

            when {
                params {
                    input_file = "${projectDir}/greetings.csv"
                }
            }

            then {
                assert file("$launchDir/results/Bonjour-output.txt").exists()
                assert file("$launchDir/results/Hello-output.txt").exists()
                assert file("$launchDir/results/Holà-output.txt").exists()
                assert file("$launchDir/results/UPPER-Bonjour-output.txt").exists()
                assert file("$launchDir/results/UPPER-Hello-output.txt").exists()
                assert file("$launchDir/results/UPPER-Holà-output.txt").exists()
            }

        }
    ```

=== "Prima"

    ```groovy title="tests/main.nf.test" linenums="1" hl_lines="14"
        test("Should run successfully with correct number of processes") {

            when {
                params {
                    input_file = "${projectDir}/greetings.csv"
                }
            }

            then {
                assert workflow.success
                assert workflow.trace.tasks().size() == 6
            }

        }
    ```

Eseguite di nuovo il test per vedere se funziona.

```bash title="nf-test pipeline pass"
nf-test test tests/main.nf.test
```

```console title="Pipeline passes with file assertions"
> nf-test test tests/main.nf.test

🚀 nf-test 0.9.3
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Workflow main.nf

  Test [f0e08a68] 'Should run successfully with correct number of processes' PASSED (8.144s)
  Test [d7e32a32] 'Should produce correct output files' PASSED (6.994s)


SUCCESS: Executed 2 tests in 15.165s
```

Successo! I test passano perché la pipeline è stata completata con successo, il numero corretto di processi è stato eseguito e i file di output sono stati creati. Questo dovrebbe anche mostrarvi quanto sia utile fornire quei nomi informativi per i vostri test.

Questa è solo la superficie; possiamo continuare a scrivere asserzioni per verificare i dettagli della pipeline, ma per ora passiamo al testing degli elementi interni della pipeline.

### Takeaway

Sapete come scrivere un nf-test per una pipeline.

### Cosa c'è dopo?

Imparate come testare un processo Nextflow.

---

## 2. Testare un processo Nextflow

Non dobbiamo scrivere test per ogni parte della pipeline, ma più test abbiamo, più possiamo essere completi riguardo alla pipeline e più possiamo essere certi che funzioni come previsto. In questa sezione testeremo entrambi i processi nella pipeline come unità individuali.

### 2.1. Testare il processo `sayHello`

Iniziamo con il processo `sayHello`.

Usiamo di nuovo il comando `nf-test generate` per generare test per il processo.

```bash
nf-test generate process main.nf
```

```console title="Output"
> nf-test generate process main.nf

Load source file '/workspaces/training/side-quests/nf-test/main.nf'
Wrote process test file '/workspaces/training/side-quests/nf-test/tests/main.sayhello.nf.test
Wrote process test file '/workspaces/training/side-quests/nf-test/tests/main.converttoupper.nf.test

SUCCESS: Generated 2 test files.
```

Concentriamoci per ora sul processo `sayhello` nel file `main.sayhello.nf.test`.

Apriamo il file e diamo un'occhiata al contenuto.

```groovy title="tests/main.sayhello.nf.test"
nextflow_process {

    name "Test Process sayHello"
    script "main.nf"
    process "sayHello"

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

Come prima, iniziamo con i dettagli del test, seguiti dai blocchi `when` e `then`. Tuttavia, abbiamo anche un blocco `process` aggiuntivo che ci permette di definire gli input del processo.

Eseguiamo il test per vedere se funziona.

```bash title="nf-test pipeline pass"
nf-test test tests/main.sayhello.nf.test
```

```console title="Process test fails"
> nf-test test tests/main.sayhello.nf.test

🚀 nf-test 0.9.3
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Process sayHello

  Test [1eaad118] 'Should run without failures' FAILED (4.876s)

  Assertion failed:

  assert process.success
         |       |
         |       false
         sayHello

  Nextflow stdout:

  Process `sayHello` declares 1 input but was called with 0 arguments
  Nextflow stderr:

FAILURE: Executed 1 tests in 4.884s (1 failed)
```

Il test fallisce perché il processo `sayHello` dichiara 1 input ma è stato chiamato con 0 argomenti. Correggiamo questo aggiungendo un input al processo. Ricordate da [Hello Workflow](../hello_nextflow/03_hello_workflow.md) (e dalla sezione di riscaldamento sopra) che il nostro processo `sayHello` prende un singolo input di tipo val, che dovremo fornire. Dovremmo anche correggere il nome del test per riflettere meglio ciò che stiamo testando.

=== "Dopo"

    ```groovy title="tests/main.sayhello.nf.test" linenums="1" hl_lines="1 10"
        test("Should run without failures and produce correct output") {

            when {
                params {
                    // define parameters here. Example:
                    // outdir = "tests/results"
                }
                process {
                    """
                    input[0] = "hello"
                    """
                }
            }

            then {
                assert process.success
                assert snapshot(process.out).match()
            }

        }
    ```

=== "Prima"

    ```groovy title="tests/main.sayhello.nf.test" linenums="1" hl_lines="1 10 11"
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
    ```

Eseguiamo di nuovo il test per vedere se funziona.

```console title="nf-test pipeline pass"
> nf-test test tests/main.sayhello.nf.test

🚀 nf-test 0.9.3
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Process sayHello

  Test [f91a1bcd] 'Should run without failures and produce correct output' PASSED (1.604s)
  Snapshots:
    1 created [Should run without failures and produce correct output]


Snapshot Summary:
  1 created

SUCCESS: Executed 1 tests in 1.611s
```

Successo! Il test passa perché il processo `sayHello` è stato eseguito con successo e l'output è stato creato.

### 2.2. Esaminare lo snapshot creato dal test

Se guardiamo il file `tests/main.sayhello.nf.test`, possiamo vedere che usa un metodo `snapshot()` nel blocco delle asserzioni:

```groovy title="tests/main.sayhello.nf.test"
assert snapshot(process.out).match()
```

Questo dice a nf-test di creare uno snapshot dell'output del processo `sayHello`. Diamo un'occhiata al contenuto del file snapshot.

```console title="Snapshot file contents"
code tests/main.sayhello.nf.test.snap
```

Non lo stamperemo qui, ma dovreste vedere un file JSON contenente i dettagli del processo e degli output del processo. In particolare, possiamo vedere una riga che assomiglia a questa:

```json title="Snapshot file contents"
"0": [
    "hello-output.txt:md5,b1946ac92492d2347c6235b4d2611184"
]
```

Questo rappresenta gli output creati dal processo `sayHello`, che stiamo testando esplicitamente. Se rieseguiamo il test, il programma verificherà che il nuovo output corrisponda all'output originariamente registrato. Questo è un modo rapido e semplice per testare che gli output del processo non cambino, motivo per cui nf-test lo fornisce come impostazione predefinita.

!!!warning "Attenzione"

    Ciò significa che dobbiamo essere sicuri che l'output che registriamo nell'esecuzione originale sia corretto!

Se, nel corso dello sviluppo futuro, qualcosa nel codice cambia causando un output diverso, il test fallirà e dovremo determinare se il cambiamento è atteso o meno.

- Se risulta che qualcosa nel codice si è rotto, dovremo correggerlo, con l'aspettativa che il codice corretto superi il test.
- Se è un cambiamento atteso (es. lo strumento è stato migliorato e i risultati sono migliori), dovremo aggiornare lo snapshot per accettare il nuovo output come riferimento da confrontare. nf-test ha un parametro `--update-snapshot` a questo scopo.

Possiamo eseguire di nuovo il test e vedere che dovrebbe passare:

```console title="nf-test process pass with snapshot"
> nf-test test tests/main.sayhello.nf.test

🚀 nf-test 0.9.3
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Process sayHello

  Test [f91a1bcd] 'Should run without failures and produce correct output' PASSED (1.675s)


SUCCESS: Executed 1 tests in 1.685s
```

Successo! Il test passa perché il processo `sayHello` è stato eseguito con successo e l'output corrisponde allo snapshot.

### 2.3. Alternativa agli snapshot: asserzioni dirette sul contenuto

Mentre gli snapshot sono ottimi per rilevare qualsiasi cambiamento nell'output, a volte si vuole verificare contenuto specifico senza essere così rigidi riguardo alla corrispondenza dell'intero file. Per esempio:

- Quando parti dell'output potrebbero cambiare (timestamp, ID casuali, ecc.) ma certi contenuti chiave devono essere presenti
- Quando si vuole verificare la presenza di pattern o valori specifici nell'output
- Quando si vuole rendere il test più esplicito riguardo a cosa costituisce il successo

Ecco come potremmo modificare il nostro test per verificare contenuto specifico:

=== "Dopo"

    ```groovy title="tests/main.sayhello.nf.test" linenums="1" hl_lines="1 5 16 17"
        test("Should run without failures and contain expected greeting") {

            when {
                params {
                    // define parameters here
                }
                process {
                    """
                    input[0] = "hello"
                    """
                }
            }

            then {
                assert process.success
                assert path(process.out[0][0]).readLines().contains('hello')
                assert !path(process.out[0][0]).readLines().contains('HELLO')
            }

        }
    ```

=== "Prima"

    ```groovy title="tests/main.sayhello.nf.test" linenums="1" hl_lines="1 5 6 17"
        test("Should run without failures and produce correct output") {

            when {
                params {
                    // define parameters here. Example:
                    // outdir = "tests/results"
                }
                process {
                    """
                    input[0] = "hello"
                    """
                }
            }

            then {
                assert process.success
                assert snapshot(process.out).match()
            }

        }
    ```

Notate che nf-test vede gli output del processo come una lista di liste, quindi `process.out[0][0]` sta recuperando la prima parte del primo elemento del canale (o 'emissione') da questo processo.

Questo approccio:

- Rende chiaro esattamente cosa ci aspettiamo nell'output
- È più resiliente ai cambiamenti irrilevanti nell'output
- Fornisce messaggi di errore migliori quando i test falliscono
- Permette validazioni più complesse (pattern regex, confronti numerici, ecc.)

Eseguiamo il test per vedere se funziona.

```bash title="nf-test pipeline pass"
nf-test test tests/main.sayhello.nf.test
```

```console title="Process test fails"
> nf-test test tests/main.sayhello.nf.test

🚀 nf-test 0.9.3
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Process sayHello

  Test [58df4e4b] 'Should run without failures and contain expected greeting' PASSED (7.196s)


SUCCESS: Executed 1 tests in 7.208s
```

### 2.4. Testare il processo `convertToUpper`

Apriamo il file `tests/main.converttoupper.nf.test` e diamo un'occhiata al contenuto:

```groovy title="tests/main.converttoupper.nf.test"
nextflow_process {

    name "Test Process convertToUpper"
    script "main.nf"
    process "convertToUpper"

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

Questo è un test simile a quello del processo `sayHello`, ma sta testando il processo `convertToUpper`. Sappiamo che questo fallirà perché, proprio come con `sayHello`, il processo `convertToUpper` prende un singolo input di tipo path, ma non ne abbiamo specificato uno.

Ora dobbiamo fornire un singolo file di input al processo convertToUpper, che include del testo che vogliamo convertire in maiuscolo. Ci sono molti modi in cui potremmo farlo:

- Potremmo creare un file dedicato per il test
- Potremmo riutilizzare il file data/greetings.csv esistente
- Potremmo crearlo al volo all'interno del test

Per ora, riutilizziamo il file data/greetings.csv esistente usando l'esempio che abbiamo usato con il test a livello di pipeline. Come prima, possiamo nominare il test per riflettere meglio ciò che stiamo testando, ma questa volta lasciamo che 'snapshot' il contenuto piuttosto che verificare stringhe specifiche (come abbiamo fatto nell'altro processo).

=== "Dopo"

    ```groovy title="tests/main.converttoupper.nf.test" linenums="1" hl_lines="1 10"
        test("Should run without failures and produce correct output") {

            when {
                params {
                    // define parameters here. Example:
                    // outdir = "tests/results"
                }
                process {
                    """
                    input[0] = "${projectDir}/greetings.csv"
                    """
                }
            }

            then {
                assert process.success
                assert snapshot(process.out).match()
            }

        }
    ```

=== "Prima"

    ```groovy title="tests/main.converttoupper.nf.test" linenums="1" hl_lines="1 10 11"
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
    ```

Ed eseguiamo il test!

```bash title="nf-test pipeline pass"
nf-test test tests/main.converttoupper.nf.test
```

```console title="nf-test process convertToUpper pass"
> nf-test test tests/main.converttoupper.nf.test

🚀 nf-test 0.9.3
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Process convertToUpper

  Test [c59b6044] 'Should run without failures and produce correct output' PASSED (1.755s)
  Snapshots:
    1 created [Should run without failures and produce correct output]


Snapshot Summary:
  1 created

SUCCESS: Executed 1 tests in 1.764s
```

Notate che abbiamo creato un file snapshot per il processo `convertToUpper` in `tests/main.converttoupper.nf.test.snap`. Se eseguiamo di nuovo il test, dovremmo vedere che nf-test passa di nuovo.

```bash title="nf-test process convertToUpper pass"
nf-test test tests/main.converttoupper.nf.test
```

```console title="nf-test process convertToUpper pass"
> nf-test test tests/main.converttoupper.nf.test

🚀 nf-test 0.9.3
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Process convertToUpper

  Test [c59b6044] 'Should run without failures and produce correct output' PASSED (1.798s)


SUCCESS: Executed 1 tests in 1.811s
```

### Takeaway

Sapete come scrivere test per un processo Nextflow ed eseguirli.

### Cosa c'è dopo?

Imparate come eseguire i test per tutto in una volta sola!

## 3. Eseguire i test per l'intero repository

Eseguire nf-test su ogni componente va bene, ma è laborioso e soggetto a errori. Non possiamo semplicemente testare tutto in una volta?

Sì, possiamo!

Eseguiamo nf-test sull'intero repository.

### 3.1. Eseguire nf-test sull'intero repository

Possiamo eseguire nf-test sull'intero repository eseguendo il comando `nf-test test`.

```bash
nf-test test .
```

Notate che stiamo usando semplicemente `.` per eseguire tutto dalla nostra directory corrente. Questo includerà ogni test!

```console title="nf-test repo pass"
> nf-test test .

🚀 nf-test 0.9.3
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Process convertToUpper

  Test [3d26d9af] 'Should run without failures and produce correct output' PASSED (4.155s)

Test Workflow main.nf

  Test [f183df37] 'Should run successfully with correct number of processes' PASSED (3.33s)
  Test [d7e32a32] 'Should produce correct output files' PASSED (3.102s)

Test Process sayHello

  Test [58df4e4b] 'Should run without failures and contain expected greeting' PASSED (2.614s)


SUCCESS: Executed 4 tests in 13.481s
```

Guardate un po'! Abbiamo eseguito 4 test, 1 per ogni processo e 2 per l'intera pipeline con un singolo comando. Immaginate quanto sia potente su una codebase di grandi dimensioni!

---

## Riepilogo

In questa side quest, avete imparato a sfruttare le funzionalità di nf-test per creare ed eseguire test per i singoli processi, nonché test end-to-end per l'intera pipeline.
Ora siete a conoscenza dei due principali approcci alla validazione dell'output, gli snapshot e le asserzioni dirette sul contenuto, e di quando usare l'uno o l'altro.
Sapete anche come eseguire i test uno per uno o per un intero progetto.

Applicare queste tecniche nel vostro lavoro vi permetterà di garantire che:

- Il vostro codice funzioni come previsto
- Le modifiche non rompano le funzionalità esistenti
- Altri sviluppatori possano contribuire con fiducia
- I problemi possano essere identificati e risolti rapidamente
- Il contenuto dell'output corrisponda alle aspettative

### Pattern chiave

<!-- TODO: Can we add snippets of code below to illustrate? -->

1. Test a livello di pipeline:
   - Test di successo di base
   - Verifica del conteggio dei processi
   - Controlli sull'esistenza dei file di output
2. Test a livello di processo
3. Due approcci alla validazione dell'output:
   - Uso degli snapshot per la verifica completa dell'output
   - Uso delle asserzioni dirette sul contenuto per controlli su contenuto specifico
4. Esecuzione di tutti i test in un repository con un singolo comando

### Risorse aggiuntive

Consultate la [documentazione di nf-test](https://www.nf-test.com/) per funzionalità di test più avanzate e best practice. Potreste voler:

- Aggiungere asserzioni più complete ai vostri test
- Scrivere test per casi limite e condizioni di errore
- Configurare l'integrazione continua per eseguire i test automaticamente
- Imparare altri tipi di test come i test di workflow e di modulo
- Esplorare tecniche di validazione del contenuto più avanzate

**Ricordate:** I test sono documentazione vivente di come il vostro codice dovrebbe comportarsi. Più test scrivete, e più specifiche sono le vostre asserzioni, più potete essere certi dell'affidabilità della vostra pipeline.

---

## Cosa c'è dopo?

Tornate al [menu delle Side Quest](../) o cliccate il pulsante in basso a destra della pagina per passare all'argomento successivo nell'elenco.
