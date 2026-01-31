# Testing con nf-test

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduzione assistita da IA - [scopri di più e suggerisci miglioramenti](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Essere in grado di testare sistematicamente che ogni parte del workflow faccia ciò che dovrebbe fare è fondamentale per la riproducibilità e la manutenzione a lungo termine, e può essere di grande aiuto durante il processo di sviluppo.

Prendiamoci un momento per parlare del perché il testing è così importante. Se state sviluppando un workflow, una delle prime cose che farete è prendere alcuni dati di test che sapete essere validi e che dovrebbero produrre un risultato. Aggiungete il primo processo alla pipeline e lo collegate ai vostri input per farlo funzionare. Poi, per verificare che tutto funzioni, lo eseguite sui dati di test. Supponendo che funzioni, passate al processo successivo ed eseguite di nuovo i dati di test. Ripetete questo processo fino a quando avete una pipeline di cui siete soddisfatti.

Poi, magari aggiungete un semplice parametro vero o falso come `--skip_process`. Ora dovete eseguire la pipeline due volte, una volta con ciascun parametro per assicurarvi che funzioni come previsto. Ma aspettate, come verifichiamo se `--skip_process` salta effettivamente il processo? Dobbiamo scavare negli output o controllare i file di log! Questo è fastidioso e soggetto a errori.

Man mano che sviluppate la vostra pipeline, diventerà rapidamente così complessa che testare manualmente ogni iterazione sarà lento e soggetto a errori. Inoltre, se trovate un errore, sarà molto difficile individuare esattamente da dove nella vostra pipeline proviene l'errore. È qui che entra in gioco il testing.

Il testing consente di verificare sistematicamente che ogni parte della pipeline funzioni come previsto. I vantaggi per uno sviluppatore di test ben scritti sono enormi:

- **Fiducia**: Poiché i test coprono l'intera pipeline, potete essere sicuri che la modifica di qualcosa non influisca su nient'altro
- **Affidabilità**: Quando più sviluppatori lavorano sulla pipeline, sanno che gli altri sviluppatori non hanno rotto la pipeline e ogni componente.
- **Trasparenza**: I test mostrano dove una pipeline sta fallendo e rendono più facile rintracciare il problema. Funzionano anche come una forma di documentazione, mostrando come eseguire un processo o un workflow.
- **Velocità**: Poiché i test sono automatizzati, possono essere eseguiti molto rapidamente e ripetutamente. Potete iterare rapidamente con meno paura di introdurre nuovi bug.

Ci sono molti tipi diversi di test che possiamo scrivere:

1. **Test a livello di modulo**: Per processi individuali
2. **Test a livello di workflow**: Per un singolo workflow
3. **Test a livello di pipeline**: Per la pipeline nel suo complesso
4. **Test di performance**: Per la velocità e l'efficienza della pipeline
5. **Stress test**: Valutazione delle prestazioni della pipeline in condizioni estreme per determinarne i limiti

Testare i processi individuali è analogo agli unit test in altri linguaggi. Testare il workflow o l'intera pipeline è analogo a quello che viene chiamato integration test in altri linguaggi, dove testiamo le interazioni dei componenti.

[**nf-test**](https://www.nf-test.com/) è uno strumento che consente di scrivere test a livello di modulo, workflow e pipeline. In breve, consente di verificare sistematicamente che ogni singola parte della pipeline funzioni come previsto, _in isolamento_.

### Obiettivi di apprendimento

In questa side quest, imparerete a utilizzare nf-test per scrivere un test a livello di workflow per la pipeline nonché test a livello di modulo per i tre processi che chiama.

Alla fine di questa side quest, sarete in grado di utilizzare efficacemente le seguenti tecniche:

- Inizializzare nf-test nel vostro progetto
- Generare test a livello di modulo e di workflow
- Aggiungere tipi comuni di asserzioni
- Comprendere quando utilizzare snapshot vs. asserzioni di contenuto
- Eseguire test per un intero progetto

Queste competenze vi aiuteranno a implementare una strategia di testing completa nei vostri progetti di pipeline, assicurando che siano più robusti e manutenibili.

### Prerequisiti

Prima di intraprendere questa side quest, dovreste:

- Aver completato il tutorial [Hello Nextflow](../hello_nextflow/README.md) o un corso per principianti equivalente.
- Essere a vostro agio con i concetti e i meccanismi di base di Nextflow (processi, canali, operatori, lavorare con file, metadati)

---

## 0. Iniziare

#### Aprire il codespace di formazione

Se non l'avete ancora fatto, assicuratevi di aprire l'ambiente di formazione come descritto in [Configurazione dell'Ambiente](../envsetup/index.md).

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### Spostarsi nella directory del progetto

Spostiamoci nella directory dove si trovano i file per questo tutorial.

```bash
cd side-quests/nf-test
```

Potete impostare VSCode per concentrarsi su questa directory:

```bash
code .
```

#### Rivedere i materiali

Troverete un file di workflow principale e un file CSV chiamato `greetings.csv` che contiene l'input alla pipeline.

```console title="Contenuto della directory"
.
├── greetings.csv
└── main.nf
```

Per una descrizione dettagliata dei file, vedere il [warmup da Hello Nextflow](../hello_nextflow/00_orientation.md).

Il workflow che testeremo è un sottoinsieme del workflow Hello costruito in [Hello Workflow](../hello_nextflow/03_hello_workflow.md).

??? example "Cosa fa il workflow Hello Nextflow?"

    Se non avete seguito la formazione [Hello Nextflow](../hello_nextflow/index.md), ecco una rapida panoramica di ciò che fa questo semplice workflow.

    Il workflow prende un file CSV contenente saluti, esegue quattro passaggi di trasformazione consecutivi su di essi e produce un singolo file di testo contenente un'immagine ASCII di un personaggio divertente che dice i saluti.

    I quattro passaggi sono implementati come processi Nextflow (`sayHello`, `convertToUpper`, `collectGreetings`, e `cowpy`) memorizzati in file di modulo separati.

    1. **`sayHello`:** Scrive ogni saluto nel proprio file di output (es. "Hello-output.txt")
    2. **`convertToUpper`:** Converte ogni saluto in maiuscolo (es. "HELLO")
    3. **`collectGreetings`:** Raccoglie tutti i saluti in maiuscolo in un singolo file batch
    4. **`cowpy`:** Genera arte ASCII utilizzando lo strumento `cowpy`

    I risultati vengono pubblicati in una directory chiamata `results/`, e l'output finale della pipeline (quando eseguita con parametri predefiniti) è un file di testo normale contenente arte ASCII di un personaggio che dice i saluti in maiuscolo.

    In questa side quest, utilizziamo una forma intermedia del workflow Hello che contiene solo i primi due processi. <!-- TODO: change this to use the full finished workflow as suggested in https://github.com/nextflow-io/training/issues/735 -->

Il sottoinsieme con cui lavoreremo è composto da due processi: `sayHello` e `convertToUpper`.
Potete vedere il codice completo del workflow qui sotto.

??? example "Codice del workflow"

    ```groovy title="main.nf"
    /*
    * Parametri della pipeline
    */
    params.input_file = "greetings.csv"

    /*
    * Usa echo per stampare 'Hello World!' sullo standard out
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
    * Usa un'utilità di sostituzione di testo per convertire il saluto in maiuscolo
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

#### Eseguire il workflow

Eseguiamo il workflow per assicurarci che funzioni come previsto.

```bash
nextflow run main.nf
```

```console title="Risultato dell'esecuzione del workflow"
 N E X T F L O W   ~  version 24.10.2

Launching `main.nf` [soggy_linnaeus] DSL2 - revision: bbf79d5c31

executor >  local (6)
[f7/c3be66] sayHello (3)       | 3 of 3 ✔
[cd/e15303] convertToUpper (3) | 3 of 3 ✔
```

CONGRATULAZIONI! Avete appena eseguito un test!

"Aspetta, cosa? Ho appena eseguito il workflow e ha funzionato! Come può essere un test?"

Bella domanda!

Analizziamo cosa è appena successo.

Avete eseguito il workflow con i parametri predefiniti, avete confermato che ha funzionato e siete soddisfatti dei risultati. Questa è l'essenza del testing. Se avete seguito il corso di formazione Hello Nextflow, avrete notato che abbiamo sempre iniziato ogni sezione eseguendo il workflow che stavamo usando come punto di partenza, per confermare che tutto fosse configurato correttamente.

Il testing del software essenzialmente fa questo processo per noi.

#### Rivedere l'incarico

La vostra sfida è aggiungere test standardizzati a questo workflow utilizzando nf-test, al fine di rendere facile verificare che ogni parte continui a funzionare come previsto nel caso vengano apportate ulteriori modifiche.

<!-- TODO: give a bit more details, similar to how it's done in the Metadata side quest -->

#### Checklist di preparazione

Pensate di essere pronti per iniziare?

- [ ] Comprendo l'obiettivo di questo corso e i suoi prerequisiti
- [ ] Il mio codespace è attivo e funzionante
- [ ] Ho impostato la mia directory di lavoro in modo appropriato
- [ ] Ho eseguito con successo il workflow
- [ ] Comprendo l'incarico

Se potete spuntare tutte le caselle, siete pronti per iniziare.

---

## 1. Inizializzare `nf-test`

Il pacchetto `nf-test` fornisce un comando di inizializzazione che configura alcune cose per consentirci di iniziare a sviluppare test per il nostro progetto.

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

Crea anche una directory `tests` contenente uno stub di file di configurazione.

### 1.1. Generare un nf-test

`nf-test` viene fornito con un set di strumenti per costruire file nf-test, risparmiandoci la maggior parte del lavoro. Questi rientrano nel sottocomando `generate`. Generiamo un test per la pipeline:

```bash
nf-test generate pipeline main.nf
```

```console title="Output"
> nf-test generate pipeline main.nf

Load source file '/workspaces/training/side-quests/nf-test/main.nf'
Wrote pipeline test file '/workspaces/training/side-quests/nf-test/tests/main.nf.test

SUCCESS: Generated 1 test files.
```

Questo creerà un file `main.nf.test` all'interno della directory `tests`. Questo è il nostro file di test a livello di pipeline. Se eseguite `tree tests/` dovreste vedere qualcosa del genere:

```console title="Contenuto della directory test"
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

- `when`: Le condizioni in cui il test dovrebbe essere eseguito. Questo include i parametri che verranno utilizzati per eseguire la pipeline.
- `then`: Le asserzioni che dovrebbero essere fatte. Questo include i risultati attesi della pipeline.

In parole semplici, la logica del test si legge come segue:
"**Quando** questi _parametri_ vengono forniti a questa _pipeline_, **allora** ci aspettiamo di vedere questi risultati."

Questo non è un test funzionale, dimostreremo come trasformarlo in uno nella prossima sezione.

### Una nota sui nomi dei test

Nell'esempio sopra, abbiamo utilizzato il nome predefinito "Should run without failures" che è appropriato per un test di base che verifica solo se la pipeline viene eseguita con successo. Tuttavia, man mano che aggiungiamo casi di test più specifici, dovremmo usare nomi più descrittivi che indicano cosa stiamo effettivamente testando. Per esempio:

- "Should convert input to uppercase" - quando testiamo funzionalità specifiche
- "Should handle empty input gracefully" - quando testiamo casi limite
- "Should respect max memory parameter" - quando testiamo vincoli di risorse
- "Should create expected output files" - quando testiamo la generazione di file

I buoni nomi dei test dovrebbero:

1. Iniziare con "Should" per rendere chiaro quale sia il comportamento atteso
2. Descrivere la funzionalità specifica o lo scenario che viene testato
3. Essere abbastanza chiari che se il test fallisce, sapete quale funzionalità è rotta

Man mano che aggiungiamo più asserzioni e casi di test specifici in seguito, useremo questi nomi più descrittivi per rendere chiaro cosa sta verificando ogni test.

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

1. nf-test ha tentato di eseguire la pipeline così com'è, utilizzando le impostazioni nel blocco `when`:

```groovy title="tests/main.nf.test"
when {
    params {
        // define parameters here. Example:
        // outdir = "tests/results"
    }
}
```

2. nf-test ha controllato lo stato della pipeline e l'ha confrontato con il blocco `when`:

```groovy title="tests/main.nf.test"
then {
    assert workflow.success
}
```

Notate come nf-test ha riportato che la pipeline è fallita e ha fornito il messaggio di errore da Nextflow:

```console title="Errore"
ERROR ~ No such file or directory: /workspaces/training/side-quests/nf-test/.nf-test/tests/693ba951a20fec36a5a9292ed1cc8a9f/greetings.csv
```

Quindi qual era il problema? Ricordate che la pipeline ha un file greetings.csv nella directory del progetto. Quando nf-test esegue la pipeline, cercherà questo file, ma non riesce a trovarlo. Il file è lì, cosa sta succedendo? Bene, se guardiamo il percorso possiamo vedere che il test sta avvenendo nel percorso `./nf-test/tests/longHashString/`. Proprio come Nextflow, nf-test crea una nuova directory per ogni test per mantenere tutto isolato. Il file di dati non si trova lì quindi dobbiamo correggere il percorso del file nel test originale.

Potreste chiedervi come faremo a puntare alla radice della pipeline nel test. Poiché questa è una situazione comune, nf-test ha una serie di variabili globali che possiamo usare per facilitarci la vita. Potete trovare l'elenco completo [qui](https://www.nf-test.com/docs/testcases/global_variables/) ma nel frattempo useremo la variabile `projectDir`, che indica la radice del progetto della pipeline.

_Prima:_

```groovy title="tests/main.nf.test" linenums="1" hl_lines="3 4"
when {
    params {
        // define parameters here. Example:
        // outdir = "tests/results"
    }
}
```

_Dopo:_

```groovy title="tests/main.nf.test" linenums="1" hl_lines="3"
when {
    params {
        input_file = "${projectDir}/greetings.csv"
    }
}
```

Eseguiamo di nuovo il test per vedere se funziona.

```bash title="nf-test pipeline pass"
nf-test test tests/main.nf.test
```

```console title="La pipeline passa"
> nf-test test tests/main.nf.test

🚀 nf-test 0.9.3
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Workflow main.nf

  Test [1d4aaf12] 'Should run without failures' PASSED (1.619s)


SUCCESS: Executed 1 tests in 1.626s
```

Successo! La pipeline viene eseguita con successo e il test passa. Eseguitelo quante volte volete e otterrete sempre lo stesso risultato!

Per impostazione predefinita, l'output di Nextflow è nascosto, ma per convincervi che nf-test stia effettivamente eseguendo il workflow, potete usare il flag `--verbose`:

```bash
nf-test test tests/main.nf.test --verbose
```

```console title="La pipeline esegue tutti i processi"
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

Un controllo semplice è assicurarsi che la nostra pipeline esegua tutti i processi che ci aspettiamo e non ne salti nessuno silenziosamente. Ricordate che la nostra pipeline esegue 6 processi, uno chiamato `sayHello` e uno chiamato `convertToUpper` per ciascuno dei 3 saluti.

Aggiungiamo un'asserzione al nostro test per verificare che la pipeline eseguite il numero atteso di processi. Aggiorneremo anche il nome del nostro test per riflettere meglio ciò che stiamo testando.

**Prima:**

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

**Dopo:**

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

Il nome del test ora riflette meglio ciò che stiamo effettivamente verificando - non solo che la pipeline viene eseguita senza fallire, ma che esegue il numero atteso di processi.

Eseguiamo di nuovo il test per vedere se funziona.

```bash title="nf-test pipeline pass"
nf-test test tests/main.nf.test
```

```console title="La pipeline passa con asserzioni"
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

**Prima:**

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

**Dopo:**

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

Eseguire di nuovo il test per vedere se funziona.

```bash title="nf-test pipeline pass"
nf-test test tests/main.nf.test
```

```console title="La pipeline passa con asserzioni sui file"
> nf-test test tests/main.nf.test

🚀 nf-test 0.9.3
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Workflow main.nf

  Test [f0e08a68] 'Should run successfully with correct number of processes' PASSED (8.144s)
  Test [d7e32a32] 'Should produce correct output files' PASSED (6.994s)


SUCCESS: Executed 2 tests in 15.165s
```

Successo! I test passano perché la pipeline è stata completata con successo, è stato eseguito il numero corretto di processi e i file di output sono stati creati. Questo dovrebbe anche mostrarvi quanto sia utile fornire quei nomi informativi per i vostri test.

Questa è solo la superficie, possiamo continuare a scrivere asserzioni per verificare i dettagli della pipeline, ma per ora passiamo a testare gli elementi interni della pipeline.

### Riepilogo

Sapete come scrivere un nf-test per una pipeline.

### Qual è il prossimo passo?

Imparate come testare un processo Nextflow.

---

## 2. Testare un processo Nextflow

Non dobbiamo scrivere test per ogni parte della pipeline, ma più test abbiamo, più possiamo essere completi sulla pipeline e più possiamo essere sicuri che funzioni come previsto. In questa sezione testeremo entrambi i processi nella pipeline come unità individuali.

### 2.1. Testare il processo `sayHello`

Iniziamo con il processo `sayHello`.

Utilizziamo di nuovo il comando `nf-test generate` per generare test per il processo.

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

Come prima, iniziamo con i dettagli del test, seguiti dai blocchi `when` e `then`. Tuttavia, abbiamo anche un blocco `process` aggiuntivo che ci consente di definire gli input al processo.

Eseguiamo il test per vedere se funziona.

```bash title="nf-test pipeline pass"
nf-test test tests/main.sayhello.nf.test
```

```console title="Il test del processo fallisce"
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

Il test fallisce perché il processo `sayHello` dichiara 1 input ma è stato chiamato con 0 argomenti. Risolviamo questo problema aggiungendo un input al processo. Ricordate da [Hello Workflow](../hello_nextflow/03_hello_workflow.md) (e dalla sezione di warmup sopra) che il nostro processo `sayHello` prende un singolo input di valore, che dovremo fornire. Dovremmo anche correggere il nome del test per riflettere meglio ciò che stiamo testando.

**Prima:**

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

**Dopo:**

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

### 2.2. Verificare lo snapshot creato dal test

Se guardiamo il file `tests/main.sayhello.nf.test`, possiamo vedere che utilizza un metodo `snapshot()` nel blocco di asserzione:

```groovy title="tests/main.sayhello.nf.test"
assert snapshot(process.out).match()
```

Questo sta dicendo a nf-test di creare uno snapshot dell'output del processo `sayHello`. Diamo un'occhiata al contenuto del file di snapshot.

```console title="Contenuto del file snapshot"
code tests/main.sayhello.nf.test.snap
```

Non lo stamperemo qui, ma dovreste vedere un file JSON contenente i dettagli del processo e degli output del processo. In particolare, possiamo vedere una riga che assomiglia a questa:

```json title="Contenuto del file snapshot"
"0": [
    "hello-output.txt:md5,b1946ac92492d2347c6235b4d2611184"
]
```

Questo rappresenta gli output creati dal processo `sayHello`, che stiamo testando esplicitamente. Se rieseguiamo il test, il programma verificherà che il nuovo output corrisponda all'output che è stato originariamente registrato. Questo è un modo rapido e semplice per testare che gli output del processo non cambino, ed è per questo che nf-test lo fornisce come predefinito.

!!!warning "Avviso"

    Ciò significa che dobbiamo essere sicuri che l'output che registriamo nell'esecuzione originale sia corretto!

Se, nel corso dello sviluppo futuro, qualcosa nel codice cambia e causa un output diverso, il test fallirà e dovremo determinare se il cambiamento è previsto o meno.

- Se si scopre che qualcosa nel codice si è rotto, dovremo risolverlo, con l'aspettativa che il codice corretto passerà il test.
- Se è un cambiamento previsto (ad esempio, lo strumento è stato migliorato e i risultati sono migliori) allora dovremo aggiornare lo snapshot per accettare il nuovo output come riferimento da abbinare. nf-test ha un parametro `--update-snapshot` a questo scopo.

Possiamo eseguire di nuovo il test e vedere che il test dovrebbe passare:

```console title="nf-test process pass con snapshot"
> nf-test test tests/main.sayhello.nf.test

🚀 nf-test 0.9.3
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Process sayHello

  Test [f91a1bcd] 'Should run without failures and produce correct output' PASSED (1.675s)


SUCCESS: Executed 1 tests in 1.685s
```

Successo! Il test passa perché il processo `sayHello` è stato eseguito con successo e l'output corrisponde allo snapshot.

### 2.3. Alternativa agli Snapshot: Asserzioni Dirette sul Contenuto

Sebbene gli snapshot siano ottimi per rilevare qualsiasi modifica nell'output, a volte si desidera verificare contenuti specifici senza essere così rigorosi sull'intera corrispondenza del file. Per esempio:

- Quando parti dell'output potrebbero cambiare (timestamp, ID casuali, ecc.) ma determinati contenuti chiave devono essere presenti
- Quando si desidera verificare modelli o valori specifici nell'output
- Quando si desidera rendere il test più esplicito su cosa costituisce il successo

Ecco come potremmo modificare il nostro test per verificare contenuti specifici:

**Prima:**

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

**Dopo:**

```groovy title="tests/main.sayhello.nf.test" linenums="1" hl_lines="1 5 16 17"
    test("Should run without failures and contain expected greeting") {

        when {
            params {
                // definire parametri qui
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

Notate che nf-test vede gli output del processo come una lista di liste, quindi `process.out[0][0]` sta recuperando la prima parte del primo elemento del canale (o 'emissione') da questo processo.

Questo approccio:

- Rende chiaro esattamente cosa ci aspettiamo nell'output
- È più resiliente a modifiche irrilevanti nell'output
- Fornisce messaggi di errore migliori quando i test falliscono
- Consente validazioni più complesse (modelli regex, confronti numerici, ecc.)

Eseguiamo il test per vedere se funziona.

```bash title="nf-test pipeline pass"
nf-test test tests/main.sayhello.nf.test
```

```console title="Il test del processo fallisce"
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

Questo è un test simile al processo `sayHello`, ma sta testando il processo `convertToUpper`. Sappiamo che questo fallirà perché proprio come con `sayHello`, il processo `convertToUpper` prende un singolo input path, ma non ne abbiamo specificato uno.

Ora dobbiamo fornire un singolo file di input al processo convertToUpper, che include del testo che vogliamo convertire in maiuscolo. Ci sono molti modi per farlo:

- Potremmo creare un file dedicato da testare
- Potremmo riutilizzare il file data/greetings.csv esistente
- Potremmo crearlo al volo all'interno del test

Per ora, riutilizziamo il file data/greetings.csv esistente utilizzando l'esempio che abbiamo usato con il test a livello di pipeline. Come prima, possiamo nominare il test per riflettere meglio ciò che stiamo testando, ma questa volta lasciamo che crei uno 'snapshot' del contenuto piuttosto che verificare stringhe specifiche (come abbiamo fatto nell'altro processo).

**Prima:**

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

**Dopo:**

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

Ed eseguire il test!

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

Notate, abbiamo creato un file di snapshot per il processo `convertToUpper` in `tests/main.converttoupper.nf.test.snap`. Se eseguiamo di nuovo il test, dovremmo vedere che nf-test passa di nuovo.

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

### Riepilogo

Sapete come scrivere test per un processo Nextflow ed eseguirli.

### Qual è il prossimo passo?

Imparate come eseguire test per tutto in una volta!

## 3. Eseguire test per l'intero repository

Eseguire nf-test su ogni componente va bene, ma è laborioso e soggetto a errori. Non possiamo testare tutto in una volta?

Sì, possiamo!

Eseguiamo nf-test sull'intero repo.

### 3.1. Eseguire nf-test sull'intero repo

Possiamo eseguire nf-test sull'intero repo eseguendo il comando `nf-test test`.

```bash
nf-test test .
```

Notate, stiamo solo usando il `.` per eseguire tutto dalla nostra directory corrente. Questo includerà ogni test!

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

Guardate! Abbiamo eseguito 4 test, 1 per ogni processo e 2 per l'intera pipeline con un singolo comando. Immaginate quanto sia potente questo su una grande base di codice!

---

## Riepilogo

In questa side quest, avete imparato a sfruttare le funzionalità di nf-test per creare ed eseguire test per processi individuali nonché test end-to-end per l'intera pipeline.
Ora siete consapevoli dei due principali approcci alla validazione dell'output, snapshot e asserzioni dirette sul contenuto, e quando utilizzare l'uno o l'altro.
Sapete anche come eseguire test uno per uno o per un intero progetto.

L'applicazione di queste tecniche nel vostro lavoro vi consentirà di garantire che:

- Il vostro codice funzioni come previsto
- Le modifiche non rompano la funzionalità esistente
- Altri sviluppatori possano contribuire con fiducia
- I problemi possano essere identificati e risolti rapidamente
- Il contenuto dell'output corrisponda alle aspettative

### Pattern chiave

<!-- TODO: Can we add snippets of code below to illustrate? -->

1. Test a livello di pipeline:
   - Testing di successo di base
   - Verifica del conteggio dei processi
   - Controlli di esistenza dei file di output
2. Test a livello di processo
3. Due approcci alla validazione dell'output:
   - Utilizzo di snapshot per la verifica completa dell'output
   - Utilizzo di asserzioni dirette sul contenuto per controlli di contenuto specifici
4. Esecuzione di tutti i test in un repository con un singolo comando

### Risorse aggiuntive

Consultate la [documentazione di nf-test](https://www.nf-test.com/) per funzionalità di testing più avanzate e best practice. Potreste voler:

- Aggiungere asserzioni più complete ai vostri test
- Scrivere test per casi limite e condizioni di errore
- Impostare l'integrazione continua per eseguire i test automaticamente
- Imparare su altri tipi di test come test di workflow e moduli
- Esplorare tecniche di validazione del contenuto più avanzate

**Ricordate:** I test sono documentazione viva di come il vostro codice dovrebbe comportarsi. Più test scrivete, e più specifiche sono le vostre asserzioni, più potete essere sicuri dell'affidabilità della vostra pipeline.

---

## Qual è il prossimo passo?

Tornate al [menu delle Side Quest](./index.md) o cliccate il pulsante in basso a destra della pagina per passare al prossimo argomento nell'elenco.
