# Elaborazione di file in input

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduzione assistita da IA - [scopri di più e suggerisci miglioramenti](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

I workflow di analisi scientifica spesso comportano l'elaborazione di un gran numero di file.
Nextflow fornisce strumenti potenti per gestire i file in modo efficiente, aiutandola ad organizzare ed elaborare i suoi dati con codice minimo.

### Obiettivi di apprendimento

In questa side quest, esploreremo come Nextflow gestisce i file, dalle operazioni di base a tecniche più avanzate per lavorare con collezioni di file.
Imparerà ad estrarre metadati dai nomi dei file, che è un requisito comune nei pipeline di analisi scientifica.

Al termine di questa side quest, sarà in grado di:

- Creare oggetti Path da stringhe di percorsi file utilizzando il metodo `file()` di Nextflow
- Accedere agli attributi dei file come nome, estensione e directory principale
- Gestire trasparentemente sia file locali che remoti utilizzando URI
- Utilizzare i canali per automatizzare la gestione dei file con `channel.fromPath()` e `channel.fromFilePairs()`
- Estrarre e strutturare metadati dai nomi dei file utilizzando la manipolazione di stringhe
- Raggruppare file correlati utilizzando pattern matching ed espressioni glob
- Integrare le operazioni sui file nei processi Nextflow con una gestione appropriata degli input
- Organizzare gli output dei processi utilizzando strutture di directory basate sui metadati

Queste competenze la aiuteranno a costruire workflow che possono gestire diversi tipi di input di file con grande flessibilità.

### Prerequisiti

Prima di intraprendere questa side quest, dovrebbe:

- Aver completato il tutorial [Hello Nextflow](../../hello_nextflow/) o un corso equivalente per principianti.
- Essere a suo agio nell'utilizzo di concetti e meccanismi base di Nextflow (processi, canali, operatori)

<!-- I removed the suggestion to do the metamaps SQ first because that works more naturally after -->

---

## 0. Iniziare

#### Aprire il codespace per la formazione

Se non l'ha ancora fatto, assicuratevi di aprire l'ambiente di formazione come descritto in [Configurazione dell'Ambiente](../envsetup/index.md).

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### Spostarsi nella directory del progetto

Spostiamoci nella directory dove si trovano i file per questo tutorial.

```bash
cd side-quests/working_with_files
```

Può configurare VSCode per focalizzarsi su questa directory:

```bash
code .
```

#### Esaminare i materiali

Troverà un semplice file di workflow chiamato `main.nf`, una directory `modules` contenente due file modulo, e una directory `data` contenente alcuni file di dati di esempio.

??? abstract "Contenuti della directory"

    ```console
    .
    ├── data
    │   ├── patientA_rep1_normal_R1_001.fastq.gz
    │   ├── patientA_rep1_normal_R2_001.fastq.gz
    │   ├── patientA_rep1_tumor_R1_001.fastq.gz
    │   ├── patientA_rep1_tumor_R2_001.fastq.gz
    │   ├── patientA_rep2_normal_R1_001.fastq.gz
    │   ├── patientA_rep2_normal_R2_001.fastq.gz
    │   ├── patientA_rep2_tumor_R1_001.fastq.gz
    │   ├── patientA_rep2_tumor_R2_001.fastq.gz
    │   ├── patientB_rep1_normal_R1_001.fastq.gz
    │   ├── patientB_rep1_normal_R2_001.fastq.gz
    │   ├── patientB_rep1_tumor_R1_001.fastq.gz
    │   ├── patientB_rep1_tumor_R2_001.fastq.gz
    │   ├── patientC_rep1_normal_R1_001.fastq.gz
    │   ├── patientC_rep1_normal_R2_001.fastq.gz
    │   ├── patientC_rep1_tumor_R1_001.fastq.gz
    │   └── patientC_rep1_tumor_R2_001.fastq.gz
    ├── main.nf
    └── modules
        ├── analyze_reads.nf
        └── count_lines.nf
    ```

Questa directory contiene dati di sequenziamento paired-end da tre pazienti (A, B, C).

Per ciascun paziente, abbiamo campioni di tipo `tumor` (tipicamente provenienti da biopsie tumorali) o `normal` (prelevati da tessuto sano o sangue).
Se non ha familiarità con l'analisi oncologica, sappia semplicemente che questo corrisponde ad un modello sperimentale che utilizza campioni tumore/normale appaiati per eseguire analisi contrastive.

Per il paziente A in particolare, abbiamo due set di repliche tecniche (ripetizioni).

I file di dati di sequenziamento sono denominati con una convenzione tipica `_R1_` e `_R2_` per quelli che sono conosciuti come 'forward reads' e 'reverse reads'.

_Non si preoccupi se non ha familiarità con questo disegno sperimentale, non è critico per comprendere questo tutorial._

#### Esaminare l'incarico

La sua sfida è scrivere un workflow Nextflow che:

1. **Carichi** i file di input utilizzando i metodi di gestione file di Nextflow
2. **Estragga** i metadati (ID paziente, replica, tipo di campione) dalla struttura del nome del file
3. **Raggruppi** i file appaiati (R1/R2) insieme utilizzando `channel.fromFilePairs()`
4. **Elabori** i file con un modulo di analisi fornito
5. **Organizzi** gli output in una struttura di directory basata sui metadati estratti

#### Checklist di preparazione

Pensa di essere pronto per iniziare?

- [ ] Comprendo l'obiettivo di questo corso e i suoi prerequisiti
- [ ] Il mio codespace è attivo e funzionante
- [ ] Ho impostato la mia directory di lavoro in modo appropriato
- [ ] Comprendo l'incarico

Se può spuntare tutte le caselle, è pronto per iniziare.

---

## 1. Operazioni di base sui file

### 1.1. Identificare il tipo di un oggetto con `.class`

Dia un'occhiata al file workflow `main.nf`:

```groovy title="main.nf" linenums="1"
#!/usr/bin/env nextflow

workflow {

    // Create a Path object from a string path
    myFile = 'data/patientA_rep1_normal_R1_001.fastq.gz'

    println "${myFile} is of class ${myFile.class}"
}
```

Questo è un mini-workflow (senza processi) che fa riferimento ad un singolo percorso di file nel suo workflow, poi lo stampa sulla console, insieme alla sua classe.

??? info "Cos'è `.class`?"

    In Nextflow, `.class` ci dice che tipo di oggetto stiamo utilizzando. È come chiedere "che tipo di cosa è questa?" per scoprire se è una stringa, un numero, un file, o qualcos'altro.
    Questo ci aiuterà ad illustrare la differenza tra una semplice stringa e un oggetto Path nelle prossime sezioni.

Eseguiamo il workflow:

```bash
nextflow run main.nf
```

??? success "Output del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [romantic_chandrasekhar] DSL2 - revision: 5a4a89bc3a

    data/patientA_rep1_normal_R1_001.fastq.gz is of class java.lang.String
    ```

Come potete vedere, Nextflow ha stampato il percorso stringa esattamente come l'abbiamo scritto.

Questo è solo output di testo; Nextflow non ha ancora fatto nulla di speciale con esso.
Abbiamo anche confermato che per quanto riguarda Nextflow, questa è solo una stringa (di classe `java.lang.String`).
Ha senso, dal momento che non abbiamo ancora detto a Nextflow che corrisponde ad un file.

### 1.2. Creare un oggetto Path con file()

Possiamo dire a Nextflow come gestire i file creando [oggetti Path](https://www.nextflow.io/docs/latest/reference/stdlib-types.html#path) da stringhe di percorsi.

Nel nostro workflow, possiamo convertire il percorso stringa `data/patientA_rep1_normal_R1_001.fastq.gz` in un oggetto Path utilizzando il metodo `file()`, che fornisce accesso alle proprietà e operazioni sui file.

Modifichi `main.nf` per avvolgere la stringa con `file()` come segue:

=== "Dopo"

    ```groovy title="main.nf" linenums="5" hl_lines="2"
        // Create a Path object from a string path
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        println "${myFile} is of class ${myFile.class}"
    ```

=== "Prima"

    ```groovy title="main.nf" linenums="5" hl_lines="2"
        // Create a Path object from a string path
        myFile = 'data/patientA_rep1_normal_R1_001.fastq.gz'

        println "${myFile} is of class ${myFile.class}"
    ```

Ora eseguite nuovamente il workflow:

```bash
nextflow run main.nf
```

??? success "Output del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [kickass_coulomb] DSL2 - revision: 5af44b1b59

    /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz is of class class sun.nio.fs.UnixPath
    ```

Questa volta, vede il percorso assoluto completo invece del percorso relativo che abbiamo fornito come input.

Nextflow ha convertito la nostra stringa in un oggetto Path e l'ha risolto nella posizione effettiva del file sul sistema.
Il percorso del file sarà ora assoluto, come in `/workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz`.

Noti anche che la classe dell'oggetto Path è `sun.nio.fs.UnixPath`: questo è il modo di Nextflow di rappresentare i file locali.
Come vedremo più avanti, i file remoti avranno nomi di classi diversi (come `nextflow.file.http.XPath` per i file HTTP), ma funzionano tutti esattamente allo stesso modo e possono essere utilizzati in modo identico nei suoi workflow.

!!! tip

    **La differenza chiave:**

    - **Stringa di percorso**: Solo testo che Nextflow tratta come caratteri
    - **Oggetto Path**: Un riferimento intelligente al file con cui Nextflow può lavorare

    Pensate a questo: una stringa di percorso è come scrivere un indirizzo su carta, mentre un oggetto Path è come avere l'indirizzo caricato in un dispositivo GPS che sa come navigare lì e può dirle dettagli sul viaggio.

### 1.3. Accedere agli attributi del file

Perché questo è utile? Beh, ora che Nextflow comprende che `myFile` è un oggetto Path e non solo una stringa, possiamo accedere ai vari attributi dell'oggetto Path.

Aggiorniamo il nostro workflow per stampare gli attributi file integrati:

=== "Dopo"

    ```groovy title="main.nf" linenums="5" hl_lines="4-9"
        // Create a Path object from a string path
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        // Print file attributes
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
    ```

=== "Prima"

    ```groovy title="main.nf" linenums="5" hl_lines="4"
        // Create a Path object from a string path
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        println "${myFile} is of class ${myFile.class}"
    ```

Esegua il workflow:

```bash
nextflow run main.nf
```

??? success "Output del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [ecstatic_ampere] DSL2 - revision: f3fa3dcb48

    File object class: sun.nio.fs.UnixPath
    File name: patientA_rep1_normal_R1_001.fastq.gz
    Simple name: patientA_rep1_normal_R1_001
    Extension: gz
    Parent directory: /workspaces/training/side-quests/working_with_files/data
    ```

Vede i vari attributi del file stampati sulla console sopra.

### 1.4. Fornire il file ad un processo

La differenza tra stringhe e oggetti Path diventa critica quando si inizia a costruire workflow effettivi con processi.
Finora abbiamo verificato che Nextflow sta ora trattando il nostro file di input come un file, ma vediamo se possiamo effettivamente eseguire qualcosa su quel file in un processo.

#### 1.4.1. Importare il processo ed esaminare il codice

Le forniamo un modulo processo pre-scritto chiamato `COUNT_LINES` che prende un input file e conta quante righe contiene.

Per utilizzare il processo nel workflow, deve solo aggiungere un'istruzione include prima del blocco workflow:

=== "Dopo"

    ```groovy title="main.nf" linenums="1" hl_lines="3"
    #!/usr/bin/env nextflow

    include { COUNT_LINES } from './modules/count_lines.nf'

    workflow {
    ```

=== "Prima"

    ```groovy title="main.nf" linenums="1"
    #!/usr/bin/env nextflow

    workflow {
    ```

Può aprire il file modulo per esaminarne il codice:

```groovy title="modules/count_lines.nf" linenums="1"
#!/usr/bin/env nextflow

process COUNT_LINES {
    debug true

    input:
    path input_file

    script:
    """
    set -o pipefail
    echo "Processing file: $input_file"
    gzip -dc $input_file | wc -l
    """
}
```

Come potete vedere, è uno script abbastanza semplice che decomprime il file e conta quante righe contiene.

??? info "Cosa fa `debug true`?"

    La direttiva `debug true` nella definizione del processo fa sì che Nextflow stampi l'output dal suo script (come il conteggio righe "40") direttamente nel log di esecuzione.
    Senza questo, vedrebbe solo lo stato di esecuzione del processo ma non l'output effettivo dal suo script.

    Per maggiori informazioni sul debugging dei workflow Nextflow, consulti la side quest [Debugging Nextflow Workflows](debugging.md).

#### 1.4.2. Aggiungere una chiamata a `COUNT_LINES`

Ora che il processo è disponibile per il workflow, possiamo aggiungere una chiamata al processo `COUNT_LINES` per eseguirlo sul file di input.

Effettui le seguenti modifiche al workflow:

=== "Dopo"

    ```groovy title="main.nf" linenums="7" hl_lines="11-12"
        // Create a Path object from a string path
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        // Print file attributes
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"

        // Count the lines in the file
        COUNT_LINES(myFile)
    ```

=== "Prima"

    ```groovy title="main.nf" linenums="7" hl_lines="4-9"
        // Create a Path object from a string path
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        // Print file attributes
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
    ```

E ora eseguite il workflow:

```bash
nextflow run main.nf
```

??? success "Output del comando"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [cheeky_hypatia] DSL2 - revision: 281d13c414

    File object class: class sun.nio.fs.UnixPath
    File name: patientA_rep1_normal_R1_001.fastq.gz
    Simple name: patientA_rep1_normal_R1_001
    Extension: gz
    Parent directory: /workspaces/training/side-quests/working_with_files/data
    executor >  local (1)
    [e9/341c05] COUNT_LINES [100%] 1 of 1 ✔
    Processing file: /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz
    40
    ```

Questo mostra che siamo in grado di operare sul file in modo appropriato all'interno di un processo.

Nello specifico, Nextflow ha eseguito con successo le seguenti operazioni:

- Preparato il file nella directory di lavoro
- Decompresso il file .gz
- Contato le righe (40 righe in questo caso)
- Completato senza errori

La chiave di questa operazione senza intoppi è che stiamo esplicitamente dicendo a Nextflow che il nostro input è un file e dovrebbe essere trattato come tale.

### 1.5. Risolvere errori base di input file

Questo spesso fa inciampare i nuovi utenti di Nextflow, quindi dedichiamo qualche minuto a vedere cosa succede quando lo si fa in modo errato.

Ci sono due posti principali dove si può sbagliare la gestione dei file: a livello del workflow e a livello del processo.

#### 1.5.1. Errore a livello di workflow

Vediamo cosa succede se torniamo a trattare il file come una stringa quando specifichiamo l'input nel blocco workflow.

Effettui le seguenti modifiche al workflow, assicurandosi di commentare le istruzioni di stampa specifiche per i percorsi:

=== "Dopo"

    ```groovy title="main.nf" linenums="7" hl_lines="2 6-11"
        // Create a Path object from a string path
        myFile = 'data/patientA_rep1_normal_R1_001.fastq.gz'

        // Print file attributes
        println "File object class: ${myFile.class}"
        /*
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
        */

        // Count the lines in the file
        COUNT_LINES(myFile)
    ```

=== "Prima"

    ```groovy title="main.nf" linenums="7" hl_lines="4-9"
        // Create a Path object from a string path
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        // Print file attributes
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"

        // Count the lines in the file
        COUNT_LINES(myFile)
    ```

E ora eseguite il workflow:

```bash
nextflow run main.nf
```

??? failure "Output del comando"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [friendly_goodall] DSL2 - revision: ae50609b20

    [-        ] COUNT_LINES -
    ERROR ~ Error executing process > 'COUNT_LINES'

    Caused by:
      Not a valid path value: 'data/patientA_rep1_normal_R1_001.fastq.gz'



    Tip: view the complete command output by changing to the process work dir and entering the command `cat .command.out`

    -- Check '.nextflow.log' file for details
    ```

Questa è la parte importante:

```console
Not a valid path value: 'data/patientA_rep1_normal_R1_001.fastq.gz'
```

Quando specifica un input `path`, Nextflow valida che stia passando effettivi riferimenti a file, non solo stringhe.
Questo errore le sta dicendo che `'data/patientA_rep1_normal_R1_001.fastq.gz'` non è un valore di percorso valido perché è una stringa, non un oggetto Path.

Nextflow ha immediatamente rilevato il problema e si è fermato prima ancora di avviare il processo.

#### 1.5.2. Errore a livello di processo

L'altro posto dove potremmo dimenticare di specificare che vogliamo che Nextflow tratti l'input come un file è nella definizione del processo.

!!! warning "Mantenga l'errore del workflow da 1.5.1"

    Affinché questo test funzioni correttamente, mantenga il workflow nel suo stato errato (usando una semplice stringa invece di `file()`).
    Quando combinato con `val` nel processo, questo produce l'errore mostrato sotto.

Effettui la seguente modifica al modulo:

=== "Dopo"

    ```groovy title="modules/count_lines.nf" linenums="3" hl_lines="5"
    process COUNT_LINES {
        debug true

        input:
        val input_file
    ```

=== "Prima"

    ```groovy title="modules/count_lines.nf" linenums="3" hl_lines="5"
    process COUNT_LINES {
        debug true

        input:
        path input_file
    ```

E ora eseguite nuovamente il workflow:

```bash
nextflow run main.nf
```

??? failure "Output del comando"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [soggy_golick] DSL2 - revision: ae50609b20

    executor >  local (1)
    [b3/b3023c] COUNT_LINES [  0%] 0 of 1 ✘
    ERROR ~ Error executing process > 'COUNT_LINES'

    Caused by:
      Process `COUNT_LINES` terminated with an error exit status (1)


    Command executed:

      set -o pipefail
      echo "Processing file: data/patientA_rep1_normal_R1_001.fastq.gz"
      gzip -dc data/patientA_rep1_normal_R1_001.fastq.gz | wc -l

    Command exit status:
      1

    Command output:
      Processing file: data/patientA_rep1_normal_R1_001.fastq.gz
      0

    Command error:
      Processing file: data/patientA_rep1_normal_R1_001.fastq.gz
      gzip: data/patientA_rep1_normal_R1_001.fastq.gz: No such file or directory
      0

    Work dir:
      /workspaces/training/side-quests/working_with_files/work/b3/b3023cb2ccb986851301d8e369e79f

    Tip: when you have fixed the problem you can continue the execution adding the option `-resume` to the run command line

    -- Check '.nextflow.log' file for details
    ```

Questo mostra molti dettagli sull'errore perché il processo è impostato per produrre informazioni di debugging, come notato sopra.

Queste sono le sezioni più rilevanti:

```console
Command executed:

  set -o pipefail
  echo "Processing file: data/patientA_rep1_normal_R1_001.fastq.gz"
  gzip -dc data/patientA_rep1_normal_R1_001.fastq.gz | wc -l
```

```console
Command error:
  Processing file: data/patientA_rep1_normal_R1_001.fastq.gz
  gzip: data/patientA_rep1_normal_R1_001.fastq.gz: No such file or directory
  0
```

Questo dice che il sistema non ha potuto trovare il file; tuttavia se cerca il percorso, c'è un file con quel nome in quella posizione.

Quando abbiamo eseguito questo, Nextflow ha passato il valore stringa allo script, ma non ha _preparato_ il file effettivo nella directory di lavoro.
Quindi il processo ha provato ad utilizzare la stringa relativa, `data/patientA_rep1_normal_R1_001.fastq.gz`, ma quel file non esiste all'interno della directory di lavoro del processo.

Presi insieme, questi due esempi le mostrano quanto sia importante dire a Nextflow se un input dovrebbe essere gestito come un file.

!!! note

    Assicuratevi di tornare indietro e correggere entrambi gli errori intenzionali prima di continuare alla sezione successiva.

### Conclusioni

- Stringhe di percorso vs oggetti Path: Le stringhe sono solo testo, gli oggetti Path sono riferimenti intelligenti a file
- Il metodo `file()` converte un percorso stringa in un oggetto Path con cui Nextflow può lavorare
- È possibile accedere alle proprietà dei file come `name`, `simpleName`, `extension` e `parent` [utilizzando gli attributi file](https://www.nextflow.io/docs/latest/working-with-files.html#getting-file-attributes)
- Utilizzare oggetti Path invece di stringhe consente a Nextflow di gestire correttamente i file nel suo workflow
- Risultati Input Processo: Una corretta gestione dei file richiede oggetti Path, non stringhe, per garantire che i file siano correttamente preparati e accessibili per l'uso da parte dei processi.

---

## 2. Utilizzo di file remoti

Una delle caratteristiche chiave di Nextflow è la capacità di passare senza problemi tra file locali (sulla stessa macchina) e file remoti accessibili via internet.

Se lo fa correttamente, non dovrebbe mai aver bisogno di cambiare la logica del suo workflow per accogliere file provenienti da diverse posizioni.
Tutto ciò che deve fare per usare un file remoto è specificare il prefisso appropriato nel percorso del file quando lo fornisce al workflow.

Ad esempio, `/path/to/data` non ha prefisso, indicando che è un percorso file locale 'normale', mentre `s3://path/to/data` include il prefisso `s3://`, indicando che si trova nell'object storage S3 di Amazon.

Sono supportati molti protocolli diversi:

- HTTP(S)/FTP (http://, https://, ftp://)
- Amazon S3 (s3://)
- Azure Blob Storage (az://)
- Google Cloud Storage (gs://)

Per utilizzare uno qualsiasi di questi, specifichi semplicemente il prefisso rilevante nella stringa, che è poi tecnicamente chiamata Uniform Resource Identifier (URI) invece di percorso file.
Nextflow gestirà l'autenticazione e la preparazione dei file nel posto giusto, scaricando o caricando e tutte le altre operazioni sui file che ci si aspetterebbe.

La forza chiave di questo sistema è che ci consente di passare tra ambienti senza cambiare alcuna logica del pipeline.
Ad esempio, può sviluppare con un piccolo set di test locale prima di passare ad un set di test su larga scala situato in storage remoto semplicemente cambiando l'URI.

### 2.1. Utilizzare un file da internet

Testiamo questo sostituendo il percorso locale che stiamo fornendo al nostro workflow con un percorso HTTPS che punta ad una copia degli stessi dati memorizzata in Github.

!!! warning

    Questo funzionerà solo se avete una connessione internet attiva.

Apra `main.nf` di nuovo e cambi il percorso input come segue:

=== "Dopo"

    ```groovy title="main.nf" linenums="2" hl_lines="2"
        // Using a remote file from the internet
        myFile = file('https://raw.github.com/nextflow-io/training/master/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz')

        // Print file attributes
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
    ```

=== "Prima"

    ```groovy title="main.nf" linenums="2" hl_lines="2"
        // Create a Path object from a string path
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        // Print file attributes
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
    ```

Eseguiamo il workflow:

```bash
nextflow run main.nf
```

??? success "Output del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [insane_swartz] DSL2 - revision: fff18abe6d

    File object class: class nextflow.file.http.XPath
    File name: patientA_rep1_normal_R1_001.fastq.gz
    Simple name: patientA_rep1_normal_R1_001
    Extension: gz
    Parent directory: /nextflow-io/training/master/side-quests/working_with_files/data
    executor >  local (1)
    [8a/2ab7ca] COUNT_LINES [100%] 1 of 1 ✔
    Processing file: patientA_rep1_normal_R1_001.fastq.gz
    40
    ```

Funziona! Può vedere che è cambiato molto poco.

L'unica differenza nell'output della console è che la classe dell'oggetto percorso è ora `nextflow.file.http.XPath`, mentre per il percorso locale la classe era `sun.nio.fs.UnixPath`.
Non ha bisogno di ricordare queste classi; menzioniamo solo questo per dimostrare che Nextflow identifica e gestisce in modo appropriato le diverse posizioni.

Dietro le quinte, Nextflow ha scaricato il file in una directory di staging situata all'interno della directory di lavoro.
Quel file preparato può quindi essere trattato come un file locale e collegato simbolicamente nella directory del processo rilevante.

Può verificare che ciò sia accaduto qui guardando i contenuti della directory di lavoro situata al valore hash del processo.

??? abstract "Contenuti della directory di lavoro"

    Se l'hash del processo era `8a/2ab7ca`, potrebbe esplorare la directory di lavoro:

    ```console
    $ ls -la work/8a/2ab7ca*/
    total 16
    drwxr-xr-x  6 user  staff   192 Jan 28 10:00 .
    drwxr-xr-x  3 user  staff    96 Jan 28 10:00 ..
    -rw-r--r--  1 user  staff     0 Jan 28 10:00 .command.begin
    -rw-r--r--  1 user  staff   127 Jan 28 10:00 .command.sh
    lrwxr-xr-x  1 user  staff    89 Jan 28 10:00 patientA_rep1_normal_R1_001.fastq.gz -> /path/to/work/stage/.../patientA_rep1_normal_R1_001.fastq.gz
    ```

    Il link simbolico punta ad una copia preparata del file remoto che Nextflow ha scaricato automaticamente.

Noti che per file più grandi, il passaggio di download richiederà tempo extra rispetto all'esecuzione su file locali.
Tuttavia, Nextflow verifica se avete già una copia preparata per evitare download non necessari.
Quindi se esegue nuovamente sullo stesso file e non ha eliminato il file preparato, Nextflow utilizzerà la copia preparata.

Questo mostra quanto sia facile passare tra dati locali e remoti usando Nextflow, che è una caratteristica chiave di Nextflow.

!!! note

    L'unica importante eccezione a questo principio è che non potete utilizzare pattern glob o percorsi di directory con HTTPS perché HTTPS non può elencare più file, quindi deve specificare URL esatti di file.
    Tuttavia, altri protocolli di storage come blob storage (`s3://`, `az://`, `gs://`) possono utilizzare sia glob che percorsi di directory.

    Ecco come potrebbe utilizzare pattern glob con cloud storage:

    ```groovy title="Esempi di cloud storage (non eseguibili in questo ambiente)"
    // S3 con pattern glob - corrisponderebbe a più file
    ch_s3_files = channel.fromPath('s3://my-bucket/data/*.fastq.gz')

    // Azure Blob Storage con pattern glob
    ch_azure_files = channel.fromPath('az://container/data/patient*_R{1,2}.fastq.gz')

    // Google Cloud Storage con pattern glob
    ch_gcs_files = channel.fromPath('gs://bucket/data/sample_*.fastq.gz')
    ```

    Le mostreremo come lavorare con i glob in pratica nella prossima sezione.

### 2.2. Tornare al file locale

Torneremo ad utilizzare i nostri file di esempio locali per il resto di questa side quest, quindi cambiamo l'input del workflow tornando al file originale:

=== "Dopo"

    ```groovy title="main.nf" linenums="2" hl_lines="2"
        // Create a Path object from a string path
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        // Print file attributes
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
    ```

=== "Prima"

    ```groovy title="main.nf" linenums="2" hl_lines="2"
        // Create a Path object from a string path
        myFile = file('https://raw.github.com/nextflow-io/training/master/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz')

        // Print file attributes
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
    ```

### Conclusioni

- L'accesso ai dati remoti avviene utilizzando un URI (HTTP, FTP, S3, Azure, Google Cloud)
- Nextflow scaricherà e preparerà automaticamente i dati nel posto giusto, purché questi percorsi vengano forniti ai processi
- Non scriva logica per scaricare o caricare file remoti!
- I file locali e remoti producono tipi di oggetti diversi ma funzionano in modo identico
- **Importante**: HTTP/HTTPS funzionano solo con singoli file (nessun pattern glob)
- Il cloud storage (S3, Azure, GCS) supporta sia singoli file che pattern glob
- Può passare senza problemi tra sorgenti di dati locali e remote senza cambiare la logica del codice (purché il protocollo supporti le operazioni richieste)

---

## 3. Utilizzo del channel factory `fromPath()`

Finora abbiamo lavorato con un singolo file alla volta, ma in Nextflow, tipicamente vorremo creare un canale di input con più file di input da elaborare.

Un modo ingenuo per farlo sarebbe combinare il metodo `file()` con [`channel.of()`](https://www.nextflow.io/docs/latest/reference/channel.html#of) così:

```groovy title="Esempio di sintassi"
ch_files = channel.of([file('data/patientA_rep1_normal_R1_001.fastq.gz')],
                      [file('data/patientA_rep1_normal_R1_001.fastq.gz')])
```

Funziona, ma è goffo.

!!! tip "Quando usare `file()` vs `channel.fromPath()`"

    - Usi `file()` quando ha bisogno di un singolo oggetto Path per manipolazione diretta (controllare se un file esiste, leggere i suoi attributi, o passarlo ad una singola invocazione di processo)
    - Usi `channel.fromPath()` quando ha bisogno di un canale che può contenere più file, specialmente con pattern glob, o quando i file fluiranno attraverso più processi

È qui che entra in gioco [`channel.fromPath()`](https://www.nextflow.io/docs/latest/reference/channel.html#frompath): un comodo channel factory che raggruppa tutta la funzionalità di cui abbiamo bisogno per generare un canale da una o più stringhe di file statici così come pattern glob.

### 3.1. Aggiungere il channel factory

Aggiorniamo il nostro workflow per usare `channel.fromPath`.

=== "Dopo"

    ```groovy title="main.nf" linenums="7" hl_lines="1-3"
        // Load files with channel.fromPath
        ch_files = channel.fromPath('data/patientA_rep1_normal_R1_001.fastq.gz')
        ch_files.view { myFile -> "Found file: $myFile" }

        // Print file attributes
        /* Comment these out for now, we'll come back to them!
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"
        */

        // Count the lines in the file
        // COUNT_LINES(myFile)
    ```

=== "Prima"

    ```groovy title="main.nf" linenums="7" hl_lines="1-2"
        // Create a Path object from a string path
        myFile = file('data/patientA_rep1_normal_R1_001.fastq.gz')

        // Print file attributes
        println "File object class: ${myFile.class}"
        println "File name: ${myFile.name}"
        println "Simple name: ${myFile.simpleName}"
        println "Extension: ${myFile.extension}"
        println "Parent directory: ${myFile.parent}"

        // Count the lines in the file
        COUNT_LINES(myFile)
    ```

Abbiamo anche commentato il codice che stampa gli attributi per ora, e aggiunto un'istruzione `.view` per stampare solo il nome del file.

Esegua il workflow:

```bash
nextflow run main.nf
```

??? success "Output del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [grave_meucci] DSL2 - revision: b09964a583

    Found file: /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz
    ```

Come potete vedere, il percorso del file viene caricato come un oggetto di tipo `Path` nel canale.
Questo è simile a quello che avrebbe fatto `file()`, eccetto che ora abbiamo un canale in cui possiamo caricare più file se vogliamo.

Usare `channel.fromPath()` è un modo conveniente di creare un nuovo canale popolato da una lista di file.

### 3.2. Visualizzare gli attributi dei file nel canale

Nel nostro primo tentativo di utilizzare il channel factory, abbiamo semplificato il codice e stampato solo il nome del file.

Torniamo a stampare gli attributi completi del file:

=== "Dopo"

    ```groovy title="main.nf" linenums="7" hl_lines="3-9 12"
        // Load files with channel.fromPath
        ch_files = channel.fromPath('data/patientA_rep1_normal_R1_001.fastq.gz')
        ch_files.view { myFile ->
            println "File object class: ${myFile.class}"
            println "File name: ${myFile.name}"
            println "Simple name: ${myFile.simpleName}"
            println "Extension: ${myFile.extension}"
            println "Parent directory: ${myFile.parent}"
        }

        // Count the lines in the file
        COUNT_LINES(ch_files)
    ```

=== "Prima"

    ```groovy title="main.nf" linenums="7" hl_lines="3"
        // Load files with channel.fromPath
        ch_files = channel.fromPath('data/patientA_rep1_normal_R1_001.fastq.gz')
        ch_files.view { myFile -> "Found file: $myFile" }

        // Count the lines in the file
        // COUNT_LINES(ch_files)
    ```

Stiamo anche riabilitando la chiamata al processo `COUNT_LINES` per verificare che l'elaborazione dei file funzioni ancora correttamente con il nostro approccio basato sui canali.

Esegua il workflow:

```bash
nextflow run main.nf
```

??? success "Output del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [furious_swanson] DSL2 - revision: c35c34950d

    executor >  local (1)
    [9d/6701a6] COUNT_LINES (1) [100%] 1 of 1 ✔
    File object class: sun.nio.fs.UnixPath
    File name: patientA_rep1_normal_R1_001.fastq.gz
    Simple name: patientA_rep1_normal_R1_001
    Extension: gz
    Parent directory: /workspaces/training/side-quests/working_with_files/data
    Processing file: patientA_rep1_normal_R1_001.fastq.gz
    40
    ```

Ed eccolo, stessi risultati di prima ma ora abbiamo il file in un canale, quindi possiamo aggiungerne altri.

### 3.3. Utilizzare un glob per corrispondere a più file

Ci sono diversi modi in cui potremmo caricare più file nel canale.
Qui le mostreremo come utilizzare i pattern glob, che sono un modo conveniente per corrispondere e recuperare nomi di file e directory basati su caratteri jolly.
Il processo di corrispondenza di questi pattern è chiamato "globbing" o "espansione del nome file".

!!! note

    Come notato precedentemente, Nextflow supporta il globbing per gestire file di input e output nella maggior parte dei casi, eccetto con i percorsi file HTTPS perché HTTPS non può elencare più file.

Diciamo che vogliamo recuperare entrambi i file in una coppia di file associati ad un dato paziente, `patientA`:

```console
patientA_rep1_normal_R1_001.fastq.gz
patientA_rep1_normal_R2_001.fastq.gz
```

Poiché l'unica differenza tra i nomi dei file è il numero della replica, _ovvero_ il numero dopo `R`, possiamo usare il carattere jolly `*` per rappresentare il numero come segue:

```console
patientA_rep1_normal_R*_001.fastq.gz
```

Questo è il pattern glob di cui abbiamo bisogno.

Ora tutto ciò che dobbiamo fare è aggiornare il percorso del file nel channel factory per utilizzare quel pattern glob come segue:

=== "Dopo"

    ```groovy title="main.nf" linenums="7"
      // Load files with channel.fromPath
        ch_files = channel.fromPath('data/patientA_rep1_normal_R*_001.fastq.gz')
    ```

=== "Prima"

    ```groovy title="main.nf" linenums="7"
      // Load files with channel.fromPath
        ch_files = channel.fromPath('data/patientA_rep1_normal_R1_001.fastq.gz')
    ```

Nextflow riconoscerà automaticamente che questo è un pattern glob e lo gestirà in modo appropriato.

Esegua il workflow per testarlo:

```bash
nextflow run main.nf
```

??? success "Output del comando"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [boring_sammet] DSL2 - revision: d2aa789c9a

    executor >  local (2)
    [3c/a65de5] COUNT_LINES (2) [100%] 2 of 2 ✔
    File object class: class sun.nio.fs.UnixPath
    File name: patientA_rep1_normal_R1_001.fastq.gz
    Simple name: patientA_rep1_normal_R1_001
    Extension: gz
    Parent directory: /workspaces/training/side-quests/working_with_files/data
    File object class: class sun.nio.fs.UnixPath
    File name: patientA_rep1_normal_R2_001.fastq.gz
    Simple name: patientA_rep1_normal_R2_001
    Extension: gz
    Parent directory: /workspaces/training/side-quests/working_with_files/data
    Processing file: patientA_rep1_normal_R1_001.fastq.gz
    40

    Processing file: patientA_rep1_normal_R2_001.fastq.gz
    40
    ```

Come potete vedere, ora abbiamo due oggetti Path nel nostro canale, il che mostra che Nextflow ha fatto correttamente l'espansione del nome file, e ha caricato ed elaborato entrambi i file come previsto.

Utilizzando questo metodo, possiamo recuperare quanti file vogliamo semplicemente cambiando il pattern glob. Se lo rendessimo più generoso, ad esempio sostituendo tutte le parti variabili dei nomi file con `*` (_ad es._ `data/patient*_rep*_*_R*_001.fastq.gz`) potremmo prendere tutti i file di esempio nella directory `data`.

### Conclusioni

- `channel.fromPath()` crea un canale con file che corrispondono ad un pattern
- Ogni file viene emesso come elemento separato nel canale
- Possiamo usare un pattern glob per corrispondere a più file
- I file vengono automaticamente convertiti in oggetti Path con attributi completi
- Il metodo `.view()` consente l'ispezione dei contenuti del canale

---

## 4. Estrazione di metadati di base dai nomi dei file

Nella maggior parte dei domini scientifici, è molto comune avere metadati codificati nei nomi dei file che contengono i dati.
Ad esempio, in bioinformatica, i file contenenti dati di sequenziamento sono spesso denominati in modo da codificare informazioni sul campione, condizione, replica e numero di lettura.

Se i nomi dei file sono costruiti secondo una convenzione consistente, può estrarre quei metadati in modo standardizzato e utilizzarli nel corso della sua analisi.
Questo è un grande 'se', ovviamente, e dovrebbe essere molto cauto ogni volta che si basa sulla struttura del nome del file; ma la realtà è che questo approccio è molto ampiamente utilizzato, quindi diamo un'occhiata a come si fa in Nextflow.

Nel caso dei nostri dati di esempio, sappiamo che i nomi dei file includono metadati strutturati in modo consistente.
Ad esempio, il nome del file `patientA_rep1_normal_R2_001` codifica quanto segue:

- ID paziente: `patientA`
- ID replica: `rep1`
- tipo di campione: `normal` (in contrapposizione a `tumor`)
- set di letture: `R1` (in contrapposizione a `R2`)

Modificheremo il nostro workflow per recuperare queste informazioni in tre passaggi:

1. Recuperare il `simpleName` del file, che include i metadati
2. Separare i metadati usando un metodo chiamato `tokenize()`
3. Utilizzare una map per organizzare i metadati

!!! warning

    Non dovrebbe mai codificare informazioni sensibili nei nomi dei file, come nomi di pazienti o altre caratteristiche identificative, poiché ciò può compromettere la privacy dei pazienti o altre restrizioni di sicurezza rilevanti.

### 4.1. Recuperare il `simpleName`

Il `simpleName` è un attributo file che corrisponde al nome del file privato del suo percorso e dell'estensione.

Effettui le seguenti modifiche al workflow:

=== "Dopo"

    ```groovy title="main.nf" linenums="7" hl_lines="3-6"
        // Load files with channel.fromPath
        ch_files = channel.fromPath('data/patientA_rep1_normal_R*_001.fastq.gz')
        ch_files.map { myFile ->
            [ myFile.simpleName, myFile ]
        }
        .view()
    ```

=== "Prima"

    ```groovy title="main.nf" linenums="7" hl_lines="3-9"
        // Load files with channel.fromPath
        ch_files = channel.fromPath('data/patientA_rep1_normal_R*_001.fastq.gz')
        ch_files.view { myFile ->
            println "File object class: ${myFile.class}"
            println "File name: ${myFile.name}"
            println "Simple name: ${myFile.simpleName}"
            println "Extension: ${myFile.extension}"
            println "Parent directory: ${myFile.parent}"
        }
    ```

Questo recupera il `simpleName` e lo associa all'oggetto file completo usando un'operazione `map()`.

Esegua il workflow per verificare che funzioni:

```bash
nextflow run main.nf
```

??? success "Output del comando"

    ```console hl_lines="7-8"
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [suspicious_mahavira] DSL2 - revision: ae8edc4e48

    executor >  local (2)
    [e9/55774b] COUNT_LINES (2) [100%] 2 of 2 ✔
    [patientA_rep1_normal_R2_001, /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R2_001.fastq.gz]
    [patientA_rep1_normal_R1_001, /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz]
    Processing file: patientA_rep1_normal_R1_001.fastq.gz
    40

    Processing file: patientA_rep1_normal_R2_001.fastq.gz
    40
    ```

Ogni elemento nel canale è ora una tupla contenente il `simpleName` e l'oggetto file originale.

### 4.2. Estrarre i metadati dal `simplename`

A questo punto, i metadati che vogliamo sono incorporati nel `simplename`, ma non possiamo accedere direttamente ai singoli elementi.
Quindi dobbiamo dividere il `simplename` nelle sue componenti.
Fortunatamente, quelle componenti sono semplicemente separate da underscore nel nome file originale, quindi possiamo applicare un metodo comune di Nextflow chiamato `tokenize()` che è perfetto per questo compito.

Effettui le seguenti modifiche al workflow:

=== "Dopo"

    ```groovy title="main.nf" linenums="7" hl_lines="4"
        // Load files with channel.fromPath
        ch_files = channel.fromPath('data/patientA_rep1_normal_R*_001.fastq.gz')
        ch_files.map { myFile ->
            [ myFile.simpleName.tokenize('_'), myFile ]
        }
    ```

=== "Prima"

    ```groovy title="main.nf" linenums="7" hl_lines="4"
        // Load files with channel.fromPath
        ch_files = channel.fromPath('data/patientA_rep1_normal_R*_001.fastq.gz')
        ch_files.map { myFile ->
            [ myFile.simpleName, myFile ]
        }
    ```

Il metodo `tokenize()` dividerà la stringa `simpleName` ovunque trovi underscore, e restituirà una lista contenente le sottostringhe.

Esegua il workflow:

```bash
nextflow run main.nf
```

??? success "Output del comando"

    ```console hl_lines="7-8"
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [gigantic_gauss] DSL2 - revision: a39baabb57

    executor >  local (2)
    [e7/da2f4b] COUNT_LINES (2) [100%] 2 of 2 ✔
    [[patientA, rep1, normal, R2, 001], /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R2_001.fastq.gz]
    [[patientA, rep1, normal, R1, 001], /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz]
    Processing file: patientA_rep1_normal_R2_001.fastq.gz
    40

    Processing file: patientA_rep1_normal_R1_001.fastq.gz
    40
    ```

Ora la tupla per ogni elemento nel nostro canale contiene la lista di metadati (_ad es._ `[patientA, rep1, normal, R1, 001]`) e l'oggetto file originale.

Fantastico!
Abbiamo suddiviso le nostre informazioni sul paziente da una singola stringa in una lista di stringhe.
Possiamo ora gestire ogni parte delle informazioni sul paziente separatamente.

### 4.3. Utilizzare una map per organizzare i metadati

I nostri metadati sono solo una lista piatta al momento.
È abbastanza facile da usare ma difficile da leggere.

```console
[patientA, rep1, normal, R1, 001]
```

Qual è l'elemento all'indice 3? Può dirlo senza fare riferimento alla spiegazione originale della struttura dei metadati?

Questa è una grande opportunità per usare un archivio chiave-valore, dove ogni elemento ha un set di chiavi e i loro valori associati, quindi può facilmente fare riferimento ad ogni chiave per ottenere il valore corrispondente.

Nel nostro esempio, ciò significa passare da questa organizzazione:

```groovy
data = [patientA, 1, normal, R1]

println data[3]
```

A questa:

```groovy
data = [id: patientA, replicate: 1, type: normal, readNum: 1]

println data.readNum
```

In Nextflow, questo è chiamato [map](https://nextflow.io/docs/latest/script.html#maps).

Convertiamo ora la nostra lista piatta in una map.
Effettui le seguenti modifiche al workflow:

=== "Dopo"

    ```groovy title="main.nf" linenums="7" hl_lines="4-13"
        // Load files with channel.fromPath
        ch_files = channel.fromPath('data/patientA_rep1_normal_R*_001.fastq.gz')
        ch_files.map { myFile ->
            def (patient, replicate, type, readNum) = myFile.simpleName.tokenize('_')
            [
              [
                id: patient,
                replicate: replicate.replace('rep', ''),
                type: type,
                readNum: readNum.replace('R', ''),
              ],
              myFile
            ]
        }
    ```

=== "Prima"

    ```groovy title="main.nf" linenums="7" hl_lines="4"
        // Load files with channel.fromPath
        ch_files = channel.fromPath('data/patientA_rep1_normal_R*_001.fastq.gz')
        ch_files.map { myFile ->
            [ myFile.simpleName.tokenize('_'), myFile ]
        }
    ```

I cambiamenti chiave qui sono:

- **Assegnazione destrutturante**: `def (patient, replicate, type, readNum) = ...` estrae i valori tokenizzati in variabili nominate in una sola riga
- **Sintassi letterale map**: `[id: patient, replicate: ...]` crea una map dove ogni chiave (come `id`) è associata ad un valore (come `patient`)
- **Struttura annidata**: La lista esterna `[..., myFile]` accoppia la map dei metadati con l'oggetto file originale

Abbiamo anche semplificato un paio di stringhe di metadati usando un metodo di sostituzione di stringhe chiamato `replace()` per rimuovere alcuni caratteri non necessari (_ad es._ `replicate.replace('rep', '')` per mantenere solo il numero dagli ID replica).

Eseguiamo nuovamente il workflow:

```bash
nextflow run main.nf
```

??? success "Output del comando"

    ```console hl_lines="7-8"
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [infallible_swartz] DSL2 - revision: 7f4e68c0cb

    executor >  local (2)
    [1b/e7fb27] COUNT_LINES (1) [100%] 2 of 2 ✔
    [[id:patientA, replicate:1, type:normal, readNum:2], /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R2_001.fastq.gz]
    [[id:patientA, replicate:1, type:normal, readNum:1], /workspaces/training/side-quests/working_with_files/data/patientA_rep1_normal_R1_001.fastq.gz]
    Processing file: patientA_rep1_normal_R2_001.fastq.gz
    40

    Processing file: patientA_rep1_normal_R1_001.fastq.gz
    40
    ```

Ora i metadati sono etichettati in modo ordinato (_ad es._ `[id:patientA, replicate:1, type:normal, readNum:2]`) quindi è molto più facile distinguere cosa è cosa.

Sarà anche molto più facile effettivamente utilizzare elementi di metadati nel workflow, e renderà il nostro codice più facile da leggere e più manutenibile.

### Conclusioni

- Possiamo gestire i nomi dei file in Nextflow con la potenza di un linguaggio di programmazione completo
- Possiamo trattare i nomi dei file come stringhe per estrarre informazioni rilevanti
- L'uso di metodi come `tokenize()` e `replace()` ci consente di manipolare le stringhe nel nome del file
- L'operazione `.map()` trasforma gli elementi del canale preservando la struttura
- Metadati strutturati (maps) rendono il codice più leggibile e manutenibile rispetto alle liste posizionali

Prossimamente, vedremo come gestire file di dati appaiati.

---

## 5. Gestione di file di dati appaiati

Molti disegni sperimentali producono file di dati appaiati che beneficiano di essere gestiti in modo esplicitamente appaiato.
Ad esempio, in bioinformatica, i dati di sequenziamento sono spesso generati sotto forma di letture appaiate, ovvero stringhe di sequenza che originano dallo stesso frammento di DNA (spesso chiamate 'forward' e 'reverse' perché vengono lette da estremità opposte).

Questo è il caso dei nostri dati di esempio, dove R1 e R2 si riferiscono ai due set di letture.

```console
data/patientA_rep1_normal_R1_001.fastq.gz
data/patientA_rep1_normal_R2_001.fastq.gz
```

Nextflow fornisce un channel factory specializzato per lavorare con file appaiati come questo chiamato `channel.fromFilePairs
