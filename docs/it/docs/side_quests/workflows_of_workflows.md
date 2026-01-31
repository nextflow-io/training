# Workflows di Workflows

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traduzione assistita da IA - [scopri di più e suggerisci miglioramenti](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Quando si sviluppa una pipeline, ci si ritrova spesso a creare sequenze simili di processi per diversi tipi di dati o fasi di analisi. Si potrebbe finire per copiare e incollare queste sequenze di processi, con conseguente duplicazione del codice difficile da mantenere; oppure si potrebbe creare un unico workflow massiccio difficile da comprendere e modificare.

Una delle funzionalità più potenti di Nextflow è la sua capacità di comporre pipeline complesse da moduli workflow più piccoli e riutilizzabili. Questo approccio modulare rende le pipeline più facili da sviluppare, testare e mantenere.

### Obiettivi di apprendimento

In questa missione secondaria, esploreremo come sviluppare moduli workflow che possono essere testati e utilizzati separatamente, comporre questi moduli in una pipeline più grande e gestire il flusso di dati tra moduli.

Al termine di questa missione secondaria, sarete in grado di:

- Suddividere pipeline complesse in unità logiche e riutilizzabili
- Testare ogni modulo workflow in modo indipendente
- Combinare e abbinare workflows per creare nuove pipeline
- Condividere moduli workflow comuni tra diverse pipeline
- Rendere il codice più manutenibile e più facile da comprendere

Queste competenze vi aiuteranno a costruire pipeline complesse mantenendo una struttura del codice pulita e manutenibile.

### Prerequisiti

Prima di affrontare questa missione secondaria dovreste:

- Aver completato il tutorial [Hello Nextflow](../hello_nextflow/README.md) o un corso per principianti equivalente.
- Essere a proprio agio con l'utilizzo di concetti e meccanismi di base di Nextflow (processi, canali, operatori, moduli)

---

## 0. Iniziare

#### Aprire il codespace di formazione

Se non lo avete ancora fatto, assicuratevi di aprire l'ambiente di formazione come descritto in [Configurazione dell'Ambiente](../envsetup/index.md).

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### Spostarsi nella directory del progetto

Spostiamoci nella directory dove si trovano i file per questo tutorial.

```bash
cd side-quests/workflows_of_workflows
```

Potete impostare VSCode per focalizzarsi su questa directory:

```bash
code .
```

#### Esaminare i materiali

Troverete una directory `modules` contenente diverse definizioni di processi che si basano su quanto appreso in 'Hello Nextflow':

```console title="Contenuto della directory"
modules/
├── say_hello.nf             # Crea un saluto (da Hello Nextflow)
├── say_hello_upper.nf       # Converte in maiuscolo (da Hello Nextflow)
├── timestamp_greeting.nf    # Aggiunge timestamp ai saluti
├── validate_name.nf         # Valida i nomi di input
└── reverse_text.nf          # Inverte il contenuto del testo
```

#### Esaminare l'assegnazione

La vostra sfida è assemblare questi moduli in due workflows separati che poi comporremo in un workflow principale:

- Un `GREETING_WORKFLOW` che valida i nomi, crea i saluti e aggiunge timestamp
- Un `TRANSFORM_WORKFLOW` che converte il testo in maiuscolo e lo inverte

#### Lista di controllo della preparazione

Pensate di essere pronti ad iniziare?

- [ ] Comprendo l'obiettivo di questo corso e i suoi prerequisiti
- [ ] Il mio codespace è attivo e funzionante
- [ ] Ho impostato la mia directory di lavoro appropriatamente
- [ ] Comprendo l'assegnazione

Se potete spuntare tutte le caselle, siete pronti per iniziare.

---

## 1. Creare il Greeting Workflow

Iniziamo creando un workflow che valida i nomi e genera saluti con timestamp.

### 1.1. Creare la struttura del workflow

```bash title="Creare la directory e il file del workflow"
mkdir -p workflows
touch workflows/greeting.nf
```

### 1.2. Aggiungere il codice del primo (sub)workflow

Aggiungete questo codice a `workflows/greeting.nf`:

```groovy title="workflows/greeting.nf" linenums="1"
include { VALIDATE_NAME } from '../modules/validate_name'
include { SAY_HELLO } from '../modules/say_hello'
include { TIMESTAMP_GREETING } from '../modules/timestamp_greeting'

workflow {

    names_ch = channel.of('Alice', 'Bob', 'Charlie')

    // Concatenare i processi: validate -> create greeting -> add timestamp
    validated_ch = VALIDATE_NAME(names_ch)
    greetings_ch = SAY_HELLO(validated_ch)
    timestamped_ch = TIMESTAMP_GREETING(greetings_ch)
}
```

Questo è un workflow completo, con una struttura simile a quelli visti nel tutorial 'Hello Nextflow', che possiamo testare in modo indipendente. Proviamolo ora:

```bash
nextflow run workflows/greeting.nf
```

??? success "Output del comando"

    ```console
    N E X T F L O W  ~  version 24.10.0
    Launching `workflows/greeting.nf` [peaceful_montalcini] DSL2 - revision: 90f61b7093
    executor >  local (9)
    [51/4f980f] process > VALIDATE_NAME (validating Bob)                    [100%] 3 of 3 ✔
    [2b/dd8dc2] process > SAY_HELLO (greeting Bob)                          [100%] 3 of 3 ✔
    [8e/882565] process > TIMESTAMP_GREETING (adding timestamp to greeting) [100%] 3 of 3 ✔
    ```

Funziona come previsto, ma per renderlo componibile ci sono alcune cose che dobbiamo modificare.

### 1.3. Rendere il workflow componibile

I workflows componibili presentano alcune differenze rispetto a quelli visti nel tutorial 'Hello Nextflow':

- Il blocco workflow deve essere nominato
- Gli input sono dichiarati utilizzando la parola chiave `take:`
- Il contenuto del workflow è posto all'interno del blocco `main:`
- Gli output sono dichiarati utilizzando la parola chiave `emit:`

Aggiorniamo il greeting workflow per corrispondere a questa struttura. Modificate il codice nel seguente modo:

```groovy title="workflows/greeting.nf" linenums="1" hl_lines="6 7 9 15 16 17"
include { VALIDATE_NAME } from '../modules/validate_name'
include { SAY_HELLO } from '../modules/say_hello'
include { TIMESTAMP_GREETING } from '../modules/timestamp_greeting'

workflow GREETING_WORKFLOW {
    take:
        names_ch        // Canale di input con i nomi

    main:
        // Concatenare i processi: validate -> create greeting -> add timestamp
        validated_ch = VALIDATE_NAME(names_ch)
        greetings_ch = SAY_HELLO(validated_ch)
        timestamped_ch = TIMESTAMP_GREETING(greetings_ch)

    emit:
        greetings = greetings_ch      // Saluti originali
        timestamped = timestamped_ch  // Saluti con timestamp
}
```

Potete vedere che il workflow è ora nominato e ha un blocco `take:` e `emit:`, e queste sono le connessioni che useremo per comporre un workflow di livello superiore.
Il contenuto del workflow è anche posto all'interno del blocco `main:`. Notate inoltre che abbiamo rimosso la dichiarazione del canale di input `names_ch`, poiché ora viene passato come argomento al workflow.

Testiamo di nuovo il workflow per vedere se funziona come previsto:

```bash
nextflow run workflows/greeting.nf
```

??? failure "Output del comando"

    ```console
    N E X T F L O W  ~  version 24.10.0
    Launching `workflows/greeting.nf` [high_brahmagupta] DSL2 - revision: 8f5857af25
    No entry workflow specified
    ```

Questo vi informa di un altro nuovo concetto, un 'entry workflow'. L'entry workflow è il workflow che viene chiamato quando si esegue uno script Nextflow. Per impostazione predefinita, Nextflow utilizzerà un workflow senza nome come entry workflow, quando presente, ed è ciò che avete fatto finora, con blocchi workflow che iniziano così:

```groovy title="hello.nf" linenums="1"
workflow {
```

Ma il nostro greeting workflow non ha un workflow senza nome, piuttosto abbiamo un workflow nominato:

```groovy title="workflows/greeting.nf" linenums="1"
workflow GREETING_WORKFLOW {
```

Ecco perché Nextflow ha generato un errore e non ha fatto ciò che volevamo.

Non abbiamo aggiunto la sintassi `take:`/`emit:` per poter chiamare direttamente il workflow - l'abbiamo fatto per poterlo comporre con altri workflows. La soluzione è creare uno script principale con un entry workflow senza nome che importa e chiama il nostro workflow nominato.

### 1.4. Creare e testare il workflow principale

Ora creeremo un workflow principale che importa e utilizza il workflow `greeting`.

Create `main.nf`:

```groovy title="main.nf" linenums="1"
include { GREETING_WORKFLOW } from './workflows/greeting'

workflow {
    names = channel.of('Alice', 'Bob', 'Charlie')
    GREETING_WORKFLOW(names)

    GREETING_WORKFLOW.out.greetings.view { "Original: $it" }
    GREETING_WORKFLOW.out.timestamped.view { "Timestamped: $it" }
}

```

Notate che la voce workflow in questo file è senza nome, e questo perché la useremo come entry workflow.

Eseguite questo e osservate l'output:

```bash
nextflow run main.nf
```

??? success "Output del comando"

    ```console
    N E X T F L O W  ~  version 24.10.0
    Launching `main.nf` [goofy_mayer] DSL2 - revision: 543f8742fe
    executor >  local (9)
    [05/3cc752] process > GREETING_WORKFLOW:VALIDATE_NAME (validating Char... [100%] 3 of 3 ✔
    [b1/b56ecf] process > GREETING_WORKFLOW:SAY_HELLO (greeting Charlie)      [100%] 3 of 3 ✔
    [ea/342168] process > GREETING_WORKFLOW:TIMESTAMP_GREETING (adding tim... [100%] 3 of 3 ✔
    Original: /workspaces/training/side_quests/workflows_of_workflows/work/bb/c8aff3df0ebc15a4d7d35f736db44c/Alice-output.txt
    Original: /workspaces/training/side_quests/workflows_of_workflows/work/fb/fa877776e8a5d90b537b1bcd3b6f5b/Bob-output.txt
    Original: /workspaces/training/side_quests/workflows_of_workflows/work/b1/b56ecf938fda8bcbec211847c8f0be/Charlie-output.txt
    Timestamped: /workspaces/training/side_quests/workflows_of_workflows/work/06/877bc909f140bbf8223343450cea36/timestamped_Alice-output.txt
    Timestamped: /workspaces/training/side_quests/workflows_of_workflows/work/aa/bd31b71cdb745b7c155ca7f8837b8a/timestamped_Bob-output.txt
    Timestamped: /workspaces/training/side_quests/workflows_of_workflows/work/ea/342168d4ba04cc899a89c56cbfd9b0/timestamped_Charlie-output.txt
    ```

Funziona! Abbiamo racchiuso il greeting workflow nominato in un workflow principale con un blocco `workflow` di ingresso senza nome. Il workflow principale sta utilizzando il workflow `GREETING_WORKFLOW` quasi (non proprio) come un processo, e sta passando il canale `names` come argomento.

### Conclusioni

In questa sezione, avete appreso diversi concetti importanti:

- **Workflows Nominati**: Creare un workflow nominato (`GREETING_WORKFLOW`) che può essere importato e riutilizzato
- **Interfacce di Workflow**: Definire input chiari con `take:` e output con `emit:` per creare un workflow componibile
- **Punti di Ingresso**: Comprendere che Nextflow necessita di un entry workflow senza nome per eseguire uno script
- **Composizione di Workflow**: Importare e utilizzare un workflow nominato all'interno di un altro workflow
- **Namespace di Workflow**: Accedere agli output del workflow utilizzando il namespace `.out` (`GREETING_WORKFLOW.out.greetings`)

Ora avete un greeting workflow funzionante che:

- Prende un canale di nomi come input
- Valida ogni nome
- Crea un saluto per ogni nome valido
- Aggiunge timestamp ai saluti
- Espone sia i saluti originali che quelli con timestamp come output

Questo approccio modulare vi consente di testare il greeting workflow in modo indipendente o di utilizzarlo come componente in pipeline più grandi.

---

## 2. Aggiungere il Transform Workflow

Ora creiamo un workflow che applica trasformazioni di testo ai saluti.

### 2.1. Creare il file del workflow

```bash
touch workflows/transform.nf
```

### 2.2. Aggiungere il codice del workflow

Aggiungete questo codice a `workflows/transform.nf`:

```groovy title="workflows/transform.nf" linenums="1"
include { SAY_HELLO_UPPER } from '../modules/say_hello_upper'
include { REVERSE_TEXT } from '../modules/reverse_text'

workflow TRANSFORM_WORKFLOW {
    take:
        input_ch         // Canale di input con i messaggi

    main:
        // Applicare le trasformazioni in sequenza
        upper_ch = SAY_HELLO_UPPER(input_ch)
        reversed_ch = REVERSE_TEXT(upper_ch)

    emit:
        upper = upper_ch        // Saluti in maiuscolo
        reversed = reversed_ch  // Saluti in maiuscolo invertiti
}
```

Non ripeteremo qui la spiegazione della sintassi componibile, ma notate che il workflow nominato è nuovamente dichiarato con un blocco `take:` e `emit:`, e il contenuto del workflow è posto all'interno del blocco `main:`.

### 2.3. Aggiornare il workflow principale

Aggiornate `main.nf` per utilizzare entrambi i workflows:

```groovy title="main.nf" linenums="1"
include { GREETING_WORKFLOW } from './workflows/greeting'
include { TRANSFORM_WORKFLOW } from './workflows/transform'

workflow {
    names = channel.of('Alice', 'Bob', 'Charlie')

    // Eseguire il greeting workflow
    GREETING_WORKFLOW(names)

    // Eseguire il transform workflow
    TRANSFORM_WORKFLOW(GREETING_WORKFLOW.out.timestamped)

    // Visualizzare i risultati
    TRANSFORM_WORKFLOW.out.upper.view { "Uppercase: $it" }
    TRANSFORM_WORKFLOW.out.reversed.view { "Reversed: $it" }
}
```

Eseguite la pipeline completa:

```bash
nextflow run main.nf
```

??? success "Output del comando"

    ```console
    N E X T F L O W  ~  version 24.10.0
    Launching `main.nf` [sick_kimura] DSL2 - revision: 8dc45fc6a8
    executor >  local (13)
    executor >  local (15)
    [83/1b51f4] process > GREETING_WORKFLOW:VALIDATE_NAME (validating Alice)  [100%] 3 of 3 ✔
    [68/556150] process > GREETING_WORKFLOW:SAY_HELLO (greeting Alice)        [100%] 3 of 3 ✔
    [de/511abd] process > GREETING_WORKFLOW:TIMESTAMP_GREETING (adding tim... [100%] 3 of 3 ✔
    [cd/e6a7e0] process > TRANSFORM_WORKFLOW:SAY_HELLO_UPPER (converting t... [100%] 3 of 3 ✔
    [f0/74ba4a] process > TRANSFORM_WORKFLOW:REVERSE_TEXT (reversing UPPER... [100%] 3 of 3 ✔
    Uppercase: /workspaces/training/side_quests/workflows_of_workflows/work/a0/d4f5df4d6344604498fa47a6084a11/UPPER-timestamped_Bob-output.txt
    Uppercase: /workspaces/training/side_quests/workflows_of_workflows/work/69/b5e37f6c79c2fd38adb75d0eca8f87/UPPER-timestamped_Charlie-output.txt
    Uppercase: /workspaces/training/side_quests/workflows_of_workflows/work/cd/e6a7e0b17e7d5a2f71bb8123cd53a7/UPPER-timestamped_Alice-output.txt
    Reversed: /workspaces/training/side_quests/workflows_of_workflows/work/7a/7a222f7957b35d1d121338566a24ac/REVERSED-UPPER-timestamped_Bob-output.txt
    Reversed: /workspaces/training/side_quests/workflows_of_workflows/work/46/8d19af62e33a5a6417c773496e0f90/REVERSED-UPPER-timestamped_Charlie-output.txt
    Reversed: /workspaces/training/side_quests/workflows_of_workflows/work/f0/74ba4a10d9ef5c82f829d1c154d0f6/REVERSED-UPPER-timestamped_Alice-output.txt
    ```

Se date un'occhiata a uno di quei file invertiti, vedrete che è la versione in maiuscolo del saluto invertita:

```bash
cat /workspaces/training/side_quests/workflows_of_workflows/work/f0/74ba4a10d9ef5c82f829d1c154d0f6/REVERSED-UPPER-timestamped_Alice-output.txt
```

```console title="Contenuto del file invertito"
!ECILA ,OLLEH ]04:50:71 60-30-5202[
```

### Conclusioni

Ora dovreste avere una pipeline completa che:

- Elabora i nomi attraverso il greeting workflow
- Invia i saluti con timestamp al transform workflow
- Produce sia versioni in maiuscolo che invertite dei saluti

---

## Riepilogo

In questa missione secondaria, abbiamo esplorato il potente concetto di composizione di workflow in Nextflow, che ci consente di costruire pipeline complesse da componenti più piccoli e riutilizzabili.

Questo approccio modulare offre diversi vantaggi rispetto alle pipeline monolitiche:

- Ogni workflow può essere sviluppato, testato e sottoposto a debug in modo indipendente
- I workflows possono essere riutilizzati in diverse pipeline
- La struttura complessiva della pipeline diventa più leggibile e manutenibile
- Le modifiche a un workflow non influenzano necessariamente gli altri se le interfacce rimangono coerenti
- I punti di ingresso possono essere configurati per eseguire diverse parti della pipeline secondo necessità

_È importante notare tuttavia che, sebbene chiamare i workflows sia un po' come chiamare i processi, non è effettivamente la stessa cosa. Non è possibile, ad esempio, eseguire un workflow N volte chiamandolo con un canale di dimensione N - sarebbe necessario passare un canale di dimensione N al workflow e iterare internamente._

L'applicazione di queste tecniche nel vostro lavoro vi consentirà di costruire pipeline Nextflow più sofisticate in grado di gestire attività bioinformatiche complesse rimanendo manutenibili e scalabili.

### Modelli chiave

1.  **Struttura del workflow**: Abbiamo definito input e output chiari per ogni workflow utilizzando la sintassi `take:` e `emit:`, creando interfacce ben definite tra i componenti, e racchiuso la logica del workflow all'interno del blocco `main:`.

    ```groovy
    workflow EXAMPLE_WORKFLOW {
        take:
            // I canali di input sono dichiarati qui
            input_ch

        main:
            // La logica del workflow va qui
            // Qui è dove i processi vengono chiamati e i canali vengono manipolati
            result_ch = SOME_PROCESS(input_ch)

        emit:
            // I canali di output sono dichiarati qui
            output_ch = result_ch
    }
    ```

2.  **Importazioni di workflow:** Abbiamo costruito due moduli workflow indipendenti e li abbiamo importati in una pipeline principale con istruzioni include.

    - Includere un singolo workflow

    ```groovy
    include { WORKFLOW_NAME } from './path/to/workflow'
    ```

    - Includere più workflows

    ```groovy
    include { WORKFLOW_A; WORKFLOW_B } from './path/to/workflows'
    ```

    - Includere con alias per evitare conflitti di nomi

    ```groovy
    include { WORKFLOW_A as WORKFLOW_A_ALIAS } from './path/to/workflow'
    ```

3.  **Punti di ingresso**: Nextflow richiede un entry workflow senza nome per sapere da dove iniziare l'esecuzione. Questo entry workflow chiama i vostri workflows nominati.

    - Workflow senza nome (punto di ingresso)

    ```groovy
    workflow {
        // Questo è il punto di ingresso quando lo script viene eseguito
        NAMED_WORKFLOW(input_ch)
    }
    ```

    - Workflow nominato (chiamato dall'entry workflow)

    ```groovy
    workflow NAMED_WORKFLOW {
        // Deve essere chiamato dall'entry workflow
    }
    ```

4.  **Gestione del flusso di dati:** Abbiamo imparato come accedere agli output del workflow utilizzando la notazione namespace (`WORKFLOW_NAME.out.channel_name`) e passarli ad altri workflows.

    ```nextflow
    WORKFLOW_A(input_ch)
    WORKFLOW_B(WORKFLOW_A.out.some_channel)
    ```

### Risorse aggiuntive

- [Documentazione Nextflow Workflow](https://www.nextflow.io/docs/latest/workflow.html)
- [Riferimento Channel Operators](https://www.nextflow.io/docs/latest/operator.html)
- [Documentazione Error Strategy](https://www.nextflow.io/docs/latest/process.html#errorstrategy)

---

## Cosa c'è dopo?

Tornate al [menu delle Missioni Secondarie](./index.md) o cliccate il pulsante in basso a destra della pagina per passare all'argomento successivo nell'elenco.
