# Part 2: Executar pipelines reals

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducció assistida per IA - [més informació i suggeriments](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

A la Part 1 d'aquest curs (Executar Operacions Bàsiques), vam començar amb un exemple de workflow que només tenia funcionalitats mínimes per mantenir la complexitat del codi baixa.
Per exemple, `1-hello.nf` utilitzava un paràmetre de línia de comandes (`--input`) per proporcionar un sol valor cada vegada.

No obstant això, la majoria de pipelines del món real utilitzen funcionalitats més sofisticades per permetre el processament eficient de grans quantitats de dades a escala, i apliquen múltiples passos de processament encadenats per una lògica de vegades complexa.

En aquesta part de la formació, demostrem funcionalitats clau de pipelines del món real provant versions ampliades del pipeline original Hello World.

## 1. Processar dades d'entrada des d'un fitxer

En un pipeline del món real, normalment volem processar múltiples punts de dades (o sèries de dades) continguts en un o més fitxers d'entrada.
I sempre que sigui possible, volem executar el processament de dades independents en paral·lel, per escurçar el temps d'espera per a l'anàlisi.

Per demostrar com Nextflow fa això, hem preparat un fitxer CSV anomenat `greetings.csv` que conté diverses salutacions d'entrada, imitant el tipus de dades en columnes que potser voldries processar en una anàlisi de dades real.
Tingues en compte que els números no són significatius, només hi són amb finalitats il·lustratives.

```csv title="data/greetings.csv" linenums="1"
Hello,English,123
Bonjour,French,456
Holà,Spanish,789
```

També hem escrit una versió millorada del workflow original, ara anomenada `2a-inputs.nf`, que llegirà el fitxer CSV, extraurà les salutacions i escriurà cadascuna d'elles en un fitxer separat.

<figure class="excalidraw">
--8<-- "docs/en/docs/nextflow_run/img/hello-pipeline-multi-inputs.svg"
</figure>

Executem primer el workflow, i després mirarem el codi Nextflow rellevant.

### 1.1. Executar el workflow

Executa la següent comanda al teu terminal.

```bash
nextflow run 2a-inputs.nf --input data/greetings.csv
```

??? success "Sortida de la comanda"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `2a-inputs.nf` [mighty_sammet] DSL2 - revision: 29fb5352b3

    executor >  local (3)
    [8e/0eb066] sayHello (2) [100%] 3 of 3 ✔
    ```

Emocionantment, això sembla indicar que es van fer '3 de 3' crides per al procés, cosa que és encoratjadora, ja que hi havia tres files de dades al CSV que vam proporcionar com a entrada.
Això suggereix que el procés `sayHello()` es va cridar tres vegades, una per cada fila d'entrada.

### 1.2. Trobar les sortides publicades al directori `results`

Mirem el directori 'results' per veure si el nostre workflow encara està escrivint una còpia de les nostres sortides allà.

??? abstract "Contingut del directori"

    ```console linenums="1" hl_lines="4-7"
    results
    ├── 1-hello
    |   └── output.txt
    └── 2a-inputs
        ├── Bonjour-output.txt
        ├── Hello-output.txt
        └── Holà-output.txt
    ```

Sí! Veiem un nou directori anomenat `2a-inputs` amb tres fitxers de sortida amb noms diferents, prou convenientment.

Pots obrir cadascun d'ells per assegurar-te que contenen la cadena de salutació apropiada.

??? abstract "Contingut del fitxer"

    ```console title="results/2a-inputs/Hello-output.txt"
    Hello
    ```

    ```console title="results/2a-inputs/Bonjour-output.txt"
    Bonjour
    ```

    ```console title="results/2a-inputs/Holà-output.txt"
    Holà
    ```

Això confirma que cada salutació del fitxer d'entrada s'ha processat adequadament.

### 1.3. Trobar les sortides originals i els registres

Potser has notat que la sortida de consola anterior només feia referència a un directori de tasca.
Això vol dir que les tres crides a `sayHello()` es van executar dins d'aquell únic directori de tasca?

#### 1.3.1. Examinar el directori de tasca donat al terminal

Donem una ullada dins d'aquell directori de tasca `8e/0eb066`.

??? abstract "Contingut del directori"

    ```console title="8e/0eb066"
    work/8e/0eb066071cdb4123906b7b4ea8b047/
    └── Bonjour-output.txt
    ```

Només trobem la sortida corresponent a una de les salutacions (així com els fitxers accessoris si habilitem la visualització de fitxers ocults).

Aleshores, què està passant aquí?

Per defecte, el sistema de registre ANSI escriu la informació d'estat per a totes les crides al mateix procés a la mateixa línia.
Com a resultat, només ens va mostrar un dels tres camins de directori de tasca (`8e/0eb066`) a la sortida de consola.
Hi ha dos altres que no es llisten allà.

#### 1.3.2. Fer que el terminal mostri més detalls

Podem modificar el comportament de registre per veure la llista completa de crides de procés afegint `-ansi-log false` a la comanda de la següent manera:

```bash
nextflow run 2a-inputs.nf --input data/greetings.csv -ansi-log false
```

??? success "Sortida de la comanda"

    ```console linenums="1"
    N E X T F L O W  ~  version 25.10.2
    Launching `2a-inputs.nf` [pedantic_hamilton] DSL2 - revision: 6bbc42e49f
    [ab/1a8ece] Submitted process > sayHello (1)
    [0d/2cae24] Submitted process > sayHello (2)
    [b5/0df1d6] Submitted process > sayHello (3)
    ```

Aquesta vegada veiem les tres execucions de procés i els seus subdirectoris de treball associats llistats a la sortida.
Deshabilitar el registre ANSI també va impedir que Nextflow utilitzés colors a la sortida del terminal.

Observa que la manera com es reporta l'estat és una mica diferent entre els dos modes de registre.
En el mode condensat, Nextflow informa si les crides es van completar amb èxit o no.
En aquest mode ampliat, només informa que es van enviar.

Això confirma que el procés `sayHello()` es crida tres vegades, i es crea un directori de tasca separat per a cadascuna.

Si mirem dins de cadascun dels directoris de tasca llistats allà, podem verificar que cadascun correspon a una de les salutacions.

??? abstract "Contingut del directori"

    ```console title="ab/1a8ece"
    work/ab/1a8ece307e53f03fce689dde904b64/
    └── Hello-output.txt
    ```

    ```console title="0d/2cae24"
    work/0d/2cae2481a53593bc607077c80c9466/
    └── Bonjour-output.txt
    ```

    ```console title="b5/0df1d6"
    work/b5/0df1d642353269909c2ce23fc2a8fa/
    └── Holà-output.txt
    ```

Això confirma que cada crida de procés s'executa de manera aïllada de totes les altres.
Això té molts avantatges, incloent evitar col·lisions si el procés produeix fitxers intermedis amb noms no únics.

!!! tip "Consell"

    Per a un workflow complex, o un gran nombre d'entrades, tenir la llista completa a la sortida del terminal pot ser una mica aclaparador, així que la gent normalment no utilitza `-ansi-log false` en l'ús rutinari.

### 1.4. Examinar el codi del workflow

Així doncs, aquesta versió del workflow és capaç de llegir un fitxer CSV d'entrades, processar les entrades per separat i nomenar les sortides de manera única.

Donem una ullada al que fa això possible al codi del workflow.

??? full-code "Fitxer de codi complet"

    ```groovy title="2a-inputs.nf" linenums="1" hl_lines="31-33 35"
    #!/usr/bin/env nextflow

    /*
    * Use echo to print 'Hello World!' to a file
    */
    process sayHello {

        input:
        val greeting

        output:
        path "${greeting}-output.txt"

        script:
        """
        echo '${greeting}' > '${greeting}-output.txt'
        """
    }

    /*
    * Pipeline parameters
    */
    params {
        input: Path
    }

    workflow {

        main:
        // crea un canal per a entrades des d'un fitxer CSV
        greeting_ch = channel.fromPath(params.input)
                            .splitCsv()
                            .map { line -> line[0] }
        // emet una salutacio
        sayHello(greeting_ch)

        publish:
        first_output = sayHello.out
    }

    output {
        first_output {
            path '2a-inputs'
            mode 'copy'
        }
    }
    ```

Una vegada més, no necessites memoritzar la sintaxi del codi, però és bo aprendre a reconèixer components clau del workflow que proporcionen funcionalitat important.

#### 1.4.1. Carregar les dades d'entrada des del CSV

Aquesta és la part més interessant: com vam canviar de prendre un sol valor de la línia de comandes, a prendre un fitxer CSV, analitzar-lo i processar les salutacions individuals que conté?

A Nextflow, fem això amb un [**canal**](https://nextflow.io/docs/latest/channel.html): una construcció de cua dissenyada per gestionar entrades de manera eficient i transportar-les d'un pas a un altre en workflows de múltiples passos, mentre proporciona paral·lelisme integrat i molts beneficis addicionals.

Desglossem-ho.

```groovy title="2a-inputs.nf" linenums="29" hl_lines="3-5"
    main:
    // crea un canal per a entrades des d'un fitxer CSV
    greeting_ch = channel.fromPath(params.input)
                        .splitCsv()
                        .map { line -> line[0] }
    // emet una salutacio
    sayHello(greeting_ch)
```

Aquest codi crea un canal anomenat `greeting_ch` que llegeix el fitxer CSV, l'analitza i extreu la primera columna de cada fila.
El resultat és un canal que conté `Hello`, `Bonjour` i `Holà`.

??? tip "Com funciona això?"

    Això és el que significa aquesta línia en llenguatge plà:

    - `channel.fromPath` és una **factoria de canals** que crea un canal a partir de camí(ns) de fitxer
    - `(params.input)` especifica que el camí del fitxer es proporciona amb `--input` a la línia de comandes

    En altres paraules, aquesta línia diu a Nextflow: pren el camí del fitxer donat amb `--input` i prepara't per tractar el seu contingut com a dades d'entrada.

    Després les dues línies següents apliquen **operadors** que fan l'anàlisi real del fitxer i la càrrega de les dades a l'estructura de dades apropiada:

    - `.splitCsv()` diu a Nextflow que analitzi el fitxer CSV en una matriu que representa files i columnes
    - `.map { line -> line[0] }` diu a Nextflow que prengui només l'element de la primera columna de cada fila

    Així que a la pràctica, començant pel següent fitxer CSV:

    ```csv title="greetings.csv" linenums="1"
    Hello,English,123
    Bonjour,French,456
    Holà,Spanish,789
    ```

    Hem transformat això en una matriu que es veu així:

    ```txt title="Array contents"
    [[Hello,English,123],[Bonjour,French,456],[Holà,Spanish,789]]
    ```

    I després hem pres el primer element de cadascuna de les tres files i els hem carregat en un canal Nextflow que ara conté: `Hello`, `Bonjour` i `Holà`.

    Si vols entendre els canals i operadors en profunditat, incloent com escriure'ls tu mateix, consulta [Hello Nextflow Part 2: Hello Channels](../hello_nextflow/02_hello_channels.md#4-read-input-values-from-a-csv-file).

#### 1.4.2. Cridar el procés a cada salutació

A continuació, a l'última línia del bloc `main:` del workflow, proporcionem el canal `greeting_ch` carregat com a entrada al procés `sayHello()`.

```groovy title="2a-inputs.nf" linenums="29" hl_lines="7"
    main:
    // crea un canal per a entrades des d'un fitxer CSV
    greeting_ch = channel.fromPath(params.input)
                        .splitCsv()
                        .map { line -> line[0] }
    // emet una salutacio
    sayHello(greeting_ch)
```

Això diu a Nextflow que executi el procés individualment a cada element del canal, _és a dir_ a cada salutació.
I perquè Nextflow és intel·ligent així, executarà aquestes crides de procés en paral·lel si és possible, depenent de la infraestructura informàtica disponible.

Així és com pots aconseguir un processament eficient i escalable de moltes dades (moltes mostres, o punts de dades, sigui quina sigui la teva unitat de recerca) amb relativament molt poc codi.

#### 1.4.3. Com es nomenen les sortides

Finalment, val la pena donar una ullada ràpida al codi del procés per veure com aconseguim que els fitxers de sortida es nomenin de manera única.

```groovy title="2a-inputs.nf" linenums="6" hl_lines="7 11"
process sayHello {

    input:
    val greeting

    output:
    path "${greeting}-output.txt"

    script:
    """
    echo '${greeting}' > '${greeting}-output.txt'
    """
}
```

Veus que, comparat amb la versió d'aquest procés a `1-hello.nf`, la declaració de sortida i la part rellevant de la comanda han canviat per incloure el valor de la salutació al nom del fitxer de sortida.

Aquesta és una manera d'assegurar que els noms dels fitxers de sortida no col·lisionaran quan es publiquin al directori de resultats comú.

I aquest és l'únic canvi que hem hagut de fer dins de la declaració del procés!

### Conclusió

Entens a un nivell bàsic com els canals i operadors ens permeten processar múltiples entrades de manera eficient.

### Què segueix?

Descobreix com es construeixen els workflows de múltiples passos i com operen.

---

## 2. Executar workflows de múltiples passos

La majoria de workflows del món real impliquen més d'un pas.
Construïm sobre el que acabem d'aprendre sobre canals, i mirem com Nextflow utilitza canals i operadors per connectar processos junts en un workflow de múltiples passos.

Amb aquesta finalitat, et proporcionem un exemple de workflow que encadena tres passos separats i demostra el següent:

1. Fer que les dades flueixin d'un procés al següent
2. Recollir sortides de múltiples crides de procés en una sola crida de procés

Específicament, hem fet una versió ampliada del workflow anomenada `2b-multistep.nf` que pren cada salutació d'entrada, la converteix a majúscules, després recull totes les salutacions en majúscules en un sol fitxer de sortida.

<figure class="excalidraw">
--8<-- "docs/en/docs/nextflow_run/img/hello-pipeline-multi-steps.svg"
</figure>

Com abans, executarem primer el workflow i després mirarem el codi per veure què hi ha de nou.

### 2.1. Executar el workflow

Executa la següent comanda al teu terminal:

```bash
nextflow run 2b-multistep.nf --input data/greetings.csv
```

??? success "Sortida de la comanda"

    ```console linenums="1"
    N E X T F L O W   ~  version 25.10.2

    Launching `2b-multistep.nf` [soggy_franklin] DSL2 - revision: bc8e1b2726

    [d6/cdf466] sayHello (1)       | 3 of 3 ✔
    [99/79394f] convertToUpper (2) | 3 of 3 ✔
    [1e/83586c] collectGreetings   | 1 of 1 ✔
    ```

Veus que com es va prometre, es van executar múltiples passos com a part del workflow; els dos primers (`sayHello` i `convertToUpper`) presumiblement es van executar a cada salutació individual, i el tercer (`collectGreetings`) s'haurà executat només una vegada, a les sortides de les tres crides de `convertToUpper`.

### 2.2. Trobar les sortides

Verifiquem que això és de fet el que va passar mirant al directori `results`.

??? abstract "Contingut del directori"

    ```console linenums="1" hl_lines="8-16"
    results
    ├── 1-hello
    |   └── output.txt
    ├── 2a-inputs
    |   ├── Bonjour-output.txt
    |   ├── Hello-output.txt
    |   └── Holà-output.txt
    └── 2b-multistep
        ├── COLLECTED-batch-output.txt
        ├── batch-report.txt
        └── intermediates
            ├── Bonjour-output.txt
            ├── Hello-output.txt
            ├── Holà-output.txt
            ├── UPPER-Bonjour-output.txt
            ├── UPPER-Hello-output.txt
            └── UPPER-Holà-output.txt

    ```

Com pots veure, tenim un nou directori anomenat `2b-multistep`, i conté bastants més fitxers que abans.
Alguns dels fitxers s'han agrupat en un subdirectori anomenat `intermediates`, mentre que dos fitxers es troben al nivell superior.

Aquests dos són els resultats finals del workflow de múltiples passos.
Pren-te un minut per mirar els noms dels fitxers i comprovar el seu contingut per confirmar que són el que esperes.

??? abstract "Contingut del fitxer"

    ```txt title="results/2b-multistep/COLLECTED-batch-output.txt"
    HELLO
    BONJOUR
    HOLà
    ```

    ```txt title="results/2b-multistep/batch-report.txt"
    There were 3 greetings in this batch.
    ```

El primer conté les nostres tres salutacions, en majúscules i recollides de nou en un sol fitxer com es va prometre.
El segon és un fitxer d'informe que resumeix alguna informació sobre l'execució.

### 2.3. Examinar el codi

Mirem el codi i identifiquem els patrons clau per a workflows de múltiples passos.

??? full-code "Fitxer de codi complet"

    ```groovy title="2b-multistep.nf" linenums="1" hl_lines="63 75-78 82-84"
    #!/usr/bin/env nextflow

    /*
    * Use echo to print 'Hello World!' to a file
    */
    process sayHello {

        input:
        val greeting

        output:
        path "${greeting}-output.txt"

        script:
        """
        echo '${greeting}' > '${greeting}-output.txt'
        """
    }

    /*
    * Use a text replacement tool to convert the greeting to uppercase
    */
    process convertToUpper {

        input:
        path input_file

        output:
        path "UPPER-${input_file}"

        script:
        """
        cat '${input_file}' | tr '[a-z]' '[A-Z]' > 'UPPER-${input_file}'
        """
    }

    /*
    * Collect uppercase greetings into a single output file
    */
    process collectGreetings {

        input:
        path input_files
        val batch_name

        output:
        path "COLLECTED-${batch_name}-output.txt", emit: outfile
        path "${batch_name}-report.txt", emit: report

        script:
        count_greetings = input_files.size()
        """
        cat ${input_files} > 'COLLECTED-${batch_name}-output.txt'
        echo 'There were ${count_greetings} greetings in this batch.' > '${batch_name}-report.txt'
        """
    }

    /*
    * Pipeline parameters
    */
    params {
        input: Path
        batch: String = 'batch'
    }

    workflow {

        main:
        // crea un canal per a entrades des d'un fitxer CSV
        greeting_ch = channel.fromPath(params.input)
                            .splitCsv()
                            .map { line -> line[0] }
        // emet una salutacio
        sayHello(greeting_ch)
        // converteix la salutacio a majuscules
        convertToUpper(sayHello.out)
        // recull totes les salutacions en un fitxer
        collectGreetings(convertToUpper.out.collect(), params.batch)

        publish:
        first_output = sayHello.out
        uppercased = convertToUpper.out
        collected = collectGreetings.out.outfile
        batch_report = collectGreetings.out.report
    }

    output {
        first_output {
            path '2b-multistep/intermediates'
            mode 'copy'
        }
        uppercased {
            path '2b-multistep/intermediates'
            mode 'copy'
        }
        collected {
            path '2b-multistep'
            mode 'copy'
        }
        batch_report {
            path '2b-multistep'
            mode 'copy'
        }
    }
    ```

Hi ha molt passant allà, però la diferència més òbvia comparada amb la versió anterior del workflow és que ara hi ha múltiples definicions de procés, i corresponentment, diverses crides de procés al bloc workflow.

Donem una ullada més de prop i vegem si podem identificar les peces més interessants.

#### 2.3.1. Visualitzar l'estructura del workflow

Si estàs utilitzant VSCode amb l'extensió Nextflow, pots obtenir un diagrama útil de com els processos estan connectats fent clic a l'enllaç petit `DAG preview` mostrat just sobre el bloc workflow en qualsevol script Nextflow.

<figure class="excalidraw">
--8<-- "docs/en/docs/nextflow_run/img/DAG-multistep.svg"
</figure>

Això et dóna una bona visió general de com els processos estan connectats i què produeixen.

Veus que a més del procés original `sayHello`, ara també tenim `convertToUpper` i `collectGreetings`, que coincideixen amb els noms dels processos que vam veure a la sortida de consola.
Les dues noves definicions de procés estan estructurades de la mateixa manera que el procés `sayHello`, excepte que `collectGreetings` pren un paràmetre d'entrada addicional anomenat `batch` i produeix dues sortides.

No entrarem en detall al codi de cadascun, però si tens curiositat, pots consultar els detalls a la [Part 2 de Hello Nextflow](../hello_nextflow/03_hello_workflow.md).

Per ara, aprofundim en com els processos estan connectats entre si.

#### 2.3.2. Com es connecten els processos

El que realment és interessant mirar aquí és com les crides de procés s'encadenen juntes al bloc `main:` del workflow.

```groovy title="2b-multistep.nf" linenums="68" hl_lines="9 11"
    main:
    // crea un canal per a entrades des d'un fitxer CSV
    greeting_ch = channel.fromPath(params.input)
                        .splitCsv()
                        .map { line -> line[0] }
    // emet una salutacio
    sayHello(greeting_ch)
    // converteix la salutacio a majuscules
    convertToUpper(sayHello.out)
    // recull totes les salutacions en un fitxer
    collectGreetings(convertToUpper.out.collect(), params.batch)
```

Pots veure que la primera crida de procés, `sayHello(greeting_ch)`, no ha canviat.
Després la següent crida de procés, a `convertToUpper`, fa referència a la sortida de `sayHello` com `sayHello.out`.

El patró és simple: `processName.out` fa referència al canal de sortida d'un procés, que es pot passar directament al següent procés.
Així és com transportem dades d'un pas al següent a Nextflow.

#### 2.3.3. Un procés pot prendre múltiples entrades

La tercera crida de procés, a `collectGreetings`, és una mica diferent.

```groovy title="2b-multistep.nf" linenums="77"
    // recull totes les salutacions en un fitxer
    collectGreetings(convertToUpper.out.collect(), params.batch)
```

Veus que aquesta crida rep dues entrades, `convertToUpper.out.collect()` i `params.batch`.
Ignorant la part `.collect()` per ara, podem generalitzar això com `collectGreetings(input1, input2)`.

Això coincideix amb les dues declaracions d'entrada al mòdul de procés:

```groovy title="2b-multistep.nf" linenums="40"
process collectGreetings {

    input:
    path input_files
    val batch_name
```

Quan Nextflow analitza això, assignarà la primera entrada de la crida a `path input_files`, i la segona a `val batch_name`.

Així que ara saps que un procés pot prendre múltiples entrades, i com es veu la crida al bloc workflow.

Ara donem una ullada més de prop a aquella primera entrada, `convertToUpper.out.collect()`.

#### 2.3.4. Què fa `collect()` a la crida `collectGreetings`

Per passar la sortida de `sayHello` a `convertToUpper`, simplement vam fer referència al canal de sortida de `sayHello` com `sayHello.out`. Però per al següent pas, estem veient una referència a `convertToUpper.out.collect()`.

Què és aquesta part `collect()` i què fa?

És un operador, per descomptat. Igual que els operadors `splitCsv` i `map` que vam trobar abans.
Aquesta vegada l'operador s'anomena `collect`, i s'aplica al canal de sortida produït per `convertToUpper`.

L'operador `collect` s'utilitza per recollir les sortides de múltiples crides al mateix procés i empaquetar-les en un sol element de canal.

En el context d'aquest workflow, està prenent les tres salutacions en majúscules al canal `convertToUpper.out` (que són tres elements de canal separats, i normalment serien gestionats en crides separades pel següent procés) i empaquetant-les en un sol element.
Així és com aconseguim que totes les salutacions tornin al mateix fitxer.

<figure class="excalidraw">
--8<-- "docs/en/docs/nextflow_run/img/with-collect-operator.svg"
</figure>

En contrast, si no apliquéssim `collect()` a la sortida de `convertToUpper()` abans d'alimentar-la a `collectGreetings()`, Nextflow simplement executaria `collectGreetings()` independentment a cada salutació, cosa que no aconseguiria el nostre objectiu.

<figure class="excalidraw">
--8<-- "docs/en/docs/nextflow_run/img/without-collect-operator.svg"
</figure>

Hi ha molts altres [operadors](https://nextflow.io/docs/latest/reference/operator.html) disponibles per aplicar transformacions al contingut dels canals entre crides de procés.

Això dóna als desenvolupadors de pipelines molta flexibilitat per personalitzar la lògica de flux del seu pipeline.
L'inconvenient és que de vegades pot fer més difícil desxifrar què està fent el pipeline.

#### 2.3.5. Un paràmetre d'entrada pot tenir un valor per defecte

Potser has notat que `collectGreetings` pren una segona entrada, `params.batch`:

```groovy title="2b-multistep.nf" linenums="77"
    // recull totes les salutacions en un fitxer
    collectGreetings(convertToUpper.out.collect(), params.batch)
```

Això passa un paràmetre CLI anomenat `--batch` al workflow.
No obstant això, quan vam llançar el workflow abans, no vam especificar un paràmetre `--batch`.

Què està passant aquí?
Dona una ullada al bloc `params`:

```groovy title="2b-multistep.nf" linenums="61" hl_lines="3"
params {
    input: Path
    batch: String = 'batch'
}
```

Hi ha un valor per defecte configurat al workflow, així que no hem de proporcionar-lo.
Però si en proporcionem un a la línia de comandes, el valor que especifiquem s'utilitzarà en lloc del valor per defecte.

Prova-ho:

```bash
nextflow run 2b-multistep.nf --input data/greetings.csv --batch test
```

??? success "Sortida de la comanda"

    ```console linenums="1"
    N E X T F L O W   ~  version 25.10.2

    Launching `2b-multistep.nf` [soggy_franklin] DSL2 - revision: bc8e1b2726

    [a5/cdff26] sayHello (1)       | 3 of 3 ✔
    [c5/78794f] convertToUpper (2) | 3 of 3 ✔
    [d3/b4d86c] collectGreetings   | 1 of 1 ✔
    ```

Hauries de veure noves sortides finals nomenades amb el teu nom de lot personalitzat.

??? abstract "Contingut del directori"

    ```console linenums="1" hl_lines="10 12"
    results
    ├── 1-hello
    |   └── output.txt
    ├── 2a-inputs
    |   ├── Bonjour-output.txt
    |   ├── Hello-output.txt
    |   └── Holà-output.txt
    └── 2b-multistep
        ├── COLLECTED-batch-output.txt
        ├── COLLECTED-test-output.txt
        ├── batch-report.txt
        ├── test-report.txt
        └── intermediates
            ├── Bonjour-output.txt
            ├── Hello-output.txt
            ├── Holà-output.txt
            ├── UPPER-Bonjour-output.txt
            ├── UPPER-Hello-output.txt
            └── UPPER-Holà-output.txt
    ```

Aquest és un aspecte de la configuració d'entrada, que cobrirem amb més detall a la Part 3, però per ara el que és important és saber que els paràmetres d'entrada poden tenir valors per defecte.

#### 2.3.6. Un procés pot produir múltiples sortides

A la definició del procés `collectGreetings`, veiem les següents declaracions de sortida:

```groovy title="2b-multistep.nf" linenums="46"
    output:
    path "COLLECTED-${batch_name}-output.txt", emit: outfile
    path "${batch_name}-report.txt", emit: report
```

Que després es fan referència pel nom donat amb `emit:` al bloc `publish:`:

```groovy title="2b-multistep.nf" linenums="80" hl_lines="4 5"
    publish:
    first_output = sayHello.out
    uppercased = convertToUpper.out
    collected = collectGreetings.out.outfile
    batch_report = collectGreetings.out.report
```

Això fa que sigui fàcil passar sortides específiques individualment a altres processos al workflow, en combinació amb diversos operadors.

#### 2.3.7. Les sortides publicades es poden organitzar

Al bloc `output`, hem utilitzat camins personalitzats per agrupar resultats intermedis per tal de facilitar la identificació només de les sortides finals del workflow.

```groovy title="2b-multistep.nf" linenums="87" hl_lines="3 7 11 15"
output {
    first_output {
        path '2b-multistep/intermediates'
        mode 'copy'
    }
    uppercased {
        path '2b-multistep/intermediates'
        mode 'copy'
    }
    collected {
        path '2b-multistep'
        mode 'copy'
    }
    batch_report {
        path '2b-multistep'
        mode 'copy'
    }
}
```

Hi ha maneres més sofisticades d'organitzar sortides publicades; tocarem algunes a la part sobre configuració.

!!! tip "Vols aprendre més sobre la construcció de workflows?"

    Per a una cobertura detallada de la construcció de workflows de múltiples passos, consulta [Hello Nextflow Part 3: Hello Workflow](../hello_nextflow/03_hello_workflow.md).

### Conclusió

Entens a un nivell bàsic com es construeixen els workflows de múltiples passos utilitzant canals i operadors i com operen.
També has vist que els processos poden prendre múltiples entrades i produir múltiples sortides, i que aquestes es poden publicar de manera estructurada.

### Què segueix?

Aprèn com els pipelines Nextflow es poden modularitzar per promoure la reutilització de codi i la mantenibilitat.

---

## 3. Executar pipelines modularitzats

Fins ara, tots els workflows que hem mirat han consistit en un sol fitxer de workflow que conté tot el codi rellevant.

No obstant això, els pipelines del món real normalment es beneficien de ser _modularitzats_, és a dir que el codi es divideix en diferents fitxers.
Això pot fer que el seu desenvolupament i manteniment siguin més eficients i sostenibles.

Aquí demostrarem la forma més comuna de modularitat de codi a Nextflow, que és l'ús de **mòduls**.

A Nextflow, un [**mòdul**](https://nextflow.io/docs/latest/module.html) és una sola definició de procés que està encapsulada per si mateixa en un fitxer de codi independent.
Per utilitzar un mòdul en un workflow, només afegeixes una declaració d'importació d'una sola línia al teu fitxer de codi de workflow; després pots integrar el procés al workflow de la mateixa manera que normalment faries.
Això fa possible reutilitzar definicions de procés en múltiples workflows sense produir múltiples còpies del codi.

Fins ara hem estat executant workflows que tenien tots els seus processos inclosos en un fitxer de codi monolític.
Ara veurem com es veu quan els processos s'emmagatzemen en mòduls individuals.

Per descomptat, hem preparat un cop més un workflow adequat per a fins de demostració, anomenat `2c-modules.nf`, juntament amb un conjunt de mòduls ubicats al directori `modules/`.

<figure class="excalidraw">
--8<-- "docs/en/docs/nextflow_run/img/modules.svg"
</figure>

??? abstract "Contingut del directori"

    ```console
    modules/
    ├── collectGreetings.nf
    ├── convertToUpper.nf
    ├── cowpy.nf
    └── sayHello.nf
    ```

Veus que hi ha quatre fitxers Nextflow, cadascun nomenat segons un dels processos.
Pots ignorar el fitxer `cowpy.nf` per ara; arribarem a aquest més tard.

### 3.1. Examinar el codi

Aquesta vegada mirarem primer el codi.
Comença obrint el fitxer de workflow `2c-modules.nf`.

??? full-code "Fitxer de codi complet"

    ```groovy title="2c-modules.nf" linenums="1"
    #!/usr/bin/env nextflow

    // Include modules
    include { sayHello } from './modules/sayHello.nf'
    include { convertToUpper } from './modules/convertToUpper.nf'
    include { collectGreetings } from './modules/collectGreetings.nf'

    /*
    * Pipeline parameters
    */
    params {
        input: Path
        batch: String = 'batch'
    }

    workflow {

        main:
        // crea un canal per a entrades des d'un fitxer CSV
        greeting_ch = channel.fromPath(params.input)
                            .splitCsv()
                            .map { line -> line[0] }
        // emet una salutacio
        sayHello(greeting_ch)
        // converteix la salutacio a majuscules
        convertToUpper(sayHello.out)
        // recull totes les salutacions en un fitxer
        collectGreetings(convertToUpper.out.collect(), params.batch)

        publish:
        first_output = sayHello.out
        uppercased = convertToUpper.out
        collected = collectGreetings.out.outfile
        batch_report = collectGreetings.out.report
    }

    output {
        first_output {
            path '2c-modules/intermediates'
            mode 'copy'
        }
        uppercased {
            path '2c-modules/intermediates'
            mode 'copy'
        }
        collected {
            path '2c-modules'
            mode 'copy'
        }
        batch_report {
            path '2c-modules'
            mode 'copy'
        }
    }
    ```

Veus que la lògica del workflow és exactament la mateixa que a la versió anterior del workflow.
No obstant això, el codi del procés ha desaparegut del fitxer de workflow, i en lloc d'això hi ha declaracions `include` que apunten a fitxers separats sota `modules`.

```groovy title="hello-modules.nf" linenums="3"
// Include modules
include { sayHello } from './modules/sayHello.nf'
include { convertToUpper } from './modules/convertToUpper.nf'
include { collectGreetings } from './modules/collectGreetings.nf'
```

Obre un d'aquests fitxers i trobaràs el codi per al procés corresponent.

??? full-code "Fitxer de codi complet"

    ```groovy title="modules/sayHello.nf" linenums="1"
    #!/usr/bin/env nextflow

    /*
    * Use echo to print 'Hello World!' to a file
    */
    process sayHello {

        input:
        val greeting

        output:
        path "${greeting}-output.txt"

        script:
        """
        echo '${greeting}' > '${greeting}-output.txt'
        """
    }
    ```

Com pots veure, el codi del procés no ha canviat; només s'ha copiat en un fitxer de mòdul individual en lloc d'estar al fitxer de workflow principal.
El mateix s'aplica als altres dos processos.

Així que vegem com es veu executar aquesta nova versió.

### 3.2. Executar el workflow

Executa aquesta comanda al teu terminal, amb la bandera `-resume`:

```bash
nextflow run 2c-modules.nf --input data/greetings.csv -resume
```

??? success "Sortida de la comanda"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `2c-modules.nf` [soggy_franklin] DSL2 - revision: bc8e1b2726

    [d6/cdf466] sayHello (1)       | 3 of 3, cached: 3 ✔
    [99/79394f] convertToUpper (2) | 3 of 3, cached: 3 ✔
    [1e/83586c] collectGreetings   | 1 of 1, cached: 1 ✔
    ```

Notaràs que les execucions de procés es van emmagatzemar totes amb èxit a la memòria cau, és a dir que Nextflow va reconèixer que ja havia fet el treball sol·licitat, encara que el codi s'ha dividit i el fitxer de workflow principal s'ha reanomenat.

Res d'això importa a Nextflow; el que importa és l'script de tasca que es genera una vegada tot el codi s'ha reunit i avaluat.

!!! tip "Consell"

    També és possible encapsular una secció d'un workflow com un 'subworkflow' que es pot importar en un pipeline més gran, però això està fora de l'abast d'aquest curs.

    Pots aprendre més sobre el desenvolupament de workflows composables a la Side Quest sobre [Workflows of Workflows](https://training.nextflow.io/latest/side_quests/workflows_of_workflows/).

### Conclusió

Saps com els processos es poden emmagatzemar en mòduls independents per promoure la reutilització de codi i millorar la mantenibilitat.

### Què segueix?

Aprèn a utilitzar contenidors per gestionar dependències de programari.

---

## 4. Utilitzar programari en contenidors

Fins ara els workflows que hem estat utilitzant com a exemples només necessitaven executar operacions de processament de text molt bàsiques utilitzant eines UNIX disponibles al nostre entorn.

No obstant això, els pipelines del món real normalment requereixen eines i paquets especialitzats que no estan inclosos per defecte a la majoria d'entorns.
Normalment, hauries d'instal·lar aquestes eines, gestionar les seves dependències i resoldre qualsevol conflicte.

Tot això és molt tediós i molest.
Una manera molt millor d'abordar aquest problema és utilitzar **contenidors**.

Un **contenidor** és una unitat de programari lleugera, independent i executable creada a partir d'una **imatge** de contenidor que inclou tot el necessari per executar una aplicació incloent codi, biblioteques del sistema i configuracions.

!!! Tip "Consell"

    Ensenyem això utilitzant la tecnologia [Docker](https://www.docker.com/get-started/), però Nextflow també suporta diverses altres tecnologies de contenidors.
    Pots aprendre més sobre el suport de Nextflow per a contenidors [aquí](https://nextflow.io/docs/latest/container.html).

### 4.1. Utilitzar un contenidor directament

Primer, provem d'interactuar amb un contenidor directament.
Això ajudarà a consolidar la teva comprensió del que són els contenidors abans de començar a utilitzar-los a Nextflow.

#### 4.1.1. Descarregar la imatge del contenidor

Per utilitzar un contenidor, normalment descarregues o "descarregues" una imatge de contenidor d'un registre de contenidors, i després executes la imatge del contenidor per crear una instància de contenidor.

La sintaxi general és la següent:

```bash title="Syntax"
docker pull '<container>'
```

- `docker pull` és la instrucció al sistema de contenidors per descarregar una imatge de contenidor d'un repositori.
- `'<container>'` és l'adreça URI de la imatge del contenidor.

Com a exemple, descarreguem una imatge de contenidor que conté [cowpy](https://github.com/jeffbuttars/cowpy), una implementació python d'una eina anomenada `cowsay` que genera art ASCII per mostrar entrades de text arbitràries d'una manera divertida.

Hi ha diversos repositoris on pots trobar contenidors publicats.
Vam utilitzar el servei [Seqera Containers](https://seqera.io/containers/) per generar aquesta imatge de contenidor Docker a partir del paquet Conda `cowpy`: `'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'`.

Executa la comanda de descàrrega completa:

```bash
docker pull 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
```

??? success "Sortida de la comanda"

    ```console
    Unable to find image 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273' locally
    131d6a1b707a8e65: Pulling from library/cowpy
    dafa2b0c44d2: Pull complete
    dec6b097362e: Pull complete
    f88da01cff0b: Pull complete
    4f4fb700ef54: Pull complete
    92dc97a3ef36: Pull complete
    403f74b0f85e: Pull complete
    10b8c00c10a5: Pull complete
    17dc7ea432cc: Pull complete
    bb36d6c3110d: Pull complete
    0ea1a16bbe82: Pull complete
    030a47592a0a: Pull complete
    622dd7f15040: Pull complete
    895fb5d0f4df: Pull complete
    Digest: sha256:fa50498b32534d83e0a89bb21fec0c47cc03933ac95c6b6587df82aaa9d68db3
    Status: Downloaded newer image for community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273
    community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273
    ```

Això diu al sistema que descarregui la imatge especificada.
Una vegada la descàrrega està completa, tens una còpia local de la imatge del contenidor.

#### 4.1.2. Iniciar el contenidor

Els contenidors es poden executar com una comanda única, però també pots utilitzar-los de manera interactiva, cosa que et dóna un indicador de shell dins del contenidor i et permet jugar amb la comanda.

La sintaxi general és la següent:

```bash title="Syntax"
docker run --rm '<container>' [tool command]
```

- `docker run --rm '<container>'` és la instrucció al sistema de contenidors per iniciar una instància de contenidor a partir d'una imatge de contenidor i executar una comanda en ella.
- `--rm` diu al sistema que tanqui la instància del contenidor després que la comanda s'hagi completat.

Completament assemblada, la comanda d'execució del contenidor es veu així:

```bash
docker run --rm -it 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
```

Executa aquesta comanda, i hauries de veure el teu indicador canviar a alguna cosa com `(base) root@b645838b3314:/tmp#`, que indica que ara estàs dins del contenidor.

Pots verificar això executant `ls` per llistar el contingut del directori:

```bash
ls /
```

??? success "Sortida de la comanda"

    ```console
    bin  boot  dev  etc  home  lib  lib64  media  mnt  opt  proc  root  run  sbin  srv  sys  tmp  usr  var
    ```

Veus que el sistema de fitxers dins del contenidor és diferent del sistema de fitxers al teu sistema amfitrió.

!!! Tip "Consell"

    Quan executes un contenidor, està aïllat del sistema amfitrió per defecte.
    Això significa que el contenidor no pot accedir a cap fitxer al sistema amfitrió tret que explícitament li permetis fer-ho especificant que vols muntar un volum com a part de la comanda `docker run` utilitzant la següent sintaxi:

    ```bash title="Syntax"
    -v <outside_path>:<inside_path>
    ```

    Això efectivament estableix un túnel a través de la paret del contenidor que pots utilitzar per accedir a aquella part del teu sistema de fitxers.

    Això es cobreix amb més detall a la [Part 5 de Hello Nextflow](../hello_nextflow/05_hello_containers.md).

#### 4.1.3. Executar l'eina `cowpy`

Des de dins del contenidor, pots executar la comanda `cowpy` directament.

```bash
cowpy "Hello Containers"
```

??? success "Sortida de la comanda"

    ```console
    ______________________________________________________
    < Hello Containers >
    ------------------------------------------------------
        \   ^__^
          \  (oo)\_______
            (__)\       )\/\
              ||----w |
              ||     ||
    ```

Això produeix art ASCII del personatge de vaca per defecte (o 'cowacter') amb una bombolla de parla que conté el text que vam especificar.

Ara que has provat l'ús bàsic, pots provar de donar-li alguns paràmetres.
Per exemple, la documentació de l'eina diu que podem establir el personatge amb `-c`.

```bash
cowpy "Hello Containers" -c tux
```

??? success "Sortida de la comanda"

    ```console
    __________________
    < Hello Containers >
    ------------------
      \
        \
            .--.
          |o_o |
          |:_/ |
          //   \ \
        (|     | )
        /'\_   _/`\
        \___)=(___/
    ```

Aquesta vegada la sortida d'art ASCII mostra el pingüí de Linux, Tux, perquè vam especificar el paràmetre `-c tux`.

Com que estàs dins del contenidor, pots executar la comanda cowpy tantes vegades com vulguis, variant els paràmetres d'entrada, sense haver de preocupar-te d'instal·lar cap biblioteca al teu sistema mateix.

??? tip "Altres personatges disponibles"

    Utilitza la bandera '-c' per triar un personatge diferent, incloent:

    `beavis`, `cheese`, `daemon`, `dragonandcow`, `ghostbusters`, `kitty`, `moose`, `milk`, `stegosaurus`, `turkey`, `turtle`, `tux`

Sent lliure de jugar amb això.
Quan hagis acabat, surt del contenidor utilitzant la comanda `exit`:

```bash
exit
```

Et trobaràs de nou al teu shell normal.

### 4.2. Utilitzar un contenidor en un workflow

Quan executem un pipeline, volem poder dir a Nextflow quin contenidor utilitzar a cada pas, i importantment, volem que gestioni tot aquell treball que acabem de fer: descarregar el contenidor, iniciar-lo, executar la comanda i desmuntar el contenidor quan hagi acabat.

Bones notícies: això és exactament el que Nextflow farà per nosaltres.
Només necessitem especificar un contenidor per a cada procés.

Per demostrar com funciona això, hem fet una altra versió del nostre workflow que executa `cowpy` al fitxer de salutacions recollides produït al tercer pas.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-pipeline-cowpy.svg"
</figure>

Això hauria de produir un fitxer que conté l'art ASCII amb les tres salutacions a la bombolla de parla.

#### 4.2.1. Examinar el codi

El workflow és molt similar a l'anterior, més el pas extra per executar `cowpy`.

??? full-code "Fitxer de codi complet"

    ```groovy title="2d-container.nf" linenums="1" hl_lines="7 15 32 39 59-62"
    #!/usr/bin/env nextflow

    // Include modules
    include { sayHello } from './modules/sayHello.nf'
    include { convertToUpper } from './modules/convertToUpper.nf'
    include { collectGreetings } from './modules/collectGreetings.nf'
    include { cowpy } from './modules/cowpy.nf'

    /*
    * Pipeline parameters
    */
    params {
        input: Path
        batch: String = 'batch'
        character: String
    }

    workflow {

        main:
        // crea un canal per a entrades des d'un fitxer CSV
        greeting_ch = channel.fromPath(params.input)
                            .splitCsv()
                            .map { line -> line[0] }
        // emet una salutacio
        sayHello(greeting_ch)
        // converteix la salutacio a majuscules
        convertToUpper(sayHello.out)
        // recull totes les salutacions en un fitxer
        collectGreetings(convertToUpper.out.collect(), params.batch)
        // genera art ASCII de les salutacions amb cowpy
        cowpy(collectGreetings.out.outfile, params.character)

        publish:
        first_output = sayHello.out
        uppercased = convertToUpper.out
        collected = collectGreetings.out.outfile
        batch_report = collectGreetings.out.report
        cowpy_art = cowpy.out
    }

    output {
        first_output {
            path '2d-container/intermediates'
            mode 'copy'
        }
        uppercased {
            path '2d-container/intermediates'
            mode 'copy'
        }
        collected {
            path '2d-container/intermediates'
            mode 'copy'
        }
        batch_report {
            path '2d-container'
            mode 'copy'
        }
        cowpy_art {
            path '2d-container'
            mode 'copy'
        }
    }
    ```

Veus que aquest workflow importa un procés `cowpy` d'un fitxer de mòdul, i el crida a la sortida de la crida `collectGreetings()`, més un paràmetre d'entrada anomenat `params.character`.

```groovy title="2d-container.nf" linenums="31"
// genera art ASCII de les salutacions amb cowpy
cowpy(collectGreetings.out.outfile, params.character)
```

El procés `cowpy`, que encapsula la comanda cowpy per generar art ASCII, està definit al mòdul `cowpy.nf`.

??? full-code "Fitxer de codi complet"

    ```groovy title="modules/cowpy.nf" linenums="1"
    #!/usr/bin/env nextflow

    // Generate ASCII art with cowpy (https://github.com/jeffbuttars/cowpy)
    process cowpy {

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'

        input:
        path input_file
        val character

        output:
        path "cowpy-${input_file}"

        script:
        """
        cat ${input_file} | cowpy -c "${character}" > cowpy-${input_file}
        """
    }
    ```

El procés `cowpy` requereix dues entrades: el camí a un fitxer d'entrada que conté el text per posar a la bombolla de parla (`input_file`), i un valor per a la variable character.

Importantment, també inclou la línia `container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'`, que apunta a l'URI del contenidor que vam utilitzar abans.

#### 4.2.2. Comprovar que Docker està habilitat a la configuració

Anticiparem lleugerament la Part 3 d'aquest curs de formació introduint el fitxer de configuració `nextflow.config`, que és una de les principals maneres que Nextflow ofereix per configurar l'execució del workflow.
Quan un fitxer anomenat `nextflow.config` està present al directori actual, Nextflow el carregarà automàticament i aplicarà qualsevol configuració que contingui.

Amb aquesta finalitat, hem inclòs un fitxer `nextflow.config` amb una sola línia de codi que habilita Docker.

```groovy title="nextflow.config" linenums="1"
docker.enabled = true
```

Aquesta configuració diu a Nextflow que utilitzi Docker per a qualsevol procés que especifiqui un contenidor compatible.

!!! tip "Consell"

    És tècnicament possible habilitar l'execució de Docker des de la línia de comandes, per execució, utilitzant el paràmetre `-with-docker <container>`.
    No obstant això, això només ens permet especificar un contenidor per a tot el workflow, mentre que l'enfocament que acabem de mostrar-te ens permet especificar un contenidor diferent per procés.
    Aquest últim és molt millor per a la modularitat, el manteniment del codi i la reproducibilitat.

#### 4.2.3. Executar el workflow

Només per recapitular, això és el que estem a punt d'executar:

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello_pipeline_complete.svg"
</figure>

Creus que funcionarà?

Executem el workflow amb la bandera `-resume`, i especifiquem que volem que el personatge sigui el gall dindi.

```bash
nextflow run 2d-container.nf --input data/greetings.csv --character turkey -resume
```

??? success "Sortida de la comanda"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `2d-container.nf` [elegant_brattain] DSL2 - revision: 028a841db1

    executor >  local (1)
    [95/fa0bac] sayHello (3)       | 3 of 3, cached: 3 ✔
    [92/32533f] convertToUpper (3) | 3 of 3, cached: 3 ✔
    [aa/e697a2] collectGreetings   | 1 of 1, cached: 1 ✔
    [7f/caf718] cowpy              | 1 of 1 ✔
    ```

Els tres primers passos es van emmagatzemar a la memòria cau ja que ja els hem executat abans, però el procés `cowpy` és nou així que aquest realment s'executa.

Pots trobar la sortida del pas `cowpy` al directori `results`.

??? abstract "Contingut del fitxer"

    ```console title="results/2d-container/cowpy-COLLECTED-batch-output.txt"
    _________
    / HOLà    \
    | HELLO   |
    \ BONJOUR /
    ---------
      \                                  ,+*^^*+___+++_
      \                           ,*^^^^              )
        \                       _+*                     ^**+_
        \                    +^       _ _++*+_+++_,         )
                  _+^^*+_    (     ,+*^ ^          \+_        )
                {       )  (    ,(    ,_+--+--,      ^)      ^\
                { (\@)    } f   ,(  ,+-^ __*_*_  ^^\_   ^\       )
              {:;-/    (_+*-+^^^^^+*+*<_ _++_)_    )    )      /
              ( /  (    (        ,___    ^*+_+* )   <    <      \
              U _/     )    *--<  ) ^\-----++__)   )    )       )
                (      )  _(^)^^))  )  )\^^^^^))^*+/    /       /
              (      /  (_))_^)) )  )  ))^^^^^))^^^)__/     +^^
            (     ,/    (^))^))  )  ) ))^^^^^^^))^^)       _)
              *+__+*       (_))^)  ) ) ))^^^^^^))^^^^^)____*^
              \             \_)^)_)) ))^^^^^^^^^^))^^^^)
              (_             ^\__^^^^^^^^^^^^))^^^^^^^)
                ^\___            ^\__^^^^^^))^^^^^^^^)\\
                      ^^^^^\uuu/^^\uuu/^^^^\^\^\^\^\^\^\^\
                        ___) >____) >___   ^\_\_\_\_\_\_\)
                        ^^^//\\_^^//\\_^       ^(\_\_\_\)
                          ^^^ ^^ ^^^ ^
    ```

Veus que el personatge està dient totes les salutacions, ja que es va executar al fitxer de salutacions recollides en majúscules.

Més important encara, vam poder executar això com a part del nostre pipeline sense haver de fer una instal·lació adequada de cowpy i totes les seves dependències.
I ara podem compartir el pipeline amb col·laboradors i fer que l'executin a la seva infraestructura sense que ells necessitin instal·lar res tampoc, a part de Docker o una de les seves alternatives (com Singularity/Apptainer) com s'ha esmentat anteriorment.

#### 4.2.4. Inspeccionar com Nextflow va llançar la tasca en contenidor

Com a coda final d'aquesta secció, donem una ullada al subdirectori de treball per a una de les crides de procés `cowpy` per obtenir una mica més d'informació sobre com Nextflow treballa amb contenidors sota el capó.

Comprova la sortida de la teva comanda `nextflow run` per trobar el camí al subdirectori de treball per al procés `cowpy`.
Mirant el que vam obtenir per a l'execució mostrada anteriorment, la línia de registre de consola per al procés `cowpy` comença amb `[7f/caf718]`.
Això correspon al següent camí de directori truncat: `work/7f/caf718`.

En aquell directori, trobaràs el fitxer `.command.run` que conté totes les comandes que Nextflow va executar en nom teu en el curs d'executar el pipeline.

??? abstract "Contingut del fitxer"

    ```console title="work/7f/caf71890cce1667c094d880f4b6dcc/.command.run"
    #!/bin/bash
    ### ---
    ### name: 'cowpy'
    ### container: 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
    ### outputs:
    ### - 'cowpy-COLLECTED-batch-output.txt'
    ### ...
    set -e
    set -u
    NXF_DEBUG=${NXF_DEBUG:=0}; [[ $NXF_DEBUG > 1 ]] && set -x
    NXF_ENTRY=${1:-nxf_main}


    nxf_sleep() {
      sleep $1 2>/dev/null || sleep 1;
    }

    nxf_date() {
        local ts=$(date +%s%3N);
        if [[ ${#ts} == 10 ]]; then echo ${ts}000
        elif [[ $ts == *%3N ]]; then echo ${ts/\%3N/000}
        elif [[ $ts == *3N ]]; then echo ${ts/3N/000}
        elif [[ ${#ts} == 13 ]]; then echo $ts
        else echo "Unexpected timestamp value: $ts"; exit 1
        fi
    }

    nxf_env() {
        echo '============= task environment ============='
        env | sort | sed "s/\(.*\)AWS\(.*\)=\(.\{6\}\).*/\1AWS\2=\3xxxxxxxxxxxxx/"
        echo '============= task output =================='
    }

    nxf_kill() {
        declare -a children
        while read P PP;do
            children[$PP]+=" $P"
        done < <(ps -e -o pid= -o ppid=)

        kill_all() {
            [[ $1 != $$ ]] && kill $1 2>/dev/null || true
            for i in ${children[$1]:=}; do kill_all $i; done
        }

        kill_all $1
    }

    nxf_mktemp() {
        local base=${1:-/tmp}
        mkdir -p "$base"
        if [[ $(uname) = Darwin ]]; then mktemp -d $base/nxf.XXXXXXXXXX
        else TMPDIR="$base" mktemp -d -t nxf.XXXXXXXXXX
        fi
    }

    nxf_fs_copy() {
      local source=$1
      local target=$2
      local basedir=$(dirname $1)
      mkdir -p $target/$basedir
      cp -fRL $source $target/$basedir
    }

    nxf_fs_move() {
      local source=$1
      local target=$2
      local basedir=$(dirname $1)
      mkdir -p $target/$basedir
      mv -f $source $target/$basedir
    }

    nxf_fs_rsync() {
      rsync -rRl $1 $2
    }

    nxf_fs_rclone() {
      rclone copyto $1 $2/$1
    }

    nxf_fs_fcp() {
      fcp $1 $2/$1
    }

    on_exit() {
        local last_err=$?
        local exit_status=${nxf_main_ret:=0}
        [[ ${exit_status} -eq 0 && ${nxf_unstage_ret:=0} -ne 0 ]] && exit_status=${nxf_unstage_ret:=0}
        [[ ${exit_status} -eq 0 && ${last_err} -ne 0 ]] && exit_status=${last_err}
        printf -- $exit_status > /workspaces/training/nextflow-run/work/7f/caf71890cce1667c094d880f4b6dcc/.exitcode
        set +u
        docker rm $NXF_BOXID &>/dev/null || true
        exit $exit_status
    }

    on_term() {
        set +e
        docker stop $NXF_BOXID
    }

    nxf_launch() {
        docker run -i --cpu-shares 1024 -e "NXF_TASK_WORKDIR" -v /workspaces/training/nextflow-run/work:/workspaces/training/nextflow-run/work -w "$NXF_TASK_WORKDIR" --name $NXF_BOXID community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273 /bin/bash -ue /workspaces/training/nextflow-run/work/7f/caf71890cce1667c094d880f4b6dcc/.command.sh
    }

    nxf_stage() {
        true
        # stage input files
        rm -f COLLECTED-batch-output.txt
        ln -s /workspaces/training/nextflow-run/work/7f/f435e3f2cf95979b5f3d7647ae6696/COLLECTED-batch-output.txt COLLECTED-batch-output.txt
    }

    nxf_unstage_outputs() {
        true
    }

    nxf_unstage_controls() {
        true
    }

    nxf_unstage() {
        if [[ ${nxf_main_ret:=0} == 0 ]]; then
            (set -e -o pipefail; (nxf_unstage_outputs | tee -a .command.out) 3>&1 1>&2 2>&3 | tee -a .command.err)
            nxf_unstage_ret=$?
        fi
        nxf_unstage_controls
    }

    nxf_main() {
        trap on_exit EXIT
        trap on_term TERM INT USR2
        trap '' USR1

        [[ "${NXF_CHDIR:-}" ]] && cd "$NXF_CHDIR"
        export NXF_BOXID="nxf-$(dd bs=18 count=1 if=/dev/urandom 2>/dev/null | base64 | tr +/ 0A | tr -d '\r\n')"
        NXF_SCRATCH=''
        [[ $NXF_DEBUG > 0 ]] && nxf_env
        touch /workspaces/training/nextflow-run/work/7f/caf71890cce1667c094d880f4b6dcc/.command.begin
        set +u
        set -u
        [[ $NXF_SCRATCH ]] && cd $NXF_SCRATCH
        export NXF_TASK_WORKDIR="$PWD"
        nxf_stage

        set +e
        (set -o pipefail; (nxf_launch | tee .command.out) 3>&1 1>&2 2>&3 | tee .command.err) &
        pid=$!
        wait $pid || nxf_main_ret=$?
        nxf_unstage
    }

    $NXF_ENTRY
    ```

Si cerques `nxf_launch` en aquest fitxer, hauries de veure alguna cosa així:

```console
nxf_launch() {
    docker run -i --cpu-shares 1024 -e "NXF_TASK_WORKDIR" -v /workspaces/training/nextflow-run/work:/workspaces/training/nextflow-run/work -w "$NXF_TASK_WORKDIR" --name $NXF_BOXID community.wave.seqera.io/library/pip_cowpy:131d6a1b707a8e65 /bin/bash -ue /workspaces/training/nextflow-run/work/7f/caf7189fca6c56ba627b75749edcb3/.command.sh
}
```

Aquesta comanda de llançament mostra que Nextflow està utilitzant una comanda `docker run` molt similar per llançar la crida de procés com vam fer quan la vam executar manualment.
També munta el subdirectori de treball corresponent al contenidor, estableix el directori de treball dins del contenidor en conseqüència, i executa el nostre script bash amb plantilla al fitxer `.command.sh`.

Això confirma que tot el treball dur que vam haver de fer manualment a la secció anterior ara el fa Nextflow per nosaltres!

### Conclusió

Entens quin paper juguen els contenidors en la gestió de versions d'eines de programari i assegurar la reproducibilitat.

Més generalment, tens una comprensió bàsica de quins són els components centrals dels pipelines Nextflow del món real i com estan organitzats.
Coneixes els fonaments de com Nextflow pot processar múltiples entrades de manera eficient, executar workflows compostos de múltiples passos connectats junts, aprofitar components de codi modulars, i utilitzar contenidors per a una major reproducibilitat i portabilitat.

### Què segueix?

Pren-te un altre descans! Això va ser una gran pila d'informació sobre com funcionen els pipelines Nextflow.

A l'última secció d'aquesta formació, aprofundirem més en el tema de la configuració.
Aprendràs com configurar l'execució del teu pipeline per adaptar-se a la teva infraestructura així com gestionar la configuració d'entrades i paràmetres.

---

## Qüestionari

<quiz>
Per què Nextflow crea un directori de tasca separat per a cada crida de procés?
- [ ] Per millorar la velocitat d'execució
- [ ] Per reduir l'ús de memòria
- [x] Per aïllar execucions i evitar col·lisions entre sortides
- [ ] Per habilitar la compressió de fitxers en paral·lel

Aprèn més: [1.3. Trobar les sortides originals i els registres](#13-find-the-original-outputs-and-logs)
</quiz>

<quiz>
Què fa l'opció `-ansi-log false` quan s'executa un workflow?
- [ ] Deshabilita tota la sortida de consola
- [x] Elimina el color de la sortida
- [x] Mostra tots els camins de directori de tasca en lloc de condensar-los en una línia
- [ ] Habilita el mode de depuració detallat

Aprèn més: [1.3.2. Fer que el terminal mostri més detalls](#132-make-the-terminal-show-more-details)

També pots utilitzar qualsevol de les següents variables d'entorn si prefereixes aquest estil:

```bash
export NXF_ANSI_LOG=0
# o
export NO_COLOR=1
```

</quiz>

<quiz>
Al codi `#!groovy channel.fromPath(params.input).splitCsv().map { line -> line[0] }`, què fa `#!groovy .map { line -> line[0] }`?
- [ ] Filtra les línies buides
- [ ] Ordena les línies alfabèticament
- [x] Extreu la primera columna de cada fila CSV
- [ ] Compta el nombre de línies

Aprèn més: [1.4.1. Carregar les dades d'entrada des del CSV](#141-loading-the-input-data-from-the-csv)
</quiz>

<quiz>
Per què és important incloure el valor d'entrada als noms dels fitxers de sortida (per exemple, `#!groovy "${greeting}-output.txt"`)?
- [ ] Per millorar la velocitat de processament
- [ ] Per habilitar la funcionalitat de resume
- [x] Per evitar que els fitxers de sortida se sobreescriguin entre si quan es processen múltiples entrades
- [ ] Per fer els fitxers més fàcils de comprimir

Aprèn més: [1.4.3. Com es nomenen les sortides](#143-how-the-outputs-are-named)
</quiz>

<quiz>
Quin és el propòsit de la declaració `include` en un workflow modularitzat?
- [ ] Copiar el codi del procés al fitxer de workflow
- [x] Importar una definició de procés d'un fitxer de mòdul extern
- [ ] Incloure configuracions
- [ ] Afegir comentaris de documentació

Aprèn més: [3. Executar pipelines modularitzats](#3-running-modularized-pipelines)
</quiz>

<quiz>
Quan modularitzes un workflow i l'executes amb `-resume`, què passa?
- [ ] La memòria cau es deshabilita per a processos modulars
- [ ] Totes les tasques s'han de tornar a executar
- [x] La memòria cau funciona normalment basant-se en els scripts de tasca generats
- [ ] Només el fitxer de workflow principal s'emmagatzema a la memòria cau

Aprèn més: [3.2. Executar el workflow](#32-run-the-workflow)
</quiz>

<quiz>
Què especifica la directiva `container` en una definició de procés?
- [ ] El directori de treball per al procés
- [ ] L'assignació màxima de memòria
- [x] L'URI de la imatge del contenidor a utilitzar per executar el procés
- [ ] El format del fitxer de sortida

Aprèn més: [4.2. Utilitzar un contenidor en un workflow](#42-use-a-container-in-a-workflow)
</quiz>

<quiz>
Al fitxer `.command.run`, què conté la funció `nxf_launch`?
- [ ] La informació de versió de Nextflow
- [ ] Els paràmetres del workflow
- [x] La comanda `docker run` amb muntatges de volum i configuracions del contenidor
- [ ] Les declaracions d'entrada del procés

Aprèn més: [4.2.4. Inspeccionar com Nextflow va llançar la tasca en contenidor](#424-inspect-how-nextflow-launched-the-containerized-task)
</quiz>

<quiz>
Què gestiona automàticament Nextflow quan executa un procés en contenidor? (Selecciona totes les que apliquin)
- [x] Descarregar la imatge del contenidor si és necessari
- [x] Muntar el directori de treball al contenidor
- [x] Executar l'script del procés dins del contenidor
- [x] Netejar la instància del contenidor després de l'execució

Aprèn més: [4. Utilitzar programari en contenidors](#4-using-containerized-software)
</quiz>
