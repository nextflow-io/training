# Part 1: Executar Nextflow

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducció assistida per IA - [més informació i suggeriments](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

En aquesta part, introduïm els conceptes bàsics per executar pipelines de Nextflow.
Comencem amb un workflow senzill de Hello World i progressem fins a un pipeline complet de múltiples passos que processa múltiples entrades en paral·lel utilitzant contenidors.

---

## 1. Hello World

El workflow `1-hello.nf` rep una salutació mitjançant un argument de línia de comandes i l'escriu en un fitxer.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello_world.svg"
</figure>

### 1.1. Llançar el workflow

Executeu la comanda següent al vostre terminal.

```bash
nextflow run 1-hello.nf --input 'Hello World!'
```

??? success "Sortida de la comanda"

    ```console hl_lines="6"
    N E X T F L O W   ~  version 26.04.4

    Launching `1-hello.nf` [infallible_volhard] revision: 82d40dbf94

    executor >  local (1)
    [6d/740edd] sayHello | 1 of 1 ✔

    Outputs:

      /workspaces/training/nextflow-run/results

      first_output: 1-hello/Hello World!-output.txt
    ```

La línia clau de la sortida és la línia d'estat del procés:

```console
[6d/740edd] sayHello | 1 of 1 ✔
```

Això ens indica que el procés `sayHello` s'ha executat correctament una vegada.
El prefix `[6d/740edd]` és una ruta truncada al directori de treball de la tasca — en parlarem més avall.
El bloc `Outputs:` que segueix llista tots els fitxers que el pipeline ha publicat, etiquetats d'acord amb el bloc `output` que es tracta a [1.4](#14-optional-code-walkthrough) més avall.

### 1.2. Trobar la sortida

Aquest workflow està configurat per publicar la seva sortida en un directori `results`.
Després d'executar-lo, hauríeu de trobar la sortida allà:

```console title="results/"
results
└── 1-hello
    └── Hello World!-output.txt
```

Obriu el fitxer per confirmar que conté `Hello World!`.

### 1.3. Explorar el directori `work/`

En segon pla, Nextflow crea un directori de tasca únic per a cada crida a un procés dins d'un directori anomenat `work/`.
El hash que es mostra a la sortida de la consola (`[6d/740edd]`) és la ruta a aquest directori.

```bash
ls work/6d/740edd*
```

A dins trobareu el fitxer de sortida juntament amb diversos fitxers de registre ocults:

- **`.command.sh`**: la comanda exacta que ha executat Nextflow
- **`.command.out`** / **`.command.err`**: stdout i stderr del procés
- **`.command.log`**: sortida de registre combinada
- **`.exitcode`**: el codi de sortida del procés

El fitxer `.command.sh` és especialment útil per depurar errors — mostra exactament el que s'ha executat.

### 1.4. Opcional: Revisió del codi

Entendre el codi no és essencial si simplement voleu executar pipelines, però si teniu curiositat, val la pena fer-hi un cop d'ull.

??? optional "Feu clic per explorar el codi associat a aquest exercici"

    Obrim `1-hello.nf` i mirem els seus components principals.

    ```groovy title="1-hello.nf" linenums="1"
    #!/usr/bin/env nextflow

    include { sayHello } from './modules/sayHello.nf'

    /*
     * Paràmetres del pipeline
     */
    params {
        input: String
    }

    workflow {

        main:
        // emet una salutació
        sayHello(params.input)

        publish:
        first_output = sayHello.out
    }

    output {
        first_output {
            path '1-hello'
            mode 'copy'
        }
    }
    ```

    Hi veiem el següent:

    - una instrucció `include` que apunta a un mòdul de `process`
    - un bloc `params` que defineix els paràmetres del pipeline
    - un bloc `workflow` que descriu el treball a realitzar
    - un bloc `output` que descriu què fer amb les sortides

    Mirem cadascun per separat.

    ### El mòdul `process`

    La instrucció `include` indica a Nextflow que carregui quelcom anomenat `sayHello` des d'un fitxer de codi separat.

    ```groovy title="1-hello.nf" linenums="3"
    include { sayHello } from './modules/sayHello.nf'
    ```

    En aquell fitxer, trobem la definició d'un procés anomenat `sayHello`:

    ```groovy title="modules/sayHello.nf" linenums="4"
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

    Un **process** defineix un únic pas del pipeline.
    Declara les seves entrades, sortides i l'script a executar.
    El qualificador `val` significa que l'entrada és un valor simple (string, número, etc.).
    El qualificador `path` significa que la sortida és una ruta de fitxer.

    És possible escriure la definició del procés al fitxer principal del workflow, però mantenir-los en fitxers de mòdul separats els fa reutilitzables: el mateix mòdul pot ser importat per múltiples scripts de workflow.

    ### El bloc `params`

    El bloc `params` declara els paràmetres de línia de comandes que accepta el workflow:

    ```groovy title="1-hello.nf" linenums="8"
    params {
        input: String
    }
    ```

    Qualsevol paràmetre declarat aquí queda disponible a la línia de comandes amb un doble guió (`--input`).
    Els tipus admesos inclouen `String`, `Integer`, `Float`, `Boolean` i `Path`.

    !!! tip "Consell"

        Els paràmetres del workflow sempre utilitzen dos guions (`--input`) per distingir-los dels indicadors propis de la CLI de Nextflow, que utilitzen un guió (p. ex. `-resume`).

    ### El bloc `workflow`

    El bloc **workflow** defineix la lògica del flux de dades: quins processos s'han d'executar i en quin ordre.

    ```groovy title="1-hello.nf" linenums="12"
    workflow {

        main:
        // emet una salutació
        sayHello(params.input)

        publish:
        first_output = sayHello.out
    }
    ```

    Aquí només s'està cridant un procés, de manera que és molt senzill; cobrirem exemples més realistes més endavant.

    La secció `main:` crida el procés `sayHello` amb el valor de `--input`.
    La secció `publish:` llista quines sortides s'han de copiar al directori de resultats.

    ### El bloc `output`

    El bloc `output` al final del fitxer especifica la ruta de destinació i el mode de còpia.

    ```groovy title="1-hello.nf" linenums="22"
    output {
        first_output {
            path '1-hello'
            mode 'copy'
        }
    }
    ```

    Cada entrada amb nom correspon a una etiqueta `publish:` del workflow i la mapeja a un subdirectori sota `results/`.

### Conclusió

Sabeu com executar un pipeline de Nextflow i trobar les seves sortides, i sabeu que el treball s'executa en directoris de tasca sota `work/`.

### Què segueix?

Descobriu com Nextflow gestiona múltiples entrades de manera eficient.

---

## 2. Processar múltiples entrades

Els pipelines del món real típicament processen moltes peces de dades, no només una.
El workflow `2-inputs.nf` llegeix d'un fitxer CSV i executa `sayHello` una vegada per fila, en paral·lel.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello-pipeline-multi-inputs-csv.svg"
</figure>

Primer executem el workflow i després mirem quin mecanisme utilitza Nextflow per gestionar aquestes múltiples entrades.

### 2.1. Executar el workflow

Executeu la comanda següent al vostre terminal.

```bash
nextflow run 2-inputs.nf --input data/greetings.csv
```

??? success "Sortida de la comanda"

    ```console hl_lines="6"
    N E X T F L O W   ~  version 26.04.4

    Launching `2-inputs.nf` [nauseous_babbage] revision: b90778224d

    executor >  local (3)
    [66/de7844] sayHello (3) | 3 of 3 ✔

    Outputs:

      /workspaces/training/nextflow-run/results

      first_output:
        - 2-inputs/Bonjour-output.txt
        - 2-inputs/Hello-output.txt
        - 2-inputs/Hola-output.txt
    ```

El `3 of 3` ens indica que el procés `sayHello` s'ha cridat tres vegades, una per cada fila del CSV.

Al directori `results`, ara hauríeu de veure tres fitxers de sortida, un per salutació:

```console title="results/2-inputs/"
2-inputs
├── Bonjour-output.txt
├── Hello-output.txt
└── Hola-output.txt
```

Obriu qualsevol dels fitxers de sortida per confirmar que cadascun conté una salutació.

La sortida condensada de dalt mostra una única línia de resum per a `sayHello`, però Nextflow en realitat ha llançat tres execucions de tasca separades, una per cada fila del CSV, i les ha executat en paral·lel tan aviat com la màquina ha tingut els recursos per fer-ho.

Igual que la tasca individual que heu explorat a [1.3](#13-explore-the-work-directory), cadascuna d'aquestes tres execucions obté el seu propi directori de tasca sota `work/`, completament aïllat de les altres:

```console title="work/"
work
├── 2d/276c63.../
│   ├── .command.sh
│   └── Hola-output.txt
├── ab/007682.../
│   ├── .command.sh
│   └── Bonjour-output.txt
└── d9/2476082.../
    ├── .command.sh
    └── Hello-output.txt
```

Cada `.command.sh` només conté la comanda per a aquella salutació concreta:

```console title="work/2d/276c63.../.command.sh"
#!/bin/bash -ue
echo 'Hola' > 'Hola-output.txt'
```

Aquest aïllament és el que fa que l'execució en paral·lel sigui segura: tres tasques que s'executen al mateix temps mai comparteixen un directori de treball, de manera que res del que escriu una tasca pot col·lisionar amb el que escriu una altra, fins i tot si produeixen fitxers amb el mateix nom.
També és per això que `-resume` (que es tracta a continuació) pot emmagatzemar en memòria cau i reutilitzar tasques individuals de manera independent: les entrades, sortides i registres de cada tasca viuen completament dins del seu propi directori, sense res compartit entre tasques que pugui quedar desincronitzat.

### 2.2. Executar el workflow de nou amb `-ansi-log false`

Per defecte, Nextflow condensa la sortida en una única línia de resum per procés.
Per veure cada crida a un procés llistada individualment, afegiu `-ansi-log false`:

```bash
nextflow run 2-inputs.nf --input data/greetings.csv -ansi-log false
```

??? success "Sortida de la comanda"

    ```console
    N E X T F L O W  ~  version 26.04.4
    Launching `2-inputs.nf` [extravagant_bardeen] - revision: b90778224d
    [43/0bac1c] Submitted process > sayHello (1)
    [2d/99f604] Submitted process > sayHello (2)
    [6d/7578d7] Submitted process > sayHello (3)

    Outputs:

      /workspaces/training/nextflow-run/results

      first_output:
        - 2-inputs/Hello-output.txt
        - 2-inputs/Bonjour-output.txt
        - 2-inputs/Hola-output.txt
    ```

Això mostra les tres crides al procés i el subdirectori de treball únic creat per a cadascuna.

### 2.3. Utilitzar `-resume` per ometre el treball completat

Ara canvieu al fitxer d'entrada ampliat, que afegeix dues salutacions més, i afegiu `-resume` a la línia de comandes:

```bash
nextflow run 2-inputs.nf --input data/greetings-extended.csv -resume
```

??? success "Sortida de la comanda"

    ```console hl_lines="6"
    N E X T F L O W   ~  version 26.04.4

    Launching `2-inputs.nf` [adoring_mayer] revision: b90778224d

    executor >  local (2)
    [84/2f3067] sayHello (5) | 5 of 5, cached: 3 ✔

    Outputs:

      /workspaces/training/nextflow-run/results

      first_output:
        - 2-inputs/Bonjour-output.txt
        - 2-inputs/Ciao-output.txt
        - 2-inputs/Hello-output.txt
        - 2-inputs/Hola-output.txt
        - 2-inputs/Ola-output.txt
    ```

Nextflow només ha executat les dues noves entrades.
Les tres salutacions processades en l'execució anterior s'han emmagatzemat en memòria cau i s'han reutilitzat automàticament.

Això també funciona per ometre l'execució de processos per a passos que ja s'han executat correctament en un pipeline de múltiples passos.
Per exemple, si una execució del pipeline ha estat interrompuda per un error del sistema, o si heu afegit nous passos a un pipeline en desenvolupament.

La capacitat de `-resume` és especialment valuosa en pipelines llargs on recuperar-se d'un error pot estalviar temps i recursos crítics.

### 2.4. Opcional: Revisió del codi

Entendre el codi no és essencial si simplement voleu executar pipelines, però si teniu curiositat, val la pena fer-hi un cop d'ull.

??? optional "Feu clic per explorar el codi associat a aquest exercici"

    El canvi clau a `2-inputs.nf` és a la secció `main:` del workflow:

    ```groovy title="2-inputs.nf" linenums="14" hl_lines="3 4 5"
        main:
        // crea un canal per a les entrades des d'un fitxer CSV
        greeting_ch = channel.fromPath(params.input)
                            .splitCsv()
                            .map { line -> line[0] }
        // emet una salutació
        sayHello(greeting_ch)
    ```

    El que veieu aquí s'anomena **canal**: una construcció de cua que gestiona les dades d'entrada d'una manera que facilita la paral·lelització d'operacions.

    - `channel.fromPath(params.input)` crea un canal a partir de la ruta del fitxer donada amb `--input`
    - `.splitCsv()` analitza el CSV en files
    - `#!groovy .map { line -> line[0] }` extreu la primera columna de cada fila

    El resultat és un canal que conté `Hello`, `Bonjour` i `Hola`.
    Quan es passa a `sayHello(greeting_ch)`, Nextflow crida automàticament el procés una vegada per element, executant-los en paral·lel quan els recursos ho permeten.

### Conclusió

Sabeu com processar múltiples entrades d'un fitxer CSV en paral·lel, i com utilitzar `-resume` per evitar repetir el treball completat.

### Què segueix?

Apreneu com un pipeline complet de múltiples passos encadena processos mitjançant canals, i com utilitzar contenidors per gestionar les eines d'anàlisi i les seves dependències.

---

## 3. Executar un pipeline de múltiples passos

Fins ara heu executat un únic procés i després l'heu executat múltiples vegades en paral·lel sobre un conjunt d'entrades.
Els pipelines reals solen anar més lluny: encadenen diversos processos, alimentant la sortida d'un cap al següent, i sovint depenen de més d'un programari al llarg del camí.
El workflow `main.nf` combina tots dos en un pipeline complet.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/hello_pipeline_complete.svg"
</figure>

Cada salutació d'entrada flueix a través dels quatre passos: `sayHello` l'escriu en un fitxer, `convertToUpper` converteix el text a majúscules, `collectGreetings` fusiona tots els resultats en un únic fitxer, i `cowpy` genera art ASCII a partir de la sortida fusionada utilitzant una eina en contenidor.
Nextflow connecta aquests passos amb canals: la sortida d'un procés es converteix en l'entrada del següent, de manera que tota la cadena s'executa automàticament a mesura que les dades estan disponibles, sense que hàgiu d'orquestrar cada pas manualment.

Tingueu en compte que aquest workflow utilitza mòduls: cada procés es defineix en el seu propi fitxer sota `modules/`, i `main.nf` els importa amb instruccions `include` en lloc de definir-los en línia.
Això fa que cada procés sigui reutilitzable en múltiples workflows sense duplicar codi. Per saber-ne més, consulteu la secció d'exploració del codi més avall.

### 3.1. Executar el workflow

Executeu la comanda següent al vostre terminal.

```bash
nextflow run main.nf --input data/greetings.csv
```

El paràmetre `character` té per defecte `turkey` a `nextflow.config`, de manera que l'art ASCII utilitza un gall dindi tret que el substituïu (proveu d'afegir `--character tux`).

??? success "Sortida de la comanda"

    ```console hl_lines="6 7 8 9"
    N E X T F L O W   ~  version 26.04.4

    Launching `main.nf` [nostalgic_brahmagupta] revision: ce74f81996

    executor >  local (8)
    [56/8499f6] sayHello (3)       | 3 of 3 ✔
    [cc/0ee42a] convertToUpper (3) | 3 of 3 ✔
    [eb/0f2e24] collectGreetings   | 1 of 1 ✔
    [b5/34e07f] cowpy              | 1 of 1 ✔

    Outputs:

      /workspaces/training/nextflow-run/results

      first_output:
        - batch/intermediates/Bonjour-output.txt
        - batch/intermediates/Hello-output.txt
        - batch/intermediates/Hola-output.txt

      uppercased:
        - batch/intermediates/UPPER-Bonjour-output.txt
        - batch/intermediates/UPPER-Hello-output.txt
        - batch/intermediates/UPPER-Hola-output.txt

      collected: batch/intermediates/COLLECTED-batch-output.txt

      batch_report: batch/batch-report.txt

      cowpy_art: batch/cowpy-COLLECTED-batch-output.txt
    ```

S'han executat quatre processos, però no el mateix nombre de vegades.
`sayHello` i `convertToUpper` s'han executat una vegada per entrada (3 of 3): cada salutació s'ha d'escriure i convertir a majúscules per separat.
`collectGreetings` i `cowpy` s'han executat només una vegada (1 of 1): fusionar les salutacions i generar l'art ASCII només té sentit un cop tots els resultats individuals estan disponibles.
Aquesta forma de dispersió i convergència, diverses tasques en paral·lel que alimenten un nombre menor de tasques posteriors, és habitual en els pipelines reals.

Nextflow no espera que un pas sencer acabi abans de començar el següent.
Tan aviat com una sortida de `sayHello` està llesta, la tasca `convertToUpper` corresponent pot començar, de manera que les tasques de diferents processos s'executen de manera concurrent en lloc de fer-ho en lots estrictes.
`collectGreetings` i `cowpy` sí que han d'esperar, ja que cadascun depèn que tots els resultats anteriors estiguin disponibles primer.

El directori `results` reflecteix aquesta convergència, a més del que l'autor del pipeline ha triat publicar i on: recordeu el bloc `output` de la revisió del codi a 1.4, que és el que defineix aquesta estructura.

```console title="results/"
results
└── batch
    ├── batch-report.txt
    ├── cowpy-COLLECTED-batch-output.txt
    └── intermediates
        ├── Bonjour-output.txt
        ├── COLLECTED-batch-output.txt
        ├── Hello-output.txt
        ├── Hola-output.txt
        ├── UPPER-Bonjour-output.txt
        ├── UPPER-Hello-output.txt
        └── UPPER-Hola-output.txt
```

El directori de nivell superior rep el nom del paràmetre `batch`, que per defecte és `batch`; el veureu canviar en exercicis posteriors.

Comproveu `cowpy-COLLECTED-batch-output.txt` per veure el fitxer d'art ASCII.

??? abstract "Contingut del fitxer"

    ```console title="results/batch/cowpy-COLLECTED-batch-output.txt"
     _________
    / HELLO   \
    | BONJOUR |
    \ HOLA    /
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

Igual que a [2.1](#21-run-the-workflow), cadascuna d'aquestes 8 execucions de tasca, de tots quatre processos, obté el seu propi directori sota `work/`, completament aïllat de les altres.
`collectGreetings` és un bon exemple de per què això importa: depèn de les sortides de les tres tasques `convertToUpper`, que viuen en tres directoris de tasca diferents, de manera que Nextflow crea enllaços simbòlics a aquests fitxers dins del propi directori de `collectGreetings` en lloc de llegir directament des dels directoris de les tasques anteriors:

```console title="work/eb/0f2e24.../"
COLLECTED-batch-output.txt
UPPER-Bonjour-output.txt -> ../../69/e6c057.../UPPER-Bonjour-output.txt
UPPER-Hello-output.txt   -> ../../cc/ba7d19.../UPPER-Hello-output.txt
UPPER-Hola-output.txt    -> ../../cc/0ee42a.../UPPER-Hola-output.txt
batch-report.txt
.command.sh
```

Cada tasca només veu els fitxers específics que necessita, sigui d'on sigui que provinguin, i mai el contingut intern del directori d'una altra tasca.
Al llarg de tot un pipeline, el mateix aïllament que heu vist amb un únic procés a [2.1](#21-run-the-workflow) és el que permet a Nextflow executar totes les tasques de tots els processos de manera concurrent i segura.

!!! note "Nota"

    El pas `cowpy` s'executa dins d'un contenidor Docker en lloc de dependre del programari instal·lat localment.
    Un contenidor empaqueta una aplicació juntament amb tot el que necessita per executar-se, de manera que no heu d'instal·lar ni gestionar les dependències vosaltres mateixos, i el pipeline es comporta de la mateixa manera en qualsevol màquina que pugui executar el contenidor.
    Nextflow també admet Conda com a alternativa als contenidors; consulteu la [Part 2](./02_configure_pipeline.md) per saber com canviar entre ells.

### 3.2. Opcional: Revisió del codi

Entendre el codi no és essencial si simplement voleu executar pipelines, però si teniu curiositat, val la pena fer-hi un cop d'ull.

??? optional "Feu clic per explorar el codi associat a aquest exercici"

    ### Com flueixen les dades d'un pas al següent

    Cada procés passa el seu canal de sortida al següent:

    ```groovy title="main.nf" linenums="19" hl_lines="7 8 9"
        main:
        // crea un canal per a les entrades des d'un fitxer CSV
        greeting_ch = channel.fromPath(params.input)
                            .splitCsv()
                            .map { line -> line[0] }
        sayHello(greeting_ch)
        convertToUpper(sayHello.out)
        collectGreetings(convertToUpper.out.collect(), params.batch)
        cowpy(collectGreetings.out.outfile, params.character)
    ```

    El patró `processName.out` fa referència al canal de sortida d'un procés.

    L'operador `.collect()` recull totes les sortides individuals de `convertToUpper` en un únic element de canal abans de passar-les a `collectGreetings`.

    ### Ús de mòduls de processos

    `main.nf` no defineix cap codi de procés directament.
    En canvi, importa cada procés des del seu propi fitxer sota `modules/`:

    ```groovy title="main.nf" linenums="3"
    include { sayHello } from './modules/sayHello.nf'
    include { convertToUpper } from './modules/convertToUpper.nf'
    include { collectGreetings } from './modules/collectGreetings.nf'
    include { cowpy } from './modules/cowpy.nf'
    ```

    Cada fitxer de mòdul conté una única definició de procés, estructurada de la mateixa manera que el mòdul `sayHello` a [1.4](#14-optional-code-walkthrough).
    Mantenir els processos en fitxers separats els fa reutilitzables en múltiples workflows sense duplicar codi.

    <figure class="excalidraw">
    --8<-- "docs/en/docs/hello_nextflow/img/modules.svg"
    </figure>

    ### Ús de programari en contenidor

    El procés `cowpy` s'executa dins d'un contenidor Docker especificat al seu fitxer de mòdul:

    ```groovy title="modules/cowpy.nf" linenums="2" hl_lines="3"
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

    Nextflow descarrega automàticament la imatge, executa l'script dins del contenidor i fa la neteja posterior.
    Docker està habilitat per a aquest projecte a `nextflow.config`:

    ```groovy title="nextflow.config"
    docker.enabled = true
    ```

    Aquesta única línia habilita Docker per a qualsevol procés del pipeline que tingui un contenidor especificat.

### Conclusió

Heu executat un pipeline complet de múltiples passos que processa múltiples entrades en paral·lel utilitzant una eina en contenidor.

### Què segueix?

Aneu a la [Part 2](./02_configure_pipeline.md), on aprendreu a configurar el comportament del pipeline utilitzant `nextflow.config`.

---

## Resum

En aquesta part heu après a:

- Executar un workflow de Nextflow i trobar les seves sortides
- Explorar el directori `work/` i els seus fitxers de registre
- Processar múltiples entrades d'un fitxer CSV en paral·lel
- Utilitzar `-resume` per ometre el treball completat en afegir noves entrades
- Executar un pipeline de múltiples passos que utilitza una eina en contenidor
