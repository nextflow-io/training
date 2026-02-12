# Part 3: Configuració d'execució

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducció assistida per IA - [més informació i suggeriments](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Aquesta secció explorarà com gestionar la configuració d'un pipeline de Nextflow per personalitzar el seu comportament, adaptar-lo a diferents entorns i optimitzar l'ús de recursos _sense alterar ni una sola línia del codi del workflow_.

Hi ha múltiples maneres de fer-ho, que es poden utilitzar en combinació i s'interpreten segons l'ordre de precedència descrit a la documentació de [Configuration](https://nextflow.io/docs/latest/config.html).

En aquesta part del curs, us mostrarem el mecanisme de fitxer de configuració més simple i comú, el fitxer `nextflow.config`, que ja vau trobar a la secció sobre contenidors de la Part 2.

Repassarem components essencials de la configuració de Nextflow com ara directives de procés, executors, perfils i fitxers de paràmetres.
Aprenent a utilitzar aquestes opcions de configuració de manera efectiva, podreu aprofitar al màxim la flexibilitat, escalabilitat i rendiment dels pipelines de Nextflow.

Per exercitar aquests elements de configuració, executarem una còpia nova del workflow que vam executar per última vegada al final de la Part 2 d'aquest curs de formació, reanomenat `3-main.nf`.

Si no esteu familiaritzats amb el pipeline Hello o us aniria bé un recordatori, consulteu [aquesta pàgina d'informació](../info/hello_pipeline.md).

---

## 1. Gestionar paràmetres d'entrada del workflow

??? example "Escenari"

    Heu descarregat un pipeline i voleu executar-lo repetidament amb els mateixos fitxers d'entrada i configuració, però no voleu escriure tots els paràmetres cada vegada.
    O potser esteu configurant el pipeline per a un col·lega que no se sent còmode amb arguments de línia de comandes.

Començarem amb un aspecte de la configuració que és simplement una extensió del que hem estat treballant fins ara: la gestió de paràmetres d'entrada.

Actualment, el nostre workflow està configurat per acceptar diversos valors de paràmetres via línia de comandes, declarats en un bloc `params` al propi script del workflow.
Un té un valor per defecte establert com a part de la seva declaració.

No obstant això, potser voleu establir valors per defecte per a tots ells, o sobreescriure el valor per defecte existent sense haver d'especificar paràmetres a la línia de comandes o modificar el fitxer de script original.

Hi ha múltiples maneres de fer-ho; us mostrarem tres maneres bàsiques que s'utilitzen molt comunament.

### 1.1. Configurar valors a `nextflow.config`

Aquest és l'enfocament més simple, tot i que possiblement és el menys flexible ja que el fitxer principal `nextflow.config` no és quelcom que vulgueu estar editant per a cada execució.
Però té l'avantatge de separar les preocupacions de _declarar_ els paràmetres al workflow (que definitivament hi pertany) versus proporcionar _valors per defecte_, que estan més a casa en un fitxer de configuració.

Fem-ho en dos passos.

#### 1.1.1. Crear un bloc `params` al fitxer de configuració

Feu els següents canvis de codi al fitxer `nextflow.config`:

=== "Després"

    ```groovy title="nextflow.config" linenums="1" hl_lines="3-10"
    docker.enabled = true

    /*
    * Paràmetres del pipeline
    */
    params {
        input = 'data/greetings.csv'
        batch = 'batch'
        character = 'turkey'
    }
    ```

=== "Abans"

    ```groovy title="nextflow.config" linenums="1"
    docker.enabled = true
    ```

Tingueu en compte que no hem copiat simplement el bloc `params` del workflow al fitxer de configuració.
Per al paràmetre `batch` que ja tenia un valor per defecte declarat, la sintaxi és una mica diferent.
Al fitxer del workflow, això és una declaració tipada.
A la configuració, són assignacions de valors.

Tècnicament, això és suficient per sobreescriure els valors per defecte encara especificats al fitxer del workflow.
Podríeu modificar el valor per defecte per a `batch` i executar el workflow per assegurar-vos que el valor establert al fitxer de configuració sobreescriu el del fitxer del workflow.

Però amb l'esperit de moure la configuració completament al fitxer de configuració, eliminem aquest valor per defecte del fitxer del workflow completament.

#### 1.1.2. Eliminar el valor per defecte per a `batch` al fitxer del workflow

Feu el següent canvi de codi al fitxer del workflow `3-main.nf`:

=== "Després"

    ```groovy title="3-main.nf" linenums="9" hl_lines="6"
    /*
    * Paràmetres del pipeline
    */
    params {
        input: Path
        batch: String
        character: String
    }
    ```

=== "Abans"

    ```groovy title="3-main.nf" linenums="9" hl_lines="6"
    /*
    * Paràmetres del pipeline
    */
    params {
        input: Path
        batch: String = 'batch'
        character: String
    }
    ```

Ara el fitxer del workflow en si no estableix cap valor per defecte per a aquests paràmetres.

#### 1.1.3. Executar el pipeline

Provem que funciona correctament sense especificar cap paràmetre a la línia de comandes.

```bash
nextflow run 3-main.nf
```

??? success "Sortida de la comanda"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `3-main.nf` [disturbed_einstein] DSL2 - revision: ede9037d02

    executor >  local (8)
    [f0/35723c] sayHello (2)       | 3 of 3 ✔
    [40/3efd1a] convertToUpper (3) | 3 of 3 ✔
    [17/e97d32] collectGreetings   | 1 of 1 ✔
    [98/c6b57b] cowpy              | 1 of 1 ✔
    ```

Això encara produeix la mateixa sortida que anteriorment.

La sortida final d'art ASCII es troba al directori `results/3-main/`, amb el nom `cowpy-COLLECTED-batch-output.txt`, igual que abans.

??? abstract "Contingut del fitxer"

    ```console title="results/3-main/cowpy-COLLECTED-batch-output.txt"
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

Funcionalment, aquest moviment no ha canviat res, però conceptualment és una mica més net tenir els valors per defecte establerts al fitxer de configuració.

### 1.2. Utilitzar un fitxer de configuració específic per a l'execució

??? example "Escenari"

    Voleu experimentar amb diferents configuracions sense modificar el vostre fitxer de configuració principal.

Podeu fer-ho creant un nou fitxer `nextflow.config` en un subdirectori que utilitzareu com a directori de treball per als vostres experiments.

#### 1.2.1. Crear el directori de treball amb una configuració en blanc

Comencem creant un nou directori i movent-nos-hi:

```bash
mkdir -p tux-run
cd tux-run
```

Després, creeu un fitxer de configuració en blanc en aquest directori:

```bash
touch nextflow.config
```

Això produeix un fitxer buit.

#### 1.2.2. Configurar la configuració experimental

Ara obriu el nou fitxer i afegiu els paràmetres que voleu personalitzar:

```groovy title="tux-run/nextflow.config" linenums="1"
params {
    input = '../data/greetings.csv'
    batch = 'experiment'
    character = 'tux'
}
```

Tingueu en compte que el camí al fitxer d'entrada ha de reflectir l'estructura de directoris.

#### 1.2.3. Executar el pipeline

Ara podem executar el nostre pipeline des del nostre nou directori de treball.
Assegureu-vos d'adaptar el camí en conseqüència!

```bash
nextflow run ../3-main.nf
```

??? success "Sortida de la comanda"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `../3-main.nf` [trusting_escher] DSL2 - revision: 356df0818d

    executor >  local (8)
    [59/b66913] sayHello (2)       [100%] 3 of 3 ✔
    [ad/f06364] convertToUpper (3) [100%] 3 of 3 ✔
    [10/714895] collectGreetings   [100%] 1 of 1 ✔
    [88/3ece98] cowpy              [100%] 1 of 1 ✔
    ```

Això crearà un nou conjunt de directoris sota `tux-run/` incloent `tux-run/work/` i `tux-run/results/`.

En aquesta execució, Nextflow combina el `nextflow.config` del nostre directori actual amb el `nextflow.config` del directori arrel del pipeline, i així sobreescriu el caràcter per defecte (turkey) amb el caràcter tux.

El fitxer de sortida final hauria de contenir el caràcter tux dient les salutacions.

??? abstract "Contingut del fitxer"

    ```console title="tux-run/results/3-main/cowpy-COLLECTED-experiment-output.txt"
    _________
    / HELLO   \
    | BONJOUR |
    \ HOLà    /
    ---------
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

Això és tot; ara teniu un espai per experimentar sense modificar la vostra configuració 'normal'.

!!! warning "Advertència"

    Assegureu-vos de tornar al directori anterior abans de passar a la següent secció!

    ```bash
    cd ..
    ```

Ara vegem una altra manera útil d'establir valors de paràmetres.

### 1.3. Utilitzar un fitxer de paràmetres

??? example "Escenari"

    Necessiteu compartir paràmetres d'execució exactes amb un col·laborador, o registrar-los per a una publicació.

L'enfocament del subdirectori funciona molt bé per experimentar, però implica una mica de configuració i requereix que adapteu els camins en conseqüència.
Hi ha un enfocament més simple per quan voleu executar el vostre pipeline amb un conjunt específic de valors, o permetre que algú altre ho faci amb un esforç mínim.

Nextflow ens permet especificar paràmetres via un [fitxer de paràmetres](https://nextflow.io/docs/latest/config.html#parameter-file) en format YAML o JSON, cosa que fa molt convenient gestionar i distribuir conjunts alternatius de valors per defecte, per exemple, així com valors de paràmetres específics d'execució.

#### 1.3.1. Examinar el fitxer de paràmetres d'exemple

Per demostrar això, proporcionem un fitxer de paràmetres d'exemple al directori actual, anomenat `test-params.yaml`:

```yaml title="test-params.yaml" linenums="1"
input: "data/greetings.csv"
batch: "yaml"
character: "stegosaurus"
```

Aquest fitxer de paràmetres conté un parell clau-valor per a cadascuna de les entrades que volem especificar.
Tingueu en compte l'ús de dos punts (`:`) en lloc de signes d'igual (`=`) si compareu la sintaxi amb el fitxer de configuració.
El fitxer de configuració està escrit en Groovy, mentre que el fitxer de paràmetres està escrit en YAML.

!!! info "Info"

    També proporcionem una versió JSON del fitxer de paràmetres com a exemple però no l'executarem aquí.
    Sentiu-vos lliures de provar-la pel vostre compte.

#### 1.3.2. Executar el pipeline

Per executar el workflow amb aquest fitxer de paràmetres, simplement afegiu `-params-file <nom_fitxer>` a la comanda base.

```bash
nextflow run 3-main.nf -params-file test-params.yaml
```

??? success "Sortida de la comanda"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `3-main.nf` [disturbed_sammet] DSL2 - revision: ede9037d02

    executor >  local (8)
    [2b/9a7d1e] sayHello (2)       | 3 of 3 ✔
    [5c/8f3b2a] convertToUpper (3) | 3 of 3 ✔
    [a3/29d8fb] collectGreetings   | 1 of 1 ✔
    [b7/83ef12] cowpy              | 1 of 1 ✔
    ```

El fitxer de sortida final hauria de contenir el caràcter stegosaurus dient les salutacions.

??? abstract "Contingut del fitxer"

    ```console title="results/3-main/cowpy-COLLECTED-yaml-output.txt"
    _________
    / HELLO   \
    | HOLà    |
    \ BONJOUR /
    ---------
    \                             .       .
    \                           / `.   .' "
      \                  .---.  <    > <    >  .---.
      \                 |    \  \ - ~ ~ - /  /    |
            _____          ..-~             ~-..-~
            |     |   \~~~\.'                    `./~~~/
          ---------   \__/                        \__/
          .'  O    \     /               /       \  "
        (_____,    `._.'               |         }  \/~~~/
          `----.          /       }     |        /    \__/
                `-.      |       /      |       /      `. ,~~|
                    ~-.__|      /_ - ~ ^|      /- _      `..-'
                        |     /        |     /     ~-.     `-. _  _  _
                        |_____|        |_____|         ~ - . _ _ _ _ _>
    ```

Utilitzar un fitxer de paràmetres pot semblar excessiu quan només teniu uns quants paràmetres per especificar, però alguns pipelines esperen desenes de paràmetres.
En aquests casos, utilitzar un fitxer de paràmetres ens permetrà proporcionar valors de paràmetres en temps d'execució sense haver d'escriure línies de comandes massives i sense modificar el script del workflow.

També facilita la distribució de conjunts de paràmetres a col·laboradors, o com a informació de suport per a una publicació, per exemple.
Això fa que el vostre treball sigui més reproduïble per altres.

### Conclusió

Sabeu com aprofitar les opcions de configuració clau per gestionar les entrades del workflow.

### Què segueix?

Apreneu com gestionar on i com es publiquen les sortides del vostre workflow.

---

## 2. Gestionar sortides del workflow

??? example "Escenari"

    El vostre pipeline publica sortides a un directori codificat de manera fixa, però voleu organitzar els resultats per projecte o nom d'experiment sense editar el codi del workflow cada vegada.

El workflow que hem heretat utilitza camins per a declaracions de sortida a nivell de workflow, cosa que no és terriblement flexible i implica molta repetició.

Vegem algunes maneres comunes en què podríeu configurar això per ser més flexible.

### 2.1. Personalitzar el nom del directori `outputDir`

Cada versió del workflow que hem executat fins ara ha publicat les seves sortides a un subdirectori diferent codificat a les definicions de sortida.

Vam canviar on estava aquest subdirectori a la Part 1 utilitzant la bandera CLI `-output-dir`, però això encara és només una cadena estàtica.
En canvi, configurem això en un fitxer de configuració, on podem definir camins dinàmics més complexos.
Podríem crear un paràmetre completament nou per a això, però utilitzem el paràmetre `batch` ja que està just allà.

#### 2.1.1. Establir un valor per a `outputDir` al fitxer de configuració

El camí que Nextflow utilitza per publicar sortides està controlat per l'opció `outputDir`.
Per canviar el camí per a totes les sortides, podeu establir un valor per a aquesta opció al fitxer de configuració `nextflow.config`.

Afegiu el següent codi al fitxer `nextflow.config`:

=== "Després"

    ```groovy title="nextflow.config" linenums="9" hl_lines="10-13"
    /*
    * Paràmetres del pipeline
    */
    params {
        input = 'data/greetings.csv'
        batch = 'batch'
        character = 'turkey'
    }

    /*
    * Configuració de sortida
    */
    outputDir = "results_config/${params.batch}"
    ```

=== "Abans"

    ```groovy title="nextflow.config" linenums="9"
    /*
    * Paràmetres del pipeline
    */
    params {
        input = 'data/greetings.csv'
        batch = 'batch'
        character = 'turkey'
    }
    ```

Això reemplaçarà el camí per defecte integrat, `results/`, amb `results_config/` més el valor del paràmetre `batch` com a subdirectori.

Recordeu que també podeu establir aquesta opció des de la línia de comandes utilitzant el paràmetre `-output-dir` a la vostra comanda (`-o` per abreujar), però llavors no podríeu utilitzar el valor del paràmetre `batch`.
Utilitzar la bandera CLI sobreescriurà `outputDir` a la configuració si està establert.

#### 2.1.2. Eliminar la part repetida del camí codificat

Encara tenim un subdirectori codificat a les opcions de sortida, així que eliminem-lo ara.

Feu els següents canvis de codi al fitxer del workflow:

=== "Després"

    ```groovy title="3-main.nf" linenums="42" hl_lines="3 7 11 15 19"
    output {
        first_output {
            path 'intermediates'
            mode 'copy'
        }
        uppercased {
            path 'intermediates'
            mode 'copy'
        }
        collected {
            path 'intermediates'
            mode 'copy'
        }
        batch_report {
            path ''
            mode 'copy'
        }
        cowpy_art {
            path ''
            mode 'copy'
        }
    }
    ```

=== "Abans"

    ```groovy title="3-main.nf" linenums="42" hl_lines="3 7 11 15 19"
    output {
        first_output {
            path '3-main/intermediates'
            mode 'copy'
        }
        uppercased {
            path '3-main/intermediates'
            mode 'copy'
        }
        collected {
            path '3-main/intermediates'
            mode 'copy'
        }
        batch_report {
            path '3-main'
            mode 'copy'
        }
        cowpy_art {
            path '3-main'
            mode 'copy'
        }
    }
    ```

També podríem haver afegit simplement `${params.batch}` a cada camí en lloc de modificar el valor per defecte d'`outputDir`, però això és més concís.

#### 2.1.3. Executar el pipeline

Provem que funciona correctament, establint el nom del batch a `outdir` des de la línia de comandes.

```bash
nextflow run 3-main.nf --batch outdir
```

??? success "Sortida de la comanda"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `3-main.nf` [amazing_church] DSL2 - revision: 6e18cd130e

    executor >  local (8)
    [9c/6a03ea] sayHello (2)       [100%] 3 of 3 ✔
    [11/9e58a6] convertToUpper (3) [100%] 3 of 3 ✔
    [c8/1977e5] collectGreetings   [100%] 1 of 1 ✔
    [38/f01eda] cowpy              [100%] 1 of 1 ✔
    ```

Això encara produeix la mateixa sortida que anteriorment, excepte que aquesta vegada trobem les nostres sortides sota `results_config/outdir/`.

??? abstract "Contingut del directori"

    ```console
    results_config/outdir
    ├── cowpy-COLLECTED-outdir-output.txt
    ├── intermediates
    │   ├── Bonjour-output.txt
    │   ├── COLLECTED-outdir-output.txt
    │   ├── Hello-output.txt
    │   ├── Holà-output.txt
    │   ├── UPPER-Bonjour-output.txt
    │   ├── UPPER-Hello-output.txt
    │   └── UPPER-Holà-output.txt
    └── outdir-report.txt
    ```

Podeu combinar aquest enfocament amb definicions de camí personalitzades per construir qualsevol jerarquia de directoris que vulgueu.

### 2.2. Organitzar sortides per procés

Una manera popular d'organitzar sortides encara més és fer-ho per procés, _és a dir_, crear subdirectoris per a cada procés executat al pipeline.

#### 2.2.1. Reemplaçar els camins de sortida per una referència als noms de procés

Tot el que necessiteu fer és referenciar el nom del procés com a `<procés>.name` a la declaració del camí de sortida.

Feu els següents canvis al fitxer del workflow:

=== "Després"

    ```groovy title="3-main.nf" linenums="42" hl_lines="3 7 11 15 19"
    output {
        first_output {
            path { sayHello.name }
            mode 'copy'
        }
        uppercased {
            path { convertToUpper.name }
            mode 'copy'
        }
        collected {
            path { collectGreetings.name }
            mode 'copy'
        }
        batch_report {
            path { collectGreetings.name }
            mode 'copy'
        }
        cowpy_art {
            path { cowpy.name }
            mode 'copy'
        }
    }
    ```

=== "Abans"

    ```groovy title="3-main.nf" linenums="42" hl_lines="3 7 11 15 19"
    output {
        first_output {
            path 'intermediates'
            mode 'copy'
        }
        uppercased {
            path 'intermediates'
            mode 'copy'
        }
        collected {
            path 'intermediates'
            mode 'copy'
        }
        batch_report {
            path ''
            mode 'copy'
        }
        cowpy_art {
            path ''
            mode 'copy'
        }
    }
    ```

Això elimina els elements codificats restants de la configuració del camí de sortida.

#### 2.2.2. Executar el pipeline

Provem que funciona correctament, establint el nom del batch a `pnames` des de la línia de comandes.

```bash
nextflow run 3-main.nf --batch pnames
```

??? success "Sortida de la comanda"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `3-main.nf` [jovial_mcclintock] DSL2 - revision: ede9037d02

    executor >  local (8)
    [4a/c2e6b8] sayHello (2)       | 3 of 3 ✔
    [6f/d4a172] convertToUpper (3) | 3 of 3 ✔
    [e8/4f19d7] collectGreetings   | 1 of 1 ✔
    [f2/a85c36] cowpy              | 1 of 1 ✔
    ```

Això encara produeix la mateixa sortida que anteriorment, excepte que aquesta vegada trobem les nostres sortides sota `results_config/pnames/`, i estan agrupades per procés.

??? abstract "Contingut del directori"

    ```console
    results_config/pnames/
    ├── collectGreetings
    │   ├── COLLECTED-pnames-output.txt
    │   └── pnames-report.txt
    ├── convertToUpper
    │   ├── UPPER-Bonjour-output.txt
    │   ├── UPPER-Hello-output.txt
    │   └── UPPER-Holà-output.txt
    ├── cowpy
    │   └── cowpy-COLLECTED-pnames-output.txt
    └── sayHello
        ├── Bonjour-output.txt
        ├── Hello-output.txt
        └── Holà-output.txt
    ```

!!! note "Nota"

    Tingueu en compte que aquí hem esborrat la distinció entre `intermediates` versus sortides finals que estan al nivell superior.
    Podeu barrejar i combinar aquests enfocaments i fins i tot incloure múltiples variables, per exemple establint el camí de la primera sortida com a `#!groovy "${params.batch}/intermediates/${sayHello.name}"`

### 2.3. Establir el mode de publicació a nivell de workflow

Finalment, amb l'esperit de reduir la quantitat de codi repetitiu, podem reemplaçar les declaracions de `mode` per sortida amb una sola línia a la configuració.

#### 2.3.1. Afegir `workflow.output.mode` al fitxer de configuració

Afegiu el següent codi al fitxer `nextflow.config`:

=== "Després"

    ```groovy title="nextflow.config" linenums="2" hl_lines="5"
    /*
    * Configuració de sortida
    */
    outputDir = "results_config/${params.batch}"
    workflow.output.mode = 'copy'
    ```

=== "Abans"

    ```groovy title="nextflow.config" linenums="12"
    /*
    * Configuració de sortida
    */
    outputDir = "results_config/${params.batch}"
    ```

Igual que l'opció `outputDir`, donar a `workflow.output.mode` un valor al fitxer de configuració seria suficient per sobreescriure el que està establert al fitxer del workflow, però eliminem el codi innecessari de totes maneres.

#### 2.3.2. Eliminar el mode de sortida del fitxer del workflow

Feu els següents canvis al fitxer del workflow:

=== "Després"

    ```groovy title="3-main.nf" linenums="42"
    output {
        first_output {
            path { sayHello.name }
        }
        uppercased {
            path { convertToUpper.name }
        }
        collected {
            path { collectGreetings.name }
        }
        batch_report {
            path { collectGreetings.name }
        }
        cowpy_art {
            path { cowpy.name }
        }
    }
    ```

=== "Abans"

    ```groovy title="3-main.nf" linenums="42" hl_lines="3 7 11 15 19"
    output {
        first_output {
            path { sayHello.name }
            mode 'copy'
        }
        uppercased {
            path { convertToUpper.name }
            mode 'copy'
        }
        collected {
            path { collectGreetings.name }
            mode 'copy'
        }
        batch_report {
            path { collectGreetings.name }
            mode 'copy'
        }
        cowpy_art {
            path { cowpy.name }
            mode 'copy'
        }
    }
    ```

Això és més concís, oi?

#### 2.3.3. Executar el pipeline

Provem que funciona correctament, establint el nom del batch a `outmode` des de la línia de comandes.

```bash
nextflow run 3-main.nf --batch outmode
```

??? success "Sortida de la comanda"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `3-main.nf` [rowdy_sagan] DSL2 - revision: ede9037d02

    executor >  local (8)
    [5b/d91e3c] sayHello (2)       | 3 of 3 ✔
    [8a/f6c241] convertToUpper (3) | 3 of 3 ✔
    [89/cd3a48] collectGreetings   | 1 of 1 ✔
    [9e/71fb52] cowpy              | 1 of 1 ✔
    ```

Això encara produeix la mateixa sortida que anteriorment, excepte que aquesta vegada trobem les nostres sortides sota `results_config/outmode/`.
Encara són totes còpies adequades, no enllaços simbòlics.

??? abstract "Contingut del directori"

    ```console
    results_config/outmode/
    ├── collectGreetings
    │   ├── COLLECTED-outmode-output.txt
    │   └── outmode-report.txt
    ├── convertToUpper
    │   ├── UPPER-Bonjour-output.txt
    │   ├── UPPER-Hello-output.txt
    │   └── UPPER-Holà-output.txt
    ├── cowpy
    │   └── cowpy-COLLECTED-outmode-output.txt
    └── sayHello
        ├── Bonjour-output.txt
        ├── Hello-output.txt
        └── Holà-output.txt
    ```

La raó principal per la qual encara podríeu voler utilitzar la manera per sortida d'establir el mode és si voleu barrejar i combinar dins del mateix workflow, _és a dir_, tenir algunes sortides copiades i algunes enllaçades simbòlicament.

Hi ha moltes altres opcions que podeu personalitzar d'aquesta manera, però esperem que això us doni una idea de la gamma d'opcions i com utilitzar-les eficaçment per adaptar-se a les vostres preferències.

### Conclusió

Sabeu com controlar la denominació i l'estructura dels directoris on es publiquen les vostres sortides, així com el mode de publicació de sortida del workflow.

### Què segueix?

Apreneu com adaptar la configuració del vostre workflow al vostre entorn de càlcul, començant amb la tecnologia d'empaquetament de programari.

---

## 3. Seleccionar una tecnologia d'empaquetament de programari

Fins ara hem estat mirant elements de configuració que controlen com entren les entrades i on surten les sortides. Ara és el moment de centrar-nos més específicament en adaptar la configuració del vostre workflow al vostre entorn de càlcul.

El primer pas en aquest camí és especificar d'on vindran els paquets de programari que s'executaran en cada pas.
Ja estan instal·lats a l'entorn de càlcul local?
Necessitem recuperar imatges i executar-les via un sistema de contenidors?
O necessitem recuperar paquets Conda i construir un entorn Conda local?

A la primera part d'aquest curs de formació (Parts 1-4) simplement vam utilitzar programari instal·lat localment al nostre workflow.
Després a la Part 5, vam introduir contenidors Docker i el fitxer `nextflow.config`, que vam utilitzar per habilitar l'ús de contenidors Docker.

Ara vegem com podem configurar una opció d'empaquetament de programari alternativa via el fitxer `nextflow.config`.

### 3.1. Deshabilitar Docker i habilitar Conda al fitxer de configuració

??? example "Escenari"

    Esteu movent el vostre pipeline a un clúster HPC on Docker no està permès per raons de seguretat.
    El clúster suporta Singularity i Conda, així que necessiteu canviar la vostra configuració en conseqüència.

Com s'ha assenyalat anteriorment, Nextflow suporta múltiples tecnologies de contenidors incloent Singularity (que s'utilitza més àmpliament en HPC), així com gestors de paquets de programari com Conda.

Podem canviar el nostre fitxer de configuració per utilitzar Conda en lloc de Docker.
Per fer-ho, canviem el valor de `docker.enabled` a `false`, i afegim una directiva habilitant l'ús de Conda:

=== "Després"

    ```groovy title="nextflow.config" linenums="1" hl_lines="1-2"
    docker.enabled = false
    conda.enabled = true
    ```

=== "Abans"

    ```groovy title="nextflow.config" linenums="1" hl_lines="1"
    docker.enabled = true
    ```

Això permetrà a Nextflow crear i utilitzar entorns Conda per a processos que tenen paquets Conda especificats.
El que significa que ara necessitem afegir un d'aquests al nostre procés `cowpy`!

### 3.2. Especificar un paquet Conda a la definició del procés

Ja hem recuperat l'URI per a un paquet Conda que conté l'eina `cowpy`: `conda-forge::cowpy==1.1.5`

Ara afegim l'URI a la definició del procés `cowpy` utilitzant la directiva `conda`:

=== "Després"

    ```groovy title="modules/cowpy.nf" linenums="4" hl_lines="4"
    process cowpy {

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
        conda 'conda-forge::cowpy==1.1.5'

        input:
    ```

=== "Abans"

    ```groovy title="modules/cowpy.nf" linenums="4"
    process cowpy {

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'

        input:
    ```

Per ser clars, no estem _reemplaçant_ la directiva `docker`, estem _afegint_ una opció alternativa.

!!! tip "Consell"

    Hi ha algunes maneres diferents d'obtenir l'URI per a un paquet conda donat.
    Recomanem utilitzar la consulta de cerca de [Seqera Containers](https://seqera.io/containers/), que us donarà un URI que podeu copiar i enganxar, fins i tot si no teniu previst crear un contenidor a partir d'ell.

### 3.3. Executar el workflow per verificar que pot utilitzar Conda

Provem-ho.

```bash
nextflow run 3-main.nf --batch conda
```

??? success "Sortida de la comanda"

    ```console title="Output"
    N E X T F L O W   ~  version 25.10.2

    Launching `3-main.nf` [trusting_lovelace] DSL2 - revision: 028a841db1

    executor >  local (8)
    [ee/4ca1f2] sayHello (3)       | 3 of 3 ✔
    [20/2596a7] convertToUpper (1) | 3 of 3 ✔
    [b3/e15de5] collectGreetings   | 1 of 1 ✔
    [c5/af5f88] cowpy              | 1 of 1 ✔
    ```

Això hauria de funcionar sense problemes i produir les mateixes sortides que anteriorment sota `results_config/conda`.

Entre bastidors, Nextflow ha recuperat els paquets Conda i ha creat l'entorn, cosa que normalment requereix una mica de feina; així que és agradable que no hàgim de fer res d'això nosaltres mateixos!

!!! info "Info"

    Això s'executa ràpidament perquè el paquet `cowpy` és força petit, però si esteu treballant amb paquets grans, pot trigar una mica més del normal la primera vegada, i podríeu veure la sortida de la consola quedar-se 'bloquejada' durant un minut o així abans de completar-se.
    Això és normal i es deu a la feina extra que Nextflow fa la primera vegada que utilitzeu un paquet nou.

Des del nostre punt de vista, sembla que funciona exactament igual que executar amb Docker, tot i que al backend la mecànica és una mica diferent.

Això significa que estem preparats per executar amb entorns Conda si cal.

??? info "Barrejar i combinar Docker i Conda"

    Com que aquestes directives s'assignen per procés, és possible 'barrejar i combinar', _és a dir_, configurar alguns dels processos del vostre workflow per executar-se amb Docker i altres amb Conda, per exemple, si la infraestructura de càlcul que esteu utilitzant suporta ambdós.
    En aquest cas, habilitaríeu tant Docker com Conda al vostre fitxer de configuració.
    Si ambdós estan disponibles per a un procés donat, Nextflow prioritzarà els contenidors.

    I com s'ha assenyalat anteriorment, Nextflow suporta múltiples altres tecnologies d'empaquetament de programari i contenidors, així que no esteu limitats només a aquests dos.

### Conclusió

Sabeu com configurar quin paquet de programari hauria d'utilitzar cada procés, i com canviar entre tecnologies.

### Què segueix?

Apreneu com canviar la plataforma d'execució utilitzada per Nextflow per fer realment la feina.

---

## 4. Seleccionar una plataforma d'execució

??? example "Escenari"

    Heu estat desenvolupant i provant el vostre pipeline al vostre portàtil, però ara necessiteu executar-lo en milers de mostres.
    La vostra institució té un clúster HPC amb un planificador Slurm que us agradaria utilitzar en canvi.

Fins ara, hem estat executant el nostre pipeline amb l'executor local.
Això executa cada tasca a la màquina on s'està executant Nextflow.
Quan Nextflow comença, mira les CPUs i memòria disponibles.
Si els recursos de les tasques preparades per executar-se excedeixen els recursos disponibles, Nextflow retindrà les últimes tasques de l'execució fins que una o més de les tasques anteriors hagin acabat, alliberant els recursos necessaris.

L'executor local és convenient i eficient, però està limitat a aquesta única màquina. Per a càrregues de treball molt grans, podeu descobrir que la vostra màquina local és un coll d'ampolla, ja sigui perquè teniu una sola tasca que requereix més recursos dels que teniu disponibles, o perquè teniu tantes tasques que esperar que una sola màquina les executi trigaria massa.

Nextflow suporta [molts backends d'execució diferents](https://nextflow.io/docs/latest/executor.html), incloent planificadors HPC (Slurm, LSF, SGE, PBS, Moab, OAR, Bridge, HTCondor i altres) així com backends d'execució al núvol (AWS Batch, Google Cloud Batch, Azure Batch, Kubernetes i més).

### 4.1. Apuntar a un backend diferent

L'elecció de l'executor s'estableix per una directiva de procés anomenada `executor`.
Per defecte està establert a `local`, així que la següent configuració està implícita:

```groovy title="Built-in configuration"
process {
    executor = 'local'
}
```

Per establir l'executor per apuntar a un backend diferent, simplement especificaríeu l'executor que voleu utilitzant una sintaxi similar a la descrita anteriorment per a assignacions de recursos (vegeu [Executors](https://nextflow.io/docs/latest/executor.html) per a totes les opcions).

```groovy title="nextflow.config"
process {
    executor = 'slurm'
}
```

!!! warning "Advertència"

    No podem provar això realment a l'entorn de formació perquè no està configurat per connectar-se a un HPC.

### 4.2. Tractar amb sintaxi específica del backend per a paràmetres d'execució

La majoria de plataformes de computació d'alt rendiment permeten (i de vegades requereixen) que especifiqueu certs paràmetres com ara sol·licituds i limitacions d'assignació de recursos (per exemple, nombre de CPUs i memòria) i nom de la cua de treballs a utilitzar.

Malauradament, cadascun d'aquests sistemes utilitza diferents tecnologies, sintaxis i configuracions per definir com s'ha de definir i enviar un treball al planificador rellevant.

??? abstract "Exemples"

    Per exemple, el mateix treball que requereix 8 CPUs i 4GB de RAM per ser executat a la cua "my-science-work" necessita ser expressat de les següents maneres diferents depenent del backend.

    ```bash title="Config for SLURM / submit using sbatch"
    #SBATCH -o /path/to/my/task/directory/my-task-1.log
    #SBATCH --no-requeue
    #SBATCH -c 8
    #SBATCH --mem 4096M
    #SBATCH -p my-science-work
    ```

    ```bash title="Config for PBS / submit using qsub"
    #PBS -o /path/to/my/task/directory/my-task-1.log
    #PBS -j oe
    #PBS -q my-science-work
    #PBS -l nodes=1:ppn=5
    #PBS -l mem=4gb
    ```

    ```bash title="Config for SGE / submit using qsub"
    #$ -o /path/to/my/task/directory/my-task-1.log
    #$ -j y
    #$ -terse
    #$ -notify
    #$ -q my-science-work
    #$ -l slots=5
    #$ -l h_rss=4096M,mem_free=4096M
    ```

Afortunadament, Nextflow simplifica tot això.
Proporciona una sintaxi estandarditzada perquè pugueu especificar les propietats rellevants com ara `cpus`, `memory` i `queue` només una vegada (vegeu [Directives de procés](https://nextflow.io/docs/latest/reference/process.html#process-directives) per a totes les opcions disponibles).
Després, en temps d'execució, Nextflow utilitzarà aquestes configuracions per generar els scripts específics del backend apropiats basats en la configuració de l'executor.

Cobrirem aquesta sintaxi estandarditzada a la següent secció.

### Conclusió

Ara sabeu com canviar l'executor per utilitzar diferents tipus d'infraestructura de càlcul.

### Què segueix?

Apreneu com avaluar i expressar assignacions i limitacions de recursos a Nextflow.

---

## 5. Controlar assignacions de recursos de càlcul

??? example "Escenari"

    El vostre pipeline continua fallant al clúster perquè les tasques estan sent eliminades per excedir els límits de memòria.
    O potser se us està cobrant per recursos que no esteu utilitzant i voleu optimitzar els costos.

La majoria de plataformes de computació d'alt rendiment permeten (i de vegades requereixen) que especifiqueu certs paràmetres d'assignació de recursos com ara nombre de CPUs i memòria.

Per defecte, Nextflow utilitzarà una sola CPU i 2GB de memòria per a cada procés.
Les directives de procés corresponents s'anomenen `cpus` i `memory`, així que la següent configuració està implícita:

```groovy title="Built-in configuration" linenums="1"
process {
    cpus = 1
    memory = 2.GB
}
```

Podeu modificar aquests valors, ja sigui per a tots els processos o per a processos específics amb nom, utilitzant directives de procés addicionals al vostre fitxer de configuració.
Nextflow les traduirà a les instruccions apropiades per a l'executor escollit.

Però com sabeu quins valors utilitzar?

### 5.1. Executar el workflow per generar un informe d'utilització de recursos

??? example "Escenari"

    No sabeu quanta memòria o CPU necessiten els vostres processos i voleu evitar malgastar recursos o que els treballs siguin eliminats.

Si no sabeu per endavant quanta CPU i memòria és probable que necessitin els vostres processos, podeu fer una mica de perfilat de recursos, és a dir, executar el workflow amb algunes assignacions per defecte, registrar quant va utilitzar cada procés, i a partir d'aquí, estimar com ajustar les assignacions base.

Convenientment, Nextflow inclou eines integrades per fer això, i generarà feliçment un informe per a vosaltres a petició.

Per fer-ho, afegiu `-with-report <nom_fitxer>.html` a la vostra línia de comandes.

```bash
nextflow run 3-main.nf -with-report report-config-1.html
```

L'informe és un fitxer html, que podeu descarregar i obrir al vostre navegador. També podeu fer clic dret sobre ell a l'explorador de fitxers de l'esquerra i fer clic a `Show preview` per veure'l a l'entorn de formació.

Preneu-vos uns minuts per mirar l'informe i veure si podeu identificar algunes oportunitats per ajustar recursos.
Assegureu-vos de fer clic a les pestanyes que mostren els resultats d'utilització com a percentatge del que es va assignar.

Vegeu [Reports](https://nextflow.io/docs/latest/reports.html) per a documentació sobre totes les funcionalitats disponibles.

### 5.2. Establir assignacions de recursos per a tots els processos

El perfilat mostra que els processos del nostre workflow de formació són molt lleugers, així que reduïm l'assignació de memòria per defecte a 1GB per procés.

Afegiu el següent al vostre fitxer `nextflow.config`, abans de la secció de paràmetres del pipeline:

```groovy title="nextflow.config" linenums="4"
/*
* Configuració de procés
*/
process {
    memory = 1.GB
}
```

Això ajudarà a reduir la quantitat de càlcul que consumim.

### 5.3. Establir assignacions de recursos per a un procés específic

Al mateix temps, farem veure que el procés `cowpy` requereix més recursos que els altres, només perquè puguem demostrar com ajustar assignacions per a un procés individual.

=== "Després"

    ```groovy title="nextflow.config" linenums="4" hl_lines="6-9"
    /*
    * Configuració de procés
    */
    process {
        memory = 1.GB
        withName: 'cowpy' {
            memory = 2.GB
            cpus = 2
        }
    }
    ```

=== "Abans"

    ```groovy title="nextflow.config" linenums="4"
    /*
    * Configuració de procés
    */
    process {
        memory = 1.GB
    }
    ```

Amb aquesta configuració, tots els processos sol·licitaran 1GB de memòria i una sola CPU (el valor per defecte implícit), excepte el procés `cowpy`, que sol·licitarà 2GB i 2 CPUs.

!!! info "Info"

    Si teniu una màquina amb poques CPUs i assigneu un nombre alt per procés, podríeu veure crides de procés posant-se en cua una darrere l'altra.
    Això és perquè Nextflow assegura que no sol·licitem més CPUs de les que estan disponibles.

### 5.4. Executar el workflow amb la configuració actualitzada

Provem això, proporcionant un nom de fitxer diferent per a l'informe de perfilat perquè puguem comparar el rendiment abans i després dels canvis de configuració.

```bash
nextflow run 3-main.nf -with-report report-config-2.html
```

Probablement no notareu cap diferència real ja que aquesta és una càrrega de treball tan petita, però aquest és l'enfocament que utilitzaríeu per analitzar el rendiment i els requisits de recursos d'un workflow del món real.

És molt útil quan els vostres processos tenen diferents requisits de recursos. Us permet ajustar correctament les assignacions de recursos que configureu per a cada procés basant-vos en dades reals, no en conjectures.

!!! tip "Consell"

    Això és només una petita mostra del que podeu fer per optimitzar l'ús de recursos.
    Nextflow mateix té una [lògica de reintent dinàmica](https://nextflow.io/docs/latest/process.html#dynamic-task-resources) realment interessant integrada per reintentar treballs que fallen a causa de limitacions de recursos.
    A més, la Plataforma Seqera ofereix eines impulsades per IA per optimitzar les vostres assignacions de recursos automàticament també.

### 5.5. Afegir límits de recursos

Depenent de quin executor de càlcul i infraestructura de càlcul estigueu utilitzant, pot haver-hi algunes restriccions sobre el que podeu (o heu) assignar.
Per exemple, el vostre clúster pot requerir que us mantingueu dins de certs límits.

Podeu utilitzar la directiva `resourceLimits` per establir les limitacions rellevants. La sintaxi es veu així quan està sola en un bloc de procés:

```groovy title="Syntax example"
process {
    resourceLimits = [
        memory: 750.GB,
        cpus: 200,
        time: 30.d
    ]
}
```

Nextflow traduirà aquests valors a les instruccions apropiades depenent de l'executor que hàgiu especificat.

No executarem això, ja que no tenim accés a infraestructura rellevant a l'entorn de formació.
No obstant això, si intentéssiu executar el workflow amb assignacions de recursos que excedeixen aquests límits, després miréssiu la comanda `sbatch` al fitxer de script `.command.run`, veuríeu que les sol·licituds que realment s'envien a l'executor estan limitades als valors especificats per `resourceLimits`.

??? info "Configuracions de referència institucionals"

    El projecte nf-core ha compilat una [col·lecció de fitxers de configuració](https://nf-co.re/configs/) compartits per diverses institucions arreu del món, cobrint una àmplia gamma d'executors HPC i núvol.

    Aquestes configuracions compartides són valuoses tant per a persones que hi treballen i per tant poden simplement utilitzar la configuració de la seva institució directament, com a model per a persones que busquen desenvolupar una configuració per a la seva pròpia infraestructura.

### Conclusió

Sabeu com generar un informe de perfilat per avaluar la utilització de recursos i com modificar assignacions de recursos per a tots els processos i/o per a processos individuals, així com establir limitacions de recursos per executar en HPC.

### Què segueix?

Apreneu com configurar perfils de configuració predefinits i canviar entre ells en temps d'execució.

---

## 6. Utilitzar perfils per canviar entre configuracions predefinides

??? example "Escenari"

    Canvieu regularment entre executar pipelines al vostre portàtil per a desenvolupament i al HPC de la vostra institució per a execucions de producció.
    Esteu cansats de canviar manualment la configuració cada vegada que canvieu d'entorn.

Us hem mostrat diverses maneres en què podeu personalitzar la configuració del vostre pipeline depenent del projecte en què estigueu treballant o de l'entorn de càlcul que estigueu utilitzant.

Potser voleu canviar entre configuracions alternatives depenent de quina infraestructura de càlcul estigueu utilitzant. Per exemple, podríeu voler desenvolupar i executar proves a petita escala localment al vostre portàtil, després executar càrregues de treball a escala completa en HPC o núvol.

Nextflow us permet configurar qualsevol nombre de [**perfils**](https://nextflow.io/docs/latest/config.html#profiles) que descriuen diferents configuracions, que després podeu seleccionar en temps d'execució utilitzant un argument de línia de comandes, en lloc d'haver de modificar el fitxer de configuració en si.

### 6.1. Crear perfils per canviar entre desenvolupament local i execució en HPC

Configurem dos perfils alternatius; un per executar càrregues a petita escala en un ordinador normal, on utilitzarem contenidors Docker, i un per executar en un HPC universitari amb un planificador Slurm, on utilitzarem paquets Conda.

#### 6.1.1. Configurar els perfils

Afegiu el següent al vostre fitxer `nextflow.config`, després de la secció de paràmetres del pipeline però abans de la configuració de sortida:

```groovy title="nextflow.config" linenums="24"
/*
* Perfils
*/
profiles {
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

Veieu que per al HPC universitari, també estem especificant limitacions de recursos.

#### 6.1.2. Executar el workflow amb un perfil

Per especificar un perfil a la nostra línia de comandes de Nextflow, utilitzem l'argument `-profile`.

Provem d'executar el workflow amb la configuració `my_laptop`.

```bash
nextflow run 3-main.nf -profile my_laptop
```

??? success "Sortida de la comanda"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `3-main.nf` [gigantic_brazil] DSL2 - revision: ede9037d02

    executor >  local (8)
    [58/da9437] sayHello (3)       | 3 of 3 ✔
    [35/9cbe77] convertToUpper (2) | 3 of 3 ✔
    [67/857d05] collectGreetings   | 1 of 1 ✔
    [37/7b51b5] cowpy              | 1 of 1 ✔
    ```

Com podeu veure, això ens permet alternar entre configuracions molt convenientment en temps d'execució.

!!! warning "Advertència"

    El perfil `univ_hpc` no s'executarà correctament a l'entorn de formació ja que no tenim accés a un planificador Slurm.

Si en el futur trobem altres elements de configuració que sempre co-ocorren amb aquests, podem simplement afegir-los al(s) perfil(s) corresponent(s).
També podem crear perfils addicionals si hi ha altres elements de configuració que volem agrupar junts.

### 6.2. Crear un perfil de paràmetres de prova

??? example "Escenari"

    Voleu que altres puguin provar el vostre pipeline ràpidament sense haver de reunir les seves pròpies dades d'entrada.

Els perfils no són només per a configuració d'infraestructura.
També podem utilitzar-los per establir valors per defecte per a paràmetres del workflow, per facilitar que altres provin el workflow sense haver de reunir valors d'entrada apropiats ells mateixos.
Podeu considerar això una alternativa a utilitzar un fitxer de paràmetres.

#### 6.2.1. Configurar el perfil

La sintaxi per expressar valors per defecte en aquest context es veu així, per a un perfil que anomenem `test`:

```groovy title="Syntax example"
    test {
        params.<parameter1>
        params.<parameter2>
        ...
    }
```

Si afegim un perfil de prova per al nostre workflow, el bloc `profiles` es converteix en:

```groovy title="nextflow.config" linenums="24"
/*
* Perfils
*/
profiles {
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
    test {
        params.input = 'data/greetings.csv'
        params.batch = 'test'
        params.character = 'dragonandcow'
    }
}
```

Igual que per als perfils de configuració tècnica, podeu configurar múltiples perfils diferents especificant paràmetres sota qualsevol nom arbitrari que vulgueu.

#### 6.2.2. Executar el workflow localment amb el perfil de prova

Convenientment, els perfils no són mútuament excloents, així que podem especificar múltiples perfils a la nostra línia de comandes utilitzant la següent sintaxi `-profile <perfil1>,<perfil2>` (per a qualsevol nombre de perfils).

Si combineu perfils que estableixen valors per als mateixos elements de configuració i estan descrits al mateix fitxer de configuració, Nextflow resoldrà el conflicte utilitzant el valor que hagi llegit per últim (_és a dir_, el que vingui més tard al fitxer).
Si les configuracions conflictives estan establertes en diferents fonts de configuració, s'aplica l'[ordre de precedència](https://www.nextflow.io/docs/latest/config.html) per defecte.

Provem d'afegir el perfil de prova a la nostra comanda anterior:

```bash
nextflow run 3-main.nf -profile my_laptop,test
```

??? success "Sortida de la comanda"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `3-main.nf` [jovial_coulomb] DSL2 - revision: 46a6763141

    executor >  local (8)
    [9b/687cdc] sayHello (2)       | 3 of 3 ✔
    [ca/552187] convertToUpper (3) | 3 of 3 ✔
    [e8/83e306] collectGreetings   | 1 of 1 ✔
    [fd/e84fa9] cowpy              | 1 of 1 ✔
    ```

Això utilitzarà Docker on sigui possible i produirà sortides sota `results_config/test`, i aquesta vegada el caràcter és el duo còmic `dragonandcow`.

??? abstract "Contingut del fitxer"

    ```console title="results_config/test/"
     _________
    / HOLà    \
    | HELLO   |
    \ BONJOUR /
    ---------
                \                    ^    /^
                  \                  / \  // \
                  \   |\___/|      /   \//  .\
                    \  /O  O  \__  /    //  | \ \           *----*
                      /     /  \/_/    //   |  \  \          \   |
                      \@___\@`    \/_   //    |   \   \         \/\ \
                    0/0/|       \/_ //     |    \    \         \ \
                0/0/0/0/|        \///      |     \     \       | |
              0/0/0/0/0/_|_ /   (  //       |      \     _\     |  /
          0/0/0/0/0/0/`/,_ _ _/  ) ; -.    |    _ _\.-~       /   /
                      ,-}        _      *-.|.-~-.           .~    ~
      \     \__/        `/\      /                 ~-. _ .-~      /
      \____(oo)           *.   }            {                   /
      (    (--)          .----~-.\        \-`                 .~
      //__\\  \__ Ack!   ///.----..<        \             _ -~
      //    \\               ///-._ _ _ _ _ _ _{^ - - - - ~
    ```

Això significa que mentre distribuïm qualsevol fitxer de dades de prova amb el codi del workflow, qualsevol pot provar ràpidament el workflow sense haver de proporcionar les seves pròpies entrades via línia de comandes o un fitxer de paràmetres.

!!! tip "Consell"

    Podem apuntar a URLs per a fitxers més grans que estan emmagatzemats externament.
    Nextflow els descarregarà automàticament sempre que hi hagi una connexió oberta.

    Per a més detalls, vegeu la Missió Secundària [Working with Files](../side_quests/working_with_files.md)

### 6.3. Utilitzar `nextflow config` per veure la configuració resolta

Com s'ha assenyalat anteriorment, de vegades el mateix paràmetre pot estar establert a diferents valors en perfils que voleu combinar.
I més generalment, hi ha nombrosos llocs on els elements de configuració poden estar emmagatzemats, i de vegades les mateixes propietats poden estar establertes a diferents valors en diferents llocs.

Nextflow aplica un [ordre de precedència](https://nextflow.io/docs/latest/config.html#configuration-file) establert per resoldre qualsevol conflicte, però això pot ser complicat de determinar vosaltres mateixos.
I fins i tot si res està en conflicte, pot ser tediós buscar tots els llocs possibles on les coses podrien estar configurades.

Afortunadament, Nextflow inclou una eina d'utilitat convenient anomenada `config` que pot automatitzar tot aquest procés per a vosaltres.

L'eina `config` explorarà tots els continguts del vostre directori de treball actual, recollirà qualsevol fitxer de configuració, i produirà la configuració completament resolta que Nextflow utilitzaria per executar el workflow.
Això us permet esbrinar quines configuracions s'utilitzaran sense haver de llançar res.

#### 6.3.1. Resoldre la configuració per defecte

Executeu aquesta comanda per resoldre la configuració que s'aplicaria per defecte.

```bash
nextflow config
```

??? success "Sortida de la comanda"

    ```groovy
    params {
      input = 'data/greetings.csv'
      batch = 'batch'
      character = 'turkey'
    }

    docker {
      enabled = false
    }

    conda {
      enabled = true
    }

    process {
      memory = '1 GB'
      withName:cowpy {
          memory = '2 GB'
          cpus = 2
      }
    }

    outputDir = 'results_config/batch'

    workflow {
      output {
          mode = 'copy'
      }
    }
    ```

Això us mostra la configuració base que obteniu si no especifiqueu res extra a la línia de comandes.

#### 6.3.2. Resoldre la configuració amb configuracions específiques activades

Si proporcioneu paràmetres de línia de comandes, per exemple habilitant un o més perfils o carregant un fitxer de paràmetres, la comanda addicionalment els tindrà en compte.

```bash
nextflow config -profile my_laptop,test
```

??? success "Sortida de la comanda"

    ```groovy
    params {
      input = 'data/greetings.csv'
      batch = 'test'
      character = 'dragonandcow'
    }

    docker {
      enabled = true
    }

    conda {
      enabled = true
    }

    process {
      memory = '1 GB'
      withName:cowpy {
          memory = '2 GB'
          cpus = 2
      }
      executor = 'local'
    }

    outputDir = 'results_config/test'

    workflow {
      output {
          mode = 'copy'
      }
    }
    ```

Això és especialment útil per a projectes complexos que impliquen múltiples capes de configuració.

### Conclusió

Sabeu com utilitzar perfils per seleccionar una configuració predefinida en temps d'execució amb un mínim d'inconvenients.
Més generalment, sabeu com configurar les vostres execucions de workflow per adaptar-se a diferents plataformes de càlcul i millorar la reproduïbilitat de les vostres anàlisis.

### Què segueix?

Apreneu com executar pipelines directament des de repositoris remots com GitHub.

---

## 7. Executar pipelines des de repositoris remots

??? example "Escenari"

    Voleu executar un pipeline ben establert com els de nf-core sense haver de descarregar i gestionar el codi vosaltres mateixos.

Fins ara hem estat executant scripts de workflow ubicats al directori actual.
A la pràctica, sovint voldreu executar pipelines emmagatzemats en repositoris remots, com ara GitHub.

Nextflow fa això senzill: podeu executar qualsevol pipeline directament des d'una URL de repositori Git sense descarregar-lo manualment primer.

### 7.1. Executar un pipeline des de GitHub

La sintaxi bàsica per executar un pipeline remot és `nextflow run <repositori>`, on `<repositori>` pot ser un camí de repositori GitHub com `nextflow-io/hello`, una URL completa, o un camí a GitLab, Bitbucket, o altres serveis d'allotjament Git.

Proveu d'executar el pipeline de demostració oficial "hello" de Nextflow:

```bash
nextflow run nextflow-io/hello
```

??? success "Sortida de la comanda"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Pulling nextflow-io/hello ...
     downloaded from https://github.com/nextflow-io/hello.git
    Launching `https://github.com/nextflow-io/hello` [sleepy_swanson] DSL2 - revision: 2ce0b0e294 [master]

    executor >  local (4)
    [ba/08236d] sayHello (4) [100%] 4 of 4 ✔
    Ciao world!

    Hello world!

    Bonjour world!

    Hola world!
    ```

La primera vegada que executeu un pipeline remot, Nextflow el descarrega i el guarda en memòria cau localment.
Les execucions posteriors utilitzen la versió en memòria cau tret que sol·liciteu explícitament una actualització.

### 7.2. Especificar una versió per a reproduïbilitat

Per defecte, Nextflow executa la versió més recent de la branca per defecte.
Podeu especificar una versió particular (etiqueta), branca, o commit utilitzant la bandera `-r`:

```bash
nextflow run nextflow-io/hello -r v1.3
```

??? success "Sortida de la comanda"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `https://github.com/nextflow-io/hello` [sick_carson] DSL2 - revision: 2ce0b0e294 [v1.3]

    executor >  local (4)
    [61/e11f77] sayHello (4) [100%] 4 of 4 ✔
    Ciao world!

    Bonjour world!

    Hello world!

    Hola world!
    ```

Especificar versions exactes és essencial per a la reproduïbilitat.

### Conclusió

Sabeu com executar pipelines directament des de GitHub i altres repositoris remots, i com especificar versions per a reproduïbilitat.

### Què segueix?

Doneu-vos una gran palmada a l'esquena!
Sabeu tot el que necessiteu saber per començar a executar i gestionar pipelines de Nextflow.

Això conclou aquest curs, però si esteu ansiosos per continuar aprenent, tenim dues recomanacions principals:

- Si voleu aprofundir en el desenvolupament dels vostres propis pipelines, feu una ullada a [Hello Nextflow](../hello_nextflow/index.md), un curs per a principiants que cobreix la mateixa progressió general que aquest però entra en molt més detall sobre canals i operadors.
- Si us agradaria continuar aprenent com executar pipelines de Nextflow sense aprofundir més en el codi, feu una ullada a la primera part de [Hello nf-core](../hello_nf-core/index.md), que introdueix les eines per trobar i executar pipelines del projecte [nf-core](https://nf-co.re/) enormement popular.

Que us divertiu!

---

## Qüestionari

<quiz>
Quan els valors de paràmetres s'estableixen tant al fitxer del workflow com a `nextflow.config`, quin té precedència?
- [ ] El valor del fitxer del workflow
- [x] El valor del fitxer de configuració
- [ ] El primer valor trobat
- [ ] Causa un error

Més informació: [1.1. Configurar valors a `nextflow.config`](#11-set-up-values-in-nextflowconfig)
</quiz>

<quiz>
Quina és la diferència de sintaxi entre establir un valor per defecte de paràmetre en un fitxer de workflow vs. un fitxer de configuració?
- [ ] Utilitzen la mateixa sintaxi
- [x] El workflow utilitza declaració tipada (`#!groovy param: Type = value`), la configuració utilitza assignació (`#!groovy param = value`)
- [ ] La configuració utilitza declaració tipada, el workflow utilitza assignació
- [ ] Només els fitxers de configuració poden establir valors per defecte

Més informació: [1.1. Configurar valors a `nextflow.config`](#11-set-up-values-in-nextflowconfig)
</quiz>

<quiz>
Com especifiqueu un fitxer de paràmetres quan executeu un workflow?
- [ ] `--params params.yaml`
- [ ] `-config params.yaml`
- [x] `-params-file params.yaml`
- [ ] `--input-params params.yaml`

Més informació: [1.3. Utilitzar un fitxer de paràmetres](#13-use-a-parameter-file)
</quiz>

<quiz>
Què controla l'opció de configuració `outputDir`?
- [ ] La ubicació del directori de treball
- [x] El camí base on es publiquen les sortides del workflow
- [ ] El directori per a fitxers de registre
- [ ] La ubicació dels fitxers de mòdul

Més informació: [2.1. Personalitzar el nom del directori outputDir](#21-customize-the-outputdir-directory-name)
</quiz>

<quiz>
Com referencieu un nom de procés dinàmicament a la configuració del camí de sortida?
- [ ] `#!groovy ${processName}`
- [ ] `#!groovy path "<process>.name"`
- [x] `#!groovy path { <process>.name }`
- [ ] `@processName`

Més informació: [2.2. Organitzar sortides per procés](#22-organize-outputs-by-process)
</quiz>

<quiz>
Si tant Docker com Conda estan habilitats i un procés té ambdues directives, quina es prioritza?
- [x] Docker (contenidors)
- [ ] Conda
- [ ] La primera definida al procés
- [ ] Causa un error

Més informació: [3. Seleccionar una tecnologia d'empaquetament de programari](#3-select-a-software-packaging-technology)
</quiz>

<quiz>
Quin és l'executor per defecte a Nextflow?
- [x] `local`
- [ ] `slurm`
- [ ] `kubernetes`
- [ ] `aws`

Més informació: [4. Seleccionar una plataforma d'execució](#4-select-an-execution-platform)
</quiz>

<quiz>
Quina comanda genera un informe d'utilització de recursos?
- [ ] `nextflow run workflow.nf -with-metrics`
- [ ] `nextflow run workflow.nf -with-stats`
- [x] `nextflow run workflow.nf -with-report report.html`
- [ ] `nextflow run workflow.nf -profile report`

Més informació: [5.1. Executar el workflow per generar un informe d'utilització de recursos](#51-run-the-workflow-to-generate-a-resource-utilization-report)
</quiz>

<quiz>
Com establiu requisits de recursos per a un procés específic anomenat `cowpy` al fitxer de configuració?
- [ ] `#!groovy cowpy.memory = '2.GB'`
- [ ] `#!groovy process.cowpy.memory = '2.GB'`
- [x] `#!groovy process { withName: 'cowpy' { memory = '2.GB' } }`
- [ ] `#!groovy resources.cowpy.memory = '2.GB'`

Més informació: [5.3. Establir assignacions de recursos per a un procés específic](#53-set-resource-allocations-for-a-specific-process)
</quiz>

<quiz>
Què fa la directiva `resourceLimits`?
- [ ] Estableix requisits mínims de recursos
- [ ] Assigna recursos als processos
- [x] Limita els recursos màxims que es poden sol·licitar
- [ ] Monitoritza l'ús de recursos en temps real

Més informació: [5.5. Afegir límits de recursos](#55-add-resource-limits)
</quiz>

<quiz>
Com especifiqueu múltiples perfils en una sola comanda?
- [ ] `-profile profile1 -profile profile2`
- [ ] `-profiles profile1,profile2`
- [x] `-profile profile1,profile2`
- [ ] `--profile profile1 --profile profile2`

Més informació: [6. Utilitzar perfils per canviar entre configuracions predefinides](#6-use-profiles-to-switch-between-preset-configurations)
</quiz>

<quiz>
Quina comanda mostra la configuració completament resolta que Nextflow utilitzaria?
- [ ] `nextflow show-config`
- [ ] `nextflow settings`
- [x] `nextflow config`
- [ ] `nextflow resolve`

Més informació: [6.3. Utilitzar `nextflow config` per veure la configuració resolta](#63-use-nextflow-config-to-see-the-resolved-configuration)
</quiz>

<quiz>
Per a què es poden utilitzar els perfils? (Seleccioneu totes les que apliquin)
- [x] Definir configuracions específiques d'infraestructura (executors, contenidors)
- [x] Establir límits de recursos per a diferents entorns
- [x] Proporcionar paràmetres de prova per a proves fàcils del workflow
- [ ] Definir nous processos

Més informació: [6. Utilitzar perfils per canviar entre configuracions predefinides](#6-use-profiles-to-switch-between-preset-configurations)
</quiz>
