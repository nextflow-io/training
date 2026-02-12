# Part 6: Hello Config

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducció assistida per IA - [més informació i suggeriments](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/FcZTiE25TeA?si=y8lAedhEHWaTV4zd&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

/// caption
:fontawesome-brands-youtube:{ .youtube } Vegeu [la llista de reproducció completa](https://youtube.com/playlist?list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&si=eF7cLR62goy-lc6n) al canal de YouTube de Nextflow.

:green_book: La transcripció del vídeo està disponible [aquí](./transcripts/06_hello_config.md).
///

Aquesta secció explorarà com configurar i gestionar la configuració del vostre pipeline de Nextflow perquè pugueu personalitzar el seu comportament, adaptar-lo a diferents entorns i optimitzar l'ús de recursos _sense alterar ni una sola línia del codi del workflow_.

Hi ha diverses maneres de fer-ho, que es poden utilitzar en combinació i s'interpreten segons l'[ordre de precedència](https://nextflow.io/docs/latest/config.html) descrit a la documentació de configuració.

En aquesta part del curs, us mostrarem el mecanisme de fitxer de configuració més senzill i comú, el fitxer [`nextflow.config`](https://nextflow.io/docs/latest/config.html), que ja vau trobar a la Part 5: Hello Containers.

Repassarem components essencials de la configuració de Nextflow com ara directives de procés, executors, perfils i fitxers de paràmetres.
Aprenent a utilitzar aquestes opcions de configuració de manera efectiva, podeu millorar la flexibilitat, escalabilitat i rendiment dels vostres pipelines.

??? info "Com començar des d'aquesta secció"

    Aquesta secció del curs assumeix que heu completat les Parts 1-5 del curs [Hello Nextflow](./index.md) i teniu un pipeline complet i funcional.

    Si esteu començant el curs des d'aquest punt, haureu de copiar el directori `modules` i el fitxer `nextflow.config` des de les solucions:

    ```bash
    cp -r solutions/5-hello-containers/modules .
    cp solutions/5-hello-containers/nextflow.config .
    ```

    El fitxer `nextflow.config` conté la línia `docker.enabled = true` que habilita l'ús de contenidors Docker.

    Si no esteu familiaritzats amb el pipeline Hello o us aniria bé un recordatori, consulteu [aquesta pàgina d'informació](../info/hello_pipeline.md).

---

## 0. Escalfament: Executeu `hello-config.nf`

Utilitzarem l'script de workflow `hello-config.nf` com a punt de partida.
És equivalent a l'script produït en treballar la Part 5 d'aquest curs de formació, excepte que hem canviat les destinacions de sortida:

```groovy title="hello-config.nf" linenums="37" hl_lines="3 7 11 15"
output {
    first_output {
        path 'hello_config/intermediates'
        mode 'copy'
    }
    uppercased {
        path 'hello_config/intermediates'
        mode 'copy'
    }
    collected {
        path 'hello_config/intermediates'
        mode 'copy'
    }
    batch_report {
        path 'hello_config'
        mode 'copy'
    }
    cowpy_art {
        path 'hello_config'
        mode 'copy'
    }
}
```

Només per assegurar-nos que tot funciona, executeu l'script una vegada abans de fer cap canvi:

```bash
nextflow run hello-config.nf
```

??? success "Sortida de la comanda"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-config.nf` [nice_escher] DSL2 - revision: d5dfdc9872

    executor > local (7)
    [6a/bc46a6] sayHello (2) [100%] 3 of 3 ✔
    [33/67bc48] convertToUpper (3) [100%] 3 of 3 ✔
    [b5/de03ba] collectGreetings [100%] 1 of 1 ✔
    [98/c6b57b] cowpy              | 1 of 1 ✔
    ```

Com abans, trobareu els fitxers de sortida al directori especificat al bloc `output` (`results/hello_config/`).

??? abstract "Contingut del directori"

    ```console
    results/hello_config/
    ├── cowpy-COLLECTED-batch-output.txt
    ├── intermediates
    │   ├── Bonjour-output.txt
    │   ├── COLLECTED-batch-output.txt
    │   ├── Hello-output.txt
    │   ├── Holà-output.txt
    │   ├── UPPER-Bonjour-output.txt
    │   ├── UPPER-Hello-output.txt
    │   └── UPPER-Holà-output.txt
    └── batch-report.txt
    ```

La sortida final d'art ASCII es troba al directori `results/hello_config/`, amb el nom `cowpy-COLLECTED-batch-output.txt`.

??? abstract "Contingut del fitxer"

    ```console title="results/hello_config/cowpy-COLLECTED-batch-output.txt"
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

Si això us ha funcionat, esteu preparats per aprendre a configurar els vostres pipelines.

---

## 1. Gestionar els paràmetres d'entrada del workflow

Començarem amb un aspecte de la configuració que és simplement una extensió del que hem estat treballant fins ara: la gestió dels paràmetres d'entrada.

Actualment, el nostre workflow està configurat per acceptar diversos valors de paràmetres via línia de comandes, amb valors per defecte establerts en un bloc `params` a l'script del workflow mateix.
Tanmateix, potser voleu sobreescriure aquests valors per defecte sense haver d'especificar paràmetres a la línia de comandes o modificar el fitxer d'script original.

Hi ha diverses maneres de fer-ho; us mostrarem tres maneres bàsiques que s'utilitzen molt comunament.

### 1.1. Moure els valors per defecte a `nextflow.config`

Aquest és l'enfocament més senzill, tot i que possiblement el menys flexible ja que el fitxer principal `nextflow.config` no és quelcom que vulgueu estar editant per a cada execució.
Però té l'avantatge de separar les preocupacions de _declarar_ els paràmetres al workflow (que definitivament hi pertany) versus proporcionar _valors per defecte_, que estan més a casa en un fitxer de configuració.

Fem-ho en dos passos.

#### 1.1.1. Crear un bloc `params` al fitxer de configuració

Feu els següents canvis de codi al fitxer `nextflow.config`:

=== "Després"

    ```groovy title="nextflow.config" linenums="1" hl_lines="3-10"
    docker.enabled = true

    /*
    * Pipeline parameters
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
La sintaxi és una mica diferent.
Al fitxer de workflow, són declaracions tipades.
A la configuració, són assignacions de valors.

Tècnicament, això és suficient per sobreescriure els valors per defecte encara especificats al fitxer de workflow.
Podríeu modificar el caràcter, per exemple, i executar el workflow per satisfer-vos que el valor establert al fitxer de configuració sobreescriu el del fitxer de workflow.

Però amb l'esperit de moure la configuració completament al fitxer de configuració, eliminem aquests valors del fitxer de workflow completament.

#### 1.1.2. Eliminar els valors del bloc `params` al fitxer de workflow

Feu els següents canvis de codi al fitxer de workflow `hello-config.nf`:

=== "Després"

    ```groovy title="hello-config.nf" linenums="9" hl_lines="5-7"
    /*
    * Pipeline parameters
    */
    params {
        input: Path
        batch: String
        character: String
    }
    ```

=== "Abans"

    ```groovy title="hello-config.nf" linenums="9" hl_lines="5-7"
    /*
    * Pipeline parameters
    */
    params {
        input: Path = 'data/greetings.csv'
        batch: String = 'batch'
        character: String = 'turkey'
    }
    ```

Ara el fitxer de workflow mateix no estableix cap valor per defecte per a aquests paràmetres.

#### 1.1.3. Executar el pipeline

Provem que funciona correctament.

```bash
nextflow run hello-config.nf
```

??? success "Sortida de la comanda"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-config.nf` [disturbed_einstein] DSL2 - revision: ede9037d02

    executor >  local (8)
    [f0/35723c] sayHello (2)       | 3 of 3 ✔
    [40/3efd1a] convertToUpper (3) | 3 of 3 ✔
    [17/e97d32] collectGreetings   | 1 of 1 ✔
    [98/c6b57b] cowpy              | 1 of 1 ✔
    ```

Això encara produeix la mateixa sortida que abans.

La sortida final d'art ASCII es troba al directori `results/hello_config/`, amb el nom `cowpy-COLLECTED-batch-output.txt`, igual que abans.

??? abstract "Contingut del fitxer"

    ```console title="results/hello_config/cowpy-COLLECTED-batch-output.txt"
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

Això és genial, però de vegades potser voleu executar alguns experiments temporals amb diferents valors per defecte sense tocar el fitxer de configuració principal.
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
nextflow run ../hello-config.nf
```

??? success "Sortida de la comanda"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `../hello-config.nf` [trusting_escher] DSL2 - revision: 356df0818d

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

    ```console title="tux-run/results/hello_config/cowpy-COLLECTED-experiment-output.txt"
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

L'enfocament del subdirectori funciona molt bé per experimentar, però implica una mica de configuració i requereix que adapteu els camins en conseqüència.
Hi ha un enfocament més senzill per quan voleu executar el vostre pipeline amb un conjunt específic de valors, o permetre que algú altre ho faci amb un esforç mínim.

Nextflow ens permet especificar paràmetres via un [fitxer de paràmetres](https://nextflow.io/docs/latest/config.html#params-file) en format YAML o JSON, cosa que fa molt convenient gestionar i distribuir conjunts alternatius de valors per defecte, per exemple, així com valors de paràmetres específics per a l'execució.

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
nextflow run hello-config.nf -params-file test-params.yaml
```

??? success "Sortida de la comanda"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-config.nf` [disturbed_sammet] DSL2 - revision: ede9037d02

    executor >  local (8)
    [f0/35723c] sayHello (2)       | 3 of 3 ✔
    [40/3efd1a] convertToUpper (3) | 3 of 3 ✔
    [17/e97d32] collectGreetings   | 1 of 1 ✔
    [98/c6b57b] cowpy              | 1 of 1 ✔
    ```

El fitxer de sortida final hauria de contenir el caràcter stegosaurus dient les salutacions.

??? abstract "Contingut del fitxer"

    ```console title="results/hello_config/cowpy-COLLECTED-yaml-output.txt"
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
En aquests casos, utilitzar un fitxer de paràmetres ens permetrà proporcionar valors de paràmetres en temps d'execució sense haver d'escriure línies de comandes massives i sense modificar l'script del workflow.

També facilita distribuir conjunts de paràmetres a col·laboradors, o com a informació de suport per a una publicació, per exemple.
Això fa que el vostre treball sigui més reproduïble per altres.

### Conclusió

Sabeu com aprofitar les opcions de configuració clau per gestionar les entrades del workflow.

### Què segueix?

Apreneu a gestionar on i com es publiquen les sortides del vostre workflow.

---

## 2. Gestionar les sortides del workflow

Fins ara hem estat codificant tots els camins per a les declaracions de sortida a nivell de workflow, i com vam notar quan vam començar a afegir múltiples sortides, pot haver-hi una mica de repetició implicada.

Vegem algunes maneres comunes en què podríeu configurar això per ser més flexible.

### 2.1. Personalitzar el directori de sortida amb `-output-dir`

Quan controlem com s'organitzen les nostres sortides 'publicades' tenim dues prioritats diferents:

- El directori de sortida de nivell superior
- Com s'organitzen els fitxers dins d'aquest directori

Hem estat utilitzant el directori de nivell superior per defecte fins ara: `results`.
Comencem personalitzant això, utilitzant l'opció CLI `-output-dir`.

#### 2.1.1. Executar el pipeline amb `-output-dir`

L'opció `-output-dir` (abreviatura: `-o`) sobreescriu el directori de sortida per defecte (`results/`) per a totes les sortides del workflow.
Aquesta és la manera recomanada de controlar el camí arrel on es publiquen les sortides.

```bash
nextflow run hello-config.nf -output-dir custom-outdir-cli/
```

??? success "Sortida de la comanda"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-config.nf` [prickly_kay] DSL2 - revision: 32ecc4fba2

    executor >  local (8)
    [9f/332636] sayHello (1)       [100%] 3 of 3 ✔
    [03/a55991] convertToUpper (3) [100%] 3 of 3 ✔
    [e5/ab7893] collectGreetings   [100%] 1 of 1 ✔
    [a8/97338e] cowpy              [100%] 1 of 1 ✔
    ```

Això publica les sortides a `custom-outdir-cli/` en lloc de `results/`:

??? abstract "Contingut del directori"

    ```console
    custom-outdir-cli/
    └── hello_config
        ├── batch-report.txt
        ├── cowpy-COLLECTED-batch-output.txt
        └── intermediates
            ├── Bonjour-output.txt
            ├── COLLECTED-batch-output.txt
            ├── Hello-output.txt
            ├── Holà-output.txt
            ├── UPPER-Bonjour-output.txt
            ├── UPPER-Hello-output.txt
            └── UPPER-Holà-output.txt
    ```

Tingueu en compte que encara tenim el subdirectori `hello_config` de les declaracions `path` al bloc output.
Netegem això.

#### 2.1.2. Eliminar els camins codificats del bloc output

El prefix `hello_config/` estava codificat en capítols anteriors, però com ara estem aprenent a configurar camins de sortida de manera flexible, podem eliminar aquesta codificació.
Per a sortides que no necessiten un subdirectori podem establir la directiva `path` a una cadena buida, o eliminar-la completament.

Feu els següents canvis de codi al fitxer de workflow:

=== "Després"

    ```groovy title="hello-config.nf" linenums="42" hl_lines="3 7 11 15 19"
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

    ```groovy title="hello-config.nf" linenums="42" hl_lines="3 7 11 15 19"
    output {
        first_output {
            path 'hello_config/intermediates'
            mode 'copy'
        }
        uppercased {
            path 'hello_config/intermediates'
            mode 'copy'
        }
        collected {
            path 'hello_config/intermediates'
            mode 'copy'
        }
        batch_report {
            path 'hello_config'
            mode 'copy'
        }
        cowpy_art {
            path 'hello_config'
            mode 'copy'
        }
    }
    ```

Executeu el pipeline de nou:

```bash
nextflow run hello-config.nf -output-dir custom-outdir-cli-2/
```

Ara les sortides es publiquen directament sota `custom-outdir-cli-2/`, sense el subdirectori `hello_config`:

??? abstract "Contingut del directori"

    ```console
    custom-outdir-cli-2/
    ├── batch-report.txt
    ├── cowpy-COLLECTED-batch-output.txt
    └── intermediates
        ├── Bonjour-output.txt
        ├── COLLECTED-batch-output.txt
        ├── Hello-output.txt
        ├── Holà-output.txt
        ├── UPPER-Bonjour-output.txt
        ├── UPPER-Hello-output.txt
        └── UPPER-Holà-output.txt
    ```

!!! tip "Consell"

    L'opció `-output-dir` s'utilitza per controlar _on_ van les sortides, mentre que la directiva `path` al bloc output controla l'_estructura de subdirectoris_.

### 2.2. Camins de sortida dinàmics

A més de canviar el directori de sortida via CLI, també podem establir un valor per defecte personalitzat al fitxer de configuració utilitzant `outputDir`.
Això ens permet establir el camí del directori dinàmicament - no només utilitzant cadenes estàtiques.

#### 2.2.1. Establir `outputDir` al fitxer de configuració

Afegiu el següent codi al fitxer `nextflow.config`:

=== "Després"

    ```groovy title="nextflow.config" linenums="9" hl_lines="10-13"
    /*
    * Pipeline parameters
    */
    params {
        input = 'data/greetings.csv'
        batch = 'batch'
        character = 'turkey'
    }

    /*
    * Output settings
    */
    outputDir = "custom-outdir-config/${params.batch}"
    ```

=== "Abans"

    ```groovy title="nextflow.config" linenums="9"
    /*
    * Pipeline parameters
    */
    params {
        input = 'data/greetings.csv'
        batch = 'batch'
        character = 'turkey'
    }
    ```

Això estableix el directori de sortida a `custom-outdir-config/` més el valor del parametre `batch` com a subdirectori.
Ara podeu canviar la ubicació de sortida establint el parametre `--batch`:

```bash
nextflow run hello-config.nf --batch my_run
```

Això publica les sortides a `custom-outdir-config/my_run/`.

!!! note "Nota"

    L'opció CLI `-output-dir` té precedència sobre la configuració `outputDir`.
    Si s'estableix, l'opció de configuració s'ignorarà completament.

#### 2.2.2. Subdirectoris amb noms de batch i procés

També podem establir declaracions de `path` de sortida de subdirectori dinàmicament, per sortida.

Per exemple, podem organitzar les nostres sortides per procés fent referència a `<procés>.name` a la declaració de camí de sortida:

=== "Després"

    ```groovy title="hello-config.nf" linenums="42" hl_lines="3 7 11 15 19"
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

    ```groovy title="hello-config.nf" linenums="42" hl_lines="3 7 11 15 19"
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

Podem anar més lluny i compondre camins de subdirectori més complexos.

A l'edició anterior vam esborrar la distinció entre `intermediates` versus sortides finals al nivell superior.
Tornem a posar-ho, i també posem els fitxers en un subdirectori `params.batch`.

!!! tip "Consell"

    Incloure `params.batch` al `path` del bloc output, en lloc del `outputDir` de configuració, significa que no se sobreescriurà amb `-output-dir` a la CLI.

Primer, actualitzeu el fitxer de configuració per eliminar `${params.batch}` de `outputDir` (ja que l'estem movent a les declaracions de camí):

=== "Després"

    ```groovy title="nextflow.config" linenums="12" hl_lines="4"
    /*
    * Output settings
    */
    outputDir = "custom-outdir-config/"
    ```

=== "Abans"

    ```groovy title="nextflow.config" linenums="12" hl_lines="4"
    /*
    * Output settings
    */
    outputDir = "custom-outdir-config/${params.batch}"
    ```

Després, feu els següents canvis al fitxer de workflow:

=== "Després"

    ```groovy title="hello-config.nf" linenums="42" hl_lines="3 7 11 15 19"
    output {
        first_output {
            path { "${params.batch}/intermediates/${sayHello.name}" }
            mode 'copy'
        }
        uppercased {
            path { "${params.batch}/intermediates/${convertToUpper.name}" }
            mode 'copy'
        }
        collected {
            path { "${params.batch}/intermediates/${collectGreetings.name}" }
            mode 'copy'
        }
        batch_report {
            path { "${params.batch}/${collectGreetings.name}" }
            mode 'copy'
        }
        cowpy_art {
            path { "${params.batch}/${cowpy.name}" }
            mode 'copy'
        }
    }
    ```

=== "Abans"

    ```groovy title="hello-config.nf" linenums="42" hl_lines="3 7 11 15 19"
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

#### 2.2.3. Executar el pipeline

Vegem com funciona a la pràctica, establint tant `-output-dir` (o `-o` per abreujar) a `custom-outdir-config-2` com el nom del batch a `rep2` des de la línia de comandes:

```bash
nextflow run hello-config.nf -output-dir custom-outdir-config-2 --batch rep2
```

??? success "Sortida de la comanda"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-config.nf` [mad_curry] DSL2 - revision: 668a98ccb9

    executor >  local (8)
    [9e/6095e0] sayHello (1)       [100%] 3 of 3 ✔
    [05/454d52] convertToUpper (3) [100%] 3 of 3 ✔
    [ed/e3ddfb] collectGreetings   [100%] 1 of 1 ✔
    [39/5e063a] cowpy              [100%] 1 of 1 ✔
    ```

Això publica les sortides a `custom-outdir-config-2/rep2/`, amb el camí base especificat _i_ el subdirectori del nom del batch _i_ resultats agrupats per procés:

??? abstract "Contingut del directori"

    ```console
    custom-outdir-config-2
    └── rep2
        ├── collectGreetings
        │   └── rep2-report.txt
        ├── cowpy
        │   └── cowpy-COLLECTED-rep2-output.txt
        └── intermediates
            ├── collectGreetings
            │   └── COLLECTED-rep2-output.txt
            ├── convertToUpper
            │   ├── UPPER-Bonjour-output.txt
            │   ├── UPPER-Hello-output.txt
            │   └── UPPER-Holà-output.txt
            └── sayHello
                ├── Bonjour-output.txt
                ├── Hello-output.txt
                └── Holà-output.txt
    ```

### 2.3. Establir el mode de publicació a nivell de workflow

Finalment, amb l'esperit de reduir la quantitat de codi repetitiu, podem reemplaçar les declaracions `mode` per sortida amb una sola línia a la configuració.

#### 2.3.1. Afegir `workflow.output.mode` al fitxer de configuració

Afegiu el següent codi al fitxer `nextflow.config`:

=== "Després"

    ```groovy title="nextflow.config" linenums="12" hl_lines="5"
    /*
    * Output settings
    */
    outputDir = "custom-outdir-config/"
    workflow.output.mode = 'copy'
    ```

=== "Abans"

    ```groovy title="nextflow.config" linenums="12"
    /*
    * Output settings
    */
    outputDir = "custom-outdir-config/"
    ```

Establir `workflow.output.mode` al fitxer de configuració és suficient per sobreescriure el que s'estableix al fitxer de workflow, però eliminem el codi innecessari igualment.

#### 2.3.2. Eliminar el mode de sortida del fitxer de workflow

Feu els següents canvis al fitxer de workflow:

=== "Després"

    ```groovy title="hello-config.nf" linenums="42"
    output {
        first_output {
            path { "${params.batch}/intermediates/${sayHello.name}" }
        }
        uppercased {
            path { "${params.batch}/intermediates/${convertToUpper.name}" }
        }
        collected {
            path { "${params.batch}/intermediates/${collectGreetings.name}" }
        }
        batch_report {
            path { "${params.batch}/${collectGreetings.name}" }
        }
        cowpy_art {
            path { "${params.batch}/${cowpy.name}" }
        }
    }
    ```

=== "Abans"

    ```groovy title="hello-config.nf" linenums="42" hl_lines="4 8 12 16 20"
    output {
        first_output {
            path { "${params.batch}/intermediates/${sayHello.name}" }
            mode 'copy'
        }
        uppercased {
            path { "${params.batch}/intermediates/${convertToUpper.name}" }
            mode 'copy'
        }
        collected {
            path { "${params.batch}/intermediates/${collectGreetings.name}" }
            mode 'copy'
        }
        batch_report {
            path { "${params.batch}/${collectGreetings.name}" }
            mode 'copy'
        }
        cowpy_art {
            path { "${params.batch}/${cowpy.name}" }
            mode 'copy'
        }
    }
    ```

Això és més concís, oi?

#### 2.3.3. Executar el pipeline

Provem que funciona correctament:

```bash
nextflow run hello-config.nf -output-dir config-output-mode
```

??? success "Sortida de la comanda"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-config.nf` [small_stone] DSL2 - revision: 024d6361b5

    executor >  local (8)
    [e8/a0e93e] sayHello (1)       [100%] 3 of 3 ✔
    [14/176c9d] convertToUpper (3) [100%] 3 of 3 ✔
    [23/d667ca] collectGreetings   [100%] 1 of 1 ✔
    [e6/1dc80e] cowpy              [100%] 1 of 1 ✔
    ```

Això publica les sortides a `config-output-mode/`, i encara són totes còpies adequades, no enllaços simbòlics.

??? abstract "Contingut del directori"

    ```console
    config-output-mode
    └── batch
        ├── collectGreetings
        │   └── batch-report.txt
        ├── cowpy
        │   └── cowpy-COLLECTED-batch-output.txt
        └── intermediates
            ├── collectGreetings
            │   └── COLLECTED-batch-output.txt
            ├── convertToUpper
            │   ├── UPPER-Bonjour-output.txt
            │   ├── UPPER-Hello-output.txt
            │   └── UPPER-Holà-output.txt
            └── sayHello
                ├── Bonjour-output.txt
                ├── Hello-output.txt
                └── Holà-output.txt
    ```

La raó principal per la qual encara podríeu voler utilitzar la manera per sortida d'establir el mode és si voleu barrejar i combinar dins del mateix workflow, _és a dir_ tenir algunes sortides copiades i algunes enllaçades simbòlicament.

Hi ha moltes altres opcions que podeu personalitzar d'aquesta manera, però esperem que això us doni una idea de la gamma d'opcions i com utilitzar-les eficaçment per adaptar-se a les vostres preferències.

### Conclusió

Sabeu com controlar la denominació i estructura dels directoris on es publiquen les vostres sortides, així com el mode de publicació de sortida del workflow.

### Què segueix?

Apreneu a adaptar la configuració del vostre workflow al vostre entorn de càlcul, començant amb la tecnologia d'empaquetament de programari.

---

## 3. Seleccionar una tecnologia d'empaquetament de programari

Fins ara hem estat mirant elements de configuració que controlen com entren les entrades i on surten les sortides. Ara és hora de centrar-nos més específicament en adaptar la configuració del vostre workflow al vostre entorn de càlcul.

El primer pas en aquest camí és especificar d'on vindran els paquets de programari que s'executaran en cada pas.
Ja estan instal·lats a l'entorn de càlcul local?
Necessitem recuperar imatges i executar-les via un sistema de contenidors?
O necessitem recuperar paquets Conda i construir un entorn Conda local?

A la primera part d'aquest curs de formació (Parts 1-4) simplement vam utilitzar programari instal·lat localment al nostre workflow.
Després a la Part 5, vam introduir contenidors Docker i el fitxer `nextflow.config`, que vam utilitzar per habilitar l'ús de contenidors Docker.

Ara vegem com podem configurar una opció alternativa d'empaquetament de programari via el fitxer `nextflow.config`.

### 3.1. Deshabilitar Docker i habilitar Conda al fitxer de configuració

Simulem que estem treballant en un clúster HPC i l'administrador no permet l'ús de Docker per raons de seguretat.
Afortunadament per a nosaltres, Nextflow suporta múltiples altres tecnologies de contenidors com Singularity (que s'utilitza més àmpliament en HPC), i gestors de paquets de programari com Conda.

Podem canviar el nostre fitxer de configuració per utilitzar [Conda](https://nextflow.io/docs/latest/conda.html) en lloc de Docker.
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
nextflow run hello-config.nf --batch conda
```

??? success "Sortida de la comanda"

    ```console title="Output"
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-config.nf` [friendly_lamport] DSL2 - revision: 024d6361b5

    executor >  local (8)
    [e8/91c116] sayHello (2)       [100%] 3 of 3 ✔
    [fe/6a70ce] convertToUpper (3) [100%] 3 of 3 ✔
    [99/7cc493] collectGreetings   [100%] 1 of 1 ✔
    [3c/09fb59] cowpy              [100%] 1 of 1 ✔
    ```

Això hauria de funcionar sense problemes i produir les mateixes sortides que abans sota `custom-outdir-config/conda`.

Entre bastidors, Nextflow ha recuperat els paquets Conda i ha creat l'entorn, cosa que normalment requereix una mica de feina; així que és agradable que no hàgim de fer res d'això nosaltres mateixos!

!!! note "Nota"

    Això s'executa ràpidament perquè el paquet `cowpy` és força petit, però si esteu treballant amb paquets grans, pot trigar una mica més del normal la primera vegada, i podríeu veure la sortida de la consola quedar-se 'bloquejada' durant un minut o així abans de completar-se.
    Això és normal i es deu a la feina extra que Nextflow fa la primera vegada que utilitzeu un paquet nou.

Des del nostre punt de vista, sembla que funciona exactament igual que executar amb Docker, tot i que al backend la mecànica és una mica diferent.

Això significa que estem preparats per executar amb entorns Conda si cal.

??? info "Barrejar i combinar Docker i Conda"

    Com que aquestes directives s'assignen per procés, és possible 'barrejar i combinar', _és a dir_ configurar alguns dels processos del vostre workflow per executar amb Docker i altres amb Conda, per exemple, si la infraestructura de càlcul que esteu utilitzant suporta ambdós.
    En aquest cas, habilitaríeu tant Docker com Conda al vostre fitxer de configuració.
    Si ambdós estan disponibles per a un procés donat, Nextflow prioritzarà els contenidors.

    I com s'ha notat anteriorment, Nextflow suporta múltiples altres tecnologies d'empaquetament de programari i contenidors, així que no esteu limitats només a aquests dos.

### Conclusió

Sabeu com configurar quin paquet de programari hauria d'utilitzar cada procés, i com canviar entre tecnologies.

### Què segueix?

Apreneu a canviar la plataforma d'execució utilitzada per Nextflow per fer realment la feina.

---

## 4. Seleccionar una plataforma d'execució

Fins ara, hem estat executant el nostre pipeline amb l'executor local.
Això executa cada tasca a la màquina on s'està executant Nextflow.
Quan Nextflow comença, mira les CPUs i memòria disponibles.
Si els recursos de les tasques preparades per executar-se excedeixen els recursos disponibles, Nextflow retindrà les últimes tasques de l'execució fins que una o més de les tasques anteriors hagin acabat, alliberant els recursos necessaris.

L'executor local és convenient i eficient, però està limitat a aquesta única màquina. Per a càrregues de treball molt grans, podeu descobrir que la vostra màquina local és un coll d'ampolla, ja sigui perquè teniu una sola tasca que requereix més recursos dels que teniu disponibles, o perquè teniu tantes tasques que esperar que una sola màquina les executi trigaria massa.

Nextflow suporta [molts executors diferents](https://nextflow.io/docs/latest/executor.html), incloent planificadors HPC (Slurm, LSF, SGE, PBS, Moab, OAR, Bridge, HTCondor i altres) així com backends d'execució al núvol (AWS Batch, Google Cloud Batch, Azure Batch, Kubernetes i més).

### 4.1. Apuntar a un backend diferent

L'elecció de l'executor s'estableix per una directiva de procés anomenada `executor`.
Per defecte està establerta a `local`, així que la següent configuració està implícita:

```groovy title="Built-in configuration"
process {
    executor = 'local'
}
```

Per establir l'executor per apuntar a un backend diferent, simplement especificaríeu l'executor que voleu utilitzant una sintaxi similar a la descrita anteriorment per a assignacions de recursos (vegeu la [documentació d'executors](https://nextflow.io/docs/latest/executor.html) per a totes les opcions).

```groovy title="nextflow.config"
process {
    executor = 'slurm'
}
```

!!! warning "Advertència"

    No podem provar això realment a l'entorn de formació perquè no està configurat per connectar-se a un HPC.

### 4.2. Tractar amb sintaxi específica del backend per a paràmetres d'execució

La majoria de plataformes de computació d'alt rendiment permeten (i de vegades requereixen) que especifiqueu certs paràmetres com ara sol·licituds i limitacions d'assignació de recursos (per exemple, nombre de CPUs i memòria) i nom de la cua de tasques a utilitzar.

Malauradament, cadascun d'aquests sistemes utilitza diferents tecnologies, sintaxis i configuracions per definir com s'ha de definir i enviar una tasca al planificador rellevant.

??? abstract "Exemples"

    Per exemple, la mateixa tasca que requereix 8 CPUs i 4GB de RAM per ser executada a la cua "my-science-work" necessita ser expressada de les següents maneres diferents segons el backend.

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
Proporciona una sintaxi estandarditzada perquè pugueu especificar les propietats rellevants com ara [`cpus`](https://nextflow.io/docs/latest/reference/process.html#cpus), [`memory`](https://nextflow.io/docs/latest/reference/process.html#memory) i [`queue`](https://nextflow.io/docs/latest/reference/process.html#queue) (vegeu [directives de procés](https://nextflow.io/docs/latest/reference/process.html#process-directives) per a altres propietats) només una vegada.
Després, en temps d'execució, Nextflow utilitzarà aquestes configuracions per generar els scripts específics del backend apropiats basats en la configuració de l'executor.

Cobrirem aquesta sintaxi estandarditzada a la següent secció.

### Conclusió

Ara sabeu com canviar l'executor per utilitzar diferents tipus d'infraestructura de càlcul.

### Què segueix?

Apreneu a avaluar i expressar assignacions i limitacions de recursos a Nextflow.

---

## 5. Controlar les assignacions de recursos de càlcul

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

Si no sabeu per endavant quanta CPU i memòria és probable que necessitin els vostres processos, podeu fer una mica de perfilat de recursos, és a dir, executar el workflow amb algunes assignacions per defecte, registrar quant va utilitzar cada procés, i a partir d'aquí, estimar com ajustar les assignacions base.

Convenientment, Nextflow inclou eines integrades per fer això, i generarà feliçment un informe per a vosaltres a petició.

Per fer-ho, afegiu `-with-report <nom_fitxer>.html` a la vostra línia de comandes.

```bash
nextflow run hello-config.nf -with-report report-config-1.html
```

L'informe és un fitxer html, que podeu descarregar i obrir al vostre navegador. També podeu fer clic dret sobre ell a l'explorador de fitxers de l'esquerra i fer clic a `Show preview` per veure'l a l'entorn de formació.

Preneu-vos uns minuts per mirar l'informe i veure si podeu identificar algunes oportunitats per ajustar recursos.
Assegureu-vos de fer clic a les pestanyes que mostren els resultats d'utilització com a percentatge del que es va assignar.

Vegeu [Informes](https://nextflow.io/docs/latest/reports.html) per a documentació sobre totes les funcionalitats disponibles.

### 5.2. Establir assignacions de recursos per a tots els processos

El perfilat mostra que els processos del nostre workflow de formació són molt lleugers, així que reduïm l'assignació de memòria per defecte a 1GB per procés.

Afegiu el següent al vostre fitxer `nextflow.config`, abans de la secció de paràmetres del pipeline:

=== "Després"

    ```groovy title="nextflow.config" linenums="1" hl_lines="4-9"
    docker.enabled = false
    conda.enabled = true

    /*
    * Process settings
    */
    process {
        memory = 1.GB
    }

    /*
    * Pipeline parameters
    */
    params {
        input = 'data/greetings.csv'
        batch = 'batch'
        character = 'turkey'
    }
    ```

=== "Abans"

    ```groovy title="nextflow.config" linenums="1"
    docker.enabled = false
    conda.enabled = true

    /*
    * Pipeline parameters
    */
    params {
        input = 'data/greetings.csv'
        batch = 'batch'
        character = 'turkey'
    }
    ```

Això ajudarà a reduir la quantitat de càlcul que consumim.

### 5.3. Establir assignacions de recursos per a un procés específic

Al mateix temps, simularem que el procés `cowpy` requereix més recursos que els altres, només perquè puguem demostrar com ajustar assignacions per a un procés individual.

=== "Després"

    ```groovy title="nextflow.config" linenums="4" hl_lines="6-9"
    /*
    * Process settings
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
    * Process settings
    */
    process {
        memory = 1.GB
    }
    ```

Amb aquesta configuració, tots els processos sol·licitaran 1GB de memòria i una sola CPU (el valor per defecte implícit), excepte el procés `cowpy`, que sol·licitarà 2GB i 2 CPUs.

!!! tip "Consell"

    Si teniu una màquina amb poques CPUs i assigneu un nombre alt per procés, podríeu veure crides de procés posant-se en cua una darrere l'altra.
    Això és perquè Nextflow assegura que no sol·licitem més CPUs de les que estan disponibles.

### 5.4. Executar el workflow amb la configuració actualitzada

Provem això, proporcionant un nom de fitxer diferent per a l'informe de perfilat perquè puguem comparar el rendiment abans i després dels canvis de configuració.

```bash
nextflow run hello-config.nf -with-report report-config-2.html
```

Probablement no notareu cap diferència real ja que aquesta és una càrrega de treball tan petita, però aquest és l'enfocament que utilitzaríeu per analitzar el rendiment i els requisits de recursos d'un workflow del món real.

És molt útil quan els vostres processos tenen diferents requisits de recursos. Us permet dimensionar correctament les assignacions de recursos que configureu per a cada procés basant-vos en dades reals, no conjectures.

!!! tip "Consell"

    Això és només una petita mostra del que podeu fer per optimitzar l'ús de recursos.
    Nextflow mateix té una [lògica de reintent dinàmic](https://nextflow.io/docs/latest/process.html#dynamic-task-resources) realment interessant integrada per reintentar tasques que fallen a causa de limitacions de recursos.
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
Tanmateix, si intentéssiu executar el workflow amb assignacions de recursos que excedeixen aquests límits, després miréssiu la comanda `sbatch` al fitxer d'script `.command.run`, veuríeu que les sol·licituds que realment s'envien a l'executor estan limitades als valors especificats per `resourceLimits`.

??? info "Configuracions de referència institucionals"

    El projecte nf-core ha compilat una [col·lecció de fitxers de configuració](https://nf-co.re/configs/) compartits per diverses institucions arreu del món, cobrint una àmplia gamma d'executors HPC i núvol.

    Aquestes configuracions compartides són valuoses tant per a persones que hi treballen i per tant poden utilitzar la configuració de la seva institució directament, com a model per a persones que busquen desenvolupar una configuració per a la seva pròpia infraestructura.

### Conclusió

Sabeu com generar un informe de perfilat per avaluar la utilització de recursos i com modificar assignacions de recursos per a tots els processos i/o per a processos individuals, així com establir limitacions de recursos per executar en HPC.

### Què segueix?

Apreneu a configurar perfils de configuració predefinits i canviar entre ells en temps d'execució.

---

## 6. Utilitzar perfils per canviar entre configuracions predefinides

Us hem mostrat diverses maneres en què podeu personalitzar la configuració del vostre pipeline depenent del projecte en què estigueu treballant o de l'entorn de càlcul que estigueu utilitzant.

Potser voleu canviar entre configuracions alternatives depenent de quina infraestructura de càlcul estigueu utilitzant. Per exemple, potser voleu desenvolupar i executar proves a petita escala localment al vostre portàtil, després executar càrregues de treball a gran escala en HPC o núvol.

Nextflow us permet configurar qualsevol nombre de [perfils](https://nextflow.io/docs/latest/config.html#config-profiles) que descriguin diferents configuracions, que després podeu seleccionar en temps d'execució utilitzant un argument de línia de comandes, en lloc d'haver de modificar el fitxer de configuració mateix.

### 6.1. Crear perfils per canviar entre desenvolupament local i execució en HPC

Configurem dos perfils alternatius; un per executar càrregues a petita escala en un ordinador normal, on utilitzarem contenidors Docker, i un per executar en un HPC universitari amb un planificador Slurm, on utilitzarem paquets Conda.

#### 6.1.1. Configurar els perfils

Afegiu el següent al vostre fitxer `nextflow.config`, després de la secció de paràmetres del pipeline però abans de la configuració de sortida:

=== "Després"

    ```groovy title="nextflow.config" linenums="15" hl_lines="10-27"
    /*
    * Pipeline parameters
    */
    params {
        input = 'data/greetings.csv'
        batch = 'batch'
        character = 'turkey'
    }

    /*
    * Profiles
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

    /*
    * Output settings
    */
    outputDir = "custom-outdir-config/"
    workflow.output.mode = 'copy'
    ```

=== "Abans"

    ```groovy title="nextflow.config" linenums="15"
    /*
    * Pipeline parameters
    */
    params {
        input = 'data/greetings.csv'
        batch = 'batch'
        character = 'turkey'
    }

    /*
    * Output settings
    */
    outputDir = "custom-outdir-config/"
    workflow.output.mode = 'copy'
    ```

Veieu que per a l'HPC universitari, també estem especificant limitacions de recursos.

#### 6.1.2. Executar el workflow amb un perfil

Per especificar un perfil a la nostra línia de comandes de Nextflow, utilitzem l'argument `-profile`.

Provem d'executar el workflow amb la configuració `my_laptop`.

```bash
nextflow run hello-config.nf -profile my_laptop
```

??? success "Sortida de la comanda"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-config.nf` [hungry_sanger] DSL2 - revision: 024d6361b5

    executor >  local (8)
    [b0/fb2ec9] sayHello (3)       [100%] 3 of 3 ✔
    [4a/e039f0] convertToUpper (3) [100%] 3 of 3 ✔
    [6f/408fa9] collectGreetings   [100%] 1 of 1 ✔
    [f1/fd6520] cowpy              [100%] 1 of 1 ✔
    ```

Com podeu veure, això ens permet alternar entre configuracions molt convenientment en temps d'execució.

!!! warning "Advertència"

    El perfil `univ_hpc` no s'executarà correctament a l'entorn de formació ja que no tenim accés a un planificador Slurm.

Si en el futur trobem altres elements de configuració que sempre co-ocorren amb aquests, podem simplement afegir-los al(s) perfil(s) corresponent(s).
També podem crear perfils addicionals si hi ha altres elements de configuració que volem agrupar junts.

### 6.2. Crear un perfil de paràmetres de prova

Els perfils no són només per a configuració d'infraestructura.
També podem utilitzar-los per establir valors per defecte per a paràmetres de workflow, per facilitar que altres provin el workflow sense haver de reunir valors d'entrada apropiats ells mateixos.
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

Si afegim un perfil de prova per al nostre workflow, el bloc `profiles` esdevé:

```groovy title="nextflow.config" linenums="24" hl_lines="18-22"
/*
* Profiles
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

Igual que per a perfils de configuració tècnica, podeu configurar múltiples perfils diferents especificant paràmetres sota qualsevol nom arbitrari que us agradi.

#### 6.2.2. Executar el workflow localment amb el perfil de prova

Convenientment, els perfils no són mútuament excloents, així que podem especificar múltiples perfils a la nostra línia de comandes utilitzant la següent sintaxi `-profile <perfil1>,<perfil2>` (per a qualsevol nombre de perfils).

Si combineu perfils que estableixen valors per als mateixos elements de configuració i es descriuen al mateix fitxer de configuració, Nextflow resoldrà el conflicte utilitzant el valor que hagi llegit en últim lloc (_és a dir_ el que vingui més tard al fitxer).
Si les configuracions conflictives s'estableixen en diferents fonts de configuració, s'aplica l'[ordre de precedència](https://nextflow.io/docs/latest/config.html) per defecte.

Provem d'afegir el perfil de prova a la nostra comanda anterior:

```bash
nextflow run hello-config.nf -profile my_laptop,test
```

??? success "Sortida de la comanda"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-config.nf` [modest_becquerel] DSL2 - revision: 024d6361b5

    executor >  local (8)
    [4c/fe2580] sayHello (1)       [100%] 3 of 3 ✔
    [fd/7d9017] convertToUpper (3) [100%] 3 of 3 ✔
    [13/1523bd] collectGreetings   [100%] 1 of 1 ✔
    [06/a1ee14] cowpy              [100%] 1 of 1 ✔
    ```

Això utilitzarà Docker on sigui possible i produirà sortides sota `custom-outdir-config/test`, i aquesta vegada el caràcter és el duo còmic `dragonandcow`.

??? abstract "Contingut del fitxer"

    ```console title="custom-outdir-config/test/cowpy/cowpy-COLLECTED-test-output.txt"
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

    Podem apuntar a URLs per a fitxers més grans que s'emmagatzemen externament.
    Nextflow els descarregarà automàticament sempre que hi hagi una connexió oberta.

    Per a més detalls, vegeu la Missió Secundària [Treballar amb Fitxers](../side_quests/working_with_files.md)

### 6.3. Utilitzar `nextflow config` per veure la configuració resolta

Com s'ha notat anteriorment, de vegades el mateix parametre pot estar establert a diferents valors en perfils que voleu combinar.
I més generalment, hi ha nombrosos llocs on es poden emmagatzemar elements de configuració, i de vegades les mateixes propietats poden estar establertes a diferents valors en diferents llocs.

Nextflow aplica un [ordre de precedència](https://nextflow.io/docs/latest/config.html) establert per resoldre qualsevol conflicte, però això pot ser complicat de determinar vosaltres mateixos.
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

    outputDir = 'custom-outdir-config/'

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

    outputDir = 'custom-outdir-config/'

    workflow {
      output {
          mode = 'copy'
      }
    }
    ```

Això esdevé especialment útil per a projectes complexos que impliquen múltiples capes de configuració.

### Conclusió

Sabeu com utilitzar perfils per seleccionar una configuració predefinida en temps d'execució amb un mínim d'esforç.
Més generalment, sabeu com configurar les vostres execucions de workflow per adaptar-se a diferents plataformes de càlcul i millorar la reproduïbilitat de les vostres anàlisis.

### Què segueix?

Celebreu i doneu-vos una gran palmada a l'esquena! Heu completat el vostre primer curs de desenvolupador de Nextflow.

Aneu al [resum final del curs](./next_steps.md) per revisar el que heu après i descobrir què ve després.

---

## Qüestionari

<quiz>
Quin és el nom del fitxer de configuració que Nextflow carrega automàticament?
- [ ] `config.nf`
- [ ] `pipeline.config`
- [x] `nextflow.config`
- [ ] `workflow.config`
</quiz>

<quiz>
Què té precedència quan el mateix parametre s'estableix tant al fitxer de configuració com a la línia de comandes?
- [ ] El valor del fitxer de configuració
- [x] El valor de la línia de comandes
- [ ] El primer valor trobat
- [ ] Cap; causa un error

Més informació: [1.1. Moure els valors per defecte a `nextflow.config`](#11-move-default-values-to-nextflowconfig)
</quiz>

<quiz>
Podeu tenir tant Docker com Conda habilitats a la mateixa configuració?
- [x] Sí, Nextflow pot utilitzar ambdós depenent de les directives de procés
- [ ] No, només un pot estar habilitat alhora
- [ ] Sí, però només en perfils
- [ ] No, són mútuament excloents
</quiz>

<quiz>
Si tant Docker com Conda estan habilitats i un procés té ambdues directives, quina es prioritza?
- [x] Docker (contenidors)
- [ ] Conda
- [ ] La primera definida
- [ ] Causa un error

Més informació: [3. Seleccionar una tecnologia d'empaquetament de programari](#3-select-a-software-packaging-technology)
</quiz>

<quiz>
Quina és l'assignació de memòria per defecte per als processos de Nextflow?
- [ ] 1 GB
- [x] 2 GB
- [ ] 4 GB
- [ ] Sense límit
</quiz>

<quiz>
Com establiu requisits de recursos per a un procés específic al fitxer de configuració?
- [ ] `#!groovy processName.memory = '4 GB'`
- [ ] `#!groovy process.memory.processName = '4 GB'`
- [x] `#!groovy process { withName: 'processName' { memory = '4 GB' } }`
- [ ] `#!groovy resources.processName.memory = '4 GB'`

Més informació: [5.3. Establir assignacions de recursos per a un procés específic](#53-set-resource-allocations-for-a-specific-process)
</quiz>

<quiz>
Quina opció de línia de comandes genera un informe d'utilització de recursos?
- [ ] `-with-metrics`
- [ ] `-with-stats`
- [x] `-with-report`
- [ ] `-with-profile`

Més informació: [5.1. Executar el workflow per generar un informe d'utilització de recursos](#51-run-the-workflow-to-generate-a-resource-utilization-report)
</quiz>

<quiz>
Què fa la directiva `resourceLimits`?
- [ ] Estableix requisits mínims de recursos
- [ ] Assigna recursos als processos
- [x] Limita els recursos màxims que es poden sol·licitar
- [ ] Monitoritza l'ús de recursos

Més informació: [5.5. Afegir límits de recursos](#55-add-resource-limits)
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
Com especifiqueu un fitxer de paràmetres quan executeu Nextflow?
- [ ] `--params params.json`
- [ ] `-config params.json`
- [x] `-params-file params.json`
- [ ] `--input params.json`

Més informació: [1.3. Utilitzar un fitxer de paràmetres](#13-use-a-parameter-file)
</quiz>

<quiz>
Per a què es poden utilitzar els perfils? (Seleccioneu totes les que apliquin)
- [x] Definir configuracions específiques d'infraestructura
- [x] Establir límits de recursos per a diferents entorns
- [x] Proporcionar paràmetres de prova
- [ ] Definir nous processos

Més informació: [6. Utilitzar perfils per canviar entre configuracions predefinides](#6-use-profiles-to-switch-between-preset-configurations)
</quiz>

<quiz>
Com especifiqueu múltiples perfils en una sola comanda?
- [ ] `-profile perfil1 -profile perfil2`
- [ ] `-profiles perfil1,perfil2`
- [x] `-profile perfil1,perfil2`
- [ ] `--profile perfil1 --profile perfil2`

Més informació: [6. Utilitzar perfils per canviar entre configuracions predefinides](#6-use-profiles-to-switch-between-preset-configurations)
</quiz>
