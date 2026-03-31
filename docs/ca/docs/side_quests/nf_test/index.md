# Proves amb nf-test

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducció assistida per IA - [més informació i suggeriments](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Poder comprovar sistemàticament que cada part del vostre workflow fa el que se suposa que ha de fer és fonamental per a la reproductibilitat i el manteniment a llarg termini, i pot ser d'una gran ajuda durant el procés de desenvolupament.

Dediquem un moment a parlar de per què les proves són tan importants. Si esteu desenvolupant un workflow, una de les primeres coses que fareu és agafar algunes dades de prova que sabeu que són vàlides i que haurien de produir un resultat. Afegiu el primer procés al pipeline i el connecteu a les vostres entrades per fer-lo funcionar. Després, per comprovar que tot funciona, l'executeu amb les dades de prova. Suposant que funciona, passeu al procés següent i torneu a executar les dades de prova. Repetiu aquest procés fins que teniu un pipeline amb el qual esteu satisfets.

Llavors, potser afegiu un paràmetre simple de vertader o fals com ara `--skip_process`. Ara heu d'executar el pipeline dues vegades, una amb cada paràmetre, per assegurar-vos que funciona com s'espera. Però espera, com comprovem si `--skip_process` realment omet el procés? Hem d'examinar les sortides o revisar els fitxers de registre! Això és una molèstia i és propens a errors.

A mesura que desenvolupeu el vostre pipeline, ràpidament es tornarà tan complex que provar manualment cada iteració serà lent i propens a errors. A més, si trobeu un error, serà molt difícil determinar exactament d'on prové dins del vostre pipeline. Aquí és on entren les proves.

Les proves us permeten comprovar sistemàticament que cada part del vostre pipeline funciona com s'espera. Els beneficis per a un desenvolupador de proves ben escrites són enormes:

- **Confiança**: Com que les proves cobreixen tot el pipeline, podeu estar segurs que canviar alguna cosa no afecta res més
- **Confiabilitat**: Quan diversos desenvolupadors treballen al pipeline, saben que els altres desenvolupadors no l'han trencat ni cap dels seus components.
- **Transparència**: Les proves mostren on falla un pipeline i faciliten la localització del problema. També funcionen com una forma de documentació, mostrant com executar un procés o workflow.
- **Velocitat**: Com que les proves estan automatitzades, es poden executar molt ràpidament i de manera repetida. Podeu iterar ràpidament amb menys por d'introduir nous errors.

Hi ha molts tipus diferents de proves que podem escriure:

1. **Proves a nivell de mòdul**: Per a processos individuals
2. **Proves a nivell de workflow**: Per a un workflow individual
3. **Proves a nivell de pipeline**: Per al pipeline en conjunt
4. **Proves de rendiment**: Per a la velocitat i l'eficiència del pipeline
5. **Proves d'estrès**: Per avaluar el rendiment del pipeline en condicions extremes per determinar els seus límits

Provar processos individuals és anàleg a les proves unitàries en altres llenguatges. Provar el workflow o tot el pipeline és anàleg al que s'anomenen proves d'integració en altres llenguatges, on provem les interaccions dels components.

[**nf-test**](https://www.nf-test.com/) és una eina que us permet escriure proves a nivell de mòdul, workflow i pipeline. En resum, us permet comprovar sistemàticament que cada part individual del pipeline funciona com s'espera, _de manera aïllada_.

### Objectius d'aprenentatge

En aquesta missió secundària, aprendreu a utilitzar nf-test per escriure una prova a nivell de workflow per al pipeline, així com proves a nivell de mòdul per als tres processos que invoca.

Al final d'aquesta missió secundària, podreu utilitzar les tècniques següents de manera efectiva:

- Inicialitzar nf-test al vostre projecte
- Generar proves a nivell de mòdul i de workflow
- Afegir tipus comuns d'assertions
- Entendre quan utilitzar snapshots en lloc d'assertions de contingut
- Executar proves per a tot un projecte

Aquestes habilitats us ajudaran a implementar una estratègia de proves completa als vostres projectes de pipeline, assegurant que siguin més robustos i fàcils de mantenir.

### Prerequisits

Abans d'abordar aquesta missió secundària, hauríeu de:

- Haver completat el tutorial [Hello Nextflow](../hello_nextflow/README.md) o un curs equivalent per a principiants.
- Estar còmodes amb els conceptes i mecanismes bàsics de Nextflow (processos, canals, operadors, treball amb fitxers, metadades)

---

## 0. Primers passos

#### Obriu el codespace de formació

Si encara no ho heu fet, assegureu-vos d'obrir l'entorn de formació tal com es descriu a la [Configuració de l'entorn](../envsetup/index.md).

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### Moveu-vos al directori del projecte

Movem-nos al directori on es troben els fitxers d'aquest tutorial.

```bash
cd side-quests/nf-test
```

Podeu configurar VSCode perquè es centri en aquest directori:

```bash
code .
```

#### Reviseu els materials

Trobareu un fitxer de workflow principal i un fitxer CSV anomenat `greetings.csv` que conté l'entrada del pipeline.

```console title="Directory contents"
.
├── greetings.csv
└── main.nf
```

Per a una descripció detallada dels fitxers, consulteu l'[escalfament de Hello Nextflow](../hello_nextflow/00_orientation.md).

El workflow que provarem és un subconjunt del Hello workflow construït a [Hello Workflow](../hello_nextflow/03_hello_workflow.md).

??? example "Què fa el workflow Hello Nextflow?"

    Si no heu fet la formació [Hello Nextflow](../hello_nextflow/index.md), aquí teniu una visió general ràpida del que fa aquest workflow senzill.

    El workflow pren un fitxer CSV que conté salutacions, executa quatre passos de transformació consecutius sobre elles i genera un únic fitxer de text que conté una imatge ASCII d'un personatge divertit dient les salutacions.

    Els quatre passos s'implementen com a processos Nextflow (`sayHello`, `convertToUpper`, `collectGreetings` i `cowpy`) emmagatzemats en fitxers de mòdul separats.

    1. **`sayHello`:** Escriu cada salutació al seu propi fitxer de sortida (p. ex., "Hello-output.txt")
    2. **`convertToUpper`:** Converteix cada salutació a majúscules (p. ex., "HELLO")
    3. **`collectGreetings`:** Recull totes les salutacions en majúscules en un únic fitxer per lots
    4. **`cowpy`:** Genera art ASCII utilitzant l'eina `cowpy`

    Els resultats es publiquen en un directori anomenat `results/`, i la sortida final del pipeline (quan s'executa amb els paràmetres per defecte) és un fitxer de text pla que conté art ASCII d'un personatge dient les salutacions en majúscules.

    En aquesta missió secundària, utilitzem una forma intermèdia del Hello workflow que només conté els dos primers processos. <!-- TODO: change this to use the full finished workflow as suggested in https://github.com/nextflow-io/training/issues/735 -->

El subconjunt amb el qual treballarem està compost de dos processos: `sayHello` i `convertToUpper`.
Podeu veure el codi complet del workflow a continuació.

??? example "Codi del workflow"

    ```groovy title="main.nf"
    /*
    * Paràmetres del pipeline
    */
    params.input_file = "greetings.csv"

    /*
    * Utilitza echo per imprimir 'Hello World!' a la sortida estàndard
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
    * Utilitza una utilitat de substitució de text per convertir la salutació a majúscules
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

        // crea un canal per a les entrades des d'un fitxer CSV
        greeting_ch = channel.fromPath(params.input_file).splitCsv().flatten()

        // emet una salutació
        sayHello(greeting_ch)

        // converteix la salutació a majúscules
        convertToUpper(sayHello.out)
    }
    ```

#### Executeu el workflow

Executem el workflow per assegurar-nos que funciona com s'espera.

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

FELICITATS! Acabeu d'executar una prova!

"Espera, com? Acabo d'executar el workflow i ha funcionat! Com és que això és una prova?"

Bona pregunta!

Analitzem el que acaba de passar.

Heu executat el workflow amb els paràmetres per defecte, heu confirmat que ha funcionat i esteu satisfets amb els resultats. Aquesta és l'essència de les proves. Si heu treballat amb el curs de formació Hello Nextflow, haureu notat que sempre començàvem cada secció executant el workflow que estàvem utilitzant com a punt de partida, per confirmar que tot estava configurat correctament.

Les proves de programari essencialment fan aquest procés per nosaltres.

#### Reviseu l'assignació

El vostre repte és afegir proves estandarditzades a aquest workflow utilitzant nf-test, per tal de facilitar la verificació que cada part continua funcionant com s'espera en cas que es facin més canvis.

<!-- TODO: give a bit more details, similar to how it's done in the Metadata side quest -->

#### Llista de comprovació de preparació

Creieu que esteu preparats per submergir-vos?

- [ ] Entenc l'objectiu d'aquest curs i els seus prerequisits
- [ ] El meu codespace està en funcionament
- [ ] He configurat el meu directori de treball adequadament
- [ ] He executat el workflow correctament
- [ ] Entenc l'assignació

Si podeu marcar totes les caselles, esteu a punt per començar.

---

## 1. Inicialitzeu `nf-test`

El paquet `nf-test` proporciona una comanda d'inicialització que configura algunes coses per tal que puguem començar a desenvolupar proves per al nostre projecte.

```bash
nf-test init
```

Això hauria de produir la sortida següent:

```bash
🚀 nf-test 0.9.3
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr

Project configured. Configuration is stored in nf-test.config
```

També crea un directori `tests` que conté un stub de fitxer de configuració.

### 1.1. Genereu un nf-test

`nf-test` inclou un conjunt d'eines per construir fitxers nf-test, estalviant-nos la major part de la feina. Aquestes es troben sota la subcomanda `generate`. Generem una prova per al pipeline:

```bash
nf-test generate pipeline main.nf
```

```console title="Output"
> nf-test generate pipeline main.nf

Load source file '/workspaces/training/side-quests/nf-test/main.nf'
Wrote pipeline test file '/workspaces/training/side-quests/nf-test/tests/main.nf.test

SUCCESS: Generated 1 test files.
```

Això crearà un fitxer `main.nf.test` dins del directori `tests`. Aquest és el nostre fitxer de proves a nivell de pipeline. Si executeu `tree tests/` hauríeu de veure alguna cosa com això:

```console title="Test directory contents"
tests/
├── main.nf.test
└── nextflow.config
```

El fitxer `main.nf.test` és el nostre fitxer de proves a nivell de pipeline. Obrim-lo i examinem el contingut.

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

Dediquem un moment a entendre l'estructura del fitxer de proves.

El bloc `nextflow_pipeline` és el punt d'entrada per a totes les proves a nivell de pipeline. Conté el següent:

- `name`: El nom de la prova.
- `script`: El camí al script del pipeline.

El bloc `test` és la prova en si. Conté el següent:

- `when`: Les condicions sota les quals s'ha d'executar la prova. Inclou els paràmetres que s'utilitzaran per executar el pipeline.
- `then`: Les assertions que s'han de fer. Inclou els resultats esperats del pipeline.

En llenguatge planer, la lògica de la prova es llegeix de la manera següent:
"**Quan** aquests _paràmetres_ es proporcionen a aquest _pipeline_, **llavors** esperem veure aquests resultats."

Aquesta no és una prova funcional; demostrarem com convertir-la en una a la secció següent.

### Una nota sobre els noms de les proves

A l'exemple anterior, hem utilitzat el nom per defecte "Should run without failures" que és adequat per a una prova bàsica que simplement comprova si el pipeline s'executa correctament. No obstant això, a mesura que afegim casos de prova més específics, hauríem d'utilitzar noms més descriptius que indiquin el que realment estem provant. Per exemple:

- "Should convert input to uppercase" - quan es prova una funcionalitat específica
- "Should handle empty input gracefully" - quan es proven casos límit
- "Should respect max memory parameter" - quan es proven restriccions de recursos
- "Should create expected output files" - quan es prova la generació de fitxers

Els bons noms de proves haurien de:

1. Començar amb "Should" per deixar clar quin és el comportament esperat
2. Descriure la funcionalitat o l'escenari específic que s'està provant
3. Ser prou clars perquè si la prova falla, sapigueu quina funcionalitat no funciona

A mesura que afegim més assertions i casos de prova específics més endavant, utilitzarem aquests noms més descriptius per deixar clar el que verifica cada prova.

### 1.2. Executeu la prova

Executem la prova per veure què passa.

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

La prova falla! Què ha passat?

1. nf-test ha intentat executar el pipeline tal com és, utilitzant la configuració del bloc `when`:

```groovy title="tests/main.nf.test"
when {
    params {
        // define parameters here. Example:
        // outdir = "tests/results"
    }
}
```

2. nf-test ha comprovat l'estat del pipeline i l'ha comparat amb el bloc `when`:

```groovy title="tests/main.nf.test"
then {
    assert workflow.success
}
```

Observeu com nf-test ha informat que el pipeline ha fallat i ha proporcionat el missatge d'error de Nextflow:

```console title="Error"
ERROR ~ No such file or directory: /workspaces/training/side-quests/nf-test/.nf-test/tests/693ba951a20fec36a5a9292ed1cc8a9f/greetings.csv
```

Quin era el problema? Recordeu que el pipeline té un fitxer greetings.csv al directori del projecte. Quan nf-test executa el pipeline, buscarà aquest fitxer, però no el pot trobar. El fitxer hi és, però què passa? Bé, si mirem el camí podem veure que la prova s'està executant al camí `./nf-test/tests/longHashString/`. Igual que Nextflow, nf-test crea un nou directori per a cada prova per mantenir tot aïllat. El fitxer de dades no es troba allà, de manera que hem de corregir el camí al fitxer a la prova original.

Tornem al fitxer de proves i canviem el camí al fitxer al bloc `when`.

Potser us pregunteu com apuntarem a l'arrel del pipeline a la prova. Com que aquesta és una situació habitual, nf-test té un conjunt de variables globals que podem utilitzar per facilitar-nos la vida. Podeu trobar la llista completa [aquí](https://www.nf-test.com/docs/testcases/global_variables/) però de moment utilitzarem la variable `projectDir`, que fa referència a l'arrel del projecte del pipeline.

=== "Després"

    ```groovy title="tests/main.nf.test" linenums="1" hl_lines="3"
    when {
        params {
            input_file = "${projectDir}/greetings.csv"
        }
    }
    ```

=== "Abans"

    ```groovy title="tests/main.nf.test" linenums="1" hl_lines="3 4"
    when {
        params {
            // define parameters here. Example:
            // outdir = "tests/results"
        }
    }
    ```

Executem la prova de nou per veure si funciona.

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

Èxit! El pipeline s'executa correctament i la prova passa. Executeu-lo tantes vegades com vulgueu i sempre obtindreu el mateix resultat!

Per defecte, la sortida de Nextflow està oculta, però per convèncer-vos que nf-test realment està executant el workflow, podeu utilitzar l'indicador `--verbose`:

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

### 1.3. Afegiu assertions

Una comprovació senzilla és assegurar-nos que el nostre pipeline executa tots els processos que esperem i no n'omet cap en silenci. Recordeu que el nostre pipeline executa 6 processos, un anomenat `sayHello` i un anomenat `convertToUpper` per a cadascuna de les 3 salutacions.

Afegim una assertion a la nostra prova per comprovar que el pipeline executa el nombre esperat de processos. També actualitzarem el nom de la prova per reflectir millor el que estem provant.

=== "Després"

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

=== "Abans"

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

El nom de la prova ara reflecteix millor el que realment estem verificant: no només que el pipeline s'executa sense fallar, sinó que executa el nombre esperat de processos.

Executem la prova de nou per veure si funciona.

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

Èxit! El pipeline s'executa correctament i la prova passa. Ara hem començat a provar els detalls del pipeline, a més de l'estat general.

### 1.4. Proveu la sortida

Afegim una assertion a la nostra prova per comprovar que el fitxer de sortida s'ha creat. L'afegirem com una prova separada, amb un nom informatiu, per facilitar la interpretació dels resultats.

=== "Després"

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

=== "Abans"

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

Executeu la prova de nou per veure si funciona.

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

Èxit! Les proves passen perquè el pipeline s'ha completat correctament, el nombre correcte de processos s'ha executat i els fitxers de sortida s'han creat. Això també us hauria de mostrar com d'útil és proporcionar noms informatius per a les vostres proves.

Això és només la superfície; podem continuar escrivint assertions per comprovar els detalls del pipeline, però de moment passem a provar els internals del pipeline.

### Conclusió

Sabeu com escriure un nf-test per a un pipeline.

### Què segueix?

Apreneu com provar un procés Nextflow.

---

## 2. Proveu un procés Nextflow

No hem d'escriure proves per a cada part del pipeline, però com més proves tenim, més exhaustius podem ser sobre el pipeline i més segurs podem estar que funciona com s'espera. En aquesta secció provarem els dos processos del pipeline com a unitats individuals.

### 2.1. Proveu el procés `sayHello`

Comencem amb el procés `sayHello`.

Utilitzem de nou la comanda `nf-test generate` per generar proves per al procés.

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

Centrem-nos ara en el procés `sayhello` al fitxer `main.sayhello.nf.test`.

Obrim el fitxer i examinem el contingut.

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

Com abans, comencem amb els detalls de la prova, seguits dels blocs `when` i `then`. No obstant això, també tenim un bloc `process` addicional que ens permet definir les entrades del procés.

Executem la prova per veure si funciona.

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

La prova falla perquè el procés `sayHello` declara 1 entrada però s'ha cridat amb 0 arguments. Corregim-ho afegint una entrada al procés. Recordeu de [Hello Workflow](../hello_nextflow/03_hello_workflow.md) (i de la secció d'escalfament anterior) que el nostre procés `sayHello` pren una única entrada de valor, que haurem de proporcionar. També hauríem de corregir el nom de la prova per reflectir millor el que estem provant.

=== "Després"

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

=== "Abans"

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

Executem la prova de nou per veure si funciona.

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

Èxit! La prova passa perquè el procés `sayHello` s'ha executat correctament i la sortida s'ha creat.

### 2.2. Examineu el snapshot creat per la prova

Si mirem el fitxer `tests/main.sayhello.nf.test`, podem veure que utilitza un mètode `snapshot()` al bloc d'assertions:

```groovy title="tests/main.sayhello.nf.test"
assert snapshot(process.out).match()
```

Això indica a nf-test que creï un snapshot de la sortida del procés `sayHello`. Examinem el contingut del fitxer de snapshot.

```console title="Snapshot file contents"
code tests/main.sayhello.nf.test.snap
```

No l'imprimirem aquí, però hauríeu de veure un fitxer JSON que conté detalls del procés i les sortides del procés. En particular, podem veure una línia que sembla això:

```json title="Snapshot file contents"
"0": [
    "hello-output.txt:md5,b1946ac92492d2347c6235b4d2611184"
]
```

Això representa les sortides creades pel procés `sayHello`, que estem provant explícitament. Si tornem a executar la prova, el programa comprovarà que la nova sortida coincideix amb la sortida que es va registrar originalment. Aquesta és una manera ràpida i senzilla de provar que les sortides del procés no canvien, per això nf-test la proporciona com a opció per defecte.

!!!warning "Advertència"

    Això significa que hem d'assegurar-nos que la sortida que registrem en l'execució original és correcta!

Si, en el curs del desenvolupament futur, alguna cosa al codi canvia i fa que la sortida sigui diferent, la prova fallarà i haurem de determinar si el canvi és esperat o no.

- Si resulta que alguna cosa al codi s'ha trencat, haurem de corregir-ho, amb l'expectativa que el codi corregit passarà la prova.
- Si és un canvi esperat (p. ex., l'eina ha millorat i els resultats són millors), llavors haurem d'actualitzar el snapshot per acceptar la nova sortida com a referència a comparar. nf-test té un paràmetre `--update-snapshot` per a aquest propòsit.

Podem tornar a executar la prova i veure que hauria de passar:

```console title="nf-test process pass with snapshot"
> nf-test test tests/main.sayhello.nf.test

🚀 nf-test 0.9.3
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Process sayHello

  Test [f91a1bcd] 'Should run without failures and produce correct output' PASSED (1.675s)


SUCCESS: Executed 1 tests in 1.685s
```

Èxit! La prova passa perquè el procés `sayHello` s'ha executat correctament i la sortida coincideix amb el snapshot.

### 2.3. Alternativa als snapshots: assertions de contingut directes

Tot i que els snapshots són excel·lents per detectar qualsevol canvi a la sortida, de vegades voleu verificar contingut específic sense ser tan estrictes sobre la coincidència de tot el fitxer. Per exemple:

- Quan parts de la sortida poden canviar (marques de temps, IDs aleatoris, etc.) però cert contingut clau ha d'estar present
- Quan voleu comprovar patrons o valors específics a la sortida
- Quan voleu fer la prova més explícita sobre el que constitueix l'èxit

Aquí teniu com podríem modificar la nostra prova per comprovar contingut específic:

=== "Després"

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

=== "Abans"

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

Tingueu en compte que nf-test veu les sortides del procés com una llista de llistes, de manera que `process.out[0][0]` obté la primera part del primer element del canal (o 'emissió') d'aquest procés.

Aquest enfocament:

- Deixa clar exactament el que esperem a la sortida
- És més resistent als canvis irrellevants a la sortida
- Proporciona millors missatges d'error quan les proves fallen
- Permet validacions més complexes (patrons regex, comparacions numèriques, etc.)

Executem la prova per veure si funciona.

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

### 2.4. Proveu el procés `convertToUpper`

Obrim el fitxer `tests/main.converttoupper.nf.test` i examinem el contingut:

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

Aquesta és una prova similar a la del procés `sayHello`, però prova el procés `convertToUpper`. Sabem que aquesta fallarà perquè, igual que amb `sayHello`, el procés `convertToUpper` pren una única entrada de camí, però no n'hem especificat cap.

Ara hem de proporcionar un únic fitxer d'entrada al procés convertToUpper, que inclou algun text que volem convertir a majúscules. Hi ha moltes maneres de fer-ho:

- Podríem crear un fitxer dedicat per a la prova
- Podríem reutilitzar el fitxer data/greetings.csv existent
- Podríem crear-lo sobre la marxa dins de la prova

De moment, reutilitzem el fitxer data/greetings.csv existent utilitzant l'exemple que hem usat amb la prova a nivell de pipeline. Com abans, podem anomenar la prova per reflectir millor el que estem provant, però aquesta vegada deixem que faci un 'snapshot' del contingut en lloc de comprovar cadenes específiques (com hem fet amb l'altre procés).

=== "Després"

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

=== "Abans"

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

I executeu la prova!

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

Tingueu en compte que hem creat un fitxer de snapshot per al procés `convertToUpper` a `tests/main.converttoupper.nf.test.snap`. Si tornem a executar la prova, hauríem de veure que nf-test passa de nou.

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

### Conclusió

Sabeu com escriure proves per a un procés Nextflow i executar-les.

### Què segueix?

Apreneu com executar proves per a tot alhora!

## 3. Executeu proves per a tot el repositori

Executar nf-test en cada component és correcte, però és laboriós i propens a errors. No podem simplement provar-ho tot alhora?

Sí que podem!

Executem nf-test en tot el repositori.

### 3.1. Executeu nf-test en tot el repositori

Podem executar nf-test en tot el repositori executant la comanda `nf-test test`.

```bash
nf-test test .
```

Tingueu en compte que simplement estem utilitzant el `.` per executar tot des del nostre directori actual. Això inclourà totes les proves!

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

Mireu això! Hem executat 4 proves, 1 per a cada procés i 2 per a tot el pipeline amb una única comanda. Imagineu com de poderós és això en una base de codi gran!

---

## Resum

En aquesta missió secundària, heu après a aprofitar les funcionalitats de nf-test per crear i executar proves per a processos individuals, així com proves d'extrem a extrem per a tot el pipeline.
Ara coneixeu els dos enfocaments principals per a la validació de sortides, els snapshots i les assertions de contingut directes, i quan utilitzar cadascun.
També sabeu com executar proves una per una o per a tot un projecte.

Aplicar aquestes tècniques en el vostre propi treball us permetrà assegurar que:

- El vostre codi funciona com s'espera
- Els canvis no trenquen la funcionalitat existent
- Altres desenvolupadors poden contribuir amb confiança
- Els problemes es poden identificar i corregir ràpidament
- El contingut de la sortida coincideix amb les expectatives

### Patrons clau

<!-- TODO: Can we add snippets of code below to illustrate? -->

1. Proves a nivell de pipeline:
   - Prova bàsica d'èxit
   - Verificació del nombre de processos
   - Comprovació de l'existència de fitxers de sortida
2. Proves a nivell de procés
3. Dos enfocaments per a la validació de sortides:
   - Ús de snapshots per a la verificació completa de la sortida
   - Ús d'assertions de contingut directes per a comprovacions de contingut específic
4. Execució de totes les proves d'un repositori amb una única comanda

### Recursos addicionals

Consulteu la [documentació de nf-test](https://www.nf-test.com/) per a funcionalitats de proves més avançades i bones pràctiques. Potser voldreu:

- Afegir assertions més exhaustives a les vostres proves
- Escriure proves per a casos límit i condicions d'error
- Configurar la integració contínua per executar les proves automàticament
- Aprendre sobre altres tipus de proves com les proves de workflow i de mòdul
- Explorar tècniques de validació de contingut més avançades

**Recordeu:** Les proves són documentació viva de com hauria de comportar-se el vostre codi. Com més proves escriviu, i com més específiques siguin les vostres assertions, més segurs podeu estar de la fiabilitat del vostre pipeline.

---

## Què segueix?

Torneu al [menú de missions secundàries](../) o feu clic al botó a la part inferior dreta de la pàgina per passar al tema següent de la llista.
