# Part 2: Reescriure Hello per a nf-core

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducció assistida per IA - [més informació i suggeriments](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

En aquesta segona part del curs de formació Hello nf-core, us mostrem com crear una versió compatible amb nf-core del pipeline produït pel curs per a principiants [Hello Nextflow](../hello_nextflow/index.md).

Haureu notat a la primera secció de la formació que els pipelines nf-core segueixen una estructura força elaborada amb molts fitxers accessoris.
Crear tot això des de zero seria molt tediós, així que la comunitat nf-core ha desenvolupat eines per fer-ho a partir d'una plantilla, per iniciar el procés.

Us mostrarem com utilitzar aquestes eines per crear una estructura de pipeline i després adaptar el codi de pipeline 'regular' existent a l'estructura nf-core.

Si no esteu familiaritzats amb el pipeline Hello o us caldria un recordatori, consulteu [aquesta pàgina d'informació](../info/hello_pipeline.md).

---

## 1. Crear un nou projecte de pipeline

Primer, creem l'estructura per al nou pipeline.

!!! note "Nota"

    Assegureu-vos que esteu al directori `hello-nf-core` al vostre terminal.

### 1.1. Executar l'eina de creació de pipelines basada en plantilla

Comencem creant un nou pipeline amb la comanda `nf-core pipelines create`.
Això crearà una nova estructura de pipeline utilitzant la plantilla base nf-core, personalitzada amb un nom de pipeline, descripció i autor.

```bash
nf-core pipelines create
```

Executar aquesta comanda obrirà una Interfície d'Usuari de Text (TUI) per a la creació del pipeline:

<div style="text-align: center;">
    <iframe width="560" height="315" src="https://www.youtube.com/embed/VwjXNXONHlY?si=d0HkFSISnKn76TeI" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" allowfullscreen="" data-ruffle-polyfilled=""></iframe>
</div>

Aquesta TUI us demanarà que proporcioneu informació bàsica sobre el vostre pipeline i us oferirà una selecció de funcionalitats per incloure o excloure a l'estructura del vostre pipeline.

- A la pantalla de benvinguda, feu clic a **Let's go!**.
- A la pantalla `Choose pipeline type`, feu clic a **Custom**.
- Introduïu els detalls del vostre pipeline de la següent manera (substituint `< EL VOSTRE NOM >` pel vostre propi nom), després feu clic a **Next**.

```
[ ] GitHub organisation: core
[ ] Workflow name: hello
[ ] A short description of your pipeline: A basic nf-core style version of Hello Nextflow
[ ] Name of the main author(s): < YOUR NAME >
```

- A la pantalla Template features, establiu `Toggle all features` a **off**, després **activeu** selectivament els següents. Comproveu les vostres seleccions i feu clic a **Continue**.

```
[ ] Add testing profiles
[ ] Use nf-core components
[ ] Use nf-schema
[ ] Add configuration files
[ ] Add documentation
```

- A la pantalla `Final details`, feu clic a **Finish**. Espereu que es creï el pipeline, després feu clic a **Continue**.
- A la pantalla Create GitHub repository, feu clic a **Finish without creating a repo**. Això mostrarà instruccions per crear un repositori GitHub més endavant. Ignoreu-les i feu clic a **Close**.

Un cop es tanqui la TUI, hauríeu de veure la següent sortida a la consola.

??? success "Sortida de la comanda"

    ```console
                                              ,--./,-.
              ___     __   __   __   ___     /,-._.--~\
        |\ | |__  __ /  ` /  \ |__) |__         }  {
        | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                              `._,._,'

        nf-core/tools version 3.4.1 - https://nf-co.re


    INFO     Launching interactive nf-core pipeline creation tool.
    ```

No hi ha cap confirmació explícita a la sortida de la consola que la creació del pipeline hagi funcionat, però hauríeu de veure un nou directori anomenat `core-hello`.

Visualitzeu el contingut del nou directori per veure quanta feina us heu estalviat utilitzant la plantilla.

```bash
tree core-hello
```

??? abstract "Contingut del directori"

    ```console
    core-hello/
    ├── assets
    │   ├── samplesheet.csv
    │   └── schema_input.json
    ├── conf
    │   ├── base.config
    │   ├── modules.config
    │   ├── test.config
    │   └── test_full.config
    ├── docs
    │   ├── output.md
    │   ├── README.md
    │   └── usage.md
    ├── main.nf
    ├── modules.json
    ├── nextflow.config
    ├── nextflow_schema.json
    ├── README.md
    ├── subworkflows
    │   ├── local
    │   │   └── utils_nfcore_hello_pipeline
    │   │       └── main.nf
    │   └── nf-core
    │       ├── utils_nextflow_pipeline
    │       │   ├── main.nf
    │       │   ├── meta.yml
    │       │   └── tests
    │       │       ├── main.function.nf.test
    │       │       ├── main.function.nf.test.snap
    │       │       ├── main.workflow.nf.test
    │       │       └── nextflow.config
    │       ├── utils_nfcore_pipeline
    │       │   ├── main.nf
    │       │   ├── meta.yml
    │       │   └── tests
    │       │       ├── main.function.nf.test
    │       │       ├── main.function.nf.test.snap
    │       │       ├── main.workflow.nf.test
    │       │       ├── main.workflow.nf.test.snap
    │       │       └── nextflow.config
    │       └── utils_nfschema_plugin
    │           ├── main.nf
    │           ├── meta.yml
    │           └── tests
    │               ├── main.nf.test
    │               ├── nextflow.config
    │               └── nextflow_schema.json
    └── workflows
        └── hello.nf

    14 directories, 34 files
    ```

Són molts fitxers!

Esperem que reconegueu molts d'ells com els mateixos que vam trobar quan vam explorar l'estructura del pipeline `nf-core/demo`.
Però no us preocupeu si encara us sentiu una mica perduts; recorrerem les parts importants junts durant aquest curs de formació.

!!! note "Nota"

    Una diferència important en comparació amb el pipeline `nf-core/demo` que vam examinar a la primera part d'aquesta formació és que no hi ha cap directori `modules`.
    Això és perquè no vam optar per incloure cap dels mòduls nf-core per defecte.

### 1.2. Comprovar que l'estructura és funcional

Cregueu-ho o no, tot i que encara no heu afegit cap mòdul per fer-lo fer feina real, l'estructura del pipeline es pot executar utilitzant el perfil de test, de la mateixa manera que vam executar el pipeline `nf-core/demo`.

```bash
nextflow run ./core-hello -profile docker,test --outdir core-hello-results
```

??? success "Sortida de la comanda"

    ```console
    N E X T F L O W   ~  version 25.04.3

    Launching `./core-hello/main.nf` [scruffy_marconi] DSL2 - revision: b9e9b3b8de

    Downloading plugin nf-schema@2.5.1
    Input/output options
      input                     : https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/samplesheet/samplesheet_test_illumina_amplicon.csv
      outdir                    : core-hello-results

    Institutional config options
      config_profile_name       : Test profile
      config_profile_description: Minimal test dataset to check pipeline function

    Generic options
      trace_report_suffix       : 2025-11-21_04-47-18

    Core Nextflow options
      runName                   : scruffy_marconi
      containerEngine           : docker
      launchDir                 : /workspaces/training/hello-nf-core
      workDir                   : /workspaces/training/hello-nf-core/work
      projectDir                : /workspaces/training/hello-nf-core/core-hello
      userName                  : root
      profile                   : docker,test
      configFiles               : /workspaces/training/hello-nf-core/core-hello/nextflow.config

    !! Only displaying parameters that differ from the pipeline defaults !!
    ------------------------------------------------------
    -[core/hello] Pipeline completed successfully-
    ```

Això us mostra que tot el cablejat bàsic està en el seu lloc.
Així doncs, on són les sortides? N'hi ha cap?

De fet, es va crear un nou directori de resultats anomenat `core-hello-results` que conté els informes d'execució estàndard:

```bash
tree core-hello-results
```

??? abstract "Contingut del directori"

    ```console
    core-hello-results
    └── pipeline_info
        ├── execution_report_2025-11-21_04-47-18.html
        ├── execution_timeline_2025-11-21_04-47-18.html
        ├── execution_trace_2025-11-21_04-47-18.txt
        ├── hello_software_versions.yml
        ├── params_2025-11-21_04-47-18.json
        └── pipeline_dag_2025-11-21_04-47-18.html

    1 directory, 6 files
    ```

Podeu donar una ullada als informes per veure què s'ha executat, i la resposta és: res en absolut!

![informe de línia de temps d'execució buit](./img/execution_timeline_empty.png)

Donem una ullada al que hi ha realment al codi.

### 1.3. Examinar el workflow de marcador de posició

Si mireu dins del fitxer `main.nf`, veureu que importa un workflow anomenat `HELLO` de `workflows/hello`.

Això és equivalent al workflow `workflows/demo.nf` que vam trobar a la Part 1, i serveix com a workflow de marcador de posició per al nostre workflow d'interès, amb alguna funcionalitat nf-core ja en el seu lloc.

```groovy title="core-hello/workflows/hello.nf" linenums="1" hl_lines="15 17 19 35"
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow HELLO {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    main:

    ch_versions = channel.empty()

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name:  'hello_software_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    emit:
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
```

En comparació amb un workflow bàsic de Nextflow com el desenvolupat a [Hello Nextflow](../hello_nextflow/index.md), notareu algunes coses que són noves aquí (línies destacades anteriorment):

- El bloc workflow té un nom
- Les entrades del workflow es declaren utilitzant la paraula clau `take:` i la construcció del canal es mou cap al workflow pare
- El contingut del workflow es col·loca dins d'un bloc `main:`
- Les sortides es declaren utilitzant la paraula clau `emit:`

Aquestes són funcionalitats opcionals de Nextflow que fan que el workflow sigui **composable**, és a dir, que es pot cridar des de dins d'un altre workflow.

!!! note "Workflows composables en profunditat"

    La [Workflows of Workflows](../side_quests/workflows_of_workflows.md) Side Quest explora la composició de workflows amb molt més detall, incloent com compondre múltiples workflows junts i gestionar fluxos de dades complexos entre ells. Introduïm la composabilitat aquí perquè és un requisit fonamental de l'arquitectura de plantilla nf-core, que utilitza workflows niats per organitzar la inicialització del pipeline, el workflow d'anàlisi principal i les tasques de finalització en components separats i reutilitzables.

Necessitarem connectar la lògica rellevant del nostre workflow d'interès a aquesta estructura.
El primer pas per a això és fer que el nostre workflow original sigui composable.

### Conclusió

Ara sabeu com crear una estructura de pipeline utilitzant les eines nf-core.

### Què segueix?

Apreneu com fer que un workflow simple sigui composable com a preludi per fer-lo compatible amb nf-core.

---

## 2. Fer el workflow original Hello Nextflow composable

Ara és el moment de posar-se a treballar integrant el nostre workflow a l'estructura nf-core.
Com a recordatori, estem treballant amb el workflow presentat al nostre curs de formació [Hello Nextflow](../hello_nextflow/index.md).

!!! tip "Consell"

    Si no esteu familiaritzats amb aquest pipeline o us caldria un recordatori, consulteu [El pipeline Hello](../info/hello_pipeline.md).

Us proporcionem una còpia neta i totalment funcional del workflow Hello Nextflow completat al directori `original-hello` juntament amb els seus mòduls i el fitxer CSV per defecte que espera utilitzar com a entrada.

```bash
tree original-hello/
```

??? abstract "Contingut del directori"

    ```console
    original-hello/
    ├── hello.nf
    ├── modules
    │   ├── collectGreetings.nf
    │   ├── convertToUpper.nf
    │   ├── cowpy.nf
    │   └── sayHello.nf
    └── nextflow.config

    1 directory, 6 files
    ```

Sentiu-vos lliures d'executar-lo per assegurar-vos que funciona:

```bash
nextflow run original-hello/hello.nf
```

??? success "Sortida de la comanda"

    ```console
    N E X T F L O W   ~  version 25.04.3

    Launching `original-hello/hello.nf` [goofy_babbage] DSL2 - revision: e9e72441e9

    executor >  local (8)
    [a4/081cec] sayHello (1)       | 3 of 3 ✔
    [e7/7e9058] convertToUpper (3) | 3 of 3 ✔
    [0c/17263b] collectGreetings   | 1 of 1 ✔
    [94/542280] cowpy              | 1 of 1 ✔
    ```

Obrim el fitxer de workflow `hello.nf` per inspeccionar el codi, que es mostra complet a continuació (sense comptar els processos, que estan en mòduls):

```groovy title="original-hello/hello.nf" linenums="1"
#!/usr/bin/env nextflow

/*
* Pipeline parameters
*/
params.greeting = 'greetings.csv'
params.batch = 'test-batch'
params.character = 'turkey'

// Include modules
include { sayHello } from './modules/sayHello.nf'
include { convertToUpper } from './modules/convertToUpper.nf'
include { collectGreetings } from './modules/collectGreetings.nf'
include { cowpy } from './modules/cowpy.nf'

workflow {

  // create a channel for inputs from a CSV file
  greeting_ch = channel.fromPath(params.greeting)
                      .splitCsv()
                      .map { line -> line[0] }

  // emit a greeting
  sayHello(greeting_ch)

  // convert the greeting to uppercase
  convertToUpper(sayHello.out)

  // collect all the greetings into one file
  collectGreetings(convertToUpper.out.collect(), params.batch)

  // generate ASCII art of the greetings with cowpy
  cowpy(collectGreetings.out.outfile, params.character)
}
```

Com podeu veure, aquest workflow es va escriure com un workflow simple sense nom que es pot executar per si sol.
Per tal de fer-lo executable des de dins d'un workflow pare com requereix la plantilla nf-core, necessitem fer-lo **composable**.

Recorrem els canvis necessaris un per un.

### 2.1. Donar nom al workflow

Primer, donem un nom al workflow perquè puguem referir-nos-hi des d'un workflow pare.

=== "Després"

    ```groovy title="original-hello/hello.nf" linenums="16"
    workflow HELLO {
    ```

=== "Abans"

    ```groovy title="original-hello/hello.nf" linenums="16"
    workflow {
    ```

Les mateixes convencions s'apliquen als noms de workflow que als noms de mòdul.

### 2.2. Substituir la construcció del canal per `take`

Ara, substituïu la construcció del canal per una simple declaració `take` que declara les entrades esperades.

=== "Després"

    ```groovy title="original-hello/hello.nf" linenums="18"
        take:
        // canal de salutacions
        greeting_ch
    ```

=== "Abans"

    ```groovy title="original-hello/hello.nf" linenums="18"
        // create a channel for inputs from a CSV file
        greeting_ch = channel.fromPath(params.greeting)
                            .splitCsv()
                            .map { line -> line[0] }
    ```

Això deixa els detalls de com es proporcionen les entrades al workflow pare.

Mentre hi som, també podem comentar la línia `params.greeting = 'greetings.csv'`

=== "Després"

    ```groovy title="original-hello/hello.nf" linenums="3" hl_lines="4"
        /*
        * Pipeline parameters
        */
        //params.greeting = 'greetings.csv'
        params.batch = 'test-batch'
        params.character = 'turkey'
    ```

=== "Abans"

    ```groovy title="original-hello/hello.nf" linenums="3" hl_lines="4"
        /*
        * Pipeline parameters
        */
        params.greeting = 'greetings.csv'
        params.batch = 'test-batch'
        params.character = 'turkey'
    ```

!!! note "Nota"

    Si teniu instal·lada l'extensió del servidor de llenguatge Nextflow, el verificador de sintaxi il·luminarà el vostre codi amb línies ondulades vermelles.
    Això és perquè si poseu una declaració `take:`, també heu de tenir un `main:`.

    Ho afegirem al següent pas.

### 2.3. Prefaciar les operacions del workflow amb la declaració `main`

A continuació, afegiu una declaració `main` abans de la resta d'operacions cridades al cos del workflow.

=== "Després"

    ```groovy title="original-hello/hello.nf" linenums="22" hl_lines="1"
        main:

        // emet una salutacio
        sayHello(greeting_ch)

        // converteix la salutacio a majuscules
        convertToUpper(sayHello.out)

        // recull totes les salutacions en un fitxer
        collectGreetings(convertToUpper.out.collect(), params.batch)

        // genera art ASCII de les salutacions amb cowpy
        cowpy(collectGreetings.out.outfile, params.character)
    ```

=== "Abans"

    ```groovy title="original-hello/hello.nf" linenums="21"
        // emit a greeting
        sayHello(greeting_ch)

        // convert the greeting to uppercase
        convertToUpper(sayHello.out)

        // collect all the greetings into one file
        collectGreetings(convertToUpper.out.collect(), params.batch)

        // generate ASCII art of the greetings with cowpy
        cowpy(collectGreetings.out.outfile, params.character)
    ```

Això bàsicament diu 'això és el que _fa_ aquest workflow'.

### 2.4. Afegir la declaració `emit`

Finalment, afegiu una declaració `emit` que declara quines són les sortides finals del workflow.

```groovy title="original-hello/hello.nf" linenums="35"
    emit:
    cowpy_hellos = cowpy.out
```

Aquesta és una addició completament nova al codi en comparació amb el workflow original.

### 2.5. Resum dels canvis completats

Si heu fet tots els canvis tal com es descriuen, el vostre workflow ara hauria de semblar així:

```groovy title="original-hello/hello.nf" linenums="1" hl_lines="16 18-20 22 36-37"
#!/usr/bin/env nextflow

/*
* Pipeline parameters
*/
// params.greeting = 'greetings.csv'
params.batch = 'test-batch'
params.character = 'turkey'

// Include modules
include { sayHello } from './modules/sayHello.nf'
include { convertToUpper } from './modules/convertToUpper.nf'
include { collectGreetings } from './modules/collectGreetings.nf'
include { cowpy } from './modules/cowpy.nf'

workflow HELLO {

    take:
    // canal de salutacions
    greeting_ch

    main:

    // emet una salutacio
    sayHello(greeting_ch)

    // converteix la salutacio a majuscules
    convertToUpper(sayHello.out)

    // recull totes les salutacions en un fitxer
    collectGreetings(convertToUpper.out.collect(), params.batch)

    // genera art ASCII de les salutacions amb cowpy
    cowpy(collectGreetings.out.outfile, params.character)

    emit:
    cowpy_hellos = cowpy.out
}
```

Això descriu tot el que Nextflow necessita EXCEPTE què alimentar al canal d'entrada.
Això es definirà al workflow pare, també anomenat workflow **entrypoint**.

### 2.6. Fer un workflow entrypoint fictici

Abans d'integrar el nostre workflow composable a l'estructura complexa nf-core, verifiquem que funciona correctament.
Podem fer un workflow entrypoint fictici simple per provar el workflow composable de manera aïllada.

Creeu un fitxer en blanc anomenat `main.nf` al mateix directori `original-hello`.

```bash
touch original-hello/main.nf
```

Copieu el següent codi al fitxer `main.nf`.

```groovy title="original-hello/main.nf" linenums="1"
#!/usr/bin/env nextflow

// importa el codi del workflow des del fitxer hello.nf
include { HELLO } from './hello.nf'

// declara el parametre d'entrada
params.greeting = 'greetings.csv'

workflow {
  // crea un canal per a entrades des d'un fitxer CSV
  greeting_ch = channel.fromPath(params.greeting)
                      .splitCsv()
                      .map { line -> line[0] }

  // crida el workflow importat sobre el canal de salutacions
  HELLO(greeting_ch)

  // visualitza les sortides emeses pel workflow
  HELLO.out.view { output -> "Output: $output" }
}
```

Hi ha dues observacions importants a fer aquí:

- La sintaxi per cridar el workflow importat és essencialment la mateixa que la sintaxi per cridar mòduls.
- Tot el que està relacionat amb portar les entrades al workflow (parametre d'entrada i construcció del canal) ara es declara en aquest workflow pare.

!!! note "Nota"

    Anomenar el fitxer de workflow entrypoint `main.nf` és una convenció, no un requisit.

    Si seguiu aquesta convenció, podeu ometre especificar el nom del fitxer de workflow a la vostra comanda `nextflow run`.
    Nextflow buscarà automàticament un fitxer anomenat `main.nf` al directori d'execució.

    No obstant això, podeu anomenar el fitxer de workflow entrypoint d'una altra manera si ho preferiu.
    En aquest cas, assegureu-vos d'especificar el nom del fitxer de workflow a la vostra comanda `nextflow run`.

### 2.7. Comprovar que el workflow s'executa

Finalment tenim totes les peces que necessitem per verificar que el workflow composable funciona.
Executem-lo!

```bash
nextflow run ./original-hello
```

Aquí veieu l'avantatge d'utilitzar la convenció de nomenclatura `main.nf`.
Si haguéssim anomenat el workflow entrypoint `something_else.nf`, hauríem hagut de fer `nextflow run original-hello/something_else.nf`.

Si heu fet tots els canvis correctament, això hauria de completar-se.

??? success "Sortida de la comanda"

    ```console
    N E X T F L O W   ~  version 25.04.3

    Launching `original-hello/main.nf` [friendly_wright] DSL2 - revision: 1ecd2d9c0a

    executor >  local (8)
    [24/c6c0d8] HELLO:sayHello (3)       | 3 of 3 ✔
    [dc/721042] HELLO:convertToUpper (3) | 3 of 3 ✔
    [48/5ab2df] HELLO:collectGreetings   | 1 of 1 ✔
    [e3/693b7e] HELLO:cowpy              | 1 of 1 ✔
    Output: /workspaces/training/hello-nf-core/work/e3/693b7e48dc119d0c54543e0634c2e7/cowpy-COLLECTED-test-batch-output.txt
    ```

Això significa que hem actualitzat amb èxit el nostre workflow HELLO per ser composable.

### Conclusió

Sabeu com fer un workflow composable donant-li un nom i afegint declaracions `take`, `main` i `emit`, i com cridar-lo des d'un workflow entrypoint.

### Què segueix?

Apreneu com injertar un workflow composable bàsic a l'estructura nf-core.

---

## 3. Ajustar la lògica del workflow actualitzat al workflow de marcador de posició

Ara que hem verificat que el nostre workflow composable funciona correctament, tornem a l'estructura del pipeline nf-core que vam crear a la secció 1.
Volem integrar el workflow composable que acabem de desenvolupar a l'estructura de plantilla nf-core, de manera que el resultat final hauria de semblar així.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nf-core/img/core-hello.svg"
</figure>

Així doncs, com ho fem? Donem una ullada al contingut actual del workflow `HELLO` a `core-hello/workflows/hello.nf` (l'estructura nf-core).

```groovy title="core-hello/workflows/hello.nf" linenums="1"
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow HELLO {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    main:

    ch_versions = channel.empty()

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name:  'hello_software_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    emit:
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
```

En general, aquest codi fa molt poc a part d'algunes tasques de manteniment que tenen a veure amb capturar la versió de qualsevol eina de programari que s'executi al pipeline.

Necessitem afegir el codi rellevant de la versió composable del workflow original que vam desenvolupar a la secció 2.

Abordarem això en les següents etapes:

1. Copiar els mòduls i configurar les importacions de mòduls
2. Deixar la declaració `take` tal com està
3. Afegir la lògica del workflow al bloc `main`
4. Actualitzar el bloc `emit`

!!! note "Nota"

    Ignorarem la captura de versions per a aquesta primera passada i veurem com connectar-ho en una part posterior d'aquesta formació.

### 3.1. Copiar els mòduls i configurar les importacions de mòduls

Els quatre processos del nostre workflow Hello Nextflow s'emmagatzemen com a mòduls a `original-hello/modules/`.
Necessitem copiar aquests mòduls a l'estructura del projecte nf-core (sota `core-hello/modules/local/`) i afegir declaracions d'importació al fitxer de workflow nf-core.

Primer copiem els fitxers de mòdul de `original-hello/` a `core-hello/`:

```bash
mkdir -p core-hello/modules/local/
cp original-hello/modules/* core-hello/modules/local/.
```

Ara hauríeu de veure el directori de mòduls llistat sota `core-hello/`.

```bash
tree core-hello/modules
```

??? abstract "Contingut del directori"

    ```console
    core-hello/modules
    └── local
        ├── collectGreetings.nf
        ├── convertToUpper.nf
        ├── cowpy.nf
        └── sayHello.nf

    1 directory, 4 files
    ```

Ara configurem les declaracions d'importació de mòduls.

Aquestes eren les declaracions d'importació al workflow `original-hello/hello.nf`:

```groovy title="original-hello/hello.nf" linenums="9"
// Include modules
include { sayHello } from './modules/sayHello.nf'
include { convertToUpper } from './modules/convertToUpper.nf'
include { collectGreetings } from './modules/collectGreetings.nf'
include { cowpy } from './modules/cowpy.nf'
```

Obriu el fitxer `core-hello/workflows/hello.nf` i transposeu aquestes declaracions d'importació tal com es mostra a continuació.

=== "Després"

    ```groovy title="core-hello/workflows/hello.nf" linenums="1" hl_lines="8-11"
    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    include { paramsSummaryMap       } from 'plugin/nf-schema'
    include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
    include { sayHello               } from '../modules/local/sayHello.nf'
    include { convertToUpper         } from '../modules/local/convertToUpper.nf'
    include { collectGreetings       } from '../modules/local/collectGreetings.nf'
    include { cowpy                  } from '../modules/local/cowpy.nf'
    ```

=== "Abans"

    ```groovy title="core-hello/workflows/hello.nf" linenums="1"
    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    include { paramsSummaryMap       } from 'plugin/nf-schema'
    include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
    ```

Dues observacions més interessants aquí:

- Hem adaptat el format de les declaracions d'importació per seguir la convenció d'estil nf-core.
- Hem actualitzat els camins relatius als mòduls per reflectir que ara s'emmagatzemen a un nivell diferent de niament.

### 3.2. Deixar la declaració `take` tal com està

El projecte nf-core té molta funcionalitat preconstruïda al voltant del concepte del samplesheet, que normalment és un fitxer CSV que conté dades en columnes.
Com que això és essencialment el que és el nostre fitxer `greetings.csv`, mantindrem la declaració `take` actual tal com està, i simplement actualitzarem el nom del canal d'entrada al següent pas.

```groovy title="core-hello/workflows/hello.nf" linenums="21"
    take:
    ch_samplesheet // channel: samplesheet read in from --input
```

El maneig de l'entrada es farà abans d'aquest workflow (no en aquest fitxer de codi).

### 3.3. Afegir la lògica del workflow al bloc `main`

Ara que els nostres mòduls estan disponibles per al workflow, podem connectar la lògica del workflow al bloc `main`.

Com a recordatori, aquest és el codi rellevant al workflow original, que no va canviar gaire quan el vam fer composable (només vam afegir la línia `main:`):

```groovy title="original-hello/hello.nf" linenums="22"
    main:

    // emet una salutacio
    sayHello(greeting_ch)

    // converteix la salutacio a majuscules
    convertToUpper(sayHello.out)

    // recull totes les salutacions en un fitxer
    collectGreetings(convertToUpper.out.collect(), params.batch)

    // genera art ASCII de les salutacions amb cowpy
    cowpy(collectGreetings.out.outfile, params.character)
```

Necessitem copiar el codi que ve després de `main:` a la nova versió del workflow.

Ja hi ha algun codi allà que té a veure amb capturar les versions de les eines que executa el workflow. Ho deixarem estar per ara (tractarem les versions de les eines més endavant).
Mantindrem la inicialització `ch_versions = channel.empty()` a la part superior, després inserirem la nostra lògica de workflow, mantenint el codi de recopilació de versions al final.
Aquest ordre té sentit perquè en un pipeline real, els processos emetrien informació de versió que s'afegiria al canal `ch_versions` mentre s'executa el workflow.

=== "Després"

    ```groovy title="core-hello/workflows/hello.nf" linenums="19" hl_lines="10-20"
    workflow HELLO {

        take:
        ch_samplesheet // channel: samplesheet read in from --input

        main:

        ch_versions = Channel.empty()

        // emet una salutacio
        sayHello(greeting_ch)

        // converteix la salutacio a majuscules
        convertToUpper(sayHello.out)

        // recull totes les salutacions en un fitxer
        collectGreetings(convertToUpper.out.collect(), params.batch)

        // genera art ASCII de les salutacions amb cowpy
        cowpy(collectGreetings.out.outfile, params.character)

        //
        // Collate and save software versions
        //
        softwareVersionsToYAML(ch_versions)
            .collectFile(
                storeDir: "${params.outdir}/pipeline_info",
                name:  'hello_software_'  + 'versions.yml',
                sort: true,
                newLine: true
            ).set { ch_collated_versions }


        emit:
        versions       = ch_versions                 // channel: [ path(versions.yml) ]

    }
    ```

=== "Abans"

    ```groovy title="core-hello/workflows/hello.nf" linenums="19"
    workflow HELLO {

        take:
        ch_samplesheet // channel: samplesheet read in from --input
        main:

        ch_versions = Channel.empty()

        //
        // Collate and save software versions
        //
        softwareVersionsToYAML(ch_versions)
            .collectFile(
                storeDir: "${params.outdir}/pipeline_info",
                name:  'hello_software_'  + 'versions.yml',
                sort: true,
                newLine: true
            ).set { ch_collated_versions }


        emit:
        versions       = ch_versions                 // channel: [ path(versions.yml) ]

    }
    ```

Notareu que també hem afegit una línia en blanc abans de `main:` per fer el codi més llegible.

Això es veu genial, però encara necessitem actualitzar el nom del canal que estem passant al procés `sayHello()` de `greeting_ch` a `ch_samplesheet` tal com es mostra a continuació, per coincidir amb el que està escrit sota la paraula clau `take:`.

=== "Després"

    ```groovy title="core-hello/workflows/hello.nf" linenums="26"
        // emet una salutacio (actualitzat per utilitzar la convencio nf-core per a samplesheets)
        sayHello(ch_samplesheet)
    ```

=== "Abans"

    ```groovy title="core-hello/workflows/hello.nf" linenums="26"
        // emet una salutacio
        sayHello(greeting_ch)
    ```

Ara la lògica del workflow està correctament connectada.

### 3.4. Actualitzar el bloc `emit`

Finalment, necessitem actualitzar el bloc `emit` per incloure la declaració de les sortides finals del workflow.

=== "Després"

    ```groovy title="core-hello/workflows/hello.nf" linenums="55" hl_lines="2"
        emit:
        cowpy_hellos   = cowpy.out
        versions       = ch_versions                 // channel: [ path(versions.yml) ]
    ```

=== "Abans"

    ```groovy title="core-hello/workflows/hello.nf" linenums="55"
        emit:
        versions       = ch_versions                 // channel: [ path(versions.yml) ]
    ```

Això conclou les modificacions que necessitem fer al workflow HELLO en si.
En aquest punt, hem aconseguit l'estructura general del codi que ens vam proposar implementar.

### Conclusió

Sabeu com ajustar les peces centrals d'un workflow composable a un workflow de marcador de posició nf-core.

### Què segueix?

Apreneu com adaptar com es gestionen les entrades a l'estructura del pipeline nf-core.

---

## 4. Adaptar el maneig d'entrades

Ara que hem integrat amb èxit la nostra lògica de workflow a l'estructura nf-core, necessitem abordar una peça més crítica: assegurar-nos que les nostres dades d'entrada es processen correctament.
La plantilla nf-core ve amb un maneig d'entrades sofisticat dissenyat per a conjunts de dades genòmiques complexes, així que necessitem adaptar-lo per treballar amb el nostre fitxer `greetings.csv` més simple.

### 4.1. Identificar on es gestionen les entrades

El primer pas és esbrinar on es fa el maneig d'entrades.

Potser recordeu que quan vam reescriure el workflow Hello Nextflow per ser composable, vam moure la declaració del parametre d'entrada un nivell cap amunt, al workflow entrypoint `main.nf`.
Així doncs, donem una ullada al workflow entrypoint `main.nf` de nivell superior que es va crear com a part de l'estructura del pipeline:

```groovy title="core-hello/main.nf" linenums="1"
#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    core/hello
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/core/hello
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { HELLO  } from './workflows/hello'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_hello_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_hello_pipeline'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow CORE_HELLO {

    take:
    samplesheet // channel: samplesheet read in from --input

    main:

    //
    // WORKFLOW: Run pipeline
    //
    HELLO (
        samplesheet
    )
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:
    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input
    )

    //
    // WORKFLOW: Run main workflow
    //
    CORE_HELLO (
        PIPELINE_INITIALISATION.out.samplesheet
    )
    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.outdir,
        params.monochrome_logs,
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
```

El projecte nf-core fa un ús intensiu de subworkflows niats, així que aquesta part pot ser una mica confusa a la primera aproximació.

El que importa aquí és que hi ha dos workflows definits:

- `CORE_HELLO` és un embolcall prim per executar el workflow HELLO que acabem d'acabar d'adaptar a `core-hello/workflows/hello.nf`.
- Un workflow sense nom que crida `CORE_HELLO` així com dos altres subworkflows, `PIPELINE_INITIALISATION` i `PIPELINE_COMPLETION`.

Aquí hi ha un diagrama de com es relacionen entre ells:

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nf-core/img/hello-nested-workflows.svg"
</figure>

Importantment, no podem trobar cap codi que construeixi un canal d'entrada a aquest nivell, només referències a un samplesheet proporcionat mitjançant el parametre `--input`.

Una mica d'investigació revela que el maneig d'entrades es fa pel subworkflow `PIPELINE_INITIALISATION`, apropiadament, que s'importa de `core-hello/subworkflows/local/utils_nfcore_hello_pipeline/main.nf`.

Si obrim aquest fitxer i desplacem cap avall, arribem a aquest fragment de codi:

```groovy title="core-hello/subworkflows/local/utils_nfcore_hello_pipeline/main.nf" linenums="76"
    //
    // Create channel from input file provided through params.input
    //

    channel
        .fromList(samplesheetToList(params.input, "${projectDir}/assets/schema_input.json"))
        .map {
            meta, fastq_1, fastq_2 ->
                if (!fastq_2) {
                    return [ meta.id, meta + [ single_end:true ], [ fastq_1 ] ]
                } else {
                    return [ meta.id, meta + [ single_end:false ], [ fastq_1, fastq_2 ] ]
                }
        }
        .groupTuple()
        .map { samplesheet ->
            validateInputSamplesheet(samplesheet)
        }
        .map {
            meta, fastqs ->
                return [ meta, fastqs.flatten() ]
        }
        .set { ch_samplesheet }

    emit:
    samplesheet = ch_samplesheet
    versions    = ch_versions
```

Aquesta és la factoria de canals que analitza el samplesheet i el passa en una forma que està llesta per ser consumida pel workflow HELLO.

!!! note "Nota"

    La sintaxi anterior és una mica diferent del que hem utilitzat anteriorment, però bàsicament això:

    ```groovy
    channel.<...>.set { ch_samplesheet }
    ```

    és equivalent a això:

    ```groovy
    ch_samplesheet = channel.<...>
    ```

Aquest codi implica alguns passos d'anàlisi i validació que són altament específics del samplesheet d'exemple inclòs amb la plantilla del pipeline nf-core, que en el moment d'escriure això és molt específic del domini i no és adequat per al nostre projecte de pipeline simple.

### 4.2. Substituir el codi del canal d'entrada de la plantilla

La bona notícia és que les necessitats del nostre pipeline són molt més simples, així que podem substituir tot això pel codi de construcció de canal que vam desenvolupar al workflow original Hello Nextflow.

Com a recordatori, així és com es veia la construcció del canal (tal com es veu al directori de solucions):

```groovy title="solutions/composable-hello/main.nf" linenums="10" hl_lines="4"
    // crea un canal per a entrades des d'un fitxer CSV
    greeting_ch = channel.fromPath(params.greeting)
        .splitCsv()
        .map { line -> line[0] }
```

Així que només necessitem connectar això al workflow d'inicialització, amb canvis menors: actualitzem el nom del canal de `greeting_ch` a `ch_samplesheet`, i el nom del parametre de `params.greeting` a `params.input` (vegeu la línia destacada).

=== "Després"

    ```groovy title="core-hello/subworkflows/local/utils_nfcore_hello_pipeline/main.nf" linenums="76" hl_lines="5-7"
        //
        // Create channel from input file provided through params.input
        //

        ch_samplesheet = channel.fromPath(params.input)
            .splitCsv()
            .map { line -> line[0] }

        emit:
        samplesheet = ch_samplesheet
        versions    = ch_versions
    ```

=== "Abans"

    ```groovy title="core-hello/subworkflows/local/utils_nfcore_hello_pipeline/main.nf" linenums="76" hl_lines="5-23"
        //
        // Create channel from input file provided through params.input
        //

        channel
            .fromList(samplesheetToList(params.input, "${projectDir}/assets/schema_input.json"))
            .map {
                meta, fastq_1, fastq_2 ->
                    if (!fastq_2) {
                        return [ meta.id, meta + [ single_end:true ], [ fastq_1 ] ]
                    } else {
                        return [ meta.id, meta + [ single_end:false ], [ fastq_1, fastq_2 ] ]
                    }
            }
            .groupTuple()
            .map { samplesheet ->
                validateInputSamplesheet(samplesheet)
            }
            .map {
                meta, fastqs ->
                    return [ meta, fastqs.flatten() ]
            }
            .set { ch_samplesheet }

        emit:
        samplesheet = ch_samplesheet
        versions    = ch_versions
    ```

Això completa els canvis que necessitem per fer que el processament d'entrades funcioni.

En la seva forma actual, això no ens permetrà aprofitar les capacitats integrades d'nf-core per a la validació d'esquemes, però podem afegir-ho més endavant.
Per ara, ens centrem en mantenir-ho tan simple com sigui possible per arribar a alguna cosa que puguem executar amb èxit amb dades de prova.

### 4.3. Actualitzar el perfil de test

Parlant de dades de prova i parametres, actualitzem el perfil de test per a aquest pipeline per utilitzar el mini-samplesheet `greetings.csv` en lloc del samplesheet d'exemple proporcionat a la plantilla.

Sota `core-hello/conf`, trobem dos perfils de test de plantilla: `test.config` i `test_full.config`, que estan destinats a provar una mostra de dades petita i una de mida completa.
Donat el propòsit del nostre pipeline, realment no té sentit configurar un perfil de test de mida completa, així que sentiu-vos lliures d'ignorar o eliminar `test_full.config`.
Ens centrarem en configurar `test.config` per executar-se al nostre fitxer `greetings.csv` amb alguns parametres per defecte.

#### 4.3.1. Copiar el fitxer `greetings.csv`

Primer necessitem copiar el fitxer `greetings.csv` a un lloc apropiat al nostre projecte de pipeline.
Normalment els fitxers de prova petits s'emmagatzemen al directori `assets`, així que copiem el fitxer des del nostre directori de treball.

```bash
cp greetings.csv core-hello/assets/.
```

Ara el fitxer `greetings.csv` està llest per ser utilitzat com a entrada de prova.

#### 4.3.2. Actualitzar el fitxer `test.config`

Ara podem actualitzar el fitxer `test.config` de la següent manera:

=== "Després"

    ```groovy title="core-hello/conf/test.config" linenums="21" hl_lines="6-10"
    params {
        config_profile_name        = 'Test profile'
        config_profile_description = 'Minimal test dataset to check pipeline function'

        // Input data
        input  = "${projectDir}/assets/greetings.csv"

        // Other parameters
        batch     = 'test'
        character = 'tux'
    }
    ```

=== "Abans"

    ```groovy title="core-hello/conf/test.config" linenums="21" hl_lines="6-8"
    params {
        config_profile_name        = 'Test profile'
        config_profile_description = 'Minimal test dataset to check pipeline function'

        // Input data
        // TODO nf-core: Specify the paths to your test data on nf-core/test-datasets
        // TODO nf-core: Give any required params for the test so that command line flags are not needed
        input  = params.pipelines_testdata_base_path + 'viralrecon/samplesheet/samplesheet_test_illumina_amplicon.csv'
    }
    ```

Punts clau:

- **Utilitzant `${projectDir}`**: Aquesta és una variable implícita de Nextflow que apunta al directori on es troba l'script de workflow principal (l'arrel del pipeline). Utilitzar-la assegura que el camí funcioni independentment d'on s'executi el pipeline.
- **Camins absoluts**: Utilitzant `${projectDir}`, creem un camí absolut, que és important per a dades de prova que s'envien amb el pipeline.
- **Ubicació de dades de prova**: Els pipelines nf-core normalment emmagatzemen dades de prova al directori `assets/` dins del repositori del pipeline per a fitxers de prova petits, o fan referència a conjunts de dades de prova externs per a fitxers més grans.

I mentre hi som, ajustem els límits de recursos per defecte per assegurar que això s'executarà en màquines molt bàsiques (com les VMs mínimes a Github Codespaces):

=== "Després"

    ```groovy title="core-hello/config/test.config" linenums="13" hl_lines="3-4"
    process {
        resourceLimits = [
            cpus: 2,
            memory: '4.GB',
            time: '1.h'
        ]
    }
    ```

=== "Abans"

    ```groovy title="core-hello/config/test.config" linenums="13" hl_lines="3-4"
    process {
        resourceLimits = [
            cpus: 4,
            memory: '15.GB',
            time: '1.h'
        ]
    }
    ```

Això completa les modificacions de codi que necessitem fer.

### 4.4. Executar el pipeline amb el perfil de test

Això ha estat molt, però finalment podem provar d'executar el pipeline!
Tingueu en compte que hem d'afegir `--validate_params false` a la línia de comandes perquè encara no hem configurat la validació (això vindrà més endavant).

```bash
nextflow run core-hello --outdir core-hello-results -profile test,docker --validate_params false
```

Si heu fet totes les modificacions correctament, hauria de completar-se.

??? success "Sortida de la comanda"

    ```console
     N E X T F L O W   ~  version 25.04.3

    Launching `core-hello/main.nf` [condescending_allen] DSL2 - revision: b9e9b3b8de

    Input/output options
      input                     : /workspaces/training/hello-nf-core/core-hello/assets/greetings.csv
      outdir                    : core-hello-results

    Institutional config options
      config_profile_name       : Test profile
      config_profile_description: Minimal test dataset to check pipeline function

    Generic options
      validate_params           : false
      trace_report_suffix       : 2025-11-21_07-29-37

    Core Nextflow options
      runName                   : condescending_allen
      containerEngine           : docker
      launchDir                 : /workspaces/training/hello-nf-core
      workDir                   : /workspaces/training/hello-nf-core/work
      projectDir                : /workspaces/training/hello-nf-core/core-hello
      userName                  : root
      profile                   : test,docker
      configFiles               : /workspaces/training/hello-nf-core/core-hello/nextflow.config

    !! Only displaying parameters that differ from the pipeline defaults !!
    ------------------------------------------------------
    executor >  local (1)
    [ed/727b7e] CORE_HELLO:HELLO:sayHello (3)       [100%] 3 of 3 ✔
    [45/bb6096] CORE_HELLO:HELLO:convertToUpper (3) [100%] 3 of 3 ✔
    [81/7e2e34] CORE_HELLO:HELLO:collectGreetings   [100%] 1 of 1 ✔
    [96/9442a1] CORE_HELLO:HELLO:cowpy              [100%] 1 of 1 ✔
    -[core/hello] Pipeline completed successfully-
    ```

Com podeu veure, això va produir el resum típic d'nf-core a l'inici gràcies al subworkflow d'inicialització, i les línies per a cada mòdul ara mostren els noms complets PIPELINE:WORKFLOW:module.

### 4.5. Trobar les sortides del pipeline

La pregunta ara és: on són les sortides del pipeline?
I la resposta és força interessant: ara hi ha dos llocs diferents per buscar els resultats.

Com potser recordeu d'abans, la nostra primera execució del workflow recentment creat va produir un directori anomenat `core-hello-results/` que contenia diversos informes d'execució i metadades.

```bash
tree core-hello-results
```

??? abstract "Contingut del directori"

    ```console
    core-hello-results
    └── pipeline_info
        ├── execution_report_2025-11-21_04-47-18.html
        ├── execution_report_2025-11-21_07-29-37.html
        ├── execution_timeline_2025-11-21_04-47-18.html
        ├── execution_timeline_2025-11-21_07-29-37.html
        ├── execution_trace_2025-11-21_04-47-18.txt
        ├── execution_trace_2025-11-21_07-29-37.txt
        ├── hello_software_versions.yml
        ├── params_2025-11-21_04-47-13.json
        ├── params_2025-11-21_07-29-41.json
        └── pipeline_dag_2025-11-21_04-47-18.html
        └── pipeline_dag_2025-11-21_07-29-37.html

    1 directory, 12 files
    ```

Veieu que vam obtenir un altre conjunt d'informes d'execució a més dels que vam obtenir de la primera execució, quan el workflow encara era només un marcador de posició.
Aquesta vegada veieu totes les tasques que es van executar com s'esperava.

![informe de línia de temps d'execució per al pipeline Hello](./img/execution_timeline_hello.png)

!!! note "Nota"

    Una vegada més les tasques no es van executar en paral·lel perquè estem executant en una màquina minimalista a Github Codespaces.
    Per veure-les executar-se en paral·lel, proveu d'augmentar l'assignació de CPU del vostre codespace i els límits de recursos a la configuració de test.

Això és genial, però els nostres resultats reals del pipeline no hi són!

Això és el que va passar: no vam canviar res als mòduls en si, així que les sortides gestionades per les directives `publishDir` a nivell de mòdul encara van a un directori `results` tal com s'especifica al pipeline original.

```bash
tree results
```

??? abstract "Contingut del directori"

    ```console
    results
    ├── Bonjour-output.txt
    ├── COLLECTED-test-batch-output.txt
    ├── COLLECTED-test-output.txt
    ├── cowpy-COLLECTED-test-batch-output.txt
    ├── cowpy-COLLECTED-test-output.txt
    ├── Hello-output.txt
    ├── Holà-output.txt
    ├── UPPER-Bonjour-output.txt
    ├── UPPER-Hello-output.txt
    └── UPPER-Holà-output.txt

    0 directories, 10 files
    ```

Ah, aquí estan, barrejats amb les sortides d'execucions anteriors del pipeline Hello original.

Si volem que estiguin organitzats de manera ordenada com ho estaven les sortides del pipeline demo, haurem de canviar com configurem les sortides per ser publicades.
Us mostrarem com fer-ho més endavant en aquest curs de formació.

<!-- TODO: Update this once we've updated Hello Nextflow to use workflow-level outputs -->

I aquí està! Pot semblar molta feina per aconseguir el mateix resultat que el pipeline original, però obteniu tots aquests informes encantadors generats automàticament, i ara teniu una base sòlida per aprofitar funcionalitats addicionals d'nf-core, incloent validació d'entrades i algunes capacitats de maneig de metadades interessants que cobrirem en una secció posterior.

---

### Conclusió

Sabeu com convertir un pipeline Nextflow regular en un pipeline d'estil nf-core utilitzant la plantilla nf-core.
Com a part d'això, heu après com fer un workflow composable, i com identificar els elements de la plantilla nf-core que més comunament necessiten ser adaptats quan es desenvolupa un pipeline d'estil nf-core personalitzat.

### Què segueix?

Feu una pausa, ha estat una feina dura! Quan estigueu preparats, passeu a [Part 3: Use an nf-core module](./03_use_module.md) per aprendre com aprofitar mòduls mantinguts per la comunitat del repositori nf-core/modules.
