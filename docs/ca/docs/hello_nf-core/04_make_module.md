# Part 4: Crear un mòdul nf-core

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducció assistida per IA - [més informació i suggeriments](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

En aquesta quarta part del curs de formació Hello nf-core, us mostrem com crear un mòdul nf-core aplicant les convencions clau que fan que els mòduls siguin portables i mantenibles.

El projecte nf-core proporciona una comanda (`nf-core modules create`) que genera plantilles de mòduls estructurades correctament de manera automàtica, similar al que vam utilitzar per al workflow a la Part 2.
No obstant això, amb finalitats didàctiques, començarem fent-ho manualment: transformant el mòdul local `cowpy` del vostre pipeline `core-hello` en un mòdul d'estil nf-core pas a pas.
Després d'això, us mostrarem com utilitzar la creació de mòduls basada en plantilles per treballar de manera més eficient en el futur.

??? info "Com començar des d'aquesta secció"

    Aquesta secció assumeix que heu completat la [Part 3: Utilitzar un mòdul nf-core](./03_use_module.md) i heu integrat el mòdul `CAT_CAT` al vostre pipeline.

    Si no vau completar la Part 3 o voleu començar de nou per a aquesta part, podeu utilitzar la solució `core-hello-part3` com a punt de partida.
    Executeu aquestes comandes des de dins del directori `hello-nf-core/`:

    ```bash
    cp -r solutions/core-hello-part3 core-hello
    cd core-hello
    ```

    Això us proporciona un pipeline amb el mòdul `CAT_CAT` ja integrat.
    Podeu comprovar que s'executa correctament executant la comanda següent:

    ```bash
    nextflow run . --outdir core-hello-results -profile test,docker --validate_params false
    ```

---

## 1. Transformar `cowpy` en un mòdul nf-core

En aquesta secció, aplicarem les convencions nf-core al mòdul local `cowpy` del vostre pipeline `core-hello`, transformant-lo en un mòdul que segueix els estàndards de la comunitat nf-core.

Aquest és el codi actual del mòdul de procés `cowpy`:

```groovy title="core-hello/modules/local/cowpy.nf" linenums="1"
#!/usr/bin/env nextflow

// Generate ASCII art with cowpy (https://github.com/jeffbuttars/cowpy)
process cowpy {

    publishDir 'results', mode: 'copy'

    container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
    conda 'conda-forge::cowpy==1.1.5'

    input:
    path input_file
    val character

    output:
    path "cowpy-${input_file}"

    script:
    """
    cat $input_file | cowpy -c "$character" > cowpy-${input_file}
    """
}
```

Aplicarem les següents convencions nf-core de manera incremental:

1. **Posar el nom del procés en majúscules a `COWPY`** per seguir la convenció.
2. **Actualitzar `COWPY` per utilitzar tuples de metadades** per propagar les metadades de mostra a través del workflow.
3. **Centralitzar la configuració d'arguments de l'eina amb `ext.args`** per augmentar la versatilitat del mòdul mantenint la interfície mínima.
4. **Estandarditzar el nomenament de sortides amb `ext.prefix`** per promoure la consistència.
5. **Centralitzar la configuració de publicació** per promoure la consistència.

Després de cada pas, executarem el pipeline per comprovar que tot funciona com s'espera.

!!! warning "Directori de treball"

    Assegureu-vos que esteu al directori `core-hello` (l'arrel del vostre pipeline) per a totes les edicions de fitxers i execucions de comandes d'aquesta secció.

    ```bash
    cd core-hello
    ```

### 1.1. Posar el nom del procés en majúscules

Això és purament una convenció estilística (no hi ha justificació tècnica) però com que és la norma per als mòduls nf-core, complim-la.

Hem de fer tres conjunts de canvis:

1. Actualitzar el nom del procés al mòdul
2. Actualitzar la declaració d'importació del mòdul a la capçalera del workflow
3. Actualitzar la crida al procés i la declaració emit al cos del workflow

Comencem!

#### 1.1.1. Actualitzar el nom del procés al mòdul

Obriu el fitxer del mòdul `cowpy.nf` (sota `core-hello/modules/local/`) i modifiqueu el nom del procés a majúscules:

=== "Després"

    ```groovy title="core-hello/modules/local/cowpy.nf" linenums="3" hl_lines="2"
    // Generate ASCII art with cowpy (https://github.com/jeffbuttars/cowpy)
    process COWPY {
    ```

=== "Abans"

    ```groovy title="core-hello/modules/local/cowpy.nf" linenums="3" hl_lines="2"
    // Generate ASCII art with cowpy (https://github.com/jeffbuttars/cowpy)
    process cowpy {
    ```

En aquest cas, posar-ho en majúscules és completament directe.

Si el nom del procés estigués compost per diverses paraules, per exemple si tinguéssim un procés anomenat MyCowpyTool originalment en camel case, la convenció nf-core seria utilitzar guions baixos per separar-les, donant MY_COWPY_TOOL.

#### 1.1.2. Actualitzar la declaració d'importació del mòdul

Els noms de procés distingeixen entre majúscules i minúscules, així que ara que hem canviat el nom del procés, hem d'actualitzar la declaració d'importació del mòdul en conseqüència a la capçalera del workflow de `hello.nf`:

=== "Després"

    ```groovy title="core-hello/workflows/hello.nf" linenums="1" hl_lines="10"
    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    include { paramsSummaryMap       } from 'plugin/nf-schema'
    include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
    include { sayHello               } from '../modules/local/sayHello.nf'
    include { convertToUpper         } from '../modules/local/convertToUpper.nf'
    include { COWPY                  } from '../modules/local/cowpy.nf'
    include { CAT_CAT                } from '../modules/nf-core/cat/cat/main'
    ```

=== "Abans"

    ```groovy title="core-hello/workflows/hello.nf" linenums="1" hl_lines="10"
    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    include { paramsSummaryMap       } from 'plugin/nf-schema'
    include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
    include { sayHello               } from '../modules/local/sayHello.nf'
    include { convertToUpper         } from '../modules/local/convertToUpper.nf'
    include { cowpy                  } from '../modules/local/cowpy.nf'
    include { CAT_CAT                } from '../modules/nf-core/cat/cat/main'
    ```

Podríem utilitzar un àlies a la declaració d'importació per evitar haver d'actualitzar les crides al procés, però això aniria en contra del propòsit d'adoptar la convenció de majúscules.

#### 1.1.3. Actualitzar la crida al procés i la declaració emit

Així que ara actualitzem les dues referències al procés al bloc workflow de `hello.nf`:

=== "Després"

    ```groovy title="core-hello/workflows/hello.nf" linenums="43" hl_lines="5 20"
    // extract the file from the tuple since cowpy doesn't use metadata yet
    ch_for_cowpy = CAT_CAT.out.file_out.map{ meta, file -> file }

    // generate ASCII art of the greetings with cowpy
    COWPY(ch_for_cowpy, params.character)

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
    cowpy_hellos   = COWPY.out
    versions       = ch_versions
    ```

=== "Abans"

    ```groovy title="core-hello/workflows/hello.nf" linenums="43" hl_lines="5 20"
    // extract the file from the tuple since cowpy doesn't use metadata yet
    ch_for_cowpy = CAT_CAT.out.file_out.map{ meta, file -> file }

    // generate ASCII art of the greetings with cowpy
    cowpy(ch_for_cowpy, params.character)

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
    cowpy_hellos   = cowpy.out
    versions       = ch_versions
    ```

Assegureu-vos de fer **ambdós** canvis, altrament obtindreu un error quan executeu això.

#### 1.1.4. Executar el pipeline per provar-lo

Executem el workflow per comprovar que tot funciona correctament després d'aquests canvis.

```bash
nextflow run . --outdir core-hello-results -profile test,docker --validate_params false
```

??? success "Sortida de la comanda"

    ```console
     N E X T F L O W   ~  version 25.04.3

    Launching `./main.nf` [elegant_plateau] DSL2 - revision: b9e9b3b8de

    Input/output options
      input                     : /workspaces/training/hello-nf-core/core-hello/assets/greetings.csv
      outdir                    : core-hello-results

    Institutional config options
      config_profile_name       : Test profile
      config_profile_description: Minimal test dataset to check pipeline function

    Generic options
      validate_params           : false
      trace_report_suffix       : 2026-01-06_04-51-29

    Core Nextflow options
      runName                   : elegant_plateau
      containerEngine           : docker
      launchDir                 : /workspaces/training/hello-nf-core/core-hello
      workDir                   : /workspaces/training/hello-nf-core/core-hello/work
      projectDir                : /workspaces/training/hello-nf-core/core-hello
      userName                  : root
      profile                   : test,docker
      configFiles               : /workspaces/training/hello-nf-core/core-hello/nextflow.config

    !! Only displaying parameters that differ from the pipeline defaults !!
    ------------------------------------------------------
    executor >  local (8)
    [7b/66ceb5] CORE_HELLO:HELLO:sayHello (3)       | 3 of 3 ✔
    [8e/1bafb9] CORE_HELLO:HELLO:convertToUpper (3) | 3 of 3 ✔
    [bb/203575] CORE_HELLO:HELLO:CAT_CAT (test)     | 1 of 1 ✔
    [39/715489] CORE_HELLO:HELLO:COWPY              | 1 of 1 ✔
    -[core/hello] Pipeline completed successfully-
    ```

D'acord, això funciona! Ara passem a fer canvis més substancials.

### 1.2. Actualitzar `COWPY` per utilitzar tuples de metadades

A la versió actual del pipeline `core-hello`, estem extraient el fitxer de la tupla de sortida de `CAT_CAT` per passar-lo a `COWPY`, com es mostra a la meitat superior del diagrama següent.

<figure class="excalidraw">
    --8<-- "docs/en/docs/hello_nf-core/img/cowpy-inputs.svg"
</figure>

Seria millor que `COWPY` acceptés tuples de metadades directament, permetent que les metadades flueixin a través del workflow, com es mostra a la meitat inferior del diagrama.

Per aconseguir-ho, haurem de fer els següents canvis:

1. Actualitzar les definicions d'entrada i sortida
2. Actualitzar la crida al procés al workflow
3. Actualitzar el bloc emit al workflow

Un cop hàgim fet tot això, executarem el pipeline per comprovar que tot encara funciona com abans.

#### 1.2.1. Actualitzar les definicions d'entrada i sortida

Torneu al fitxer del mòdul `cowpy.nf` i modifiqueu-lo per acceptar tuples de metadades com es mostra a continuació.

=== "Després"

    ```groovy title="core-hello/modules/local/cowpy.nf" linenums="11" hl_lines="2 6"
        input:
            tuple val(meta), path(input_file)
            val character

        output:
            tuple val(meta), path("cowpy-${input_file}"), emit: cowpy_output
    ```

=== "Abans"

    ```groovy title="core-hello/modules/local/cowpy.nf" linenums="11" hl_lines="2 6"
        input:
            path input_file
            val character

        output:
            path "cowpy-${input_file}"
    ```

Com podeu veure, hem canviat tant l'**entrada principal** com la **sortida** a una tupla que segueix el patró `tuple val(meta), path(input_file)` introduït a la Part 3 d'aquesta formació.
Per a la sortida, també hem aprofitat per afegir `emit: cowpy_output` per donar al canal de sortida un nom descriptiu.

Ara que hem canviat el que el procés espera, hem d'actualitzar el que li proporcionem a la crida del procés.

#### 1.2.2. Actualitzar la crida al procés al workflow

La bona notícia és que aquest canvi simplificarà la crida al procés.
Ara que la sortida de `CAT_CAT` i l'entrada de `COWPY` tenen la mateixa 'forma', és a dir, ambdues consisteixen en una estructura `tuple val(meta), path(input_file)`, podem simplement connectar-les directament en lloc d'haver d'extreure el fitxer explícitament de la sortida del procés `CAT_CAT`.

Obriu el fitxer del workflow `hello.nf` (sota `core-hello/workflows/`) i actualitzeu la crida a `COWPY` com es mostra a continuació.

=== "Després"

    ```groovy title="core-hello/workflows/hello.nf" linenums="43" hl_lines="2"
        // generate ASCII art of the greetings with cowpy
        COWPY(CAT_CAT.out.file_out, params.character)
    ```

=== "Abans"

    ```groovy title="core-hello/workflows/hello.nf" linenums="43" hl_lines="1-2 5"
        // extract the file from the tuple since cowpy doesn't use metadata yet
        ch_for_cowpy = CAT_CAT.out.file_out.map{ meta, file -> file }

        // generate ASCII art of the greetings with cowpy
        COWPY(ch_for_cowpy, params.character)
    ```

Ara cridem `COWPY` sobre `CAT_CAT.out.file_out` directament.

Com a resultat, ja no necessitem construir el canal `ch_for_cowpy`, així que aquesta línia (i la seva línia de comentari) es pot eliminar completament.

#### 1.2.3. Actualitzar el bloc emit al workflow

Com que `COWPY` ara emet una sortida amb nom, `cowpy_output`, podem actualitzar el bloc `emit:` del workflow `hello.nf` per utilitzar-la.

=== "Després"

    ```groovy title="core-hello/workflows/hello.nf" linenums="60" hl_lines="2"
        emit:
        cowpy_hellos   = COWPY.out.cowpy_output
        versions       = ch_versions
    ```

=== "Abans"

    ```groovy title="core-hello/workflows/hello.nf" linenums="60" hl_lines="2"
        emit:
        cowpy_hellos   = COWPY.out
        versions       = ch_versions
    ```

Això tècnicament no és necessari, però és una bona pràctica referir-se a sortides amb nom sempre que sigui possible.

#### 1.2.4. Executar el pipeline per provar-lo

Executem el workflow per comprovar que tot funciona correctament després d'aquests canvis.

```bash
nextflow run . --outdir core-hello-results -profile test,docker --validate_params false
```

??? success "Sortida de la comanda"

    ```console
     N E X T F L O W   ~  version 25.04.3

    Launching `./main.nf` [modest_saha] DSL2 - revision: b9e9b3b8de

    Downloading plugin nf-schema@2.5.1
    Input/output options
      input                     : /workspaces/training/hello-nf-core/core-hello/assets/greetings.csv
      outdir                    : core-hello-results

    Institutional config options
      config_profile_name       : Test profile
      config_profile_description: Minimal test dataset to check pipeline function

    Generic options
      validate_params           : false
      trace_report_suffix       : 2025-12-27_06-16-55

    Core Nextflow options
      runName                   : modest_saha
      containerEngine           : docker
      launchDir                 : /workspaces/training/hello-nf-core/core-hello
      workDir                   : /workspaces/training/hello-nf-core/core-hello/work
      projectDir                : /workspaces/training/hello-nf-core/core-hello
      userName                  : root
      profile                   : test,docker
      configFiles               : /workspaces/training/hello-nf-core/core-hello/nextflow.config

    !! Only displaying parameters that differ from the pipeline defaults !!
    ------------------------------------------------------
    executor >  local (8)
    [a8/447993] CORE_HELLO:HELLO:sayHello (3)       [100%] 3 of 3 ✔
    [00/1fc59c] CORE_HELLO:HELLO:convertToUpper (3) [100%] 3 of 3 ✔
    [57/ac800d] CORE_HELLO:HELLO:CAT_CAT (test)     [100%] 1 of 1 ✔
    [b7/092f2b] CORE_HELLO:HELLO:COWPY              [100%] 1 of 1 ✔
    -[core/hello] Pipeline completed successfully-
    ```

El pipeline hauria d'executar-se correctament, amb les metadades fluint ara de `CAT_CAT` a través de `COWPY`.

Això completa el que necessitàvem fer perquè `COWPY` gestioni tuples de metadades.
Ara, vegem què més podem fer per aprofitar els patrons de mòduls nf-core.

### 1.3. Centralitzar la configuració d'arguments de l'eina amb `ext.args`

En el seu estat actual, el procés `COWPY` espera rebre un valor per al paràmetre `character`.
Com a resultat, hem de proporcionar un valor cada vegada que cridem el procés, fins i tot si estaríem contents amb els valors per defecte establerts per l'eina.
Per a `COWPY` això no és un gran problema, però per a eines amb molts paràmetres opcionals, pot ser força feixuc.

El projecte nf-core recomana utilitzar una funcionalitat de Nextflow anomenada [`ext.args`](https://www.nextflow.io/docs/latest/reference/process.html#ext) per gestionar els arguments de l'eina de manera més convenient mitjançant fitxers de configuració.

En lloc de declarar entrades de procés per a cada opció de l'eina, escriviu el mòdul per referenciar `ext.args` en la construcció de la seva línia de comandes.
Llavors només cal configurar la variable `ext.args` per contenir els arguments i valors que voleu utilitzar al fitxer `modules.config`, que consolida els detalls de configuració per a tots els mòduls.
Nextflow afegirà aquests arguments amb els seus valors a la línia de comandes de l'eina en temps d'execució.

Apliquem aquest enfocament al mòdul `COWPY`.
Haurem de fer els següents canvis:

1. Actualitzar el mòdul `COWPY`
2. Configurar `ext.args` al fitxer `modules.config`
3. Actualitzar el workflow `hello.nf`

Un cop hàgim fet tot això, executarem el pipeline per comprovar que tot encara funciona com abans.

#### 1.3.1. Actualitzar el mòdul `COWPY`

Fem-ho.
Obriu el fitxer del mòdul `cowpy.nf` (sota `core-hello/modules/local/`) i modifiqueu-lo per referenciar `ext.args` com es mostra a continuació.

=== "Després"

    ```groovy title="modules/local/cowpy.nf" linenums="1" hl_lines="18 20"
    #!/usr/bin/env nextflow

    // Generate ASCII art with cowpy (https://github.com/jeffbuttars/cowpy)
    process COWPY {

        publishDir 'results', mode: 'copy'

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
        conda 'conda-forge::cowpy==1.1.5'

        input:
            tuple val(meta), path(input_file)

        output:
            tuple val(meta), path("cowpy-${input_file}"), emit: cowpy_output

        script:
        def args = task.ext.args ?: ''
        """
        cat $input_file | cowpy $args > cowpy-${input_file}
        """
    }
    ```

=== "Abans"

    ```groovy title="core-hello/modules/local/cowpy.nf" linenums="1" hl_lines="13 20"
    #!/usr/bin/env nextflow

    // Generate ASCII art with cowpy (https://github.com/jeffbuttars/cowpy)
    process COWPY {

        publishDir 'results', mode: 'copy'

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
        conda 'conda-forge::cowpy==1.1.5'

        input:
            tuple val(meta), path(input_file)
            val character

        output:
            tuple val(meta), path("cowpy-${input_file}"), emit: cowpy_output

        script:
        """
        cat $input_file | cowpy -c "$character" > cowpy-${input_file}
        """
    }
    ```

Podeu veure que hem fet tres canvis.

1. **Al bloc `input:`, hem eliminat l'entrada `val character`.**
   D'ara endavant, proporcionarem aquest argument mitjançant la configuració `ext.args` com es descriu més avall.

2. **Al bloc `script:`, hem afegit la línia `def args = task.ext.args ?: ''`.**
   Aquesta línia utilitza l'operador `?:` per determinar el valor de la variable `args`: el contingut de `task.ext.args` si no està buit, o una cadena buida si ho està.
   Tingueu en compte que tot i que generalment ens referim a `ext.args`, aquest codi ha de referenciar `task.ext.args` per extreure la configuració `ext.args` a nivell de mòdul.

3. **A la línia de comandes, hem reemplaçat `-c "$character"` amb `$args`.**
   Aquí és on Nextflow injectarà qualsevol argument d'eina establert a `ext.args` al fitxer `modules.config`.

Com a resultat, la interfície del mòdul ara és més simple: només espera les entrades essencials de metadades i fitxers.

!!! note "Nota"

    L'operador `?:` sovint s'anomena 'operador Elvis' perquè sembla una cara d'Elvis Presley de costat, amb el caràcter `?` simbolitzant l'ona del seu cabell.

#### 1.3.2. Configurar `ext.args` al fitxer `modules.config`

Ara que hem tret la declaració de `character` del mòdul, hem d'afegir-la a `ext.args` al fitxer de configuració `modules.config`.

Específicament, afegirem aquest petit fragment de codi al bloc `process {}`:

```groovy title="Code to add"
withName: 'COWPY' {
    ext.args = { "-c ${params.character}" }
}
```

La sintaxi `withName:` assigna aquesta configuració només al procés `COWPY`, i `ext.args = { "-c ${params.character}" }` simplement compon una cadena que inclourà el valor del paràmetre `character`.
Tingueu en compte l'ús de claus, que indiquen a Nextflow que avaluï el valor del paràmetre en temps d'execució.

Té sentit? Afegim-ho.

Obriu `conf/modules.config` i afegiu el codi de configuració dins del bloc `process {}` com es mostra a continuació.

=== "Després"

    ```groovy title="core-hello/conf/modules.config" linenums="13" hl_lines="8-10"
    process {
        publishDir = [
            path: { "${params.outdir}/${task.process.tokenize(':')[-1].tokenize('_')[0].toLowerCase()}" },
            mode: params.publish_dir_mode,
            saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
        ]

        withName: 'COWPY' {
            ext.args = { "-c ${params.character}" }
        }
    }
    ```

=== "Abans"

    ```groovy title="core-hello/conf/modules.config" linenums="13"
    process {
        publishDir = [
            path: { "${params.outdir}/${task.process.tokenize(':')[-1].tokenize('_')[0].toLowerCase()}" },
            mode: params.publish_dir_mode,
            saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
        ]
    }
    ```

Esperem que pugueu imaginar tenir tots els mòduls d'un pipeline amb els seus `ext.args` especificats en aquest fitxer, amb els següents beneficis:

- La **interfície del mòdul es manté simple** - Només accepta les entrades essencials de metadades i fitxers
- El **pipeline encara exposa `params.character`** - Els usuaris finals encara poden configurar-lo com abans
- El **mòdul ara és portable** - Es pot reutilitzar en altres pipelines sense esperar un nom de paràmetre específic
- La configuració està **centralitzada** a `modules.config`, mantenint la lògica del workflow neta

Utilitzant el fitxer `modules.config` com el lloc on tots els pipelines centralitzen la configuració per mòdul, fem que els nostres mòduls siguin més reutilitzables entre diferents pipelines.

#### 1.3.3. Actualitzar el workflow `hello.nf`

Com que el mòdul `COWPY` ja no requereix el paràmetre `character` com a entrada, hem d'actualitzar la crida del workflow en conseqüència.

Obriu el fitxer del workflow `hello.nf` (sota `core-hello/workflows/`) i actualitzeu la crida a `COWPY` com es mostra a continuació.

=== "Després"

    ```groovy title="core-hello/workflows/hello.nf" linenums="39" hl_lines="2"
        // generate ASCII art of the greetings with cowpy
        COWPY(CAT_CAT.out.file_out)
    ```

=== "Abans"

    ```groovy title="core-hello/workflows/hello.nf" linenums="39" hl_lines="2"
        // generate ASCII art of the greetings with cowpy
        COWPY(CAT_CAT.out.file_out, params.character)
    ```

El codi del workflow ara és més net: no necessitem passar `params.character` directament al procés.
La interfície del mòdul es manté mínima, fent-la més portable, mentre que el pipeline encara proporciona l'opció explícita mitjançant configuració.

#### 1.3.4. Executar el pipeline per provar-lo

Comprovem que el workflow encara funciona com s'espera, especificant un caràcter diferent per verificar que la configuració `ext.args` està funcionant.

Executeu aquesta comanda utilitzant `kosh`, una de les opcions més... enigmàtiques:

```bash
nextflow run . --outdir core-hello-results -profile test,docker --validate_params false --character kosh
```

??? success "Sortida de la comanda"

    ```console
     N E X T F L O W   ~  version 25.04.3

    Launching `./main.nf` [exotic_planck] DSL2 - revision: b9e9b3b8de

    Input/output options
      input                     : /workspaces/training/hello-nf-core/core-hello/assets/greetings.csv
      outdir                    : core-hello-results

    Institutional config options
      config_profile_name       : Test profile
      config_profile_description: Minimal test dataset to check pipeline function

    Generic options
      validate_params           : false
      trace_report_suffix       : 2025-12-27_06-23-13

    Core Nextflow options
      runName                   : exotic_planck
      containerEngine           : docker
      launchDir                 : /workspaces/training/hello-nf-core/core-hello
      workDir                   : /workspaces/training/hello-nf-core/core-hello/work
      projectDir                : /workspaces/training/hello-nf-core/core-hello
      userName                  : root
      profile                   : test,docker
      configFiles               : /workspaces/training/hello-nf-core/core-hello/nextflow.config

    !! Only displaying parameters that differ from the pipeline defaults !!
    ------------------------------------------------------
    executor >  local (8)
    [13/9e3c0e] CORE_HELLO:HELLO:sayHello (3)       [100%] 3 of 3 ✔
    [e2/5b0ee5] CORE_HELLO:HELLO:convertToUpper (3) [100%] 3 of 3 ✔
    [b6/4fb569] CORE_HELLO:HELLO:CAT_CAT (test)     [100%] 1 of 1 ✔
    [38/eb29ea] CORE_HELLO:HELLO:COWPY              [100%] 1 of 1 ✔
    -[core/hello] Pipeline completed successfully-
    ```

Això hauria d'executar-se correctament com abans.

Verifiquem que la configuració `ext.args` ha funcionat comprovant la sortida.
Trobeu la sortida al navegador de fitxers o utilitzeu el hash de la tasca (la part `38/eb29ea` a l'exemple anterior) per mirar el fitxer de sortida:

```bash
cat work/38/eb29ea*/cowpy-test.txt
```

??? success "Sortida de la comanda"

    ```console
    _________
    / HELLO   \
    | HOLà    |
    \ BONJOUR /
    ---------
        \
        \
          \
      ___       _____     ___
    /   \     /    /|   /   \
    |     |   /    / |  |     |
    |     |  /____/  |  |     |
    |     |  |    |  |  |     |
    |     |  | {} | /   |     |
    |     |  |____|/    |     |
    |     |    |==|     |     |
    |      \___________/      |
    |                         |
    |                         |
    ```

Hauríeu de veure l'art ASCII mostrat amb el caràcter `kosh`, confirmant que la configuració `ext.args` ha funcionat!

??? info "(Opcional) Inspeccionar el fitxer de comandes"

    Si voleu veure exactament com s'ha aplicat la configuració, podeu inspeccionar el fitxer `.command.sh`:

    ```bash
    cat work/38/eb29ea*/.command.sh
    ```

    Veureu la comanda `cowpy` amb l'argument `-c kosh`:

    ```console
    #!/usr/bin/env bash
    ...
    cat test.txt | cowpy -c kosh > cowpy-test.txt
    ```

    Això mostra que el fitxer `.command.sh` s'ha generat correctament basant-se en la configuració `ext.args`.

Preneu-vos un moment per pensar en el que hem aconseguit aquí.
Aquest enfocament manté la interfície del mòdul centrada en dades essencials (fitxers, metadades i qualsevol paràmetre obligatori per mostra), mentre que les opcions que controlen el comportament de l'eina es gestionen per separat mitjançant configuració.

Això pot semblar innecessari per a una eina simple com `cowpy`, però pot fer una gran diferència per a eines d'anàlisi de dades que tenen molts arguments opcionals.

Per resumir els beneficis d'aquest enfocament:

- **Interfície neta**: El mòdul se centra en entrades de dades essencials (metadades i fitxers)
- **Flexibilitat**: Els usuaris poden especificar arguments d'eina mitjançant configuració, incloent valors específics per mostra
- **Consistència**: Tots els mòduls nf-core segueixen aquest patró
- **Portabilitat**: Els mòduls es poden reutilitzar sense opcions d'eina codificades
- **Sense canvis al workflow**: Afegir o canviar opcions d'eina no requereix actualitzar el codi del workflow

!!! note "Nota"

    El sistema `ext.args` té capacitats addicionals potents no cobertes aquí, incloent canviar valors d'arguments dinàmicament basant-se en metadades. Consulteu les [especificacions de mòduls nf-core](https://nf-co.re/docs/guidelines/components/modules) per a més detalls.

### 1.4. Estandarditzar el nomenament de sortides amb `ext.prefix`

Ara que hem donat al procés `COWPY` accés al metamap, podem començar a aprofitar un altre patró útil d'nf-core: nomenar fitxers de sortida basant-se en metadades.

Aquí utilitzarem una funcionalitat de Nextflow anomenada `ext.prefix` que ens permetrà estandarditzar el nomenament de fitxers de sortida entre mòduls utilitzant `meta.id` (l'identificador inclòs al metamap), mentre encara podem configurar mòduls individualment si es desitja.

Això serà similar al que vam fer amb `ext.args`, amb algunes diferències que detallarem a mesura que avancem.

Apliquem aquest enfocament al mòdul `COWPY`.
Haurem de fer els següents canvis:

1. Actualitzar el mòdul `COWPY`
2. Configurar `ext.prefix` al fitxer `modules.config`

(No calen canvis al workflow.)

Un cop hàgim fet això, executarem el pipeline per comprovar que tot encara funciona com abans.

#### 1.4.1. Actualitzar el mòdul `COWPY`

Obriu el fitxer del mòdul `cowpy.nf` (sota `core-hello/modules/local/`) i modifiqueu-lo per referenciar `ext.prefix` com es mostra a continuació.

=== "Després"

    ```groovy title="core-hello/modules/local/cowpy.nf" linenums="1" hl_lines="2 6 8"
        output:
            tuple val(meta), path("${prefix}.txt"), emit: cowpy_output

        script:
        def args = task.ext.args ?: ''
        prefix = task.ext.prefix ?: "${meta.id}"
        """
        cat $input_file | cowpy $args > ${prefix}.txt
        """
    }
    ```

=== "Abans"

    ```groovy title="core-hello/modules/local/cowpy.nf" linenums="1" hl_lines="2 7"
        output:
            tuple val(meta), path("cowpy-${input_file}"), emit: cowpy_output

        script:
        def args = task.ext.args ?: ''
        """
        cat $input_file | cowpy $args > cowpy-${input_file}
        """
    }
    ```

Podeu veure que hem fet tres canvis.

1. **Al bloc `script:`, hem afegit la línia `prefix = task.ext.prefix ?: "${meta.id}"`.**
   Aquesta línia utilitza l'operador `?:` per determinar el valor de la variable `prefix`: el contingut de `task.ext.prefix` si no està buit, o l'identificador del metamap (`meta.id`) si ho està.
   Tingueu en compte que tot i que generalment ens referim a `ext.prefix`, aquest codi ha de referenciar `task.ext.prefix` per extreure la configuració `ext.prefix` a nivell de mòdul.

2. **A la línia de comandes, hem reemplaçat `cowpy-${input_file}` amb `${prefix}.txt`.**
   Aquí és on Nextflow injectarà el valor de `prefix` determinat per la línia anterior.

3. **Al bloc `output:`, hem reemplaçat `path("cowpy-${input_file}")` amb `path("${prefix}.txt")`.**
   Això simplement reitera quin serà el camí del fitxer segons el que està escrit a la línia de comandes.

Com a resultat, el nom del fitxer de sortida ara es construeix utilitzant un valor per defecte raonable (l'identificador del metamap) combinat amb l'extensió de format de fitxer apropiada.

#### 1.4.2. Configurar `ext.prefix` al fitxer `modules.config`

En aquest cas, el valor per defecte raonable no és prou expressiu per al nostre gust; volem utilitzar un patró de nomenament personalitzat que inclogui el nom de l'eina, `cowpy-<id>.txt`, com teníem abans.

Ho farem configurant `ext.prefix` a `modules.config`, igual que vam fer per al paràmetre `character` amb `ext.args`, excepte que aquesta vegada el bloc `withName: 'COWPY' {}` ja existeix, i només necessitem afegir la línia següent:

```groovy title="Code to add"
ext.prefix = { "cowpy-${meta.id}" }
```

Això compon la cadena que volem.
Tingueu en compte que un cop més utilitzem claus, aquesta vegada per indicar a Nextflow que avaluï el valor de `meta.id` en temps d'execució.

Afegim-ho.

Obriu `conf/modules.config` i afegiu el codi de configuració dins del bloc `process {}` com es mostra a continuació.

=== "Després"

    ```groovy title="core-hello/conf/modules.config" linenums="21" hl_lines="3"
        withName: 'COWPY' {
            ext.args = { "-c ${params.character}" }
            ext.prefix = { "cowpy-${meta.id}" }
        }
    ```

=== "Abans"

    ```groovy title="core-hello/conf/modules.config" linenums="21"
        withName: 'COWPY' {
            ext.args = { "-c ${params.character}" }
        }
    ```

En cas que us ho estigueu preguntant, el closure `ext.prefix` té accés a la peça correcta de metadades perquè la configuració s'avalua en el context de l'execució del procés, on les metadades estan disponibles.

#### 1.4.3. Executar el pipeline per provar-lo

Comprovem que el workflow encara funciona com s'espera.

```bash
nextflow run . --outdir core-hello-results -profile test,docker --validate_params false
```

??? success "Sortida de la comanda"

    ```console
     N E X T F L O W   ~  version 25.04.3

    Launching `./main.nf` [admiring_turing] DSL2 - revision: b9e9b3b8de

    Input/output options
      input                     : /workspaces/training/hello-nf-core/core-hello/assets/greetings.csv
      outdir                    : core-hello-results

    Institutional config options
      config_profile_name       : Test profile
      config_profile_description: Minimal test dataset to check pipeline function

    Generic options
      validate_params           : false
      trace_report_suffix       : 2025-12-27_06-29-02

    Core Nextflow options
      runName                   : admiring_turing
      containerEngine           : docker
      launchDir                 : /workspaces/training/hello-nf-core/core-hello
      workDir                   : /workspaces/training/hello-nf-core/core-hello/work
      projectDir                : /workspaces/training/hello-nf-core/core-hello
      userName                  : root
      profile                   : test,docker
      configFiles               : /workspaces/training/hello-nf-core/core-hello/nextflow.config

    !! Only displaying parameters that differ from the pipeline defaults !!
    ------------------------------------------------------
    executor >  local (8)
    [b2/e08524] CORE_HELLO:HELLO:sayHello (3)       [100%] 3 of 3 ✔
    [13/88939f] CORE_HELLO:HELLO:convertToUpper (3) [100%] 3 of 3 ✔
    [23/4554e1] CORE_HELLO:HELLO:CAT_CAT (test)     [100%] 1 of 1 ✔
    [a3/c6cbe9] CORE_HELLO:HELLO:COWPY              [100%] 1 of 1 ✔
    -[core/hello] Pipeline completed successfully-
    ```

Mireu la sortida al directori de resultats.
Hauríeu de veure el fitxer de sortida de cowpy amb el mateix nomenament que abans: `cowpy-test.txt`, basat en el nom de lot per defecte.

??? abstract "Contingut del directori"

    ```console hl_lines="3"
    results
    ├── Bonjour-output.txt
    ├── cowpy-test.txt
    ├── Hello-output.txt
    ├── Holà-output.txt
    ├── UPPER-Bonjour-output.txt
    ├── UPPER-Hello-output.txt
    └── UPPER-Holà-output.txt
    ```

Sentiu-vos lliures de canviar la configuració `ext.prefix` a `conf/modules.config` per satisfer-vos que podeu canviar el patró de nomenament sense haver de fer cap canvi al codi del mòdul o del workflow.

Alternativament, també podeu provar d'executar això de nou amb un paràmetre `--batch` diferent especificat a la línia de comandes per satisfer-vos que aquesta part encara és personalitzable sobre la marxa.

Això demostra com `ext.prefix` us permet mantenir la vostra convenció de nomenament preferida mentre manteniu la interfície del mòdul flexible.

Per resumir els beneficis d'aquest enfocament:

- **Nomenament estandarditzat**: Els fitxers de sortida normalment es nomenen utilitzant IDs de mostra de les metadades
- **Configurable**: Els usuaris poden sobreescriure el nomenament per defecte si cal
- **Consistent**: Tots els mòduls nf-core segueixen aquest patró
- **Predictible**: És fàcil saber com es diran els fitxers de sortida

Força bé, oi?
Bé, hi ha un canvi més important que hem de fer per millorar el nostre mòdul perquè s'ajusti a les directrius nf-core.

### 1.5. Centralitzar la configuració de publicació

Potser heu notat que hem estat publicant sortides a dos directoris diferents:

- **`results`** — El directori de sortida original que hem estat utilitzant des del principi per als nostres mòduls locals, establert individualment utilitzant directives `publishDir` per mòdul;
- **`core-hello-results`** — El directori de sortida establert amb `--outdir` a la línia de comandes, que ha estat rebent els registres nf-core i els resultats publicats per `CAT_CAT`.

Això és desordenat i subòptim; seria millor tenir una ubicació per a tot.
Per descomptat, podríem anar a cadascun dels nostres mòduls locals i actualitzar la directiva `publishDir` manualment per utilitzar el directori `core-hello-results`, però què passa la propera vegada que decidim canviar el directori de sortida?

Tenir mòduls individuals prenent decisions de publicació clarament no és el camí a seguir, especialment en un món on el mateix mòdul podria utilitzar-se en molts pipelines diferents, per persones que tenen necessitats o preferències diferents.
Volem poder controlar on es publiquen les sortides al nivell de la configuració del workflow.

"Ei," podríeu dir, "`CAT_CAT` està enviant les seves sortides a `--outdir`. Potser hauríem de copiar la seva directiva `publishDir`?"

Sí, és una gran idea.

Excepte que no té una directiva `publishDir`. (Endavant, mireu el codi del mòdul.)

Això és perquè els pipelines nf-core centralitzen el control al nivell del workflow configurant `publishDir` a `conf/modules.config` en lloc de fer-ho en mòduls individuals.
Específicament, la plantilla nf-core declara una directiva `publishDir` per defecte (amb una estructura de directoris predefinida) que s'aplica a tots els mòduls tret que es proporcioni una directiva que la sobreescrigui.

No sona genial? Podria ser que per aprofitar aquesta directiva per defecte, tot el que necessitem fer és eliminar la directiva `publishDir` actual dels nostres mòduls locals?

Provem-ho amb `COWPY` per veure què passa, després mirarem el codi de la configuració per defecte per entendre com funciona.

Finalment, demostrarem com sobreescriure el comportament per defecte si es desitja.

#### 1.5.1. Eliminar la directiva `publishDir` de `COWPY`

Fem això.
Obriu el fitxer del mòdul `cowpy.nf` (sota `core-hello/modules/local/`) i elimineu la directiva `publishDir` com es mostra a continuació.

=== "Després"

    ```groovy title="core-hello/modules/local/cowpy.nf (excerpt)" linenums="1"
    #!/usr/bin/env nextflow

    // Generate ASCII art with cowpy (https://github.com/jeffbuttars/cowpy)
    process COWPY {

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
        conda 'conda-forge::cowpy==1.1.5'
    ```

=== "Abans"

    ```groovy title="core-hello/modules/local/cowpy.nf (excerpt)" linenums="1" hl_lines="6"
    #!/usr/bin/env nextflow

    // Generate ASCII art with cowpy (https://github.com/jeffbuttars/cowpy)
    process COWPY {

        publishDir 'results', mode: 'copy'

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
        conda 'conda-forge::cowpy==1.1.5'

    ```

Això és tot!

#### 1.5.2. Executar el pipeline per provar-lo

Vegem què passa si executem el pipeline ara.

```bash
nextflow run . --outdir core-hello-results -profile test,docker --validate_params false
```

??? success "Sortida de la comanda"

    ```console
     N E X T F L O W   ~  version 25.04.3

    Launching `./main.nf` [silly_caravaggio] DSL2 - revision: b9e9b3b8de

    Input/output options
      input                     : /workspaces/training/hello-nf-core/core-hello/assets/greetings.csv
      outdir                    : core-hello-results

    Institutional config options
      config_profile_name       : Test profile
      config_profile_description: Minimal test dataset to check pipeline function

    Generic options
      validate_params           : false
      trace_report_suffix       : 2025-12-27_06-35-56

    Core Nextflow options
      runName                   : silly_caravaggio
      containerEngine           : docker
      launchDir                 : /workspaces/training/hello-nf-core/core-hello
      workDir                   : /workspaces/training/hello-nf-core/core-hello/work
      projectDir                : /workspaces/training/hello-nf-core/core-hello
      userName                  : root
      profile                   : test,docker
      configFiles               : /workspaces/training/hello-nf-core/core-hello/nextflow.config

    !! Only displaying parameters that differ from the pipeline defaults !!
    ------------------------------------------------------
    executor >  local (8)
    [db/39978e] CORE_HELLO:HELLO:sayHello (3)       [100%] 3 of 3 ✔
    [b5/bf6a8d] CORE_HELLO:HELLO:convertToUpper (3) [100%] 3 of 3 ✔
    [b7/c61842] CORE_HELLO:HELLO:CAT_CAT (test)     [100%] 1 of 1 ✔
    [46/5839d6] CORE_HELLO:HELLO:COWPY              [100%] 1 of 1 ✔
    -[core/hello] Pipeline completed successfully-
    ```

Mireu el vostre directori de treball actual.
Ara el `core-hello-results` també conté les sortides del mòdul `COWPY`.

??? abstract "Contingut del directori"

    ```console hl_lines="4-5"
    core-hello-results/
    ├── cat
    │   └── test.txt
    ├── cowpy
    │   └── cowpy-test.txt
    └── pipeline_info
        ├── execution_report_2025-12-27_06-16-55.html
        ├── execution_report_2025-12-27_06-23-13.html
        ├── execution_report_2025-12-27_06-29-02.html
        ├── execution_report_2025-12-27_06-35-56.html
        ├── execution_timeline_2025-12-27_06-16-55.html
        ├── execution_timeline_2025-12-27_06-23-13.html
        ├── execution_timeline_2025-12-27_06-29-02.html
        ├── execution_timeline_2025-12-27_06-35-56.html
        ├── execution_trace_2025-12-27_06-16-55.txt
        ├── execution_trace_2025-12-27_06-23-13.txt
        ├── execution_trace_2025-12-27_06-29-02.txt
        ├── execution_trace_2025-12-27_06-35-56.txt
        ├── hello_software_versions.yml
        ├── params_2025-12-27_06-17-00.json
        ├── params_2025-12-27_06-23-17.json
        ├── params_2025-12-27_06-29-07.json
        ├── params_2025-12-27_06-36-01.json
        ├── pipeline_dag_2025-12-27_06-16-55.html
        ├── pipeline_dag_2025-12-27_06-23-13.html
        ├── pipeline_dag_2025-12-27_06-29-02.html
        └── pipeline_dag_2025-12-27_06-35-56.html
    ```

Podeu veure que Nextflow ha creat aquesta jerarquia de directoris basant-se en els noms del workflow i del mòdul.

El codi responsable viu al fitxer `conf/modules.config`.
Aquesta és la configuració `publishDir` per defecte que forma part de la plantilla nf-core i s'aplica a tots els processos:

```groovy
process {
    publishDir = [
        path: { "${params.outdir}/${task.process.tokenize(':')[-1].tokenize('_')[0].toLowerCase()}" },
        mode: params.publish_dir_mode,
        saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
    ]
}
```

Això pot semblar complicat, així que vegem cadascun dels tres components:

- **`path:`** Determina el directori de sortida basat en el nom del procés.
  El nom complet d'un procés contingut a `task.process` inclou la jerarquia d'importacions de workflow i mòdul (com ara `CORE_HELLO:HELLO:CAT_CAT`).
  Les operacions `tokenize` eliminen aquesta jerarquia per obtenir només el nom del procés, després prenen la primera part abans de qualsevol guió baix (si escau), i la converteixen a minúscules.
  Això és el que determina que els resultats de `CAT_CAT` es publiquin a `${params.outdir}/cat/`.
- **`mode:`** Controla com es publiquen els fitxers (còpia, enllaç simbòlic, etc.).
  Això és configurable mitjançant el paràmetre `params.publish_dir_mode`.
- **`saveAs:`** Filtra quins fitxers publicar.
  Aquest exemple exclou fitxers `versions.yml` retornant `null` per a ells, evitant que es publiquin.

Això proporciona una lògica consistent per organitzar sortides.

La sortida es veu encara millor quan tots els mòduls d'un pipeline adopten aquesta convenció, així que sentiu-vos lliures d'anar a eliminar les directives `publishDir` dels altres mòduls del vostre pipeline.
Aquest valor per defecte s'aplicarà fins i tot als mòduls que no vam modificar explícitament per seguir les directrius nf-core.

Dit això, podeu decidir que voleu organitzar les vostres entrades de manera diferent, i la bona notícia és que és fàcil fer-ho.

#### 1.5.3. Sobreescriure el valor per defecte

Per sobreescriure la directiva `publishDir` per defecte, simplement podeu afegir les vostres pròpies directives al fitxer `conf/modules.config`.

Per exemple, podríeu sobreescriure el valor per defecte per a un sol procés utilitzant el selector `withName:`, com en aquest exemple on afegim una directiva `publishDir` personalitzada per al procés 'COWPY'.

```groovy title="core-hello/conf/modules.config" linenums="13" hl_lines="8-10"
process {
    publishDir = [
        path: { "${params.outdir}/${task.process.tokenize(':')[-1].tokenize('_')[0].toLowerCase()}" },
        mode: params.publish_dir_mode,
        saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
    ]

    withName: 'COWPY' {
        ext.args = { "-c ${params.character}" }
        publishDir = [
            path: 'my_custom_results'
        ]
    }
}
```

En realitat no farem aquest canvi, però sentiu-vos lliures de jugar amb això i veure quina lògica podeu implementar.

El punt és que aquest sistema us permet tenir el millor dels dos mons: consistència per defecte i la flexibilitat de personalitzar la configuració sota demanda.

Per resumir, obteniu:

- **Font única de veritat**: Tota la configuració de publicació viu a `modules.config`
- **Valor per defecte útil**: Els processos funcionen immediatament sense configuració per mòdul
- **Personalització fàcil**: Sobreescriviu el comportament de publicació a la configuració, no al codi del mòdul
- **Mòduls portables**: Els mòduls no codifiquen ubicacions de sortida

Això completa el conjunt de funcionalitats de mòduls nf-core que absolutament hauríeu d'aprendre a utilitzar, però n'hi ha d'altres sobre les quals podeu llegir a les [especificacions de mòduls nf-core](https://nf-co.re/docs/guidelines/components/modules).

### Conclusió

Ara sabeu com adaptar mòduls locals per seguir les convencions nf-core:

- Dissenyeu els vostres mòduls per acceptar i propagar tuples de metadades;
- Utilitzeu `ext.args` per mantenir les interfícies de mòdul mínimes i portables;
- Utilitzeu `ext.prefix` per a un nomenament de fitxers de sortida configurable i estandarditzat;
- Adopteu la directiva `publishDir` centralitzada per defecte per a una estructura de directori de resultats consistent.

### Què segueix?

Apreneu com utilitzar les eines integrades basades en plantilles d'nf-core per crear mòduls de manera fàcil.

---

## 2. Crear un mòdul amb les eines nf-core

Ara que heu après els patrons de mòduls nf-core aplicant-los manualment, vegem com crearíeu mòduls a la pràctica.

### 2.1. Generar una estructura de mòdul a partir d'una plantilla

Similar al que existeix per crear pipelines, el projecte nf-core proporciona eines per generar mòduls estructurats correctament basats en una plantilla, amb tots aquests patrons integrats des del principi.

#### 2.1.1. Executar la comanda de creació de mòdul

La comanda `nf-core modules create` genera una plantilla de mòdul que ja segueix totes les convencions que heu après.

Creem una nova versió del mòdul `COWPY` amb una plantilla mínima executant aquesta comanda:

```bash
nf-core modules create --empty-template COWPY
```

La bandera `--empty-template` crea una plantilla inicial neta sense codi extra, facilitant veure l'estructura essencial.

La comanda s'executa de manera interactiva, guiant-vos a través de la configuració.
Busca automàticament informació de l'eina des de repositoris de paquets com Bioconda i bio.tools per pre-omplir metadades.

Se us demanarà diverses opcions de configuració:

- **Informació de l'autor**: El vostre nom d'usuari de GitHub per a l'atribució
- **Etiqueta de recursos**: Un conjunt predefinit de requisits computacionals.
  El projecte nf-core proporciona etiquetes estàndard com `process_single` per a eines lleugeres i `process_high` per a les exigents.
  Aquestes etiquetes ajuden a gestionar l'assignació de recursos entre diferents entorns d'execució.
- **Requisit de metadades**: Si el mòdul necessita informació específica de mostra mitjançant un mapa `meta` (normalment sí per a mòduls de processament de dades).

L'eina gestiona la complexitat de trobar informació de paquets i configurar l'estructura, permetent-vos centrar-vos en implementar la lògica específica de l'eina.

#### 2.1.2. Examinar l'estructura del mòdul

L'eina crea una estructura de mòdul completa a `modules/local/` (o `modules/nf-core/` si esteu al repositori nf-core/modules):

??? abstract "Contingut del directori"

    ```console
    modules/local/cowpy
    ├── environment.yml
    ├── main.nf
    ├── meta.yml
    └── tests
        └── main.nf.test
    ```

Cada fitxer té un propòsit específic:

- **`main.nf`**: Definició del procés amb tots els patrons nf-core integrats
- **`meta.yml`**: Documentació del mòdul descrivint entrades, sortides i l'eina
- **`environment.yml`**: Especificació de l'entorn Conda per a dependències
- **`tests/main.nf.test`**: Casos de prova nf-test per validar que el mòdul funciona

!!! tip "Apreneu més sobre proves"

    El fitxer de prova generat utilitza nf-test, un marc de proves per a pipelines i mòduls Nextflow. Per aprendre com escriure i executar aquestes proves, consulteu la [missió secundària nf-test](../side_quests/nf-test.md).

El `main.nf` generat inclou tots els patrons que acabeu d'aprendre, més algunes funcionalitats addicionals:

```groovy title="modules/local/cowpy/main.nf" hl_lines="11 21 22"
process COWPY {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE':
        'biocontainers/YOUR-TOOL-HERE' }"

    input:
    tuple val(meta), path(input)        // Patró 1: Tuples de metadades ✓

    output:
    tuple val(meta), path("*"), emit: output
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''                // Patró 2: ext.args ✓
    def prefix = task.ext.prefix ?: "${meta.id}"  // Patró 3: ext.prefix ✓

    """
    // Add your tool command here

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        COWPY: \$(cowpy --version)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    echo $args
    touch ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        COWPY: \$(cowpy --version)
    END_VERSIONS
    """
}
```

Noteu com tots els patrons que vau aplicar manualment anteriorment ja hi són!

La plantilla també inclou diverses convencions addicionals d'nf-core.
Algunes d'aquestes funcionen immediatament, mentre que d'altres són marcadors de posició que haurem d'omplir, com es descriu a continuació.

**Funcionalitats que funcionen tal com estan:**

- **`tag "$meta.id"`**: Afegeix l'ID de mostra als noms de procés als registres per a un seguiment més fàcil
- **`label 'process_single'`**: Etiqueta de recursos per configurar requisits de CPU/memòria
- **Bloc `when:`**: Permet l'execució condicional mitjançant configuració `task.ext.when`

Aquestes funcionalitats ja són funcionals i fan que els mòduls siguin més mantenibles.

**Marcadors de posició que personalitzarem a continuació:**

- **Blocs `input:` i `output:`**: Declaracions genèriques que actualitzarem per coincidir amb la nostra eina
- **Bloc `script:`**: Conté un comentari on afegirem la comanda `cowpy`
- **Bloc `stub:`**: Plantilla que actualitzarem per produir les sortides correctes
- **Contenidor i entorn**: Marcadors de posició que omplirem amb informació de paquets

Les següents seccions expliquen com completar aquestes personalitzacions.

### 2.2. Configurar el contenidor i l'entorn conda

Les directrius nf-core requereixen que especifiquem tant un contenidor com un entorn Conda com a part del mòdul.

#### 2.2.1. Contenidor

Per al contenidor, podeu utilitzar [Seqera Containers](https://seqera.io/containers/) per construir automàticament un contenidor des de qualsevol paquet Conda, incloent paquets conda-forge.
En aquest cas estem utilitzant el mateix contenidor preconstruït que abans.

El codi per defecte ofereix alternar entre Docker i Singularity, però simplificarem aquesta línia i només especificarem el contenidor Docker que vam obtenir de Seqera Containers anteriorment.

=== "Després"

```groovy title="modules/local/cowpy/main.nf" linenums="3" hl_lines="6"
process COWPY {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273"
```

=== "Abans"

```groovy title="modules/local/cowpy/main.nf" linenums="3" hl_lines="6"
process COWPY {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE':
        'biocontainers/YOUR-TOOL-HERE' }"
```

#### 2.2.2. Entorn Conda

Per a l'entorn Conda, el codi del mòdul especifica `conda "${moduleDir}/environment.yml"` que significa que s'ha de configurar al fitxer `environment.yml`.

L'eina de creació de mòduls ens va advertir que no podia trobar el paquet `cowpy` a Bioconda (el canal principal per a eines bioinformàtiques).
No obstant això, `cowpy` està disponible a conda-forge, així que podeu completar l'`environment.yml` així:

=== "Després"

    ```yaml title="modules/local/cowpy/environment.yml"  linenums="1" hl_lines="1 3 5"
    name: COWPY
    channels:
      - conda-forge
    dependencies:
      - cowpy=1.1.5
    ```

=== "Abans"

    ```yaml title="modules/local/cowpy/environment.yml" linenums="1"
    ---
    # yaml-language-server: $schema=https://raw.githubusercontent.com/nf-core/modules/master/modules/environment-schema.json
    channels:
      - conda-forge
      - bioconda
    dependencies:
      # TODO nf-core: List required Conda package(s).
      #               Software MUST be pinned to channel (i.e. "bioconda"), version (i.e. "1.10").
      #               For Conda, the build (i.e. "h9402c20_2") must be EXCLUDED to support installation on different operating systems.
      - "YOUR-TOOL-HERE"
    ```

Per a la submissió a nf-core, hauríem de seguir els valors per defecte més de prop, però per al nostre propi ús podem simplificar el codi d'aquesta manera.

!!! tip "Paquets Bioconda vs conda-forge"

    - **Paquets Bioconda**: Obtenen automàticament BioContainers construïts, proporcionant contenidors llestos per utilitzar
    - **Paquets conda-forge**: Poden utilitzar Seqera Containers per construir contenidors sota demanda des de la recepta Conda

    La majoria d'eines bioinformàtiques estan a Bioconda, però per a eines conda-forge, Seqera Containers proporciona una solució fàcil per a la contenidorització.

### 2.3. Connectar la lògica de `COWPY`

Ara actualitzem els elements de codi que són específics del que fa el procés `COWPY`: les entrades i sortides, i el bloc script.

#### 2.3.1. Entrades i sortides

La plantilla generada inclou declaracions d'entrada i sortida genèriques que haureu de personalitzar per a la vostra eina específica.
Mirant enrere al nostre mòdul `COWPY` manual de la secció 1, podem utilitzar-lo com a guia.

Actualitzeu els blocs d'entrada i sortida:

=== "Després"

    ```groovy title="modules/local/cowpy/main.nf" linenums="8" hl_lines="2 5"
    input:
    tuple val(meta), path(input_file)

    output:
    tuple val(meta), path("${prefix}.txt"), emit: cowpy_output
    path "versions.yml"           , emit: versions
    ```

=== "Abans"

    ```groovy title="modules/local/cowpy/main.nf" linenums="8" hl_lines="2 5"
    input:
    tuple val(meta), path(input)

    output:
    tuple val(meta), path("*"), emit: output
    path "versions.yml"           , emit: versions
    ```

Això especifica:

- El nom del paràmetre del fitxer d'entrada (`input_file` en lloc de `input` genèric)
- El nom del fitxer de sortida utilitzant el patró de prefix configurable (`${prefix}.txt` en lloc del comodí `*`)
- Un nom emit descriptiu (`cowpy_output` en lloc de `output` genèric)

Si esteu utilitzant el servidor de llenguatge Nextflow per validar la sintaxi, la part `${prefix}` es marcarà com un error en aquesta etapa perquè encara no l'hem afegit al bloc script.
Fem-ho ara.

#### 2.3.2. El bloc script

La plantilla proporciona un comentari marcador de posició al bloc script on hauríeu d'afegir la comanda real de l'eina.

Basant-nos en el mòdul que vam escriure manualment abans, hauríem de fer les següents edicions:

=== "Després"

    ```groovy title="modules/local/cowpy/main.nf" linenums="15" hl_lines="3 6"
    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    cat $input_file | cowpy $args > ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        COWPY: \$(cowpy --version)
    END_VERSIONS
    """
    ```

=== "Abans"

    ```groovy title="modules/local/cowpy/main.nf" linenums="15" hl_lines="6"
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    // Add your tool command here

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        COWPY: \$(cowpy --version)
    END_VERSIONS
    """
    ```

Canvis clau:

- Canviar `def prefix` a només `prefix` (sense `def`) per fer-lo accessible al bloc output
- Reemplaçar el comentari amb la comanda `cowpy` real que utilitza tant `$args` com `${prefix}.txt`

Tingueu en compte que si no haguéssim fet ja el treball d'afegir la configuració `ext.args` i `ext.prefix` per al procés `COWPY` al fitxer `modules.config`, hauríem de fer-ho ara.

#### 2.3.3. Implementar el bloc stub

En el context de Nextflow, un bloc [stub](https://www.nextflow.io/docs/latest/process.html#stub) us permet definir un script lleuger i fictici utilitzat per a prototipatge ràpid i proves de la lògica d'un pipeline sense executar la comanda real.

No us preocupeu massa si això sembla misteriós; ho incloem per completesa però també podeu simplement eliminar la secció stub si no voleu tractar amb ella, ja que és completament opcional.

=== "Després"

    ```groovy title="modules/local/cowpy/main.nf" linenums="27" hl_lines="3 6"
    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        COWPY: \$(cowpy --version)
    END_VERSIONS
    """
    ```

=== "Abans"

    ```groovy title="modules/local/cowpy/main.nf" linenums="27" hl_lines="3 6"
    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    echo $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        COWPY: \$(cowpy --version)
    END_VERSIONS
    """
    ```

Canvis clau:

- Canviar `def prefix` a només `prefix` per coincidir amb el bloc script
- Eliminar la línia `echo $args` (que era només codi marcador de posició de plantilla)
- El stub crea un fitxer `${prefix}.txt` buit que coincideix amb el que produeix el bloc script

Això us permet provar la lògica del workflow i la gestió de fitxers sense esperar que l'eina real s'executi.

Un cop hàgiu completat la configuració de l'entorn (secció 2.2), entrades/sortides (secció 2.3.1), bloc script (secció 2.3.2) i bloc stub (secció 2.3.3), el mòdul està llest per provar!

### 2.4. Intercanviar el nou mòdul `COWPY` i executar el pipeline

Tot el que necessitem fer per provar aquesta nova versió del mòdul `COWPY` és canviar la declaració d'importació al fitxer del workflow `hello.nf` per apuntar al nou fitxer.

=== "Després"

    ```groovy title="workflows/hello.nf" linenums="1" hl_lines="10"
    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    include { paramsSummaryMap       } from 'plugin/nf-schema'
    include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
    include { sayHello               } from '../modules/local/sayHello.nf'
    include { convertToUpper         } from '../modules/local/convertToUpper.nf'
    include { COWPY                  } from '../modules/local/cowpy/main.nf'
    include { CAT_CAT                } from '../modules/nf-core/cat/cat/main'
    ```

=== "Abans"

    ```groovy title="modules/local/cowpy/main.nf" linenums="1" hl_lines="10"
    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    include { paramsSummaryMap       } from 'plugin/nf-schema'
    include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
    include { sayHello               } from '../modules/local/sayHello.nf'
    include { convertToUpper         } from '../modules/local/convertToUpper.nf'
    include { COWPY                  } from '../modules/local/cowpy.nf'
    include { CAT_CAT                } from '../modules/nf-core/cat/cat/main'
    ```

Executem el pipeline per provar-lo.

```bash
nextflow run . --outdir core-hello-results -profile test,docker --validate_params false
```

??? success "Sortida de la comanda"

    ```console hl_lines="33"
      N E X T F L O W   ~  version 25.04.3

    Launching `./main.nf` [prickly_neumann] DSL2 - revision: b9e9b3b8de

    Input/output options
      input                     : /workspaces/training/hello-nf-core/core-hello/assets/greetings.csv
      outdir                    : core-hello-results

    Institutional config options
      config_profile_name       : Test profile
      config_profile_description: Minimal test dataset to check pipeline function

    Generic options
      validate_params           : false
      trace_report_suffix       : 2025-12-27_08-23-51

    Core Nextflow options
      runName                   : prickly_neumann
      containerEngine           : docker
      launchDir                 : /workspaces/training/hello-nf-core/core-hello
      workDir                   : /workspaces/training/hello-nf-core/core-hello/work
      projectDir                : /workspaces/training/hello-nf-core/core-hello
      userName                  : root
      profile                   : test,docker
      configFiles               : /workspaces/training/hello-nf-core/core-hello/nextflow.config

    !! Only displaying parameters that differ from the pipeline defaults !!
    ------------------------------------------------------
    executor >  local (8)
    [e9/008ede] CORE_HELLO:HELLO:sayHello (3)       [100%] 3 of 3 ✔
    [f0/d70cfe] CORE_HELLO:HELLO:convertToUpper (3) [100%] 3 of 3 ✔
    [be/0ecc58] CORE_HELLO:HELLO:CAT_CAT (test)     [100%] 1 of 1 ✔
    [11/8e082f] CORE_HELLO:HELLO:COWPY (test)       [100%] 1 of 1 ✔
    -[core/hello] Pipeline completed successfully-
    ```

Això produeix els mateixos resultats que abans.

### Conclusió

Ara sabeu com utilitzar les eines integrades d'nf-core per crear mòduls de manera eficient utilitzant plantilles en lloc d'escriure tot des de zero.

### Què segueix?

Apreneu quins són els beneficis de contribuir mòduls a nf-core i quins són els principals passos i requisits implicats.

---

## 3. Contribuir mòduls de tornada a nf-core

El repositori [nf-core/modules](https://github.com/nf-core/modules) dóna la benvinguda a contribucions de mòduls ben provats i estandarditzats.

### 3.1. Per què contribuir?

Contribuir els vostres mòduls a nf-core:

- Fa que les vostres eines estiguin disponibles per a tota la comunitat nf-core a través del catàleg de mòduls a [nf-co.re/modules](https://nf-co.re/modules)
- Assegura el manteniment i millores contínues de la comunitat
- Proporciona garantia de qualitat mitjançant revisió de codi i proves automatitzades
- Dóna visibilitat i reconeixement al vostre treball

### 3.2. Llista de verificació del contribuïdor

Per contribuir un mòdul a nf-core, haureu de passar pels següents passos:

1. Comproveu si ja existeix a [nf-co.re/modules](https://nf-co.re/modules)
2. Feu un fork del repositori [nf-core/modules](https://github.com/nf-core/modules)
3. Utilitzeu `nf-core modules create` per generar la plantilla
4. Ompliu la lògica del mòdul i les proves
5. Proveu amb `nf-core modules test tool/subtool`
6. Valideu amb `nf-core modules lint tool/subtool`
7. Envieu una pull request

Per a instruccions detallades, consulteu el [tutorial de components nf-core](https://nf-co.re/docs/tutorials/nf-core_components/components).

### 3.3. Recursos

- **Tutorial de components**: [Guia completa per crear i contribuir mòduls](https://nf-co.re/docs/tutorials/nf-core_components/components)
- **Especificacions de mòduls**: [Requisits tècnics i directrius](https://nf-co.re/docs/guidelines/components/modules)
- **Suport de la comunitat**: [Slack d'nf-core](https://nf-co.re/join) - Uniu-vos al canal `#modules`

### Conclusió

Ara sabeu com crear mòduls nf-core! Heu après els quatre patrons clau que fan que els mòduls siguin portables i mantenibles:

- Les **tuples de metadades** propaguen metadades a través del workflow
- **`ext.args`** simplifica les interfícies de mòdul gestionant arguments opcionals mitjançant configuració
- **`ext.prefix`** estandarditza el nomenament de fitxers de sortida
- La **publicació centralitzada** mitjançant `publishDir` configurada a `modules.config` en lloc de codificada als mòduls

Transformant `COWPY` pas a pas, heu desenvolupat una comprensió profunda d'aquests patrons, fent-vos capaços de treballar amb, depurar i crear mòduls nf-core.
A la pràctica, utilitzareu `nf-core modules create` per generar mòduls estructurats correctament amb aquests patrons integrats des del principi.

Finalment, heu après com contribuir mòduls a la comunitat nf-core, fent que les eines estiguin disponibles per a investigadors d'arreu del món mentre us beneficieu del manteniment continu de la comunitat.

### Què segueix?

Quan estigueu preparats, continueu a la [Part 5: Validació d'entrada](./05_input_validation.md) per aprendre com afegir validació d'entrada basada en esquema al vostre pipeline.
