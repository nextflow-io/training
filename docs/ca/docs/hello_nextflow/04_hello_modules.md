# Part 4: Hello Modules

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducció assistida per IA - [més informació i suggeriments](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/43Ot-f0iOME?si=y8lAedhEHWaTV4zd&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

/// caption
:fontawesome-brands-youtube:{ .youtube } Consulteu [la llista de reproducció completa](https://youtube.com/playlist?list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&si=eF7cLR62goy-lc6n) al canal de YouTube de Nextflow.

:green_book: La transcripció del vídeo està disponible [aquí](./transcripts/04_hello_modules.md).
///

Aquesta secció cobreix com organitzar el codi del vostre workflow per fer el desenvolupament i el manteniment del vostre pipeline més eficient i sostenible.
Específicament, demostrarem com utilitzar [**mòduls**](https://nextflow.io/docs/latest/module.html).

A Nextflow, un **mòdul** és un fitxer de codi independent, sovint encapsulant una única definició de procés.
Per utilitzar un mòdul en un workflow, només cal afegir una única línia amb una declaració `include` al fitxer de codi del vostre workflow; llavors podeu integrar el procés al workflow de la mateixa manera que ho faríeu normalment.
Això fa possible reutilitzar definicions de processos en múltiples workflows sense produir múltiples còpies del codi.

Quan vam començar a desenvolupar el nostre workflow, vam escriure tot en un únic fitxer de codi.
Ara mourem els processos a mòduls individuals.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nextflow/img/modules.svg"
</figure>

Això farà que el nostre codi sigui més compartible, flexible i mantenible.

??? info "Com començar des d'aquesta secció"

    Aquesta secció del curs assumeix que heu completat les Parts 1-3 del curs [Hello Nextflow](./index.md), però si us sentiu còmodes amb els conceptes bàsics coberts en aquestes seccions, podeu començar des d'aquí sense fer res especial.

---

## 0. Escalfament: Executeu `hello-modules.nf`

Utilitzarem l'script de workflow `hello-modules.nf` com a punt de partida.
És equivalent a l'script produït en completar la Part 3 d'aquest curs de formació, excepte que hem canviat les destinacions de sortida:

```groovy title="hello-modules.nf" linenums="37" hl_lines="3 7 11 15"
output {
    first_output {
        path 'hello_modules'
        mode 'copy'
    }
    uppercased {
        path 'hello_modules'
        mode 'copy'
    }
    collected {
        path 'hello_modules'
        mode 'copy'
    }
    batch_report {
        path 'hello_modules'
        mode 'copy'
    }
}
```

Només per assegurar-nos que tot funciona, executeu l'script una vegada abans de fer cap canvi:

```bash
nextflow run hello-modules.nf
```

??? success "Sortida de la comanda"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `hello-modules.nf` [hopeful_avogadro] DSL2 - revision: b09af1237d

    executor >  local (7)
    [0f/8795c9] sayHello (3)       [100%] 3 of 3 ✔
    [6a/eb2510] convertToUpper (3) [100%] 3 of 3 ✔
    [af/479117] collectGreetings   [100%] 1 of 1 ✔
    ```

Com abans, trobareu els fitxers de sortida al directori especificat al bloc `output` (aquí, `results/hello_modules/`).

??? abstract "Contingut del directori"

    ```console
    results/hello_modules/
    ├── Bonjour-output.txt
    ├── COLLECTED-batch-output.txt
    ├── Hello-output.txt
    ├── Holà-output.txt
    ├── batch-report.txt
    ├── UPPER-Bonjour-output.txt
    ├── UPPER-Hello-output.txt
    └── UPPER-Holà-output.txt
    ```

Si això ha funcionat, esteu preparats per aprendre com modularitzar el codi del vostre workflow.

---

## 1. Creeu un directori per emmagatzemar mòduls

És una bona pràctica emmagatzemar els vostres mòduls en un directori específic.
Podeu anomenar aquest directori com vulgueu, però la convenció és anomenar-lo `modules/`.

```bash
mkdir modules
```

---

## 2. Creeu un mòdul per a `sayHello()`

En la seva forma més simple, convertir un procés existent en un mòdul és poc més que una operació de copiar i enganxar.
Crearem un fitxer base per al mòdul, copiarem el codi rellevant i després l'eliminarem del fitxer principal del workflow.

Després només caldrà afegir una declaració `include` perquè Nextflow sàpiga que ha d'incorporar el codi rellevant en temps d'execució.

### 2.1. Creeu un fitxer base per al nou mòdul

Creem un fitxer buit per al mòdul anomenat `sayHello.nf`.

```bash
touch modules/sayHello.nf
```

Això ens dóna un lloc on posar el codi del procés.

### 2.2. Moveu el codi del procés `sayHello` al fitxer del mòdul

Copieu tota la definició del procés del fitxer del workflow al fitxer del mòdul.

```groovy title="modules/sayHello.nf" linenums="1"
/*
 * Utilitza echo per imprimir 'Hello World!' a un fitxer
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

Un cop fet això, elimineu la definició del procés del fitxer del workflow.

### 2.3. Afegiu una declaració include abans del bloc workflow

La sintaxi per incloure un procés des d'un mòdul és força directa:

```groovy title="Syntax: include declaration"
include { <PROCESS_NAME> } from '<path_to_module>'
```

Inserim això sobre el bloc `params` i l'omplim adequadament.

=== "Després"

    ```groovy title="hello-modules.nf" linenums="44" hl_lines="1 2"
    // Inclou mòduls
    include { sayHello } from './modules/sayHello.nf'

    /*
    * Paràmetres del pipeline
    */
    params {
        greeting: Path = 'data/greetings.csv'
        batch: String = 'batch'
    }
    ```

=== "Abans"

    ```groovy title="hello-modules.nf" linenums="44"
    /*
    * Paràmetres del pipeline
    */
    params {
        greeting: Path = 'data/greetings.csv'
        batch: String = 'batch'
    }
    ```

Veieu que hem omplert el nom del procés, `sayHello`, i la ruta al fitxer que conté el codi del mòdul, `./modules/sayHello.nf`.

### 2.4. Executeu el workflow

Estem executant el workflow amb essencialment el mateix codi i entrades que abans, així que executem amb la bandera `-resume` i vegem què passa.

```bash
nextflow run hello-modules.nf -resume
```

??? success "Sortida de la comanda"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-modules.nf` [romantic_poisson] DSL2 - revision: 96edfa9ad3

    [f6/cc0107] sayHello (1)       | 3 of 3, cached: 3 ✔
    [3c/4058ba] convertToUpper (2) | 3 of 3, cached: 3 ✔
    [1a/bc5901] collectGreetings   | 1 of 1, cached: 1 ✔
    ```

Això hauria d'executar-se molt ràpidament perquè tot està en memòria cau.
Podeu comprovar les sortides publicades si voleu.

Nextflow ha reconegut que encara és tota la mateixa feina a fer, encara que el codi estigui dividit en múltiples fitxers.

### Conclusió

Sabeu com extreure un procés a un mòdul local i sabeu que fer això no trenca la capacitat de reprendre el workflow.

### Què segueix?

Practiqueu fent més mòduls.
Un cop n'heu fet un, podeu fer-ne un milió més...
Però fem-ne només dos més per ara.

---

## 3. Modularitzeu el procés `convertToUpper()`

### 3.1. Creeu un fitxer base per al nou mòdul

Creeu un fitxer buit per al mòdul anomenat `convertToUpper.nf`.

```bash
touch modules/convertToUpper.nf
```

### 3.2. Moveu el codi del procés `convertToUpper` al fitxer del mòdul

Copieu tota la definició del procés del fitxer del workflow al fitxer del mòdul.

```groovy title="modules/convertToUpper.nf" linenums="1"
/*
 * Utilitza una eina de substitució de text per convertir la salutació a majúscules
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
```

Un cop fet això, elimineu la definició del procés del fitxer del workflow.

### 3.3. Afegiu una declaració include abans del bloc `params`

Inseriu la declaració include sobre el bloc `params` i ompliu-la adequadament.

=== "Després"

    ```groovy title="hello-modules.nf" linenums="23" hl_lines="3"
    // Inclou mòduls
    include { sayHello } from './modules/sayHello.nf'
    include { convertToUpper } from './modules/convertToUpper.nf'

    /*
    * Paràmetres del pipeline
    */
    params {
        greeting: Path = 'data/greetings.csv'
        batch: String = 'batch'
    }
    ```

=== "Abans"

    ```groovy title="hello-modules.nf" linenums="23"
    // Inclou mòduls
    include { sayHello } from './modules/sayHello.nf'

    /*
    * Paràmetres del pipeline
    */
    params {
        greeting: Path = 'data/greetings.csv'
        batch: String = 'batch'
    }
    ```

Això hauria de començar a semblar molt familiar.

### 3.4. Executeu el workflow de nou

Executeu això amb la bandera `-resume`.

```bash
nextflow run hello-modules.nf -resume
```

??? success "Sortida de la comanda"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-modules.nf` [nauseous_heisenberg] DSL2 - revision: a04a9f2da0

    [c9/763d42] sayHello (3)       | 3 of 3, cached: 3 ✔
    [60/bc6831] convertToUpper (3) | 3 of 3, cached: 3 ✔
    [1a/bc5901] collectGreetings   | 1 of 1, cached: 1 ✔
    ```

Això encara hauria de produir la mateixa sortida que abans.

Dos fets, un més per fer!

---

## 4. Modularitzeu el procés `collectGreetings()`

### 4.1. Creeu un fitxer base per al nou mòdul

Creeu un fitxer buit per al mòdul anomenat `collectGreetings.nf`.

```bash
touch modules/collectGreetings.nf
```

### 4.2. Moveu el codi del procés `collectGreetings` al fitxer del mòdul

Copieu tota la definició del procés del fitxer del workflow al fitxer del mòdul.

```groovy title="modules/collectGreetings.nf" linenums="1"
/*
 * Recull salutacions en majúscules en un únic fitxer de sortida
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
```

Un cop fet això, elimineu la definició del procés del fitxer del workflow.

### 4.3. Afegiu una declaració include abans del bloc `params`

Inseriu la declaració include sobre el bloc `params` i ompliu-la adequadament.

=== "Després"

    ```groovy title="hello-modules.nf" linenums="3" hl_lines="4"
    // Inclou mòduls
    include { sayHello } from './modules/sayHello.nf'
    include { convertToUpper } from './modules/convertToUpper.nf'
    include { collectGreetings } from './modules/collectGreetings.nf'

    /*
    * Paràmetres del pipeline
    */
    params {
        greeting: Path = 'data/greetings.csv'
        batch: String = 'batch'
    }
    ```

=== "Abans"

    ```groovy title="hello-modules.nf" linenums="3"
    // Inclou mòduls
    include { sayHello } from './modules/sayHello.nf'
    include { convertToUpper } from './modules/convertToUpper.nf'

    /*
    * Paràmetres del pipeline
    */
    params {
        greeting: Path = 'data/greetings.csv'
        batch: String = 'batch'
    }
    ```

L'últim!

### 4.4. Executeu el workflow

Executeu això amb la bandera `-resume`.

```bash
nextflow run hello-modules.nf -resume
```

??? success "Sortida de la comanda"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `hello-modules.nf` [friendly_coulomb] DSL2 - revision: 7aa2b9bc0f

    [f6/cc0107] sayHello (1)       | 3 of 3, cached: 3 ✔
    [3c/4058ba] convertToUpper (2) | 3 of 3, cached: 3 ✔
    [1a/bc5901] collectGreetings   | 1 of 1, cached: 1 ✔
    ```

Això encara hauria de produir la mateixa sortida que abans.

### Conclusió

Sabeu com modularitzar múltiples processos en un workflow.

Felicitats, heu fet tota aquesta feina i absolutament res ha canviat en com funciona el pipeline!

Bromes a part, ara el vostre codi és més modular, i si decidiu escriure un altre pipeline que cridi un d'aquests processos, només cal que escriviu una curta declaració `include` per utilitzar el mòdul rellevant.
Això és millor que copiar i enganxar el codi, perquè si més endavant decidiu millorar el mòdul, tots els vostres pipelines heretaran les millores.

### Què segueix?

Feu una petita pausa si us ve de gust.

Quan estigueu preparats, passeu a [**Part 5: Hello Containers**](./05_hello_containers.md) per aprendre com utilitzar contenidors per gestionar dependències de programari de manera més convenient i reproduïble.

---

## Qüestionari

<quiz>
Què és un mòdul a Nextflow?
- [ ] Un fitxer de configuració
- [x] Un fitxer independent que pot contenir definicions de processos
- [ ] Una definició de workflow
- [ ] Un operador de canal

Més informació: [2. Creeu un mòdul per a `sayHello()`](#2-create-a-module-for-sayhello)
</quiz>

<quiz>
Quina convenció s'utilitza típicament per emmagatzemar fitxers de mòduls?
- [ ] Al mateix directori que el workflow
- [ ] En un directori `bin/`
- [x] En un directori `modules/`
- [ ] En un directori `lib/`

Més informació: [1. Creeu un directori per emmagatzemar mòduls](#1-create-a-directory-to-store-modules)
</quiz>

<quiz>
Quina és la sintaxi correcta per utilitzar un mòdul?

- [ ] `#!groovy import { SAYHELLO } from './modules/sayhello.nf'`
- [ ] `#!groovy require { SAYHELLO } from './modules/sayhello.nf'`
- [x] `#!groovy include { SAYHELLO } from './modules/sayhello.nf'`
- [ ] `#!groovy load { SAYHELLO } from './modules/sayhello.nf'`

Més informació: [2.3. Afegiu una declaració include](#23-add-an-include-declaration-before-the-workflow-block)
</quiz>

<quiz>
Què passa amb la funcionalitat `-resume` quan s'utilitzen mòduls?
- [ ] Ja no funciona
- [ ] Requereix configuració addicional
- [x] Funciona igual que abans
- [ ] Només funciona per a mòduls locals
</quiz>

<quiz>
Quins són els beneficis d'utilitzar mòduls? (Seleccioneu tots els que corresponguin)
- [x] Reutilització de codi entre workflows
- [x] Manteniment més fàcil
- [x] Millor organització del codi del workflow
- [ ] Velocitat d'execució més ràpida
</quiz>
