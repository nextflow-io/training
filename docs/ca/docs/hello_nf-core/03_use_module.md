# Part 3: Utilitzar un mòdul nf-core

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducció assistida per IA - [més informació i suggeriments](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

En aquesta tercera part del curs de formació Hello nf-core, us mostrem com trobar, instal·lar i utilitzar un mòdul nf-core existent al vostre pipeline.

Un dels grans avantatges de treballar amb nf-core és la capacitat d'aprofitar mòduls preconstruïts i provats del repositori [nf-core/modules](https://github.com/nf-core/modules).
En lloc d'escriure cada procés des de zero, podeu instal·lar i utilitzar mòduls mantinguts per la comunitat que segueixen les millors pràctiques.

Per demostrar com funciona això, substituirem el mòdul personalitzat `collectGreetings` pel mòdul `cat/cat` de nf-core/modules al pipeline `core-hello`.

??? info "Com començar des d'aquesta secció"

    Aquesta secció del curs assumeix que heu completat la [Part 2: Reescriure Hello per a nf-core](./02_rewrite_hello.md) i teniu un pipeline `core-hello` funcional.

    Si no heu completat la Part 2 o voleu començar de nou per a aquesta part, podeu utilitzar la solució `core-hello-part2` com a punt de partida.
    Executeu aquesta comanda des del directori `hello-nf-core/`:

    ```bash
    cp -r solutions/core-hello-part2 core-hello
    cd core-hello
    ```

    Això us proporciona un pipeline nf-core completament funcional preparat per afegir mòduls.
    Podeu comprovar que s'executa correctament executant la comanda següent:

    ```bash
    nextflow run . --outdir core-hello-results -profile test,docker --validate_params false
    ```

---

## 1. Trobar i instal·lar un mòdul nf-core adequat

Primer, aprendrem com trobar un mòdul nf-core existent i instal·lar-lo al nostre pipeline.

Volem substituir el procés `collectGreetings`, que utilitza la comanda Unix `cat` per concatenar múltiples fitxers de salutacions en un de sol.
Concatenar fitxers és una operació molt comuna, així que és raonable pensar que ja podria existir un mòdul a nf-core dissenyat per a aquest propòsit.

Comencem.

### 1.1. Explorar els mòduls disponibles al lloc web d'nf-core

El projecte nf-core manté un catàleg centralitzat de mòduls a [https://nf-co.re/modules](https://nf-co.re/modules).

Navegueu a la pàgina de mòduls al vostre navegador web i utilitzeu la barra de cerca per cercar 'concatenate'.

![resultats de cerca de mòduls](./img/module-search-results.png)

Com podeu veure, hi ha força resultats, molts d'ells mòduls dissenyats per concatenar tipus de fitxers molt específics.
Entre ells, hauríeu de veure un anomenat `cat_cat` que és de propòsit general.

!!! note "Convenció de nomenclatura de mòduls"

    El guió baix (`_`) s'utilitza com a substitut del caràcter barra (`/`) als noms de mòduls.

    Els mòduls nf-core segueixen la convenció de nomenclatura `programari/comanda` quan una eina proporciona múltiples comandes, com `samtools/view` (paquet samtools, comanda view) o `gatk/haplotypecaller` (paquet GATK, comanda HaplotypeCaller).
    Per a eines que només proporcionen una comanda principal, els mòduls utilitzen un sol nivell com `fastqc` o `multiqc`.

Feu clic a la caixa del mòdul `cat_cat` per veure la documentació del mòdul.

La pàgina del mòdul mostra:

- Una breu descripció: "A module for concatenation of gzipped or uncompressed files"
- Comanda d'instal·lació: `nf-core modules install cat/cat`
- Estructura del canal d'entrada i sortida
- Paràmetres disponibles

### 1.2. Llistar els mòduls disponibles des de la línia de comandes

Alternativament, també podeu cercar mòduls directament des de la línia de comandes utilitzant les eines nf-core.

```bash
nf-core modules list remote
```

Això mostrarà una llista de tots els mòduls disponibles al repositori nf-core/modules, tot i que és una mica menys convenient si no coneixeu ja el nom del mòdul que esteu cercant.
No obstant això, si el coneixeu, podeu canalitzar la llista a `grep` per trobar mòduls específics:

```bash
nf-core modules list remote | grep 'cat/cat'
```

??? success "Sortida de la comanda"

    ```console
    │ cat/cat
    ```

Tingueu en compte que l'enfocament amb `grep` només extraurà resultats amb el terme de cerca al seu nom, cosa que no funcionaria per a `cat_cat`.

### 1.3. Obtenir informació detallada sobre el mòdul

Per veure informació detallada sobre un mòdul específic des de la línia de comandes, utilitzeu la comanda `info`:

```bash
nf-core modules info cat/cat
```

Això mostra documentació sobre el mòdul, incloent les seves entrades, sortides i informació bàsica d'ús.

??? success "Sortida de la comanda"

    ```console

                                              ,--./,-.
              ___     __   __   __   ___     /,-._.--~\
        |\ | |__  __ /  ` /  \ |__) |__         }  {
        | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                              `._,._,'

        nf-core/tools version 3.4.1 - https://nf-co.re


    ╭─ Module: cat/cat  ─────────────────────────────────────────────────╮
    │ 🌐 Repository: https://github.com/nf-core/modules.git              │
    │ 🔧 Tools: cat                                                      │
    │ 📖 Description: A module for concatenation of gzipped or           │
    │ uncompressed files                                                 │
    ╰────────────────────────────────────────────────────────────────────╯
                      ╷                                          ╷
    📥 Inputs        │Description                               │Pattern
    ╺━━━━━━━━━━━━━━━━━┿━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┿━━━━━━━╸
    input[0]         │                                          │
    ╶─────────────────┼──────────────────────────────────────────┼───────╴
      meta  (map)     │Groovy Map containing sample information  │
                      │e.g. [ id:'test', single_end:false ]      │
    ╶─────────────────┼──────────────────────────────────────────┼───────╴
      files_in  (file)│List of compressed / uncompressed files   │      *
                      ╵                                          ╵
                          ╷                                 ╷
    📥 Outputs           │Description                      │     Pattern
    ╺━━━━━━━━━━━━━━━━━━━━━┿━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┿━━━━━━━━━━━━╸
    file_out             │                                 │
    ╶─────────────────────┼─────────────────────────────────┼────────────╴
      meta  (map)         │Groovy Map containing sample     │
                          │information                      │
    ╶─────────────────────┼─────────────────────────────────┼────────────╴
      ${prefix}  (file)   │Concatenated file. Will be       │ ${file_out}
                          │gzipped if file_out ends with    │
                          │".gz"                            │
    ╶─────────────────────┼─────────────────────────────────┼────────────╴
    versions             │                                 │
    ╶─────────────────────┼─────────────────────────────────┼────────────╴
      versions.yml  (file)│File containing software versions│versions.yml
                          ╵                                 ╵

    💻  Installation command: nf-core modules install cat/cat

    ```

Aquesta és exactament la mateixa informació que podeu trobar al lloc web.

### 1.4. Instal·lar el mòdul cat/cat

Ara que hem trobat el mòdul que volem, hem d'afegir-lo al codi font del nostre pipeline.

La bona notícia és que el projecte nf-core inclou eines per facilitar aquesta part.
Específicament, la comanda `nf-core modules install` permet automatitzar la recuperació del codi i fer-lo disponible al vostre projecte en un sol pas.

Navegueu al directori del vostre pipeline i executeu la comanda d'instal·lació:

```bash
cd core-hello
nf-core modules install cat/cat
```

L'eina pot primer demanar-vos que especifiqueu un tipus de repositori.
(Si no, salteu a "Finalment, l'eina procedirà a instal·lar el mòdul.")

??? success "Sortida de la comanda"

    ```console

                                          ,--./,-.
          ___     __   __   __   ___     /,-._.--~\
    |\ | |__  __ /  ` /  \ |__) |__         }  {
    | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                          `._,._,'

    nf-core/tools version 3.4.1 - https://nf-co.re


    WARNING  'repository_type' not defined in .nf-core.yml
    ? Is this repository a pipeline or a modules repository? (Use arrow keys)
    » Pipeline
      Modules repository
    ```

Si és així, premeu enter per acceptar la resposta per defecte (`Pipeline`) i continuar.

L'eina oferirà després modificar la configuració del vostre projecte per evitar aquesta pregunta en el futur.

??? success "Sortida de la comanda"

    ```console
        INFO     To avoid this prompt in the future, add the 'repository_type' key to your .nf-core.yml file.
        ? Would you like me to add this config now? [y/n] (y):
    ```

Val la pena aprofitar aquesta eina convenient!
Premeu enter per acceptar la resposta per defecte (sí).

Finalment, l'eina procedirà a instal·lar el mòdul.

??? success "Sortida de la comanda"

    ```console
    INFO Config added to '.nf-core.yml'
    INFO Reinstalling modules found in 'modules.json' but missing from directory:
    INFO Installing 'cat/cat'
    INFO Use the following statement to include this module:

        include { CAT_CAT } from '../modules/nf-core/cat/cat/main'
    ```

La comanda automàticament:

- Descarrega els fitxers del mòdul a `modules/nf-core/cat/cat/`
- Actualitza `modules.json` per fer seguiment del mòdul instal·lat
- Us proporciona la declaració `include` correcta per utilitzar al vostre workflow

!!! tip "Consell"

    Assegureu-vos sempre que el vostre directori de treball actual és l'arrel del vostre projecte de pipeline abans d'executar la comanda d'instal·lació del mòdul.

Comprovem que el mòdul s'ha instal·lat correctament:

```bash
tree -L 4 modules
```

??? abstract "Contingut del directori"

    ```console
    modules
    ├── local
    │   ├── collectGreetings.nf
    │   ├── convertToUpper.nf
    │   ├── cowpy.nf
    │   └── sayHello.nf
    └── nf-core
        └── cat
            └── cat
                ├── environment.yml
                ├── main.nf
                ├── meta.yml
                └── tests

    5 directories, 7 files
    ```

També podeu verificar la instal·lació demanant a la utilitat nf-core que llisti els mòduls instal·lats localment:

```bash
nf-core modules list local
```

??? success "Sortida de la comanda"

    ```console
    INFO     Repository type: pipeline
    INFO     Modules installed in '.':

    ┏━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━┓
    ┃ Module Name ┃ Repository      ┃ Version SHA ┃ Message                                ┃ Date       ┃
    ┡━━━━━━━━━━━━━╇━━━━━━━━━━━━━━━━━╇━━━━━━━━━━━━━╇━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╇━━━━━━━━━━━━┩
    │ cat/cat     │ nf-core/modules │ 41dfa3f     │ update meta.yml of all modules (#8747) │ 2025-07-07 │
    └─────────────┴─────────────────┴─────────────┴────────────────────────────────────────┴────────────┘
    ```

Això confirma que el mòdul `cat/cat` ara forma part del codi font del vostre projecte.

No obstant això, per utilitzar realment el nou mòdul, hem d'importar-lo al nostre pipeline.

### 1.5. Actualitzar les importacions de mòduls

Substituïm la declaració `include` del mòdul `collectGreetings` per la de `CAT_CAT` a la secció d'importacions del workflow `workflows/hello.nf`.

Com a recordatori, l'eina d'instal·lació de mòduls ens va proporcionar la declaració exacta a utilitzar:

```groovy title="Import statement produced by install command"
include { CAT_CAT } from '../modules/nf-core/cat/cat/main'`
```

Tingueu en compte que la convenció nf-core és utilitzar majúscules per als noms de mòduls quan s'importen.

Obriu [core-hello/workflows/hello.nf](core-hello/workflows/hello.nf) i feu la substitució següent:

=== "Després"

    ```groovy title="core-hello/workflows/hello.nf" linenums="1" hl_lines="11"
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
    include { collectGreetings       } from '../modules/local/collectGreetings.nf'
    include { cowpy                  } from '../modules/local/cowpy.nf'
    ```

Observeu com el camí per al mòdul nf-core difereix dels mòduls locals:

- **Mòdul nf-core**: `'../modules/nf-core/cat/cat/main'` (referència a `main.nf`)
- **Mòdul local**: `'../modules/local/collectGreetings.nf'` (referència a un sol fitxer)

El mòdul ara està disponible per al workflow, així que només cal que substituïm la crida a `collectGreetings` per utilitzar `CAT_CAT`. Oi?

No tan ràpid.

En aquest punt, podríeu estar temptats de començar a editar codi, però val la pena prendre's un moment per examinar acuradament què espera el nou mòdul i què produeix.

Tractarem això com una secció separada perquè implica un nou mecanisme que encara no hem cobert: els mapes de metadades.

!!! note "Nota"

    Opcionalment podeu eliminar el fitxer `collectGreetings.nf`:

    ```bash
    rm modules/local/collectGreetings.nf
    ```

    No obstant això, potser voldreu conservar-lo com a referència per entendre les diferències entre mòduls locals i nf-core.

### Conclusió

Sabeu com trobar un mòdul nf-core i fer-lo disponible al vostre projecte.

### Què segueix?

Avaluar què requereix un nou mòdul i identificar qualsevol canvi important necessari per integrar-lo a un pipeline.

---

## 2. Avaluar els requisits del nou mòdul

Específicament, hem d'examinar la **interfície** del mòdul, és a dir, les seves definicions d'entrada i sortida, i comparar-la amb la interfície del mòdul que volem substituir.
Això ens permetrà determinar si podem tractar el nou mòdul com un reemplaçament directe o si haurem d'adaptar part del cablejat.

Idealment això és quelcom que hauríeu de fer _abans_ fins i tot d'instal·lar el mòdul, però bé, més val tard que mai.
(Per cert, hi ha una comanda `uninstall` per desfer-se dels mòduls que decidiu que ja no voleu.)

!!! note "Nota"

    El procés CAT_CAT inclou una gestió força intel·ligent de diferents tipus de compressió, extensions de fitxer i altres aspectes que no són estrictament rellevants per al que intentem mostrar-vos aquí, així que ignorarem la major part i ens centrarem només en les parts que són importants.

### 2.1. Comparar les interfícies dels dos mòduls

Com a recordatori, així és com es veu la interfície del nostre mòdul `collectGreetings`:

```groovy title="modules/local/collectGreetings.nf (excerpt)" linenums="1" hl_lines="6-7 10"
process collectGreetings {

    publishDir 'results', mode: 'copy'

    input:
        path input_files
        val batch_name

    output:
        path "COLLECTED-${batch_name}-output.txt" , emit: outfile
```

El mòdul `collectGreetings` pren dues entrades:

- `input_files` conté un o més fitxers d'entrada a processar;
- `batch_name` és un valor que utilitzem per assignar un nom específic de l'execució al fitxer de sortida, que és una forma de metadades.

En completar-se, `collectGreetings` produeix un sol camí de fitxer, emès amb l'etiqueta `outfile`.

En comparació, la interfície del mòdul `cat/cat` és més complexa:

```groovy title="modules/nf-core/cat/cat/main.nf (excerpt)" linenums="1" hl_lines="11 14"
process CAT_CAT {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pigz:2.3.4' :
        'biocontainers/pigz:2.3.4' }"

    input:
    tuple val(meta), path(files_in)

    output:
    tuple val(meta), path("${prefix}"), emit: file_out
    path "versions.yml"               , emit: versions
```

El mòdul CAT_CAT pren una sola entrada, però aquesta entrada és una tupla que conté dues coses:

- `meta` és una estructura que conté metadades, anomenada metamap;
- `files_in` conté un o més fitxers d'entrada a processar, equivalent a `input_files` de `collectGreetings`.

En completar-se, CAT_CAT lliura les seves sortides en dues parts:

- Una altra tupla que conté el metamap i el fitxer de sortida concatenat, emès amb l'etiqueta `file_out`;
- Un fitxer `versions.yml` que captura informació sobre la versió del programari que s'ha utilitzat, emès amb l'etiqueta `versions`.

Tingueu en compte també que per defecte, el fitxer de sortida es nomenarà basant-se en un identificador que forma part de les metadades (codi no mostrat aquí).

Això pot semblar molt a tenir en compte només mirant el codi, així que aquí teniu un diagrama per ajudar-vos a visualitzar com encaixa tot plegat.

<figure class="excalidraw">
--8<-- "docs/en/docs/hello_nf-core/img/module_comparison.svg"
</figure>

Podeu veure que els dos mòduls tenen requisits d'entrada similars pel que fa al contingut (un conjunt de fitxers d'entrada més algunes metadades) però expectatives molt diferents sobre com s'empaqueta aquest contingut.
Ignorant el fitxer de versions per ara, la seva sortida principal també és equivalent (un fitxer concatenat), excepte que CAT_CAT també emet el metamap conjuntament amb el fitxer de sortida.

Les diferències d'empaquetament seran força fàcils de gestionar, com veureu d'aquí a poc.
No obstant això, per entendre la part del metamap, hem d'introduir-vos a algun context addicional.

### 2.2. Entendre els metamaps

Acabem de dir-vos que el mòdul CAT_CAT espera un mapa de metadades com a part de la seva tupla d'entrada.
Prenguem-nos uns minuts per examinar més de prop què és això.

El **mapa de metadades**, sovint anomenat **metamap** per abreujar, és un mapa d'estil Groovy que conté informació sobre unitats de dades.
En el context dels pipelines Nextflow, les unitats de dades poden ser qualsevol cosa que vulgueu: mostres individuals, lots de mostres o conjunts de dades complets.

Per convenció, un metamap nf-core s'anomena `meta` i conté el camp obligatori `id`, que s'utilitza per nomenar sortides i fer seguiment d'unitats de dades.

Per exemple, un mapa de metadades típic podria tenir aquest aspecte:

```groovy title="Example of sample-level metamap"
[id: 'sample1', single_end: false, strandedness: 'forward']
```

O en un cas on les metadades s'adjunten a nivell de lot:

```groovy title="Example of batch-level metamap"
[id: 'batch1', date: '25.10.01']
```

Ara posem això en el context del procés `CAT_CAT`, que espera que els fitxers d'entrada s'empaqueten en una tupla amb un metamap, i també produeix el metamap com a part de la tupla de sortida.

```groovy title="modules/nf-core/cat/cat/main.nf (excerpt)" linenums="1" hl_lines="2 5"
input:
tuple val(meta), path(files_in)

output:
tuple val(meta), path("${prefix}"), emit: file_out
```

Com a resultat, cada unitat de dades viatja pel pipeline amb les metadades rellevants adjuntes.
Els processos posteriors poden accedir fàcilment a aquestes metadades també.

Recordeu com us vam dir que el fitxer produït per `CAT_CAT` es nomenarà basant-se en un identificador que forma part de les metadades?
Aquest és el codi rellevant:

```groovy title="modules/nf-core/cat/cat/main.nf (excerpt)" linenums="35"
prefix   = task.ext.prefix ?: "${meta.id}${getFileSuffix(file_list[0])}"
```

Això es tradueix aproximadament com segueix: si es proporciona un `prefix` mitjançant el sistema de paràmetres de tasca extern (`task.ext`), utilitzeu-lo per nomenar el fitxer de sortida; en cas contrari, creeu-ne un utilitzant `${meta.id}`, que correspon al camp `id` del metamap.

Podeu imaginar el canal d'entrada que arriba a aquest mòdul amb continguts com aquest:

```groovy title="Example input channel contents"
ch_input = [[[id: 'batch1', date: '25.10.01'], ['file1A.txt', 'file1B.txt']],
            [[id: 'batch2', date: '25.10.26'], ['file2A.txt', 'file2B.txt']],
            [[id: 'batch3', date: '25.11.14'], ['file3A.txt', 'file3B.txt']]]
```

Llavors el contingut del canal de sortida que surt seria així:

```groovy title="Example output channel contents"
ch_input = [[[id: 'batch1', date: '25.10.01'], 'batch1.txt'],
            [[id: 'batch2', date: '25.10.26'], 'batch2.txt'],
            [[id: 'batch3', date: '25.11.14'], 'batch3.txt']]
```

Com s'ha esmentat anteriorment, la configuració d'entrada `tuple val(meta), path(files_in)` és un patró estàndard utilitzat a tots els mòduls nf-core.

Esperem que pugueu començar a veure com d'útil pot ser això.
No només us permet nomenar sortides basant-vos en metadades, sinó que també podeu fer coses com utilitzar-les per aplicar diferents valors de paràmetres, i en combinació amb operadors específics, fins i tot podeu agrupar, ordenar o filtrar dades mentre flueixen pel pipeline.

!!! note "Més informació sobre metadades"

    Per a una introducció completa sobre com treballar amb metadades als workflows Nextflow, incloent com llegir metadades des de fulls de mostres i utilitzar-les per personalitzar el processament, consulteu la missió secundària [Metadades als workflows](../side_quests/metadata).

### 2.3. Resumir els canvis a fer

Basant-nos en el que hem revisat, aquests són els canvis principals que hem de fer al nostre pipeline per utilitzar el mòdul `cat/cat`:

- Crear un metamap que contingui el nom del lot;
- Empaquetar el metamap en una tupla amb el conjunt de fitxers d'entrada a concatenar (provinents de `convertToUpper`);
- Canviar la crida de `collectGreetings()` a `CAT_CAT`;
- Extreure el fitxer de sortida de la tupla produïda pel procés `CAT_CAT` abans de passar-lo a `cowpy`.

Això hauria de ser suficient! Ara que tenim un pla, estem preparats per començar.

### Conclusió

Sabeu com avaluar la interfície d'entrada i sortida d'un nou mòdul per identificar els seus requisits, i heu après com els metamaps són utilitzats pels pipelines nf-core per mantenir les metadades estretament associades amb les dades mentre flueixen pel pipeline.

### Què segueix?

Integrar el nou mòdul en un workflow.

---

## 3. Integrar CAT_CAT al workflow `hello.nf`

Ara que sabeu tot sobre els metamaps (o prou per als propòsits d'aquest curs, almenys), és hora d'implementar realment els canvis que hem descrit anteriorment.

Per claredat, dividirem això i cobrirem cada pas per separat.

!!! note "Nota"

    Tots els canvis mostrats a continuació es fan a la lògica del workflow al bloc `main` del fitxer de workflow `core-hello/workflows/hello.nf`.

### 3.1. Crear un mapa de metadades

Primer, hem de crear un mapa de metadades per a `CAT_CAT`, tenint en compte que els mòduls nf-core requereixen que el metamap tingui almenys un camp `id`.

Com que no necessitem cap altra metadada, podem mantenir-ho simple i utilitzar quelcom així:

```groovy title="Syntax example"
def cat_meta = [id: 'test']
```

Excepte que no volem codificar el valor `id` de forma fixa; volem utilitzar el valor del parametre `params.batch`.
Així que el codi es converteix en:

```groovy title="Syntax example"
def cat_meta = [id: params.batch]
```

Sí, és literalment així de simple crear un metamap bàsic.

Afegim aquestes línies després de la crida a `convertToUpper`, eliminant la crida a `collectGreetings`:

=== "Després"

    ```groovy title="core-hello/workflows/hello.nf" linenums="26" hl_lines="7-8"
        // emet una salutacio
        sayHello(ch_samplesheet)

        // converteix la salutacio a majuscules
        convertToUpper(sayHello.out)

        // crea un mapa de metadades amb el nom del lot com a ID
        def cat_meta = [ id: params.batch ]

        // genera art ASCII de les salutacions amb cowpy
        cowpy(collectGreetings.out.outfile, params.character)
    ```

=== "Abans"

    ```groovy title="core-hello/workflows/hello.nf" linenums="26" hl_lines="7-8"
        // emet una salutacio
        sayHello(ch_samplesheet)

        // converteix la salutacio a majuscules
        convertToUpper(sayHello.out)

        // recull totes les salutacions en un fitxer
        collectGreetings(convertToUpper.out.collect(), params.batch)

        // genera art ASCII de les salutacions amb cowpy
        cowpy(collectGreetings.out.outfile, params.character)
    ```

Això crea un mapa de metadades simple on l'`id` s'estableix al nom del nostre lot (que serà `test` quan s'utilitzi el perfil de test).

### 3.2. Crear un canal amb tuples de metadades

A continuació, transformeu el canal de fitxers en un canal de tuples que continguin metadades i fitxers:

=== "Després"

    ```groovy title="core-hello/workflows/hello.nf" linenums="26" hl_lines="10-11"
        // emet una salutacio
        sayHello(ch_samplesheet)

        // converteix la salutacio a majuscules
        convertToUpper(sayHello.out)

        // crea un mapa de metadades amb el nom del lot com a ID
        def cat_meta = [ id: params.batch ]

        // crea un canal amb metadades i fitxers en format tupla
        ch_for_cat = convertToUpper.out.collect().map { files -> tuple(cat_meta, files) }

        // genera art ASCII de les salutacions amb cowpy
        cowpy(collectGreetings.out.outfile, params.character)
    ```

=== "Abans"

    ```groovy title="core-hello/workflows/hello.nf" linenums="26"
        // emet una salutacio
        sayHello(ch_samplesheet)

        // converteix la salutacio a majuscules
        convertToUpper(sayHello.out)

        // crea un mapa de metadades amb el nom del lot com a ID
        def cat_meta = [ id: params.batch ]

        // genera art ASCII de les salutacions amb cowpy
        cowpy(collectGreetings.out.outfile, params.character)
    ```

La línia que hem afegit aconsegueix dues coses:

- `.collect()` recull tots els fitxers de la sortida de `convertToUpper` en una sola llista
- `.map { files -> tuple(cat_meta, files) }` crea una tupla de `[metadades, fitxers]` en el format que `CAT_CAT` espera

Això és tot el que necessitem fer per configurar la tupla d'entrada per a `CAT_CAT`.

### 3.3. Cridar el mòdul CAT_CAT

Ara cridem `CAT_CAT` sobre el canal acabat de crear:

=== "Després"

    ```groovy title="core-hello/workflows/hello.nf" linenums="26" hl_lines="13-14"
        // emet una salutacio
        sayHello(ch_samplesheet)

        // converteix la salutacio a majuscules
        convertToUpper(sayHello.out)

        // crea un mapa de metadades amb el nom del lot com a ID
        def cat_meta = [ id: params.batch ]

        // crea un canal amb metadades i fitxers en format tupla
        ch_for_cat = convertToUpper.out.collect().map { files -> tuple(cat_meta, files) }

        // concatena els fitxers utilitzant el modul nf-core cat/cat
        CAT_CAT(ch_for_cat)

        // genera art ASCII de les salutacions amb cowpy
        cowpy(collectGreetings.out.outfile, params.character)
    ```

=== "Abans"

    ```groovy title="core-hello/workflows/hello.nf" linenums="26"
        // emet una salutacio
        sayHello(ch_samplesheet)

        // converteix la salutacio a majuscules
        convertToUpper(sayHello.out)

        // crea un mapa de metadades amb el nom del lot com a ID
        def cat_meta = [ id: params.batch ]

        // crea un canal amb metadades i fitxers en format tupla
        ch_for_cat = convertToUpper.out.collect().map { files -> tuple(cat_meta, files) }

        // genera art ASCII de les salutacions amb cowpy
        cowpy(collectGreetings.out.outfile, params.character)
    ```

Això completa la part més complicada d'aquesta substitució, però encara no hem acabat del tot: encara hem d'actualitzar com passem la sortida concatenada al procés `cowpy`.

### 3.4. Extreure el fitxer de sortida de la tupla per a `cowpy`

Anteriorment, el procés `collectGreetings` només produïa un fitxer que podíem passar directament a `cowpy`.
No obstant això, el procés `CAT_CAT` produeix una tupla que inclou el metamap a més del fitxer de sortida.

Com que `cowpy` encara no accepta tuples de metadades (ho arreglarem a la següent part del curs), hem d'extreure el fitxer de sortida de la tupla produïda per `CAT_CAT` abans de passar-lo a `cowpy`:

=== "Després"

    ```groovy title="core-hello/workflows/hello.nf" linenums="26" hl_lines="16-17 20"
        // emet una salutacio
        sayHello(ch_samplesheet)

        // converteix la salutacio a majuscules
        convertToUpper(sayHello.out)

        // crea un mapa de metadades amb el nom del lot com a ID
        def cat_meta = [ id: params.batch ]

        // crea un canal amb metadades i fitxers en format tupla
        ch_for_cat = convertToUpper.out.collect().map { files -> tuple(cat_meta, files) }

        // concatena les salutacions
        CAT_CAT(ch_for_cat)

        // extreu el fitxer de la tupla ja que cowpy encara no utilitza metadades
        ch_for_cowpy = CAT_CAT.out.file_out.map{ meta, file -> file }

        // genera art ASCII de les salutacions amb cowpy
        cowpy(ch_for_cowpy, params.character)
    ```

=== "Abans"

    ```groovy title="core-hello/workflows/hello.nf" linenums="26" hl_lines="17"
        // emet una salutacio
        sayHello(ch_samplesheet)

        // converteix la salutacio a majuscules
        convertToUpper(sayHello.out)

        // crea un mapa de metadades amb el nom del lot com a ID
        def cat_meta = [ id: params.batch ]

        // crea un canal amb metadades i fitxers en format tupla
        ch_for_cat = convertToUpper.out.collect().map { files -> tuple(cat_meta, files) }

        // concatena les salutacions
        CAT_CAT(ch_for_cat)

        // genera art ASCII de les salutacions amb cowpy
        cowpy(collectGreetings.out.outfile, params.character)
    ```

L'operació `.map{ meta, file -> file }` extreu el fitxer de la tupla `[metadades, fitxer]` produïda per `CAT_CAT` en un nou canal, `ch_for_cowpy`.

Llavors només cal passar `ch_for_cowpy` a `cowpy` en lloc de `collectGreetings.out.outfile` en aquesta última línia.

!!! note "Nota"

    A la següent part del curs, actualitzarem `cowpy` per treballar directament amb tuples de metadades, així que aquest pas d'extracció ja no serà necessari.

### 3.5. Provar el workflow

Provem que el workflow funciona amb el mòdul `cat/cat` acabat d'integrar:

```bash
nextflow run . --outdir core-hello-results -profile test,docker --validate_params false
```

Això hauria d'executar-se raonablement ràpid.

??? success "Sortida de la comanda"

    ```console
    N E X T F L O W ~ version 25.04.3

        Launching `./main.nf` [evil_pike] DSL2 - revision: b9e9b3b8de

        Input/output options
          input                     : /workspaces/training/hello-nf-core/core-hello/assets/greetings.csv
          outdir                    : core-hello-results

        Institutional config options
          config_profile_name       : Test profile
          config_profile_description: Minimal test dataset to check pipeline function

        Generic options
          validate_params           : false
          trace_report_suffix       : 2025-10-30_18-50-58

        Core Nextflow options
          runName                   : evil_pike
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
        [b3/f005fd] CORE_HELLO:HELLO:sayHello (3)       [100%] 3 of 3 ✔
        [08/f923d0] CORE_HELLO:HELLO:convertToUpper (3) [100%] 3 of 3 ✔
        [34/3729a9] CORE_HELLO:HELLO:CAT_CAT (test)     [100%] 1 of 1 ✔
        [24/df918a] CORE_HELLO:HELLO:cowpy              [100%] 1 of 1 ✔
        -[core/hello] Pipeline completed successfully-
    ```

Observeu que `CAT_CAT` ara apareix a la llista d'execució de processos en lloc de `collectGreetings`.

I això és tot! Ara estem utilitzant un mòdul robust mantingut per la comunitat en lloc de codi personalitzat de qualitat prototip per a aquest pas del pipeline.

### Conclusió

Ara sabeu com:

- Trobar i instal·lar mòduls nf-core
- Avaluar els requisits d'un mòdul nf-core
- Crear un mapa de metadades simple per utilitzar amb un mòdul nf-core
- Integrar un mòdul nf-core al vostre workflow

### Què segueix?

Aprendre a adaptar els vostres mòduls locals per seguir les convencions nf-core.
També us mostrarem com crear nous mòduls nf-core a partir d'una plantilla utilitzant les eines nf-core.
