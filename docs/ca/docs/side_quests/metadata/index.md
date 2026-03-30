# Metadades i mapes de metadades

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducció assistida per IA - [més informació i suggeriments](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

En qualsevol anàlisi científica, rarament treballem només amb els fitxers de dades en brut.
Cada fitxer porta la seva pròpia informació addicional: què és, d'on prové i què el fa especial.
Aquesta informació extra és el que anomenem metadades.

Les metadades són dades que descriuen altres dades.
Les metadades fan un seguiment de detalls importants sobre els fitxers i les condicions experimentals, i ajuden a adaptar les anàlisis a les característiques úniques de cada conjunt de dades.

Penseu-hi com en un catàleg de biblioteca: mentre que els llibres contenen el contingut real (les dades en brut), les fitxes del catàleg proporcionen informació essencial sobre cada llibre—quan es va publicar, qui el va escriure, on trobar-lo (metadades).
En els pipelines de Nextflow, les metadades es poden utilitzar per a:

- Fer un seguiment de la informació específica de cada fitxer al llarg del workflow
- Configurar processos en funció de les característiques dels fitxers
- Agrupar fitxers relacionats per a una anàlisi conjunta

### Objectius d'aprenentatge

En aquesta missió secundària, explorarem com gestionar metadades en workflows.
Partint d'un full de dades senzill (sovint anomenat samplesheet en bioinformàtica) que conté informació bàsica sobre els fitxers, aprendreu com:

- Llegir i analitzar metadades de fitxers CSV
- Crear i manipular mapes de metadades
- Afegir nous camps de metadades durant l'execució del workflow
- Utilitzar metadades per personalitzar el comportament dels processos

Aquestes habilitats us ajudaran a construir pipelines més robustos i flexibles que puguin gestionar relacions complexes entre fitxers i requisits de processament.

### Prerequisits

Abans d'emprendre aquesta missió secundària, hauríeu de:

- Haver completat el tutorial [Hello Nextflow](../hello_nextflow/README.md) o un curs equivalent per a principiants.
- Estar còmodes amb els conceptes i mecanismes bàsics de Nextflow (processos, canals, operadors)

---

## 0. Primers passos

#### Obriu el codespace de formació

Si encara no ho heu fet, assegureu-vos d'obrir l'entorn de formació tal com es descriu a la [Configuració de l'entorn](../envsetup/index.md).

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### Moveu-vos al directori del projecte

Anem al directori on es troben els fitxers d'aquest tutorial.

```bash
cd side-quests/metadata
```

Podeu configurar VSCode perquè es centri en aquest directori:

```bash
code .
```

#### Reviseu els materials

Trobareu un fitxer de workflow principal i un directori `data` que conté un full de dades i uns quants fitxers de dades.

??? abstract "Contingut del directori"

    ```console
    .
    ├── data
    │   ├── bonjour.txt
    │   ├── ciao.txt
    │   ├── guten_tag.txt
    │   ├── hallo.txt
    │   ├── hello.txt
    │   ├── hola.txt
    │   ├── salut.txt
    │   └── datasheet.csv
    ├── main.nf
    └── nextflow.config
    ```

El workflow del fitxer `main.nf` és un esquelet que anireu ampliant gradualment fins a convertir-lo en un workflow completament funcional.

El full de dades llista les rutes als fitxers de dades i algunes metadades associades, organitzades en 3 columnes:

- `id`: autoexplicatiu, un identificador assignat al fitxer
- `character`: un nom de personatge, que utilitzarem més endavant per dibuixar diferents criatures
- `data`: rutes als fitxers `.txt` que contenen salutacions en diferents idiomes

```console title="datasheet.csv"
id,character,recording
sampleA,squirrel,/workspaces/training/side-quests/metadata/data/bonjour.txt
sampleB,tux,/workspaces/training/side-quests/metadata/data/guten_tag.txt
sampleC,sheep,/workspaces/training/side-quests/metadata/data/hallo.txt
sampleD,turkey,/workspaces/training/side-quests/metadata/data/hello.txt
sampleE,stegosaurus,/workspaces/training/side-quests/metadata/data/hola.txt
sampleF,moose,/workspaces/training/side-quests/metadata/data/salut.txt
sampleG,turtle,/workspaces/training/side-quests/metadata/data/ciao.txt
```

Cada fitxer de dades conté text de salutació en un dels cinc idiomes (fr: francès, de: alemany, es: espanyol, it: italià, en: anglès).

També us proporcionarem una eina d'anàlisi de llenguatge en contenidor anomenada `langid`.

#### Reviseu l'assignació

El vostre repte és escriure un workflow de Nextflow que:

1. **Identifiqui** automàticament l'idioma de cada fitxer
2. **Agrupe** els fitxers per família lingüística (llengües germàniques vs. romàniques)
3. **Personalitzi** el processament de cada fitxer en funció del seu idioma i metadades
4. **Organitzi** les sortides per grup lingüístic

Això representa un patró típic de workflow on les metadades específiques de cada fitxer guien les decisions de processament; exactament el tipus de problema que els mapes de metadades resolen de manera elegant.

#### Llista de comprovació de preparació

Creieu que esteu a punt per submergir-vos?

- [ ] Entenc l'objectiu d'aquest curs i els seus prerequisits
- [ ] El meu codespace està en funcionament
- [ ] He configurat el meu directori de treball adequadament
- [ ] Entenc l'assignació

Si podeu marcar totes les caselles, esteu a punt per començar.

---

## 1. Carregar metadades des d'un full de dades

Obriu el fitxer de workflow `main.nf` per examinar l'esquelet del workflow que us donem com a punt de partida.

```groovy title="main.nf" linenums="1"
#!/usr/bin/env nextflow

workflow  {

    ch_datasheet = channel.fromPath("./data/datasheet.csv")

}
```

Podeu veure que hem configurat una factoria de canals bàsica per carregar el full de dades d'exemple com a fitxer, però que encara no llegirà el contingut del fitxer.
Comencem afegint-ho.

### 1.1. Llegir el contingut amb `splitCsv`

Hem de triar un operador que analitzi el contingut del fitxer de manera adequada amb el mínim esforç per la nostra part.
Com que el nostre full de dades és en format CSV, aquesta és una feina per a l'operador [`splitCsv`](https://www.nextflow.io/docs/latest/reference/operator.html#splitcsv), que carrega cada fila del fitxer com un element del canal.

Feu els canvis següents per afegir una operació `splitCsv()` al codi de construcció del canal, més una operació `view()` per comprovar que el contingut del fitxer s'està carregant correctament al canal.

=== "Després"

    ```groovy title="main.nf" linenums="3" hl_lines="4-5"
    workflow  {

        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .view()

    }
    ```

=== "Abans"

    ```groovy title="main.nf" linenums="3"
    workflow  {

        ch_datasheet = channel.fromPath("./data/datasheet.csv")

    }
    ```

Tingueu en compte que estem utilitzant l'opció `header: true` per indicar a Nextflow que llegeixi la primera fila del fitxer CSV com a fila de capçalera.

Vegem què en surt, oi?
Executeu el workflow:

```bash
nextflow run main.nf
```

??? success "Sortida de la comanda"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [exotic_albattani] DSL2 - revision: c0d03cec83

    [id:sampleA, character:squirrel, recording:/workspaces/training/side-quests/metadata/data/bonjour.txt]
    [id:sampleB, character:tux, recording:/workspaces/training/side-quests/metadata/data/guten_tag.txt]
    [id:sampleC, character:sheep, recording:/workspaces/training/side-quests/metadata/data/hallo.txt]
    [id:sampleD, character:turkey, recording:/workspaces/training/side-quests/metadata/data/hello.txt]
    [id:sampleE, character:stegosaurus, recording:/workspaces/training/side-quests/metadata/data/hola.txt]
    [id:sampleF, character:moose, recording:/workspaces/training/side-quests/metadata/data/salut.txt]
    [id:sampleG, character:turtle, recording:/workspaces/training/side-quests/metadata/data/ciao.txt]
    ```

Podem veure que l'operador ha construït un mapa de parells clau-valor per a cada fila del fitxer CSV, amb les capçaleres de columna com a claus per als valors corresponents.

Cada entrada del mapa correspon a una columna del nostre full de dades:

- `id`
- `character`
- `recording`

Molt bé! Això facilita l'accés a camps específics de cada fitxer.
Per exemple, podríem accedir a l'identificador del fitxer amb `id` o a la ruta del fitxer txt amb `recording`.

??? info "(Opcional) Més informació sobre els mapes"

    En Groovy, el llenguatge de programació sobre el qual es construeix Nextflow, un mapa és una estructura de dades de clau-valor similar als diccionaris en Python, els objectes en JavaScript o els hashes en Ruby.

    Aquí teniu un script executable que mostra com podeu definir un mapa i accedir al seu contingut en la pràctica:

    ```groovy title="examples/map_demo.nf"
    #!/usr/bin/env nextflow

    // Crea un mapa senzill
    def my_map = [id:'sampleA', character:'squirrel']

    // Imprimeix el mapa sencer
    println "map: ${my_map}"

    // Accedeix als valors individuals amb notació de punt
    println "id: ${my_map.id}"
    println "character: ${my_map.character}"
    ```

    Tot i que no té un bloc `workflow` pròpiament dit, Nextflow pot executar-lo com si fos un workflow:

    ```bash
    nextflow run examples/map_demo.nf
    ```

    I aquí teniu el que podeu esperar veure a la sortida:

    ```console title="Output"
     N E X T F L O W   ~  version 25.10.2

    Launching `map_demo.nf` [cheesy_plateau] DSL2 - revision: fae5b8496e

    map: [id:sampleA, character:squirrel]
    id: sampleA
    character: squirrel
    ```

### 1.2. Seleccionar camps específics amb `map`

Suposem que volem accedir a la columna `character` del full de dades i imprimir-la.
Podem utilitzar l'operador `map` de Nextflow per iterar sobre cada element del nostre canal i seleccionar específicament l'entrada `character` de l'objecte mapa.

Feu les edicions següents al workflow:

=== "Després"

    ```groovy title="main.nf" linenums="3" hl_lines="5-7"
    workflow  {

        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map{ row ->
                row.character
            }
            .view()

    }
    ```

=== "Abans"

    ```groovy title="main.nf" linenums="3"
    workflow  {

        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .view()

    }
    ```

Ara executeu el workflow de nou:

```bash
nextflow run main.nf
```

??? success "Sortida de la comanda"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [exotic_albattani] DSL2 - revision: c0d03cec83

    squirrel
    tux
    sheep
    turkey
    stegosaurus
    moose
    turtle
    ```

Èxit! Hem aprofitat l'estructura de mapa derivada del nostre full de dades per accedir als valors de columnes individuals per a cada fila.

Ara que hem llegit correctament el full de dades i tenim accés a les dades de cada fila, podem començar a implementar la lògica del nostre pipeline.

### 1.3. Organitzar les metadades en un 'meta map'

En l'estat actual del workflow, els fitxers d'entrada (sota la clau `recording`) i les metadades associades (`id`, `character`) estan tots al mateix nivell, com si estiguessin tots en una gran bossa.
La conseqüència pràctica és que cada procés que consumeixi aquest canal hauria de ser configurat tenint en compte aquesta estructura:

```groovy
    input:
    tuple val(id), val(character), file(recording)
```

Això és correcte sempre que el nombre de columnes del full de dades no canviï.
Tanmateix, si afegiu fins i tot una sola columna al full de dades, la forma del canal ja no coincidirà amb el que espera el procés, i el workflow produirà errors.
També fa que el procés sigui difícil de compartir amb altres persones que puguin tenir dades d'entrada lleugerament diferents, i podríeu acabar havent de codificar variables al procés que no són necessàries per al bloc script.

Per evitar aquest problema, hem de trobar una manera de mantenir l'estructura del canal consistent independentment del nombre de columnes que contingui el full de dades.

Podem fer-ho recollint totes les metadades en un element dins de la tupla, que anomenarem el mapa de metadades, o més simplement 'meta map'.

Feu les edicions següents a l'operació `map`:

=== "Després"

    ```groovy title="main.nf" linenums="5" hl_lines="4"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                    [ [id: row.id, character: row.character], row.recording ]
            }
            .view()
    ```

=== "Abans"

    ```groovy title="main.nf" linenums="5" hl_lines="4"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map{ row ->
                row.character
            }
            .view()
    ```

Hem reestructurat els elements del canal en una tupla formada per dos elements, el meta map i l'objecte fitxer corresponent.

Executem el workflow:

```bash
nextflow run main.nf
```

??? success "Sortida de la comanda"

    ```console title="View meta map"
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [lethal_booth] DSL2 - revision: 0d8f844c07

    [[id:sampleA, character:squirrel], /workspaces/training/side-quests/metadata/data/bonjour.txt]
    [[id:sampleB, character:tux], /workspaces/training/side-quests/metadata/data/guten_tag.txt]
    [[id:sampleC, character:sheep], /workspaces/training/side-quests/metadata/data/hallo.txt]
    [[id:sampleD, character:turkey], /workspaces/training/side-quests/metadata/data/hello.txt]
    [[id:sampleE, character:stegosaurus], /workspaces/training/side-quests/metadata/data/hola.txt]
    [[id:sampleF, character:moose], /workspaces/training/side-quests/metadata/data/salut.txt]
    [[id:sampleG, character:turtle], /workspaces/training/side-quests/metadata/data/ciao.txt]
    ```

Ara, cada element del canal conté primer el mapa de metadades i segon l'objecte fitxer corresponent:

```console title="Example output structure"
[
  [id:sampleA, character:squirrel],
  /workspaces/training/side-quests/metadata/data/bonjour.txt
]
```

Com a resultat, afegir més columnes al full de dades farà que hi hagi més metadades disponibles al mapa `meta`, però no canviarà la forma del canal.
Això ens permet escriure processos que consumeixin el canal sense haver de codificar els elements de metadades a l'especificació d'entrada:

```groovy title="Syntax example"
    input:
    tuple val(meta), file(recording)
```

Aquesta és una convenció àmpliament utilitzada per organitzar metadades en workflows de Nextflow.

### Conclusió

En aquesta secció, heu après:

- **Per què les metadades són importants:** Mantenir les metadades amb les vostres dades preserva informació important sobre els fitxers al llarg del workflow.
- **Com llegir fulls de dades:** Utilitzant `splitCsv` per llegir fitxers CSV amb informació de capçalera i transformar les files en dades estructurades
- **Com crear un meta map:** Separant les metadades de les dades dels fitxers utilitzant l'estructura de tupla `[ [id:valor, ...], fitxer ]`

---

## 2. Manipulació de metadades

Ara que tenim les metadades carregades, fem-ne alguna cosa!

Utilitzarem una eina anomenada [`langid`](https://github.com/saffsd/langid.py) per identificar l'idioma contingut al fitxer de gravació de cada criatura.
L'eina ve preentrenada en un conjunt d'idiomes i, donat un fragment de text, produirà una predicció d'idioma i una puntuació de probabilitat associada, tots dos a `stdout`.

### 2.1. Importar el procés i examinar el codi

Us proporcionem un mòdul de procés preescrit anomenat `IDENTIFY_LANGUAGE` que encapsula l'eina `langid`, de manera que només cal afegir una instrucció include abans del bloc workflow.

Feu la següent edició al workflow:

=== "Després"

    ```groovy title="main.nf" linenums="1" hl_lines="3"
    #!/usr/bin/env nextflow

    include { IDENTIFY_LANGUAGE } from './modules/langid.nf'

    workflow {
    ```

=== "Abans"

    ```groovy title="main.nf" linenums="1"
    #!/usr/bin/env nextflow

    workflow {
    ```

Podeu obrir el fitxer del mòdul per examinar el seu codi:

```groovy title="modules/langid.nf" linenums="1" hl_lines="9 12"
#!/usr/bin/env nextflow

// Utilitza langid per predir l'idioma de cada fitxer d'entrada
process IDENTIFY_LANGUAGE {

    container 'community.wave.seqera.io/library/pip_langid:b2269f456a5629ff'

    input:
    tuple val(meta), path(file)

    output:
    tuple val(meta), path(file), stdout

    script:
    """
    langid < ${file} -l en,de,fr,es,it | sed -E "s/.*\\('([a-z]+)'.*/\\1/" | tr -d '\\n'
    """
}
```

Com podeu veure, la definició d'entrada utilitza la mateixa estructura `tuple val(meta), path(file)` que acabem d'aplicar al nostre canal d'entrada.

La definició de sortida està estructurada com una tupla amb una estructura similar a la de l'entrada, excepte que també conté `stdout` com a tercer element.
Aquest patró `tuple val(meta), path(file), <sortida>` manté les metadades associades tant amb les dades d'entrada com amb les sortides a mesura que flueixen pel pipeline.

Tingueu en compte que estem utilitzant el qualificador de sortida [`stdout`](https://www.nextflow.io/docs/latest/process.html#outputs) de Nextflow perquè l'eina imprimeix la seva sortida directament a la consola en lloc d'escriure un fitxer; i utilitzem `sed` a la línia de comandes per eliminar la puntuació de probabilitat, netejar la cadena eliminant els caràcters de nova línia i retornar només la predicció d'idioma.

### 2.2. Afegir una crida a `IDENTIFY_LANGUAGE`

Ara que el procés està disponible per al workflow, podem afegir una crida al procés `IDENTIFY_LANGUAGE` per executar-lo sobre el canal de dades.

Feu les edicions següents al workflow:

=== "Després"

    ```groovy title="main.nf" linenums="7" hl_lines="7-9"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                [[id: row.id, character: row.character], row.recording]
            }

        // Executa langid per identificar l'idioma de cada salutació
        IDENTIFY_LANGUAGE(ch_datasheet)
        IDENTIFY_LANGUAGE.out.view()
    ```

=== "Abans"

    ```groovy title="main.nf" linenums="7" hl_lines="6"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                [[id: row.id, character: row.character], row.recording]
            }
            .view()
    ```

Tingueu en compte que hem eliminat l'operació `.view()` original de la construcció del canal.

Ara podem executar el workflow:

```bash
nextflow run main.nf
```

??? success "Sortida de la comanda"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [voluminous_mcnulty] DSL2 - revision: f9bcfebabb

    executor >  local (7)
    [4e/f722fe] IDENTIFY_LANGUAGE (7) [100%] 7 of 7 ✔
    [[id:sampleA, character:squirrel], /workspaces/training/side-quests/metadata/work/eb/f7148ebdd898fbe1136bec6a714acb/bonjour.txt, fr]
    [[id:sampleB, character:tux], /workspaces/training/side-quests/metadata/work/16/71d72410952c22cd0086d9bca03680/guten_tag.txt, de]
    [[id:sampleD, character:turkey], /workspaces/training/side-quests/metadata/work/c4/b7562adddc1cc0b7d414ec45d436eb/hello.txt, en]
    [[id:sampleC, character:sheep], /workspaces/training/side-quests/metadata/work/ea/04f5d979429e4455e14b9242fb3b45/hallo.txt, de]
    [[id:sampleF, character:moose], /workspaces/training/side-quests/metadata/work/5a/6c2b84bf8fadb98e28e216426be079/salut.txt, fr]
    [[id:sampleE, character:stegosaurus], /workspaces/training/side-quests/metadata/work/af/ee7c69bcab891c40d0529305f6b9e7/hola.txt, es]
    [[id:sampleG, character:turtle], /workspaces/training/side-quests/metadata/work/4e/f722fe47271ba7ebcd69afa42964ca/ciao.txt, it]
    ```

Excel·lent! Ara tenim una predicció de quin idioma parla cada personatge.

I tal com s'ha assenyalat anteriorment, també hem inclòs el fitxer d'entrada i el meta map a la sortida, la qual cosa significa que tots dos romanen associats amb la nova informació que acabem de produir.
Això serà útil al pas següent.

!!! note "Nota"

    De manera més general, aquest patró de mantenir el meta map associat amb els resultats facilita l'associació de resultats relacionats que comparteixen els mateixos identificadors.

    Com ja haureu après, no podeu confiar en l'ordre dels elements dels canals per fer coincidir els resultats entre ells.
    En canvi, heu d'utilitzar claus per associar les dades correctament, i els meta maps proporcionen una estructura ideal per a aquest propòsit.

    Explorem aquest cas d'ús en detall a la missió secundària [Splitting & Grouping](../splitting_and_grouping/).

### 2.3. Augmentar les metadades amb les sortides del procés

Donat que els resultats que acabem de produir són en si mateixos una forma de metadades sobre el contingut dels fitxers, seria útil afegir-los al nostre meta map.

Tanmateix, no volem modificar el meta map existent in situ.
Des d'un punt de vista tècnic, és _possible_ fer-ho, però no és segur.

Per tant, en lloc d'això, crearem un nou meta map que contingui el contingut del meta map existent més un nou parell clau-valor `lang: lang_id` que conté la nova informació, utilitzant l'operador `+` (una característica de Groovy).
I ho combinarem amb una operació [`map`](https://www.nextflow.io/docs/latest/operator.html#map) per substituir el mapa antic pel nou.

Aquí teniu les edicions que heu de fer al workflow:

=== "Després"

    ```groovy title="main.nf" linenums="13" hl_lines="3-7"
        // Executa langid per identificar l'idioma de cada salutació
        IDENTIFY_LANGUAGE(ch_datasheet)
        IDENTIFY_LANGUAGE.out
            .map { meta, file, lang_id ->
                [meta + [lang: lang_id], file]
            }
            .view()
    ```

=== "Abans"

    ```groovy title="main.nf" linenums="13" hl_lines="3"
        // Executa langid per identificar l'idioma de cada salutació
        IDENTIFY_LANGUAGE(ch_datasheet)
        IDENTIFY_LANGUAGE.out.view()
    ```

Si encara no esteu familiaritzats amb l'operador `+`, o si això us sembla confús, dediqueu uns minuts a llegir l'explicació detallada que trobareu a continuació.

??? info "Creació del nou meta map utilitzant l'operador `+`"

    **Primer, heu de saber que podem fusionar el contingut de dos mapes utilitzant l'operador Groovy `+`.**

    Suposem que tenim els mapes següents:

    ```groovy
    map1 = [id: 'sampleA', character: 'squirrel']
    map2 = [lang: 'fr']
    ```

    Podem fusionar-los així:

    ```groovy
    new_map = map1 + map2
    ```

    El contingut de `new_map` serà:

    ```groovy
    [id: 'sampleA', character: 'squirrel', lang: 'fr']
    ```

    Molt bé!

    **Però, i si heu d'afegir un camp que encara no forma part d'un mapa?**

    Suposem que torneu a partir de `map1`, però la predicció d'idioma no és al seu propi mapa (no hi ha cap `map2`).
    En canvi, es troba en una variable anomenada `lang_id`, i sabeu que voleu emmagatzemar el seu valor (`'fr'`) amb la clau `lang`.

    En realitat podeu fer el següent:

    ```groovy
    new_map = [map1 + [lang: lang_id]]
    ```

    Aquí, `[lang: new_info]` crea un nou mapa sense nom al vol, i `map1 + ` fusiona `map1` amb el nou mapa sense nom, produint el mateix contingut de `new_map` que abans.

    Elegant, oi?

    **Ara transposem-ho al context d'una operació `channel.map()` de Nextflow.**

    El codi es converteix en:

    ```groovy
    .map { map1, lang_id ->
        [map1 + [lang: lang_id]]
    }
    ```

    Això fa el següent:

    - `map1, lang_id ->` pren els dos elements de la tupla
    - `[map1 + [lang: lang_id]]` crea el nou mapa tal com s'ha detallat anteriorment

    La sortida és un únic mapa sense nom amb el mateix contingut que `new_map` en el nostre exemple anterior.
    Per tant, hem transformat efectivament:

    ```groovy
    [id: 'sampleA', character: 'squirrel'], 'fr'
    ```

    en:

    ```groovy
    [id: 'sampleA', character: 'squirrel', lang: 'fr']
    ```

    Esperem que pugueu veure que si canviem `map1` per `meta`, això és bàsicament tot el que necessitem per afegir la predicció d'idioma al nostre meta map en el nostre workflow.

    Excepte per una cosa!

    En el cas del nostre workflow, **també hem de tenir en compte la presència de l'objecte `file` a la tupla**, que es compon de `meta, file, lang_id`.

    Per tant, el codi aquí es convertiria en:

    ```groovy
    .map { meta, file, lang_id ->
        [meta + [lang: lang_id], file]
    }
    ```

    Si us costa seguir per què el `file` sembla que es mou a l'operació `map`, imagineu que en lloc de `[meta + [lang: lang_id], file]`, aquella línia llegís `[new_map, file]`.
    Això hauria de deixar més clar que simplement estem deixant el `file` al seu lloc original en la segona posició de la tupla. Simplement hem pres el valor `new_info` i l'hem incorporat al mapa que es troba en la primera posició.

    **I això ens porta de nou a l'estructura de canal `tuple val(meta), path(file)`!**

Un cop estigueu segurs que enteneu el que fa aquest codi, executeu el workflow per veure si ha funcionat:

```bash
nextflow run main.nf -resume
```

??? success "Sortida de la comanda"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [cheeky_fermat] DSL2 - revision: d096281ee4

    [4e/f722fe] IDENTIFY_LANGUAGE (7) [100%] 7 of 7, cached: 7 ✔
    [[id:sampleA, character:squirrel, lang:fr], /workspaces/training/side-quests/metadata/work/eb/f7148ebdd898fbe1136bec6a714acb/bonjour.txt]
    [[id:sampleB, character:tux, lang:de], /workspaces/training/side-quests/metadata/work/16/71d72410952c22cd0086d9bca03680/guten_tag.txt]
    [[id:sampleC, character:sheep, lang:de], /workspaces/training/side-quests/metadata/work/ea/04f5d979429e4455e14b9242fb3b45/hallo.txt]
    [[id:sampleD, character:turkey, lang:en], /workspaces/training/side-quests/metadata/work/c4/b7562adddc1cc0b7d414ec45d436eb/hello.txt]
    [[id:sampleF, character:moose, lang:fr], /workspaces/training/side-quests/metadata/work/5a/6c2b84bf8fadb98e28e216426be079/salut.txt]
    [[id:sampleE, character:stegosaurus, lang:es], /workspaces/training/side-quests/metadata/work/af/ee7c69bcab891c40d0529305f6b9e7/hola.txt]
    [[id:sampleG, character:turtle, lang:it], /workspaces/training/side-quests/metadata/work/4e/f722fe47271ba7ebcd69afa42964ca/ciao.txt]
    ```

Sí, això és correcte!
Hem reorganitzat ordenadament la sortida del procés de `meta, file, lang_id` de manera que `lang_id` és ara una de les claus del meta map, i les tuples del canal s'ajusten de nou al model `meta, file`.

<!-- TODO (futur) Hauríem de mostrar també com eliminar una clau amb subMap?! O indicar on trobar-ho. -->

### 2.4. Assignar un grup lingüístic mitjançant condicionals

Ara que tenim les nostres prediccions d'idioma, utilitzem la informació per assignar nous agrupaments.

En les nostres dades d'exemple, els idiomes que utilitzen els nostres personatges es poden agrupar en llengües germàniques (anglès, alemany) i llengües romàniques (francès, espanyol, italià).
Podria ser útil tenir aquesta classificació disponible en algun punt posterior del pipeline, de manera que afegim aquesta informació al meta map.

I, bona notícia, aquest és un altre cas que es presta perfectament a l'ús de l'operador `map`!

Concretament, definirem una variable anomenada `lang_group`, utilitzarem una lògica condicional senzilla per determinar quin valor assignar a `lang_group` per a cada peça de dades.

La sintaxi general tindrà aquest aspecte:

```groovy
.map { meta, file ->

    // la lògica condicional que defineix lang_group va aquí

    [meta + [lang_group: lang_group], file]
}
```

Podeu veure que és molt similar a l'operació de fusió de mapes al vol que vam utilitzar al pas anterior.
Només cal que escrivim les instruccions condicionals.

Aquí teniu la lògica condicional que volem aplicar:

- Definiu una variable anomenada `lang_group` amb el valor per defecte `'unknown'`.
- Si `lang` és alemany (`'de'`) o anglès (`'en'`), canvieu `lang_group` a `germanic`.
- En cas contrari, si `lang` es troba en una llista que conté francès (`'fr'`), espanyol (`'es'`) i italià (`'it'`), canvieu `lang_group` a `romance`.

Intenteu escriure-ho vosaltres mateixos si ja sabeu com escriure instruccions condicionals en Nextflow.

!!! tip "Consell"

    Podeu accedir al valor de `lang` dins de l'operació map amb `meta.lang`.

Hauríeu d'acabar fent els canvis següents al workflow:

=== "Després"

    ```groovy title="main.nf" linenums="13" hl_lines="7-19"
        // Executa langid per identificar l'idioma de cada salutació
        IDENTIFY_LANGUAGE(ch_datasheet)
        IDENTIFY_LANGUAGE.out
            .map { meta, file, lang_id ->
                [meta + [lang: lang_id], file]
            }
            .map { meta, file ->

                def lang_group = "unknown"
                if (meta.lang.equals("de") || meta.lang.equals('en')) {
                    lang_group = "germanic"
                }
                else if (meta.lang in ["fr", "es", "it"]) {
                    lang_group = "romance"
                }

                [meta + [lang_group: lang_group], file]
            }
            .view()
    ```

=== "Abans"

    ```groovy title="main.nf" linenums="13" hl_lines="7"
        // Executa langid per identificar l'idioma de cada salutació
        IDENTIFY_LANGUAGE(ch_datasheet)
        IDENTIFY_LANGUAGE.out
            .map { meta, file, lang_id ->
                [meta + [lang: lang_id], file]
            }
            .view()
    ```

Aquí teniu els punts clau:

- Utilitzem `def lang_group = "unknown"` per crear la variable `lang_group` amb el valor per defecte `unknown`.
- Utilitzem una estructura `if {} else if {}` per a la lògica condicional, amb proves `.equals()` alternatives per als dos idiomes germànics, i una prova d'existència en una llista per als tres idiomes romànics.
- Utilitzem l'operació de fusió `meta + [lang_group:lang_group]` com anteriorment per generar el meta map actualitzat.

<!-- TODO (futur) Afegir nota/enllaços a la documentació rellevant a la secció de recursos addicionals -->

Un cop tot tingui sentit, executeu el workflow de nou per veure el resultat:

```bash
nextflow run main.nf -resume
```

??? success "Sortida de la comanda"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [wise_almeida] DSL2 - revision: 46778c3cd0

    [da/652cc6] IDENTIFY_LANGUAGE (7) [100%] 7 of 7, cached: 7 ✔
    [[id:sampleA, character:squirrel, lang:fr, lang_group:romance], /workspaces/training/side-quests/metadata/data/bonjour.txt]
    [[id:sampleB, character:tux, lang:de, lang_group:germanic], /workspaces/training/side-quests/metadata/data/guten_tag.txt]
    [[id:sampleC, character:sheep, lang:de, lang_group:germanic], /workspaces/training/side-quests/metadata/data/hallo.txt]
    [[id:sampleD, character:turkey, lang:en, lang_group:germanic], /workspaces/training/side-quests/metadata/data/hello.txt]
    [[id:sampleE, character:stegosaurus, lang:es, lang_group:romance], /workspaces/training/side-quests/metadata/data/hola.txt]
    [[id:sampleF, character:moose, lang:fr, lang_group:romance], /workspaces/training/side-quests/metadata/data/salut.txt]
    [[id:sampleG, character:turtle, lang:it, lang_group:romance], /workspaces/training/side-quests/metadata/data/ciao.txt]
    ```

Com podeu veure, els elements del canal mantenen la seva estructura `[meta, file]`, però el meta map ara inclou aquesta nova classificació.

### Conclusió

En aquesta secció, heu après com:

- **Aplicar metadades d'entrada als canals de sortida**: Copiar metadades d'aquesta manera ens permet associar resultats posteriorment basant-nos en el contingut de les metadades.
- **Crear claus personalitzades**: Heu creat dues noves claus al vostre meta map, fusionant-les amb `meta + [nova_clau:valor]` al meta map existent. Una basada en un valor calculat d'un procés, i una altra basada en una condició que heu establert a l'operador `map`.

Aquestes us permeten associar metadades noves i existents amb fitxers a mesura que avanceu pel vostre pipeline.
Fins i tot si no esteu utilitzant metadades com a part d'un procés, mantenir el meta map associat amb les dades d'aquesta manera facilita tenir tota la informació rellevant junta.

---

## 3. Ús de la informació del meta map en un procés

Ara que sabeu com crear i actualitzar el meta map, podem arribar a la part realment divertida: utilitzar les metadades en un procés.

Més concretament, afegirem un segon pas al nostre workflow per dibuixar cada animal com a art ASCII i fer-lo dir el text gravat en un globus de diàleg.
Ho farem utilitzant una eina anomenada [`cowpy`](https://github.com/jeffbuttars/cowpy).

??? info "Què fa `cowpy`?"

    `cowpy` és una eina de línia de comandes que genera art ASCII per mostrar entrades de text arbitràries d'una manera divertida.
    És una implementació en Python de l'eina clàssica [cowsay](https://en.wikipedia.org/wiki/Cowsay) de Tony Monroe.

    ```console
    cowpy "Hello Nextflow"
    ```

    ```console
    ______________________________________________________
    < Hello Nextflow >
    ------------------------------------------------------
        \   ^__^
          \  (oo)\_______
            (__)\       )\/\
              ||----w |
              ||     ||
    ```

    Opcionalment, podeu seleccionar un personatge (o 'cowacter') per utilitzar en lloc de la vaca per defecte.

    ```console
    cowpy "Hello Nextflow" -c tux
    ```

    ```console
    __________________
    < Hello Nextflow >
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

Si heu treballat el curs Hello Nextflow, ja heu vist aquesta eina en acció.
Si no, no us preocupeu; cobrirem tot el que necessiteu saber a mesura que avancem.

### 3.1. Importar el procés i examinar el codi

Us proporcionem un mòdul de procés preescrit anomenat `COWPY` que encapsula l'eina `cowpy`, de manera que només cal afegir una instrucció include abans del bloc workflow.

Feu la següent edició al workflow:

=== "Després"

    ```groovy title="main.nf" linenums="1" hl_lines="4"
    #!/usr/bin/env nextflow

    include { IDENTIFY_LANGUAGE } from './modules/langid.nf'
    include { COWPY } from './modules/cowpy.nf'

    workflow {
    ```

=== "Abans"

    ```groovy title="main.nf" linenums="1"
    #!/usr/bin/env nextflow

    include { IDENTIFY_LANGUAGE } from './modules/langid.nf'

    workflow {
    ```

Podeu obrir el fitxer del mòdul per examinar el seu codi:

```groovy title="modules/cowpy.nf" linenums="1"
#!/usr/bin/env nextflow

// Genera art ASCII amb cowpy
process COWPY {

    publishDir "results/", mode: 'copy'

    container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'

    input:
    path input_file
    val character

    output:
    path "cowpy-${input_file}"

    script:
    """
    cat ${input_file} | cowpy -c ${character} > cowpy-${input_file}
    """
}
```

Com podeu veure, aquest procés està dissenyat actualment per prendre un fitxer d'entrada (que conté el text a mostrar) i un valor que especifica el personatge que s'ha de dibuixar en art ASCII, normalment proporcionat a nivell de workflow per un parametre de línia de comandes.

### 3.2. Passar un camp del meta map com a entrada

Quan vam utilitzar l'eina `cowpy` al curs Hello Nextflow, vam utilitzar un parametre de línia de comandes per determinar quin personatge utilitzar per dibuixar la imatge final.
Tenia sentit, perquè només generàvem una imatge per execució del pipeline.

Tanmateix, en aquest tutorial, volem generar una imatge adequada per a cada subjecte que estem processant, de manera que utilitzar un parametre de línia de comandes seria massa limitant.

Bona notícia: tenim una columna `character` al nostre full de dades i, per tant, al nostre meta map.
Utilitzem-la per establir el personatge que el procés hauria d'utilitzar per a cada entrada.

Per a això, haurem de fer tres coses:

1. Donar un nom al canal de sortida que prové del procés anterior per poder operar-hi de manera més còmoda.
2. Determinar com accedir a la informació d'interès
3. Afegir una crida al segon procés i alimentar-lo amb la informació adequada.

Comencem.

#### 3.2.1. Donar nom al canal de sortida anterior

Hem aplicat les manipulacions anteriors directament sobre el canal de sortida del primer procés, `IDENTIFY_LANGUAGE.out`.
Per tal d'alimentar el contingut d'aquest canal al procés següent (i fer-ho d'una manera clara i fàcil de llegir) volem donar-li el seu propi nom, `ch_languages`.

Podem fer-ho utilitzant l'operador [`set`](https://www.nextflow.io/docs/latest/reference/operator.html#set).

Al workflow principal, substituïu l'operador `.view()` per `.set { ch_languages }`, i afegiu una línia per comprovar que podem fer referència al canal pel seu nom.

=== "Després"

    ```groovy title="main.nf" linenums="14" hl_lines="19 21 22"
        // Executa langid per identificar l'idioma de cada salutació
        IDENTIFY_LANGUAGE(ch_datasheet)
        IDENTIFY_LANGUAGE.out
            .map { meta, file, lang_id ->
                [meta + [lang: lang_id], file]
            }
            .map { meta, file ->

                def lang_group = "unknown"
                if (meta.lang.equals("de") || meta.lang.equals('en')) {
                    lang_group = "germanic"
                }
                else if (meta.lang in ["fr", "es", "it"]) {
                    lang_group = "romance"
                }

                [meta + [lang_group: lang_group], file]
            }
            .set { ch_languages }

        // Temporal: examinar ch_languages
        ch_languages.view()
    ```

=== "Abans"

    ```groovy title="main.nf" linenums="14" hl_lines="19"
        // Executa langid per identificar l'idioma de cada salutació
        IDENTIFY_LANGUAGE(ch_datasheet)
        IDENTIFY_LANGUAGE.out
            .map { meta, file, lang_id ->
                [meta + [lang: lang_id], file]
            }
            .map { meta, file ->

                def lang_group = "unknown"
                if (meta.lang.equals("de") || meta.lang.equals('en')) {
                    lang_group = "germanic"
                }
                else if (meta.lang in ["fr", "es", "it"]) {
                    lang_group = "romance"
                }

                [meta + [lang_group: lang_group], file]
            }
            .view()
    ```

Executem-ho:

```bash
nextflow run main.nf
```

??? success "Sortida de la comanda"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `./main.nf` [friendly_austin] DSL2 - revision: 3dbe460fd6

    [36/cca6a7] IDENTIFY_LANGUAGE (7) | 7 of 7 ✔
    [[id:sampleB, character:tux, lang:de, lang_group:germanic], /workspaces/training/side-quests/metadata/work/e2/6db2402d83cf72081bcd2d11784714/guten_tag.txt]
    [[id:sampleA, character:squirrel, lang:fr, lang_group:romance], /workspaces/training/side-quests/metadata/work/6c/114c818317d169457d6e7336d5d55b/bonjour.txt]
    [[id:sampleC, character:sheep, lang:de, lang_group:germanic], /workspaces/training/side-quests/metadata/work/55/68c69c5efb527f3604ddb3daab8057/hallo.txt]
    [[id:sampleD, character:turkey, lang:en, lang_group:germanic], /workspaces/training/side-quests/metadata/work/2a/4752055ccb5d1370b0ef9da41d3993/hello.txt]
    [[id:sampleE, character:stegosaurus, lang:es, lang_group:romance], /workspaces/training/side-quests/metadata/work/f4/fcd3186dc666d5d239ffa6c37d125d/hola.txt]
    [[id:sampleF, character:moose, lang:fr, lang_group:romance], /workspaces/training/side-quests/metadata/work/c3/3b2627f733f278a7088332a5806108/salut.txt]
    [[id:sampleG, character:turtle, lang:it, lang_group:romance], /workspaces/training/side-quests/metadata/work/36/cca6a7dbfa26ac24f9329787a32e9d/ciao.txt]
    ```

Això confirma que ara podem fer referència al canal pel seu nom.

#### 3.2.2. Accedir al fitxer i a les metadades del personatge

Sabem, mirant el codi del mòdul, que el procés `COWPY` espera rebre un fitxer de text i un valor `character`.
Per escriure la crida al procés `COWPY`, només cal saber com extreure l'objecte fitxer i les metadades corresponents de cada element del canal.

Com sol ser el cas, la manera més senzilla de fer-ho és utilitzar una operació `map`.

El nostre canal conté tuples estructurades com `[meta, file]`, de manera que podem accedir directament a l'objecte `file`, i podem accedir al valor `character` emmagatzemat dins del meta map fent referència a ell com `meta.character`.

Al workflow principal, feu els canvis de codi següents:

=== "Després"

    ```groovy title="main.nf" linenums="34"
        // Temporal: accedir al fitxer i al personatge
        ch_languages.map { meta, file -> file }.view { file -> "File: " + file }
        ch_languages.map { meta, file -> meta.character }.view { character -> "Character: " + character }
    ```

=== "Abans"

    ```groovy title="main.nf" linenums="34"
        // Temporal: examinar ch_languages
        ch_languages.view()
    ```

Tingueu en compte que estem utilitzant closures (com ara `{ file -> "File: " + file }`) per fer que la sortida de les operacions `.view` sigui més llegible.

Executem-ho:

```bash
nextflow run main.nf -resume
```

??? success "Sortida de la comanda"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `./main.nf` [cheesy_cantor] DSL2 - revision: 15af9c1ec7

    [43/05df08] IDENTIFY_LANGUAGE (7) [100%] 7 of 7, cached: 7 ✔
    Character: squirrel
    File: /workspaces/training/side-quests/metadata/work/8d/4b9498bbccb7a74f04e41877cdc3e5/bonjour.txt
    File: /workspaces/training/side-quests/metadata/work/d3/604274985406e40d79021dea658e60/guten_tag.txt
    Character: tux
    Character: turkey
    File: /workspaces/training/side-quests/metadata/work/d4/fafcc9415b61d2b0fea872e6a05e8a/hello.txt
    File: /workspaces/training/side-quests/metadata/work/02/468ac9efb27f636715e8144b37e9a7/hallo.txt
    Character: sheep
    Character: moose
    Character: stegosaurus
    File: /workspaces/training/side-quests/metadata/work/d4/61a7e1188b4f2742bc72004e226eca/salut.txt
    File: /workspaces/training/side-quests/metadata/work/ae/68364be238c11149c588bf6fc858b1/hola.txt
    File: /workspaces/training/side-quests/metadata/work/43/05df081af5d879ab52e5828fa0357e/ciao.txt
    Character: turtle
    ```

_Les rutes dels fitxers i els valors dels personatges poden aparèixer en un ordre diferent a la vostra sortida._

Això confirma que podem accedir al fitxer i al personatge per a cada element del canal.

#### 3.2.3. Cridar el procés `COWPY`

Ara posem-ho tot junt i cridem realment el procés `COWPY` sobre el canal `ch_languages`.

Al workflow principal, feu els canvis de codi següents:

=== "Després"

    ```groovy title="main.nf" linenums="34"
        // Executa cowpy per generar art ASCII
        COWPY(
            ch_languages.map { meta, file -> file },
            ch_languages.map { meta, file -> meta.character }
        )
    ```

=== "Abans"

    ```groovy title="main.nf" linenums="34"
        // Temporal: accedir al fitxer i al personatge
        ch_languages.map { meta, file -> [file, meta.character] }
            .view()
    ```

Veieu que simplement copiem les dues operacions map (sense les instruccions `.view()`) com a entrades a la crida del procés.
Assegureu-vos de no oblidar la coma entre elles!

És una mica incòmode, però veurem com millorar-ho a la secció següent.

Executem-ho:

```bash
nextflow run main.nf -resume
```

??? success "Sortida de la comanda"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [suspicious_crick] DSL2 - revision: 25541014c5

    executor >  local (7)
    [43/05df08] IDENTIFY_LANGUAGE (7) [100%] 7 of 7, cached: 7 ✔
    [e7/317c18] COWPY (6)             [100%] 7 of 7 ✔
    ```

Si mireu al directori de resultats, hauríeu de veure els fitxers individuals que contenen l'art ASCII de cada salutació pronunciada pel personatge corresponent.

??? abstract "Contingut del directori i exemple de fitxer"

    ```console
    results/
    ├── cowpy-bonjour.txt
    ├── cowpy-ciao.txt
    ├── cowpy-guten_tag.txt
    ├── cowpy-hallo.txt
    ├── cowpy-hello.txt
    ├── cowpy-hola.txt
    └── cowpy-salut.txt
    ```

    ```text title="results/cowpy-bonjour.txt"
     _________________
    / Bonjour         \
    \ Salut, à demain /
    -----------------
      \
        \
                      _ _
          | \__/|  .~    ~.
          /oo `./      .'
          {o__,   \    {
            / .  . )    \
            `-` '-' \    }
          .(   _(   )_.'
          '---.~_ _ _|
    ```

Això demostra que vam poder utilitzar la informació del meta map per parametritzar la comanda al segon pas del pipeline.

Tanmateix, com s'ha assenyalat anteriorment, part del codi implicat era una mica incòmode, ja que vam haver de desempaquetar les metadades mentre encara estàvem en el context del cos del workflow.
Aquest enfocament funciona bé per utilitzar un nombre petit de camps del meta map, però escalaria malament si volguéssim utilitzar-ne molts més.

Hi ha un altre operador anomenat `multiMap()` que ens permet agilitzar-ho una mica, però fins i tot llavors no és ideal.

??? info "(Opcional) Versió alternativa amb `multiMap()`"

    Per si us ho pregunteu, no podíem simplement escriure una única operació `map()` que produís tant el `file` com el `character`, perquè això els retornaria com una tupla.
    Vam haver d'escriure dues operacions `map()` separades per alimentar els elements `file` i `character` al procés per separat.

    Tècnicament hi ha una altra manera de fer-ho mitjançant una única operació de mapeig, utilitzant l'operador `multiMap()`, que és capaç d'emetre múltiples canals.
    Per exemple, podríeu substituir la crida a `COWPY` anterior pel codi següent:

    === "Després"

        ```groovy title="main.nf" linenums="34"
            // Executa cowpy per generar art ASCII
            COWPY(
                ch_languages.multiMap { meta, file ->
                    file: file
                    character: meta.character
                }
            )
        ```

    === "Abans"

        ```groovy title="main.nf" linenums="34"
            // Executa cowpy per generar art ASCII
            COWPY(
                ch_languages.map { meta, file -> file },
                ch_languages.map { meta, file -> meta.character }
            )
        ```

    Això produeix exactament el mateix resultat.

En qualsevol cas, és incòmode haver de fer algun desempaquetament a nivell de workflow.

Seria millor si poguéssim alimentar el meta map sencer al procés i seleccionar el que necessitem un cop allà.

### 3.3. Passar i utilitzar el meta map sencer

El propòsit del meta map és, al cap i a la fi, passar totes les metadades juntes com un paquet.
L'única raó per la qual no podíem fer-ho anteriorment és que el procés no estava configurat per acceptar un meta map.
Però com que controlem el codi del procés, podem canviar-ho.

Modifiquem el procés `COWPY` per acceptar l'estructura de tupla `[meta, file]` que vam utilitzar al primer procés per poder agilitzar el workflow.

Per a això, haurem de fer tres coses:

1. Modificar les definicions d'entrada del mòdul del procés `COWPY`
2. Actualitzar la comanda del procés per utilitzar el meta map
3. Actualitzar la crida al procés al cos del workflow

Preparats? Comencem!

#### 3.3.1. Modificar l'entrada del mòdul `COWPY`

Feu les edicions següents al fitxer del mòdul `cowpy.nf`:

=== "Després"

    ```groovy title="cowpy.nf" linenums="10" hl_lines="2"
    input:
    tuple val(meta), path(input_file)
    ```

=== "Abans"

    ```groovy title="cowpy.nf" linenums="10" hl_lines="2-3"
    input:
    path(input_file)
    val character
    ```

Això ens permet utilitzar l'estructura de tupla `[meta, file]` que hem vist anteriorment al tutorial.

Tingueu en compte que no hem actualitzat la definició de sortida del procés per produir el meta map, per tal de mantenir el tutorial simplificat, però podeu fer-ho vosaltres mateixos com a exercici seguint el model del procés `IDENTIFY_LANGUAGE`.

#### 3.3.2. Actualitzar la comanda per utilitzar el camp del meta map

El meta map sencer ara està disponible dins del procés, de manera que podem fer referència a la informació que conté directament des de dins del bloc de comandes.

Feu les edicions següents al fitxer del mòdul `cowpy.nf`:

=== "Després"

    ```groovy title="cowpy.nf" linenums="16" hl_lines="3"
    script:
    """
    cat ${input_file} | cowpy -c ${meta.character} > cowpy-${input_file}
    """
    ```

=== "Abans"

    ```groovy title="cowpy.nf" linenums="16" hl_lines="3"
    script:
    """
    cat ${input_file} | cowpy -c ${character} > cowpy-${input_file}
    """
    ```

Hem substituït la referència al valor `character` que anteriorment es passava com a entrada independent pel valor contingut al meta map, al qual fem referència utilitzant `meta.character`.

Ara actualitzem la crida al procés en conseqüència.

#### 3.3.3. Actualitzar la crida al procés i executar-lo

El procés ara espera que la seva entrada utilitzi l'estructura de tupla `[meta, file]`, que és el que produeix el procés anterior, de manera que simplement podem passar el canal `ch_languages` sencer al procés `COWPY`.

Feu les edicions següents al workflow principal:

=== "Després"

    ```groovy title="main.nf" linenums="34" hl_lines="2"
    // Executa cowpy per generar art ASCII
    COWPY(ch_languages)
    ```

=== "Abans"

    ```groovy title="main.nf" linenums="34" hl_lines="3-4"
    // Executa cowpy per generar art ASCII
    COWPY(
        ch_languages.map { meta, file -> file },
        ch_languages.map { meta, file -> meta.character }
    )
    ```

Això simplifica considerablement la crida!

Esborrem els resultats de l'execució anterior i executem-ho:

```bash
rm -r results
nextflow run main.nf
```

??? success "Sortida de la comanda"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [wise_sammet] DSL2 - revision: 99797b1e92

    executor >  local (14)
    [5d/dffd4e] process > IDENTIFY_LANGUAGE (7) [100%] 7 of 7 ✔
    [25/9243df] process > COWPY (7)             [100%] 7 of 7 ✔
    ```

Si mireu al directori de resultats, hauríeu de veure les mateixes sortides que abans, _és a dir_, fitxers individuals que contenen l'art ASCII de cada salutació pronunciada pel personatge corresponent.

??? abstract "Contingut del directori"

    ```console
    ./results/
    ├── cowpy-bonjour.txt
    ├── cowpy-ciao.txt
    ├── cowpy-guten_tag.txt
    ├── cowpy-hallo.txt
    ├── cowpy-hello.txt
    ├── cowpy-hola.txt
    └── cowpy-salut.txt
    ```

Per tant, això produeix els mateixos resultats que abans amb un codi més senzill.

Per descomptat, això assumeix que podeu modificar el codi del procés.
En alguns casos, potser haureu de dependre de processos existents que no teniu llibertat per modificar, la qual cosa limita les vostres opcions.
La bona notícia, si teniu previst utilitzar mòduls del projecte [nf-core](https://nf-co.re/), és que els mòduls nf-core estan tots configurats per utilitzar l'estructura de tupla `[meta, file]` com a estàndard.

### 3.4. Resolució de problemes amb entrades requerides que falten

El valor `character` és necessari perquè el procés `COWPY` s'executi correctament.
Si no establim un valor per defecte en un fitxer de configuració, HEM de proporcionar un valor per a ell al full de dades.

**Què passa si no ho fem?**
Depèn del que contingui el full de dades d'entrada i de quina versió del workflow estem executant.

#### 3.4.1. La columna character existeix però està buida

Suposem que eliminem el valor del personatge per a una de les entrades del nostre full de dades per simular un error de recollida de dades:

```csv title="datasheet.csv" linenums="1" hl_lines="2"
id,character,recording
sampleA,,/workspaces/training/side-quests/metadata/data/bonjour.txt
sampleB,tux,/workspaces/training/side-quests/metadata/data/guten_tag.txt
sampleC,sheep,/workspaces/training/side-quests/metadata/data/hallo.txt
sampleD,turkey,/workspaces/training/side-quests/metadata/data/hello.txt
sampleE,stegosaurus,/workspaces/training/side-quests/metadata/data/hola.txt
sampleF,moose,/workspaces/training/side-quests/metadata/data/salut.txt
sampleG,turtle,/workspaces/training/side-quests/metadata/data/ciao.txt
```

Per a qualsevol de les versions del workflow que hem utilitzat anteriorment, la clau `character` es crearà per a totes les entrades quan es llegeixi el full de dades, però per a `sampleA` el valor serà una cadena buida.

Això causarà un error.

??? failure "Sortida de la comanda"

    ```console hl_lines="8 11 16 28"
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [marvelous_hirsch] DSL2 - revision: 0dfeee3cc1

    executor >  local (9)
    [c1/c5dd4f] process > IDENTIFY_LANGUAGE (7) [ 85%] 6 of 7
    [d3/b7c415] process > COWPY (2)             [  0%] 0 of 6
    ERROR ~ Error executing process > 'COWPY (1)'

    Caused by:
      Process `COWPY (1)` terminated with an error exit status (2)


    Command executed:

      cat bonjour.txt | cowpy -c  > cowpy-bonjour.txt

    Command exit status:
      2

    Command output:
      (empty)

    Command error:
      usage: cowpy [-h] [-l] [-L] [-t] [-u] [-e EYES] [-c COWACTER] [-E] [-r] [-x]
                  [-C]
                  [msg ...]
      cowpy: error: argument -c/--cowacter: expected one argument

    Work dir:
      /workspaces/training/side-quests/metadata/work/ca/9d49796612a54dec5ed466063c809b

    Container:
      community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273

    Tip: you can try to figure out what's wrong by changing to the process work dir and showing the script file named `.command.sh`

    -- Check '.nextflow.log' file for details
    ```

Quan Nextflow executa la línia de comandes `cowpy` per a aquesta mostra, `${meta.character}` s'omple amb una cadena buida a la línia de comandes `cowpy`, de manera que l'eina `cowpy` llança un error dient que no s'ha proporcionat cap valor per a l'argument `-c`.

#### 3.4.2. La columna character no existeix al full de dades

Ara suposem que eliminem la columna `character` completament del nostre full de dades:

```csv title="datasheet.csv" linenums="1"
id,recording
sampleA,/workspaces/training/side-quests/metadata/data/bonjour.txt
sampleB,/workspaces/training/side-quests/metadata/data/guten_tag.txt
sampleC,/workspaces/training/side-quests/metadata/data/hallo.txt
sampleD,/workspaces/training/side-quests/metadata/data/hello.txt
sampleE,/workspaces/training/side-quests/metadata/data/hola.txt
sampleF,/workspaces/training/side-quests/metadata/data/salut.txt
sampleG,/workspaces/training/side-quests/metadata/data/ciao.txt
```

En aquest cas, la clau `character` no es crearà en absolut quan es llegeixi el full de dades.

##### 3.4.2.1. Valor accedit a nivell de workflow

Si estem utilitzant la versió del codi que vam escriure a la secció 3.2, Nextflow intentarà accedir a la clau `character` al meta map ABANS de cridar el procés `COWPY`.

No trobarà cap element que coincideixi amb la instrucció, de manera que no executarà `COWPY` en absolut.

??? success "Sortida de la comanda"

    ```console hl_lines="7"
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [desperate_montalcini] DSL2 - revision: 0dfeee3cc1

    executor >  local (7)
    [1a/df2544] process > IDENTIFY_LANGUAGE (7) [100%] 7 of 7 ✔
    [-        ] process > COWPY                 -
    ```

Pel que fa a Nextflow, aquest workflow s'ha executat correctament!
Tanmateix, no es produirà cap de les sortides que volem.

##### 3.4.2.2. Valor accedit a nivell de procés

Si estem utilitzant la versió de la secció 3.3, Nextflow passarà el meta map sencer al procés `COWPY` i intentarà executar la comanda.

Això causarà un error, però diferent del primer cas.

??? failure "Sortida de la comanda"

    ```console hl_lines="8 11 16"
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [jovial_bohr] DSL2 - revision: eaaf375827

    executor >  local (9)
    [0d/ada9db] process > IDENTIFY_LANGUAGE (5) [ 85%] 6 of 7
    [06/28065f] process > COWPY (2)             [  0%] 0 of 6
    ERROR ~ Error executing process > 'COWPY (2)'

    Caused by:
      Process `COWPY (2)` terminated with an error exit status (1)


    Command executed:

      cat guten_tag.txt | cowpy -c null > cowpy-guten_tag.txt

    Command exit status:
      1

    Command output:
      (empty)

    Command error:
      Traceback (most recent call last):
        File "/opt/conda/bin/cowpy", line 10, in <module>
          sys.exit(main())
                  ~~~~^^
        File "/opt/conda/lib/python3.13/site-packages/cowpy/cow.py", line 1215, in main
          print(cow(eyes=args.eyes,
                ~~~^^^^^^^^^^^^^^^^
                tongue=args.tongue,
                ^^^^^^^^^^^^^^^^^^^
                thoughts=args.thoughts
                ^^^^^^^^^^^^^^^^^^^^^^
                    ).milk(msg)
                    ^
      TypeError: 'str' object is not callable

    Work dir:
      /workspaces/training/side-quests/metadata/work/06/28065f7d9fd7d22bba084aa941b6d6

    Container:
      community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273

    Tip: you can replicate the issue by changing to the process work dir and entering the command `bash .command.run`

    -- Check '.nextflow.log' file for details
    ```

Això passa perquè `meta.character` no existeix, de manera que el nostre intent d'accedir-hi retorna `null`. Com a resultat, Nextflow literalment insereix `null` a la línia de comandes, que per descomptat no és reconegut per l'eina `cowpy`.

#### 3.4.3. Solucions

A part de proporcionar un valor per defecte com a part de la configuració del workflow, hi ha dues coses que podem fer per gestionar-ho de manera més robusta:

1. Implementar validació d'entrada al vostre workflow per assegurar-vos que el full de dades conté tota la informació requerida. Podeu trobar una [introducció a la validació d'entrada](../hello_nf-core/05_input_validation.md) al curs de formació Hello nf-core. <!-- TODO (futur) pendent d'una missió secundària de Validació pròpia -->

2. Si voleu assegurar-vos que qualsevol persona que utilitzi el vostre mòdul de procés pugui identificar immediatament les entrades requerides, també podeu fer que la propietat de metadades requerida sigui una entrada explícita.

Aquí teniu un exemple de com funcionaria.

Primer, a nivell de procés, actualitzeu la definició d'entrada de la manera següent:

=== "Després"

    ```groovy title="cowpy.nf" linenums="12" hl_lines="2"
        input:
        tuple val(meta), val(character), path(input_file)
    ```

=== "Abans"

    ```groovy title="cowpy.nf" linenums="12" hl_lines="2"
        input:
        tuple val(meta), path(input_file)
    ```

Després, a nivell de workflow, utilitzeu una operació de mapeig per extreure la propietat `character` de les metadades i convertir-la en un component explícit de la tupla d'entrada:

=== "Després"

    ```groovy title="main.nf" linenums="37" hl_lines="1"
        COWPY(ch_languages.map{meta, file -> [meta, meta.character, file]})
    ```

=== "Abans"

    ```groovy title="main.nf" linenums="37" hl_lines="1"
        COWPY(ch_languages)
    ```

Aquest enfocament té l'avantatge de mostrar explícitament que `character` és necessari, i fa que el procés sigui més fàcil de reutilitzar en altres contextos.

Això posa de manifest un principi de disseny important:

**Utilitzeu el meta map per a informació opcional i descriptiva, però extraieu els valors requerits com a entrades explícites.**

El meta map és excel·lent per mantenir les estructures de canal netes i evitar estructures de canal arbitràries, però per als elements obligatoris que es referencien directament en un procés, extreure'ls com a entrades explícites crea un codi més robust i mantenible.

### Conclusió

En aquesta secció, heu après com utilitzar metadades per personalitzar l'execució d'un procés, accedint-hi ja sigui a nivell de workflow o a nivell de procés.

---

## Exercici suplementari

Si voleu practicar l'ús de la informació del meta map des de dins d'un procés, proveu d'utilitzar altres peces d'informació del meta map com ara `lang` i `lang_group` per personalitzar com s'anomenen i/o s'organitzen les sortides.

Per exemple, proveu de modificar el codi per produir aquest resultat:

```console title="Results directory contents"
results/
├── germanic
│   ├── de-guten_tag.txt
│   ├── de-hallo.txt
│   └── en-hello.txt
└── romance
    ├── es-hola.txt
    ├── fr-bonjour.txt
    ├── fr-salut.txt
    └── it-ciao.txt
```

<!-- TODO (futur) Proporcionar una solució elaborada -->
<!-- el canvi de nom hauria d'utilitzar el meta dins del procés -->
<!-- l'organització de la sortida hauria d'utilitzar el meta a les sortides del workflow -->

---

## Resum

En aquesta missió secundària, heu explorat com treballar eficaçment amb metadades en workflows de Nextflow.

Aquest patró de mantenir les metadades explícites i adjuntes a les dades és una pràctica recomanada fonamental en Nextflow, que ofereix diversos avantatges respecte a la codificació de la informació dels fitxers:

- Les metadades dels fitxers romanen associades amb els fitxers al llarg del workflow
- El comportament dels processos es pot personalitzar per fitxer
- L'organització de les sortides pot reflectir les metadades dels fitxers
- La informació dels fitxers es pot ampliar durant l'execució del pipeline

Aplicar aquest patró en el vostre propi treball us permetrà construir workflows de bioinformàtica robustos i mantenibles.

### Patrons clau

1.  **Lectura i estructuració de metadades:** Llegir fitxers CSV i crear mapes de metadades organitzats que romanguin associats amb els vostres fitxers de dades.

    ```groovy
    channel.fromPath('datasheet.csv')
      .splitCsv(header: true)
      .map { row ->
          [ [id:row.id, character:row.character], row.recording ]
      }
    ```

2.  **Ampliació de metadades durant el workflow:** Afegir nova informació a les vostres metadades a mesura que el vostre pipeline avança, afegint sortides de processos i derivant valors mitjançant lògica condicional.

    - Afegir noves claus basades en la sortida del procés

    ```groovy
    .map { meta, file, lang ->
      [ meta + [lang:lang], file ]
    }
    ```

    - Afegir noves claus utilitzant una clàusula condicional

    ```groovy
    .map{ meta, file ->
        if ( meta.lang.equals("de") || meta.lang.equals('en') ){
            lang_group = "germanic"
        } else if ( meta.lang in ["fr", "es", "it"] ) {
            lang_group = "romance"
        } else {
            lang_group = "unknown"
        }
    }
    ```

3.  **Personalització del comportament del procés:** Utilitzar metadades dins del procés.

    ```groovy
    cat $input_file | cowpy -c ${meta.character} > cowpy-${input_file}
    ```

### Recursos addicionals

- [map](https://www.nextflow.io/docs/latest/operator.html#map)
- [stdout](https://www.nextflow.io/docs/latest/process.html#outputs)

---

## Què segueix?

Torneu al [menú de missions secundàries](../) o feu clic al botó a la part inferior dreta de la pàgina per continuar amb el tema següent de la llista.
