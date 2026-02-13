# Metadades i mapes de metadades

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducció assistida per IA - [més informació i suggeriments](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

En qualsevol anàlisi científica, rarament treballem només amb els fitxers de dades en brut.
Cada fitxer ve amb la seva pròpia informació addicional: què és, d'on prové i què el fa especial.
Aquesta informació extra és el que anomenem metadades.

Les metadades són dades que descriuen altres dades.
Les metadades fan un seguiment de detalls importants sobre els fitxers i les condicions experimentals, i ajuden a adaptar les anàlisis a les característiques úniques de cada conjunt de dades.

Penseu-hi com un catàleg de biblioteca: mentre que els llibres contenen el contingut real (dades en brut), les fitxes del catàleg proporcionen informació essencial sobre cada llibre—quan es va publicar, qui el va escriure, on trobar-lo (metadades).
En pipelines de Nextflow, les metadades es poden utilitzar per:

- Fer un seguiment de la informació específica dels fitxers al llarg del workflow
- Configurar processos basant-se en les característiques dels fitxers
- Agrupar fitxers relacionats per a una anàlisi conjunta

### Objectius d'aprenentatge

En aquesta missió secundària, explorarem com gestionar metadades en workflows.
Començant amb un full de dades senzill (sovint anomenat samplesheet en bioinformàtica) que conté informació bàsica dels fitxers, aprendreu com:

- Llegir i analitzar metadades de fitxers des de fitxers CSV
- Crear i manipular mapes de metadades
- Afegir nous camps de metadades durant l'execució del workflow
- Utilitzar metadades per personalitzar el comportament dels processos

Aquestes habilitats us ajudaran a construir pipelines més robustos i flexibles que poden gestionar relacions complexes entre fitxers i requisits de processament.

### Prerequisits

Abans d'emprendre aquesta missió secundària, hauríeu de:

- Haver completat el tutorial [Hello Nextflow](../hello_nextflow/README.md) o un curs equivalent per a principiants.
- Sentir-vos còmodes utilitzant conceptes i mecanismes bàsics de Nextflow (processos, canals, operadors)

---

## 0. Primers passos

#### Obriu l'espai de codi de formació

Si encara no ho heu fet, assegureu-vos d'obrir l'entorn de formació tal com es descriu a [Configuració de l'entorn](../envsetup/index.md).

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### Moveu-vos al directori del projecte

Anem al directori on es troben els fitxers per a aquest tutorial.

```bash
cd side-quests/metadata
```

Podeu configurar VSCode perquè se centri en aquest directori:

```bash
code .
```

#### Reviseu els materials

Trobareu un fitxer de workflow principal i un directori `data` que conté un full de dades i un grapat de fitxers de dades.

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

El workflow al fitxer `main.nf` és un esbós que expandireu gradualment fins a convertir-lo en un workflow completament funcional.

El full de dades llista els camins als fitxers de dades i algunes metadades associades, organitzades en 3 columnes:

- `id`: autoexplicatiu, un ID donat al fitxer
- `character`: un nom de personatge, que utilitzarem més endavant per dibuixar diferents criatures
- `data`: camins a fitxers `.txt` que contenen salutacions en diferents idiomes

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

Cada fitxer de dades conté text de salutació en un de cinc idiomes (fr: francès, de: alemany, es: espanyol, it: italià, en: anglès).

També us proporcionarem una eina d'anàlisi de llenguatge contenidoritzada anomenada `langid`.

#### Reviseu l'assignació

El vostre repte és escriure un workflow de Nextflow que:

1. **Identifiqui** l'idioma de cada fitxer automàticament
2. **Agrupí** els fitxers per família lingüística (llengües germàniques vs llengües romàniques)
3. **Personalitzi** el processament de cada fitxer basant-se en el seu idioma i metadades
4. **Organitzi** les sortides per grup lingüístic

Això representa un patró de workflow típic on les metadades específiques dels fitxers impulsen les decisions de processament; exactament el tipus de problema que els mapes de metadades resolen elegantment.

#### Llista de verificació de preparació

Creieu que esteu preparats per començar?

- [ ] Entenc l'objectiu d'aquest curs i els seus prerequisits
- [ ] El meu espai de codi està en funcionament
- [ ] He establert el meu directori de treball adequadament
- [ ] Entenc l'assignació

Si podeu marcar totes les caselles, esteu preparats per començar.

---

## 1. Carregar metadades des d'un full de dades

Obriu el fitxer de workflow `main.nf` per examinar l'esbós de workflow que us donem com a punt de partida.

```groovy title="main.nf" linenums="1"
#!/usr/bin/env nextflow

workflow  {

    ch_datasheet = channel.fromPath("./data/datasheet.csv")

}
```

Podeu veure que hem configurat una factoria de canals bàsica per carregar el full de dades d'exemple com a fitxer, però això encara no llegirà el contingut del fitxer.
Comencem afegint això.

### 1.1. Llegir el contingut amb `splitCsv`

Necessitem triar un operador que analitzi el contingut del fitxer adequadament amb un mínim esforç per part nostra.
Com que el nostre full de dades està en format CSV, aquesta és una feina per a l'operador [`splitCsv`](https://www.nextflow.io/docs/latest/reference/operator.html#splitcsv), que carrega cada fila del fitxer com un element al canal.

Feu els següents canvis per afegir una operació `splitCsv()` al codi de construcció del canal, més una operació `view()` per comprovar que el contingut del fitxer s'està carregant correctament al canal.

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

Vegem què surt d'això, d'acord?
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

Això és genial! Facilita l'accés a camps específics de cada fitxer.
Per exemple, podríem accedir a l'ID del fitxer amb `id` o al camí del fitxer txt amb `recording`.

??? info "(Opcional) Més sobre mapes"

    A Groovy, el llenguatge de programació sobre el qual es construeix Nextflow, un mapa és una estructura de dades clau-valor similar als diccionaris en Python, objectes en JavaScript o hashes en Ruby.

    Aquí teniu un script executable que mostra com podeu definir un mapa i accedir al seu contingut a la pràctica:

    ```groovy title="examples/map_demo.nf"
    #!/usr/bin/env nextflow

    // Crea un mapa senzill
    def my_map = [id:'sampleA', character:'squirrel']

    // Imprimeix tot el mapa
    println "map: ${my_map}"

    // Accedeix a valors individuals utilitzant notació de punt
    println "id: ${my_map.id}"
    println "character: ${my_map.character}"
    ```

    Tot i que no té un bloc `workflow` adequat, Nextflow pot executar-lo com si fos un workflow:

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

Diguem que volem accedir a la columna `character` del full de dades i imprimir-la.
Podem utilitzar l'operador `map` de Nextflow per iterar sobre cada element del nostre canal i seleccionar específicament l'entrada `character` de l'objecte mapa.

Feu les següents edicions al workflow:

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

### 1.3. Organitzar les metadades en un 'mapa de metadades'

En l'estat actual del workflow, els fitxers d'entrada (sota la clau `recording`) i les metadades associades (`id`, `character`) estan tots al mateix nivell, com si estiguessin tots en una gran bossa.
La conseqüència pràctica és que cada procés que consumeixi aquest canal hauria de ser configurat amb aquesta estructura en ment:

```groovy
    input:
    tuple val(id), val(character), file(recording)
```

Això està bé mentre el nombre de columnes al full de dades no canviï.
No obstant això, si afegiu fins i tot només una columna al full de dades, la forma del canal ja no coincidirà amb el que el procés espera, i el workflow produirà errors.
També fa que el procés sigui difícil de compartir amb altres que podrien tenir dades d'entrada lleugerament diferents, i podríeu acabar havent de codificar variables al procés que no són necessàries pel bloc de script.

Per evitar aquest problema, necessitem trobar una manera de mantenir l'estructura del canal consistent independentment de quantes columnes contingui el full de dades.

Podem fer-ho recollint totes les metadades en un element dins de la tupla, que anomenarem el mapa de metadades, o més simplement 'mapa meta'.

Feu les següents edicions a l'operació `map`:

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

Hem reestructurat els elements del nostre canal en una tupla que consisteix en dos elements, el mapa meta i l'objecte fitxer corresponent.

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

Com a resultat, afegir més columnes al full de dades farà que més metadades estiguin disponibles al mapa `meta`, però no canviarà la forma del canal.
Això ens permet escriure processos que consumeixen el canal sense haver de codificar els elements de metadades a l'especificació d'entrada:

```groovy title="Syntax example"
    input:
    tuple val(meta), file(recording)
```

Aquesta és una convenció àmpliament utilitzada per organitzar metadades en workflows de Nextflow.

### Conclusió

En aquesta secció, heu après:

- **Per què les metadades són importants:** Mantenir les metadades amb les vostres dades preserva informació important dels fitxers al llarg del workflow.
- **Com llegir fulls de dades:** Utilitzar `splitCsv` per llegir fitxers CSV amb informació de capçalera i transformar files en dades estructurades
- **Com crear un mapa meta:** Separar metadades de dades de fitxers utilitzant l'estructura de tupla `[ [id:value, ...], file ]`

---

## 2. Manipular metadades

Ara que tenim les nostres metadades carregades, fem-ne alguna cosa!

Utilitzarem una eina anomenada [`langid`](https://github.com/saffsd/langid.py) per identificar l'idioma contingut a cada fitxer de gravació de criatura.
L'eina ve pre-entrenada en un conjunt d'idiomes, i donat un fragment de text, produirà una predicció d'idioma i una puntuació de probabilitat associada, ambdues a `stdout`.

### 2.1. Importar el procés i examinar el codi

Us proporcionem un mòdul de procés pre-escrit anomenat `IDENTIFY_LANGUAGE` que encapsula l'eina `langid`, així que només necessiteu afegir una declaració d'inclusió abans del bloc de workflow.

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
Aquest patró `tuple val(meta), path(file), <output>` manté les metadades associades tant amb les dades d'entrada com amb les sortides mentre flueix pel pipeline.

Tingueu en compte que estem utilitzant el qualificador de sortida [`stdout`](https://www.nextflow.io/docs/latest/process.html#outputs) de Nextflow aquí perquè l'eina imprimeix la seva sortida directament a la consola en lloc d'escriure un fitxer; i utilitzem `sed` a la línia de comandes per eliminar la puntuació de probabilitat, netejar la cadena eliminant caràcters de nova línia, i retornar només la predicció d'idioma.

### 2.2. Afegir una crida a `IDENTIFY_LANGUAGE`

Ara que el procés està disponible per al workflow, podem afegir una crida al procés `IDENTIFY_LANGUAGE` per executar-lo sobre el canal de dades.

Feu les següents edicions al workflow:

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

Tingueu en compte que hem eliminat l'operació `.view()` original a la construcció del canal.

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

I com s'ha indicat anteriorment, també hem inclòs el fitxer d'entrada i el mapa meta a la sortida, el que significa que ambdós romanen associats amb la nova informació que acabem de produir.
Això serà útil al següent pas.

!!! note "Nota"

    Més generalment, aquest patró de mantenir el mapa meta associat amb els resultats facilita l'associació de resultats relacionats que comparteixen els mateixos identificadors.

    Com ja haureu après, no podeu confiar en l'ordre dels elements als canals per fer coincidir resultats entre ells.
    En canvi, heu d'utilitzar claus per associar dades correctament, i els mapes meta proporcionen una estructura ideal per a aquest propòsit.

    Explorem aquest cas d'ús en detall a la missió secundària [Splitting & Grouping](./splitting_and_grouping.md).

### 2.3. Augmentar metadades amb sortides de processos

Atès que els resultats que acabem de produir són en si mateixos una forma de metadades sobre el contingut dels fitxers, seria útil afegir-los al nostre mapa meta.

No obstant això, no volem modificar el mapa meta existent in situ.
Des d'un punt de vista tècnic, és _possible_ fer-ho, però no és segur.

Així que en canvi, crearem un nou mapa meta que contingui el contingut del mapa meta existent més un nou parell clau-valor `lang: lang_id` que conté la nova informació, utilitzant l'operador `+` (una característica de Groovy).
I combinarem això amb una operació [`map`](https://www.nextflow.io/docs/latest/operator.html#map) per reemplaçar el mapa antic amb el nou.

Aquí teniu les edicions que necessiteu fer al workflow:

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

Si encara no esteu familiaritzats amb l'operador `+`, o si això sembla confús, preneu-vos uns minuts per repassar l'explicació detallada a continuació.

??? info "Creació del nou mapa meta utilitzant l'operador `+`"

    **Primer, heu de saber que podem fusionar el contingut de dos mapes utilitzant l'operador Groovy `+`.**

    Diguem que tenim els següents mapes:

    ```groovy
    map1 = [id: 'sampleA', character: 'squirrel']
    map2 = [lang: 'fr']
    ```

    Els podem fusionar així:

    ```groovy
    new_map = map1 + map2
    ```

    El contingut de `new_map` serà:

    ```groovy
    [id: 'sampleA', character: 'squirrel', lang: 'fr']
    ```

    Genial!

    **Però què passa si necessiteu afegir un camp que encara no forma part d'un mapa?**

    Diguem que torneu a començar des de `map1`, però la predicció d'idioma no està al seu propi mapa (no hi ha `map2`).
    En canvi, es manté en una variable anomenada `lang_id`, i sabeu que voleu emmagatzemar el seu valor (`'fr'`) amb la clau `lang`.

    De fet, podeu fer el següent:

    ```groovy
    new_map = [map1 + [lang: lang_id]]
    ```

    Aquí, `[lang: new_info]` crea un nou mapa sense nom sobre la marxa, i `map1 + ` fusiona `map1` amb el nou mapa sense nom, produint el mateix contingut de `new_map` que abans.

    Genial, oi?

    **Ara transposem això al context d'una operació `channel.map()` de Nextflow.**

    El codi es converteix en:

    ```groovy
    .map { map1, lang_id ->
        [map1 + [lang: lang_id]]
    }
    ```

    Això fa el següent:

    - `map1, lang_id ->` pren els dos elements de la tupla
    - `[map1 + [lang: lang_id]]` crea el nou mapa tal com es detalla anteriorment

    La sortida és un únic mapa sense nom amb el mateix contingut que `new_map` al nostre exemple anterior.
    Així que hem transformat efectivament:

    ```groovy
    [id: 'sampleA', character: 'squirrel'], 'fr'
    ```

    en:

    ```groovy
    [id: 'sampleA', character: 'squirrel', lang: 'fr']
    ```

    Esperem que pugueu veure que si canviem `map1` a `meta`, això és bàsicament tot el que necessitem per afegir la predicció d'idioma al nostre mapa meta al nostre workflow.

    Excepte per una cosa!

    En el cas del nostre workflow, **també necessitem tenir en compte la presència de l'objecte `file` a la tupla**, que està composada de `meta, file, lang_id`.

    Així que el codi aquí es convertiria en:

    ```groovy
    .map { meta, file, lang_id ->
        [meta + [lang: lang_id], file]
    }
    ```

    Si teniu dificultats per seguir per què el `file` sembla estar movent-se a l'operació `map`, imagineu que en lloc de `[meta + [lang: lang_id], file]`, aquesta línia llegeix `[new_map, file]`.
    Això hauria de deixar més clar que simplement estem deixant el `file` al seu lloc original en segona posició a la tupla. Simplement hem pres el valor `new_info` i l'hem incorporat al mapa que està en primera posició.

    **I això ens porta de nou a l'estructura de canal `tuple val(meta), path(file)`!**

Un cop estigueu segurs que enteneu què fa aquest codi, executeu el workflow per veure si ha funcionat:

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

Sí, això quadra!
Hem reorganitzat acuradament la sortida del procés de `meta, file, lang_id` de manera que `lang_id` ara és una de les claus al mapa meta, i les tuples del canal s'ajusten de nou al model `meta, file`.

### 2.4. Assignar un grup lingüístic utilitzant condicionals

Ara que tenim les nostres prediccions d'idioma, utilitzem la informació per assignar algunes agrupacions noves.

A les nostres dades d'exemple, els idiomes utilitzats pels nostres personatges es poden agrupar en llengües germàniques (anglès, alemany) i llengües romàniques (francès, espanyol, italià).
Podria ser útil tenir aquesta classificació fàcilment disponible en algun lloc més endavant al pipeline, així que afegim aquesta informació al mapa meta.

I, bones notícies, aquest és un altre cas que es presta perfectament a utilitzar l'operador `map`!

Específicament, definirem una variable anomenada `lang_group`, utilitzarem una lògica condicional simple per determinar quin valor assignar a `lang_group` per a cada peça de dades.

La sintaxi general tindrà aquest aspecte:

```groovy
.map { meta, file ->

    // la lògica condicional que defineix lang_group va aquí

    [meta + [lang_group: lang_group], file]
}
```

Podeu veure que això és molt similar a l'operació de fusió de mapes sobre la marxa que vam utilitzar al pas anterior.
Només necessitem escriure les declaracions condicionals.

Aquí teniu la lògica condicional que volem aplicar:

- Definir una variable anomenada `lang_group` amb valor per defecte `'unknown'`.
- Si `lang` és alemany (`'de'`) o anglès (`'en'`), canviar `lang_group` a `germanic`.
- Si no, si `lang` està inclòs en una llista que conté francès (`'fr'`), espanyol (`'es'`) i italià (`'it'`), canviar `lang_group` a `romance`.

Proveu d'escriure-ho vosaltres mateixos si ja sabeu com escriure declaracions condicionals a Nextflow.

!!! tip "Consell"

    Podeu accedir al valor de `lang` dins de l'operació map amb `meta.lang`.

Hauríeu d'acabar fent els següents canvis al workflow:

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

- Utilitzem `def lang_group = "unknown"` per crear la variable `lang_group` amb valor per defecte establert a `unknown`.
- Utilitzem una estructura `if {} else if {}` per a la lògica condicional, amb proves alternatives `.equals()` per a les dues llengües germàniques, i una prova d'existència en una llista per a les tres llengües romàniques.
- Utilitzem l'operació de fusió `meta + [lang_group:lang_group]` com anteriorment per generar el mapa meta actualitzat.

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

Com podeu veure, els elements del canal mantenen la seva estructura `[meta, file]`, però el mapa meta ara inclou aquesta nova classificació.

### Conclusió

En aquesta secció, heu après com:

- **Aplicar metadades d'entrada a canals de sortida**: Copiar metadades d'aquesta manera ens permet associar resultats més endavant basant-nos en el contingut de les metadades.
- **Crear claus personalitzades**: Heu creat dues claus noves al vostre mapa meta, fusionant-les amb `meta + [new_key:value]` al mapa meta existent. Una basada en un valor calculat d'un procés, i una basada en una condició que heu establert a l'operador `map`.

Aquestes us permeten associar metadades noves i existents amb fitxers a mesura que progresseu pel vostre pipeline.
Fins i tot si no utilitzeu metadades com a part d'un procés, mantenir el mapa meta associat amb les dades d'aquesta manera facilita mantenir tota la informació rellevant junta.

---

## 3. Utilitzar informació del mapa meta en un procés

Ara que sabeu com crear i actualitzar el mapa meta, podem arribar a la part realment divertida: utilitzar realment les metadades en un procés.

Més específicament, afegirem un segon pas al nostre workflow per dibuixar cada animal com a art ASCII i fer-lo dir el text gravat en una bombolla de diàleg.
Ho farem utilitzant una eina anomenada [`cowpy`](https://github.com/jeffbuttars/cowpy).

??? info "Què fa `cowpy`?"

    `cowpy` és una eina de línia de comandes que genera art ASCII per mostrar entrades de text arbitràries d'una manera divertida.
    És una implementació en python de la clàssica eina [cowsay](https://en.wikipedia.org/wiki/Cowsay) de Tony Monroe.

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

Si vau treballar pel curs Hello Nextflow, ja heu vist aquesta eina en acció.
Si no, no us preocupeu; cobrirem tot el que necessiteu saber a mesura que avancem.

### 3.1. Importar el procés i examinar el codi

Us proporcionem un mòdul de procés pre-escrit anomenat `COWPY` que encapsula l'eina `cowpy`, així que només necessiteu afegir una declaració d'inclusió abans del bloc de workflow.

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

Com podeu veure, aquest procés està actualment dissenyat per prendre un fitxer d'entrada (que conté el text a mostrar) i un valor que especifica el personatge que s'hauria de dibuixar en art ASCII, normalment proporcionat al nivell de workflow per un paràmetre de línia de comandes.

### 3.2. Passar un camp del mapa meta com a entrada

Quan vam utilitzar l'eina `cowpy` al curs Hello Nextflow, vam utilitzar un paràmetre de línia de comandes per determinar quin personatge utilitzar per dibuixar la imatge final.
Això tenia sentit, perquè només estàvem generant una imatge per execució del pipeline.

No obstant això, en aquest tutorial, volem generar una imatge apropiada per a cada subjecte que estem processant, així que utilitzar un paràmetre de línia de comandes seria massa limitant.

Bones notícies: tenim una columna `character` al nostre full de dades i, per tant, al nostre mapa meta.
Utilitzem això per establir el personatge que el procés hauria d'utilitzar per a cada entrada.

Per a això, haurem de fer tres coses:

1. Donar un nom al canal de sortida que surt del procés anterior perquè puguem operar-hi més convenientment.
2. Determinar com accedir a la informació d'interès
3. Afegir una crida al segon procés i alimentar la informació adequadament.

Comencem.

#### 3.2.1. Nomenar el canal de sortida anterior

Vam aplicar les manipulacions anteriors directament sobre el canal de sortida del primer procés, `IDENTIFY_LANGUAGE.out`.
Per alimentar el contingut d'aquest canal al següent procés (i fer-ho d'una manera que sigui clara i fàcil de llegir) volem donar-li el seu propi nom, `ch_languages`.

Podem fer-ho utilitzant l'operador [`set`](https://www.nextflow.io/docs/latest/reference/operator.html#set).

Al workflow principal, reemplaceu l'operador `.view()` amb `.set { ch_languages }`, i afegiu una línia provant que podem referir-nos al canal pel nom.

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

        // Temporal: mirar dins de ch_languages
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

Executem això:

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

Això confirma que ara podem referir-nos al canal pel nom.

#### 3.2.2. Accedir al fitxer i les metadades del personatge

Sabem per mirar el codi del mòdul que el procés `COWPY` espera rebre un fitxer de text i un valor `character`.
Per escriure la crida al procés `COWPY`, només necessitem saber com extreure l'objecte fitxer corresponent i les metadades de cada element al canal.

Com sovint passa, la manera més senzilla de fer-ho és utilitzar una operació `map`.

El nostre canal conté tuples estructurades com `[meta, file]`, així que podem accedir a l'objecte `file` directament, i podem accedir al valor `character` emmagatzemat dins del mapa meta referint-nos-hi com `meta.character`.

Al workflow principal, feu els següents canvis de codi:

=== "Després"

    ```groovy title="main.nf" linenums="34"
        // Temporal: accedir al fitxer i al personatge
        ch_languages.map { meta, file -> file }.view { file -> "File: " + file }
        ch_languages.map { meta, file -> meta.character }.view { character -> "Character: " + character }
    ```

=== "Abans"

    ```groovy title="main.nf" linenums="34"
        // Temporal: mirar dins de ch_languages
        ch_languages.view()
    ```

Tingueu en compte que estem utilitzant closures (com ara `{ file -> "File: " + file }`) per fer la sortida de les operacions `.view` més llegible.

Executem això:

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

_Els camins de fitxer i els valors de personatge poden sortir en un ordre diferent a la vostra sortida._

Això confirma que som capaços d'accedir al fitxer i al personatge per a cada element al canal.

#### 3.2.3. Cridar el procés `COWPY`

Ara posem-ho tot plegat i cridem realment el procés `COWPY` sobre el canal `ch_languages`.

Al workflow principal, feu els següents canvis de codi:

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

Veieu que simplement copiem les dues operacions map (menys les declaracions `.view()`) com a entrades a la crida del procés.
Només assegureu-vos de no oblidar la coma entre elles!

És una mica feixuc, però veurem com fer-ho millor a la següent secció.

Executem això:

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

Si mireu al directori de resultats, hauríeu de veure els fitxers individuals que contenen l'art ASCII de cada salutació dita pel personatge corresponent.

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

Això mostra que vam poder utilitzar la informació al mapa meta per parametritzar la comanda al segon pas del pipeline.

No obstant això, com s'ha indicat anteriorment, part del codi implicat era una mica feixuc, ja que vam haver de desempaquetar metadades mentre encara estàvem al context del cos del workflow.
Aquest enfocament funciona bé per utilitzar un petit nombre de camps del mapa meta, però escalaria malament si volguéssim utilitzar molts més.

Hi ha un altre operador anomenat `multiMap()` que ens permet racionalitzar això una mica, però fins i tot així no és ideal.

??? info "(Opcional) Versió alternativa amb `multiMap()`"

    En cas que us ho estigueu preguntant, no podíem simplement escriure una única operació `map()` que produís tant el `file` com el `character`, perquè això els retornaria com una tupla.
    Vam haver d'escriure dues operacions `map()` separades per alimentar els elements `file` i `character` al procés per separat.

    Tècnicament hi ha una altra manera de fer-ho mitjançant una única operació de mapatge, utilitzant l'operador `multiMap()`, que és capaç d'emetre múltiples canals.
    Per exemple, podríeu reemplaçar la crida a `COWPY` anterior amb el següent codi:

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

En qualsevol cas, és incòmode que hàgim de fer algun desempaquetament al nivell de workflow.

Seria millor si poguéssim alimentar tot el mapa meta al procés i triar el que necessitem un cop allà.

### 3.3. Passar i utilitzar tot el mapa meta

El punt del mapa meta és, després de tot, passar totes les metadades juntes com un paquet.
L'única raó per la qual no podíem fer-ho anteriorment és que el procés no està configurat per acceptar un mapa meta.
Però com que controlem el codi del procés, podem canviar això.

Modifiquem el procés `COWPY` per acceptar l'estructura de tupla `[meta, file]` que vam utilitzar al primer procés perquè puguem racionalitzar el workflow.

Per a això, haurem de fer tres coses:

1. Modificar les definicions d'entrada del mòdul del procés `COWPY`
2. Actualitzar la comanda del procés per utilitzar el mapa meta
3. Actualitzar la crida del procés al cos del workflow

Preparats? Anem-hi!

#### 3.3.1. Modificar l'entrada del mòdul `COWPY`

Feu les següents edicions al fitxer del mòdul `cowpy.nf`:

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

Això ens permet utilitzar l'estructura de tupla `[meta, file]` que vam cobrir anteriorment al tutorial.

Tingueu en compte que no vam actualitzar la definició de sortida del procés per produir el mapa meta, per mantenir el tutorial racionalitzat, però sentiu-vos lliures de fer-ho vosaltres mateixos com a exercici seguint el model del procés `IDENTIFY_LANGUAGE`.

#### 3.3.2. Actualitzar la comanda per utilitzar el camp del mapa meta

Tot el mapa meta està ara disponible dins del procés, així que podem referir-nos a la informació que conté directament des de dins del bloc de comandes.

Feu les següents edicions al fitxer del mòdul `cowpy.nf`:

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

Hem reemplaçat la referència al valor `character` prèviament passat com a entrada independent amb el valor mantingut al mapa meta, al qual ens referim utilitzant `meta.character`.

Ara actualitzem la crida del procés en conseqüència.

#### 3.3.3. Actualitzar la crida del procés i executar-la

El procés ara espera que la seva entrada utilitzi l'estructura de tupla `[meta, file]`, que és el que el procés anterior produeix, així que podem simplement passar tot el canal `ch_languages` al procés `COWPY`.

Feu les següents edicions al workflow principal:

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

Això simplifica significativament la crida!

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

Si mireu al directori de resultats, hauríeu de veure les mateixes sortides que anteriorment, _és a dir_ fitxers individuals que contenen l'art ASCII de cada salutació dita pel personatge corresponent.

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

Així que això produeix els mateixos resultats que abans amb codi més senzill.

Per descomptat, això assumeix que podeu modificar el codi del procés.
En alguns casos, podeu haver de confiar en processos existents que no teniu llibertat de modificar, el que limita les vostres opcions.
La bona notícia, si teniu previst utilitzar mòduls del projecte [nf-core](https://nf-co.re/), és que els mòduls nf-core estan tots configurats per utilitzar l'estructura de tupla `[meta, file]` com a estàndard.

### 3.4. Resolució de problemes d'entrades requerides que falten

El valor `character` és necessari perquè el procés `COWPY` s'executi amb èxit.
Si no establim un valor per defecte per a ell en un fitxer de configuració, HEM de proporcionar un valor per a ell al full de dades.

**Què passa si no ho fem?**
Depèn del que contingui el full de dades d'entrada i quina versió del workflow estem executant.

#### 3.4.1. La columna character existeix però està buida

Diguem que esborrem el valor de character per a una de les entrades al nostre full de dades per simular un error de recollida de dades:

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

Per a qualsevol versió del workflow que hem utilitzat anteriorment, la clau `character` es crearà per a totes les entrades quan es llegeixi el full de dades, però per a `sampleA` el valor serà una cadena buida.

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

Quan Nextflow executa la línia de comandes `cowpy` per a aquesta mostra, `${meta.character}` s'omple amb una cadena buida a la línia de comandes `cowpy`, així que l'eina `cowpy` llança un error dient que no es va proporcionar cap valor per a l'argument `-c`.

#### 3.4.2. La columna character no existeix al full de dades

Ara diguem que esborrem la columna `character` completament del nostre full de dades:

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

##### 3.4.2.1. Valor accedit al nivell de workflow

Si estem utilitzant la versió del codi que vam escriure a la secció 3.2, Nextflow intentarà accedir a la clau `character` al mapa meta ABANS de cridar el procés `COWPY`.

No trobarà cap element que coincideixi amb la instrucció, així que no executarà `COWPY` en absolut.

??? success "Sortida de la comanda"

    ```console hl_lines="7"
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [desperate_montalcini] DSL2 - revision: 0dfeee3cc1

    executor >  local (7)
    [1a/df2544] process > IDENTIFY_LANGUAGE (7) [100%] 7 of 7 ✔
    [-        ] process > COWPY                 -
    ```

Pel que fa a Nextflow, aquest workflow s'ha executat amb èxit!
No obstant això, no es produirà cap de les sortides que volem.

##### 3.4.2.2. Valor accedit al nivell de procés

Si estem utilitzant la versió de la secció 3.3, Nextflow passarà tot el mapa meta al procés `COWPY` i intentarà executar la comanda.

Això causarà un error, però un de diferent comparat amb el primer cas.

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

Això passa perquè `meta.character` no existeix, així que el nostre intent d'accedir-hi retorna `null`. Com a resultat, Nextflow literalment connecta `null` a la línia de comandes, que per descomptat no és reconegut per l'eina `cowpy`.

#### 3.4.3. Solucions

A part de proporcionar un valor per defecte com a part de la configuració del workflow, hi ha dues coses que podem fer per gestionar això de manera més robusta:

1. Implementar validació d'entrada al vostre workflow per assegurar que el full de dades conté tota la informació requerida. Podeu trobar una [introducció a la validació d'entrada](../hello_nf-core/05_input_validation.md) al curs de formació Hello nf-core.

2. Si voleu assegurar-vos que qualsevol que utilitzi el vostre mòdul de procés pugui identificar immediatament les entrades requerides, també podeu fer que la propietat de metadades requerida sigui una entrada explícita.

Aquí teniu un exemple de com funcionaria això.

Primer, al nivell de procés, actualitzeu la definició d'entrada de la següent manera:

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

Després, al nivell de workflow, utilitzeu una operació de mapatge per extreure la propietat `character` de les metadades i fer-la un component explícit de la tupla d'entrada:

=== "Després"

    ```groovy title="main.nf" linenums="37" hl_lines="1"
        COWPY(ch_languages.map{meta, file -> [meta, meta.character, file]})
    ```

=== "Abans"

    ```groovy title="main.nf" linenums="37" hl_lines="1"
        COWPY(ch_languages)
    ```

Aquest enfocament té l'avantatge de mostrar explícitament que `character` és necessari, i fa que el procés sigui més fàcil de redesplegar en altres contextos.

Això destaca un principi de disseny important:

**Utilitzeu el mapa meta per a informació opcional i descriptiva, però extraieu valors requerits com a entrades explícites.**

El mapa meta és excel·lent per mantenir les estructures de canal netes i prevenir estructures de canal arbitràries, però per a elements obligatoris que es referencien directament en un procés, extreure'ls com a entrades explícites crea codi més robust i mantenible.

### Conclusió

En aquesta secció, heu après com utilitzar metadades per personalitzar l'execució d'un procés, accedint-hi ja sigui al nivell de workflow o al nivell de procés.

---

## Exercici suplementari

Si voleu practicar l'ús d'informació del mapa meta des de dins d'un procés, proveu d'utilitzar altres peces d'informació del mapa meta com ara `lang` i `lang_group` per personalitzar com es nomenen i/o organitzen les sortides.

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

---

## Resum

En aquesta missió secundària, heu explorat com treballar efectivament amb metadades en workflows de Nextflow.

Aquest patró de mantenir les metadades explícites i adjuntes a les dades és una pràctica recomanada central a Nextflow, oferint diversos avantatges sobre codificar informació de fitxers:

- Les metadades dels fitxers romanen associades amb els fitxers al llarg del workflow
- El comportament del procés es pot personalitzar per fitxer
- L'organització de sortides pot reflectir les metadades dels fitxers
- La informació dels fitxers es pot expandir durant l'execució del pipeline

Aplicar aquest patró al vostre propi treball us permetrà construir workflows de bioinformàtica robustos i mantenibles.

### Patrons clau

1.  **Llegir i Estructurar Metadades:** Llegir fitxers CSV i crear mapes de metadades organitzats que romanen associats amb els vostres fitxers de dades.

    ```groovy
    channel.fromPath('datasheet.csv')
      .splitCsv(header: true)
      .map { row ->
          [ [id:row.id, character:row.character], row.recording ]
      }
    ```

2.  **Expandir Metadades Durant el Workflow** Afegir nova informació a les vostres metadades a mesura que el vostre pipeline progressa afegint sortides de processos i derivant valors mitjançant lògica condicional.

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

3.  **Personalitzar el Comportament del Procés:** Utilitzar metadades dins del procés.

    ```groovy
    cat $input_file | cowpy -c ${meta.character} > cowpy-${input_file}
    ```

### Recursos addicionals

- [map](https://www.nextflow.io/docs/latest/operator.html#map)
- [stdout](https://www.nextflow.io/docs/latest/process.html#outputs)

---

## Què segueix?

Torneu al [menú de Missions Secundàries](./index.md) o feu clic al botó a la part inferior dreta de la pàgina per passar al següent tema de la llista.
