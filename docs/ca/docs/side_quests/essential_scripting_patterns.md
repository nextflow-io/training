# Patrons Essencials de Scripting en Nextflow

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducció assistida per IA - [més informació i suggeriments](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Nextflow és un llenguatge de programació que s'executa a la Màquina Virtual de Java. Tot i que Nextflow està construït sobre [Groovy](http://groovy-lang.org/) i comparteix gran part de la seva sintaxi, Nextflow és més que "Groovy amb extensions" -- és un llenguatge independent amb una [sintaxi](https://nextflow.io/docs/latest/reference/syntax.html) i una [biblioteca estàndard](https://nextflow.io/docs/latest/reference/stdlib.html) completament especificades.

Podeu escriure molt codi Nextflow sense aventurar-vos més enllà de la sintaxi bàsica per a variables, mapes i llistes. La majoria de tutorials de Nextflow se centren en l'orquestració de workflows (canals, processos i flux de dades), i podeu arribar sorprenentment lluny només amb això.

No obstant això, quan necessiteu manipular dades, analitzar noms de fitxers complexos, implementar lògica condicional o construir workflows de producció robustos, ajuda pensar en dos aspectes diferents del vostre codi: **flux de dades** (canals, operadors, processos i workflows) i **scripting** (el codi dins de closures, funcions i scripts de processos). Tot i que aquesta distinció és una mica arbitrària—tot és codi Nextflow—proporciona un model mental útil per entendre quan esteu orquestrant el vostre pipeline versus quan esteu manipulant dades. Dominar ambdós aspectes millora dràsticament la vostra capacitat per escriure workflows clars i mantenibles.

### Objectius d'aprenentatge

Aquesta missió secundària us porta en un viatge pràctic des de conceptes bàsics fins a patrons preparats per a producció.
Transformarem un workflow simple de lectura de CSV en un pipeline bioinformàtic sofisticat, evolucionant-lo pas a pas a través de reptes realistes:

- **Comprendre els límits:** Distingir entre operacions de flux de dades i scripting, i entendre com treballen junts
- **Manipulació de dades:** Extreure, transformar i subconjuntar mapes i col·leccions utilitzant operadors potents
- **Processament de cadenes:** Analitzar esquemes complexos de nomenclatura de fitxers amb patrons regex i dominar la interpolació de variables
- **Funcions reutilitzables:** Extreure lògica complexa en funcions amb nom per a workflows més nets i mantenibles
- **Lògica dinàmica:** Construir processos que s'adaptin a diferents tipus d'entrada i utilitzar closures per a l'assignació dinàmica de recursos
- **Enrutament condicional:** Enrutar intel·ligentment mostres a través de diferents processos basant-se en les seves característiques de metadades
- **Operacions segures:** Gestionar dades absents amb gràcia amb operadors segurs contra nuls i validar entrades amb missatges d'error clars
- **Gestors basats en configuració:** Utilitzar gestors d'esdeveniments de workflow per a registre, notificacions i gestió del cicle de vida

### Prerequisits

Abans d'emprendre aquesta missió secundària, hauríeu de:

- Haver completat el tutorial [Hello Nextflow](../hello_nextflow/README.md) o un curs equivalent per a principiants.
- Sentir-vos còmodes utilitzant conceptes i mecanismes bàsics de Nextflow (processos, canals, operadors, treballar amb fitxers, metadades)
- Tenir familiaritat bàsica amb construccions de programació comunes (variables, mapes, llistes)

Aquest tutorial explicarà conceptes de programació a mesura que els trobem, així que no necessiteu una experiència extensa en programació.
Començarem amb conceptes fonamentals i construirem fins a patrons avançats.

---

## 0. Primers passos

#### Obrir el codespace de formació

Si encara no ho heu fet, assegureu-vos d'obrir l'entorn de formació tal com es descriu a [Configuració de l'entorn](../envsetup/index.md).

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### Moure's al directori del projecte

Movem-nos al directori on es troben els fitxers per a aquest tutorial.

```bash
cd side-quests/essential_scripting_patterns
```

#### Revisar els materials

Trobareu un fitxer de workflow principal i un directori `data` que conté fitxers de dades d'exemple.

```console title="Contingut del directori"
.
├── collect.nf
├── data
│   ├── samples.csv
│   └── sequences
│       ├── SAMPLE_001_S1_L001_R1_001.fastq
│       ├── SAMPLE_002_S2_L001_R1_001.fastq
│       └── SAMPLE_003_S3_L001_R1_001.fastq
├── main.nf
├── modules
│   ├── fastp.nf
│   ├── generate_report.nf
│   └── trimgalore.nf
└── nextflow.config
```

El nostre CSV de mostra conté informació sobre mostres biològiques que necessiten diferents processaments basant-se en les seves característiques:

```console title="samples.csv"
sample_id,organism,tissue_type,sequencing_depth,file_path,quality_score
SAMPLE_001,human,liver,30000000,data/sequences/SAMPLE_001_S1_L001_R1_001.fastq,38.5
SAMPLE_002,mouse,brain,25000000,data/sequences/SAMPLE_002_S2_L001_R1_001.fastq,35.2
SAMPLE_003,human,kidney,45000000,data/sequences/SAMPLE_003_S3_L001_R1_001.fastq,42.1
```

Utilitzarem aquest conjunt de dades realista per explorar tècniques de programació pràctiques que trobareu en workflows bioinformàtics reals.

#### Llista de verificació de preparació

Creieu que esteu preparats per començar?

- [ ] Entenc l'objectiu d'aquest curs i els seus prerequisits
- [ ] El meu codespace està en funcionament
- [ ] He establert el meu directori de treball adequadament

Si podeu marcar totes les caselles, esteu preparats per començar.

---

## 1. Flux de Dades vs Scripting: Comprendre els Límits

### 1.1. Identificar Què és Què

Quan escriviu workflows de Nextflow, és important distingir entre **flux de dades** (com es mouen les dades a través de canals i processos) i **scripting** (el codi que manipula dades i pren decisions). Construïm un workflow que demostra com treballen junts.

#### 1.1.1. Workflow Bàsic de Nextflow

Comencem amb un workflow simple que només llegeix el fitxer CSV (ja ho hem fet per vosaltres a `main.nf`):

```groovy title="main.nf" linenums="1"
workflow {
    ch_samples = channel.fromPath("./data/samples.csv")
        .splitCsv(header: true)
        .view()
}
```

El bloc `workflow` defineix l'estructura del nostre pipeline, mentre que `channel.fromPath()` crea un canal des d'un camí de fitxer. L'operador `.splitCsv()` processa el fitxer CSV i converteix cada fila en una estructura de dades de tipus mapa.

Executeu aquest workflow per veure les dades CSV en brut:

```bash
nextflow run main.nf
```

??? success "Sortida de la comanda"

    ```console
    Launching `main.nf` [marvelous_tuckerman] DSL2 - revision: 6113e05c17

    [sample_id:SAMPLE_001, organism:human, tissue_type:liver, sequencing_depth:30000000, file_path:data/sequences/SAMPLE_001_S1_L001_R1_001.fastq, quality_score:38.5]
    [sample_id:SAMPLE_002, organism:mouse, tissue_type:brain, sequencing_depth:25000000, file_path:data/sequences/SAMPLE_002_S2_L001_R1_001.fastq, quality_score:35.2]
    [sample_id:SAMPLE_003, organism:human, tissue_type:kidney, sequencing_depth:45000000, file_path:data/sequences/SAMPLE_003_S3_L001_R1_001.fastq, quality_score:42.1]
    ```

#### 1.1.2. Afegir l'Operador Map

Ara afegirem scripting per transformar les dades, utilitzant l'operador `.map()` que probablement ja coneixeu. Aquest operador pren una 'closure' on podem escriure codi per transformar cada element.

!!! note "Nota"

    Una **closure** és un bloc de codi que es pot passar i executar més tard. Penseu-hi com una funció que definiu en línia. Les closures s'escriuen amb claus `{ }` i poden prendre paràmetres. Són fonamentals per a com funcionen els operadors de Nextflow i si heu estat escrivint Nextflow durant un temps, potser ja les heu estat utilitzant sense adonar-vos-en!

Així és com es veu aquesta operació map:

=== "Després"

    ```groovy title="main.nf" linenums="2" hl_lines="3-6"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row ->
                return row
            }
            .view()
    ```

=== "Abans"

    ```groovy title="main.nf" linenums="2" hl_lines="3"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .view()
    ```

Aquesta és la nostra primera **closure** - una funció anònima que podeu passar com a argument (similar a lambdas en Python o funcions fletxa en JavaScript). Les closures són essencials per treballar amb operadors de Nextflow.

La closure `{ row -> return row }` pren un paràmetre `row` (podria ser qualsevol nom: `item`, `sample`, etc.).

Quan l'operador `.map()` processa cada element del canal, passa aquest element a la vostra closure. Aquí, `row` conté una fila CSV a la vegada.

Apliqueu aquest canvi i executeu el workflow:

```bash
nextflow run main.nf
```

Veureu la mateixa sortida que abans, perquè simplement estem retornant l'entrada sense canvis. Això confirma que l'operador map està funcionant correctament. Ara comencem a transformar les dades.

#### 1.1.3. Crear una Estructura de Dades de Tipus Mapa

Ara escriurem lògica de **scripting** dins de la nostra closure per transformar cada fila de dades. Aquí és on processem elements de dades individuals en lloc d'orquestrar el flux de dades.

=== "Després"

    ```groovy title="main.nf" linenums="2" hl_lines="4-12"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row ->
                // Scripting per a transformació de dades
                def sample_meta = [
                    id: row.sample_id.toLowerCase(),
                    organism: row.organism,
                    tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                    depth: row.sequencing_depth.toInteger(),
                    quality: row.quality_score.toDouble()
                ]
                return sample_meta
            }
            .view()
    ```

=== "Abans"

    ```groovy title="main.nf" linenums="2" hl_lines="4"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row ->
                return row
            }
            .view()
    ```

El mapa `sample_meta` és una estructura de dades clau-valor (com diccionaris en Python, objectes en JavaScript o hashes en Ruby) que emmagatzema informació relacionada: ID de mostra, organisme, tipus de teixit, profunditat de seqüenciació i puntuació de qualitat.

Utilitzem mètodes de manipulació de cadenes com `.toLowerCase()` i `.replaceAll()` per netejar les nostres dades, i mètodes de conversió de tipus com `.toInteger()` i `.toDouble()` per convertir dades de cadena del CSV en els tipus numèrics apropiats.

Apliqueu aquest canvi i executeu el workflow:

```bash
nextflow run main.nf
```

??? success "Sortida de la comanda"

    ```console
    [id:sample_001, organism:human, tissue:liver, depth:30000000, quality:38.5]
    [id:sample_002, organism:mouse, tissue:brain, depth:25000000, quality:35.2]
    [id:sample_003, organism:human, tissue:kidney, depth:45000000, quality:42.1]
    ```

#### 1.1.4. Afegir Lògica Condicional

Ara afegirem més scripting - aquesta vegada utilitzant un operador ternari per prendre decisions basades en valors de dades.

Feu el següent canvi:

=== "Després"

    ```groovy title="main.nf" linenums="2" hl_lines="11-12"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row ->
                def sample_meta = [
                    id: row.sample_id.toLowerCase(),
                    organism: row.organism,
                    tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                    depth: row.sequencing_depth.toInteger(),
                    quality: row.quality_score.toDouble()
                ]
                def priority = sample_meta.quality > 40 ? 'high' : 'normal'
                return sample_meta + [priority: priority]
            }
            .view()
    ```

=== "Abans"

    ```groovy title="main.nf" linenums="2" hl_lines="11"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row ->
                def sample_meta = [
                    id: row.sample_id.toLowerCase(),
                    organism: row.organism,
                    tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                    depth: row.sequencing_depth.toInteger(),
                    quality: row.quality_score.toDouble()
                ]
                return sample_meta
            }
            .view()
    ```

L'operador ternari és una forma abreujada d'una sentència if/else que segueix el patró `condició ? valor_si_cert : valor_si_fals`. Aquesta línia significa: "Si la qualitat és superior a 40, utilitza 'high', altrament utilitza 'normal'". El seu cosí, l'**operador Elvis** (`?:`), proporciona valors per defecte quan alguna cosa és nul·la o buida - explorarem aquest patró més endavant en aquest tutorial.

L'operador d'addició de mapes `+` crea un **nou mapa** en lloc de modificar l'existent. Aquesta línia crea un nou mapa que conté tots els parells clau-valor de `sample_meta` més la nova clau `priority`.

!!! Note "Nota"

    Mai modifiqueu mapes passats a closures - sempre creeu-ne de nous utilitzant `+` (per exemple). A Nextflow, les mateixes dades sovint flueixen a través de múltiples operacions simultàniament. Modificar un mapa in situ pot causar efectes secundaris impredictibles quan altres operacions referencien aquest mateix objecte. Crear nous mapes assegura que cada operació té la seva pròpia còpia neta.

Executeu el workflow modificat:

```bash
nextflow run main.nf
```

??? success "Sortida de la comanda"

    ```console
    [id:sample_001, organism:human, tissue:liver, depth:30000000, quality:38.5, priority:normal]
    [id:sample_002, organism:mouse, tissue:brain, depth:25000000, quality:35.2, priority:normal]
    [id:sample_003, organism:human, tissue:kidney, depth:45000000, quality:42.1, priority:high]
    ```

Hem afegit amb èxit lògica condicional per enriquir les nostres metadades amb un nivell de prioritat basat en puntuacions de qualitat.

#### 1.1.5. Subconjuntar Mapes amb `.subMap()`

Mentre que l'operador `+` afegeix claus a un mapa, de vegades necessiteu fer el contrari - extreure només claus específiques. El mètode `.subMap()` és perfecte per a això.

Afegim una línia per crear una versió simplificada de les nostres metadades que només conté camps d'identificació:

=== "Després"

    ```groovy title="main.nf" linenums="2" hl_lines="12-15"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row ->
                // Scripting per a transformació de dades
                def sample_meta = [
                    id: row.sample_id.toLowerCase(),
                    organism: row.organism,
                    tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                    depth: row.sequencing_depth.toInteger(),
                    quality: row.quality_score.toDouble()
                ]
                def id_only = sample_meta.subMap(['id', 'organism', 'tissue'])
                println "Només camps d'ID: ${id_only}"

                def priority = sample_meta.quality > 40 ? 'high' : 'normal'
                return sample_meta + [priority: priority]
            }
            .view()
    ```

=== "Abans"

    ```groovy title="main.nf" linenums="2" hl_lines="12"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row ->
                // Scripting per a transformació de dades
                def sample_meta = [
                    id: row.sample_id.toLowerCase(),
                    organism: row.organism,
                    tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                    depth: row.sequencing_depth.toInteger(),
                    quality: row.quality_score.toDouble()
                ]
                def priority = sample_meta.quality > 40 ? 'high' : 'normal'
                return sample_meta + [priority: priority]
            }
            .view()
    ```

Executeu el workflow modificat:

```bash
nextflow run main.nf
```

??? success "Sortida de la comanda"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [peaceful_cori] DSL2 - revision: 4cc4a8340f

    Només camps d'ID: [id:sample_001, organism:human, tissue:liver]
    Només camps d'ID: [id:sample_002, organism:mouse, tissue:brain]
    Només camps d'ID: [id:sample_003, organism:human, tissue:kidney]
    [id:sample_001, organism:human, tissue:liver, depth:30000000, quality:38.5, priority:normal]
    [id:sample_002, organism:mouse, tissue:brain, depth:25000000, quality:35.2, priority:normal]
    [id:sample_003, organism:human, tissue:kidney, depth:45000000, quality:42.1, priority:high]
    ```

Això mostra tant les metadades completes mostrades per l'operació `view()` com el subconjunt extret que hem imprès amb `println`.

El mètode `.subMap()` pren una llista de claus i retorna un nou mapa que conté només aquestes claus. Si una clau no existeix al mapa original, simplement no s'inclou al resultat.

Això és particularment útil quan necessiteu crear diferents versions de metadades per a diferents processos - alguns poden necessitar metadades completes mentre que altres necessiten només camps d'identificació mínims.

Ara elimineu aquestes sentències println per restaurar el vostre workflow al seu estat anterior, ja que no les necessitem d'ara endavant.

!!! tip "Resum d'Operacions amb Mapes"

    - **Afegir claus**: `map1 + [new_key: value]` - Crea un nou mapa amb claus addicionals
    - **Extreure claus**: `map1.subMap(['key1', 'key2'])` - Crea un nou mapa amb només les claus especificades
    - **Ambdues operacions creen nous mapes** - Els mapes originals romanen sense canvis

#### 1.1.6. Combinar Mapes i Retornar Resultats

Fins ara, només hem estat retornant el que la comunitat Nextflow anomena el 'mapa meta', i hem estat ignorant els fitxers als quals aquestes metadades es refereixen. Però si esteu escrivint workflows de Nextflow, probablement voleu fer alguna cosa amb aquests fitxers.

Generem una estructura de canal que comprèn una tupla de 2 elements: el mapa de metadades enriquit i el camí de fitxer corresponent. Aquest és un patró comú a Nextflow per passar dades a processos.

=== "Després"

    ```groovy title="main.nf" linenums="2" hl_lines="12"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row ->
                def sample_meta = [
                    id: row.sample_id.toLowerCase(),
                    organism: row.organism,
                    tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                    depth: row.sequencing_depth.toInteger(),
                    quality: row.quality_score.toDouble()
                ]
                def priority = sample_meta.quality > 40 ? 'high' : 'normal'
                return tuple( sample_meta + [priority: priority], file(row.file_path) )
            }
            .view()
    ```

=== "Abans"

    ```groovy title="main.nf" linenums="2" hl_lines="12"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row ->
                def sample_meta = [
                    id: row.sample_id.toLowerCase(),
                    organism: row.organism,
                    tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                    depth: row.sequencing_depth.toInteger(),
                    quality: row.quality_score.toDouble()
                ]
                def priority = sample_meta.quality > 40 ? 'high' : 'normal'
                return sample_meta + [priority: priority]
            }
            .view()
    ```

Apliqueu aquest canvi i executeu el workflow:

```bash
nextflow run main.nf
```

??? success "Sortida de la comanda"

    ```console
    [[id:sample_001, organism:human, tissue:liver, depth:30000000, quality:38.5, priority:normal], /workspaces/training/side-quests/essential_scripting_patterns/data/sequences/SAMPLE_001_S1_L001_R1_001.fastq]
    [[id:sample_002, organism:mouse, tissue:brain, depth:25000000, quality:35.2, priority:normal], /workspaces/training/side-quests/essential_scripting_patterns/data/sequences/SAMPLE_002_S2_L001_R1_001.fastq]
    [[id:sample_003, organism:human, tissue:kidney, depth:45000000, quality:42.1, priority:high], /workspaces/training/side-quests/essential_scripting_patterns/data/sequences/SAMPLE_003_S3_L001_R1_001.fastq]
    ```

Aquesta estructura de tupla `[meta, file]` és un patró comú a Nextflow per passar tant metadades com fitxers associats a processos.

!!! note "Nota"

    **Mapes i Metadades**: Els mapes són fonamentals per treballar amb metadades a Nextflow. Per a una explicació més detallada sobre com treballar amb mapes de metadades, consulteu la missió secundària [Treballar amb metadades](./metadata.md).

El nostre workflow demostra el patró central: les **operacions de flux de dades** (`workflow`, `channel.fromPath()`, `.splitCsv()`, `.map()`, `.view()`) orquestren com es mouen les dades a través del pipeline, mentre que el **scripting** (mapes `[key: value]`, mètodes de cadenes, conversions de tipus, operadors ternaris) dins de la closure `.map()` gestiona la transformació d'elements de dades individuals.

### 1.2. Comprendre Diferents Tipus: Channel vs List

Fins ara, molt bé, podem distingir entre operacions de flux de dades i scripting. Però què passa quan el mateix nom de mètode existeix en ambdós contextos?

Un exemple perfecte és el mètode `collect`, que existeix tant per a tipus de canal com per a tipus List a la biblioteca estàndard de Nextflow. El mètode `collect()` en una List transforma cada element, mentre que l'operador `collect()` en un canal agrupa totes les emissions del canal en un canal d'un sol element.

Demostrem això amb algunes dades de mostra, començant per refrescar-nos sobre què fa l'operador `collect()` del canal. Consulteu `collect.nf`:

```groovy title="collect.nf" linenums="1"
def sample_ids = ['sample_001', 'sample_002', 'sample_003']

// channel.collect() - agrupa múltiples emissions de canal en una
ch_input = channel.fromList(sample_ids)
ch_input.view { sample -> "Element individual del canal: ${sample}" }
ch_collected = ch_input.collect()
ch_collected.view { list -> "Resultat de channel.collect(): ${list} (${list.size()} elements agrupats en 1)" }
```

Passos:

- Definir una List d'IDs de mostra
- Crear un canal amb `fromList()` que emet cada ID de mostra per separat
- Imprimir cada element amb `view()` a mesura que flueix
- Agrupar tots els elements en una sola llista amb l'operador `collect()` del canal
- Imprimir el resultat agrupat (un sol element que conté tots els IDs de mostra) amb un segon `view()`

Hem canviat l'estructura del canal, però no hem canviat les dades en si.

Executeu el workflow per confirmar-ho:

```bash
nextflow run collect.nf
```

??? success "Sortida de la comanda"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `collect.nf` [loving_mendel] DSL2 - revision: e8d054a46e

    Element individual del canal: sample_001
    Element individual del canal: sample_002
    Element individual del canal: sample_003
    Resultat de channel.collect(): [sample_001, sample_002, sample_003] (3 elements agrupats en 1)
    ```

`view()` retorna una sortida per cada emissió del canal, així que sabem que aquesta única sortida conté tots els 3 elements originals agrupats en una llista.

Ara vegem el mètode `collect` en una List en acció. Modifiqueu `collect.nf` per aplicar el mètode `collect` de la List a la llista original d'IDs de mostra:

=== "Després"

    ```groovy title="main.nf" linenums="1" hl_lines="9-13"
    def sample_ids = ['sample_001', 'sample_002', 'sample_003']

    // channel.collect() - agrupa múltiples emissions de canal en una
    ch_input = channel.fromList(sample_ids)
    ch_input.view { sample -> "Element individual del canal: ${sample}" }
    ch_collected = ch_input.collect()
    ch_collected.view { list -> "Resultat de channel.collect(): ${list} (${list.size()} elements agrupats en 1)" }

    // List.collect() - transforma cada element, preserva l'estructura
    def formatted_ids = sample_ids.collect { id ->
        id.toUpperCase().replace('SAMPLE_', 'SPECIMEN_')
    }
    println "Resultat de List.collect(): ${formatted_ids} (${sample_ids.size()} elements transformats en ${formatted_ids.size()})"
    ```

=== "Abans"

    ```groovy title="main.nf" linenums="1"
    def sample_ids = ['sample_001', 'sample_002', 'sample_003']

    // channel.collect() - agrupa múltiples emissions de canal en una
    ch_input = channel.fromList(sample_ids)
    ch_input.view { sample -> "Element individual del canal: ${sample}" }
    ch_collected = ch_input.collect()
    ch_collected.view { list -> "Resultat de channel.collect(): ${list} (${list.size()} elements agrupats en 1)" }
    ```

En aquest nou fragment:

- Definim una nova variable `formatted_ids` que utilitza el mètode `collect` de la List per transformar cada ID de mostra a la llista original
- Imprimim el resultat utilitzant `println`

Executeu el workflow modificat:

```bash
nextflow run collect.nf
```

??? success "Sortida de la comanda"

    ```console hl_lines="5"
    N E X T F L O W   ~  version 25.10.2

    Launching `collect.nf` [cheeky_stonebraker] DSL2 - revision: 2d5039fb47

    Resultat de List.collect(): [SPECIMEN_001, SPECIMEN_002, SPECIMEN_003] (3 elements transformats en 3)
    Element individual del canal: sample_001
    Element individual del canal: sample_002
    Element individual del canal: sample_003
    Resultat de channel.collect(): [sample_001, sample_002, sample_003] (3 elements agrupats en 1)
    ```

Aquesta vegada, NO hem canviat l'estructura de les dades, encara tenim 3 elements a la llista, però SÍ hem transformat cada element utilitzant el mètode `collect` de la List per produir una nova llista amb valors modificats. Això és similar a utilitzar l'operador `map` en un canal, però està operant sobre una estructura de dades List en lloc d'un canal.

`collect` és un cas extrem que estem utilitzant aquí per fer un punt. La lliçó clau és que quan escriviu workflows, sempre distingiu entre **estructures de dades** (Lists, Maps, etc.) i **canals** (construccions de flux de dades). Les operacions poden compartir noms però comportar-se completament diferent depenent del tipus sobre el qual es criden.

### 1.3. L'Operador Spread (`*.`) - Forma Abreujada per a Extracció de Propietats

Relacionat amb el mètode `collect` de la List està l'operador spread (`*.`), que proporciona una manera concisa d'extreure propietats de col·leccions. És essencialment sucre sintàctic per a un patró `collect` comú.

Afegim una demostració al nostre fitxer `collect.nf`:

=== "Després"

    ```groovy title="collect.nf" linenums="1" hl_lines="15-18"
    def sample_ids = ['sample_001', 'sample_002', 'sample_003']

    // channel.collect() - agrupa múltiples emissions de canal en una
    ch_input = channel.fromList(sample_ids)
    ch_input.view { sample -> "Element individual del canal: ${sample}" }
    ch_collected = ch_input.collect()
    ch_collected.view { list -> "Resultat de channel.collect(): ${list} (${list.size()} elements agrupats en 1)" }

    // List.collect() - transforma cada element, preserva l'estructura
    def formatted_ids = sample_ids.collect { id ->
        id.toUpperCase().replace('SAMPLE_', 'SPECIMEN_')
    }
    println "Resultat de List.collect(): ${formatted_ids} (${sample_ids.size()} elements transformats en ${formatted_ids.size()})"

    // Operador spread - accés concís a propietats
    def sample_data = [[id: 's1', quality: 38.5], [id: 's2', quality: 42.1], [id: 's3', quality: 35.2]]
    def all_ids = sample_data*.id
    println "Resultat de l'operador spread: ${all_ids}"
    ```

=== "Abans"

    ```groovy title="collect.nf" linenums="1"
    def sample_ids = ['sample_001', 'sample_002', 'sample_003']

    // channel.collect() - agrupa múltiples emissions de canal en una
    ch_input = channel.fromList(sample_ids)
    ch_input.view { sample -> "Element individual del canal: ${sample}" }
    ch_collected = ch_input.collect()
    ch_collected.view { list -> "Resultat de channel.collect(): ${list} (${list.size()} elements agrupats en 1)" }

    // List.collect() - transforma cada element, preserva l'estructura
    def formatted_ids = sample_ids.collect { id ->
        id.toUpperCase().replace('SAMPLE_', 'SPECIMEN_')
    }
    println "Resultat de List.collect(): ${formatted_ids} (${sample_ids.size()} elements transformats en ${formatted_ids.size()})"
    ```

Executeu el workflow actualitzat:

```bash title="Provar l'operador spread"
nextflow run collect.nf
```

??? success "Sortida de la comanda"

    ```console hl_lines="6"
    N E X T F L O W   ~  version 25.10.2

    Launching `collect.nf` [cranky_galileo] DSL2 - revision: 5f3c8b2a91

    Resultat de List.collect(): [SPECIMEN_001, SPECIMEN_002, SPECIMEN_003] (3 elements transformats en 3)
    Resultat de l'operador spread: [s1, s2, s3]
    Element individual del canal: sample_001
    Element individual del canal: sample_002
    Element individual del canal: sample_003
    Resultat de channel.collect(): [sample_001, sample_002, sample_003] (3 elements agrupats en 1)
    ```

L'operador spread `*.` és una forma abreujada per a un patró collect comú:

```groovy
// Aquests són equivalents:
def ids = samples*.id
def ids = samples.collect { it.id }

// També funciona amb crides a mètodes:
def names = files*.getName()
def names = files.collect { it.getName() }
```

L'operador spread és particularment útil quan necessiteu extreure una sola propietat d'una llista d'objectes - és més llegible que escriure la closure `collect` completa.

!!! tip "Quan Utilitzar Spread vs Collect"

    - **Utilitzeu spread (`*.`)** per a accés simple a propietats: `samples*.id`, `files*.name`
    - **Utilitzeu collect** per a transformacions o lògica complexa: `samples.collect { it.id.toUpperCase() }`, `samples.collect { [it.id, it.quality > 40] }`

### Conclusió

En aquesta secció, heu après:

- **Flux de dades vs scripting**: Els operadors de canal orquestren com flueixen les dades a través del vostre pipeline, mentre que el scripting transforma elements de dades individuals
- **Comprendre tipus**: El mateix nom de mètode (com `collect`) pot comportar-se diferent depenent del tipus sobre el qual es crida (Channel vs List)
- **El context importa**: Sempre sigueu conscients de si esteu treballant amb canals (flux de dades) o estructures de dades (scripting)

Comprendre aquests límits és essencial per a la depuració, documentació i escriptura de workflows mantenibles.

A continuació, aprofundirem en les capacitats de processament de cadenes, que són essencials per gestionar dades del món real.

---

## 2. Processament de Cadenes i Generació Dinàmica de Scripts

Dominar el processament de cadenes separa els workflows fràgils dels pipelines robustos. Aquesta secció cobreix l'anàlisi de noms de fitxers complexos, la generació dinàmica de scripts i la interpolació de variables.

### 2.1. Coincidència de Patrons i Expressions Regulars

Els fitxers bioinformàtics sovint tenen convencions de nomenclatura complexes que codifiquen metadades. Extraiem això automàticament utilitzant coincidència de patrons amb expressions regulars.

Tornarem al nostre workflow `main.nf` i afegirem alguna lògica de coincidència de patrons per extreure informació addicional de mostra dels noms de fitxers. Els fitxers FASTQ al nostre conjunt de dades segueixen convencions de nomenclatura d'estil Illumina amb noms com `SAMPLE_001_S1_L001_R1_001.fastq.gz`. Aquests poden semblar críptics, però en realitat codifiquen metadades útils com l'ID de mostra, el número de carril i la direcció de lectura. Utilitzarem capacitats regex per analitzar aquests noms.

Feu el següent canvi al vostre workflow `main.nf` existent:

=== "Després"

    ```groovy title="main.nf" linenums="4" hl_lines="10-21"
            .map { row ->
                // Scripting per a transformació de dades
                def sample_meta = [
                    id: row.sample_id.toLowerCase(),
                    organism: row.organism,
                    tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                    depth: row.sequencing_depth.toInteger(),
                    quality: row.quality_score.toDouble()
                ]
                def fastq_path = file(row.file_path)

                def m = (fastq_path.name =~ /^(.+)_S(\d+)_L(\d{3})_(R[12])_(\d{3})\.fastq(?:\.gz)?$/)
                def file_meta = m ? [
                    sample_num: m[0][2].toInteger(),
                    lane: m[0][3],
                    read: m[0][4],
                    chunk: m[0][5]
                ] : [:]

                def priority = sample_meta.quality > 40 ? 'high' : 'normal'
                return tuple(sample_meta + file_meta + [priority: priority], fastq_path)
            }
    ```

=== "Abans"

    ```groovy title="main.nf" linenums="4" hl_lines="10-11"
            .map { row ->
                // Scripting per a transformació de dades
                def sample_meta = [
                    id: row.sample_id.toLowerCase(),
                    organism: row.organism,
                    tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                    depth: row.sequencing_depth.toInteger(),
                    quality: row.quality_score.toDouble()
                ]
                def priority = sample_meta.quality > 40 ? 'high' : 'normal'
                return tuple(sample_meta + [priority: priority], file(row.file_path))
            }
    ```

Això demostra **conceptes clau de processament de cadenes**:

1. **Literals d'expressions regulars** utilitzant la sintaxi `~/patró/` - això crea un patró regex sense necessitat d'escapar barres invertides
2. **Coincidència de patrons** amb l'operador `=~` - això intenta fer coincidir una cadena amb un patró regex
3. **Objectes Matcher** que capturen grups amb `[0][1]`, `[0][2]`, etc. - `[0]` es refereix a la coincidència completa, `[1]`, `[2]`, etc. es refereixen a grups capturats entre parèntesis

Desglossem el patró regex `^(.+)_S(\d+)_L(\d{3})_(R[12])_(\d{3})\.fastq(?:\.gz)?$`:

| Patró               | Coincideix                                | Captura                           |
| ------------------- | ----------------------------------------- | --------------------------------- |
| `^(.+)`             | Nom de mostra des de l'inici              | Grup 1: nom de mostra             |
| `_S(\d+)`           | Número de mostra `_S1`, `_S2`, etc.       | Grup 2: número de mostra          |
| `_L(\d{3})`         | Número de carril `_L001`                  | Grup 3: carril (3 dígits)         |
| `_(R[12])`          | Direcció de lectura `_R1` o `_R2`         | Grup 4: direcció de lectura       |
| `_(\d{3})`          | Número de fragment `_001`                 | Grup 5: fragment (3 dígits)       |
| `\.fastq(?:\.gz)?$` | Extensió de fitxer `.fastq` o `.fastq.gz` | No capturat (?: és no capturador) |

Això analitza convencions de nomenclatura d'estil Illumina per extreure metadades automàticament.

Executeu el workflow modificat:

```bash title="Provar coincidència de patrons"
nextflow run main.nf
```

??? success "Sortida de la comanda"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [clever_pauling] DSL2 - revision: 605d2058b4

    [[id:sample_001, organism:human, tissue:liver, depth:30000000, quality:38.5, sample_num:1, lane:001, read:R1, chunk:001, priority:normal], /workspaces/training/side-quests/essential_scripting_patterns/data/sequences/SAMPLE_001_S1_L001_R1_001.fastq]
    [[id:sample_002, organism:mouse, tissue:brain, depth:25000000, quality:35.2, sample_num:2, lane:001, read:R1, chunk:001, priority:normal], /workspaces/training/side-quests/essential_scripting_patterns/data/sequences/SAMPLE_002_S2_L001_R1_001.fastq]
    [[id:sample_003, organism:human, tissue:kidney, depth:45000000, quality:42.1, sample_num:3, lane:001, read:R1, chunk:001, priority:high], /workspaces/training/side-quests/essential_scripting_patterns/data/sequences/SAMPLE_003_S3_L001_R1_001.fastq]
    ```

Això mostra les metadades enriquides dels noms de fitxers.

### 2.2. Generació Dinàmica de Scripts en Processos

Els blocs script dels processos són essencialment cadenes multilínia que es passen al shell. Podeu utilitzar **lògica condicional** (if/else, operadors ternaris) per generar dinàmicament diferents cadenes de script basant-se en característiques d'entrada. Això és essencial per gestionar diversos tipus d'entrada—com lectures single-end vs paired-end—sense duplicar definicions de processos.

Afegim un procés al nostre workflow que demostra aquest patró. Obriu `modules/fastp.nf` i doneu-hi un cop d'ull:

```groovy title="modules/fastp.nf" linenums="1"
process FASTP {
    container 'community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*_trimmed*.fastq.gz"), emit: reads

    script:
    """
    fastp \\
        --in1 ${reads[0]} \\
        --in2 ${reads[1]} \\
        --out1 ${meta.id}_trimmed_R1.fastq.gz \\
        --out2 ${meta.id}_trimmed_R2.fastq.gz \\
        --json ${meta.id}.fastp.json \\
        --html ${meta.id}.fastp.html \\
        --thread $task.cpus
    """
}
```

El procés pren fitxers FASTQ com a entrada i executa l'eina `fastp` per retallar adaptadors i filtrar lectures de baixa qualitat. Malauradament, la persona que va escriure aquest procés no va permetre les lectures single-end que tenim al nostre conjunt de dades d'exemple. Afegim-lo al nostre workflow i vegem què passa:

Primer, incloeu el mòdul a la primera línia del vostre workflow `main.nf`:

```groovy title="main.nf" linenums="1"
include { FASTP } from './modules/fastp.nf'
```

Després modifiqueu el bloc `workflow` per connectar el canal `ch_samples` al procés `FASTP`:

=== "Després"

    ```groovy title="main.nf" linenums="25" hl_lines="27"
    workflow {

        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row ->
                def sample_meta = [
                    id: row.sample_id.toLowerCase(),
                    organism: row.organism,
                    tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                    depth: row.sequencing_depth.toInteger(),
                    quality: row.quality_score.toDouble()
                ]
                def fastq_path = file(row.file_path)

                def m = (fastq_path.name =~ /^(.+)_S(\d+)_L(\d{3})_(R[12])_(\d{3})\.fastq(?:\.gz)?$/)
                def file_meta = m ? [
                    sample_num: m[0][2].toInteger(),
                    lane: m[0][3],
                    read: m[0][4],
                    chunk: m[0][5]
                ] : [:]

                def priority = sample_meta.quality > 40 ? 'high' : 'normal'
                return tuple(sample_meta + file_meta + [priority: priority], fastq_path)
            }

        ch_fastp = FASTP(ch_samples)
    }
    ```

=== "Abans"

    ```groovy title="main.nf" linenums="25" hl_lines="26"
    workflow {

        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row ->
                def sample_meta = [
                    id: row.sample_id.toLowerCase(),
                    organism: row.organism,
                    tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                    depth: row.sequencing_depth.toInteger(),
                    quality: row.quality_score.toDouble()
                ]
                def fastq_path = file(row.file_path)

                def m = (fastq_path.name =~ /^(.+)_S(\d+)_L(\d{3})_(R[12])_(\d{3})\.fastq(?:\.gz)?$/)
                def file_meta = m ? [
                    sample_num: m[0][2].toInteger(),
                    lane: m[0][3],
                    read: m[0][4],
                    chunk: m[0][5]
                ] : [:]

                def priority = sample_meta.quality > 40 ? 'high' : 'normal'
                return [sample_meta + file_meta + [priority: priority], file(row.file_path)]
            }
            .view()
    }
    ```

Executeu aquest workflow modificat:

```bash
nextflow run main.nf
```

??? failure "Sortida de la comanda"

    ```console
    ERROR ~ Error executing process > 'FASTP (3)'

    Caused by:
      Process `FASTP (3)` terminated with an error exit status (255)


    Command executed:

      fastp \
          --in1 SAMPLE_003_S3_L001_R1_001.fastq \
          --in2 null \
          --out1 sample_003_trimmed_R1.fastq.gz \
          --out2 sample_003_trimmed_R2.fastq.gz \
          --json sample_003.fastp.json \
          --html sample_003.fastp.html \
          --thread 2

    Command exit status:
      255

    Command output:
      (empty)
    ```

Podeu veure que el procés està intentant executar `fastp` amb un valor `null` per al segon fitxer d'entrada, cosa que està causant que falli. Això és perquè el nostre conjunt de dades conté lectures single-end, però el procés està codificat per esperar lectures paired-end (dos fitxers d'entrada a la vegada).

Arregleu això afegint lògica condicional al bloc `script:` del procés `FASTP`. Una sentència if/else comprova el recompte de fitxers de lectura i ajusta la comanda en conseqüència.

=== "Després"

    ```groovy title="main.nf" linenums="10" hl_lines="3-27"
        script:
        // Detecció simple de single-end vs paired-end
        def is_single = reads instanceof List ? reads.size() == 1 : true

        if (is_single) {
            def input_file = reads instanceof List ? reads[0] : reads
            """
            fastp \\
                --in1 ${input_file} \\
                --out1 ${meta.id}_trimmed.fastq.gz \\
                --json ${meta.id}.fastp.json \\
                --html ${meta.id}.fastp.html \\
                --thread $task.cpus
            """
        } else {
            """
            fastp \\
                --in1 ${reads[0]} \\
                --in2 ${reads[1]} \\
                --out1 ${meta.id}_trimmed_R1.fastq.gz \\
                --out2 ${meta.id}_trimmed_R2.fastq.gz \\
                --json ${meta.id}.fastp.json \\
                --html ${meta.id}.fastp.html \\
                --thread $task.cpus
            """
        }
    ```

=== "Abans"

    ```groovy title="main.nf" linenums="10" hl_lines="2-11"
            script:
            """
            fastp \\
                --in1 ${reads[0]} \\
                --in2 ${reads[1]} \\
                --out1 ${meta.id}_trimmed_R1.fastq.gz \\
                --out2 ${meta.id}_trimmed_R2.fastq.gz \\
                --json ${meta.id}.fastp.json \\
                --html ${meta.id}.fastp.html \\
                --thread $task.cpus
            """
        }
    ```

Ara el workflow pot gestionar tant lectures single-end com paired-end amb gràcia. La lògica condicional comprova el nombre de fitxers d'entrada i construeix la comanda apropiada per a `fastp`. Vegem si funciona:

```bash
nextflow run main.nf
```

??? success "Sortida de la comanda"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [adoring_rosalind] DSL2 - revision: 04b1cd93e9

    executor >  local (3)
    [31/a8ad4d] process > FASTP (3) [100%] 3 of 3 ✔
    ```

Es veu bé! Si comprovem les comandes reals que es van executar (personalitzeu per al vostre hash de tasca):

```console title="Comprovar comandes executades"
cat work/31/a8ad4d95749e685a6d842d3007957f/.command.sh
```

Podem veure que Nextflow va triar correctament la comanda adequada per a lectures single-end:

```bash title=".command.sh"
#!/bin/bash -ue
fastp \
    --in1 SAMPLE_003_S3_L001_R1_001.fastq \
    --out1 sample_003_trimmed.fastq.gz \
    --json sample_003.fastp.json \
    --html sample_003.fastp.html \
    --thread 2
```

Un altre ús comú de la lògica dinàmica de scripts es pot veure a [el mòdul Nextflow for Science Genomics](../../nf4science/genomics/02_joint_calling). En aquest mòdul, el procés GATK que es crida pot prendre múltiples fitxers d'entrada, però cadascun ha de tenir el prefix `-V` per formar una línia de comandes correcta. El procés utilitza scripting per transformar una col·lecció de fitxers d'entrada (`all_gvcfs`) en els arguments de comanda correctes:

```groovy title="manipulació de línia de comandes per a GATK" linenums="1" hl_lines="2 5"
    script:
    def gvcfs_line = all_gvcfs.collect { gvcf -> "-V ${gvcf}" }.join(' ')
    """
    gatk GenomicsDBImport \
        ${gvcfs_line} \
        -L ${interval_list} \
        --genomicsdb-workspace-path ${cohort_name}_gdb
    """
```

Aquests patrons d'utilitzar scripting en blocs script de processos són extremadament potents i es poden aplicar en molts escenaris - des de gestionar tipus d'entrada variables fins a construir arguments de línia de comandes complexos a partir de col·leccions de fitxers, fent que els vostres processos siguin veritablement adaptables als requisits diversos de dades del món real.

### 2.3. Interpolació de Variables: Variables de Nextflow i Shell

Els scripts de processos barregen variables de Nextflow, variables de shell i substitucions de comandes, cadascuna amb sintaxi d'interpolació diferent. Utilitzar la sintaxi incorrecta causa errors. Explorem-les amb un procés que crea un informe de processament.

Doneu un cop d'ull al fitxer de mòdul `modules/generate_report.nf`:

```groovy title="modules/generate_report.nf" linenums="1"
process GENERATE_REPORT {

    publishDir 'results/reports', mode: 'copy'

    input:
    tuple val(meta), path(reads)

    output:
    path "${meta.id}_report.txt"

    script:
    """
    echo "Processant ${reads}" > ${meta.id}_report.txt
    echo "Mostra: ${meta.id}" >> ${meta.id}_report.txt
    """
}
```

Aquest procés escriu un informe simple amb l'ID de mostra i el nom del fitxer. Ara executem-lo per veure què passa quan necessitem barrejar diferents tipus de variables.

Incloeu el procés al vostre `main.nf` i afegiu-lo al workflow:

=== "Després"

    ```groovy title="main.nf" linenums="1" hl_lines="2 30"
    include { FASTP } from './modules/fastp.nf'
    include { GENERATE_REPORT } from './modules/generate_report.nf'

    workflow {
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row ->
                def sample_meta = [
                    id: row.sample_id.toLowerCase(),
                    organism: row.organism,
                    tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                    depth: row.sequencing_depth.toInteger(),
                    quality: row.quality_score.toDouble()
                ]
                def fastq_path = file(row.file_path)

                def m = (fastq_path.name =~ /^(.+)_S(\d+)_L(\d{3})_(R[12])_(\d{3})\.fastq(?:\.gz)?$/)
                def file_meta = m ? [
                    sample_num: m[0][2].toInteger(),
                    lane: m[0][3],
                    read: m[0][4],
                    chunk: m[0][5]
                ] : [:]

                def priority = sample_meta.quality > 40 ? 'high' : 'normal'
                return tuple(sample_meta + file_meta + [priority: priority], fastq_path)
            }

        ch_fastp = FASTP(ch_samples)
        GENERATE_REPORT(ch_samples)
    }
    ```

=== "Abans"

    ```groovy title="main.nf" linenums="1" hl_lines="1 10-29"
    include { FASTP } from './modules/fastp.nf'

    workflow {
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row ->
                def sample_meta = [
                    id: row.sample_id.toLowerCase(),
                    organism: row.organism,
                    tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                    depth: row.sequencing_depth.toInteger(),
                    quality: row.quality_score.toDouble()
                ]
                def fastq_path = file(row.file_path)

                def m = (fastq_path.name =~ /^(.+)_S(\d+)_L(\d{3})_(R[12])_(\d{3})\.fastq(?:\.gz)?$/)
                def file_meta = m ? [
                    sample_num: m[0][2].toInteger(),
                    lane: m[0][3],
                    read: m[0][4],
                    chunk: m[0][5]
                ] : [:]

                def priority = sample_meta.quality > 40 ? 'high' : 'normal'
                return tuple(sample_meta + file_meta + [priority: priority], fastq_path)
            }

        ch_fastp = FASTP(ch_samples)
    }
    ```

Ara executeu el workflow i comproveu els informes generats a `results/reports/`. Haurien de contenir informació bàsica sobre cada mostra.

Però què passa si volem afegir informació sobre quan i on es va produir el processament? Modifiquem el procés per utilitzar variables de **shell** i una mica de substitució de comandes per incloure l'usuari actual, el nom de l'amfitrió i la data a l'informe:

=== "Després"

    ```groovy title="modules/generate_report.nf" linenums="10" hl_lines="5-7"
        script:
        """
        echo "Processant ${reads}" > ${meta.id}_report.txt
        echo "Mostra: ${meta.id}" >> ${meta.id}_report.txt
        echo "Processat per: ${USER}" >> ${meta.id}_report.txt
        echo "Nom d'amfitrió: $(hostname)" >> ${meta.id}_report.txt
        echo "Data: $(date)" >> ${meta.id}_report.txt
        """
    ```

=== "Abans"

    ```groovy title="modules/generate_report.nf" linenums="10"
        script:
        """
        echo "Processant ${reads}" > ${meta.id}_report.txt
        echo "Mostra: ${meta.id}" >> ${meta.id}_report.txt
        """
    ```

Si executeu això, notareu un error - Nextflow intenta interpretar `${USER}` com una variable de Nextflow que no existeix.

??? failure "Sortida de la comanda"

    ```console
    Error modules/generate_report.nf:15:27: `USER` is not defined
    │  15 |     echo "Processat per: ${USER}" >> ${meta.id}_report.txt
    ╰     |                           ^^^^

    ERROR ~ Script compilation failed
    ```

Necessitem escapar-la perquè Bash pugui gestionar-la en lloc de Nextflow.

Arregleu això escapant les variables de shell i les substitucions de comandes amb una barra invertida (`\`):

=== "Després"

    ```groovy title="modules/generate_report.nf" linenums="10" hl_lines="5-7"
        script:
        """
        echo "Processant ${reads}" > ${meta.id}_report.txt
        echo "Mostra: ${meta.id}" >> ${meta.id}_report.txt
        echo "Processat per: \${USER}" >> ${meta.id}_report.txt
        echo "Nom d'amfitrió: \$(hostname)" >> ${meta.id}_report.txt
        echo "Data: \$(date)" >> ${meta.id}_report.txt
        """
    ```

=== "Abans"

    ```groovy title="modules/generate_report.nf" linenums="10"
        script:
        """
        echo "Processant ${reads}" > ${meta.id}_report.txt
        echo "Mostra: ${meta.id}" >> ${meta.id}_report.txt
        echo "Processat per: ${USER}" >> ${meta.id}_report.txt
        echo "Nom d'amfitrió: $(hostname)" >> ${meta.id}_report.txt
        echo "Data: $(date)" >> ${meta.id}_report.txt
        """
    ```

Ara funciona! La barra invertida (`\`) diu a Nextflow "no interpretis això, passa-ho a Bash."

### Conclusió

En aquesta secció, heu après tècniques de **processament de cadenes**:

- **Expressions regulars per a anàlisi de fitxers**: Utilitzar l'operador `=~` i patrons regex (`~/patró/`) per extreure metadades de convencions de nomenclatura de fitxers complexes
- **Generació dinàmica de scripts**: Utilitzar lògica condicional (if/else, operadors ternaris) per generar diferents cadenes de script basant-se en característiques d'entrada
- **Interpolació de variables**: Comprendre quan Nextflow interpreta cadenes vs quan ho fa el shell
  - `${var}` - Variables de Nextflow (interpolades per Nextflow en temps de compilació del workflow)
  - `\${var}` - Variables d'entorn del shell (escapades, passades a bash en temps d'execució)
  - `\$(cmd)` - Substitució de comandes del shell (escapada, executada per bash en temps d'execució)

Aquests patrons de processament i generació de cadenes són essencials per gestionar els diversos formats de fitxers i convencions de nomenclatura que trobareu en workflows bioinformàtics del món real.

---

## 3. Crear Funcions Reutilitzables

La lògica complexa de workflow en línia en operadors de canal o definicions de processos redueix la llegibilitat i mantenibilitat. Les **funcions** us permeten extreure aquesta lògica en components amb nom i reutilitzables.

La nostra operació map ha crescut llarga i complexa. Extraiem-la en una funció reutilitzable utilitzant la paraula clau `def`.

Per il·lustrar com es veu això amb el nostre workflow existent, feu la modificació següent, utilitzant `def` per definir una funció reutilitzable anomenada `separateMetadata`:

=== "Després"

    ```groovy title="main.nf" linenums="1" hl_lines="4-24 29"
    include { FASTP } from './modules/fastp.nf'
    include { GENERATE_REPORT } from './modules/generate_report.nf'

    def separateMetadata(row) {
        def sample_meta = [
            id: row.sample_id.toLowerCase(),
            organism: row.organism,
            tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
            depth: row.sequencing_depth.toInteger(),
            quality: row.quality_score.toDouble()
        ]
        def fastq_path = file(row.file_path)

        def m = (fastq_path.name =~ /^(.+)_S(\d+)_L(\d{3})_(R[12])_(\d{3})\.fastq(?:\.gz)?$/)
        def file_meta = m ? [
            sample_num: m[0][2].toInteger(),
            lane: m[0][3],
            read: m[0][4],
            chunk: m[0][5]
        ] : [:]

        def priority = sample_meta.quality > 40 ? 'high' : 'normal'
        return tuple(sample_meta + file_meta + [priority: priority], fastq_path)
    }

    workflow {
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map{ row -> separateMetadata(row) }

        ch_fastp = FASTP(ch_samples)
        GENERATE_REPORT(ch_samples)
    }
    ```

=== "Abans"

    ```groovy title="main.nf" linenums="1" hl_lines="7-27"
    include { FASTP } from './modules/fastp.nf'
    include { GENERATE_REPORT } from './modules/generate_report.nf'

    workflow {
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row ->
                def sample_meta = [
                    id: row.sample_id.toLowerCase(),
                    organism: row.organism,
                    tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
                    depth: row.sequencing_depth.toInteger(),
                    quality: row.quality_score.toDouble()
                ]
                def fastq_path = file(row.file_path)

                def m = (fastq_path.name =~ /^(.+)_S(\d+)_L(\d{3})_(R[12])_(\d{3})\.fastq(?:\.gz)?$/)
                def file_meta = m ? [
                    sample_num: m[0][2].toInteger(),
                    lane: m[0][3],
                    read: m[0][4],
                    chunk: m[0][5]
                ] : [:]

                def priority = sample_meta.quality > 40 ? 'high' : 'normal'
                return tuple(sample_meta + file_meta + [priority: priority], fastq_path)
            }

        ch_fastp = FASTP(ch_samples)
        GENERATE_REPORT(ch_samples)
    }
    ```

En extreure aquesta lògica en una funció, hem reduït la lògica real del workflow a alguna cosa molt més neta:

```groovy title="workflow mínim"
    ch_samples = channel.fromPath("./data/samples.csv")
        .splitCsv(header: true)
        .map{ row -> separateMetadata(row) }

    ch_fastp = FASTP(ch_samples)
    GENERATE_REPORT(ch_samples)
```

Això fa que la lògica del workflow sigui molt més fàcil de llegir i entendre d'una ullada. La funció `separateMetadata` encapsula tota la lògica complexa per analitzar i enriquir metadades, fent-la reutilitzable i testable.

Executeu el workflow per assegurar-vos que encara funciona:

```bash
nextflow run main.nf
```

??? success "Sortida de la comanda"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [admiring_panini] DSL2 - revision: 8cc832e32f

    executor >  local (6)
    [8c/2e3f91] process > FASTP (3)           [100%] 3 of 3 ✔
    [7a/1b4c92] process > GENERATE_REPORT (3) [100%] 3 of 3 ✔
    ```

La sortida hauria de mostrar ambdós processos completant-se amb èxit. El workflow ara és molt més net i fàcil de mantenir, amb tota la lògica complexa de processament de metadades encapsulada a la funció `separateMetadata`.

### Conclusió

En aquesta secció, heu après **creació de funcions**:

- **Definir funcions amb `def`**: La paraula clau per crear funcions amb nom (com `def` en Python o `function` en JavaScript)
- **Àmbit de funcions**: Les funcions definides a nivell de script són accessibles a tot el vostre workflow de Nextflow
- **Valors de retorn**: Les funcions retornen automàticament l'última expressió, o utilitzen `return` explícit
- **Codi més net**: Extreure lògica complexa en funcions és una pràctica fonamental d'enginyeria de programari en qualsevol llenguatge

A continuació, explorarem com utilitzar closures en directives de processos per a l'assignació dinàmica de recursos.

---

## 4. Directives de Recursos Dinàmiques amb Closures

Fins ara hem utilitzat scripting al bloc `script` dels processos. Però les **closures** (introduïdes a la Secció 1.1) també són increïblement útils en directives de processos, especialment per a l'assignació dinàmica de recursos. Afegim directives de recursos al nostre procés FASTP que s'adaptin basant-se en les característiques de la mostra.

### 4.1. Assignació de recursos específica per mostra

Actualment, el nostre procés FASTP utilitza recursos per defecte. Fem-lo més intel·ligent assignant més CPUs per a mostres d'alta profunditat. Editeu `modules/fastp.nf` per incloure una directiva `cpus` dinàmica i una directiva `memory` estàtica:

=== "Després"

    ```groovy title="modules/fastp.nf" linenums="1" hl_lines="4-5"
    process FASTP {
        container 'community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690'

        cpus { meta.depth > 40000000 ? 2 : 1 }
        memory 2.GB

        input:
        tuple val(meta), path(reads)
    ```

=== "Abans"

    ```groovy title="modules/fastp.nf" linenums="1"
    process FASTP {
        container 'community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690'

        input:
        tuple val(meta), path(reads)
    ```

La closure `{ meta.depth > 40000000 ? 2 : 1 }` utilitza l'**operador ternari** (cobert a la Secció 1.1) i s'avalua per a cada tasca, permetent l'assignació de recursos per mostra. Les mostres d'alta profunditat (>40M lectures) obtenen 2 CPUs, mentre que les altres obtenen 1 CPU.

!!! note "Accedir a Variables d'Entrada en Directives"

    La closure pot accedir a qualsevol variable d'entrada (com `meta` aquí) perquè Nextflow avalua aquestes closures en el context de cada execució de tasca.

Executeu el workflow de nou amb l'opció `-ansi-log false` per facilitar veure els hashs de tasca.

```bash
nextflow run main.nf -ansi-log false
```

??? success "Sortida de la comanda"

    ```console
    N E X T F L O W  ~  version 25.10.2
    Launching `main.nf` [fervent_albattani] DSL2 - revision: fa8f249759
    [bd/ff3d41] Submitted process > FASTP (2)
    [a4/a3aab2] Submitted process > FASTP (1)
    [48/6db0c9] Submitted process > FASTP (3)
    [ec/83439d] Submitted process > GENERATE_REPORT (3)
    [bd/15d7cc] Submitted process > GENERATE_REPORT (2)
    [42/699357] Submitted process > GENERATE_REPORT (1)
    ```

Podeu comprovar la comanda `docker` exacta que es va executar per veure l'assignació de CPU per a qualsevol tasca donada:

```console title="Comprovar comanda docker"
cat work/48/6db0c9e9d8aa65e4bb4936cd3bd59e/.command.run | grep "docker run"
```

Hauríeu de veure alguna cosa com:

```bash title="comanda docker"
    docker run -i --cpu-shares 4096 --memory 2048m -e "NXF_TASK_WORKDIR" -v /workspaces/training/side-quests/essential_scripting_patterns:/workspaces/training/side-quests/essential_scripting_patterns -w "$NXF_TASK_WORKDIR" --name $NXF_BOXID community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690 /bin/bash -ue /workspaces/training/side-quests/essential_scripting_patterns/work/48/6db0c9e9d8aa65e4bb4936cd3bd59e/.command.sh
```

En aquest exemple hem triat un exemple que va sol·licitar 2 CPUs (`--cpu-shares 2048`), perquè era una mostra d'alta profunditat, però hauríeu de veure diferents assignacions de CPU depenent de la profunditat de la mostra. Proveu això per a les altres tasques també.

### 4.2. Estratègies de reintent

Un altre patró potent és utilitzar `task.attempt` per a estratègies de reintent. Per mostrar per què això és útil, començarem reduint l'assignació de memòria a FASTP a menys del que necessita. Canvieu la directiva `memory` a `modules/fastp.nf` a `1.GB`:

=== "Després"

    ```groovy title="modules/fastp.nf" linenums="1" hl_lines="5"
    process FASTP {
        container 'community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690'

        cpus { meta.depth > 40000000 ? 4 : 2 }
        memory 1.GB

        input:
        tuple val(meta), path(reads)
    ```

=== "Abans"

    ```groovy title="modules/fastp.nf" linenums="1" hl_lines="5"
    process FASTP {
        container 'community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690'

        cpus { meta.depth > 40000000 ? 4 : 2 }
        memory 2.GB

        input:
        tuple val(meta), path(reads)
    ```

... i executeu el workflow de nou:

```bash
nextflow run main.nf
```

??? failure "Sortida de la comanda"

    ```console hl_lines="2 11"
    Command exit status:
      137

    Command output:
      (empty)

    Command error:
      Detecting adapter sequence for read1...
      No adapter detected for read1

      .command.sh: line 7:   101 Killed                  fastp --in1 SAMPLE_002_S2_L001_R1_001.fastq --out1 sample_002_trimmed.fastq.gz --json sample_002.fastp.json --html sample_002.fastp.html --thread 2
    ```

Això indica que el procés va ser eliminat per excedir els límits de memòria.

Aquest és un escenari molt comú en workflows del món real - de vegades simplement no sabeu quanta memòria necessitarà una tasca fins que l'executeu.

Per fer el nostre workflow més robust, podem implementar una estratègia de reintent que augmenta l'assignació de memòria en cada intent, una vegada més utilitzant una closure de Groovy. Modifiqueu la directiva `memory` per multiplicar la memòria base per `task.attempt`, i afegiu directives `errorStrategy 'retry'` i `maxRetries 2`:

=== "Després"

    ```groovy title="modules/fastp.nf" linenums="1" hl_lines="5-7"
    process FASTP {
        container 'community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690'

        cpus { meta.depth > 40000000 ? 4 : 2 }
        memory { 1.GB * task.attempt }
        errorStrategy 'retry'
        maxRetries 2

        input:
        tuple val(meta), path(reads)
    ```

=== "Abans"

    ```groovy title="modules/fastp.nf" linenums="1" hl_lines="5"
    process FASTP {
        container 'community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690'

        cpus { meta.depth > 40000000 ? 4 : 2 }
        memory 2.GB

        input:
        tuple val(meta), path(reads)
    ```

Ara si el procés falla per memòria insuficient, Nextflow reintentarà amb més memòria:

- Primer intent: 1 GB (task.attempt = 1)
- Segon intent: 2.GB (task.attempt = 2)

... i així successivament, fins al límit `maxRetries`.

### Conclusió

Les directives dinàmiques amb closures us permeten:

- Assignar recursos basant-se en característiques d'entrada
- Implementar estratègies de reintent automàtic amb recursos creixents
- Combinar múltiples factors (metadades, número d'intent, prioritats)
- Utilitzar lògica condicional per a càlculs de recursos complexos

Això fa que els vostres workflows siguin tant més eficients (no sobre-assignant) com més robustos (reintent automàtic amb més recursos).

---

## 5. Lògica Condicional i Control de Processos

Anteriorment, vam utilitzar `.map()` amb scripting per transformar dades de canal. Ara utilitzarem lògica condicional per controlar quins processos s'executen basant-se en dades—essencial per a workflows flexibles que s'adapten a diferents tipus de mostres.

Els [operadors de flux de dades](https://www.nextflow.io/docs/latest/reference/operator.html) de Nextflow prenen closures avaluades en temps d'execució, permetent que la lògica condicional impulsi decisions de workflow basades en el contingut del canal.

### 5.1. Enrutament amb `.branch()`

Per exemple, suposem que les nostres mostres de seqüenciació necessiten ser retallades amb FASTP només si són mostres humanes amb una cobertura per sobre d'un cert llindar. Les mostres de ratolí o les mostres de baixa cobertura haurien de ser executades amb Trimgalore en lloc d'això (aquest és un exemple artificial, però il·lustra el punt).

Hem proporcionat un procés Trimgalore simple a `modules/trimgalore.nf`, doneu-hi un cop d'ull si voleu, però els detalls no són importants per a aquest exercici. El punt clau és que volem enrutar mostres basant-nos en les seves metadades.

Incloeu el nou mòdul de `modules/trimgalore.nf`:

=== "Després"

    ```groovy title="main.nf" linenums="1" hl_lines="2"
    include { FASTP } from './modules/fastp.nf'
    include { TRIMGALORE } from './modules/trimgalore.nf'
    ```

=== "Abans"

    ```groovy title="main.nf" linenums="1"
    include { FASTP } from './modules/fastp.nf'
    ```

... i després modifiqueu el vostre workflow `main.nf` per ramificar mostres basant-se en les seves metadades i enrutar-les a través del procés de retallat apropiat, així:

=== "Després"

    ```groovy title="main.nf" linenums="28" hl_lines="5-12"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row -> separateMetadata(row) }

        trim_branches = ch_samples
            .branch { meta, reads ->
                fastp: meta.organism == 'human' && meta.depth >= 30000000
                trimgalore: true
            }

        ch_fastp = FASTP(trim_branches.fastp)
        ch_trimgalore = TRIMGALORE(trim_branches.trimgalore)
        GENERATE_REPORT(ch_samples)
    ```

=== "Abans"

    ```groovy title="main.nf" linenums="28" hl_lines="5"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row -> separateMetadata(row) }

        ch_fastp = FASTP(ch_samples)
        GENERATE_REPORT(ch_samples)
    ```

Executeu aquest workflow modificat:

```bash
nextflow run main.nf
```

??? success "Sortida de la comanda"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [adoring_galileo] DSL2 - revision: c9e83aaef1

    executor >  local (6)
    [1d/0747ac] process > FASTP (2)           [100%] 2 of 2 ✔
    [cc/c44caf] process > TRIMGALORE (1)      [100%] 1 of 1 ✔
    [34/bd5a9f] process > GENERATE_REPORT (1) [100%] 3 of 3 ✔
    ```

Aquí, hem utilitzat expressions condicionals petites però potents dins de l'operador `.branch{}` per enrutar mostres basant-se en les seves metadades. Les mostres humanes amb alta cobertura passen per `FASTP`, mentre que totes les altres mostres passen per `TRIMGALORE`.

### 5.2. Utilitzar `.filter()` amb Truthiness

Un altre patró potent per controlar l'execució del workflow és l'operador `.filter()`, que utilitza una closure per determinar quins elements haurien de continuar pel pipeline. Dins de la closure de filtre, escriureu **expressions booleanes** que decideixen quins elements passen.

Nextflow (com molts llenguatges dinàmics) té un concepte de **"truthiness"** que determina quins valors s'avaluen com a `true` o `false` en contextos booleans:

- **Truthy**: Valors no nuls, cadenes no buides, números no zero, col·leccions no buides
- **Falsy**: `null`, cadenes buides `""`, zero `0`, col·leccions buides `[]` o `[:]`, `false`

Això significa que `meta.id` sol (sense `!= null` explícit) comprova si l'ID existeix i no està buit. Utilitzem això per filtrar mostres que no compleixen els nostres requisits de qualitat.

Afegiu el següent abans de l'operació branch:

=== "Després"

    ```groovy title="main.nf" linenums="28" hl_lines="5-11"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row -> separateMetadata(row) }

        // Filtrar mostres invàlides o de baixa qualitat
        ch_valid_samples = ch_samples
            .filter { meta, reads ->
                meta.id && meta.organism && meta.depth >= 25000000
            }

        trim_branches = ch_valid_samples
            .branch { meta, reads ->
                fastp: meta.organism == 'human' && meta.depth >= 30000000
                trimgalore: true
            }
    ```

=== "Abans"

    ```groovy title="main.nf" linenums="28" hl_lines="5"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map { row -> separateMetadata(row) }

        trim_branches = ch_samples
            .branch { meta, reads ->
                fastp: meta.organism == 'human' && meta.depth >= 30000000
                trimgalore: true
            }
    ```

Executeu el workflow de nou:

```bash
nextflow run main.nf
```

??? success "Sortida de la comanda"

    ```console
    N E X T F L O W  ~  version 25.10.2
    Launching `main.nf` [lonely_williams] DSL2 - revision: d0b3f121ec
    [94/b48eac] Submitted process > FASTP (2)
    [2c/d2b28f] Submitted process > GENERATE_REPORT (2)
    [65/2e3be4] Submitted process > GENERATE_REPORT (1)
    [94/b48eac] NOTE: Process `FASTP (2)` terminated with an error exit status (137) -- Execution is retried (1)
    [3e/0d8664] Submitted process > TRIMGALORE (1)
    [6a/9137b0] Submitted process > FASTP (1)
    [6a/9137b0] NOTE: Process `FASTP (1)` terminated with an error exit status (137) -- Execution is retried (1)
    [83/577ac0] Submitted process > GENERATE_REPORT (3)
    [a2/5117de] Re-submitted process > FASTP (1)
    [1f/a1a4ca] Re-submitted process > FASTP (2)
    ```

Com que hem triat un filtre que exclou algunes mostres, es van executar menys tasques.

L'expressió de filtre `meta.id && meta.organism && meta.depth >= 25000000` combina truthiness amb comparacions explícites:

- `meta.id && meta.organism` comprova que ambdós camps existeixen i no estan buits (utilitzant truthiness)
- `meta.depth >= 25000000` assegura profunditat de seqüenciació suficient amb una comparació explícita

!!! note "Truthiness a la Pràctica"

    L'expressió `meta.id && meta.organism` és més concisa que escriure:
    ```groovy
    meta.id != null && meta.id != '' && meta.organism != null && meta.organism != ''
    ```

    Això fa que la lògica de filtratge sigui molt més neta i fàcil de llegir.

### Conclusió

En aquesta secció, heu après a utilitzar lògica condicional per controlar l'execució del workflow utilitzant les interfícies de closure d'operadors de Nextflow com `.branch{}` i `.filter{}`, aprofitant truthiness per escriure expressions condicionals concises.

El nostre pipeline ara enruta intel·ligentment mostres a través de processos apropiats, però els workflows de producció necessiten gestionar dades invàlides amb gràcia. Fem que el nostre workflow sigui robust contra valors absents o nuls.

---

## 6. Operadors de Navegació Segura i Elvis

La nostra funció `separateMetadata` actualment assumeix que tots els camps CSV estan presents i són vàlids. Però què passa amb dades incompletes? Descobrim-ho.

### 6.1. El Problema: Accedir a Propietats que No Existeixen

Suposem que volem afegir suport per a informació opcional d'execució de seqüenciació. En alguns laboratoris, les mostres poden tenir un camp addicional per a l'ID d'execució de seqüenciació o número de lot, però el nostre CSV actual no té aquesta columna. Intentem accedir-hi igualment.

Modifiqueu la funció `separateMetadata` per incloure un camp run_id:

=== "Després"

    ```groovy title="main.nf" linenums="5" hl_lines="9"
    def separateMetadata(row) {
        def sample_meta = [
            id: row.sample_id.toLowerCase(),
            organism: row.organism,
            tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
            depth: row.sequencing_depth.toInteger(),
            quality: row.quality_score.toDouble()
        ]
        def run_id = row.run_id.toUpperCase()
    ```

=== "Abans"

    ```groovy title="main.nf" linenums="5"
    def separateMetadata(row) {
        def sample_meta = [
            id: row.sample_id.toLowerCase(),
            organism: row.organism,
            tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
            depth: row.sequencing_depth.toInteger(),
            quality: row.quality_score.toDouble()
        ]
    ```

Ara executeu el workflow:

```bash
nextflow run main.nf
```

??? failure "Sortida de la comanda"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [trusting_torvalds] DSL2 - revision: b56fbfbce2

    ERROR ~ Cannot invoke method toUpperCase() on null object

    -- Check script 'main.nf' at line: 13 or see '.nextflow.log' file for more details
    ```

Això falla amb una NullPointerException.

El problema és que `row.run_id` retorna `null` perquè la columna `run_id` no existeix al nostre CSV. Quan intentem cridar `.toUpperCase()` sobre `null`, falla. Aquí és on l'operador de navegació segura salva el dia.

### 6.2. Operador de Navegació Segura (`?.`)

L'operador de navegació segura (`?.`) retorna `null` en lloc de llançar una excepció quan es crida sobre un valor `null`. Si l'objecte abans de `?.` és `null`, tota l'expressió s'avalua com a `null` sense executar el mètode.

Actualitzeu la funció per utilitzar navegació segura:

=== "Després"

    ```groovy title="main.nf" linenums="4" hl_lines="9"
    def separateMetadata(row) {
        def sample_meta = [
            id: row.sample_id.toLowerCase(),
            organism: row.organism,
            tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
            depth: row.sequencing_depth.toInteger(),
            quality: row.quality_score.toDouble()
        ]
        def run_id = row.run_id?.toUpperCase()
    ```

=== "Abans"

    ```groovy title="main.nf" linenums="4" hl_lines="9"
    def separateMetadata(row) {
        def sample_meta = [
            id: row.sample_id.toLowerCase(),
            organism: row.organism,
            tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
            depth: row.sequencing_depth.toInteger(),
            quality: row.quality_score.toDouble()
        ]
        def run_id = row.run_id.toUpperCase()
    ```

Executeu de nou:

```bash
nextflow run main.nf
```

No falla! El workflow ara gestiona el camp absent amb gràcia. Quan `row.run_id` és `null`, l'operador `?.` evita la crida `.toUpperCase()`, i `run_id` es converteix en `null` en lloc de causar una excepció.

### 6.3. Operador Elvis (`?:`) per a Valors per Defecte

L'operador Elvis (`?:`) proporciona valors per defecte quan el costat esquerre és "falsy" (com s'ha explicat anteriorment). Rep el nom d'Elvis Presley perquè `?:` sembla el seu famós cabell i ulls quan es veu de costat!

Ara que estem utilitzant navegació segura, `run_id` serà `null` per a mostres sense aquest camp. Utilitzem l'operador Elvis per proporcionar un valor per defecte i afegir-lo al nostre mapa `sample_meta`:

=== "Després"

    ```groovy title="main.nf" linenums="5" hl_lines="9-10"
    def separateMetadata(row) {
        def sample_meta = [
            id: row.sample_id.toLowerCase(),
            organism: row.organism,
            tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
            depth: row.sequencing_depth.toInteger(),
            quality: row.quality_score.toDouble()
        ]
        def run_id = row.run_id?.toUpperCase() ?: 'UNSPECIFIED'
        sample_meta.run = run_id
    ```

=== "Abans"

    ```groovy title="main.nf" linenums="5" hl_lines="9"
    def separateMetadata(row) {
        def sample_meta = [
            id: row.sample_id.toLowerCase(),
            organism: row.organism,
            tissue: row.tissue_type.replaceAll('_', ' ').toLowerCase(),
            depth: row.sequencing_depth.toInteger(),
            quality: row.quality_score.toDouble()
        ]
        def run_id = row.run_id?.toUpperCase()
    ```

També afegiu un operador `view()` al workflow per veure els resultats:

=== "Després"

    ```groovy title="main.nf" linenums="30" hl_lines="4"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map{ row -> separateMetadata(row) }
            .view()
    ```

=== "Abans"

    ```groovy title="main.nf" linenums="30"
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map{ row -> separateMetadata(row) }
    ```

i executeu el workflow:

```bash
nextflow run main.nf
```

??? success "Sortida de la comanda"

    ```console
    [[id:sample_001, organism:human, tissue:liver, depth:30000000, quality:38.5, run:UNSPECIFIED, sample_num:1, lane:001, read:R1, chunk:001, priority:normal], /workspaces/training/side-quests/essential_scripting_patterns/data/sequences/SAMPLE_001_S1_L001_R1_001.fastq]
    [[id:sample_002, organism:mouse, tissue:brain, depth:25000000, quality:35.2, run:UNSPECIFIED, sample_num:2, lane:001, read:R1, chunk:001, priority:normal], /workspaces/training/side-quests/essential_scripting_patterns/data/sequences/SAMPLE_002_S2_L001_R1_001.fastq]
    [[id:sample_003, organism:human, tissue:kidney, depth:45000000, quality:42.1, run:UNSPECIFIED, sample_num:3, lane:001, read:R1, chunk:001, priority:high], /workspaces/training/side-quests/essential_scripting_patterns/data/sequences/SAMPLE_003_S3_L001_R1_001.fastq]
    ```

Perfecte! Ara totes les mostres tenen un camp `run` amb el seu ID d'execució real (en majúscules) o el valor per defecte 'UNSPECIFIED'. La combinació de `?.` i `?:` proporciona tant seguretat (sense fallades) com valors per defecte sensats.

Traieu l'operador `.view()` ara que hem confirmat que funciona.

!!! tip "Combinar Navegació Segura i Elvis"

    El patró `value?.method() ?: 'default'` és comú en workflows de producció:

    - `value?.method()` - Crida el mètode de forma segura, retorna `null` si `value` és `null`
    - `?: 'default'` - Proporciona un valor de reserva si el resultat és `null`

    Aquest patró gestiona dades absents/incompletes amb gràcia.

Utilitzeu aquests operadors de forma consistent en funcions, closures d'operadors (`.map{}`, `.filter{}`), scripts de processos i fitxers de configuració. Eviten fallades quan es gestionen dades del món real.

### Conclusió

- **Navegació segura (`?.`)**: Evita fallades amb valors nuls - retorna null en lloc de llançar una excepció
- **Operador Elvis (`?:`)**: Proporciona valors per defecte - `value ?: 'default'`
- **Combinant**: `value?.method() ?: 'default'` és el patró comú

Aquests operadors fan que els workflows siguin resilients a dades incompletes - essencial per a treball del món real.

---

## 7. Validació amb `error()` i `log.warn`

De vegades necessiteu aturar el workflow immediatament si els paràmetres d'entrada són invàlids. A Nextflow, podeu utilitzar funcions integrades com `error()` i `log.warn`, així com construccions de programació estàndard com sentències `if` i lògica booleana, per implementar lògica de validació. Afegim validació al nostre workflow.

Creeu una funció de validació abans del vostre bloc workflow, crideu-la des del workflow, i canvieu la creació del canal per utilitzar un paràmetre per al camí del fitxer CSV. Si el paràmetre falta o el fitxer no existeix, crideu `error()` per aturar l'execució amb un missatge clar.

=== "Després"

    ```groovy title="main.nf" linenums="1" hl_lines="5-20 23-24"
    include { FASTP } from './modules/fastp.nf'
    include { TRIMGALORE } from './modules/trimgalore.nf'
    include { GENERATE_REPORT } from './modules/generate_report.nf'

    def validateInputs() {
        // Comprovar que el paràmetre d'entrada està proporcionat
        if (!params.input) {
            error("Camí del fitxer CSV d'entrada no proporcionat. Si us plau, especifiqueu --input <file.csv>")
        }

        // Comprovar que el fitxer CSV existeix
        if (!file(params.input).exists()) {
            error("Fitxer CSV d'entrada no trobat: ${params.input}")
        }
    }
    ...
    workflow {
        validateInputs()
        ch_samples = channel.fromPath(params.input)
    ```

=== "Abans"

    ```groovy title="main.nf" linenums="1"
    include { FASTP } from './modules/fastp.nf'
    include { TRIMGALORE } from './modules/trimgalore.nf'
    include { GENERATE_REPORT } from './modules/generate_report.nf'

    ...
    workflow {
        ch_samples = channel.fromPath("./data/samples.csv")
    ```

Ara proveu d'executar sense el fitxer CSV:

```bash
nextflow run main.nf
```

??? failure "Sortida de la comanda"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [confident_coulomb] DSL2 - revision: 07059399ed

    WARN: Access to undefined parameter `input` -- Initialise it to a default value eg. `params.input = some_value`
    Camí del fitxer CSV d'entrada no proporcionat. Si us plau, especifiqueu --input <file.csv>
    ```

El workflow s'atura immediatament amb un missatge d'error clar en lloc de fallar misteriosament més tard

Ara executeu-lo amb un fitxer inexistent:

```bash
nextflow run main.nf --input ./data/nonexistent.csv
```

??? failure "Sortida de la comanda"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [cranky_gates] DSL2 - revision: 26839ae3eb

    Fitxer CSV d'entrada no trobat: ./data/nonexistent.csv
    ```

Finalment, executeu-lo amb el fitxer correcte:

```bash
nextflow run main.nf --input ./data/samples.csv
```

També podeu afegir validació dins de la funció `separateMetadata`. Utilitzem el no fatal `log.warn` per emetre advertències per a mostres amb baixa profunditat de seqüenciació, però encara permetre que el workflow continuï:

=== "Després"

    ```groovy title="main.nf" linenums="1" hl_lines="3-6"
        def priority = sample_meta.quality > 40 ? 'high' : 'normal'

        // Validar que les dades tenen sentit
        if (sample_meta.depth < 30000000) {
            log.warn "Baixa profunditat de seqüenciació per a ${sample_meta.id}: ${sample_meta.depth}"
        }

        return tuple(sample_meta + file_meta + [priority: priority], fastq_path)
    }
    ```

=== "Abans"

    ```groovy title="main.nf" linenums="1"
        def priority = sample_meta.quality > 40 ? 'high' : 'normal'

        return tuple(sample_meta + file_meta + [priority: priority], fastq_path)
    }
    ```

Executeu el workflow de nou amb el CSV original:

```bash
nextflow run main.nf --input ./data/samples.csv
```

??? warning "Sortida de la comanda"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [awesome_goldwasser] DSL2 - revision: a31662a7c1

    executor >  local (5)
    [ce/df5eeb] process > FASTP (2)           [100%] 2 of 2 ✔
    [-        ] process > TRIMGALORE          -
    [d1/7d2b4b] process > GENERATE_REPORT (3) [100%] 3 of 3 ✔
    WARN: Baixa profunditat de seqüenciació per a sample_002: 25000000
    ```

Veiem una advertència sobre baixa profunditat de seqüenciació per a una de les mostres.

### Conclusió

- **`error()`**: Atura el workflow immediatament amb un missatge clar
- **`log.warn`**: Emet advertències sense aturar el workflow
- **Validació primerenca**: Comproveu les entrades abans de processar per fallar ràpidament amb errors útils
- **Funcions de validació**: Creeu lògica de validació reutilitzable que es pot cridar a l'inici del workflow

La validació adequada fa que els workflows siguin més robustos i fàcils d'utilitzar capturant problemes primerencament amb missatges d'error clars.

---

## 8. Gestors d'Esdeveniments de Workflow

Fins ara, hem estat escrivint codi als nostres scripts de workflow i definicions de processos. Però hi ha una altra característica important que hauríeu de conèixer: els gestors d'esdeveniments de workflow.

Els gestors d'esdeveniments són closures que s'executen en punts específics del cicle de vida del vostre workflow. Són perfectes per afegir registre, notificacions o operacions de neteja. Aquests gestors s'han de definir al vostre script de workflow juntament amb la vostra definició de workflow.

### 8.1. El Gestor `onComplete`

El gestor d'esdeveniments més utilitzat és `onComplete`, que s'executa quan el vostre workflow finalitza (tant si ha tingut èxit com si ha fallat). Afegim-ne un per resumir els resultats del nostre pipeline.

Afegiu el gestor d'esdeveniments al vostre fitxer `main.nf`, dins de la vostra definició de workflow:

=== "Després"

    ```groovy title="main.nf" linenums="66" hl_lines="5-16"
        ch_fastp = FASTP(trim_branches.fastp)
        ch_trimgalore = TRIMGALORE(trim_branches.trimgalore)
        GENERATE_REPORT(ch_samples)

        workflow.onComplete = {
            println ""
            println "Resum d'execució del pipeline:"
            println "=========================="
            println "Completat a: ${workflow.complete}"
            println "Durada      : ${workflow.duration}"
            println "Èxit        : ${workflow.success}"
            println "workDir     : ${workflow.workDir}"
            println "estat sortida : ${workflow.exitStatus}"
            println ""
        }
    }
    ```

=== "Abans"

    ```groovy title="main.nf" linenums="66" hl_lines="4"
        ch_fastp = FASTP(trim_branches.fastp)
        ch_trimgalore = TRIMGALORE(trim_branches.trimgalore)
        GENERATE_REPORT(ch_samples)
    }
    ```

Aquesta closure s'executa quan el workflow es completa. Dins, teniu accés a l'objecte `workflow` que proporciona propietats útils sobre l'execució.

Executeu el vostre workflow i veureu aquest resum aparèixer al final!

```bash
nextflow run main.nf --input ./data/samples.csv -ansi-log false
```

??? success "Sortida de la comanda"

    ```console
    N E X T F L O W  ~  version 25.10.2
    Launching `main.nf` [marvelous_boltzmann] DSL2 - revision: a31662a7c1
    WARN: Baixa profunditat de seqüenciació per a sample_002: 25000000
    [9b/d48e40] Submitted process > FASTP (2)
    [6a/73867a] Submitted process > GENERATE_REPORT (2)
    [79/ad0ac5] Submitted process > GENERATE_REPORT (1)
    [f3/bda6cb] Submitted process > FASTP (1)
    [34/d5b52f] Submitted process > GENERATE_REPORT (3)

    Resum d'execució del pipeline:
    ==========================
    Completat a: 2025-10-10T12:14:24.885384+01:00
    Durada      : 2.9s
    Èxit        : true
    workDir     : /workspaces/training/side-quests/essential_scripting_patterns/work
    estat sortida : 0
    ```

Fem-lo més útil afegint lògica condicional:

=== "Després"

    ```groovy title="main.nf" linenums="66" hl_lines="5-22"
        ch_fastp = FASTP(trim_branches.fastp)
        ch_trimgalore = TRIMGALORE(trim_branches.trimgalore)
        GENERATE_REPORT(ch_samples)

        workflow.onComplete = {
            println ""
            println "Resum d'execució del pipeline:"
            println "=========================="
            println "Completat a: ${workflow.complete}"
            println "Durada      : ${workflow.duration}"
            println "Èxit        : ${workflow.success}"
            println "workDir     : ${workflow.workDir}"
            println "estat sortida : ${workflow.exitStatus}"
            println ""

            if (workflow.success) {
                println "✅ Pipeline completat amb èxit!"
            } else {
                println "❌ Pipeline ha fallat!"
                println "Error: ${workflow.errorMessage}"
            }
        }
    }
    ```

=== "Abans"

    ```groovy title="main.nf" linenums="66" hl_lines="5-16"
        ch_fastp = FASTP(trim_branches.fastp)
        ch_trimgalore = TRIMGALORE(trim_branches.trimgalore)
        GENERATE_REPORT(ch_samples)

        workflow.onComplete = {
            println ""
            println "Resum d'execució del pipeline:"
            println "=========================="
            println "Completat a: ${workflow.complete}"
            println "Durada      : ${workflow.duration}"
            println "Èxit        : ${workflow.success}"
            println "workDir     : ${workflow.workDir}"
            println "estat sortida : ${workflow.exitStatus}"
            println ""
        }
    }
    ```

Ara obtenim un resum encara més informatiu, incloent un missatge d'èxit/fallada i el directori de sortida si s'especifica:

??? success "Sortida de la comanda"

    ```console
    N E X T F L O W  ~  version 25.10.2
    Launching `main.nf` [boring_linnaeus] DSL2 - revision: a31662a7c1
    WARN: Baixa profunditat de seqüenciació per a sample_002: 25000000
    [e5/242efc] Submitted process > FASTP (2)
    [3b/74047c] Submitted process > GENERATE_REPORT (3)
    [8a/7a57e6] Submitted process > GENERATE_REPORT (1)
    [a8/b1a31f] Submitted process > GENERATE_REPORT (2)
    [40/648429] Submitted process > FASTP (1)

    Resum d'execució del pipeline:
    ==========================
    Completat a: 2025-10-10T12:16:00.522569+01:00
    Durada      : 3.6s
    Èxit        : true
    workDir     : /workspaces/training/side-quests/essential_scripting_patterns/work
    estat sortida : 0

    ✅ Pipeline completat amb èxit!
    ```

També podeu escriure el resum a un fitxer utilitzant operacions de fitxers:

```groovy title="main.nf - Escriure resum a fitxer"
workflow {
    // ... el vostre codi de workflow ...

    workflow.onComplete = {
        def summary = """
        Resum d'Execució del Pipeline
        ===========================
        Completat: ${workflow.complete}
        Durada : ${workflow.duration}
        Èxit  : ${workflow.success}
        Comanda  : ${workflow.commandLine}
        """

        println summary

        // Escriure a un fitxer de registre
        def log_file = file("${workflow.launchDir}/pipeline_summary.txt")
        log_file.text = summary
    }
}
```

### 8.2. El Gestor `onError`

A més de `onComplete`, hi ha un altre gestor d'esdeveniments que podeu utilitzar: `onError`, que s'executa només si el workflow falla:

```groovy title="main.nf - gestor onError"
workflow {
    // ... el vostre codi de workflow ...

    workflow.onError = {
        println "="* 50
        println "L'execució del pipeline ha fallat!"
        println "Missatge d'error: ${workflow.errorMessage}"
        println "="* 50

        // Escriure registre d'error detallat
        def error_file = file("${workflow.launchDir}/error.log")
        error_file.text = """
        Informe d'Error del Workflow
        =====================
        Hora: ${new Date()}
        Error: ${workflow.errorMessage}
        Informe d'error: ${workflow.errorReport ?: 'No hi ha informe detallat disponible'}
        """

        println "Detalls de l'error escrits a: ${error_file}"
    }
}
```

Podeu utilitzar múltiples gestors junts al vostre script de workflow:

```groovy title="main.nf - Gestors combinats"
workflow {
    // ... el vostre codi de workflow ...

    workflow.onError = {
        println "Workflow ha fallat: ${workflow.errorMessage}"
    }

    workflow.onComplete = {
        def duration_mins = workflow.duration.toMinutes().round(2)
        def status = workflow.success ? "ÈXIT ✅" : "FALLAT ❌"

        println """
        Pipeline finalitzat: ${status}
        Durada: ${duration_mins} minuts
        """
    }
}
```

### Conclusió

En aquesta secció, heu après:

- **Closures de gestors d'esdeveniments**: Closures al vostre script de workflow que s'executen en diferents punts del cicle de vida
- **Gestor `onComplete`**: Per a resums d'execució i informes de resultats
- **Gestor `onError`**: Per a gestió d'errors i registre de fallades
- **Propietats de l'objecte workflow**: Accedir a `workflow.success`, `workflow.duration`, `workflow.errorMessage`, etc.

Els gestors d'esdeveniments mostren com podeu utilitzar tot el poder del llenguatge Nextflow dins dels vostres scripts de workflow per afegir capacitats sofisticades de registre i notificació.

---

## Resum

Felicitats, ho heu aconseguit!

Al llarg d'aquesta missió secundària, heu construït un pipeline complet de processament de mostres que ha evolucionat des de la gestió bàsica de metadades fins a un workflow sofisticat i preparat per a producció.
Cada secció s'ha basat en l'anterior, demostrant com les construccions de programació transformen workflows simples en sistemes potents de processament de dades, amb els següents beneficis:

- **Codi més clar**: Comprendre flux de dades vs scripting us ajuda a escriure workflows més organitzats
- **Gestió robusta**: La navegació segura i els operadors Elvis fan que els workflows siguin resilients a dades absents
- **Processament flexible**: La lògica condicional permet que els vostres workflows processin diferents tipus de mostres adequadament
- **Recursos adaptatius**: Les directives dinàmiques optimitzen l'ús de recursos basant-se en característiques d'entrada

Aquesta progressió reflecteix l'evolució del món real dels pipelines bioinformàtics, des de prototips de recerca que gestionen unes poques mostres fins a sistemes de producció que processen milers de mostres a través de laboratoris i institucions.
Cada repte que heu resolt i patró que heu après reflecteix problemes reals que els desenvolupadors afronten quan escalen workflows de Nextflow.

Aplicar aquests patrons al vostre propi treball us permetrà construir workflows robustos i preparats per a producció.

### Patrons clau

1.  **Flux de dades vs Scripting:** Heu après a distingir entre operacions de flux de dades (orquestració de canals) i scripting (codi que manipula dades), incloent les diferències crucials entre operacions sobre diferents tipus com `collect` en Channel vs List.

    - Flux de dades: orquestració de canals

    ```groovy
    channel.fromPath('*.fastq').splitCsv(header: true)
    ```

    - Scripting: processament de dades en col·leccions

    ```groovy
    sample_data.collect { it.toUpperCase() }
    ```

2.  **Processament Avançat de Cadenes**: Heu dominat expressions regulars per analitzar noms de fitxers, generació dinàmica de scripts en processos, i interpolació de variables (Nextflow vs Bash vs Shell).

    - Coincidència de patrons

    ```groovy
    filename =~ ~/^(\w+)_(\w+)_(\d+)\.fastq$/
    ```

    - Funció amb retorn condicional

    ```groovy
    def parseSample(filename) {
        def matcher = filename =~ pattern
        return matcher ? [valid: true, data: matcher[0]] : [valid: false]
    }
    ```

    - Col·lecció de fitxers a arguments de comanda (en bloc script de procés)

    ```groovy
    script:
    def file_args = input_files.collect { file -> "--input ${file}" }.join(' ')
    """
    analysis_tool ${file_args} --output results.txt
    """
    ```

3.  **Crear Funcions Reutilitzables**: Heu après a extreure lògica complexa en funcions amb nom que es poden cridar des d'operadors de canal, fent que els workflows siguin més llegibles i mantenibles.

    - Definir una funció amb nom

    ```groovy
    def separateMetadata(row) {
        def sample_meta = [ /* codi ocult per brevetat */ ]
        def fastq_path = file(row.file_path)
        def m = (fastq_path.name =~ /^(.+)_S(\d+)_L(\d{3})_(R[12])_(\d{3})\.fastq(?:\.gz)?$/)
        def file_meta = m ? [ /* codi ocult per brevetat */ ] : [:]
        def priority = sample_meta.quality > 40 ? 'high' : 'normal'

        return tuple(sample_meta + file_meta + [priority: priority], fastq_path)
    }
    ```

    - Cridar la funció amb nom en un workflow

    ```groovy
    workflow {
        ch_samples = channel.fromPath("./data/samples.csv")
            .splitCsv(header: true)
            .map{ row -> separateMetadata(row) }

        ch_fastp = FASTP(ch_samples)
    }
    ```

4.  **Directives de Recursos Dinàmiques amb Closures**: Heu explorat l'ús de closures en directives de processos per a l'assignació adaptativa de recursos basant-se en característiques d'entrada.

    - Closures amb nom i composició

    ```groovy
    def enrichData = normalizeId >> addQualityCategory >> addFlags
    def processor = generalFunction.curry(fixedParam)
    ```

    - Closures amb accés a l'àmbit

    ```groovy
    def collectStats = { data -> stats.count++; return data }
    ```

5.  **Lògica Condicional i Control de Processos**: Heu afegit enrutament intel·ligent utilitzant operadors `.branch()` i `.filter()`, aprofitant truthiness per a expressions condicionals concises.

    - Utilitzar `.branch()` per enrutar dades a través de diferents branques de workflow

    ```groovy
    trim_branches = ch_samples
    .branch { meta, reads ->
        fastp: meta.organism == 'human' && meta.depth >= 30000000
        trimgalore: true
    }

    ch_fastp = FASTP(trim_branches.fastp)
    ch_trimgalore = TRIMGALORE(trim_branches.trimgalore)
    ```

    - Avaluació booleana amb Groovy Truth

    ```groovy
    if (sample.files) println "Té fitxers"
    ```

    - Utilitzar `filter()` per subconjuntar dades amb 'truthiness'

    ```groovy
    ch_valid_samples = ch_samples
        .filter { meta, reads ->
            meta.id && meta.organism && meta.depth >= 25000000
        }
    ```

6.  **Operadors de Navegació Segura i Elvis**: Heu fet el pipeline robust contra dades absents utilitzant `?.` per a accés segur a propietats contra nuls i `?:` per proporcionar valors per defecte.

    ```groovy
    def id = data?.sample?.id ?: 'unknown'
    ```

7.  **Validació amb error() i log.warn**: Heu après a validar entrades primerencament i fallar ràpidament amb missatges d'error clars.

    ```groovy
    try {
        def errors = validateSample(sample)
        if (errors) throw new RuntimeException("Invàlid: ${errors.join(', ')}")
    } catch (Exception e) {
        println "Error: ${e.message}"
    }
    ```

8.  **Gestors d'Esdeveniments de Configuració**: Heu après a utilitzar gestors d'esdeveniments de workflow (`onComplete` i `onError`) per a registre, notificacions i gestió del cicle de vida.

    - Utilitzar `onComplete` per registrar i notificar

    ```groovy
    workflow.onComplete = {
        println "Èxit        : ${workflow.success}"
        println "estat sortida : ${workflow.exitStatus}"

        if (workflow.success) {
            println "✅ Pipeline completat amb èxit!"
        } else {
            println "❌ Pipeline ha fallat!"
            println "Error: ${workflow.errorMessage}"
        }
    }
    ```

    - Utilitzar `onError` per prendre accions específicament en cas de fallada

    ```groovy
    workflow.onError = {
        // Escriure registre d'error detallat
        def error_file = file("${workflow.launchDir}/error.log")
        error_file.text = """
        Hora: ${new Date()}
        Error: ${workflow.errorMessage}
        Informe d'error: ${workflow.errorReport ?: 'No hi ha informe detallat disponible'}
        """

        println "Detalls de l'error escrits a: ${error_file}"
    }
    ```

### Recursos addicionals

- [Referència del Llenguatge Nextflow](https://nextflow.io/docs/latest/reference/syntax.html)
- [Operadors de Nextflow](https://www.nextflow.io/docs/latest/operator.html)
- [Sintaxi de Scripts de Nextflow](https://www.nextflow.io/docs/latest/script.html)
- [Biblioteca Estàndard de Nextflow](https://nextflow.io/docs/latest/reference/stdlib.html)

Assegureu-vos de consultar aquests recursos quan necessiteu explorar característiques més avançades.

Us beneficiareu de practicar i expandir les vostres habilitats per:

- Escriure workflows més nets amb separació adequada entre flux de dades i scripting
- Dominar la interpolació de variables per evitar errors comuns amb variables de Nextflow, Bash i shell
- Utilitzar directives de recursos dinàmiques per a workflows eficients i adaptatius
- Transformar col·leccions de fitxers en arguments de línia de comandes correctament formatats
- Gestionar diferents convencions de nomenclatura de fitxers i formats d'entrada amb gràcia utilitzant regex i processament de cadenes
- Construir codi reutilitzable i mantenible utilitzant patrons avançats de closures i programació funcional
- Processar i organitzar conjunts de dades complexos utilitzant operacions de col·leccions
- Afegir validació, gestió d'errors i registre per fer que els vostres workflows estiguin preparats per a producció
- Implementar gestió del cicle de vida del workflow amb gestors d'esdeveniments

---

## Què segueix?

Torneu al [menú de Missions Secundàries](./index.md) o feu clic al botó a la part inferiordreta de la pàgina per passar al següent tema de la llista.
