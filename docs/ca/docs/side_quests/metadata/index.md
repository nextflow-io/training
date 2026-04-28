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
- Entendre per què la interfície "meta map + fitxer de dades" és una convenció àmpliament utilitzada
- Afegir nous camps de metadades durant l'execució del workflow
- Utilitzar metadades per personalitzar el comportament dels processos i organitzar les sortides

Aquestes habilitats us ajudaran a construir pipelines més robustos i flexibles que puguin gestionar relacions complexes entre fitxers i requisits de processament.

### Prerequisits

Abans d'emprendre aquesta missió secundària, hauríeu de:

- Haver completat el tutorial [Hello Nextflow](../../hello_nextflow/index.md) o un curs equivalent per a principiants.
- Estar còmodes amb els conceptes i mecanismes bàsics de Nextflow (processos, canals, operadors)

---

## 0. Primers passos

#### Obriu el codespace de formació

Si encara no ho heu fet, assegureu-vos d'obrir l'entorn de formació tal com es descriu a la [Configuració de l'entorn](../../envsetup/index.md).

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

L'editor s'obre amb el directori del projecte en focus.

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

Utilitzarem una eina anomenada [`COWPY`](https://github.com/jeffbuttars/cowpy) per generar art ASCII de cada personatge dient la seva salutació gravada.

??? info "Què fa `COWPY`?"

    `COWPY` és una eina de línia de comandes que genera art ASCII per mostrar entrades de text arbitràries d'una manera divertida.
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

A més, utilitzarem una eina d'anàlisi de llenguatge anomenada `langid` per identificar en quin idioma parla cada personatge i organitzar les sortides del pipeline en conseqüència.

#### Reviseu l'assignació

El vostre repte és escriure un workflow de Nextflow que:

1. **Generi art ASCII** de cada personatge
2. **Organitzi** les sortides per família lingüística (llengües germàniques vs. romàniques)

Això representa un patró típic de workflow on les metadades específiques de cada fitxer guien les decisions de processament; exactament el tipus de problema que els mapes de metadades resolen de manera elegant.

#### Llista de comprovació de preparació

Creieu que esteu a punt per submergir-vos?

- [ ] Entenc l'objectiu d'aquest curs i els seus prerequisits
- [ ] El meu codespace està en funcionament
- [ ] He configurat el meu directori de treball adequadament
- [ ] Entenc l'assignació

Si podeu marcar totes les caselles, esteu a punt per començar.

---

## 1. Opcions bàsiques per carregar i utilitzar metadades

Obriu el fitxer de workflow `main.nf` per examinar l'esquelet del workflow que us donem com a punt de partida.

```groovy title="main.nf" linenums="1"
#!/usr/bin/env nextflow

workflow  {
    main:
    ch_datasheet = channel.fromPath("./data/datasheet.csv")
        .splitCsv(header: true)
        .view()

    publish:
    cowpy_art = channel.empty()
}

output {
    cowpy_art {
    }
}
```

L'operador [`splitCsv`](https://www.nextflow.io/docs/latest/reference/operator.html#splitcsv) llegeix cada fila del fitxer com un element del canal.
Aquest és el mateix enfocament que utilitzem per carregar dades CSV a Hello Nextflow, el nostre curs per a principiants.
Consulteu [aquesta secció](../../hello_nextflow/02_hello_channels.md#4-read-input-values-from-a-csv-file) si necessiteu un recordatori de com funciona.

Amb `header: true`, la primera fila es tracta com a capçaleres de columna, de manera que cada element es converteix en un mapa de parells clau-valor indexats pel nom de columna.

Tingueu en compte que com que encara no estem executant cap procés sobre les dades, els blocs `publish` i `output` són simplement esquelets.

### 1.1. Executeu el workflow

Executeu el workflow per veure com s'estructura el contingut del canal un cop tot s'ha carregat:

```bash
nextflow run main.nf
```

??? success "Sortida de la comanda"

    ```console
     N E X T F L O W   ~  version 25.10.4

    Launching `main.nf` [exotic_albattani] DSL2 - revision: c0d03cec83

    [id:sampleA, character:squirrel, recording:/workspaces/training/side-quests/metadata/data/bonjour.txt]
    [id:sampleB, character:tux, recording:/workspaces/training/side-quests/metadata/data/guten_tag.txt]
    [id:sampleC, character:sheep, recording:/workspaces/training/side-quests/metadata/data/hallo.txt]
    [id:sampleD, character:turkey, recording:/workspaces/training/side-quests/metadata/data/hello.txt]
    [id:sampleE, character:stegosaurus, recording:/workspaces/training/side-quests/metadata/data/hola.txt]
    [id:sampleF, character:moose, recording:/workspaces/training/side-quests/metadata/data/salut.txt]
    [id:sampleG, character:turtle, recording:/workspaces/training/side-quests/metadata/data/ciao.txt]
    ```

Com podeu veure, l'operador ha construït un mapa de parells clau-valor per a cada fila del fitxer CSV, amb les capçaleres de columna com a claus per als valors corresponents.

Cada entrada del mapa correspon a una columna del nostre full de dades:

- `id`
- `character`
- `recording`

Això facilita l'accés a camps específics de cada fila.
Per exemple, podríem accedir a l'identificador del fitxer amb `id` o a la ruta del fitxer txt amb `recording`.

??? info "(Opcional) Més informació sobre els mapes de Groovy"

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
     N E X T F L O W   ~  version 25.10.4

    Launching `map_demo.nf` [cheesy_plateau] DSL2 - revision: fae5b8496e

    map: [id:sampleA, character:squirrel]
    id: sampleA
    character: squirrel
    ```

### 1.2. Seleccionar un camp específic amb `map`

Utilitzarem l'operador `map` per iterar sobre cada element del canal i seleccionar només el camp `character`, al qual podem accedir pel nom utilitzant la notació de punt.

#### 1.2.1. Afegiu l'operació map

Per accedir a la columna `character`, afegiu l'operació `map` abans de l'operació `.view()` de la manera següent:

=== "Després"

    ```groovy title="main.nf" linenums="5" hl_lines="3-5"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                row.character
            }
            .view()
    ```

=== "Abans"

    ```groovy title="main.nf" linenums="5"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .view()
    ```

Aquesta manera d'accedir a un camp específic s'explica amb més detall a [aquesta secció](../../hello_nextflow/02_hello_channels.md#43-use-the-map-operator-to-extract-the-greetings) de Hello Nextflow, si necessiteu un recordatori de com funciona.

#### 1.2.2. Executeu el workflow

Executeu el workflow per verificar que podeu veure els noms dels personatges extrets.

```bash
nextflow run main.nf
```

??? success "Sortida de la comanda"

    ```console
     N E X T F L O W   ~  version 25.10.4

    Launching `main.nf` [exotic_albattani] DSL2 - revision: c0d03cec83

    squirrel
    tux
    sheep
    turkey
    stegosaurus
    moose
    turtle
    ```

Això mostra que podem accedir als valors de la columna `character` per a cada fila.

Ara fem alguna cosa amb aquestes dades: utilitzem els camps `character` i `recording` junts per generar art ASCII amb `COWPY`.

### 1.3. Emetre sub-canals amb `multiMap`

Us proporcionem un mòdul de procés preescrit per a `COWPY`, de manera que primer heu d'examinar els requisits d'entrada del procés.

Podeu obrir el fitxer per veure com és el procés:

```groovy title="modules/cowpy.nf" linenums="1" hl_lines="7 8"
// Genera art ASCII amb cowpy
process COWPY {

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

Com podeu veure, el procés pren dues entrades separades: un fitxer de gravació i un nom de personatge.
Importantment, tenim valors per a tots dos, però actualment estan agrupats dins de cada element del canal.

Una manera d'extreure múltiples camps en canals separats és l'operador [`multiMap`](https://www.nextflow.io/docs/latest/reference/operator.html#multimap), que divideix un canal en múltiples sub-canals amb nom en una sola operació.

#### 1.3.1. Afegiu l'operació multiMap

Substituïu l'operació `map` per `multiMap`:

=== "Després"

    ```groovy title="main.nf" linenums="5" hl_lines="3-6"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .multiMap { row ->
                file: row.recording
                character: row.character
            }
    ```

=== "Abans"

    ```groovy title="main.nf" linenums="5" hl_lines="3-6"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                row.character
            }
            .view()
    ```

El bloc `multiMap` defineix dos sub-canals amb nom (`file` i `character`) a partir de cada fila, als quals podem accedir com a `ch_datasheet.file` i `ch_datasheet.character`.

#### 1.3.2. Crideu COWPY sobre els sub-canals

Ara, incloeu el procés `COWPY` i passeu-li cada sub-canal com a argument separat:

=== "Després"

    ```groovy title="main.nf" linenums="1" hl_lines="3 14"
    #!/usr/bin/env nextflow

    include { COWPY } from './modules/cowpy.nf'

    workflow {
        main:
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .multiMap { row ->
                file: row.recording
                character: row.character
            }

        COWPY(ch_datasheet.file, ch_datasheet.character)

        publish:
        cowpy_art = channel.empty()
    }

    output {
        cowpy_art {
        }
    }
    ```

=== "Abans"

    ```groovy title="main.nf" linenums="1"
    #!/usr/bin/env nextflow

    workflow {
        main:
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .multiMap { row ->
                file: row.recording
                character: row.character
            }

        publish:
        cowpy_art = channel.empty()
    }

    output {
        cowpy_art {
        }
    }
    ```

Això ens permet passar els dos camps per separat tal com requereix `COWPY`.

#### 1.3.3. Configureu la publicació de la sortida

Finalment, afegiu la sortida de `COWPY` al bloc `publish:`:

=== "Després"

    ```groovy title="main.nf" linenums="14" hl_lines="4"
        COWPY(ch_datasheet.file, ch_datasheet.character)

        publish:
        cowpy_art = COWPY.out
    ```

=== "Abans"

    ```groovy title="main.nf" linenums="14" hl_lines="4"
        COWPY(ch_datasheet.file, ch_datasheet.character)

        publish:
        cowpy_art = channel.empty()
    ```

Això ens permetrà veure fàcilment les sortides produïdes pel workflow.

#### 1.3.4. Executeu el workflow

Executeu el workflow per comprovar que `COWPY` s'executa sobre les entrades que hem proporcionat:

```bash
nextflow run main.nf
```

??? success "Sortida de la comanda"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [clever_dijkstra] DSL2 - revision: a1b2c3d4e5

    executor >  local (7)
    [3a/f1c290] COWPY (7) [100%] 7 of 7 ✔
    ```

Com podeu veure, `COWPY` s'ha executat sobre cada fitxer utilitzant el personatge correcte per a cadascun.

??? abstract "Contingut del directori de resultats"

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

??? example "Contingut de results/cowpy-guten_tag.txt"

    ```console
    $ cat results/cowpy-guten_tag.txt
     _____________________________
    / Guten Tag, wie geht es dir? \
    \ Auf Wiedersehen, bis morgen /
     -----------------------------
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

Aquest enfocament funciona, però té una limitació: hem hagut de dividir el canal en dos sub-canals separats.
Si volguéssim passar més camps al procés, hauríem de dividir-los en més sub-canals.
Això podria resultar molest i desordenat.

Bona notícia: hi ha una manera més senzilla de fer-ho.

### 1.4. Agrupar-ho tot com una sola entrada al procés

En lloc de dividir els camps en canals separats, podem actualitzar el procés per rebre totes les entrades com una sola tupla, la qual cosa simplifica la crida al procés.

#### 1.4.1. Actualitzeu el procés COWPY

Actualitzeu `COWPY` per acceptar una tupla corresponent als tres elements de cada fila:

=== "Després"

    ```groovy title="modules/cowpy.nf" linenums="1" hl_lines="7 10 14"
    // Genera art ASCII amb cowpy
    process COWPY {

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'

        input:
        tuple val(id), val(character), path(recording)

        output:
        path "cowpy-${recording}"

        script:
        """
        cat ${recording} | cowpy -c ${character} > cowpy-${recording}
        """
    }
    ```

=== "Abans"

    ```groovy title="modules/cowpy.nf" linenums="1" hl_lines="7-8 11 15"
    // Genera art ASCII amb cowpy
    process COWPY {

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

Ara el procés pren una sola entrada que conté tot el que li volem passar.

#### 1.4.2. Utilitzeu `map()` per crear la tupla d'entrada

Encara hem d'utilitzar una operació de mapeig per enumerar els elements que volem passar a la tupla del procés:

=== "Després"

    ```groovy title="main.nf" linenums="5" hl_lines="3-5"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                [row.id, row.character, row.recording]
            }
    ```

=== "Abans"

    ```groovy title="main.nf" linenums="5" hl_lines="3-6"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .multiMap { row ->
                file: row.recording
                character: row.character
            }
    ```

Potser us pregunteu per què no podem simplement passar el mapa de Groovy sencer que prové de `splitCsv` tal com és.
És perquè hem de dir-li a Nextflow explícitament que el fitxer de gravació s'ha de gestionar com a ruta (és a dir, s'ha de preparar correctament).
Això passa a nivell de la interfície d'entrada de `COWPY`, on l'element `recording` es designa explícitament com a `path`.

#### 1.4.3. Actualitzeu la crida al procés

Finalment, substituïm les dues entrades separades de la crida al procés per la tupla única que acabem de crear:

=== "Després"

    ```groovy title="main.nf" linenums="5" hl_lines="7"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                [row.id, row.character, row.recording]
            }

        COWPY(ch_datasheet)
    ```

=== "Abans"

    ```groovy title="main.nf" linenums="5" hl_lines="7"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                [row.id, row.character, row.recording]
            }

        COWPY(ch_datasheet.file, ch_datasheet.character)
    ```

Això simplifica una mica la crida al procés.

#### 1.4.4. Executeu el workflow

Executeu el workflow per verificar que `COWPY` encara pot processar les dades correctament:

```bash
nextflow run main.nf
```

??? success "Sortida de la comanda"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [pedantic_lovelace] DSL2 - revision: b2c3d4e5f6

    executor >  local (7)
    [5e/2a1b34] COWPY (7) [100%] 7 of 7 ✔
    ```

La sortida és els mateixos set fitxers `cowpy-*.txt` que abans, ara produïts amb una crida més senzilla a `COWPY`.

??? abstract "Contingut del directori de resultats"

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

??? example "Contingut de results/cowpy-guten_tag.txt"

    ```console
    $ cat results/cowpy-guten_tag.txt
     _____________________________
    / Guten Tag, wie geht es dir? \
    \ Auf Wiedersehen, bis morgen /
     -----------------------------
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

Això és una lleugera millora respecte a l'enfocament amb `multiMap`.
Però encara hem hagut de desempaquetar el mapa de Groovy original per crear la tupla d'entrada, i hi ha un acoblament estret entre el procés i el full de dades: la definició d'entrada de `COWPY` ara fa referència directament als noms de columna `id`, `character` i `recording`.

```groovy
input:
tuple val(id), val(character), path(recording)
```

Si un col·laborador utilitza un full de dades amb una estructura diferent —amb columnes addicionals, o columnes en un ordre diferent— aquest procés no funcionarà sense modificacions.
Això fa que el procés sigui fràgil, perquè la seva estructura d'entrada està lligada a la composició exacta del full de dades.

Per resoldre-ho, necessitem una manera de passar totes les metadades com un paquet sense codificar la seva estructura exacta a la interfície del procés.

### 1.5. Utilitzar una interfície de meta map + fitxer

La solució és separar dues preocupacions diferents al canal: les **metadades sobre una mostra** i el **fitxer de dades** en si.
Agrupant totes les metadades en un sol mapa —el "meta map"— obtenim una tupla consistent de dos elements independentment de quantes columnes de metadades contingui el full de dades:

```groovy title="Syntax example"
input:
tuple val(meta), path(file)
```

Afegir o eliminar columnes del full de dades canvia el contingut de `meta`, però la forma de la tupla `[meta, file]` es manté constant.
Els processos que accepten aquesta estructura no necessiten saber ni preocupar-se de quants camps de metadades existeixen.

#### 1.5.1. Reorganitzeu el contingut de la tupla en un meta map

Reestructurem l'operació `map` per produir una tupla `[meta, file]`:

=== "Després"

    ```groovy title="main.nf" linenums="5" hl_lines="4 6 8 11"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                [[id: row.id, character: row.character], row.recording]
            }
            .view()

        // COWPY(ch_datasheet)  // S'actualitzarà al pas següent

        publish:
        cowpy_art = channel.empty() // COWPY.out
    ```

=== "Abans"

    ```groovy title="main.nf" linenums="5" hl_lines="4 7"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                [row.id, row.character, row.recording]
            }

        COWPY(ch_datasheet)

        publish:
        cowpy_art = COWPY.out
    ```

Notareu que també hem afegit una instrucció `view()`, hem comentat la crida a `COWPY` i hem substituït `COWPY.out` per `channel.empty()` perquè la definició d'entrada del procés encara no coincideix amb la nova estructura.

#### 1.5.2. Executeu el workflow per inspeccionar el contingut reorganitzat

Executeu el workflow per veure la nova forma del canal:

```bash
nextflow run main.nf
```

??? success "Sortida de la comanda"

    ```console title="View meta map"
     N E X T F L O W   ~  version 25.10.4

    Launching `main.nf` [lethal_booth] DSL2 - revision: 0d8f844c07

    [[id:sampleA, character:squirrel], /workspaces/training/side-quests/metadata/data/bonjour.txt]
    [[id:sampleB, character:tux], /workspaces/training/side-quests/metadata/data/guten_tag.txt]
    [[id:sampleC, character:sheep], /workspaces/training/side-quests/metadata/data/hallo.txt]
    [[id:sampleD, character:turkey], /workspaces/training/side-quests/metadata/data/hello.txt]
    [[id:sampleE, character:stegosaurus], /workspaces/training/side-quests/metadata/data/hola.txt]
    [[id:sampleF, character:moose], /workspaces/training/side-quests/metadata/data/salut.txt]
    [[id:sampleG, character:turtle], /workspaces/training/side-quests/metadata/data/ciao.txt]
    ```

Ara, cada element del canal és una tupla de dos elements: primer el meta map i segon el fitxer.

```console title="Example element structure"
[
  [id:sampleA, character:squirrel],
  /workspaces/training/side-quests/metadata/data/bonjour.txt
]
```

Si més endavant afegim una columna `language` al full de dades, estarà disponible com a `meta.language` sense necessitat de fer cap canvi a la definició d'entrada del procés.

#### 1.5.3. Actualitzeu el procés `COWPY` per utilitzar el meta map

Actualitzeu `COWPY` per acceptar l'estructura de tupla `[meta, file]`:

=== "Després"

    ```groovy title="modules/cowpy.nf" linenums="1" hl_lines="7 10 14"
    // Genera art ASCII amb cowpy
    process COWPY {

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'

        input:
        tuple val(meta), path(input_file)

        output:
        path "cowpy-${input_file}"

        script:
        """
        cat ${input_file} | cowpy -c ${meta.character} > cowpy-${input_file}
        """
    }
    ```

=== "Abans"

    ```groovy title="modules/cowpy.nf" linenums="1" hl_lines="7 10 14"
    // Genera art ASCII amb cowpy
    process COWPY {

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'

        input:
        tuple val(id), val(character), path(recording)

        output:
        path "cowpy-${recording}"

        script:
        """
        cat ${recording} | cowpy -c ${character} > cowpy-${recording}
        """
    }
    ```

Dins del bloc script, `meta.character` accedeix al camp `character` del meta map.
Qualsevol camp del meta map és accessible de la mateixa manera.

#### 1.5.4. Actualitzeu la crida al procés

Restaureu la crida a `COWPY` i connecteu la seva sortida per a la publicació:

=== "Després"

    ```groovy title="main.nf" linenums="5" hl_lines="7 10"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                [[id: row.id, character: row.character], row.recording]
            }

        COWPY(ch_datasheet)

        publish:
        cowpy_art = COWPY.out
    ```

=== "Abans"

    ```groovy title="main.nf" linenums="5" hl_lines="6 8 11"
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                [[id: row.id, character: row.character], row.recording]
            }
            .view()

        // COWPY(ch_datasheet)  // S'actualitzarà al pas següent

        publish:
        cowpy_art = channel.empty() // COWPY.out
    ```

També hem restaurat la publicació de la sortida.

#### 1.5.5. Executeu el workflow

Executeu el workflow per comprovar que tot funciona:

```bash
nextflow run main.nf
```

??? success "Sortida de la comanda"

    ```console
     N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [wise_sammet] DSL2 - revision: 99797b1e92

    executor >  local (7)
    [5d/dffd4e] COWPY (7) [100%] 7 of 7 ✔
    ```

El directori de resultats ara conté els fitxers d'art ASCII.

??? abstract "Contingut del directori"

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

??? example "Contingut de results/cowpy-guten_tag.txt"

    ```console
    $ cat results/cowpy-guten_tag.txt
     _____________________________
    / Guten Tag, wie geht es dir? \
    \ Auf Wiedersehen, bis morgen /
     -----------------------------
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

El procés ara rep totes les metadades com un paquet a través de `meta`, utilitza el que necessita (`meta.character`) i ignora la resta.

Aquesta és la interfície estàndard que utilitzen tots els mòduls de [nf-core](https://nf-co.re/).
El patró `tuple val(meta), path(file)` apareix de manera consistent a tota la biblioteca de mòduls nf-core, per la qual cosa els workflows que adopten aquesta convenció poden incorporar mòduls nf-core amb una fricció mínima.

### Conclusió

En aquesta secció, heu après:

- **Com llegir fulls de dades:** Utilitzant `splitCsv` per analitzar fitxers CSV amb informació de capçalera
- **Per què existeix la convenció del meta map:** Separar les metadades dels fitxers de dades en tuples `[meta, file]` manté l'estructura del canal estable a mesura que el full de dades evoluciona
- **Com utilitzar els camps del meta map dins d'un procés:** Qualsevol camp del meta map és accessible mitjançant la notació de punt al bloc script

---

## 2. Manipulacions addicionals de metadades

Ara que la interfície del meta map està en funcionament, podem enriquir-la a mesura que les dades flueixen pel pipeline.

Utilitzarem una eina anomenada [`langid`](https://github.com/saffsd/langid.py) per identificar l'idioma de cada fitxer de gravació.
Donat un fragment de text, produeix una predicció d'idioma i una puntuació de probabilitat a `stdout`.

### 2.1. Afegir un pas d'identificació de l'idioma

Us proporcionem un mòdul de procés preescrit anomenat `IDENTIFY_LANGUAGE` que encapsula l'eina `langid`.

Obriu el fitxer del mòdul per examinar el seu codi:

```groovy title="modules/langid.nf" linenums="1" hl_lines="7 10"
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

La definició d'entrada utilitza la mateixa estructura `tuple val(meta), path(file)` que acabem de construir a la secció 1, de manera que `ch_datasheet` pot alimentar directament aquest procés sense cap adaptació.

La sortida afegeix `stdout` com a tercer element: això captura la predicció d'idioma que `langid` imprimeix a la consola.
La comanda `sed` elimina la puntuació de probabilitat i el salt de línia final, deixant només el codi d'idioma de dues lletres.

#### 2.1.1. Afegiu una crida a `IDENTIFY_LANGUAGE`

Incloeu el mòdul de procés `IDENTIFY_LANGUAGE` i crideu-lo sobre el canal del full de dades:

=== "Després"

    ```groovy title="main.nf" linenums="1" hl_lines="4 14-16"
    #!/usr/bin/env nextflow

    include { COWPY } from './modules/cowpy.nf'
    include { IDENTIFY_LANGUAGE } from './modules/langid.nf'

    workflow {
        main:
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                [[id: row.id, character: row.character], row.recording]
            }

        // Executa langid per identificar l'idioma de cada salutació
        IDENTIFY_LANGUAGE(ch_datasheet)
        IDENTIFY_LANGUAGE.out.view()

        COWPY(ch_datasheet)
    ```

=== "Abans"

    ```groovy title="main.nf" linenums="1"
    #!/usr/bin/env nextflow

    include { COWPY } from './modules/cowpy.nf'

    workflow {
        main:
        ch_datasheet = channel.fromPath("./data/datasheet.csv")
            .splitCsv(header: true)
            .map { row ->
                [[id: row.id, character: row.character], row.recording]
            }

        COWPY(ch_datasheet)
    ```

La sortida principal d'aquest procés és simplement una cadena de text, de manera que no hi ha fitxers de sortida per publicar.
En canvi, utilitzem `IDENTIFY_LANGUAGE.out.view()` per veure els resultats de l'operació.

#### 2.1.2. Executeu el workflow

Executeu el workflow per produir la identificació de l'idioma, utilitzant `-resume` per evitar tornar a executar les tasques de `COWPY`:

```bash
nextflow run main.nf -resume
```

??? success "Sortida de la comanda"

    ```console
     N E X T F L O W   ~  version 25.10.4

    Launching `main.nf` [voluminous_mcnulty] DSL2 - revision: f9bcfebabb

    executor >  local (14)
    [5d/dffd4e] COWPY (7)             [100%] 7 of 7, cached: 7 ✔
    [4e/f722fe] IDENTIFY_LANGUAGE (7) [100%] 7 of 7 ✔
    [[id:sampleA, character:squirrel], /workspaces/training/side-quests/metadata/work/eb/f7148ebdd898fbe1136bec6a714acb/bonjour.txt, fr]
    [[id:sampleB, character:tux], /workspaces/training/side-quests/metadata/work/16/71d72410952c22cd0086d9bca03680/guten_tag.txt, de]
    [[id:sampleD, character:turkey], /workspaces/training/side-quests/metadata/work/c4/b7562adddc1cc0b7d414ec45d436eb/hello.txt, en]
    [[id:sampleC, character:sheep], /workspaces/training/side-quests/metadata/work/ea/04f5d979429e4455e14b9242fb3b45/hallo.txt, de]
    [[id:sampleF, character:moose], /workspaces/training/side-quests/metadata/work/5a/6c2b84bf8fadb98e28e216426be079/salut.txt, fr]
    [[id:sampleE, character:stegosaurus], /workspaces/training/side-quests/metadata/work/af/ee7c69bcab891c40d0529305f6b9e7/hola.txt, es]
    [[id:sampleG, character:turtle], /workspaces/training/side-quests/metadata/work/4e/f722fe47271ba7ebcd69afa42964ca/ciao.txt, it]
    ```

Ara tenim una predicció d'idioma per a cada fitxer del conjunt de dades.

Tingueu en compte que la tupla de sortida es compon de `[meta, file, lang_id]`, la qual cosa significa que el meta map i el fitxer es mantenen associats amb el nou resultat.

!!! note "Nota"

    Aquest patró de mantenir el meta map associat amb els resultats facilita l'associació de resultats entre canals més endavant.
    No podeu confiar en l'ordre dels elements dels canals per associar les dades correctament.
    En canvi, heu d'utilitzar claus.
    Els meta maps proporcionen una estructura ideal per a aquest propòsit.

    Aquest cas d'ús s'explora en detall a la missió secundària [Splitting & Grouping](../splitting_and_grouping/index.md).

### 2.2. Augmentar les metadades amb les sortides del procés

La predicció d'idioma és en si mateixa una metadada sobre les dades del fitxer.
En lloc de mantenir-la com un element separat, incorporem-la de nou al meta map.

#### 2.2.1. Creeu un meta map nou i ampliat

Podem crear un nou meta map per substituir l'original utilitzant l'operador `+` de Groovy:

=== "Després"

    ```groovy title="main.nf" linenums="14" hl_lines="3-7"
        // Executa langid per identificar l'idioma de cada salutació
        IDENTIFY_LANGUAGE(ch_datasheet)
        IDENTIFY_LANGUAGE.out
            .map { meta, file, lang_id ->
                [meta + [lang: lang_id], file]
            }
            .view()
    ```

=== "Abans"

    ```groovy title="main.nf" linenums="14" hl_lines="3"
        // Executa langid per identificar l'idioma de cada salutació
        IDENTIFY_LANGUAGE(ch_datasheet)
        IDENTIFY_LANGUAGE.out.view()
    ```

El nucli d'aquesta operació és `#!groovy meta + [lang: lang_id]`.

Aquest codi essencialment crea un mapa temporal amb un sol parell clau-valor que conté el codi d'idioma (`[lang: lang_id]`), i després utilitza l'operador `+` de Groovy per combinar-lo amb el mapa `meta` original que conté les metadades preexistents, produint un meta map nou i ampliat.

Per a una explicació més detallada, consulteu el quadre a continuació.

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
    new_map = map1 + [lang: lang_id]
    ```

    Aquí, `[lang: lang_id]` crea un nou mapa sense nom al vol, i `map1 + ` fusiona `map1` amb el nou mapa sense nom, produint el mateix contingut de `new_map` que abans.

    Elegant, oi?

    **Ara transposem-ho al context d'una operació `channel.map()` de Nextflow.**

    El codi es converteix en:

    ```groovy
    .map { map1, lang_id ->
        map1 + [lang: lang_id]
    }
    ```

    Això fa el següent:

    - `#!groovy map1, lang_id ->` pren els dos elements de la tupla
    - `#!groovy map1 + [lang: lang_id]` crea el nou mapa tal com s'ha detallat anteriorment

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

    Si us costa seguir per què el `file` sembla que es mou a l'operació `map`, imagineu que en lloc de `#!groovy [meta + [lang: lang_id], file]`, aquella línia llegís `[new_map, file]`.
    Això hauria de deixar més clar que simplement estem deixant el `file` al seu lloc original en la segona posició de la tupla. Simplement hem pres el valor `new_info` i l'hem incorporat al mapa que es troba en la primera posició.

    **I això ens porta de nou a l'estructura de canal `tuple val(meta), path(file)`!**

#### 2.2.2. Executeu el workflow

Un cop estigueu segurs que enteneu el que fa aquest codi, executeu el workflow per veure si ha funcionat:

```bash
nextflow run main.nf -resume
```

??? success "Sortida de la comanda"

    ```console
     N E X T F L O W   ~  version 25.10.4

    Launching `main.nf` [cheeky_fermat] DSL2 - revision: d096281ee4

    [5d/dffd4e] COWPY (7)             [100%] 7 of 7, cached: 7 ✔
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

!!! tip "Consell: Eliminar claus d'un meta map"

    Podeu eliminar una clau d'un meta map utilitzant el mètode [`subMap`](https://docs.groovy-lang.org/latest/html/groovy-jdk/java/util/Map.html#subMap(java.util.Collection)) de Groovy, que retorna un nou mapa que conté només les claus que especifiqueu:

    ```groovy
    meta.subMap(['id', 'character'])  // retorna un mapa amb només 'id' i 'character'
    ```

    Això és útil quan un procés o mòdul posterior no necessita tots els camps que s'han acumulat al meta map.

### 2.3. Assignar un grup lingüístic mitjançant condicionals

Amb la predicció d'idioma al meta map, podem derivar-ne més metadades.
Els idiomes del nostre conjunt de dades es divideixen en dues famílies: germàniques (anglès, alemany) i romàniques (francès, espanyol, italià).
Afegir un camp `lang_group` farà que aquesta classificació estigui disponible més endavant al pipeline.

#### 2.3.1. Afegiu una operació `map` amb la lògica condicional

Utilitzarem una segona operació `map` amb lògica condicional per assignar la família lingüística:

```groovy
.map { meta, file ->

    // la lògica condicional que defineix lang_group va aquí

    [meta + [lang_group: lang_group], file]
}
```

Aquí teniu la lògica a aplicar:

- Comenceu amb `lang_group = 'unknown'` com a valor per defecte.
- Si `meta.lang` és `'de'` o `'en'`, establiu `lang_group` a `'germanic'`.
- En cas contrari, si `meta.lang` es troba a `['fr', 'es', 'it']`, establiu `lang_group` a `'romance'`.

!!! tip "Consell"

    Podeu accedir al valor de `lang` dins de l'operació map amb `meta.lang`.

Feu els canvis següents al workflow:

=== "Després"

    ```groovy title="main.nf" linenums="14" hl_lines="7-19 21"
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

        ch_languages.view()
    ```

=== "Abans"

    ```groovy title="main.nf" linenums="14" hl_lines="7"
        // Executa langid per identificar l'idioma de cada salutació
        IDENTIFY_LANGUAGE(ch_datasheet)
        IDENTIFY_LANGUAGE.out
            .map { meta, file, lang_id ->
                [meta + [lang: lang_id], file]
            }
            .view()
    ```

Punts clau:

- `def lang_group = "unknown"` inicialitza la variable amb un valor per defecte segur.
- L'estructura `if / else if` gestiona les dues famílies lingüístiques; qualsevol altra cosa es manté com a `'unknown'`.
- `#!groovy .set { ch_languages }` dóna un nom al canal resultant per utilitzar-lo al pas següent.

<!-- TODO (futur) Afegir nota/enllaços a la documentació rellevant a la secció de recursos addicionals -->

#### 2.3.2. Executeu el workflow:

Executeu el workflow per verificar que funciona:

```bash
nextflow run main.nf -resume
```

??? success "Sortida de la comanda"

    ```console
     N E X T F L O W   ~  version 25.10.4

    Launching `main.nf` [wise_almeida] DSL2 - revision: 46778c3cd0

    [5d/dffd4e] COWPY (7)             [100%] 7 of 7, cached: 7 ✔
    [da/652cc6] IDENTIFY_LANGUAGE (7) [100%] 7 of 7, cached: 7 ✔
    [[id:sampleA, character:squirrel, lang:fr, lang_group:romance], /workspaces/training/side-quests/metadata/data/bonjour.txt]
    [[id:sampleB, character:tux, lang:de, lang_group:germanic], /workspaces/training/side-quests/metadata/data/guten_tag.txt]
    [[id:sampleC, character:sheep, lang:de, lang_group:germanic], /workspaces/training/side-quests/metadata/data/hallo.txt]
    [[id:sampleD, character:turkey, lang:en, lang_group:germanic], /workspaces/training/side-quests/metadata/data/hello.txt]
    [[id:sampleE, character:stegosaurus, lang:es, lang_group:romance], /workspaces/training/side-quests/metadata/data/hola.txt]
    [[id:sampleF, character:moose, lang:fr, lang_group:romance], /workspaces/training/side-quests/metadata/data/salut.txt]
    [[id:sampleG, character:turtle, lang:it, lang_group:romance], /workspaces/training/side-quests/metadata/data/ciao.txt]
    ```

El meta map ara conté quatre camps: `id`, `character`, `lang` i `lang_group`.
L'estructura del canal continua sent `[meta, file]`.

### 2.4. Utilitzar metadades per anomenar i organitzar les sortides

Amb `lang` i `lang_group` ara disponibles al meta map, podem utilitzar-los per afegir un codi d'idioma als noms dels fitxers de sortida i organitzar-los en subdirectoris per família lingüística.

Això requereix tres canvis: actualitzar el procés `COWPY` per canviar el nom de la seva sortida i incloure `meta` en el que emet, actualitzar la crida a `COWPY` per executar-se sobre `ch_languages`, i actualitzar el bloc de sortida per especificar la ruta del subdirectori.

#### 2.4.1. Actualitzeu el procés `COWPY`

Canvieu el nom del fitxer de sortida utilitzant el codi d'idioma del meta map, i afegiu `meta` a la sortida perquè el bloc de sortida pugui accedir a `lang_group` per a l'enrutament al subdirectori:

=== "Després"

    ```groovy title="modules/cowpy.nf" linenums="9" hl_lines="2 6"
        output:
        tuple val(meta), path("${meta.lang}-${input_file}")

        script:
        """
        cat ${input_file} | cowpy -c ${meta.character} > ${meta.lang}-${input_file}
        """
    ```

=== "Abans"

    ```groovy title="modules/cowpy.nf" linenums="9" hl_lines="2 6"
        output:
        path "cowpy-${input_file}"

        script:
        """
        cat ${input_file} | cowpy -c ${meta.character} > cowpy-${input_file}
        """
    ```

Això mostra com podem aprofitar altres camps de metadades per personalitzar el comportament d'un procés, sense haver de modificar la definició d'entrada en absolut.

#### 2.4.2. Actualitzeu la crida a `COWPY` per executar-se sobre `ch_languages`

Substituïu `COWPY(ch_datasheet)` per `COWPY(ch_languages)`:

=== "Després"

    ```groovy title="main.nf" linenums="32" hl_lines="3"
        .set { ch_languages }

        COWPY(ch_languages)

        publish:
        cowpy_art = COWPY.out
    }
    ```

=== "Abans"

    ```groovy title="main.nf" linenums="32" hl_lines="3 5"
        .set { ch_languages }

        ch_languages.view()

        COWPY(ch_datasheet)

        publish:
        cowpy_art = COWPY.out
    }
    ```

També eliminem la línia `ch_languages.view()` ja que no necessitem inspeccionar el contingut del canal.

#### 2.4.3. Actualitzeu el bloc de sortida

Afegiu una closure `path` al bloc `output {}` per enrutar cada fitxer al subdirectori del seu grup lingüístic:

=== "Després"

    ```groovy title="main.nf" linenums="40" hl_lines="3"
    output {
        cowpy_art {
            path { meta, file -> meta.lang_group }
        }
    }
    ```

=== "Abans"

    ```groovy title="main.nf" linenums="40" hl_lines="2 3"
    output {
        cowpy_art {
        }
    }
    ```

Això mostra com podem utilitzar metadades per organitzar les sortides amb gran flexibilitat.

#### 2.4.4. Executeu el pipeline complet

Esborreu els resultats anteriors i executeu el pipeline complet:

```bash
rm -r results
nextflow run main.nf
```

??? success "Sortida de la comanda"

    ```console
     N E X T F L O W   ~  version 25.10.4

    Launching `main.nf` [suspicious_crick] DSL2 - revision: 25541014c5

    executor >  local (14)
    [5d/dffd4e] IDENTIFY_LANGUAGE (7) [100%] 7 of 7 ✔
    [e7/317c18] COWPY (7)             [100%] 7 of 7 ✔
    ```

El directori de resultats ara està organitzat per família lingüística, amb cada fitxer anomenat segons l'idioma detectat:

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

La closure `path` al bloc `output {}` rep cada tupla `[meta, file]` i retorna `meta.lang_group` com a nom del subdirectori.
El nom del fitxer en si prové del que el procés produeix (`#!groovy "${meta.lang}-${input_file}"`).
Tots dos elements de metadades (codi d'idioma i grup lingüístic) provenen del meta map enriquit construït en aquesta secció.

### Conclusió

En aquesta secció, heu après:

- **Com augmentar el meta map amb les sortides del procés:** Afegir noves claus amb `#!groovy meta + [clau: valor]` manté l'estructura del canal `[meta, file]` intacta mentre s'enriqueixen les metadades.
- **Com derivar metadades a partir de metadades:** La lògica condicional dins d'una operació `map` pot calcular nous camps a partir dels existents.
- **Com utilitzar metadades per organitzar les sortides:** La closure `path` al bloc `output {}` pot llegir del meta map per enrutar fitxers a subdirectoris.

---

## 3. Consideracions de robustesa

Quan els valors de les metadades guien el comportament del procés, les dades que falten o estan incompletes poden causar problemes difícils de diagnosticar.
Aquí teniu el que cal esperar i com gestionar-ho.

### 3.1. Què passa quan falta un camp de metadades requerit

El valor `character` és necessari perquè el procés `COWPY` produeixi un resultat vàlid.
El mode de fallada depèn de si la columna existeix al full de dades però està buida, o si és absent del tot.

#### 3.1.1. La columna existeix però un valor està buit

Suposem que una entrada del full de dades té el camp `character` en blanc:

```csv title="datasheet.csv" linenums="1" hl_lines="2"
id,character,recording
sampleA,,/workspaces/training/side-quests/metadata/data/bonjour.txt
sampleB,tux,/workspaces/training/side-quests/metadata/data/guten_tag.txt
...
```

La clau `character` es crea per a totes les entrades quan s'analitza el full de dades, però `meta.character` per a `sampleA` serà una cadena buida.
Quan Nextflow substitueix `#!groovy ${meta.character}` a la comanda, l'eina `COWPY` rep un argument buit per a `-c` i falla:

??? failure "Sortida de la comanda"

    ```console hl_lines="8 11 16 28"
     N E X T F L O W   ~  version 25.10.4

    Launching `main.nf` [marvelous_hirsch] DSL2 - revision: 0dfeee3cc1

    executor >  local (9)
    [c1/c5dd4f] process > IDENTIFY_LANGUAGE (7) [ 85%] 6 of 7
    [d3/b7c415] process > COWPY (2)             [  0%] 0 of 6
    ERROR ~ Error executing process > 'COWPY (1)'

    Caused by:
      Process `COWPY (1)` terminated with an error exit status (2)


    Command executed:

      cat bonjour.txt | cowpy -c  > fr-bonjour.txt

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

El missatge d'error (`expected one argument`) apunta a l'indicador `-c` buit.
Comprovar el fitxer `.command.sh` al directori de treball confirma que la comanda s'ha executat amb un valor buit.

#### 3.1.2. La columna no existeix al full de dades

Si la columna `character` és absent del tot:

```csv title="datasheet.csv" linenums="1"
id,recording
sampleA,/workspaces/training/side-quests/metadata/data/bonjour.txt
sampleB,/workspaces/training/side-quests/metadata/data/guten_tag.txt
...
```

La clau `character` no es crea mai al meta map.
Quan el script del procés avalua `#!groovy ${meta.character}`, la clau absent retorna `null`, i Nextflow literalment substitueix la cadena `null` a la comanda:

??? failure "Sortida de la comanda"

    ```console hl_lines="8 11 16"
     N E X T F L O W   ~  version 25.10.4

    Launching `main.nf` [jovial_bohr] DSL2 - revision: eaaf375827

    executor >  local (9)
    [0d/ada9db] process > IDENTIFY_LANGUAGE (5) [ 85%] 6 of 7
    [06/28065f] process > COWPY (2)             [  0%] 0 of 6
    ERROR ~ Error executing process > 'COWPY (2)'

    Caused by:
      Process `COWPY (2)` terminated with an error exit status (1)


    Command executed:

      cat guten_tag.txt | cowpy -c null > de-guten_tag.txt

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

El `cowpy -c null` a la comanda executada és la pista de diagnòstic.

### 3.2. Estratègies per gestionar metadades que falten

Hi ha dos enfocaments complementaris per fer els workflows més robustos davant de metadades que falten.

**1. Validació d'entrada**

La solució més fiable és validar el full de dades abans que comenci qualsevol processament, de manera que els problemes es detectin aviat amb un missatge d'error clar en lloc de manifestar-se com una fallada críptica del procés a mig camí.
La formació de [Hello nf-core](../../hello_nf-core/05_input_validation.md) explica com afegir validació d'entrada utilitzant el connector nf-schema. <!-- TODO (futur) pendent d'una missió secundària de Validació pròpia -->

**2. Entrades explícites del procés per als valors requerits**

Si voleu que la interfície del procés comuniqui que un valor determinat és obligatori, considereu extreure'l del meta map com a entrada explícita:

=== "Definició del procés"

    ```groovy title="modules/cowpy.nf" linenums="6"
    input:
    tuple val(meta), val(character), path(input_file)
    ```

=== "Crida al workflow"

    ```groovy title="main.nf"
    COWPY(ch_languages.map { meta, file -> [meta, meta.character, file] })
    ```

Aquest enfocament fa que `character` sigui una part visible i requerida del contracte del procés.
Qualsevol persona que llegeixi el mòdul pot veure immediatament que s'ha de proporcionar un valor de personatge.
Si el camp és absent, el workflow falla clarament a nivell del canal abans que el procés s'executi.

Això posa de manifest un principi de disseny útil:

**Utilitzeu el meta map per a informació opcional o descriptiva; extraieu els valors requerits com a entrades explícites.**

El meta map manté les estructures del canal netes i estables, però per als valors que un procés realment necessita, exposar-los com a entrades amb nom millora la claredat i fa que el mòdul sigui més fàcil d'utilitzar correctament en altres contextos.

### Conclusió

En aquesta secció, heu vist:

- **Com es manifesten les metadades que falten:** Un camp buit produeix un argument buit; un camp absent produeix `null` substituït literalment a la comanda.
- **Dues estratègies complementàries:** Validació d'entrada per detectar problemes aviat, i entrades explícites del procés per comunicar els requisits clarament.

---

## Resum

En aquesta missió secundària, heu explorat com treballar eficaçment amb metadades en workflows de Nextflow.

El patró de tupla "meta map + fitxer de dades" és una convenció fonamental en Nextflow, que ofereix diversos avantatges respecte a passar metadades com a valors individuals:

- L'estructura del canal es manté estable a mesura que el full de dades evoluciona
- El comportament del procés es pot personalitzar per mostra sense codificar noms de camps
- Les metadades estan disponibles al llarg del pipeline per a l'anomenament, l'agrupació i l'organització de les sortides
- Els mòduls escrits per a aquesta interfície són intercanviables, inclosos els mòduls nf-core

### Patrons clau

1.  **Lectura i estructuració de metadades:** Analitzar un full de dades CSV i crear un meta map.

    ```groovy
    channel.fromPath('datasheet.csv')
        .splitCsv(header: true)
        .map { row ->
            [ [id: row.id, character: row.character], row.recording ]
        }
    ```

2.  **Ampliació de metadades durant el workflow:** Afegir noves claus a partir de les sortides del procés o de lògica derivada.

    ```groovy
    // A partir d'una sortida del procés
    .map { meta, file, lang ->
        [ meta + [lang: lang], file ]
    }

    // A partir de lògica condicional
    .map { meta, file ->
        def lang_group = "unknown"
        if (meta.lang in ["de", "en"]) { lang_group = "germanic" }
        else if (meta.lang in ["fr", "es", "it"]) { lang_group = "romance" }
        [ meta + [lang_group: lang_group], file ]
    }
    ```

3.  **Ús de metadades dins d'un procés:** Accediu a qualsevol camp mitjançant la notació de punt al bloc script.

    ```groovy
    cat ${input_file} | cowpy -c ${meta.character} > ${meta.lang}-${input_file}
    ```

4.  **Organització de les sortides per valor de metadades:** Utilitzeu la closure `path` al bloc `output {}`.

    ```groovy
    output {
        cowpy_art {
            path { meta, file -> meta.lang_group }
        }
    }
    ```

### Recursos addicionals

- [operador map](https://www.nextflow.io/docs/latest/operator.html#map)
- [operador multiMap](https://www.nextflow.io/docs/latest/reference/operator.html#multimap)
- [qualificador de sortida stdout](https://www.nextflow.io/docs/latest/process.html#outputs)

---

## Què segueix?

Torneu al [menú de missions secundàries](../index.md) o feu clic al botó a la part inferior dreta de la pàgina per continuar amb el tema següent de la llista.
