# Divisió i Agrupació

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducció assistida per IA - [més informació i suggeriments](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Nextflow proporciona eines potents per treballar amb dades de manera flexible. Una capacitat clau és dividir les dades en diferents fluxos i després agrupar els elements relacionats de nou. Això és especialment valuós en workflows de bioinformàtica on necessiteu processar diferents tipus de mostres per separat abans de combinar els resultats per a l'anàlisi.

Penseu-ho com classificar el correu: separeu les cartes per destinació, processeu cada pila de manera diferent i després recombineu els elements que van a la mateixa persona. Nextflow utilitza operadors especials per aconseguir això amb dades científiques. Aquest enfocament també es coneix comunament com el patró **scatter/gather** en computació distribuïda i workflows de bioinformàtica.

El sistema de canals de Nextflow és el cor d'aquesta flexibilitat. Els canals connecten diferents parts del vostre workflow, permetent que les dades flueixin a través de la vostra anàlisi. Podeu crear múltiples canals des d'una única font de dades, processar cada canal de manera diferent i després fusionar els canals de nou quan sigui necessari. Aquest enfocament us permet dissenyar workflows que reflecteixen naturalment els camins de ramificació i convergència d'anàlisis bioinformàtiques complexes.

### Objectius d'aprenentatge

En aquesta missió secundària, aprendreu a dividir i agrupar dades utilitzant els operadors de canal de Nextflow.
Començarem amb un fitxer CSV que conté informació de mostres i fitxers de dades associats, després manipularem i reorganitzarem aquestes dades.

Al final d'aquesta missió secundària, podreu separar i combinar fluxos de dades de manera efectiva, utilitzant les tècniques següents:

- Llegir dades de fitxers utilitzant `splitCsv`
- Filtrar i transformar dades amb `filter` i `map`
- Combinar dades relacionades utilitzant `join` i `groupTuple`
- Crear combinacions de dades amb `combine` per a processament paral·lel
- Optimitzar l'estructura de dades utilitzant `subMap` i estratègies de deduplicació
- Construir funcions reutilitzables amb closures amb nom per ajudar-vos a manipular estructures de canal

Aquestes habilitats us ajudaran a construir workflows que poden gestionar múltiples fitxers d'entrada i diferents tipus de dades de manera eficient, mentre manteniu una estructura de codi neta i mantenible.

### Prerequisits

Abans d'emprendre aquesta missió secundària, hauríeu de:

- Haver completat el tutorial [Hello Nextflow](../hello_nextflow/README.md) o un curs equivalent per a principiants.
- Sentir-vos còmodes utilitzant conceptes i mecanismes bàsics de Nextflow (processos, canals, operadors, treballar amb fitxers, metadades)

**Opcional:** Recomanem completar primer la missió secundària [Metadades en workflows](./metadata.md).
Això cobreix els fonaments de llegir fitxers CSV amb `splitCsv` i crear mapes de metadades, que utilitzarem intensament aquí.

---

## 0. Primers passos

#### Obriu l'espai de codi de formació

Si encara no ho heu fet, assegureu-vos d'obrir l'entorn de formació tal com es descriu a [Configuració de l'entorn](../envsetup/index.md).

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### Moveu-vos al directori del projecte

Movem-nos al directori on es troben els fitxers per a aquest tutorial.

```bash
cd side-quests/splitting_and_grouping
```

Podeu configurar VSCode per centrar-vos en aquest directori:

```bash
code .
```

#### Reviseu els materials

Trobareu un fitxer de workflow principal i un directori `data` que conté un full de mostres anomenat `samplesheet.csv`.

```console title="Contingut del directori"
.
├── data
│   └── samplesheet.csv
└── main.nf
```

El full de mostres conté informació sobre mostres de diferents pacients, incloent l'ID del pacient, el número de repetició de la mostra, el tipus (normal o tumor) i els camins a fitxers de dades hipotètics (que en realitat no existeixen, però farem veure que sí).

```console title="samplesheet.csv"
id,repeat,type,bam
patientA,1,normal,patientA_rep1_normal.bam
patientA,1,tumor,patientA_rep1_tumor.bam
patientA,2,normal,patientA_rep2_normal.bam
patientA,2,tumor,patientA_rep2_tumor.bam
patientB,1,normal,patientB_rep1_normal.bam
patientB,1,tumor,patientB_rep1_tumor.bam
patientC,1,normal,patientC_rep1_normal.bam
patientC,1,tumor,patientC_rep1_tumor.bam
```

Aquest full de mostres llista vuit mostres de tres pacients (A, B, C).

Per a cada pacient, tenim mostres que són de tipus `tumor` (típicament originades de biòpsies de tumor) o `normal` (preses de teixit sa o sang).
Si no esteu familiaritzats amb l'anàlisi de càncer, només cal saber que això correspon a un model experimental que utilitza mostres aparellades tumor/normal per realitzar anàlisis contrastives.

Per al pacient A específicament, tenim dos conjunts de rèpliques tècniques (repeticions).

!!! note "Nota"

    No us preocupeu si no esteu familiaritzats amb aquest disseny experimental, no és crític per entendre aquest tutorial.

#### Reviseu l'assignació

El vostre repte és escriure un workflow de Nextflow que:

1. **Llegeixi** dades de mostres d'un fitxer CSV i les estructuri amb mapes de metadades
2. **Separi** mostres en diferents canals segons el tipus (normal vs tumor)
3. **Uneixi** parelles aparellades tumor/normal per ID de pacient i número de rèplica
4. **Distribueixi** mostres a través d'intervals genòmics per a processament paral·lel
5. **Agrupí** mostres relacionades de nou per a anàlisi posterior

Això representa un patró comú de bioinformàtica on necessiteu dividir dades per a processament independent, després recombinar elements relacionats per a anàlisi comparativa.

#### Llista de verificació de preparació

Creieu que esteu preparats per començar?

- [ ] Entenc l'objectiu d'aquest curs i els seus prerequisits
- [ ] El meu espai de codi està en funcionament
- [ ] He establert el meu directori de treball adequadament
- [ ] Entenc l'assignació

Si podeu marcar totes les caselles, esteu preparats per començar.

---

## 1. Llegir dades de mostres

### 1.1. Llegir dades de mostres amb `splitCsv` i crear mapes de metadades

Comencem llegint les dades de mostres amb `splitCsv` i organitzant-les en el patró de mapa de metadades. Al `main.nf`, veureu que ja hem començat el workflow.

```groovy title="main.nf" linenums="1" hl_lines="2"
workflow {
    ch_samplesheet = channel.fromPath("./data/samplesheet.csv")
}
```

!!! note "Nota"

    Al llarg d'aquest tutorial, utilitzarem el prefix `ch_` per a totes les variables de canal per indicar clarament que són canals de Nextflow.

Si heu completat la missió secundària [Metadades en workflows](./metadata.md), reconeixereu aquest patró. Utilitzarem `splitCsv` per llegir el CSV i estructurar immediatament les dades amb un mapa de metadades per separar les metadades dels camins de fitxers.

!!! info "Info"

    Trobarem dos conceptes diferents anomenats `map` en aquesta formació:

    - **Estructura de dades**: El mapa de Groovy (equivalent a diccionaris/hashs en altres llenguatges) que emmagatzema parells clau-valor
    - **Operador de canal**: L'operador `.map()` que transforma elements en un canal

    Aclararem quin volem dir en context, però aquesta distinció és important per entendre quan es treballa amb Nextflow.

Apliqueu aquests canvis a `main.nf`:

=== "Després"

    ```groovy title="main.nf" linenums="2" hl_lines="2-6"
        ch_samples = channel.fromPath("./data/samplesheet.csv")
            .splitCsv(header: true)
            .map{ row ->
              [[id:row.id, repeat:row.repeat, type:row.type], row.bam]
            }
            .view()
    ```

=== "Abans"

    ```groovy title="main.nf" linenums="2" hl_lines="1"
        ch_samplesheet = channel.fromPath("./data/samplesheet.csv")
    ```

Això combina l'operació `splitCsv` (llegir el CSV amb capçaleres) i l'operació `map` (estructurar dades com a tuples `[meta, file]`) en un sol pas. Apliqueu aquest canvi i executeu el pipeline:

```bash
nextflow run main.nf
```

??? success "Sortida de la comanda"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [deadly_mercator] DSL2 - revision: bd6b0224e9

    [[id:patientA, repeat:1, type:normal], patientA_rep1_normal.bam]
    [[id:patientA, repeat:1, type:tumor], patientA_rep1_tumor.bam]
    [[id:patientA, repeat:2, type:normal], patientA_rep2_normal.bam]
    [[id:patientA, repeat:2, type:tumor], patientA_rep2_tumor.bam]
    [[id:patientB, repeat:1, type:normal], patientB_rep1_normal.bam]
    [[id:patientB, repeat:1, type:tumor], patientB_rep1_tumor.bam]
    [[id:patientC, repeat:1, type:normal], patientC_rep1_normal.bam]
    [[id:patientC, repeat:1, type:tumor], patientC_rep1_tumor.bam]
    ```

Ara tenim un canal on cada element és una tupla `[meta, file]` - metadades separades dels camins de fitxers. Aquesta estructura ens permet dividir i agrupar la nostra càrrega de treball basant-nos en camps de metadades.

---

## 2. Filtrar i transformar dades

### 2.1. Filtrar dades amb `filter`

Podem utilitzar l'[operador `filter`](https://www.nextflow.io/docs/latest/operator.html#filter) per filtrar les dades segons una condició. Diguem que només volem processar mostres normals. Podem fer això filtrant les dades segons el camp `type`. Inserim això abans de l'operador `view`.

=== "Després"

    ```groovy title="main.nf" linenums="2" hl_lines="6"
        ch_samples = channel.fromPath("./data/samplesheet.csv")
            .splitCsv(header: true)
            .map{ row ->
              [[id:row.id, repeat:row.repeat, type:row.type], row.bam]
            }
            .filter { meta, file -> meta.type == 'normal' }
            .view()
    ```

=== "Abans"

    ```groovy title="main.nf" linenums="2"
        ch_samples = channel.fromPath("./data/samplesheet.csv")
            .splitCsv(header: true)
            .map{ row ->
              [[id:row.id, repeat:row.repeat, type:row.type], row.bam]
            }
            .view()
    ```

Executeu el workflow de nou per veure el resultat filtrat:

```bash
nextflow run main.nf
```

??? success "Sortida de la comanda"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [admiring_brown] DSL2 - revision: 194d61704d

    [[id:patientA, repeat:1, type:normal], patientA_rep1_normal.bam]
    [[id:patientA, repeat:2, type:normal], patientA_rep2_normal.bam]
    [[id:patientB, repeat:1, type:normal], patientB_rep1_normal.bam]
    [[id:patientC, repeat:1, type:normal], patientC_rep1_normal.bam]
    ```

Hem filtrat amb èxit les dades per incloure només mostres normals. Recapitulem com funciona això.

L'operador `filter` pren una closure que s'aplica a cada element del canal. Si la closure retorna `true`, l'element s'inclou; si retorna `false`, l'element s'exclou.

En el nostre cas, volem mantenir només mostres on `meta.type == 'normal'`. La closure utilitza la tupla `meta,file` per referir-se a cada mostra, accedeix al tipus de mostra amb `meta.type` i comprova si és igual a `'normal'`.

Això s'aconsegueix amb la única closure que hem introduït anteriorment:

```groovy title="main.nf" linenums="7"
    .filter { meta, file -> meta.type == 'normal' }
```

### 2.2. Crear canals filtrats separats

Actualment estem aplicant el filtre al canal creat directament des del CSV, però volem filtrar això de més d'una manera, així que reescrivim la lògica per crear un canal filtrat separat per a mostres normals:

=== "Després"

    ```groovy title="main.nf" linenums="2" hl_lines="6 8"
        ch_samples = channel.fromPath("./data/samplesheet.csv")
            .splitCsv(header: true)
            .map{ row ->
                [[id:row.id, repeat:row.repeat, type:row.type], row.bam]
            }
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
        ch_normal_samples
            .view()
    ```

=== "Abans"

    ```groovy title="main.nf" linenums="2"
        ch_samples = channel.fromPath("./data/samplesheet.csv")
            .splitCsv(header: true)
            .map{ row ->
              [[id:row.id, repeat:row.repeat, type:row.type], row.bam]
            }
            .filter { meta, file -> meta.type == 'normal' }
            .view()
    ```

Executeu el pipeline per veure els resultats:

```bash
nextflow run main.nf
```

??? success "Sortida de la comanda"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [trusting_poisson] DSL2 - revision: 639186ee74

    [[id:patientA, repeat:1, type:normal], patientA_rep1_normal.bam]
    [[id:patientA, repeat:2, type:normal], patientA_rep2_normal.bam]
    [[id:patientB, repeat:1, type:normal], patientB_rep1_normal.bam]
    [[id:patientC, repeat:1, type:normal], patientC_rep1_normal.bam]
    ```

Hem filtrat amb èxit les dades i creat un canal separat per a mostres normals.

Creem també un canal filtrat per a les mostres de tumor:

=== "Després"

    ```groovy title="main.nf" linenums="7" hl_lines="3-8"
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
        ch_tumor_samples = ch_samples
            .filter { meta, file -> meta.type == 'tumor' }
        ch_normal_samples
            .view{'Normal sample: ' + it}
        ch_tumor_samples
            .view{'Tumor sample: ' + it}
    ```

=== "Abans"

    ```groovy title="main.nf" linenums="7" hl_lines="3 4"
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
        ch_normal_samples
            .view()
    ```

```bash
nextflow run main.nf
```

??? success "Sortida de la comanda"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [maniac_boltzmann] DSL2 - revision: 3636b6576b

    Tumor sample: [[id:patientA, repeat:1, type:tumor], patientA_rep1_tumor.bam]
    Tumor sample: [[id:patientA, repeat:2, type:tumor], patientA_rep2_tumor.bam]
    Normal sample: [[id:patientA, repeat:1, type:normal], patientA_rep1_normal.bam]
    Normal sample: [[id:patientA, repeat:2, type:normal], patientA_rep2_normal.bam]
    Normal sample: [[id:patientB, repeat:1, type:normal], patientB_rep1_normal.bam]
    Normal sample: [[id:patientC, repeat:1, type:normal], patientC_rep1_normal.bam]
    Tumor sample: [[id:patientB, repeat:1, type:tumor], patientB_rep1_tumor.bam]
    Tumor sample: [[id:patientC, repeat:1, type:tumor], patientC_rep1_tumor.bam]
    ```

Hem separat les mostres normals i de tumor en dos canals diferents, i hem utilitzat una closure proporcionada a `view()` per etiquetar-les de manera diferent a la sortida: `ch_tumor_samples.view{'Tumor sample: ' + it}`.

### Conclusió

En aquesta secció, heu après:

- **Filtrar dades**: Com filtrar dades amb `filter`
- **Dividir dades**: Com dividir dades en diferents canals segons una condició
- **Visualitzar dades**: Com utilitzar `view` per imprimir les dades i etiquetar la sortida de diferents canals

Ara hem separat les mostres normals i de tumor en dos canals diferents. A continuació, unirem les mostres normals i de tumor pel camp `id`.

---

## 3. Unir canals per identificadors

A la secció anterior, hem separat les mostres normals i de tumor en dos canals diferents. Aquestes podrien processar-se independentment utilitzant processos o workflows específics segons el seu tipus. Però què passa quan volem comparar les mostres normals i de tumor del mateix pacient? En aquest punt, necessitem unir-les de nou assegurant-nos de fer coincidir les mostres segons el seu camp `id`.

Nextflow inclou molts mètodes per combinar canals, però en aquest cas l'operador més apropiat és [`join`](https://www.nextflow.io/docs/latest/operator.html#join). Si esteu familiaritzats amb SQL, actua com l'operació `JOIN`, on especifiquem la clau per unir i el tipus d'unió a realitzar.

### 3.1. Utilitzar `map` i `join` per combinar segons l'ID del pacient

#### 3.1.1. Comprovar l'estructura de dades

Si no teniu la sortida de consola encara disponible, executem el pipeline per comprovar la nostra estructura de dades i veure com necessitem modificar-la per unir pel camp `id`.

```bash
nextflow run main.nf
```

??? success "Sortida de la comanda"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [maniac_boltzmann] DSL2 - revision: 3636b6576b

    Tumor sample: [[id:patientA, repeat:1, type:tumor], patientA_rep1_tumor.bam]
    Tumor sample: [[id:patientA, repeat:2, type:tumor], patientA_rep2_tumor.bam]
    Normal sample: [[id:patientA, repeat:1, type:normal], patientA_rep1_normal.bam]
    Normal sample: [[id:patientA, repeat:2, type:normal], patientA_rep2_normal.bam]
    Normal sample: [[id:patientB, repeat:1, type:normal], patientB_rep1_normal.bam]
    Normal sample: [[id:patientC, repeat:1, type:normal], patientC_rep1_normal.bam]
    Tumor sample: [[id:patientB, repeat:1, type:tumor], patientB_rep1_tumor.bam]
    Tumor sample: [[id:patientC, repeat:1, type:tumor], patientC_rep1_tumor.bam]
    ```

Podem veure que el camp `id` és el primer element de cada mapa de metadades. Perquè `join` funcioni, hauríem d'aïllar el camp `id` en cada tupla. Després d'això, podem simplement utilitzar l'operador `join` per combinar els dos canals.

#### 3.1.2. Aïllar el camp `id`

Per aïllar el camp `id`, podem utilitzar l'[operador `map`](https://www.nextflow.io/docs/latest/operator.html#map) per crear una nova tupla amb el camp `id` com a primer element.

=== "Després"

    ```groovy title="main.nf" linenums="7" hl_lines="3 6"
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
            .map { meta, file -> [meta.id, meta, file] }
        ch_tumor_samples = ch_samples
            .filter { meta, file -> meta.type == 'tumor' }
            .map { meta, file -> [meta.id, meta, file] }
        ch_normal_samples
            .view{'Normal sample: ' + it}
        ch_tumor_samples
            .view{'Tumor sample: ' + it}
    ```

=== "Abans"

    ```groovy title="main.nf" linenums="7"
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
        ch_tumor_samples = ch_samples
            .filter { meta, file -> meta.type == 'tumor' }
        ch_normal_samples
            .view{'Normal sample: ' + it}
        ch_tumor_samples
            .view{'Tumor sample: ' + it}
    ```

```bash
nextflow run main.nf
```

??? success "Sortida de la comanda"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [mad_lagrange] DSL2 - revision: 9940b3f23d

    Tumor sample: [patientA, [id:patientA, repeat:1, type:tumor], patientA_rep1_tumor.bam]
    Tumor sample: [patientA, [id:patientA, repeat:2, type:tumor], patientA_rep2_tumor.bam]
    Normal sample: [patientA, [id:patientA, repeat:1, type:normal], patientA_rep1_normal.bam]
    Normal sample: [patientA, [id:patientA, repeat:2, type:normal], patientA_rep2_normal.bam]
    Tumor sample: [patientB, [id:patientB, repeat:1, type:tumor], patientB_rep1_tumor.bam]
    Tumor sample: [patientC, [id:patientC, repeat:1, type:tumor], patientC_rep1_tumor.bam]
    Normal sample: [patientB, [id:patientB, repeat:1, type:normal], patientB_rep1_normal.bam]
    Normal sample: [patientC, [id:patientC, repeat:1, type:normal], patientC_rep1_normal.bam]
    ```

Pot ser subtil, però hauríeu de poder veure que el primer element de cada tupla és el camp `id`.

#### 3.1.3. Combinar els dos canals

Ara podem utilitzar l'operador `join` per combinar els dos canals segons el camp `id`.

Un cop més, utilitzarem `view` per imprimir les sortides unides.

=== "Després"

    ```groovy title="main.nf" linenums="7" hl_lines="7-9"
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
            .map { meta, file -> [meta.id, meta, file] }
        ch_tumor_samples = ch_samples
            .filter { meta, file -> meta.type == 'tumor' }
            .map { meta, file -> [meta.id, meta, file] }
        ch_joined_samples = ch_normal_samples
            .join(ch_tumor_samples)
        ch_joined_samples.view()
    ```

=== "Abans"

    ```groovy title="main.nf" linenums="7" hl_lines="7-10"
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
            .map { meta, file -> [meta.id, meta, file] }
        ch_tumor_samples = ch_samples
            .filter { meta, file -> meta.type == 'tumor' }
            .map { meta, file -> [meta.id, meta, file] }
        ch_normal_samples
            .view{'Normal sample: ' + it}
        ch_tumor_samples
            .view{'Tumor sample: ' + it}
    ```

```bash
nextflow run main.nf
```

??? success "Sortida de la comanda"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [soggy_wiles] DSL2 - revision: 3bc1979889

    [patientA, [id:patientA, repeat:1, type:normal], patientA_rep1_normal.bam, [id:patientA, repeat:1, type:tumor], patientA_rep1_tumor.bam]
    [patientA, [id:patientA, repeat:2, type:normal], patientA_rep2_normal.bam, [id:patientA, repeat:2, type:tumor], patientA_rep2_tumor.bam]
    [patientB, [id:patientB, repeat:1, type:normal], patientB_rep1_normal.bam, [id:patientB, repeat:1, type:tumor], patientB_rep1_tumor.bam]
    [patientC, [id:patientC, repeat:1, type:normal], patientC_rep1_normal.bam, [id:patientC, repeat:1, type:tumor], patientC_rep1_tumor.bam]
    ```

És una mica difícil de veure perquè és molt ample, però hauríeu de poder veure que les mostres s'han unit pel camp `id`. Cada tupla ara té el format:

- `id`: L'ID de la mostra
- `normal_meta_map`: Les metadades de la mostra normal incloent tipus, rèplica i camí al fitxer bam
- `normal_sample_file`: El fitxer de mostra normal
- `tumor_meta_map`: Les metadades de la mostra de tumor incloent tipus, rèplica i camí al fitxer bam
- `tumor_sample`: La mostra de tumor incloent tipus, rèplica i camí al fitxer bam

!!! warning "Advertència"

    L'operador `join` descartarà qualsevol tupla no coincident. En aquest exemple, ens hem assegurat que totes les mostres coincidissin per a tumor i normal, però si això no és cert heu d'utilitzar el paràmetre `remainder: true` per mantenir les tuples no coincidents. Consulteu la [documentació](https://www.nextflow.io/docs/latest/operator.html#join) per a més detalls.

Així que ara sabeu com utilitzar `map` per aïllar un camp en una tupla, i com utilitzar `join` per combinar tuples segons el primer camp.
Amb aquest coneixement, podem combinar amb èxit canals segons un camp compartit.

A continuació, considerarem la situació on voleu unir per múltiples camps.

### 3.2. Unir per múltiples camps

Tenim 2 rèpliques per a sampleA, però només 1 per a sampleB i sampleC. En aquest cas vam poder unir-les efectivament utilitzant el camp `id`, però què passaria si estiguessin desincronitzades? Podríem barrejar les mostres normals i de tumor de diferents rèpliques!

Per evitar això, podem unir per múltiples camps. En realitat hi ha múltiples maneres d'aconseguir això, però ens centrarem en crear una nova clau d'unió que inclogui tant l'`id` de la mostra com el número de `replicate`.

Comencem creant una nova clau d'unió. Podem fer això de la mateixa manera que abans utilitzant l'[operador `map`](https://www.nextflow.io/docs/latest/operator.html#map) per crear una nova tupla amb els camps `id` i `repeat` com a primer element.

=== "Després"

    ```groovy title="main.nf" linenums="7" hl_lines="3 6"
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
            .map { meta, file -> [[meta.id, meta.repeat], meta, file] }
        ch_tumor_samples = ch_samples
            .filter { meta, file -> meta.type == 'tumor' }
            .map { meta, file -> [[meta.id, meta.repeat], meta, file] }
    ```

=== "Abans"

    ```groovy title="main.nf" linenums="7" hl_lines="3 6"
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
            .map { meta, file -> [meta.id, meta, file] }
        ch_tumor_samples = ch_samples
            .filter { meta, file -> meta.type == 'tumor' }
            .map { meta, file -> [meta.id, meta, file] }
    ```

Ara hauríem de veure que la unió s'està produint però utilitzant tant els camps `id` com `repeat`. Executeu el workflow:

```bash
nextflow run main.nf
```

??? success "Sortida de la comanda"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [prickly_wing] DSL2 - revision: 3bebf22dee

    [[patientA, 1], [id:patientA, repeat:1, type:normal], patientA_rep1_normal.bam, [id:patientA, repeat:1, type:tumor], patientA_rep1_tumor.bam]
    [[patientA, 2], [id:patientA, repeat:2, type:normal], patientA_rep2_normal.bam, [id:patientA, repeat:2, type:tumor], patientA_rep2_tumor.bam]
    [[patientB, 1], [id:patientB, repeat:1, type:normal], patientB_rep1_normal.bam, [id:patientB, repeat:1, type:tumor], patientB_rep1_tumor.bam]
    [[patientC, 1], [id:patientC, repeat:1, type:normal], patientC_rep1_normal.bam, [id:patientC, repeat:1, type:tumor], patientC_rep1_tumor.bam]
    ```

Noteu com tenim una tupla de dos elements (camps `id` i `repeat`) com a primer element de cada resultat unit. Això demostra com elements complexos poden utilitzar-se com a clau d'unió, permetent coincidències força intrincades entre mostres de les mateixes condicions.

Si voleu explorar més maneres d'unir per diferents claus, consulteu la [documentació de l'operador join](https://www.nextflow.io/docs/latest/operator.html#join) per a opcions i exemples addicionals.

### 3.3. Utilitzar `subMap` per crear una nova clau d'unió

L'enfocament anterior perd els noms de camp de la nostra clau d'unió - els camps `id` i `repeat` es converteixen només en una llista de valors. Per retenir els noms de camp per a accés posterior, podem utilitzar el [mètode `subMap`](<https://docs.groovy-lang.org/latest/html/groovy-jdk/java/util/Map.html#subMap(java.util.Collection)>).

El mètode `subMap` extreu només els parells clau-valor especificats d'un mapa. Aquí extraurem només els camps `id` i `repeat` per crear la nostra clau d'unió.

=== "Després"

    ```groovy title="main.nf" linenums="7" hl_lines="3 6"
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
            .map { meta, file -> [meta.subMap(['id', 'repeat']), meta, file] }
        ch_tumor_samples = ch_samples
            .filter { meta, file -> meta.type == 'tumor' }
            .map { meta, file -> [meta.subMap(['id', 'repeat']), meta, file] }
    ```

=== "Abans"

    ```groovy title="main.nf" linenums="7" hl_lines="3 6"
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
            .map { meta, file -> [[meta.id, meta.repeat], meta, file] }
        ch_tumor_samples = ch_samples
            .filter { meta, file -> meta.type == 'tumor' }
            .map { meta, file -> [[meta.id, meta.repeat], meta, file] }
    ```

```bash
nextflow run main.nf
```

??? success "Sortida de la comanda"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [reverent_wing] DSL2 - revision: 847016c3b7

    [[id:patientA, repeat:1], [id:patientA, repeat:1, type:normal], patientA_rep1_normal.bam, [id:patientA, repeat:1, type:tumor], patientA_rep1_tumor.bam]
    [[id:patientA, repeat:2], [id:patientA, repeat:2, type:normal], patientA_rep2_normal.bam, [id:patientA, repeat:2, type:tumor], patientA_rep2_tumor.bam]
    [[id:patientB, repeat:1], [id:patientB, repeat:1, type:normal], patientB_rep1_normal.bam, [id:patientB, repeat:1, type:tumor], patientB_rep1_tumor.bam]
    [[id:patientC, repeat:1], [id:patientC, repeat:1, type:normal], patientC_rep1_normal.bam, [id:patientC, repeat:1, type:tumor], patientC_rep1_tumor.bam]
    ```

Ara tenim una nova clau d'unió que no només inclou els camps `id` i `repeat` sinó que també conserva els noms de camp perquè puguem accedir-hi més tard per nom, per exemple `meta.id` i `meta.repeat`.

### 3.4. Utilitzar una closure amb nom en map

Per evitar duplicació i reduir errors, podem utilitzar una closure amb nom. Una closure amb nom ens permet crear una funció reutilitzable que podem cridar en múltiples llocs.

Per fer-ho, primer definim la closure com una nova variable:

=== "Després"

    ```groovy title="main.nf" linenums="2" hl_lines="7"
        ch_samples = channel.fromPath("./data/samplesheet.csv")
            .splitCsv(header: true)
            .map{ row ->
                [[id:row.id, repeat:row.repeat, type:row.type], row.bam]
            }

        getSampleIdAndReplicate = { meta, bam -> [ meta.subMap(['id', 'repeat']), meta, file(bam) ] }

        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
    ```

=== "Abans"

    ```groovy title="main.nf" linenums="2"
        ch_samples = channel.fromPath("./data/samplesheet.csv")
            .splitCsv(header: true)
            .map{ row ->
                [[id:row.id, repeat:row.repeat, type:row.type], row.bam]
            }
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
    ```

Hem definit la transformació map com una variable amb nom que podem reutilitzar.

Noteu que també convertim el camí del fitxer a un objecte Path utilitzant `file()` perquè qualsevol procés que rebi aquest canal pugui gestionar el fitxer correctament (per a més informació vegeu [Treballar amb fitxers](./working_with_files.md)).

Implementem la closure al nostre workflow:

=== "Després"

    ```groovy title="main.nf" linenums="10" hl_lines="3 6"
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
             .map ( getSampleIdAndReplicate )
        ch_tumor_samples = ch_samples
            .filter { meta, file -> meta.type == 'tumor' }
             .map ( getSampleIdAndReplicate )

    ```

=== "Abans"

    ```groovy title="main.nf" linenums="10" hl_lines="3 6"
        ch_normal_samples = ch_samples
            .filter { meta, file -> meta.type == 'normal' }
            .map { meta, file -> [meta.subMap(['id', 'repeat']), meta, file] }
        ch_tumor_samples = ch_samples
            .filter { meta, file -> meta.type == 'tumor' }
            .map { meta, file -> [meta.subMap(['id', 'repeat']), meta, file] }
    ```

!!! note "Nota"

    L'operador `map` ha canviat d'utilitzar `{ }` a utilitzar `( )` per passar la closure com a argument. Això és perquè l'operador `map` espera una closure com a argument i `{ }` s'utilitza per definir una closure anònima. Quan es crida una closure amb nom, utilitzeu la sintaxi `( )`.

Executeu el workflow un cop més per comprovar que tot segueix funcionant:

```bash
nextflow run main.nf
```

??? success "Sortida de la comanda"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [angry_meninsky] DSL2 - revision: 2edc226b1d

    [[id:patientA, repeat:1], [id:patientA, repeat:1, type:normal], patientA_rep1_normal.bam, [id:patientA, repeat:1, type:tumor], patientA_rep1_tumor.bam]
    [[id:patientA, repeat:2], [id:patientA, repeat:2, type:normal], patientA_rep2_normal.bam, [id:patientA, repeat:2, type:tumor], patientA_rep2_tumor.bam]
    [[id:patientB, repeat:1], [id:patientB, repeat:1, type:normal], patientB_rep1_normal.bam, [id:patientB, repeat:1, type:tumor], patientB_rep1_tumor.bam]
    [[id:patientC, repeat:1], [id:patientC, repeat:1, type:normal], patientC_rep1_normal.bam, [id:patientC, repeat:1, type:tumor], patientC_rep1_tumor.bam]
    ```

Utilitzar una closure amb nom ens permet reutilitzar la mateixa transformació en múltiples llocs, reduint el risc d'errors i fent el codi més llegible i mantenible.

### 3.5. Reduir la duplicació de dades

Tenim moltes dades duplicades al nostre workflow. Cada element de les mostres unides repeteix els camps `id` i `repeat`. Com que aquesta informació ja està disponible a la clau d'agrupació, podem evitar aquesta redundància. Com a recordatori, la nostra estructura de dades actual té aquest aspecte:

```groovy
[
  [
    "id": "sampleC",
    "repeat": "1",
  ],
  [
    "id": "sampleC",
    "repeat": "1",
    "type": "normal",
  ],
  "sampleC_rep1_normal.bam"
  [
    "id": "sampleC",
    "repeat": "1",
    "type": "tumor",
  ],
  "sampleC_rep1_tumor.bam"
]
```

Com que els camps `id` i `repeat` estan disponibles a la clau d'agrupació, eliminem-los de la resta de cada element del canal per evitar duplicació. Podem fer això utilitzant el mètode `subMap` per crear un nou mapa amb només el camp `type`. Aquest enfocament ens permet mantenir tota la informació necessària mentre eliminem la redundància a la nostra estructura de dades.

=== "Després"

    ```groovy title="main.nf" linenums="8" hl_lines="1"
        getSampleIdAndReplicate = { meta, bam -> [ meta.subMap(['id', 'repeat']), meta.subMap(['type']), file(bam) ] }
    ```

=== "Abans"

    ```groovy title="main.nf" linenums="8" hl_lines="1"
        getSampleIdAndReplicate = { meta, bam -> [ meta.subMap(['id', 'repeat']), meta, file(bam) ] }
    ```

Ara la closure retorna una tupla on el primer element conté els camps `id` i `repeat`, i el segon element conté només el camp `type`. Això elimina la redundància emmagatzemant la informació d'`id` i `repeat` un cop a la clau d'agrupació, mentre es manté tota la informació necessària.

Executeu el workflow per veure com es veu això:

```bash
nextflow run main.nf
```

??? success "Sortida de la comanda"

    ```console
    [[id:patientA, repeat:1], [type:normal], /workspaces/training/side-quests/splitting_and_grouping/patientA_rep1_normal.bam, [type:tumor], /workspaces/training/side-quests/splitting_and_grouping/patientA_rep1_tumor.bam]
    [[id:patientA, repeat:2], [type:normal], /workspaces/training/side-quests/splitting_and_grouping/patientA_rep2_normal.bam, [type:tumor], /workspaces/training/side-quests/splitting_and_grouping/patientA_rep2_tumor.bam]
    [[id:patientB, repeat:1], [type:normal], /workspaces/training/side-quests/splitting_and_grouping/patientB_rep1_normal.bam, [type:tumor], /workspaces/training/side-quests/splitting_and_grouping/patientB_rep1_tumor.bam]
    [[id:patientC, repeat:1], [type:normal], /workspaces/training/side-quests/splitting_and_grouping/patientC_rep1_normal.bam, [type:tumor], /workspaces/training/side-quests/splitting_and_grouping/patientC_rep1_tumor.bam]
    ```

Podem veure que només declarem els camps `id` i `repeat` un cop a la clau d'agrupació i tenim el camp `type` a les dades de mostra. No hem perdut cap informació però hem aconseguit fer els continguts del nostre canal més concisos.

### 3.6. Eliminar informació redundant

Hem eliminat informació duplicada anteriorment, però encara tenim alguna altra informació redundant als nostres canals.

Al principi, vam separar les mostres normals i de tumor utilitzant `filter`, després les vam unir segons les claus `id` i `repeat`. L'operador `join` preserva l'ordre en què es fusionen les tuples, així que en el nostre cas, amb mostres normals al costat esquerre i mostres de tumor a la dreta, el canal resultant manté aquesta estructura: `id, <elements normals>, <elements tumor>`.

Com que coneixem la posició de cada element al nostre canal, podem simplificar encara més l'estructura eliminant les metadades `[type:normal]` i `[type:tumor]`.

=== "Després"

    ```groovy title="main.nf" linenums="8" hl_lines="1"
        getSampleIdAndReplicate = { meta, file -> [ meta.subMap(['id', 'repeat']), file ] }
    ```

=== "Abans"

    ```groovy title="main.nf" linenums="8" hl_lines="1"
        getSampleIdAndReplicate = { meta, file -> [ meta.subMap(['id', 'repeat']), meta.subMap(['type']), file ] }
    ```

Executeu de nou per veure el resultat:

```bash
nextflow run main.nf
```

??? success "Sortida de la comanda"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [confident_leavitt] DSL2 - revision: a2303895bd

    [[id:patientA, repeat:1], patientA_rep1_normal.bam, patientA_rep1_tumor.bam]
    [[id:patientA, repeat:2], patientA_rep2_normal.bam, patientA_rep2_tumor.bam]
    [[id:patientB, repeat:1], patientB_rep1_normal.bam, patientB_rep1_tumor.bam]
    [[id:patientC, repeat:1], patientC_rep1_normal.bam, patientC_rep1_tumor.bam]
    ```

### Conclusió

En aquesta secció, heu après:

- **Manipular Tuples**: Com utilitzar `map` per aïllar un camp en una tupla
- **Unir Tuples**: Com utilitzar `join` per combinar tuples segons el primer camp
- **Crear Claus d'Unió**: Com utilitzar `subMap` per crear una nova clau d'unió
- **Closures amb Nom**: Com utilitzar una closure amb nom en map
- **Unió per Múltiples Camps**: Com unir per múltiples camps per a coincidències més precises
- **Optimització d'Estructura de Dades**: Com racionalitzar l'estructura del canal eliminant informació redundant

Ara teniu un workflow que pot dividir un full de mostres, filtrar les mostres normals i de tumor, unir-les per ID de mostra i número de rèplica, i després imprimir els resultats.

Aquest és un patró comú en workflows de bioinformàtica on necessiteu fer coincidir mostres o altres tipus de dades després de processar-les independentment, així que és una habilitat útil. A continuació, veurem com repetir una mostra múltiples vegades.

## 4. Distribuir mostres per intervals

Un patró clau en workflows de bioinformàtica és distribuir l'anàlisi a través de regions genòmiques. Per exemple, la crida de variants pot paral·lelitzar-se dividint el genoma en intervals (com cromosomes o regions més petites). Aquesta estratègia de paral·lelització millora significativament l'eficiència del pipeline distribuint la càrrega computacional a través de múltiples nuclis o nodes, reduint el temps d'execució global.

A la secció següent, demostrarem com distribuir les nostres dades de mostres a través de múltiples intervals genòmics. Aparellarem cada mostra amb cada interval, permetent el processament paral·lel de diferents regions genòmiques. Això multiplicarà la mida del nostre conjunt de dades pel nombre d'intervals, creant múltiples unitats d'anàlisi independents que poden reunir-se més tard.

### 4.1. Distribuir mostres per intervals utilitzant `combine`

Comencem creant un canal d'intervals. Per mantenir la vida simple, utilitzarem només 3 intervals que definirem manualment. En un workflow real, podríeu llegir-los d'una entrada de fitxer o fins i tot crear un canal amb molts fitxers d'intervals.

=== "Després"

    ```groovy title="main.nf" linenums="17" hl_lines="2"
            .join(ch_tumor_samples)
        ch_intervals = channel.of('chr1', 'chr2', 'chr3')
    ```

=== "Abans"

    ```groovy title="main.nf" linenums="17" hl_lines="2"
            .join(ch_tumor_samples)
        ch_joined_samples.view()
    ```

Ara recordeu, volem repetir cada mostra per a cada interval. Això de vegades es coneix com el producte cartesià de les mostres i els intervals. Podem aconseguir això utilitzant l'[operador `combine`](https://www.nextflow.io/docs/latest/operator.html#combine). Això prendrà cada element del canal 1 i el repetirà per a cada element del canal 2. Afegim un operador combine al nostre workflow:

=== "Després"

    ```groovy title="main.nf" linenums="18" hl_lines="3-5"
        ch_intervals = channel.of('chr1', 'chr2', 'chr3')

        ch_combined_samples = ch_joined_samples
            .combine(ch_intervals)
            .view()
    ```

=== "Abans"

    ```groovy title="main.nf" linenums="18"
        ch_intervals = channel.of('chr1', 'chr2', 'chr3')
    ```

Ara executem-ho i vegem què passa:

```bash
nextflow run main.nf
```

??? success "Sortida de la comanda"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [mighty_tesla] DSL2 - revision: ae013ab70b

    [[id:patientA, repeat:1], patientA_rep1_normal.bam, patientA_rep1_tumor.bam, chr1]
    [[id:patientA, repeat:1], patientA_rep1_normal.bam, patientA_rep1_tumor.bam, chr2]
    [[id:patientA, repeat:1], patientA_rep1_normal.bam, patientA_rep1_tumor.bam, chr3]
    [[id:patientA, repeat:2], patientA_rep2_normal.bam, patientA_rep2_tumor.bam, chr1]
    [[id:patientA, repeat:2], patientA_rep2_normal.bam, patientA_rep2_tumor.bam, chr2]
    [[id:patientA, repeat:2], patientA_rep2_normal.bam, patientA_rep2_tumor.bam, chr3]
    [[id:patientB, repeat:1], patientB_rep1_normal.bam, patientB_rep1_tumor.bam, chr1]
    [[id:patientB, repeat:1], patientB_rep1_normal.bam, patientB_rep1_tumor.bam, chr2]
    [[id:patientB, repeat:1], patientB_rep1_normal.bam, patientB_rep1_tumor.bam, chr3]
    [[id:patientC, repeat:1], patientC_rep1_normal.bam, patientC_rep1_tumor.bam, chr1]
    [[id:patientC, repeat:1], patientC_rep1_normal.bam, patientC_rep1_tumor.bam, chr2]
    [[id:patientC, repeat:1], patientC_rep1_normal.bam, patientC_rep1_tumor.bam, chr3]
    ```

Èxit! Hem repetit cada mostra per a cada interval de la nostra llista de 3 intervals. Hem triplicat efectivament el nombre d'elements al nostre canal.

És una mica difícil de llegir, però a la següent secció ho endreçarem.

### 4.2. Organitzar el canal

Podem utilitzar l'operador `map` per endreçar i refactoritzar les nostres dades de mostra perquè siguin més fàcils d'entendre. Movem la cadena d'intervals al mapa d'unió al primer element.

=== "Després"

    ```groovy title="main.nf" linenums="20" hl_lines="3-9"
        ch_combined_samples = ch_joined_samples
            .combine(ch_intervals)
            .map { grouping_key, normal, tumor, interval ->
                [
                    grouping_key + [interval: interval],
                    normal,
                    tumor
                ]
            }
            .view()
    ```

=== "Abans"

    ```groovy title="main.nf" linenums="20"
        ch_combined_samples = ch_joined_samples
            .combine(ch_intervals)
            .view()
    ```

Desglossem què fa aquesta operació map pas a pas.

Primer, utilitzem paràmetres amb nom per fer el codi més llegible. Utilitzant els noms `grouping_key`, `normal`, `tumor` i `interval`, podem referir-nos als elements de la tupla per nom en lloc de per índex:

```groovy
        .map { grouping_key, normal, tumor, interval ->
```

A continuació, combinem el `grouping_key` amb el camp `interval`. El `grouping_key` és un mapa que conté els camps `id` i `repeat`. Creem un nou mapa amb l'`interval` i els fusionem utilitzant l'addició de mapes de Groovy (`+`):

```groovy
                grouping_key + [interval: interval],
```

Finalment, retornem això com una tupla amb tres elements: el mapa de metadades combinat, el fitxer de mostra normal i el fitxer de mostra de tumor:

```groovy
            [
                grouping_key + [interval: interval],
                normal,
                tumor
            ]
```

Executem-ho de nou i comprovem els continguts del canal:

```bash
nextflow run main.nf
```

??? success "Sortida de la comanda"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [sad_hawking] DSL2 - revision: 1f6f6250cd

    [[id:patientA, repeat:1, interval:chr1], patientA_rep1_normal.bam, patientA_rep1_tumor.bam]
    [[id:patientA, repeat:1, interval:chr2], patientA_rep1_normal.bam, patientA_rep1_tumor.bam]
    [[id:patientA, repeat:1, interval:chr3], patientA_rep1_normal.bam, patientA_rep1_tumor.bam]
    [[id:patientA, repeat:2, interval:chr1], patientA_rep2_normal.bam, patientA_rep2_tumor.bam]
    [[id:patientA, repeat:2, interval:chr2], patientA_rep2_normal.bam, patientA_rep2_tumor.bam]
    [[id:patientA, repeat:2, interval:chr3], patientA_rep2_normal.bam, patientA_rep2_tumor.bam]
    [[id:patientB, repeat:1, interval:chr1], patientB_rep1_normal.bam, patientB_rep1_tumor.bam]
    [[id:patientB, repeat:1, interval:chr2], patientB_rep1_normal.bam, patientB_rep1_tumor.bam]
    [[id:patientB, repeat:1, interval:chr3], patientB_rep1_normal.bam, patientB_rep1_tumor.bam]
    [[id:patientC, repeat:1, interval:chr1], patientC_rep1_normal.bam, patientC_rep1_tumor.bam]
    [[id:patientC, repeat:1, interval:chr2], patientC_rep1_normal.bam, patientC_rep1_tumor.bam]
    [[id:patientC, repeat:1, interval:chr3], patientC_rep1_normal.bam, patientC_rep1_tumor.bam]
    ```

Utilitzar `map` per coaccionar les vostres dades a l'estructura correcta pot ser complicat, però és crucial per a la manipulació efectiva de dades.

Ara tenim cada mostra repetida a través de tots els intervals genòmics, creant múltiples unitats d'anàlisi independents que poden processar-se en paral·lel. Però què passa si volem reunir mostres relacionades? A la següent secció, aprendrem com agrupar mostres que comparteixen atributs comuns.

### Conclusió

En aquesta secció, heu après:

- **Distribuir mostres per intervals**: Com utilitzar `combine` per repetir mostres per intervals
- **Crear productes cartesians**: Com generar totes les combinacions de mostres i intervals
- **Organitzar l'estructura del canal**: Com utilitzar `map` per reestructurar dades per a millor llegibilitat
- **Preparació per a processament paral·lel**: Com configurar dades per a anàlisi distribuïda

## 5. Agregar mostres utilitzant `groupTuple`

A les seccions anteriors, hem après com dividir dades d'un fitxer d'entrada i filtrar per camps específics (en el nostre cas mostres normals i de tumor). Però això només cobreix un únic tipus d'unió. Què passa si volem agrupar mostres per un atribut específic? Per exemple, en lloc d'unir parelles aparellades normal-tumor, podríem voler processar totes les mostres de "sampleA" juntes independentment del seu tipus. Aquest patró és comú en workflows de bioinformàtica on podeu voler processar mostres relacionades per separat per raons d'eficiència abans de comparar o combinar els resultats al final.

Nextflow inclou mètodes integrats per fer això, el principal que veurem és `groupTuple`.

Comencem agrupant totes les nostres mostres que tenen els mateixos camps `id` i `interval`, això seria típic d'una anàlisi on volguéssim agrupar rèpliques tècniques però mantenir mostres significativament diferents separades.

Per fer això, hauríem de separar les nostres variables d'agrupació perquè puguem utilitzar-les de manera aïllada.

El primer pas és similar al que vam fer a la secció anterior. Hem d'aïllar la nostra variable d'agrupació com el primer element de la tupla. Recordeu, el nostre primer element és actualment un mapa dels camps `id`, `repeat` i `interval`:

```groovy title="main.nf" linenums="1"
{
  "id": "sampleA",
  "repeat": "1",
  "interval": "chr1"
}
```

Podem reutilitzar el mètode `subMap` d'abans per aïllar els nostres camps `id` i `interval` del mapa. Com abans, utilitzarem l'operador `map` per aplicar el mètode `subMap` al primer element de la tupla per a cada mostra.

=== "Després"

    ```groovy title="main.nf" linenums="20" hl_lines="11-19"
        ch_combined_samples = ch_joined_samples
            .combine(ch_intervals)
            .map { grouping_key, normal, tumor, interval ->
                [
                    grouping_key + [interval: interval],
                    normal,
                    tumor
                ]
            }

        ch_grouped_samples = ch_combined_samples
            .map { grouping_key, normal, tumor ->
                [
                    grouping_key.subMap('id', 'interval'),
                    normal,
                    tumor
                ]
              }
              .view()
    ```

=== "Abans"

    ```groovy title="main.nf" linenums="20" hl_lines="10"
        ch_combined_samples = ch_joined_samples
            .combine(ch_intervals)
            .map { grouping_key, normal, tumor, interval ->
                [
                    grouping_key + [interval: interval],
                    normal,
                    tumor
                ]
            }
            .view()
    ```

Executem-ho de nou i comprovem els continguts del canal:

```bash
nextflow run main.nf
```

??? success "Sortida de la comanda"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [hopeful_brenner] DSL2 - revision: 7f4f7fea76

    [[id:patientA, interval:chr1], patientA_rep1_normal.bam, patientA_rep1_tumor.bam]
    [[id:patientA, interval:chr2], patientA_rep1_normal.bam, patientA_rep1_tumor.bam]
    [[id:patientA, interval:chr3], patientA_rep1_normal.bam, patientA_rep1_tumor.bam]
    [[id:patientA, interval:chr1], patientA_rep2_normal.bam, patientA_rep2_tumor.bam]
    [[id:patientA, interval:chr2], patientA_rep2_normal.bam, patientA_rep2_tumor.bam]
    [[id:patientA, interval:chr3], patientA_rep2_normal.bam, patientA_rep2_tumor.bam]
    [[id:patientB, interval:chr1], patientB_rep1_normal.bam, patientB_rep1_tumor.bam]
    [[id:patientB, interval:chr2], patientB_rep1_normal.bam, patientB_rep1_tumor.bam]
    [[id:patientB, interval:chr3], patientB_rep1_normal.bam, patientB_rep1_tumor.bam]
    [[id:patientC, interval:chr1], patientC_rep1_normal.bam, patientC_rep1_tumor.bam]
    [[id:patientC, interval:chr2], patientC_rep1_normal.bam, patientC_rep1_tumor.bam]
    [[id:patientC, interval:chr3], patientC_rep1_normal.bam, patientC_rep1_tumor.bam]
    ```

Podem veure que hem aïllat amb èxit els camps `id` i `interval`, però encara no hem agrupat les mostres.

!!! note "Nota"

    Estem descartant el camp `replicate` aquí. Això és perquè no el necessitem per a processament posterior. Després de completar aquest tutorial, vegeu si podeu incloure'l sense afectar l'agrupació posterior!

Ara agrupem les mostres per aquest nou element d'agrupació, utilitzant l'[operador `groupTuple`](https://www.nextflow.io/docs/latest/operator.html#grouptuple).

=== "Després"

    ```groovy title="main.nf" linenums="30" hl_lines="9"
        ch_grouped_samples = ch_combined_samples
            .map { grouping_key, normal, tumor ->
                [
                    grouping_key.subMap('id', 'interval'),
                    normal,
                    tumor
                ]
              }
              .groupTuple()
              .view()
    ```

=== "Abans"

    ```groovy title="main.nf" linenums="30"
        ch_grouped_samples = ch_combined_samples
            .map { grouping_key, normal, tumor ->
                [
                    grouping_key.subMap('id', 'interval'),
                    normal,
                    tumor
                ]
              }
              .view()
    ```

Això és tot! Només hem afegit una única línia de codi. Vegem què passa quan l'executem:

```bash
nextflow run main.nf
```

??? success "Sortida de la comanda"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `main.nf` [friendly_jang] DSL2 - revision: a1bee1c55d

    [[id:patientA, interval:chr1], [patientA_rep1_normal.bam, patientA_rep2_normal.bam], [patientA_rep1_tumor.bam, patientA_rep2_tumor.bam]]
    [[id:patientA, interval:chr2], [patientA_rep1_normal.bam, patientA_rep2_normal.bam], [patientA_rep1_tumor.bam, patientA_rep2_tumor.bam]]
    [[id:patientA, interval:chr3], [patientA_rep1_normal.bam, patientA_rep2_normal.bam], [patientA_rep1_tumor.bam, patientA_rep2_tumor.bam]]
    [[id:patientB, interval:chr1], [patientB_rep1_normal.bam], [patientB_rep1_tumor.bam]]
    [[id:patientB, interval:chr2], [patientB_rep1_normal.bam], [patientB_rep1_tumor.bam]]
    [[id:patientB, interval:chr3], [patientB_rep1_normal.bam], [patientB_rep1_tumor.bam]]
    [[id:patientC, interval:chr1], [patientC_rep1_normal.bam], [patientC_rep1_tumor.bam]]
    [[id:patientC, interval:chr2], [patientC_rep1_normal.bam], [patientC_rep1_tumor.bam]]
    [[id:patientC, interval:chr3], [patientC_rep1_normal.bam], [patientC_rep1_tumor.bam]]
    ```

Noteu que les nostres dades han canviat d'estructura i dins de cada element del canal els fitxers ara estan continguts en tuples com `[patientA_rep1_normal.bam, patientA_rep2_normal.bam]`. Això és perquè quan utilitzem `groupTuple`, Nextflow combina els fitxers individuals per a cada mostra d'un grup. Això és important recordar quan s'intenta gestionar les dades posteriorment.

!!! note "Nota"

    [`transpose`](https://www.nextflow.io/docs/latest/reference/operator.html#transpose) és l'oposat de groupTuple. Desempaqueta els elements d'un canal i els aplana. Proveu d'afegir `transpose` i desfer l'agrupació que hem realitzat anteriorment!

### Conclusió

En aquesta secció, heu après:

- **Agrupar mostres relacionades**: Com utilitzar `groupTuple` per agregar mostres per atributs comuns
- **Aïllar claus d'agrupació**: Com utilitzar `subMap` per extreure camps específics per a l'agrupació
- **Gestionar estructures de dades agrupades**: Com treballar amb l'estructura niuada creada per `groupTuple`
- **Gestió de rèpliques tècniques**: Com agrupar mostres que comparteixen les mateixes condicions experimentals

---

## Resum

En aquesta missió secundària, heu après com dividir i agrupar dades utilitzant canals.

Modificant les dades a mesura que flueixen pel pipeline, podeu construir un pipeline escalable sense utilitzar bucles o sentències while, oferint diversos avantatges sobre enfocaments més tradicionals:

- Podem escalar a tantes o tan poques entrades com vulguem sense codi addicional
- Ens centrem en gestionar el flux de dades a través del pipeline, en lloc de la iteració
- Podem ser tan complexos o simples com sigui necessari
- El pipeline es torna més declaratiu, centrant-se en què hauria de passar en lloc de com hauria de passar
- Nextflow optimitzarà l'execució per a nosaltres executant operacions independents en paral·lel

Dominar aquestes operacions de canal us permetrà construir pipelines flexibles i escalables que gestionen relacions de dades complexes sense recórrer a bucles o programació iterativa, permetent a Nextflow optimitzar l'execució i paral·lelitzar operacions independents automàticament.

### Patrons clau

1.  **Crear dades d'entrada estructurades:** Començant des d'un fitxer CSV amb mapes de metadades (basant-se en patrons de [Metadades en workflows](./metadata.md))

    ```groovy
    ch_samples = channel.fromPath("./data/samplesheet.csv")
        .splitCsv(header: true)
        .map{ row ->
          [[id:row.id, repeat:row.repeat, type:row.type], row.bam]
        }
    ```

2.  **Dividir dades en canals separats:** Hem utilitzat `filter` per dividir dades en fluxos independents segons el camp `type`

    ```groovy
    channel.filter { it.type == 'tumor' }
    ```

3.  **Unir mostres aparellades:** Hem utilitzat `join` per recombinar mostres relacionades segons els camps `id` i `repeat`

    - Unir dos canals per clau (primer element de la tupla)

    ```groovy
    tumor_ch.join(normal_ch)
    ```

    - Extreure clau d'unió i unir per aquest valor

    ```groovy
    tumor_ch.map { meta, file -> [meta.id, meta, file] }
        .join(
          normal_ch.map { meta, file -> [meta.id, meta, file] }
        )
    ```

    - Unir per múltiples camps utilitzant subMap

    ```groovy
    tumor_ch.map { meta, file -> [meta.subMap(['id', 'repeat']), meta, file] }
        .join(
          normal_ch.map { meta, file -> [meta.subMap(['id', 'repeat']), meta, file] }
        )
    ```

4.  **Distribuir a través d'intervals:** Hem utilitzat `combine` per crear productes cartesians de mostres amb intervals genòmics per a processament paral·lel.

    ```groovy
    samples_ch.combine(intervals_ch)
    ```

5.  **Agregar per claus d'agrupació:** Hem utilitzat `groupTuple` per agrupar pel primer element de cada tupla, recollint així mostres que comparteixen els camps `id` i `interval` i fusionant rèpliques tècniques.

    ```groovy
    channel.groupTuple()
    ```

6.  **Optimitzar l'estructura de dades:** Hem utilitzat `subMap` per extreure camps específics i hem creat una closure amb nom per fer transformacions reutilitzables.

    - Extreure camps específics d'un mapa

    ```groovy
    meta.subMap(['id', 'repeat'])
    ```

    - Utilitzar closure amb nom per a transformacions reutilitzables

    ```groovy
    getSampleIdAndReplicate = { meta, file -> [meta.subMap(['id', 'repeat']), file] }
    channel.map(getSampleIdAndReplicate)
    ```

### Recursos addicionals

- [filter](https://www.nextflow.io/docs/latest/operator.html#filter)
- [map](https://www.nextflow.io/docs/latest/operator.html#map)
- [join](https://www.nextflow.io/docs/latest/operator.html#join)
- [groupTuple](https://www.nextflow.io/docs/latest/operator.html#grouptuple)
- [combine](https://www.nextflow.io/docs/latest/operator.html#combine)

---

## Què segueix?

Torneu al [menú de Missions Secundàries](./index.md) o feu clic al botó a la part inferior dreta de la pàgina per passar al següent tema de la llista.
