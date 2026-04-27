---

# Entorn de Desenvolupament

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducció assistida per IA - [més informació i suggeriments](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Els Entorns de Desenvolupament Integrats (IDE) moderns poden transformar radicalment la vostra experiència de desenvolupament amb Nextflow. Aquesta missió secundària se centra específicament en aprofitar VS Code i la seva extensió de Nextflow per escriure codi més ràpidament, detectar errors d'hora i navegar per workflows complexos de manera eficient.

!!! note "Això no és un tutorial tradicional"

    A diferència d'altres mòduls de formació, aquesta guia està organitzada com una col·lecció de consells ràpids, suggeriments i exemples pràctics en lloc d'un tutorial pas a pas. Cada secció es pot explorar de manera independent segons els vostres interessos i les vostres necessitats de desenvolupament actuals. Sentiu-vos lliures de saltar entre seccions i centrar-vos en les funcionalitats que siguin més útils immediatament per al vostre desenvolupament de workflows.

## Què hauríeu de saber primer

Aquesta guia assumeix que heu completat el curs de formació [Hello Nextflow](../hello_nextflow/) i que esteu familiaritzats amb els conceptes fonamentals de Nextflow, incloent-hi:

- **Estructura bàsica del workflow**: Comprensió dels processos, workflows i com es connecten entre si
- **Operacions amb canals**: Creació de canals, pas de dades entre processos i ús d'operadors bàsics
- **Mòduls i organització**: Creació de mòduls reutilitzables i ús de sentències include
- **Conceptes bàsics de configuració**: Ús de `nextflow.config` per a paràmetres, directives de procés i perfils

## Què aprendreu aquí

Aquesta guia se centra en les **funcionalitats de productivitat de l'IDE** que us convertiran en un desenvolupador de Nextflow més eficient:

- **Ressaltat de sintaxi avançat**: Comprensió del que VS Code us mostra sobre l'estructura del vostre codi
- **Autocompleció intel·ligent**: Aprofitament de suggeriments contextuals per escriure codi més ràpidament
- **Detecció d'errors i diagnòstics**: Detecció d'errors de sintaxi abans d'executar el vostre workflow
- **Navegació pel codi**: Moviment ràpid entre processos, mòduls i definicions
- **Formatació i organització**: Manteniment d'un estil de codi consistent i llegible
- **Desenvolupament assistit per IA** (opcional): Ús d'eines d'IA modernes integrades amb el vostre IDE

!!! info "Per què les funcionalitats de l'IDE ara?"

    Probablement ja heu estat utilitzant VS Code durant el curs [Hello Nextflow](../hello_nextflow/), però hem mantingut el focus en aprendre els fonaments de Nextflow en lloc de les funcionalitats de l'IDE. Ara que esteu familiaritzats amb els conceptes bàsics de Nextflow com processos, workflows, canals i mòduls, esteu preparats per aprofitar les sofisticades funcionalitats de l'IDE que us convertiran en un desenvolupador més eficient.

    Penseu en això com a "pujar de nivell" el vostre entorn de desenvolupament: el mateix editor que heu estat utilitzant té capacitats molt més potents que es tornen veritablement valuoses un cop enteneu per a què us estan ajudant.

---

## 0. Configuració i Escalfament

Configurem un espai de treball específicament per explorar les funcionalitats de l'IDE:

```bash title="Navigate to the IDE features directory"
cd side-quests/ide_features
```

Obriu aquest directori a VS Code:

```bash title="Open VS Code in current directory"
code .
```

El directori `ide_features` conté workflows d'exemple que demostren diverses funcionalitats de l'IDE:

```bash title="Show directory structure"
tree .
```

```console title="Project structure"
tree .
.
├── basic_workflow.nf
├── complex_workflow.nf
├── data
│   ├── sample_001.fastq.gz
│   ├── sample_002.fastq.gz
│   ├── sample_003.fastq.gz
│   ├── sample_004.fastq.gz
│   ├── sample_005.fastq.gz
│   └── sample_data.csv
├── modules
│   ├── fastqc.nf
│   ├── star.nf
│   └── utils.nf
└── nextflow.config

3 directories, 12 files
```

!!! note "Sobre els fitxers d'exemple"

    - `basic_workflow.nf` és un workflow bàsic funcional que podeu executar i modificar
    - `complex_workflow.nf` està dissenyat únicament per a il·lustració per demostrar les funcionalitats de navegació; pot ser que no s'executi correctament, però mostra una estructura de workflow multifitxer realista

### Dreceres de teclat

Algunes de les funcionalitats d'aquesta guia utilitzen dreceres de teclat opcionals. És possible que accediu a aquest material a través de GitHub Codespaces al navegador, i en aquest cas de vegades les dreceres no funcionaran com s'espera perquè s'utilitzen per a altres coses al vostre sistema.

Si executeu VS Code localment, com probablement fareu quan estigueu escrivint workflows, les dreceres funcionaran tal com es descriu.

Si utilitzeu un Mac, algunes (no totes) les dreceres de teclat utilitzaran "cmd" en lloc de "ctrl", i ho indicarem al text com `Ctrl/Cmd`.

### 0.1. Instal·lació de l'extensió de Nextflow

!!! note "Ja utilitzeu Devcontainers?"

    Si esteu treballant a **GitHub Codespaces** o utilitzant un **devcontainer local**, l'extensió de Nextflow probablement ja està instal·lada i configurada per a vosaltres. Podeu ometre els passos d'instal·lació manual que hi ha a continuació i procedir directament a explorar les funcionalitats de l'extensió.

Per instal·lar l'extensió manualment:

1. Obriu VS Code
2. Aneu a la vista d'Extensions fent clic a la icona d'extensions a l'esquerra: ![icona d'extensions](img/extensions_icon.png) (drecera `Ctrl/Cmd+Shift+X` si executeu VSCode localment)
3. Cerqueu "Nextflow"
4. Instal·leu l'extensió oficial de Nextflow

![Instal·la l'extensió de Nextflow](img/install_extension.png)

### 0.2. Disposició de l'espai de treball

Com que heu estat utilitzant VS Code durant tot Hello Nextflow, ja esteu familiaritzats amb els conceptes bàsics. Aquí teniu com organitzar el vostre espai de treball de manera eficient per a aquesta sessió:

- **Àrea de l'editor**: Per visualitzar i editar fitxers. Podeu dividir-la en múltiples panells per comparar fitxers un al costat de l'altre.
- **Explorador de fitxers** feu clic a (![icona de l'explorador de fitxers](img/files_icon.png)) (`Ctrl/Cmd+Shift+E`): Els fitxers i carpetes locals del vostre sistema. Manteniu-lo obert a l'esquerra per navegar entre fitxers.
- **Terminal integrat** (`Ctrl+Shift+` accent greu tant per a Windows com per a MacOS): Un terminal per interactuar amb l'ordinador a la part inferior. Utilitzeu-lo per executar Nextflow o altres comandes.
- **Panell de problemes** (`Ctrl+Shift+M`): VS Code mostrarà aquí qualsevol error i problema que detecti. Això és útil per destacar problemes d'un cop d'ull.

Podeu arrossegar els panells o amagar-los (`Ctrl/Cmd+B` per alternar la barra lateral) per personalitzar la vostra disposició mentre treballem pels exemples.

### Conclusió

Teniu VS Code configurat amb l'extensió de Nextflow i enteneu la disposició de l'espai de treball per a un desenvolupament eficient.

### Què segueix?

Apreneu com el ressaltat de sintaxi us ajuda a entendre l'estructura del codi Nextflow d'un cop d'ull.

---

## 1. Ressaltat de Sintaxi i Estructura del Codi

Ara que el vostre espai de treball està configurat, explorem com el ressaltat de sintaxi de VS Code us ajuda a llegir i escriure codi Nextflow de manera més efectiva.

### 1.1. Elements de sintaxi de Nextflow

Obriu `basic_workflow.nf` per veure el ressaltat de sintaxi en acció:

![Mostra de sintaxi](img/syntax_showcase.png)

Observeu com VS Code ressalta:

- **Paraules clau** (`process`, `workflow`, `input`, `output`, `script`) amb colors distintius
- **Literals de string** i **paràmetres** amb estils diferents
- **Comentaris** amb un color apagat
- **Variables** i **crides a funcions** amb l'èmfasi adequat
- **Blocs de codi** amb guies d'indentació adequades

!!! note "Colors depenents del tema"

    Els colors específics que veureu dependran del vostre tema de VS Code (mode fosc/clar), la configuració de colors i qualsevol personalització que hàgiu fet. El important és que els diferents elements de sintaxi es distingeixen visualment entre si, fent que l'estructura del codi sigui més fàcil d'entendre independentment de l'esquema de colors que hàgiu triat.

### 1.2. Comprensió de l'estructura del codi

El ressaltat de sintaxi us ajuda a identificar ràpidament:

- **Límits dels processos**: Distinció clara entre diferents processos
- **Blocs d'entrada/sortida**: Fàcil de localitzar les definicions del flux de dades
- **Blocs de script**: Les comandes reals que s'estan executant
- **Operacions amb canals**: Passos de transformació de dades
- **Directives de configuració**: Configuració específica del procés

Aquesta organització visual es torna inestimable quan es treballa amb workflows complexos que contenen múltiples processos i fluxos de dades intricats.

### Conclusió

Enteneu com el ressaltat de sintaxi de VS Code us ajuda a llegir l'estructura del codi Nextflow i identificar els diferents elements del llenguatge per a un desenvolupament més ràpid.

### Què segueix?

Apreneu com l'autocompleció intel·ligent accelera l'escriptura de codi amb suggeriments contextuals.

---

## 2. Autocompleció Intel·ligent

Les funcionalitats d'autocompleció de VS Code us ajuden a escriure codi més ràpidament i amb menys errors suggerint opcions adequades basades en el context.

### 2.1. Suggeriments contextuals

Les opcions d'autocompleció varien depenent d'on us trobeu al vostre codi:

#### Operacions amb canals

Obriu `basic_workflow.nf` de nou i proveu d'escriure `channel.` al bloc del workflow:

![Autocompleció de canal](img/autocomplete_channel.png)

Veureu suggeriments per a:

- `fromPath()` - Crea un canal a partir de rutes de fitxers
- `fromFilePairs()` - Crea un canal a partir de fitxers aparellats
- `of()` - Crea un canal a partir de valors
- `fromSRA()` - Crea un canal a partir d'accessions SRA
- I molts més...

Això us ajuda a trobar ràpidament la fàbrica de canals adequada sense necessitat de recordar els noms exactes dels mètodes.

També podeu descobrir els operadors disponibles per aplicar als canals. Per exemple, escriviu `FASTQC.out.html.` per veure les operacions disponibles:

![Autocompleció d'operadors de canal](img/autocomplete_operators.png)

#### Directives de procés

Dins d'un bloc script d'un procés, escriviu `task.` per veure les propietats d'execució disponibles:

![Autocompleció de propietats de tasca](img/autocomplete_task.png)

#### Configuració

Obriu nextflow.config i escriviu `process.` en qualsevol lloc per veure les directives de procés disponibles:

![Autocompleció de configuració](img/autocomplete_config.png)

Veureu suggeriments per a:

- `executor`
- `memory`
- `cpus`

Això estalvia temps en configurar processos i funciona en diferents àmbits de configuració. Per exemple, proveu d'escriure `docker.` per veure les opcions de configuració específiques de Docker.

### Conclusió

Podeu utilitzar l'autocompleció intel·ligent de VS Code per descobrir les operacions de canal disponibles, les directives de procés i les opcions de configuració sense memoritzar la sintaxi.

### Què segueix?

Apreneu com la detecció d'errors en temps real us ajuda a detectar problemes abans d'executar el vostre workflow, simplement llegint el codi.

## 3. Detecció d'Errors i Diagnòstics

La detecció d'errors en temps real de VS Code us ajuda a detectar problemes abans d'executar el vostre workflow.

### 3.1. Detecció d'errors de sintaxi

Creem un error deliberat per veure la detecció en acció. Obriu `basic_workflow.nf` i canvieu el nom del procés de `FASTQC` a `FASTQ` (o qualsevol altre nom no vàlid). VS Code ressaltarà immediatament l'error al bloc del workflow amb un subratllat vermell ondulat:

![Subratllat d'error](img/error_underline.png)

### 3.2. Panell de problemes

Més enllà del ressaltat individual d'errors, VS Code proporciona un panell de Problemes centralitzat que agrega tots els errors, advertències i missatges d'informació de tot el vostre espai de treball. Obriu-lo amb `Ctrl/Cmd+Shift+M` i utilitzeu la icona de filtre per mostrar només els errors rellevants per al fitxer actual:

![Filtra el panell de problemes](img/active_file.png)

Feu clic en qualsevol problema per saltar directament a la línia problemàtica

![Panell de problemes](img/problems_panel.png)

Corregiu l'error canviant el nom del procés de nou a `FASTQC`.

### 3.3. Patrons d'errors comuns

Els errors comuns en la sintaxi de Nextflow inclouen:

- **Claudàtors que falten**: `{` o `}` sense parella
- **Blocs incomplets**: Seccions requerides que falten als processos
- **Sintaxi no vàlida**: DSL de Nextflow mal format
- **Errors tipogràfics en paraules clau**: Directives de procés mal escrites
- **Incompatibilitats de canals**: Incompatibilitats de tipus

El servidor de llenguatge de Nextflow ressalta aquests problemes al panell de Problemes. Podeu revisar-los d'hora per evitar errors de sintaxi mentre executeu un pipeline.

### Conclusió

Podeu utilitzar la detecció d'errors de VS Code i el panell de Problemes per detectar errors de sintaxi i problemes abans d'executar el vostre workflow, estalviant temps i evitant frustracions.

### Què segueix?

Apreneu a navegar eficientment entre processos, mòduls i definicions en workflows complexos.

---

## 4. Navegació pel Codi i Gestió de Símbols

La navegació eficient és crucial quan es treballa amb workflows complexos que abasten múltiples fitxers. Per entendre això, substituïu la definició del procés a `basic_workflow.nf` per una importació del mòdul que us hem proporcionat:

=== "Després"

    ```groovy title="basic_workflow.nf" linenums="3"
    include { FASTQC } from './modules/fastqc.nf'
    ```

=== "Abans"

    ```groovy title="basic_workflow.nf" linenums="3"
    process FASTQC {
        tag "${sample_id}"
        publishDir "${params.output_dir}/fastqc", mode: 'copy'

        input:
        tuple val(sample_id), path(reads)

        output:
        tuple val(sample_id), path("*.html"), emit: html
        tuple val(sample_id), path("*.zip"), emit: zip

        script:
        def args = task.ext.args ?: ''
        """
        fastqc \\
            ${args} \\
            --threads ${task.cpus} \\
            ${reads}
        """
    }
    ```

### 4.1. Anar a la definició

Si passeu el ratolí per sobre del nom d'un procés com `FASTQC`, veureu una finestra emergent amb la interfície del mòdul (entrades i sortides):

![Anar a la definició](img/syntax.png)

Aquesta funcionalitat és particularment valuosa quan s'escriuen workflows, ja que us permet entendre la interfície del mòdul sense obrir el fitxer del mòdul directament.

Podeu navegar ràpidament a qualsevol definició de procés, mòdul o variable utilitzant **Ctrl/Cmd-clic**. Passeu el ratolí per sobre de l'enllaç al fitxer del mòdul a la part superior del script i seguiu l'enllaç tal com se suggereix:

![Segueix l'enllaç](img/follow_link.png)

El mateix funciona per als noms de processos. Torneu a `basic_workflow.nf` i proveu-ho amb el nom del procés `FASTQC` al bloc del workflow. Això us porta directament al nom del procés (que és el mateix que el fitxer del mòdul en aquest exemple, però podria estar a mig camí d'un fitxer molt més gran).

Per tornar on éreu, utilitzeu **Alt+←** (o **Ctrl+-** al Mac). Aquesta és una manera potent d'explorar el codi sense perdre el vostre lloc.

Ara explorem la navegació en un workflow més complex utilitzant `complex_workflow.nf` (el fitxer només per a il·lustració esmentat anteriorment). Aquest workflow conté múltiples processos definits en fitxers de mòdul separats, així com alguns en línia. Tot i que les estructures multifitxer complexes poden ser difícils de navegar manualment, la capacitat de saltar a les definicions fa que l'exploració sigui molt més manejable.

1. Obriu `complex_workflow.nf`
2. Navegueu a les definicions dels mòduls
3. Utilitzeu **Alt+←** (o **Ctrl+-**) per navegar enrere
4. Navegueu al nom del procés `FASTQC` al bloc del workflow. Això us porta directament al nom del procés (que és el mateix que el fitxer del mòdul en aquest exemple, però podria estar a mig camí d'un fitxer molt més gran).
5. Navegueu enrere de nou
6. Navegueu al procés `TRIM_GALORE` al bloc del workflow. Està definit en línia, de manera que no us portarà a un fitxer separat, però encara us mostrarà la definició del procés i podreu navegar enrere on éreu.

### 4.2. Navegació per símbols

Amb `complex_workflow.nf` encara obert, podeu obtenir una visió general de tots els símbols del fitxer escrivint `@` a la barra de cerca a la part superior de VSCode (la drecera de teclat és `Ctrl/Cmd+Shift+O`, però pot no funcionar a Codespaces). Això obre el panell de navegació per símbols, que llista tots els símbols del fitxer actual:

![Navegació per símbols](img/symbols.png)

Això mostra:

- Totes les definicions de processos
- Definicions de workflows (hi ha dos workflows definits en aquest fitxer)
- Definicions de funcions

Comenceu a escriure per filtrar els resultats.

### 4.3. Trobar totes les referències

Entendre on s'utilitza un procés o variable a tot el vostre codi pot ser molt útil. Per exemple, si voleu trobar totes les referències al procés `FASTQC`, comenceu navegant a la seva definició. Podeu fer-ho obrint `modules/fastqc.nf` directament, o utilitzant la funcionalitat de navegació ràpida de VS Code amb `Ctrl/Cmd-clic` com hem fet anteriorment. Un cop a la definició del procés, feu clic dret sobre el nom del procés `FASTQC` i seleccioneu "Find All References" del menú contextual per veure totes les instàncies on s'utilitza.

![Trobar referències](img/references.png)

Aquesta funcionalitat mostra totes les instàncies on es fa referència a `FASTQC` dins del vostre espai de treball, incloent-hi el seu ús als dos workflows diferents. Aquesta informació és crucial per avaluar l'impacte potencial de les modificacions al procés `FASTQC`.

### 4.4. Panell d'esquema

El panell d'Esquema, situat a la barra lateral de l'Explorador (feu clic a ![Icona de l'Explorador](img/files_icon.png)), proporciona una visió general convenient de tots els símbols del vostre fitxer actual. Aquesta funcionalitat us permet navegar ràpidament i gestionar l'estructura del vostre codi mostrant funcions, variables i altres elements clau en una vista jeràrquica.

![Panell d'esquema](img/outline.png)

Utilitzeu el panell d'Esquema per navegar ràpidament a diferents parts del vostre codi sense utilitzar el navegador de fitxers.

### 4.5. Visualització del DAG

L'extensió de Nextflow de VS Code pot visualitzar el vostre workflow com un Graf Acíclic Dirigit (DAG). Això us ajuda a entendre el flux de dades i les dependències entre processos. Obriu `complex_workflow.nf` i feu clic al botó "Preview DAG" que apareix sobre `workflow {` (el segon bloc `workflow` d'aquest fitxer):

![Previsualització del DAG](img/dag_preview.png)

Això és només el workflow d'entrada, però també podeu previsualitzar el DAG per als workflows interns fent clic al botó "Preview DAG" sobre el workflow `RNASEQ_PIPELINE {` més amunt:

![Previsualització del DAG del workflow intern](img/dag_preview_inner.png)

Per a aquest workflow, podeu utilitzar els nodes del DAG per navegar a les definicions de processos corresponents al codi. Feu clic en un node i us portarà a la definició del procés rellevant a l'editor. Especialment quan un workflow creix fins a una mida gran, això pot ajudar-vos realment a navegar pel codi i entendre com estan connectats els processos.

### Conclusió

Podeu navegar per workflows complexos de manera eficient utilitzant anar a la definició, cerca de símbols, trobar referències i visualització del DAG per entendre l'estructura del codi i les dependències.

### Què segueix?

Apreneu a treballar eficaçment en múltiples fitxers interconnectats en projectes Nextflow més grans.

## 5. Treballar amb Múltiples Fitxers

El desenvolupament real amb Nextflow implica treballar amb múltiples fitxers interconnectats. Explorem com VS Code us ajuda a gestionar projectes complexos de manera eficient.

### 5.1. Navegació ràpida entre fitxers

Amb `complex_workflow.nf` obert, observareu que importa diversos mòduls. Practiquem la navegació ràpida entre ells.

Premeu **Ctrl+P** (o **Cmd+P**) i comenceu a escriure "fast":

VS Code us mostrarà els fitxers coincidents. Seleccioneu `modules/fastqc.nf` per saltar-hi instantàniament. Això és molt més ràpid que fer clic a través de l'explorador de fitxers quan sabeu aproximadament quin fitxer esteu buscant.

Proveu-ho amb altres patrons:

- Escriviu "star" per trobar el fitxer del mòdul d'alineament STAR (`star.nf`)
- Escriviu "utils" per trobar el fitxer de funcions utilitàries (`utils.nf`)
- Escriviu "config" per saltar als fitxers de configuració (`nextflow.config`)

### 5.2. Editor dividit per al desenvolupament multifitxer

Quan treballeu amb mòduls, sovint necessiteu veure tant el workflow principal com les definicions dels mòduls simultàniament. Configurem-ho:

1. Obriu `complex_workflow.nf`
2. Obriu `modules/fastqc.nf` en una pestanya nova
3. Feu clic dret a la pestanya `modules/fastqc.nf` i seleccioneu "Split Right"
4. Ara podeu veure els dos fitxers un al costat de l'altre

![Editor dividit](img/split_editor.png)

Això és inestimable quan:

- Comproveu les interfícies dels mòduls mentre escriviu crides al workflow, i la previsualització no és suficient
- Compareu processos similars en diferents mòduls
- Depureu el flux de dades entre el workflow i els mòduls

### 5.3. Cerca a tot el projecte

De vegades necessiteu trobar on s'utilitzen patrons específics a tot el vostre projecte. Premeu `Ctrl/Cmd+Shift+F` per obrir el panell de cerca.

Proveu de cercar `publishDir` a tot l'espai de treball:

![Cerca al projecte](img/project_search.png)

Això us mostra tots els fitxers que utilitzen directoris de publicació, ajudant-vos a:

- Entendre els patrons d'organització de sortides
- Trobar exemples de directives específiques
- Garantir la consistència entre mòduls

### Conclusió

Podeu gestionar projectes complexos multifitxer utilitzant la navegació ràpida entre fitxers, editors dividits i la cerca a tot el projecte per treballar eficientment entre workflows i mòduls.

### Què segueix?

Apreneu com les funcionalitats de formatació i manteniment del codi mantenen els vostres workflows organitzats i llegibles.

---

## 6. Formatació i Manteniment del Codi

La formatació adequada del codi és essencial no només per a l'estètica, sinó també per millorar la llegibilitat, la comprensió i la facilitat d'actualitzar workflows complexos.

### 6.1. Formatació automàtica en acció

Obriu `basic_workflow.nf` i desordeneu deliberadament la formatació:

- Elimineu part de la indentació: Seleccioneu tot el document i premeu `shift+tab` moltes vegades per eliminar tantes indentacions com sigui possible.
- Afegiu espais addicionals en llocs aleatoris: a la sentència `channel.fromPath`, afegiu 30 espais després del `(`.
- Trenqueu algunes línies de manera incòmoda: Afegiu una nova línia entre l'operador `.view {` i la cadena `Processing sample:` però no afegiu un salt de línia corresponent abans del parèntesi de tancament `}`.

Ara premeu `Shift+Alt+F` (o `Shift+Option+F` al MacOS) per formatar automàticament:

VS Code immediatament:

- Corregeix la indentació per mostrar l'estructura del procés clarament
- Alinea elements similars de manera consistent
- Elimina els espais en blanc innecessaris
- Manté salts de línia llegibles

Tingueu en compte que la formatació automàtica pot no resoldre tots els problemes d'estil del codi. El servidor de llenguatge de Nextflow intenta mantenir el vostre codi ordenat, però també respecta les vostres preferències personals en certes àrees. Per exemple, si elimineu la indentació dins del bloc `script` d'un procés, el formatador ho deixarà tal com està, ja que potser preferiu intencionadament aquest estil.

Actualment, no hi ha una aplicació estricta d'estil per a Nextflow, de manera que el servidor de llenguatge ofereix una certa flexibilitat. No obstant això, aplicarà de manera consistent les regles de formatació al voltant de les definicions de mètodes i funcions per mantenir la claredat.

### 6.2. Funcionalitats d'organització del codi

#### Comentaris ràpids

Seleccioneu un bloc de codi al vostre workflow i premeu **Ctrl+/** (o **Cmd+/**) per comentar-lo:

```groovy
// workflow {
//     ch_input = channel.fromPath(params.input)
//         .splitCsv(header: true)
//         .map { row -> [row.sample_id, file(row.fastq_path)] }
//
//     FASTQC(ch_input)
// }
```

Això és perfecte per a:

- Desactivar temporalment parts dels workflows durant el desenvolupament
- Afegir comentaris explicatius a operacions de canal complexes
- Documentar seccions del workflow

Utilitzeu **Ctrl+/** (o **Cmd+/**) de nou per descomentar el codi.

#### Plegament de codi per a una visió general

A `complex_workflow.nf`, observeu les petites fletxes al costat de les definicions de processos. Feu-hi clic per plegar (col·lapsar) els processos:

![Plegament de codi](img/code_folding.png)

Això us dóna una visió general d'alt nivell de l'estructura del vostre workflow sense perdre-us en els detalls d'implementació.

#### Coincidència de claudàtors

Col·loqueu el cursor al costat de qualsevol claudàtor `{` o `}` i VS Code ressaltarà el claudàtor coincident. Utilitzeu **Ctrl+Shift+\\** (o **Cmd+Shift+\\**) per saltar entre claudàtors coincidents.

Això és crucial per a:

- Entendre els límits dels processos
- Trobar claudàtors que falten o sobren
- Navegar per estructures de workflow niades

#### Selecció i edició multilínia

Per editar múltiples línies simultàniament, VS Code ofereix potents capacitats de multicursor:

- **Selecció multilínia**: Manteniu premut **Ctrl+Alt** (o **Cmd+Option** al MacOS) i utilitzeu les tecles de fletxa per seleccionar múltiples línies
- **Indentació multilínia**: Seleccioneu múltiples línies i utilitzeu **Tab** per indentar o **Shift+Tab** per desdentar blocs sencers

Això és particularment útil per a:

- Indentar blocs de processos sencers de manera consistent
- Afegir comentaris a múltiples línies alhora
- Editar definicions de paràmetres similars en múltiples processos

### Conclusió

Podeu mantenir un codi net i llegible utilitzant la formatació automàtica, les funcionalitats de comentaris, el plegament de codi, la coincidència de claudàtors i l'edició multilínia per organitzar workflows complexos de manera eficient.

### Què segueix?

Apreneu com VS Code s'integra amb el vostre flux de treball de desenvolupament més ampli més enllà de simplement editar codi.

---

## 7. Integració del Flux de Treball de Desenvolupament

VS Code s'integra bé amb el vostre flux de treball de desenvolupament més enllà de simplement editar codi.

### 7.1. Integració amb el control de versions

!!! note "Codespaces i integració amb Git"

    Si esteu treballant a **GitHub Codespaces**, algunes funcionalitats d'integració amb Git poden no funcionar com s'espera, especialment les dreceres de teclat per al Control de Codi Font. És possible que també hàgiu declinat obrir el directori com a repositori Git durant la configuració inicial, cosa que és correcta per a propòsits de formació.

Si el vostre projecte és un repositori git (com ho és aquest), VS Code mostra:

- Fitxers modificats amb indicadors de color
- Estat de Git a la barra d'estat
- Vistes de diferències en línia
- Capacitats de commit i push

Obriu el panell de Control de Codi Font utilitzant el botó de control de codi font (![Icona de control de codi font](img/source_control_icon.png)) (`Ctrl+Shift+G` o `Cmd+Shift+G` si esteu treballant amb VSCode localment) per veure els canvis de git i fer commits directament a l'editor.

![Panell de Control de Codi Font](img/source_control.png)

### 7.2. Execució i inspecció de workflows

Executem un workflow i després inspeccionem els resultats. Al terminal integrat (`Ctrl+Shift+` accent greu tant per a Windows com per a MacOS), executeu el workflow bàsic:

```bash title="Run the basic workflow"
nextflow run basic_workflow.nf --input data/sample_data.csv --output_dir results
```

Mentre s'executa el workflow, veureu la sortida en temps real al terminal. Després de la finalització, podeu utilitzar VS Code per inspeccionar els resultats sense sortir de l'editor:

1. **Navegueu als directoris de treball**: Utilitzeu l'explorador de fitxers o el terminal per navegar per `.nextflow/work`
2. **Obriu fitxers de registre**: Feu clic a les rutes dels fitxers de registre a la sortida del terminal per obrir-los directament a VS Code
3. **Inspeccioneu les sortides**: Navegueu pels directoris de resultats publicats a l'explorador de fitxers
4. **Visualitzeu informes d'execució**: Obriu informes HTML directament a VS Code o al vostre navegador

Això manté tot en un sol lloc en lloc de canviar entre múltiples aplicacions.

### Conclusió

Podeu integrar VS Code amb el control de versions i l'execució de workflows per gestionar tot el vostre procés de desenvolupament des d'una única interfície.

### Què segueix?

Veieu com totes aquestes funcionalitats de l'IDE funcionen juntes en el vostre flux de treball de desenvolupament diari.

---

## 8. Resum i notes ràpides

Aquí teniu algunes notes ràpides sobre cadascuna de les funcionalitats de l'IDE comentades anteriorment:

### 8.1. Iniciar una nova funcionalitat

1. **Obertura ràpida de fitxers** (`Ctrl+P` o `Cmd+P`) per trobar mòduls existents rellevants
2. **Editor dividit** per veure processos similars un al costat de l'altre
3. **Navegació per símbols** (`Ctrl+Shift+O` o `Cmd+Shift+O`) per entendre l'estructura del fitxer
4. **Autocompleció** per escriure codi nou ràpidament

### 8.2. Depurar problemes

1. **Panell de problemes** (`Ctrl+Shift+M` o `Cmd+Shift+M`) per veure tots els errors alhora
2. **Anar a la definició** (`Ctrl-clic` o `Cmd-clic`) per entendre les interfícies dels processos
3. **Trobar totes les referències** per veure com s'utilitzen els processos
4. **Cerca a tot el projecte** per trobar patrons o problemes similars

### 8.3. Refactorització i millora

1. **Cerca a tot el projecte** (`Ctrl+Shift+F` o `Cmd+Shift+F`) per trobar patrons
2. **Formatació automàtica** (`Shift+Alt+F` o `Shift+Option+F`) per mantenir la consistència
3. **Plegament de codi** per centrar-se en l'estructura
4. **Integració amb Git** per fer un seguiment dels canvis

---

## Resum

Ara heu fet un recorregut ràpid per les funcionalitats de l'IDE de VS Code per al desenvolupament amb Nextflow. Aquestes eines us faran significativament més productius:

- **Reduint errors** mitjançant la comprovació de sintaxi en temps real
- **Accelerant el desenvolupament** amb l'autocompleció intel·ligent
- **Millorant la navegació** en workflows complexos multifitxer
- **Mantenint la qualitat** mitjançant una formatació consistent
- **Millorant la comprensió** mitjançant el ressaltat avançat i la visualització de l'estructura

No esperem que recordeu tot, però ara que sabeu que existeixen aquestes funcionalitats, podreu trobar-les quan les necessiteu. A mesura que continueu desenvolupant workflows amb Nextflow, aquestes funcionalitats de l'IDE es tornaran una segona naturalesa, permetent-vos centrar-vos en escriure codi d'alta qualitat en lloc de lluitar amb la sintaxi i l'estructura.

### Què segueix?

Apliqueu aquestes habilitats de l'IDE mentre treballeu en altres mòduls de formació, per exemple:

- **[nf-test](nf-test.md)**: Creeu suites de proves exhaustives per als vostres workflows
- **[Hello nf-core](../../hello_nf-core/)**: Construïu pipelines de qualitat de producció amb estàndards de la comunitat

El veritable poder d'aquestes funcionalitats de l'IDE emergeix quan treballeu en projectes més grans i complexos. Comenceu a incorporar-les al vostre flux de treball gradualment: en poques sessions, es tornaran una segona naturalesa i transformaran la manera com abordeu el desenvolupament amb Nextflow.

Des de detectar errors abans que us alentissin fins a navegar per bases de codi complexes amb facilitat, aquestes eines us convertiran en un desenvolupador més segur i eficient.

Bon codi!
