# Entorn de Desenvolupament

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducció assistida per IA - [més informació i suggeriments](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Els Entorns de Desenvolupament Integrats (IDE) moderns poden transformar dramàticament la vostra experiència de desenvolupament amb Nextflow. Aquesta missió secundària se centra específicament en aprofitar VS Code i la seva extensió de Nextflow per escriure codi més ràpidament, detectar errors de manera primerenca i navegar per workflows complexos de manera eficient.

!!! note "Això no és un tutorial tradicional"

    A diferència d'altres mòduls de formació, aquesta guia està organitzada com una col·lecció de consells ràpids, trucs i exemples pràctics en lloc d'un tutorial pas a pas. Cada secció es pot explorar de manera independent segons els vostres interessos i necessitats actuals de desenvolupament. Sentiu-vos lliures de saltar d'una secció a una altra i centrar-vos en les funcionalitats que seran més útils immediatament per al vostre desenvolupament de workflows.

## Què hauríeu de saber primer

Aquesta guia assumeix que heu completat el curs de formació [Hello Nextflow](../hello_nextflow/) i que us sentiu còmodes amb els conceptes fonamentals de Nextflow, incloent:

- **Estructura bàsica de workflow**: Comprendre processos, workflows i com es connecten entre ells
- **Operacions de canal**: Crear canals, passar dades entre processos i utilitzar operadors bàsics
- **Mòduls i organització**: Crear mòduls reutilitzables i utilitzar declaracions include
- **Conceptes bàsics de configuració**: Utilitzar `nextflow.config` per a paràmetres, directives de procés i perfils

## Què aprendreu aquí

Aquesta guia se centra en **funcionalitats de productivitat de l'IDE** que us faran un desenvolupador de Nextflow més eficient:

- **Ressaltat de sintaxi avançat**: Entendre què us mostra VS Code sobre l'estructura del vostre codi
- **Autocompletat intel·ligent**: Aprofitar suggeriments contextuals per escriure codi més ràpidament
- **Detecció d'errors i diagnòstics**: Detectar errors de sintaxi abans d'executar el vostre workflow
- **Navegació de codi**: Moure's ràpidament entre processos, mòduls i definicions
- **Format i organització**: Mantenir un estil de codi consistent i llegible
- **Desenvolupament assistit per IA** (opcional): Utilitzar eines d'IA modernes integrades amb el vostre IDE

!!! info "Per què funcionalitats d'IDE ara?"

    Probablement ja heu estat utilitzant VS Code durant el curs [Hello Nextflow](../hello_nextflow/), però vam mantenir el focus en aprendre els fonaments de Nextflow en lloc de les funcionalitats de l'IDE. Ara que us sentiu còmodes amb conceptes bàsics de Nextflow com processos, workflows, canals i mòduls, esteu preparats per aprofitar les funcionalitats sofisticades de l'IDE que us faran un desenvolupador més eficient.

    Penseu en això com a "pujar de nivell" el vostre entorn de desenvolupament - el mateix editor que heu estat utilitzant té capacitats molt més potents que esdevenen realment valuoses un cop enteneu amb què us estan ajudant.

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

!!! note "Sobre els Fitxers d'Exemple"

    - `basic_workflow.nf` és un workflow bàsic funcional que podeu executar i modificar
    - `complex_workflow.nf` està dissenyat només per a il·lustració per demostrar funcionalitats de navegació - pot ser que no s'executi correctament però mostra una estructura de workflow multi-fitxer realista

### Dreceres de Teclat

Algunes de les funcionalitats d'aquesta guia utilitzaran dreceres de teclat opcionals. És possible que estigueu accedint a aquest material via GitHub Codespaces al navegador, i en aquest cas de vegades les dreceres no funcionaran com s'espera perquè s'utilitzen per a altres coses al vostre sistema.

Si esteu executant VS Code localment, com probablement fareu quan estigueu escrivint workflows realment, les dreceres funcionaran com es descriu.

Si utilitzeu un Mac, algunes (no totes) dreceres de teclat utilitzaran "cmd" en lloc de "ctrl", i ho indicarem al text com `Ctrl/Cmd`.

### 0.1. Instal·lació de l'Extensió de Nextflow

!!! note "Ja Utilitzeu Devcontainers?"

    Si esteu treballant a **GitHub Codespaces** o utilitzant un **devcontainer local**, l'extensió de Nextflow probablement ja està instal·lada i configurada per a vosaltres. Podeu saltar els passos d'instal·lació manual a continuació i procedir directament a explorar les funcionalitats de l'extensió.

Per instal·lar l'extensió manualment:

1. Obriu VS Code
2. Aneu a la vista d'Extensions clicant la icona d'extensions a l'esquerra: ![icona d'extensions](img/extensions_icon.png) (drecera `Ctrl/Cmd+Shift+X` si esteu executant VSCode localment)
3. Cerqueu "Nextflow"
4. Instal·leu l'extensió oficial de Nextflow

![Instal·lar l'Extensió de Nextflow](img/install_extension.png)

### 0.2. Disposició de l'Espai de Treball

Com que heu estat utilitzant VS Code durant Hello Nextflow, ja esteu familiaritzats amb els conceptes bàsics. Així és com organitzar el vostre espai de treball de manera eficient per a aquesta sessió:

- **Àrea d'Editor**: Per visualitzar i editar fitxers. Podeu dividir-la en múltiples panells per comparar fitxers costat a costat.
- **Explorador de Fitxers** cliqueu (![icona d'explorador de fitxers](img/files_icon.png)) (`Ctrl/Cmd+Shift+E`): Els fitxers i carpetes locals al vostre sistema. Manteniu-lo obert a l'esquerra per navegar entre fitxers
- **Terminal Integrat** (`Ctrl+Shift+` accent greu tant per a Windows com per a MacOS): Un terminal per interactuar amb l'ordinador a la part inferior. Utilitzeu-lo per executar Nextflow o altres comandes.
- **Panell de Problemes** (`Ctrl+Shift+M`): VS Code mostrarà aquí qualsevol error i problema que detecti. Això és útil per ressaltar problemes d'una ullada.

Podeu arrossegar panells o amagar-los (`Ctrl/Cmd+B` per alternar la barra lateral) per personalitzar la vostra disposició mentre treballem amb els exemples.

### Conclusió

Teniu VS Code configurat amb l'extensió de Nextflow i enteneu la disposició de l'espai de treball per a un desenvolupament eficient.

### Què segueix?

Apreneu com el ressaltat de sintaxi us ajuda a entendre l'estructura del codi Nextflow d'una ullada.

---

## 1. Ressaltat de Sintaxi i Estructura del Codi

Ara que el vostre espai de treball està configurat, explorem com el ressaltat de sintaxi de VS Code us ajuda a llegir i escriure codi Nextflow de manera més efectiva.

### 1.1. Elements de Sintaxi de Nextflow

Obriu `basic_workflow.nf` per veure el ressaltat de sintaxi en acció:

![Mostra de Sintaxi](img/syntax_showcase.png)

Observeu com VS Code ressalta:

- **Paraules clau** (`process`, `workflow`, `input`, `output`, `script`) amb colors diferents
- **Literals de cadena** i **paràmetres** amb estils diferents
- **Comentaris** amb un color apagat
- **Variables** i **crides a funcions** amb l'èmfasi apropiat
- **Blocs de codi** amb guies d'indentació adequades

!!! note "Colors Dependents del Tema"

    Els colors específics que veieu dependran del vostre tema de VS Code (mode fosc/clar), configuració de colors i qualsevol personalització que hàgiu fet. L'important és que diferents elements de sintaxi es distingeixen visualment entre ells, fent que l'estructura del codi sigui més fàcil d'entendre independentment de l'esquema de colors escollit.

### 1.2. Comprensió de l'Estructura del Codi

El ressaltat de sintaxi us ajuda a identificar ràpidament:

- **Límits de procés**: Distinció clara entre diferents processos
- **Blocs d'entrada/sortida**: Fàcil de detectar definicions de flux de dades
- **Blocs de script**: Les comandes reals que s'executen
- **Operacions de canal**: Passos de transformació de dades
- **Directives de configuració**: Configuracions específiques de procés

Aquesta organització visual esdevé invaluable quan es treballa amb workflows complexos que contenen múltiples processos i fluxos de dades intrincats.

### Conclusió

Enteneu com el ressaltat de sintaxi de VS Code us ajuda a llegir l'estructura del codi Nextflow i identificar diferents elements del llenguatge per a un desenvolupament més ràpid.

### Què segueix?

Apreneu com l'autocompletat intel·ligent accelera l'escriptura de codi amb suggeriments contextuals.

---

## 2. Autocompletat Intel·ligent

Les funcionalitats d'autocompletat de VS Code us ajuden a escriure codi més ràpidament i amb menys errors suggerint opcions apropiades segons el context.

### 2.1. Suggeriments Contextuals

Les opcions d'autocompletat varien segons on us trobeu al vostre codi:

#### Operacions de Canal

Obriu `basic_workflow.nf` de nou i proveu d'escriure `channel.` al bloc de workflow:

![Autocompletat de canal](img/autocomplete_channel.png)

Veureu suggeriments per a:

- `fromPath()` - Crear canal des de rutes de fitxer
- `fromFilePairs()` - Crear canal des de fitxers aparellats
- `of()` - Crear canal des de valors
- `fromSRA()` - Crear canal des d'accessos SRA
- I molts més...

Això us ajuda a trobar ràpidament la factory de canal correcta a utilitzar sense necessitat de recordar noms de mètode exactes.

També podeu descobrir els operadors disponibles per aplicar als canals. Per exemple, escriviu `FASTQC.out.html.` per veure les operacions disponibles:

![Autocompletat d'operacions de canal](img/autocomplete_operators.png)

#### Directives de Procés

Dins d'un bloc de script de procés, escriviu `task.` per veure les propietats d'execució disponibles:

![Autocompletat de propietats de tasca](img/autocomplete_task.png)

#### Configuració

Obriu nextflow.config i escriviu `process.` a qualsevol lloc per veure les directives de procés disponibles:

![Autocompletat de configuració](img/autocomplete_config.png)

Veureu suggeriments per a:

- `executor`
- `memory`
- `cpus`

Això estalvia temps quan es configuren processos i funciona a través de diferents àmbits de configuració. Per exemple, proveu d'escriure `docker.` per veure opcions de configuració específiques de Docker.

### Conclusió

Podeu utilitzar l'autocompletat intel·ligent de VS Code per descobrir operacions de canal disponibles, directives de procés i opcions de configuració sense memoritzar sintaxi.

### Què segueix?

Apreneu com la detecció d'errors en temps real us ajuda a detectar problemes abans d'executar el vostre workflow, simplement llegint el codi.

## 3. Detecció d'Errors i Diagnòstics

La detecció d'errors en temps real de VS Code us ajuda a detectar problemes abans d'executar el vostre workflow.

### 3.1. Detecció d'Errors de Sintaxi

Creem un error deliberat per veure la detecció en acció. Obriu `basic_workflow.nf` i canvieu el nom del procés de `FASTQC` a `FASTQ` (o qualsevol altre nom invàlid). VS Code ressaltarà immediatament l'error al bloc de workflow amb una línia ondulada vermella:

![Subratllat d'error](img/error_underline.png)

### 3.2. Panell de Problemes

Més enllà del ressaltat d'errors individuals, VS Code proporciona un Panell de Problemes centralitzat que agrega tots els errors, advertències i missatges d'informació a través del vostre espai de treball. Obriu-lo amb `Ctrl/Cmd+Shift+M` i utilitzeu la icona de filtre per mostrar només errors rellevants al fitxer actual:

![Filtrar el panell de problemes](img/active_file.png)

Cliqueu a qualsevol problema per saltar directament a la línia problemàtica

![Panell de Problemes](img/problems_panel.png)

Corregiu l'error canviant el nom del procés de nou a `FASTQC`.

### 3.3. Patrons d'Error Comuns

Els errors comuns a la sintaxi de Nextflow inclouen:

- **Claudàtors que falten**: `{` o `}` sense parella
- **Blocs incomplets**: Seccions requerides que falten als processos
- **Sintaxi invàlida**: DSL de Nextflow mal format
- **Errors tipogràfics a paraules clau**: Directives de procés mal escrites
- **Desajustos de canal**: Incompatibilitats de tipus

El servidor de llenguatge de Nextflow ressalta aquests problemes al Panell de Problemes. Podeu revisar-los de manera primerenca per evitar errors de sintaxi mentre executeu un pipeline.

### Conclusió

Podeu utilitzar la detecció d'errors de VS Code i el Panell de Problemes per detectar errors de sintaxi i problemes abans d'executar el vostre workflow, estalviant temps i evitant frustracions.

### Què segueix?

Apreneu com navegar eficientment entre processos, mòduls i definicions en workflows complexos.

---

## 4. Navegació de Codi i Gestió de Símbols

La navegació eficient és crucial quan es treballa amb workflows complexos que abasten múltiples fitxers. Per entendre això, substituïu la definició del procés a `basic_workflow.nf` amb una importació del mòdul que us hem proporcionat:

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

### 4.1. Anar a la Definició

Si passeu el ratolí per sobre d'un nom de procés com `FASTQC`, veureu una finestra emergent amb la interfície del mòdul (entrades i sortides):

![Anar a la definició](img/syntax.png)

Aquesta funcionalitat és particularment valuosa quan s'escriuen workflows, ja que us permet entendre la interfície del mòdul sense obrir directament el fitxer del mòdul.

Podeu navegar ràpidament a qualsevol definició de procés, mòdul o variable utilitzant **Ctrl/Cmd-clic**. Passeu el ratolí per sobre de l'enllaç al fitxer del mòdul a la part superior de l'script, i seguiu l'enllaç com es suggereix:

![Seguir enllaç](img/follow_link.png)

El mateix funciona per als noms de procés. Torneu a `basic_workflow.nf` i proveu-ho amb el nom del procés `FASTQC` al bloc de workflow. Això us enllaça directament al nom del procés (que és el mateix que el fitxer del mòdul en aquest exemple, però podria estar a mig camí d'un fitxer molt més gran).

Per tornar a on éreu, utilitzeu **Alt+←** (o **Ctrl+-** a Mac). Aquesta és una manera potent d'explorar codi sense perdre la vostra posició.

Ara explorem la navegació en un workflow més complex utilitzant `complex_workflow.nf` (el fitxer només d'il·lustració mencionat anteriorment). Aquest workflow conté múltiples processos definits en fitxers de mòdul separats, així com alguns en línia. Tot i que les estructures multi-fitxer complexes poden ser difícils de navegar manualment, la capacitat de saltar a definicions fa l'exploració molt més manejable.

1. Obriu `complex_workflow.nf`
2. Navegueu a definicions de mòdul
3. Utilitzeu **Alt+←** (o **Ctrl+-**) per navegar enrere
4. Navegueu al nom del procés `FASTQC` al bloc de workflow. Això us enllaça directament al nom del procés (que és el mateix que el fitxer del mòdul en aquest exemple, però podria estar a mig camí d'un fitxer molt més gran).
5. Navegueu enrere de nou
6. Navegueu al procés `TRIM_GALORE` al bloc de workflow. Aquest està definit en línia, així que no us portarà a un fitxer separat, però encara us mostrarà la definició del procés, i encara podreu navegar enrere a on éreu.

### 4.2. Navegació de Símbols

Amb `complex_workflow.nf` encara obert, podeu obtenir una visió general de tots els símbols al fitxer escrivint `@` a la barra de cerca a la part superior de VSCode (la drecera de teclat és `Ctrl/Cmd+Shift+O`, però pot ser que no funcioni a Codespaces). Això obre el panell de navegació de símbols, que llista tots els símbols al fitxer actual:

![Navegació de símbols](img/symbols.png)

Això mostra:

- Totes les definicions de procés
- Definicions de workflow (hi ha dos workflows definits en aquest fitxer)
- Definicions de funció

Comenceu a escriure per filtrar resultats.

### 4.3. Trobar Totes les Referències

Entendre on s'utilitza un procés o variable a través del vostre codi pot ser molt útil. Per exemple, si voleu trobar totes les referències al procés `FASTQC`, comenceu navegant a la seva definició. Podeu fer-ho obrint `modules/fastqc.nf` directament, o utilitzant la funcionalitat de navegació ràpida de VS Code amb `Ctrl/Cmd-clic` com vam fer anteriorment. Un cop a la definició del procés, feu clic dret sobre el nom del procés `FASTQC` i seleccioneu "Find All References" del menú contextual per veure totes les instàncies on s'utilitza.

![Trobar referències](img/references.png)

Aquesta funcionalitat mostra totes les instàncies on es fa referència a `FASTQC` dins del vostre espai de treball, incloent el seu ús als dos workflows diferents. Aquesta visió és crucial per avaluar l'impacte potencial de modificacions al procés `FASTQC`.

### 4.4. Panell d'Esquema

El panell d'Esquema, ubicat a la barra lateral de l'Explorador (cliqueu ![icona d'Explorador](img/files_icon.png)), proporciona una visió general convenient de tots els símbols al vostre fitxer actual. Aquesta funcionalitat us permet navegar ràpidament i gestionar l'estructura del vostre codi mostrant funcions, variables i altres elements clau en una vista jeràrquica.

![Panell d'esquema](img/outline.png)

Utilitzeu el panell d'Esquema per navegar ràpidament a diferents parts del vostre codi sense utilitzar el navegador de fitxers.

### 4.5. Visualització DAG

L'extensió de Nextflow de VS Code pot visualitzar el vostre workflow com un Graf Acíclic Dirigit (DAG). Això us ajuda a entendre el flux de dades i les dependències entre processos. Obriu `complex_workflow.nf` i cliqueu el botó "Preview DAG" sobre `workflow {` (el segon bloc `workflow` en aquest fitxer):

![Vista prèvia DAG](img/dag_preview.png)

Aquest és només el workflow 'entry', però també podeu previsualitzar el DAG per als workflows interns clicant el botó "Preview DAG" sobre el workflow `RNASEQ_PIPELINE {` més amunt:

![Vista prèvia DAG workflow intern](img/dag_preview_inner.png)

Per a aquest workflow, podeu utilitzar els nodes del DAG per navegar a les definicions de procés corresponents al codi. Cliqueu a un node, i us portarà a la definició de procés rellevant a l'editor. Particularment quan un workflow creix a una mida gran, això pot ajudar-vos realment a navegar pel codi i entendre com estan connectats els processos.

### Conclusió

Podeu navegar workflows complexos eficientment utilitzant anar-a-definició, cerca de símbols, trobar referències i visualització DAG per entendre l'estructura del codi i les dependències.

### Què segueix?

Apreneu com treballar efectivament a través de múltiples fitxers interconnectats en projectes Nextflow més grans.

## 5. Treballar a Través de Múltiples Fitxers

El desenvolupament real de Nextflow implica treballar amb múltiples fitxers interconnectats. Explorem com VS Code us ajuda a gestionar projectes complexos de manera eficient.

### 5.1. Navegació Ràpida de Fitxers

Amb `complex_workflow.nf` obert, notareu que importa diversos mòduls. Practiquem la navegació ràpida entre ells.

Premeu **Ctrl+P** (o **Cmd+P**) i comenceu a escriure "fast":

VS Code us mostrarà fitxers coincidents. Seleccioneu `modules/fastqc.nf` per saltar-hi instantàniament. Això és molt més ràpid que clicar a través de l'explorador de fitxers quan sabeu aproximadament quin fitxer esteu buscant.

Proveu això amb altres patrons:

- Escriviu "star" per trobar el fitxer del mòdul d'alineament STAR (`star.nf`)
- Escriviu "utils" per trobar el fitxer de funcions d'utilitat (`utils.nf`)
- Escriviu "config" per saltar a fitxers de configuració (`nextflow.config`)

### 5.2. Editor Dividit per a Desenvolupament Multi-fitxer

Quan es treballa amb mòduls, sovint necessiteu veure tant el workflow principal com les definicions de mòdul simultàniament. Configurem això:

1. Obriu `complex_workflow.nf`
2. Obriu `modules/fastqc.nf` en una pestanya nova
3. Feu clic dret a la pestanya `modules/fastqc.nf` i seleccioneu "Split Right"
4. Ara podeu veure ambdós fitxers costat a costat

![Editor dividit](img/split_editor.png)

Això és invaluable quan:

- Es comproven interfícies de mòdul mentre s'escriuen crides de workflow, i la vista prèvia no és suficient
- Es comparen processos similars a través de diferents mòduls
- Es depura el flux de dades entre workflow i mòduls

### 5.3. Cerca a Tot el Projecte

De vegades necessiteu trobar on s'utilitzen patrons específics a través de tot el vostre projecte. Premeu `Ctrl/Cmd+Shift+F` per obrir el panell de cerca.

Proveu de cercar `publishDir` a través de l'espai de treball:

![Cerca de projecte](img/project_search.png)

Això us mostra cada fitxer que utilitza directoris de publicació, ajudant-vos a:

- Entendre patrons d'organització de sortida
- Trobar exemples de directives específiques
- Assegurar consistència a través de mòduls

### Conclusió

Podeu gestionar projectes multi-fitxer complexos utilitzant navegació ràpida de fitxers, editors dividits i cerca a tot el projecte per treballar eficientment a través de workflows i mòduls.

### Què segueix?

Apreneu com les funcionalitats de format i manteniment de codi mantenen els vostres workflows organitzats i llegibles.

---

## 6. Format i Manteniment de Codi

El format adequat del codi és essencial no només per a l'estètica sinó també per millorar la llegibilitat, la comprensió i la facilitat d'actualització de workflows complexos.

### 6.1. Format Automàtic en Acció

Obriu `basic_workflow.nf` i desorganitzeu deliberadament el format:

- Elimineu alguna indentació: Ressalteu tot el document i premeu `shift+tab` moltes vegades per eliminar tantes indentacions com sigui possible.
- Afegiu espais extra en llocs aleatoris: a la declaració `channel.fromPath`, afegiu 30 espais després del `(`.
- Trenqueu algunes línies de manera estranya: Afegiu una línia nova entre l'operador `.view {` i la cadena `Processing sample:` però no afegiu una línia nova corresponent abans del parèntesi de tancament `}`.

Ara premeu `Shift+Alt+F` (o `Shift+Option+F` a MacOS) per formatar automàticament:

VS Code immediatament:

- Corregeix la indentació per mostrar l'estructura del procés clarament
- Alinea elements similars de manera consistent
- Elimina espais en blanc innecessaris
- Manté salts de línia llegibles

Tingueu en compte que el format automàtic pot no resoldre tots els problemes d'estil de codi. El servidor de llenguatge de Nextflow pretén mantenir el vostre codi endreçat, però també respecta les vostres preferències personals en certes àrees. Per exemple, si elimineu la indentació dins del bloc `script` d'un procés, el formatador ho deixarà tal com està, ja que podríeu preferir intencionadament aquest estil.

Actualment, no hi ha una aplicació estricta d'estil per a Nextflow, així que el servidor de llenguatge ofereix certa flexibilitat. No obstant això, aplicarà consistentment regles de format al voltant de definicions de mètodes i funcions per mantenir la claredat.

### 6.2. Funcionalitats d'Organització de Codi

#### Comentaris Ràpids

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

- Desactivar temporalment parts de workflows durant el desenvolupament
- Afegir comentaris explicatius a operacions de canal complexes
- Documentar seccions de workflow

Utilitzeu **Ctrl+/** (o **Cmd+/**) de nou per descomentar el codi.

#### Plegat de Codi per a Visió General

A `complex_workflow.nf`, observeu les petites fletxes al costat de les definicions de procés. Cliqueu-les per plegar (col·lapsar) processos:

![Plegat de codi](img/code_folding.png)

Això us dóna una visió general d'alt nivell de l'estructura del vostre workflow sense perdre's en detalls d'implementació.

#### Coincidència de Claudàtors

Col·loqueu el cursor al costat de qualsevol claudàtor `{` o `}` i VS Code ressalta el claudàtor coincident. Utilitzeu **Ctrl+Shift+\\** (o **Cmd+Shift+\\**) per saltar entre claudàtors coincidents.

Això és crucial per a:

- Entendre límits de procés
- Trobar claudàtors que falten o sobren
- Navegar estructures de workflow niuades

#### Selecció i Edició Multi-línia

Per editar múltiples línies simultàniament, VS Code ofereix capacitats multi-cursor potents:

- **Selecció multi-línia**: Manteniu **Ctrl+Alt** (o **Cmd+Option** per a MacOS) i utilitzeu les tecles de fletxa per seleccionar múltiples línies
- **Indentació multi-línia**: Seleccioneu múltiples línies i utilitzeu **Tab** per indentar o **Shift+Tab** per desindentar blocs sencers

Això és particularment útil per a:

- Indentar blocs de procés sencers de manera consistent
- Afegir comentaris a múltiples línies alhora
- Editar definicions de paràmetres similars a través de múltiples processos

### Conclusió

Podeu mantenir codi net i llegible utilitzant format automàtic, funcionalitats de comentaris, plegat de codi, coincidència de claudàtors i edició multi-línia per organitzar workflows complexos de manera eficient.

### Què segueix?

Apreneu com VS Code s'integra amb el vostre workflow de desenvolupament més ampli més enllà de simplement editar codi.

---

## 7. Integració del Workflow de Desenvolupament

VS Code s'integra bé amb el vostre workflow de desenvolupament més enllà de simplement editar codi.

### 7.1. Integració de Control de Versions

!!! note "Codespaces i Integració Git"

    Si esteu treballant a **GitHub Codespaces**, algunes funcionalitats d'integració Git poden no funcionar com s'espera, particularment dreceres de teclat per a Control de Codi. També podeu haver declinat obrir el directori com a repositori Git durant la configuració inicial, que està bé per a propòsits de formació.

Si el vostre projecte és un repositori git (com aquest ho és), VS Code mostra:

- Fitxers modificats amb indicadors de color
- Estat de Git a la barra d'estat
- Vistes de diferències en línia
- Capacitats de commit i push

Obriu el panell de Control de Codi utilitzant el botó de control de codi (![icona de Control de codi](img/source_control_icon.png)) (`Ctrl+Shift+G` o `Cmd+Shift+G` si esteu treballant amb VSCode localment) per veure canvis de git i preparar commits directament a l'editor.

![Panell de Control de Codi](img/source_control.png)

### 7.2. Execució i Inspecció de Workflows

Executem un workflow i després inspeccionem els resultats. Al terminal integrat (`Ctrl+Shift+` accent greu tant a Windows com a MacOS), executeu el workflow bàsic:

```bash title="Run the basic workflow"
nextflow run basic_workflow.nf --input data/sample_data.csv --output_dir results
```

Mentre el workflow s'executa, veureu sortida en temps real al terminal. Després de completar-se, podeu utilitzar VS Code per inspeccionar resultats sense sortir del vostre editor:

1. **Navegueu a directoris de treball**: Utilitzeu l'explorador de fitxers o el terminal per navegar `.nextflow/work`
2. **Obriu fitxers de registre**: Cliqueu a rutes de fitxers de registre a la sortida del terminal per obrir-los directament a VS Code
3. **Inspeccioneu sortides**: Navegueu directoris de resultats publicats a l'explorador de fitxers
4. **Visualitzeu informes d'execució**: Obriu informes HTML directament a VS Code o al vostre navegador

Això manté tot en un sol lloc en lloc de canviar entre múltiples aplicacions.

### Conclusió

Podeu integrar VS Code amb control de versions i execució de workflows per gestionar tot el vostre procés de desenvolupament des d'una única interfície.

### Què segueix?

Veieu com totes aquestes funcionalitats d'IDE treballen juntes al vostre workflow de desenvolupament diari.

---

## 8. Recapitulació i notes ràpides

Aquí teniu algunes notes ràpides sobre cadascuna de les funcionalitats d'IDE discutides anteriorment:

### 8.1. Començar una Nova Funcionalitat

1. **Obertura ràpida de fitxers** (`Ctrl+P` o `Cmd+P`) per trobar mòduls existents rellevants
2. **Editor dividit** per veure processos similars costat a costat
3. **Navegació de símbols** (`Ctrl+Shift+O` o `Cmd+Shift+O`) per entendre l'estructura del fitxer
4. **Autocompletat** per escriure codi nou ràpidament

### 8.2. Depuració de Problemes

1. **Panell de problemes** (`Ctrl+Shift+M` o `Cmd+Shift+M`) per veure tots els errors alhora
2. **Anar a la definició** (`Ctrl-clic` o `Cmd-clic`) per entendre interfícies de procés
3. **Trobar totes les referències** per veure com s'utilitzen els processos
4. **Cerca a tot el projecte** per trobar patrons o problemes similars

### 8.3. Refactorització i Millora

1. **Cerca a tot el projecte** (`Ctrl+Shift+F` o `Cmd+Shift+F`) per trobar patrons
2. **Format automàtic** (`Shift+Alt+F` o `Shift+Option+F`) per mantenir consistència
3. **Plegat de codi** per centrar-se en l'estructura
4. **Integració Git** per fer seguiment de canvis

---

## Resum

Ara heu tingut un recorregut ràpid per les funcionalitats d'IDE de VS Code per al desenvolupament de Nextflow. Aquestes eines us faran significativament més productius:

- **Reduint errors** mitjançant comprovació de sintaxi en temps real
- **Accelerant el desenvolupament** amb autocompletat intel·ligent
- **Millorant la navegació** en workflows multi-fitxer complexos
- **Mantenint la qualitat** mitjançant format consistent
- **Millorant la comprensió** mitjançant ressaltat avançat i visualització d'estructura

No esperem que recordeu tot, però ara sabeu que aquestes funcionalitats existeixen i podreu trobar-les quan les necessiteu. A mesura que continueu desenvolupant workflows de Nextflow, aquestes funcionalitats d'IDE esdevindran una segona naturalesa, permetent-vos centrar-vos en escriure codi d'alta qualitat en lloc de lluitar amb la sintaxi i l'estructura.

### Què segueix?

Apliqueu aquestes habilitats d'IDE mentre treballeu amb altres mòduls de formació, per exemple:

- **[nf-test](nf-test.md)**: Creeu suites de proves completes per als vostres workflows
- **[Hello nf-core](../../hello_nf-core/)**: Construïu pipelines de qualitat de producció amb estàndards de la comunitat

El veritable poder d'aquestes funcionalitats d'IDE emergeix a mesura que treballeu en projectes més grans i complexos. Comenceu a incorporar-les al vostre workflow gradualment—en poques sessions, esdevindran una segona naturalesa i transformaran com abordeu el desenvolupament de Nextflow.

Des de detectar errors abans que us alenteixen fins a navegar per codis complexos amb facilitat, aquestes eines us faran un desenvolupador més confiat i eficient.

Feliç programació!
